"""Tests for ``viralconseq rerun``: config loading, --set overrides, validation
before Snakemake, workflow selection and the kwargs handed to run_workflow."""

import os
import tempfile
import unittest
from unittest.mock import patch

import yaml
from click.testing import CliRunner

from viralconseq.rerun_cli import apply_overrides, args_from_config, rerun


def _config(**overrides):
    config = {
        "samples": {"sample-a": ["/r/a.fastq.gz"]},
        "data": "nanopore",
        "output": "/out/run1/",
        "threads": 2,
        "threads_total": 4,
        "max_memory_mb": 4096,
        "reference": "/ref.fa",
        "scheme": "NA",
        "minimum_depth": 20,
        "minimum_length": 50,
        "af_threshold": 0.51,
        "chunk_size": 10000,
        "clair3_model": "r941_prom_hac_g360+g422",
        "variant_quality": 20,
        "variant_depth": 10,
        "minimum_map_quality": 30,
        "infer_consensus_sequence_ram": 2,
        "run_viralqc": False,
    }
    config.update(overrides)
    return config


class Test_Overrides(unittest.TestCase):
    def test_values_are_parsed_as_yaml(self):
        config = apply_overrides(
            _config(), ["minimum_depth=30", "af_threshold=0.7", "run_viralqc=true", "scheme='x y'"]
        )
        self.assertEqual(config["minimum_depth"], 30)
        self.assertEqual(config["af_threshold"], 0.7)
        self.assertIs(config["run_viralqc"], True)
        self.assertEqual(config["scheme"], "x y")

    def test_unknown_key_and_malformed_entry_rejected(self):
        from click import BadParameter

        with self.assertRaises(BadParameter):
            apply_overrides(_config(), ["minimum_dept=30"])
        with self.assertRaises(BadParameter):
            apply_overrides(_config(), ["nonsense"])

    def test_args_from_config_targets_the_run_directory(self):
        args = args_from_config(_config(), "c.yml", "/envs")
        self.assertEqual(args["output"], "/out")
        self.assertEqual(args["run_name"], "run1")
        self.assertEqual(args["threads_total"], 4)
        self.assertEqual(args["max_memory_mb"], 4096)
        self.assertEqual(args["data_type"], "nanopore")
        self.assertTrue(os.path.isabs(args["config_file"]))


class Test_RerunCommand(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()
        self._tmp = tempfile.TemporaryDirectory()
        self.config_path = os.path.join(self._tmp.name, "config.yml")

    def tearDown(self):
        self._tmp.cleanup()

    def _write(self, config):
        with open(self.config_path, "w") as fh:
            yaml.safe_dump(config, fh)

    def _invoke(self, *extra, ok=True):
        captured = {}

        def fake_run(workflow, args, **kwargs):
            captured["workflow"] = workflow
            captured["args"] = args
            captured.update(kwargs)
            return ok

        with patch("viralconseq.rerun_cli._orchestrator.run_workflow", side_effect=fake_run):
            result = self.runner.invoke(rerun, [self.config_path, *extra])
        return result, captured

    def test_runs_the_matching_workflow_with_mode_flags(self):
        self._write(_config())
        result, captured = self._invoke("--dry-run", "--keep-going")
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertTrue(captured["workflow"].endswith("consensus_nanopore.smk"))
        self.assertEqual(
            (captured["dryrun"], captured["unlock"], captured["keepgoing"]), (True, False, True)
        )
        self.assertEqual(captured["args"]["threads_total"], 4)
        self.assertFalse(os.path.exists(self.config_path + ".bak"))

    def test_segmented_reference_picks_segmented_workflow(self):
        self._write(_config(reference={"S": "/S.fa", "L": "/L.fa"}))
        _, captured = self._invoke("--dry-run")
        self.assertTrue(captured["workflow"].endswith("consensus_nanopore_segmented.smk"))

    def test_set_rewrites_the_config_and_keeps_a_backup(self):
        self._write(_config())
        result, captured = self._invoke("--set", "minimum_depth=30", "--dry-run")
        self.assertEqual(result.exit_code, 0, result.output)
        with open(self.config_path) as fh:
            self.assertEqual(yaml.safe_load(fh)["minimum_depth"], 30)
        with open(self.config_path + ".bak") as fh:
            self.assertEqual(yaml.safe_load(fh)["minimum_depth"], 20)
        self.assertIn("# --- clair3 ---", open(self.config_path).read())

    def test_invalid_override_rejected_before_snakemake(self):
        self._write(_config())
        result, captured = self._invoke("--set", "af_threshold=abc")
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("configuration_error", result.output)
        self.assertNotIn("workflow", captured)
        with open(self.config_path) as fh:
            self.assertEqual(yaml.safe_load(fh)["af_threshold"], 0.51)  # untouched

    def test_missing_required_key_is_a_configuration_error(self):
        config = _config()
        del config["chunk_size"]
        self._write(config)
        result, captured = self._invoke()
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("chunk_size", result.output)
        self.assertNotIn("workflow", captured)

    def test_failed_workflow_exits_one(self):
        self._write(_config())
        result, _ = self._invoke(ok=False)
        self.assertEqual(result.exit_code, 1)


if __name__ == "__main__":
    unittest.main()
