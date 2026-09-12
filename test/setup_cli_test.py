"""Tests for ``viralconseq setup``."""

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml
from click.testing import CliRunner

from viralconseq.config_generator import ConfigGenerator
from viralconseq.setup_cli import (
    _ALL_PIPELINES,
    _PIPELINE_TO_WORKFLOW,
    _collect_env_yamls,
    _scripts_dir,
    setup,
)


class Test_ConfigGeneratorSkeleton(unittest.TestCase):
    """``ConfigGenerator.write_skeleton`` writes loadable YAML with the
    required keys for each pipeline x data_type combination."""

    def _write_and_load(self, pipeline: str, data_type: str) -> dict:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "skel.yaml"
            placeholder = Path(tmp) / "ph"
            placeholder.mkdir()
            ConfigGenerator.write_skeleton(pipeline, data_type, str(path), str(placeholder))
            with open(path) as fh:
                return yaml.safe_load(fh)

    def test_consensus_illumina_has_required_keys(self):
        cfg = self._write_and_load("consensus", "illumina")
        for k in ("samples", "data", "output", "threads", "reference", "scheme", "adapters"):
            self.assertIn(k, cfg, f"missing key: {k}")
        self.assertEqual(cfg["data"], "illumina")

    def test_consensus_nanopore_has_required_keys(self):
        cfg = self._write_and_load("consensus", "nanopore")
        for k in ("samples", "data", "output", "threads", "reference", "clair3_model"):
            self.assertIn(k, cfg, f"missing key: {k}")
        self.assertEqual(cfg["data"], "nanopore")

    def test_unknown_pipeline_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "skel.yaml"
            placeholder = Path(tmp) / "ph"
            placeholder.mkdir()
            with self.assertRaises(ValueError):
                ConfigGenerator.write_skeleton("bogus", "illumina", str(path), str(placeholder))

    def test_optional_features_enabled_consensus_illumina(self):
        """``--run-isnv`` must be on so the LoFreq branch is in the DAG and
        ``viralconseq setup`` materializes ``envs/consensus.yaml``."""
        cfg = self._write_and_load("consensus", "illumina")
        self.assertTrue(cfg["run_isnv"])


class Test_CollectEnvYamls(unittest.TestCase):
    """``_collect_env_yamls`` enumerates the conda envs declared in a
    workflow and every rule module it includes."""

    def test_consensus_illumina_finds_known_envs(self):
        smk = _scripts_dir() / "consensus_illumina.smk"
        yamls = _collect_env_yamls(smk)
        # The package ships five env YAMLs; consensus_illumina pulls a subset.
        self.assertTrue(yamls, "no envs found")
        # qc.yaml is the env that failed first in the original bug report.
        self.assertIn("qc.yaml", yamls)

    def test_all_pipelines_resolve_to_existing_smk(self):
        scripts = _scripts_dir()
        for name, (_, _, smk_name) in _PIPELINE_TO_WORKFLOW.items():
            with self.subTest(pipeline=name):
                self.assertTrue(
                    (scripts / smk_name).is_file(),
                    f"workflow not found for {name}: {smk_name}",
                )


class Test_SetupCli(unittest.TestCase):
    """``viralconseq setup`` honours its flags and forwards them to
    Snakemake correctly."""

    def setUp(self):
        self.runner = CliRunner()

    def test_dry_run_lists_envs(self):
        result = self.runner.invoke(
            setup,
            ["--pipelines", "consensus-illumina", "--dry-run"],
            catch_exceptions=False,
        )
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertIn("consensus-illumina", result.output)
        self.assertIn("qc.yaml", result.output)
        # Regression: consensus.yaml (LoFreq) must appear for consensus-illumina.
        # Previously the skeleton config left run_isnv=False, so the
        # detect_isnv rule was pruned out of the DAG and consensus.yaml was
        # never materialized by `viralconseq setup`. Users who later passed
        # --run-isnv hit dynamic env creation on the hot path.
        self.assertIn("consensus.yaml", result.output)

    def test_pipeline_choices_are_the_two_consensus_workflows(self):
        self.assertEqual(_ALL_PIPELINES, ["consensus-illumina", "consensus-nanopore"])
        self.assertEqual(
            sorted(_PIPELINE_TO_WORKFLOW),
            ["consensus-illumina", "consensus-nanopore"],
        )

    def test_metagenomics_pipeline_is_rejected(self):
        result = self.runner.invoke(setup, ["--pipelines", "meta-illumina", "--dry-run"])
        self.assertEqual(result.exit_code, 2, result.output)
        self.assertIn("meta-illumina", result.output)

    def test_dry_run_all_pipelines(self):
        result = self.runner.invoke(setup, ["--dry-run"], catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        for name in _ALL_PIPELINES:
            self.assertIn(name, result.output)

    def test_conda_prefix_forwarded_to_snakemake(self):
        with (
            patch("viralconseq.setup_cli.snakemake", return_value=True) as mock_snake,
            tempfile.TemporaryDirectory() as tmp,
        ):
            result = self.runner.invoke(
                setup,
                [
                    "--pipelines",
                    "consensus-illumina",
                    "--conda-prefix",
                    tmp,
                ],
                catch_exceptions=False,
            )
        self.assertEqual(result.exit_code, 0, result.output)
        mock_snake.assert_called_once()
        kwargs = mock_snake.call_args.kwargs
        self.assertEqual(kwargs["conda_prefix"], tmp)
        self.assertTrue(kwargs["use_conda"])
        self.assertTrue(kwargs["conda_create_envs_only"])

    def test_default_conda_prefix(self):
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop("VIRALCONSEQ_CONDA_PREFIX", None)
            with patch("viralconseq.setup_cli.snakemake", return_value=True) as mock_snake:
                result = self.runner.invoke(
                    setup, ["--pipelines", "consensus-illumina"], catch_exceptions=False
                )
        self.assertEqual(result.exit_code, 0, result.output)
        expected = str(Path.home() / ".cache" / "viralconseq" / "conda-envs")
        self.assertEqual(mock_snake.call_args.kwargs["conda_prefix"], expected)

    def test_env_var_conda_prefix(self):
        with (
            tempfile.TemporaryDirectory() as tmp,
            patch.dict(os.environ, {"VIRALCONSEQ_CONDA_PREFIX": tmp}),
        ):
            with patch("viralconseq.setup_cli.snakemake", return_value=True) as mock_snake:
                result = self.runner.invoke(
                    setup, ["--pipelines", "consensus-illumina"], catch_exceptions=False
                )
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(mock_snake.call_args.kwargs["conda_prefix"], tmp)

    def test_failure_exits_nonzero(self):
        """If Snakemake reports an env-build failure, setup exits non-zero."""
        with (
            patch("viralconseq.setup_cli.snakemake", return_value=False),
            tempfile.TemporaryDirectory() as tmp,
        ):
            result = self.runner.invoke(
                setup,
                ["--pipelines", "consensus-illumina", "--conda-prefix", tmp],
                catch_exceptions=False,
            )
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("FAILED", result.output)


if __name__ == "__main__":
    unittest.main()
