"""Tests for viralconseq consensus CLI (click-based)."""

import os
import unittest
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from viralconseq.consensus_cli import consensus
from viralconseq.constants import detect_cores


class Test_ConsensusIlluminaCommand(unittest.TestCase):
    """Tests for `viralconseq consensus illumina`."""

    def setUp(self):
        self.runner = CliRunner()
        self._required = [
            "illumina",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ]

    def _invoke(self, extra_args=None):
        args = self._required + (extra_args or [])
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0):
            return self.runner.invoke(consensus, args, catch_exceptions=False)

    def test_required_args_missing_causes_error(self):
        """Missing required args should exit with non-zero code."""
        result = self.runner.invoke(consensus, ["illumina"])
        self.assertNotEqual(result.exit_code, 0)

    def test_required_args_success_with_reference(self):
        result = self._invoke()
        self.assertEqual(result.exit_code, 0, result.output)

    def test_required_args_success_without_reference(self):
        """Reference flags are optional at parse time; validated by core logic."""
        args = [
            "illumina",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
        ]
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0):
            result = self.runner.invoke(consensus, args, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)

    def test_required_args_success_with_segmented_reference(self):
        args = [
            "illumina",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--segmented-reference",
            "S=/path/to/S.fasta",
            "--segmented-reference",
            "L=/path/to/L.fasta",
        ]
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0) as mock_main:
            result = self.runner.invoke(consensus, args, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        called_args = mock_main.call_args[0][0]
        self.assertEqual(
            called_args["segmented_reference"],
            {
                "S": "/path/to/S.fasta",
                "L": "/path/to/L.fasta",
            },
        )

    def test_single_reference_defaults_to_false(self):
        """--single-reference is off unless explicitly passed."""
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0) as mock_main:
            result = self.runner.invoke(consensus, self._required, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertFalse(mock_main.call_args[0][0]["single_reference"])

    def test_single_reference_flag_threads_into_args(self):
        """--single-reference sets args["single_reference"] to True."""
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0) as mock_main:
            result = self.runner.invoke(
                consensus, self._required + ["--single-reference"], catch_exceptions=False
            )
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertTrue(mock_main.call_args[0][0]["single_reference"])

    def test_default_values_optional_args(self):
        """Check that all optional args have correct defaults for illumina."""
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0) as mock_main:
            result = self.runner.invoke(consensus, self._required, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        args = mock_main.call_args[0][0]
        self.assertEqual(args["data_type"], "illumina")
        self.assertEqual(args["run_name"], "undefined")
        self.assertIsNone(args["adapters"])
        self.assertEqual(args["trim_head"], 0)
        self.assertEqual(args["trim_tail"], 0)
        self.assertEqual(args["cut_front_mean_quality"], 10)
        self.assertEqual(args["cut_tail_mean_quality"], 10)
        self.assertEqual(args["cut_right_window_size"], 4)
        self.assertEqual(args["cut_right_mean_quality"], 15)
        self.assertEqual(args["af_threshold"], 0.51)
        self.assertEqual(args["af_isnv_threshold"], 0.0)
        self.assertFalse(args["run_isnv"])
        self.assertEqual(args["minimum_coverage"], 20)
        self.assertEqual(args["minimum_read_length"], 50)
        self.assertEqual(args["threads"], 1)
        self.assertEqual(args["threads_total"], detect_cores())
        self.assertFalse(args["create_config_only"])


class Test_RemovedOptionsRejected(unittest.TestCase):
    """Options dropped in the split from ViralUnity (HTML report, gene
    annotation) must be rejected by click rather than silently ignored."""

    REMOVED = [
        ["--gene-annotation", "genes.gff3"],
        ["--segmented-gene-annotation", "S=genes.gff3"],
        ["--generate-html-report"],
        ["--no-generate-html-report"],
    ]

    def setUp(self):
        self.runner = CliRunner()

    def _required(self, data_type):
        return [
            data_type,
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ]

    def test_removed_options_are_unknown(self):
        for data_type in ("illumina", "nanopore"):
            for extra in self.REMOVED:
                with self.subTest(data_type=data_type, option=extra[0]):
                    with patch("viralconseq.consensus_cli.consensus_main", return_value=0):
                        result = self.runner.invoke(consensus, self._required(data_type) + extra)
                    self.assertEqual(result.exit_code, 2, result.output)
                    self.assertIn("No such option", result.output)


class Test_ThreadOptionsRequireAtLeastOne(unittest.TestCase):
    """``--threads``, ``--threads-total`` and every ``--<rule>-cpus/-ram`` are
    ``IntRange(min=1)``: a zero would make Snakemake plan with no cores or a
    rule with no memory, which fails late and confusingly."""

    def setUp(self):
        self.runner = CliRunner()

    def _required(self, data_type):
        return [
            data_type,
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ]

    def test_zero_or_negative_values_rejected(self):
        cases = [
            ("illumina", ["--threads", "0"]),
            ("illumina", ["--threads-total", "-1"]),
            ("illumina", ["--map-reads-cpus", "0"]),
            ("nanopore", ["--infer-consensus-sequence-ram", "0"]),
        ]
        for data_type, extra in cases:
            with self.subTest(option=extra[0]):
                with patch("viralconseq.consensus_cli.consensus_main", return_value=0):
                    result = self.runner.invoke(consensus, self._required(data_type) + extra)
                self.assertEqual(result.exit_code, 2, result.output)
                self.assertIn("is not in the range", result.output)


class Test_Clair3ModelOptions(unittest.TestCase):
    def _invoke(self, extra, env=None):
        args = [
            "nanopore",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ] + extra
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0) as mock_main:
            result = CliRunner().invoke(consensus, args, env=env, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        return mock_main.call_args[0][0]

    def test_explicit_model_and_dir_thread_into_args(self):
        args = self._invoke(
            ["--clair3-model", "r941_prom_hac_g360+g422", "--clair3-model-dir", "/m"]
        )
        self.assertEqual(args["clair3_model"], "r941_prom_hac_g360+g422")
        self.assertEqual(args["clair3_model_dir"], "/m")

    def test_model_dir_env_var(self):
        args = self._invoke([], env={"VIRALCONSEQ_CLAIR3_MODELS": "/from/env"})
        self.assertEqual(args["clair3_model_dir"], "/from/env")

    def test_illumina_has_no_clair3_options(self):
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0):
            result = CliRunner().invoke(
                consensus,
                [
                    "illumina",
                    "--sample-sheet",
                    "s.csv",
                    "--config-file",
                    "c.yml",
                    "--output",
                    "o",
                    "--reference",
                    "r.fa",
                    "--clair3-model-dir",
                    "/m",
                ],
            )
        self.assertEqual(result.exit_code, 2)


class Test_ResourceOptions(unittest.TestCase):
    """Per-rule ``--<rule>-cpus`` / ``--<rule>-ram`` options are routed through
    ``**kwargs`` into the args dict under their snake_case key."""

    def setUp(self):
        self.runner = CliRunner()

    def _invoke(self, data_type, extra):
        args = [
            data_type,
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ] + extra
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0) as mock_main:
            result = self.runner.invoke(consensus, args, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        return mock_main.call_args[0][0]

    def test_illumina_resource_options_land_in_args(self):
        args = self._invoke("illumina", ["--map-reads-cpus", "8", "--run-viralqc-ram", "3"])
        self.assertEqual(args["map_reads_cpus"], 8)
        self.assertEqual(args["run_viralqc_ram"], 3)
        # Untouched rules carry None: the config then omits <rule>_cpus and the
        # rule falls back to --threads.
        self.assertIsNone(args["detect_isnv_cpus"])
        # Only memory-declaring rules have a --ram option.
        self.assertNotIn("trim_primer_sequences_ram", args)
        self.assertNotIn("perform_qc_ram", args)

    def test_ram_option_exists_only_for_memory_rules(self):
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0):
            result = self.runner.invoke(
                consensus, self._invoke_args("illumina", ["--perform-qc-ram", "16"])
            )
        self.assertEqual(result.exit_code, 2, result.output)
        self.assertIn("No such option", result.output)

    def _invoke_args(self, data_type, extra):
        return [
            data_type,
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ] + extra

    def test_nanopore_resource_options_land_in_args(self):
        args = self._invoke("nanopore", ["--infer-consensus-sequence-cpus", "6"])
        self.assertEqual(args["infer_consensus_sequence_cpus"], 6)
        self.assertNotIn("perform_qc_cpus", args)

    def test_max_memory_and_threads_total(self):
        args = self._invoke("nanopore", ["--max-memory", "12"])
        self.assertEqual(args["max_memory"], 12)
        # No explicit --max-memory -> None (detected later, during validation).
        self.assertIsNone(self._invoke("nanopore", [])["max_memory"])
        with patch("viralconseq.consensus_cli.detect_cores", return_value=7):
            self.assertEqual(self._invoke("illumina", [])["threads_total"], 7)
        self.assertEqual(self._invoke("illumina", ["--threads-total", "3"])["threads_total"], 3)


class Test_ConsensusNanoporeCommand(unittest.TestCase):
    """Tests for `viralconseq consensus nanopore`."""

    def setUp(self):
        self.runner = CliRunner()
        self._required = [
            "nanopore",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ]

    def test_required_args_success(self):
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0):
            result = self.runner.invoke(consensus, self._required, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)

    def test_default_values_optional_args(self):
        """Check nanopore-specific defaults."""
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0) as mock_main:
            result = self.runner.invoke(consensus, self._required, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        args = mock_main.call_args[0][0]
        self.assertEqual(args["data_type"], "nanopore")
        self.assertEqual(args["af_threshold"], 0.51)
        self.assertEqual(args["chunk_size"], 10000)
        self.assertEqual(args["clair3_model"], "auto")
        self.assertTrue(
            args["clair3_model_dir"].endswith(
                os.path.join(".cache", "viralconseq", "clair3-models")
            )
        )
        self.assertEqual(args["variant_quality"], 20)
        self.assertEqual(args["variant_depth"], 10)
        self.assertEqual(args["minimum_map_quality"], 30)


class Test_ConsensusCondaPrefix(unittest.TestCase):
    """``--conda-prefix`` (and the ``$VIRALCONSEQ_CONDA_PREFIX`` env var) must
    land in the args dict passed to ``consensus_main`` so the orchestrator
    can forward it to Snakemake."""

    def setUp(self):
        self.runner = CliRunner()
        self._required_illumina = [
            "illumina",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ]
        self._required_nanopore = [
            "nanopore",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ]

    def _invoke(self, cli_args, env=None):
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0) as mock_main:
            result = self.runner.invoke(consensus, cli_args, env=env or {}, catch_exceptions=False)
        return result, mock_main

    def test_explicit_conda_prefix_illumina(self):
        result, mock_main = self._invoke(self._required_illumina + ["--conda-prefix", "/tmp/foo"])
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(mock_main.call_args[0][0]["conda_prefix"], "/tmp/foo")

    def test_explicit_conda_prefix_nanopore(self):
        result, mock_main = self._invoke(self._required_nanopore + ["--conda-prefix", "/tmp/foo"])
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(mock_main.call_args[0][0]["conda_prefix"], "/tmp/foo")

    def test_default_conda_prefix(self):
        # Ensure env var is unset so we hit the fallback.
        env = {"VIRALCONSEQ_CONDA_PREFIX": ""}
        env.pop("VIRALCONSEQ_CONDA_PREFIX")
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop("VIRALCONSEQ_CONDA_PREFIX", None)
            result, mock_main = self._invoke(self._required_illumina)
        self.assertEqual(result.exit_code, 0, result.output)
        expected = str(Path.home() / ".cache" / "viralconseq" / "conda-envs")
        self.assertEqual(mock_main.call_args[0][0]["conda_prefix"], expected)

    def test_env_var_conda_prefix(self):
        with patch.dict(os.environ, {"VIRALCONSEQ_CONDA_PREFIX": "/srv/shared/envs"}):
            result, mock_main = self._invoke(self._required_illumina)
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(mock_main.call_args[0][0]["conda_prefix"], "/srv/shared/envs")


class Test_ConsensusViralQCOptions(unittest.TestCase):
    """``--run-viralqc/--no-run-viralqc`` and ``--viralqc-db`` (plus the
    ``$VIRALCONSEQ_VIRALQC_DB`` env var) must land in the args dict."""

    def setUp(self):
        self.runner = CliRunner()

    def _required(self, data_type):
        return [
            data_type,
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--reference",
            "reference.fasta",
        ]

    def _invoke(self, cli_args, env=None):
        with patch("viralconseq.consensus_cli.consensus_main", return_value=0) as mock_main:
            result = self.runner.invoke(consensus, cli_args, env=env or {}, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        return mock_main.call_args[0][0]

    def test_run_viralqc_defaults_true(self):
        for dt in ("illumina", "nanopore"):
            with self.subTest(data_type=dt):
                args = self._invoke(self._required(dt))
                self.assertTrue(args["run_viralqc"])

    def test_no_run_viralqc_flag_threads_into_args(self):
        for dt in ("illumina", "nanopore"):
            with self.subTest(data_type=dt):
                args = self._invoke(self._required(dt) + ["--no-run-viralqc"])
                self.assertFalse(args["run_viralqc"])

    def test_explicit_viralqc_db_threads_into_args(self):
        args = self._invoke(self._required("illumina") + ["--viralqc-db", "/srv/vqc"])
        self.assertEqual(args["viralqc_db"], "/srv/vqc")

    def test_default_viralqc_db_is_cache_dir(self):
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop("VIRALCONSEQ_VIRALQC_DB", None)
            args = self._invoke(self._required("nanopore"))
        self.assertEqual(
            args["viralqc_db"], str(Path.home() / ".cache" / "viralconseq" / "viralqc-db")
        )

    def test_env_var_viralqc_db(self):
        args = self._invoke(
            self._required("illumina"), env={"VIRALCONSEQ_VIRALQC_DB": "/shared/vqc"}
        )
        self.assertEqual(args["viralqc_db"], "/shared/vqc")

    def test_viralqc_db_tilde_expanded(self):
        args = self._invoke(self._required("illumina") + ["--viralqc-db", "~/vqc"])
        self.assertEqual(args["viralqc_db"], os.path.join(str(Path.home()), "vqc"))

    def test_run_viralqc_resource_options_do_not_clash_with_flag(self):
        args = self._invoke(
            self._required("illumina")
            + ["--no-run-viralqc", "--run-viralqc-cpus", "8", "--run-viralqc-ram", "16"]
        )
        self.assertFalse(args["run_viralqc"])
        self.assertEqual(args["run_viralqc_cpus"], 8)
        self.assertEqual(args["run_viralqc_ram"], 16)
        nano = self._invoke(self._required("nanopore"))
        self.assertIsNone(nano["run_viralqc_cpus"])


if __name__ == "__main__":
    unittest.main()
