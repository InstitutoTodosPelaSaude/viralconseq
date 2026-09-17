"""Tests for path-argument resolution in the viralconseq CLI.

An earlier CLI passed ``workdir=os.path.dirname(<config_file>)`` to Snakemake,
which silently re-rooted every relative path the user provided on the command
line. Two behaviours prevent that:

1. ``consensus.main`` resolves every path-typed CLI argument to an absolute
   path against ``os.getcwd()`` *before* validation runs.
2. ``consensus.run_snakemake_workflow`` does not pass ``workdir`` to
   Snakemake, so Snakemake uses the shell's cwd (which is what we want)
   rather than the config file's parent directory.

These tests pin down both behaviours.
"""

import os
import tempfile
import unittest
from unittest.mock import patch

from viralconseq import consensus
from viralconseq.validators import (
    CONSENSUS_PATH_ARG_KEYS,
    _is_path_sentinel,
    resolve_path_args,
)

# ---------------------------------------------------------------------------
# _is_path_sentinel
# ---------------------------------------------------------------------------


class TestIsPathSentinel(unittest.TestCase):
    def test_none_is_sentinel(self):
        self.assertTrue(_is_path_sentinel(None))

    def test_empty_string_is_sentinel(self):
        self.assertTrue(_is_path_sentinel(""))
        self.assertTrue(_is_path_sentinel("   "))

    def test_NA_is_sentinel(self):
        self.assertTrue(_is_path_sentinel("NA"))
        self.assertTrue(_is_path_sentinel("  NA  "))

    def test_lowercase_na_is_not_sentinel(self):
        """Matches the case-sensitive ``"NA"`` check used by the validators."""
        self.assertFalse(_is_path_sentinel("na"))

    def test_non_string_is_sentinel(self):
        self.assertTrue(_is_path_sentinel(0))
        self.assertTrue(_is_path_sentinel(False))
        self.assertTrue(_is_path_sentinel(123))

    def test_real_paths_are_not_sentinel(self):
        self.assertFalse(_is_path_sentinel("/abs/path"))
        self.assertFalse(_is_path_sentinel("relative/path"))
        self.assertFalse(_is_path_sentinel("file.fasta"))


# ---------------------------------------------------------------------------
# resolve_path_args
# ---------------------------------------------------------------------------


class TestResolvePathArgs(unittest.TestCase):
    def test_relative_paths_become_absolute_against_base_dir(self):
        args = {
            "reference": "refs/ref.fasta",
            "primer_scheme": "refs/scheme.bed",
        }
        resolve_path_args(args, CONSENSUS_PATH_ARG_KEYS, base_dir="/work/proj")
        self.assertEqual(args["reference"], "/work/proj/refs/ref.fasta")
        self.assertEqual(args["primer_scheme"], "/work/proj/refs/scheme.bed")

    def test_absolute_paths_are_left_untouched(self):
        args = {"reference": "/some/abs/path/ref.fasta"}
        resolve_path_args(args, CONSENSUS_PATH_ARG_KEYS, base_dir="/work/proj")
        self.assertEqual(args["reference"], "/some/abs/path/ref.fasta")

    def test_sentinels_are_left_untouched(self):
        args = {
            "primer_scheme": "NA",
            "reference": "  NA  ",
            "adapters": None,
            "output": "",
        }
        resolve_path_args(args, CONSENSUS_PATH_ARG_KEYS, base_dir="/work/proj")
        self.assertEqual(args["primer_scheme"], "NA")
        self.assertEqual(args["reference"], "  NA  ")
        self.assertIsNone(args["adapters"])
        self.assertEqual(args["output"], "")

    def test_missing_keys_are_ignored(self):
        args = {"reference": "refs/ref.fasta"}
        resolve_path_args(args, CONSENSUS_PATH_ARG_KEYS, base_dir="/work/proj")
        self.assertEqual(args["reference"], "/work/proj/refs/ref.fasta")
        for k in CONSENSUS_PATH_ARG_KEYS:
            if k == "reference":
                continue
            self.assertNotIn(k, args)

    def test_segmented_reference_dict_values_are_resolved(self):
        args = {
            "reference": {
                "L": "refs/L.fasta",
                "S": "/abs/refs/S.fasta",
            }
        }
        resolve_path_args(args, ("reference",), base_dir="/work/proj")
        self.assertEqual(
            args["reference"],
            {"L": "/work/proj/refs/L.fasta", "S": "/abs/refs/S.fasta"},
        )

    def test_defaults_to_cwd_when_base_dir_omitted(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_real = os.path.realpath(tmp)
            old_cwd = os.getcwd()
            try:
                os.chdir(tmp_real)
                args = {"output": "results/foo"}
                resolve_path_args(args, ("output",))
                self.assertEqual(args["output"], os.path.join(tmp_real, "results/foo"))
            finally:
                os.chdir(old_cwd)

    def test_viralqc_db_is_a_path_key(self):
        self.assertIn("viralqc_db", CONSENSUS_PATH_ARG_KEYS)
        args = {"viralqc_db": "dbs/vqc"}
        resolve_path_args(args, CONSENSUS_PATH_ARG_KEYS, base_dir="/work/proj")
        self.assertEqual(args["viralqc_db"], "/work/proj/dbs/vqc")
        args = {"viralqc_db": "NA"}
        resolve_path_args(args, CONSENSUS_PATH_ARG_KEYS, base_dir="/work/proj")
        self.assertEqual(args["viralqc_db"], "NA")

    def test_returns_args_dict_for_chaining(self):
        args = {"output": "results"}
        returned = resolve_path_args(args, ("output",), base_dir="/x")
        self.assertIs(returned, args)


# ---------------------------------------------------------------------------
# run_snakemake_workflow runs Snakemake inside the run directory
# ---------------------------------------------------------------------------


class Test_WorkdirIsRunDir(unittest.TestCase):
    """``run_snakemake_workflow`` must run Snakemake with the run directory
    (``<output>/<run_name>``) as its working directory, never the config file's
    parent directory (the historical bug) and never the caller's cwd (which
    scattered ``.snakemake/`` wherever the CLI was launched from).

    Every path in the generated config is absolute by then, so only
    ``.snakemake/`` moves.
    """

    def _captured_kwargs(self, args):
        captured = {}

        def fake_snakemake(*_args, **kwargs):
            captured.update(kwargs)
            return True

        with patch("viralconseq._orchestrator.snakemake", side_effect=fake_snakemake):
            consensus.run_snakemake_workflow(args)
        return captured

    def test_workdir_is_output_slash_run_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_real = os.path.realpath(tmp)
            captured = self._captured_kwargs(
                {
                    "data_type": "nanopore",
                    "reference": os.path.join(tmp_real, "ref.fasta"),
                    # The config lives somewhere else entirely: it must not
                    # influence the working directory.
                    "config_file": os.path.join(tmp_real, "scratch", "run.yml"),
                    "output": os.path.join(tmp_real, "results"),
                    "run_name": "run1",
                    "threads_total": 1,
                }
            )
        self.assertEqual(captured["workdir"], os.path.join(tmp_real, "results", "run1"))
        self.assertNotEqual(captured["workdir"], os.path.join(tmp_real, "scratch"))
        self.assertTrue(captured["force_incomplete"])
        self.assertTrue(os.path.isabs(captured["configfiles"][0]))

    def test_workdir_omitted_without_output(self):
        """Callers that pass no output/run_name (tests, embedding) keep Snakemake's default."""
        captured = self._captured_kwargs(
            {
                "data_type": "illumina",
                "reference": "ref.fasta",
                "config_file": "c.yml",
                "threads_total": 1,
            }
        )
        self.assertIsNone(captured["workdir"])


# ---------------------------------------------------------------------------
# Integration: consensus.main calls resolve_path_args
# ---------------------------------------------------------------------------


class TestConsensusMainResolvesPaths(unittest.TestCase):
    def test_main_resolves_paths_before_validate_args(self):
        captured = {}

        def fake_validate(args):
            captured.update(args)
            return {"sample1": ["/abs/sample1.fastq"]}

        with (
            patch("viralconseq.consensus.validate_args", side_effect=fake_validate),
            patch("viralconseq.consensus.generate_config_file"),
            patch("viralconseq.consensus.run_snakemake_workflow", return_value=True),
            tempfile.TemporaryDirectory() as tmp,
        ):
            tmp_real = os.path.realpath(tmp)
            old_cwd = os.getcwd()
            os.chdir(tmp_real)
            try:
                exit_code = consensus.main(
                    {
                        "data_type": "illumina",
                        "config_file": "scratch/cons.yml",
                        "output": "out/cons",
                        "reference": "refs/ref.fasta",
                        "primer_scheme": "schemes/scheme.bed",
                        "adapters": "adapters/illumina.fa",
                        "create_config_only": True,
                        "threads_total": 1,
                    }
                )
            finally:
                os.chdir(old_cwd)

        self.assertEqual(exit_code, 0)
        for key in (
            "config_file",
            "output",
            "reference",
            "primer_scheme",
            "adapters",
        ):
            self.assertTrue(
                os.path.isabs(captured[key]),
                f"{key} should be absolute, got {captured[key]!r}",
            )
            self.assertTrue(captured[key].startswith(tmp_real))


if __name__ == "__main__":
    unittest.main()
