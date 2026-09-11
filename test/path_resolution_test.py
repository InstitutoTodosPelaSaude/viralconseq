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

    def test_returns_args_dict_for_chaining(self):
        args = {"output": "results"}
        returned = resolve_path_args(args, ("output",), base_dir="/x")
        self.assertIs(returned, args)


# ---------------------------------------------------------------------------
# run_snakemake_workflow must not pin workdir to the config directory
# ---------------------------------------------------------------------------


class TestWorkdirUsesCwd(unittest.TestCase):
    """``run_snakemake_workflow`` must not pin the working directory to the
    config file's parent directory (this was the bug).

    The fix relies on Snakemake's default behaviour: when ``workdir`` is not
    passed to ``snakemake()``, Snakemake does not ``chdir`` and the shell's
    cwd is used as-is. This test asserts that ``workdir`` is either omitted
    or set to cwd — never the config file's directory.
    """

    def test_workdir_is_not_config_dirname(self):
        captured = {}

        def fake_snakemake(*args, **kwargs):
            captured.update(kwargs)
            return True

        with (
            patch("viralconseq._orchestrator.snakemake", side_effect=fake_snakemake),
            tempfile.TemporaryDirectory() as tmp,
        ):
            tmp_real = os.path.realpath(tmp)
            old_cwd = os.getcwd()
            os.chdir(tmp_real)
            try:
                consensus.run_snakemake_workflow(
                    {
                        "data_type": "nanopore",
                        "reference": os.path.join(tmp_real, "ref.fasta"),
                        # Pretend the user put the config in a sub-directory
                        # so the old code would have set workdir=that subdir.
                        "config_file": os.path.join(tmp_real, "scratch", "run.yml"),
                        "threads_total": 1,
                    }
                )
            finally:
                os.chdir(old_cwd)

        config_dir = os.path.join(tmp_real, "scratch")
        passed_workdir = captured.get("workdir")
        self.assertNotEqual(
            passed_workdir,
            config_dir,
            "workdir must NOT be the config file's parent directory",
        )
        # Either workdir was omitted (Snakemake defaults to cwd) or it was
        # explicitly set to cwd. Both are acceptable.
        if passed_workdir is not None:
            self.assertEqual(passed_workdir, tmp_real)


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
