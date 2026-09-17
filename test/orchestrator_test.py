"""Tests for viralconseq._orchestrator.run_workflow.

Focused on the Snakemake-call kwargs the orchestrator threads through —
notably ``conda_prefix``, which Snakemake uses to share per-rule conda
envs across working directories.
"""

import os
import tempfile
import unittest
from unittest.mock import patch

from viralconseq import _orchestrator


def _write_dummy_workflow(path: str) -> None:
    with open(path, "w") as fh:
        fh.write("rule all:\n    input: []\n")


class Test_RunWorkflowForwardsCondaPrefix(unittest.TestCase):
    """``run_workflow`` must forward ``args['conda_prefix']`` to ``snakemake()``."""

    def _run_with_args(self, args: dict) -> dict:
        """Invoke run_workflow against a real on-disk dummy .smk and return the
        kwargs the patched ``snakemake`` received."""
        captured = {}

        def fake_snakemake(*_args, **kwargs):
            captured.update(kwargs)
            return True

        with tempfile.TemporaryDirectory() as tmp:
            workflow_path = os.path.join(tmp, "dummy.smk")
            _write_dummy_workflow(workflow_path)
            config_path = os.path.join(tmp, "config.yaml")
            with open(config_path, "w") as fh:
                fh.write("samples: {}\n")
            args = {**args, "config_file": config_path, "threads_total": 1}
            with patch("viralconseq._orchestrator.snakemake", side_effect=fake_snakemake):
                _orchestrator.run_workflow(workflow_path, args)
        return captured

    def test_conda_prefix_forwarded_when_set(self):
        captured = self._run_with_args({"conda_prefix": "/tmp/viralconseq-envs"})
        self.assertEqual(captured.get("conda_prefix"), "/tmp/viralconseq-envs")
        self.assertTrue(captured.get("use_conda"))

    def test_conda_prefix_none_when_missing(self):
        """Absent ``conda_prefix`` -> ``None`` (preserves pre-fix per-workdir behaviour)."""
        captured = self._run_with_args({})
        self.assertIsNone(captured.get("conda_prefix"))
        self.assertTrue(captured.get("use_conda"))

    def test_relative_conda_prefix_is_made_absolute(self):
        """Snakemake resolves a relative prefix against workdir, so pass it absolute."""
        captured = self._run_with_args({"conda_prefix": "envs-cache"})
        self.assertEqual(captured.get("conda_prefix"), os.path.abspath("envs-cache"))


class Test_RunWorkflowRunsInsideTheRunDirectory(unittest.TestCase):
    """``.snakemake/`` must live in ``<output>/<run_name>``, and interrupted
    runs must resume rather than stop on IncompleteFilesException."""

    def _run_with_args(self, args: dict) -> dict:
        captured = {}

        def fake_snakemake(*_args, **kwargs):
            captured.update(kwargs)
            return True

        with tempfile.TemporaryDirectory() as tmp:
            workflow_path = os.path.join(tmp, "dummy.smk")
            _write_dummy_workflow(workflow_path)
            config_path = os.path.join(tmp, "config.yaml")
            with open(config_path, "w") as fh:
                fh.write("samples: {}\n")
            args = {**args, "config_file": config_path, "threads_total": 1}
            with patch("viralconseq._orchestrator.snakemake", side_effect=fake_snakemake):
                _orchestrator.run_workflow(workflow_path, args)
        return captured

    def test_workdir_and_force_incomplete_forwarded(self):
        captured = self._run_with_args({"output": "/data/out", "run_name": "runA"})
        self.assertEqual(captured.get("workdir"), os.path.join("/data/out", "runA"))
        self.assertTrue(captured.get("force_incomplete"))

    def test_workdir_none_without_output(self):
        captured = self._run_with_args({})
        self.assertIsNone(captured.get("workdir"))
        self.assertTrue(captured.get("force_incomplete"))

    def test_mode_kwargs_forwarded(self):
        captured = self._run_with_args({})
        self.assertEqual(
            (captured["dryrun"], captured["unlock"], captured["keepgoing"]), (False, False, False)
        )

    def test_memory_budget_forwarded_as_resources(self):
        captured = self._run_with_args({"max_memory_mb": 4096})
        self.assertEqual(captured.get("resources"), {"mem_mb": 4096})

    def test_no_resources_kwarg_when_budget_is_off(self):
        captured = self._run_with_args({"max_memory_mb": 0})
        self.assertNotIn("resources", captured)
        self.assertNotIn("resources", self._run_with_args({}))

    def test_describe_resources(self):
        rules = ["map_reads", "infer_consensus_sequence", "run_viralqc"]
        off = _orchestrator.describe_resources(
            {"threads_total": 4, "max_memory_mb": 0, "memory_budget_source": "--max-memory"}, rules
        )
        self.assertIn("4 core(s)", off)
        self.assertIn("memory budget off", off)
        on = _orchestrator.describe_resources(
            {
                "threads_total": 8,
                "max_memory_mb": 8192,
                "memory_detected_mb": 16384,
                "memory_budget_source": "detected 16.0 GB, 10% headroom",
            },
            rules,
        )
        self.assertIn("memory budget 8.0 GB", on)
        self.assertIn("infer_consensus_sequence 2 GB -> at most 4 at once", on)
        self.assertIn("run_viralqc 1 GB -> at most 8 at once", on)
        self.assertEqual(on.count("detected"), 1)  # not repeated after the source text
        explicit = _orchestrator.describe_resources(
            {
                "threads_total": 8,
                "max_memory_mb": 4096,
                "memory_detected_mb": 16384,
                "memory_budget_source": "--max-memory",
            },
            rules,
        )
        self.assertIn("(--max-memory, detected 16.0 GB)", explicit)

    def test_run_dir_for(self):
        self.assertIsNone(_orchestrator.run_dir_for({"output": "x"}))
        self.assertIsNone(_orchestrator.run_dir_for({"run_name": "x"}))
        self.assertEqual(
            _orchestrator.run_dir_for({"output": "out", "run_name": "r"}),
            os.path.abspath(os.path.join("out", "r")),
        )


if __name__ == "__main__":
    unittest.main()
