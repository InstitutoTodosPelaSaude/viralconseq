"""Tests for scripts/python/tool_versions.py, driven exactly as the workflow
runs it (``main(argv)``), with the probe subprocesses patched."""

import os
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from viralconseq.scripts.python import tool_versions

# Real first lines printed by the tools in the pinned environments.
OUTPUTS = {
    "minimap2": ("2.30-r1287\n", ""),
    "samtools": ("samtools 1.23.1\nUsing htslib 1.24\n", ""),
    "bedtools": ("bedtools v2.31.1\n", ""),
    "GSAlign": ("\nGenAlign v1.0.22\nUsage: GSAlign [-i IndexFile Prefix]\n", ""),
    "gofasta": ("gofasta version 1.2.3\n", ""),
    "fastp": ("", "fastp 1.1.0\n"),
    "multiqc": ("multiqc, version 1.21\n", ""),
    "lofreq": ("version: 2.1.5\ncommit: unknown\n", ""),
    "bcftools": ("bcftools 1.23.1\n", ""),
    "clair3": ("Clair3 v1.2.0\n", ""),
    "python": ("Python 3.11.16\n", ""),
    "nextclade": ("nextclade 3.15.3\n", ""),
    "blastn": ("blastn: 2.16.0+\n Package: blast 2.16.0, build Mar 28 2025\n", ""),
}
EXPECTED = {
    "minimap2": "2.30-r1287",
    "samtools": "1.23.1",
    "bedtools": "2.31.1",
    "GSAlign": "1.0.22",
    "gofasta": "1.2.3",
    "fastp": "1.1.0",
    "multiqc": "1.21",
    "lofreq": "2.1.5",
    "bcftools": "1.23.1",
    "clair3": "1.2.0",
    "python": "3.11.16",
    "nextclade": "3.15.3",
    "blastn": "2.16.0+",
}


def _fake_run(command, **_kwargs):
    tool = next(t for t, probe in tool_versions.PROBES.items() if probe == list(command))
    out, err = OUTPUTS[tool]
    return subprocess.CompletedProcess(command, 0, stdout=out, stderr=err)


class Test_ParseVersion(unittest.TestCase):
    def test_real_outputs_parse(self):
        for tool, (out, err) in OUTPUTS.items():
            with self.subTest(tool=tool):
                self.assertEqual(tool_versions.parse_version(out + "\n" + err), EXPECTED[tool])

    def test_no_version_token(self):
        self.assertIsNone(tool_versions.parse_version("usage: thing [options]"))


class Test_Main(unittest.TestCase):
    def _run(self, argv):
        with (
            patch.object(tool_versions.shutil, "which", return_value="/env/bin/x"),
            patch.object(tool_versions.subprocess, "run", side_effect=_fake_run),
            patch.object(tool_versions.importlib.metadata, "version", return_value="9.9.9"),
        ):
            return tool_versions.main(argv)

    def test_writes_tsv_with_tools_and_dists(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "v.tsv")
            rc = self._run(
                ["--tools", "minimap2", "GSAlign", "--python-dists", "pandas", "--output", out]
            )
            self.assertEqual(rc, 0)
            with open(out) as fh:
                rows = [line.rstrip("\n").split("\t") for line in fh]
        self.assertEqual(
            rows,
            [
                ["tool", "version"],
                ["minimap2", "2.30-r1287"],
                ["GSAlign", "1.0.22"],
                ["pandas", "9.9.9"],
            ],
        )

    def test_missing_binary_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "v.tsv")
            with patch.object(tool_versions.shutil, "which", return_value=None):
                rc = tool_versions.main(["--tools", "minimap2", "--output", out])
            self.assertEqual(rc, 1)
            self.assertFalse(os.path.exists(out))

    def test_unparseable_output_recorded_verbatim(self):
        def usage_only(command, **_kwargs):
            return subprocess.CompletedProcess(command, 2, stdout="", stderr="usage: tool [opts]\n")

        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "v.tsv")
            with (
                patch.object(tool_versions.shutil, "which", return_value="/env/bin/x"),
                patch.object(tool_versions.subprocess, "run", side_effect=usage_only),
            ):
                rc = tool_versions.main(["--tools", "gofasta", "--output", out])
            self.assertEqual(rc, 0)
            with open(out) as fh:
                self.assertIn("gofasta\tunparsed: usage: tool [opts]", fh.read())

    def test_unknown_tool_exits(self):
        with self.assertRaises(SystemExit):
            tool_versions.probe_tool("no-such-tool")

    def test_missing_dist_reported(self):
        self.assertEqual(tool_versions.probe_dist("definitely-not-installed-xyz"), "not installed")


if __name__ == "__main__":
    unittest.main()
