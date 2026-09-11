"""Tests for the run-provenance manifest."""

import json
import os
import tempfile
import unittest

from viralconseq import __version__
from viralconseq.provenance import (
    MANIFEST_FILENAME,
    build_run_manifest,
    record_run_completion,
    write_run_manifest,
)


class Test_RunManifest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.r1 = os.path.join(self.tmp, "s_R1.fastq.gz")
        self.r2 = os.path.join(self.tmp, "s_R2.fastq.gz")
        with open(self.r1, "wb") as fh:
            fh.write(b"read1")
        with open(self.r2, "wb") as fh:
            fh.write(b"read2")
        self.args = {
            "output": os.path.join(self.tmp, "out"),
            "run_name": "run1",
            "config_file": os.path.join(self.tmp, "config.yml"),
            "data_type": "illumina",
        }
        self.samples = {"s": [self.r1, self.r2]}

    def tearDown(self):
        self._tmp.cleanup()

    def test_build_manifest_records_version_and_checksums(self):
        manifest = build_run_manifest(
            self.args, self.samples, timestamp="2026-01-01T00:00:00+00:00"
        )
        self.assertEqual(manifest["viralconseq_version"], __version__)
        self.assertEqual(manifest["sample_count"], 1)
        recs = manifest["samples"]["s"]
        self.assertEqual(len(recs), 2)
        self.assertEqual(recs[0]["size_bytes"], 5)
        self.assertEqual(len(recs[0]["sha256"]), 64)
        self.assertNotEqual(recs[0]["sha256"], recs[1]["sha256"])  # r1 != r2

    def test_missing_input_flagged_not_crash(self):
        samples = {"s": [os.path.join(self.tmp, "nope.fastq.gz")]}
        manifest = build_run_manifest(self.args, samples)
        self.assertTrue(manifest["samples"]["s"][0]["missing"])

    def test_write_manifest_creates_json_in_run_dir(self):
        path = write_run_manifest(self.args, self.samples)
        self.assertTrue(path.endswith(os.path.join("out", "run1", MANIFEST_FILENAME)))
        with open(path) as fh:
            loaded = json.load(fh)
        self.assertEqual(loaded["run_name"], "run1")
        self.assertEqual(loaded["data_type"], "illumina")

    def test_config_sha256_none_when_config_absent(self):
        manifest = build_run_manifest(self.args, self.samples)
        self.assertIsNone(manifest["config_sha256"])

    def test_config_sha256_recorded_when_config_present(self):
        with open(self.args["config_file"], "w") as fh:
            fh.write("samples: {}\n")
        manifest = build_run_manifest(self.args, self.samples)
        self.assertEqual(len(manifest["config_sha256"]), 64)

    def test_record_run_completion_updates_manifest(self):
        path = write_run_manifest(self.args, self.samples)
        record_run_completion(path, status="success", timestamp="2026-01-02T00:00:00+00:00")
        with open(path) as fh:
            loaded = json.load(fh)
        self.assertEqual(loaded["status"], "success")
        self.assertEqual(loaded["finished_utc"], "2026-01-02T00:00:00+00:00")


if __name__ == "__main__":
    unittest.main()
