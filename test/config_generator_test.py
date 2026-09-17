"""Tests for ConfigGenerator.add_samples serialization.

Focus: R1/R2 are stored as a YAML list, not a space-joined string, so a file
path containing a space is not silently corrupted by the workflow's split on
the sample value.
"""

import os
import tempfile
import unittest
from unittest.mock import patch

import yaml

from viralconseq.config_generator import ConfigGenerator
from viralconseq.constants import ConfigKeys, DataType
from viralconseq.exceptions import ConfigurationError


class Test_AddSamples(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def _gen(self):
        return ConfigGenerator(os.path.join(self.tmp, "config.yml"))

    def test_illumina_samples_stored_as_list(self):
        gen = self._gen()
        gen.add_samples({"s1": ["/x/R1.fastq.gz", "/x/R2.fastq.gz"]}, DataType.ILLUMINA)
        self.assertEqual(
            gen.config[ConfigKeys.SAMPLES]["sample-s1"],
            ["/x/R1.fastq.gz", "/x/R2.fastq.gz"],
        )

    def test_nanopore_samples_stored_as_list(self):
        gen = self._gen()
        gen.add_samples({"s1": ["/x/reads.fastq.gz"]}, DataType.NANOPORE)
        self.assertEqual(gen.config[ConfigKeys.SAMPLES]["sample-s1"], ["/x/reads.fastq.gz"])

    def test_space_in_path_survives_yaml_roundtrip(self):
        """The whole point of the list change: a space in a path is preserved."""
        gen = self._gen()
        gen.add_samples(
            {"s1": ["/data dir/R1.fastq.gz", "/data dir/R2.fastq.gz"]}, DataType.ILLUMINA
        )
        gen.save()
        with open(gen.config_path) as fh:
            loaded = yaml.safe_load(fh)
        self.assertEqual(
            loaded["samples"]["sample-s1"],
            ["/data dir/R1.fastq.gz", "/data dir/R2.fastq.gz"],
        )

    def test_illumina_wrong_file_count_raises(self):
        gen = self._gen()
        with self.assertRaises(ConfigurationError):
            gen.add_samples({"s1": ["only_one.fastq.gz"]}, DataType.ILLUMINA)

    def test_nanopore_wrong_file_count_raises(self):
        gen = self._gen()
        with self.assertRaises(ConfigurationError):
            gen.add_samples({"s1": ["a.fastq.gz", "b.fastq.gz"]}, DataType.NANOPORE)


class Test_Save(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def _full_generator(self):
        gen = ConfigGenerator(os.path.join(self.tmp, "nested", "config.yml"))
        gen.add_samples({"s1": ["/x/reads.fastq.gz"]}, DataType.NANOPORE)
        gen.add_output("/out", "run")
        gen.add_threads(2)
        gen.add_consensus_settings("/ref.fa", "NA", 20)
        gen.add_consensus_nanopore_settings(50, 0.51, 10000, "model", 20, 10, 30)
        gen.add_viralqc_settings(run_viralqc=True, viralqc_db="/db")
        gen.add_resource_settings({}, ["map_reads"])
        return gen

    def test_sections_written_in_order_with_comments(self):
        gen = self._full_generator()
        gen.save()
        text = open(gen.config_path).read()
        headers = [line for line in text.splitlines() if line.startswith("# --- ")]
        self.assertEqual(
            headers,
            [
                "# --- run ---",
                "# --- consensus ---",
                "# --- read_qc ---",
                "# --- clair3 ---",
                "# --- viralqc ---",
                "# --- resources ---",
            ],
        )
        # Every section header is followed by at least one explanatory comment.
        lines = text.splitlines()
        for header in headers:
            nxt = lines[lines.index(header) + 1]
            self.assertTrue(nxt.startswith("# ") and not nxt.startswith("# ---"), nxt)

    def test_saved_yaml_roundtrips_and_leaves_no_tempfile(self):
        gen = self._full_generator()
        gen.save()
        with open(gen.config_path) as fh:
            loaded = yaml.safe_load(fh)
        self.assertEqual(loaded, gen.config)
        leftovers = [f for f in os.listdir(os.path.dirname(gen.config_path)) if f != "config.yml"]
        self.assertEqual(leftovers, [])

    def test_write_failure_raises_configuration_error_and_cleans_up(self):
        gen = ConfigGenerator(os.path.join(self.tmp, "config.yml"))
        gen.add_threads(1)
        with patch("viralconseq.config_generator.os.replace", side_effect=OSError("disk full")):
            with self.assertRaises(ConfigurationError):
                gen.save()
        self.assertEqual(os.listdir(self.tmp), [])


if __name__ == "__main__":
    unittest.main()
