"""Tests for ConfigGenerator.add_samples serialization.

Focus: R1/R2 are stored as a YAML list, not a space-joined string, so a file
path containing a space is not silently corrupted by the workflow's split on
the sample value.
"""

import os
import tempfile
import unittest

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


if __name__ == "__main__":
    unittest.main()
