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
        gen.add_report_settings(run_report=False)
        gen.add_resource_settings({"map_reads_cpus": 4}, ["map_reads", "run_viralqc"])
        return gen

    def test_resource_settings_contract(self):
        gen = ConfigGenerator(os.path.join(self.tmp, "c.yml"))
        gen.add_resource_settings(
            {"map_reads_cpus": 4, "run_viralqc_ram": 3}, ["map_reads", "perform_qc", "run_viralqc"]
        )
        # cpus only when given; ram only for memory rules (default when not given)
        self.assertEqual(gen.config, {"map_reads_cpus": 4, "run_viralqc_ram": 3})
        gen2 = ConfigGenerator(os.path.join(self.tmp, "d.yml"))
        gen2.add_resource_settings({}, ["infer_consensus_sequence", "map_reads"])
        self.assertEqual(gen2.config, {"infer_consensus_sequence_ram": 2})

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
                "# --- report ---",
                "# --- resources ---",
            ],
        )
        # Every section header is followed by at least one explanatory comment.
        lines = text.splitlines()
        for header in headers:
            nxt = lines[lines.index(header) + 1]
            self.assertTrue(nxt.startswith("# ") and not nxt.startswith("# ---"), nxt)

    def test_add_provenance_lands_in_its_own_section(self):
        gen = self._full_generator()
        gen.add_provenance("1.2.3")
        gen.save()
        text = open(gen.config_path).read()
        self.assertIn("# --- provenance ---", text)
        self.assertEqual(yaml.safe_load(text)["viralconseq_version"], "1.2.3")

    def test_saved_yaml_roundtrips_and_leaves_no_tempfile(self):
        gen = self._full_generator()
        gen.save()
        with open(gen.config_path) as fh:
            loaded = yaml.safe_load(fh)
        self.assertEqual(loaded, gen.config)
        leftovers = [f for f in os.listdir(os.path.dirname(gen.config_path)) if f != "config.yml"]
        self.assertEqual(leftovers, [])

    def test_per_sample_clair3_models_are_rekeyed_like_samples(self):
        gen = ConfigGenerator(os.path.join(self.tmp, "c.yml"))
        gen.add_consensus_nanopore_settings(
            50, 0.51, 10000, {"a": "m1", "b": "m2"}, 20, 10, 30, clair3_model_dir="/models"
        )
        self.assertEqual(gen.config["clair3_model"], {"sample-a": "m1", "sample-b": "m2"})
        self.assertEqual(gen.config["clair3_model_dir"], "/models")

    def test_from_dict_places_keys_in_sections_and_backup_keeps_previous(self):
        path = os.path.join(self.tmp, "config.yml")
        with open(path, "w") as fh:
            fh.write("data: nanopore\n")
        gen = ConfigGenerator.from_dict(
            path,
            {
                "samples": {},
                "data": "nanopore",
                "clair3_model": "m",
                "run_viralqc": True,
                "map_reads_cpus": 4,
                "max_memory_mb": 0,
                "viralconseq_version": "1",
                "odd_key": 1,
            },
        )
        self.assertEqual(gen.section_for("clair3_model"), "clair3")
        self.assertEqual(gen.section_for("map_reads_cpus"), "resources")
        self.assertEqual(gen.section_for("max_memory_mb"), "resources")
        self.assertEqual(gen.section_for("viralconseq_version"), "provenance")
        self.assertEqual(gen.section_for("odd_key"), "run")
        gen.save(backup=True)
        with open(path + ".bak") as fh:
            self.assertEqual(fh.read(), "data: nanopore\n")
        with open(path) as fh:
            self.assertEqual(yaml.safe_load(fh)["clair3_model"], "m")

    def test_write_failure_raises_configuration_error_and_cleans_up(self):
        gen = ConfigGenerator(os.path.join(self.tmp, "config.yml"))
        gen.add_threads(1)
        with patch("viralconseq.config_generator.os.replace", side_effect=OSError("disk full")):
            with self.assertRaises(ConfigurationError):
                gen.save()
        self.assertEqual(os.listdir(self.tmp), [])


if __name__ == "__main__":
    unittest.main()
