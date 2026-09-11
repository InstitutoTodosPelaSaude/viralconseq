"""Tests for viralconseq create-samplesheet CLI command (click-based)."""

import csv
import os
import tempfile
import unittest
from unittest.mock import patch

from click.testing import CliRunner

from viralconseq.create_samplesheet import (
    create_samplesheet,
    extract_sample_name,
    find_samples_level_0,
    find_samples_level_1,
    generate_sample_sheet,
    validate_args,
    write_sample_sheet,
)
from viralconseq.exceptions import ValidationError


class Test_SampleSheetWritingSafety(unittest.TestCase):
    def test_extract_sample_name_rejects_unsafe_name(self):
        # A name with a space would corrupt the sheet / become a wildcard hazard;
        # reject it at the source instead of downstream.
        with self.assertRaises(ValidationError):
            extract_sample_name("bad name_R1.fastq", "_")

    def test_write_sample_sheet_quotes_paths_with_commas(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "s.csv")
            write_sample_sheet({"S1": ["a,b.fastq", "c.fastq"]}, out)
            with open(out) as f:
                rows = list(csv.reader(f))
            # A comma in a path must stay inside one field, not split the row.
            self.assertEqual(rows[0], ["S1", "a,b.fastq", "c.fastq"])


class Test_CreateSamplesheeetCommand(unittest.TestCase):
    """Tests for `viralconseq create-samplesheet`."""

    def setUp(self):
        self.runner = CliRunner()

    def _invoke(self, extra_args=None):
        args = ["--input", "input/dir", "--output", "output.file"] + (extra_args or [])
        with (
            patch("viralconseq.create_samplesheet.validate_args"),
            patch("viralconseq.create_samplesheet.generate_sample_sheet"),
        ):
            return self.runner.invoke(create_samplesheet, args, catch_exceptions=False)

    def test_get_args_required(self):
        """Missing required args (--input, --output) should exit non-zero."""
        result = self.runner.invoke(create_samplesheet, ["--level", "1"])
        self.assertNotEqual(result.exit_code, 0)

    def test_required_args_success_when_only_required_set(self):
        result = self._invoke()
        self.assertEqual(result.exit_code, 0, result.output)

    def test_default_values_optional_args(self):
        """Check that defaults match the expected values."""
        with (
            patch("viralconseq.create_samplesheet.validate_args") as mock_validate,
            patch("viralconseq.create_samplesheet.generate_sample_sheet"),
        ):
            result = self.runner.invoke(
                create_samplesheet,
                ["--input", "input/dir", "--output", "output.file"],
                catch_exceptions=False,
            )
        self.assertEqual(result.exit_code, 0, result.output)
        called_args = mock_validate.call_args[0][0]
        self.assertEqual(called_args["input"], "input/dir")
        self.assertEqual(called_args["output"], "output.file")
        self.assertEqual(called_args["level"], 1)
        self.assertEqual(called_args["pattern"], "R1")
        self.assertEqual(called_args["separator"], "-")


class Test_ValidateArgs(unittest.TestCase):
    def setUp(self):
        self.args = {
            "input": "input/dir",
            "output": "output.file",
            "level": 1,
        }

    @patch("os.path.isdir", return_value=True)
    @patch("os.path.isfile", return_value=False)
    def test_validate_args_success(self, mock_isfile, mock_isdir):
        validated_args = validate_args(self.args)
        self.assertEqual(validated_args, None)

    @patch("os.path.isdir", return_value=False)
    def test_validate_args_input_not_exist(self, mock_isdir):
        with self.assertRaises(Exception):
            validate_args(self.args)

    @patch("os.path.isdir", return_value=True)
    @patch("os.path.isfile", return_value=True)
    def test_validate_args_output_exists(self, mock_isfile, mock_isdir):
        with self.assertRaises(Exception):
            validate_args(self.args)


class Test_GenerateSamplesheet(unittest.TestCase):
    def setUp(self):
        self.args = {
            "input": "input/dir",
            "output": "output.file",
            "separator": "_",
            "pattern": "R1",
        }

    @patch(
        "glob.glob",
        side_effect=[
            ["input/dir/1", "input/dir/2"],
            ["input/dir/1/R1_sample1.fastq", "input/dir/1/R1_sample2.fastq"],
            ["input/dir/2/R2_sample1.fastq", "input/dir/2/R2_sample2.fastq"],
        ],
    )
    @patch("os.path.isdir", return_value=True)
    @patch("os.path.isfile", return_value=True)
    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    def test_generate_samplesheet_level_1(self, mock_open, mock_isfile, mock_isdir, mock_glob):
        self.args["level"] = 1
        generate_sample_sheet(self.args)
        mock_open.assert_called_with("output.file", "w", newline="")
        handle = mock_open()
        handle.write.assert_any_call(
            "1,input/dir/1/R1_sample1.fastq,input/dir/1/R1_sample2.fastq\n"
        )
        handle.write.assert_any_call(
            "2,input/dir/2/R2_sample1.fastq,input/dir/2/R2_sample2.fastq\n"
        )

    @patch(
        "glob.glob",
        side_effect=[
            ["input/dir/1/R1_sample1.fastq", "input/dir/1/R1_sample2.fastq"],
            ["input/dir/1/R1_sample1.fastq", "input/dir/1/R1_sample2.fastq"],
            ["input/dir/1/R1_sample1.fastq", "input/dir/1/R1_sample2.fastq"],
        ],
    )
    @patch("os.path.isfile", return_value=True)
    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    def test_generate_samplesheet_level_0(self, mock_open, mock_isfile, mock_glob):
        self.args["level"] = 0
        generate_sample_sheet(self.args)
        mock_open.assert_called_with("output.file", "w", newline="")
        handle = mock_open()
        handle.write.assert_any_call(
            "R1,input/dir/1/R1_sample1.fastq,input/dir/1/R1_sample2.fastq\n"
        )


class Test_SampleGroupingIntegrity(unittest.TestCase):
    """Regressions for silent sample collisions in create-samplesheet."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def _touch(self, *relparts):
        path = os.path.join(self.tmp, *relparts)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        open(path, "w").close()
        return path

    def test_level0_prefix_no_substring_collision(self):
        """'s1' must not swallow 's10' files (old substring glob bug)."""
        self._touch("s1_R1.fastq.gz")
        self._touch("s1_R2.fastq.gz")
        self._touch("s10_R1.fastq.gz")
        self._touch("s10_R2.fastq.gz")

        samples = find_samples_level_0(self.tmp, separator="_", pattern="R1")

        self.assertEqual(set(samples), {"s1", "s10"})
        self.assertEqual(len(samples["s1"]), 2)
        self.assertEqual(len(samples["s10"]), 2)
        self.assertTrue(all(os.path.basename(p).startswith("s1_") for p in samples["s1"]))
        self.assertTrue(all(os.path.basename(p).startswith("s10_") for p in samples["s10"]))

    def test_level1_duplicate_sample_names_rejected(self):
        """Two subdirectories that reduce to the same sample name must error."""
        self._touch("sampleA-1", "reads.fastq.gz")
        self._touch("sampleA-2", "reads.fastq.gz")

        with self.assertRaises(ValidationError):
            find_samples_level_1(self.tmp, separator="-")


if __name__ == "__main__":
    unittest.main()
