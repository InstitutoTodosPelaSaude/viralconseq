"""Tests for ``viralconseq create-report`` (report_cli.py)."""

import os
import tempfile
import unittest
from unittest.mock import patch

from click.testing import CliRunner

from viralconseq.cli import cli
from viralconseq.exceptions import ReportError
from viralconseq.report_cli import create_report_command


class Test_CreateReportCommand(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()
        self._tmp = tempfile.TemporaryDirectory()
        self.run_dir = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def test_registered_on_the_top_level_cli(self):
        self.assertIn("create-report", cli.commands)
        result = self.runner.invoke(cli, ["create-report", "--help"])
        self.assertEqual(result.exit_code, 0)
        self.assertIn("RUN_DIR", result.output)

    def test_passes_arguments_through_and_prints_the_path(self):
        with patch("viralconseq.report_cli.create_report", return_value="/x/report.html") as create:
            result = self.runner.invoke(
                create_report_command,
                [self.run_dir, "-o", "/x/report.html", "--label", "Run 1"],
            )
        self.assertEqual(result.exit_code, 0, result.output)
        create.assert_called_once_with(self.run_dir, output="/x/report.html", label="Run 1")
        self.assertIn("/x/report.html", result.output)

    def test_defaults_are_none(self):
        with patch("viralconseq.report_cli.create_report", return_value="p") as create:
            result = self.runner.invoke(create_report_command, [self.run_dir])
        self.assertEqual(result.exit_code, 0, result.output)
        create.assert_called_once_with(self.run_dir, output=None, label=None)

    def test_report_error_becomes_a_clean_exit(self):
        with patch(
            "viralconseq.report_cli.create_report", side_effect=ReportError("no summary.tsv")
        ):
            result = self.runner.invoke(create_report_command, [self.run_dir])
        self.assertEqual(result.exit_code, 1)
        self.assertIn("[report_error] no summary.tsv", result.output)

    def test_missing_directory_is_a_usage_error(self):
        result = self.runner.invoke(create_report_command, [os.path.join(self.run_dir, "nope")])
        self.assertEqual(result.exit_code, 2)


if __name__ == "__main__":
    unittest.main()
