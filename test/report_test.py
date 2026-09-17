"""Tests for viralconseq.report (the create-report engine): run-dir checks,
config loading, the argv handed to build_report.py and its agreement with the
workflow's rule report and the script's own parser."""

import os
import re
import tempfile
import unittest

import yaml

from viralconseq import report
from viralconseq.exceptions import ConfigurationError, ReportError
from viralconseq.scripts.python import build_report, build_summary

RULES_DIR = os.path.join(os.path.dirname(report.__file__), "scripts", "rules")


def _write_summary(run_dir, rows, segmented=False):
    columns = [c for c in build_summary.COLUMNS if segmented or c != "segment"]
    with open(os.path.join(run_dir, "summary.tsv"), "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            values = dict(zip(columns, ["NA"] * len(columns)))
            values.update(row)
            fh.write("\t".join(str(values[c]) for c in columns) + "\n")


def _ok_row(sample, cov=95.0):
    return {
        "sample_id": sample,
        "status": "ok",
        "total_reads": 1000,
        "qc_passed_reads": 900,
        "mapped_reads": 800,
        "pct_mapped": 80.0,
        "mean_depth": 100.0,
        "median_depth": 90.0,
        "coverage_min_depth": cov,
        "min_depth": 20,
        "consensus_length": 100,
        "n_count": 5,
        "n_pct": 5.0,
    }


class _RunDir(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.run_dir = os.path.join(self._tmp.name, "run1")
        os.makedirs(self.run_dir)

    def tearDown(self):
        self._tmp.cleanup()

    def write_config(self, **config):
        with open(os.path.join(self.run_dir, "config.yml"), "w") as fh:
            yaml.safe_dump(config, fh)


class Test_ValidateRunDir(_RunDir):
    def test_missing_summary_is_a_report_error(self):
        with self.assertRaises(ReportError) as ctx:
            report.validate_run_dir(self.run_dir)
        self.assertIn("summary.tsv", str(ctx.exception))
        self.assertEqual(ctx.exception.code, "report_error")

    def test_not_a_directory(self):
        with self.assertRaises(ReportError):
            report.validate_run_dir(os.path.join(self.run_dir, "nope"))

    def test_warnings_list_every_missing_companion_file(self):
        _write_summary(self.run_dir, [_ok_row("sample-a")])
        warnings = report.validate_run_dir(self.run_dir)
        self.assertEqual(len(warnings), 2)
        self.assertTrue(any("config.yml" in w for w in warnings))
        self.assertTrue(any("versions.tsv" in w for w in warnings))


class Test_LoadRunConfig(_RunDir):
    def test_absent_is_empty(self):
        self.assertEqual(report.load_run_config(self.run_dir), {})

    def test_mapping_is_returned(self):
        self.write_config(data="nanopore", minimum_depth=30)
        self.assertEqual(report.load_run_config(self.run_dir)["minimum_depth"], 30)

    def test_non_mapping_is_a_configuration_error(self):
        with open(os.path.join(self.run_dir, "config.yml"), "w") as fh:
            fh.write("- just\n- a list\n")
        with self.assertRaises(ConfigurationError):
            report.load_run_config(self.run_dir)


class Test_ReportArgv(_RunDir):
    def _argv(self, config, label=None):
        return report.report_argv(self.run_dir, config, "/out/report.html", label)

    def _value(self, argv, flag):
        return argv[argv.index(flag) + 1]

    def test_full_nanopore_config(self):
        scheme = os.path.join(self._tmp.name, "primers.bed")
        open(scheme, "w").close()
        argv = self._argv(
            {
                "data": "nanopore",
                "minimum_depth": 30,
                "consensus_coverage_threshold": 80,
                "reference": {"S": "/s.fa", "L": "/l.fa"},
                "scheme": scheme,
                "af_threshold": 0.6,
                "minimum_length": 200,
                "minimum_map_quality": 30,
                "clair3_model": "r941_prom_hac_g360+g422",
            }
        )
        self.assertEqual(self._value(argv, "--data-type"), "nanopore")
        self.assertEqual(self._value(argv, "--min-depth"), "30")
        self.assertEqual(self._value(argv, "--consensus-coverage-threshold"), "80.0")
        i = argv.index("--segments")
        self.assertEqual(argv[i + 1 : i + 3], ["S", "L"])
        self.assertEqual(self._value(argv, "--scheme"), scheme)
        self.assertEqual(self._value(argv, "--clair3-model"), "r941_prom_hac_g360+g422")
        self.assertEqual(self._value(argv, "--label"), "run1")
        self.assertEqual(self._value(argv, "--template"), report.TEMPLATE_PATH)
        self.assertTrue(os.path.isfile(report.TEMPLATE_PATH))

    def test_defaults_without_config_and_label_override(self):
        argv = self._argv({}, label="My run")
        self.assertEqual(self._value(argv, "--data-type"), "illumina")
        self.assertEqual(self._value(argv, "--min-depth"), "20")
        self.assertEqual(self._value(argv, "--consensus-coverage-threshold"), "70.0")
        self.assertEqual(self._value(argv, "--label"), "My run")
        for flag in ("--segments", "--scheme", "--clair3-model", "--minimum-map-quality"):
            self.assertNotIn(flag, argv)

    def test_data_type_inferred_from_clair3_layout(self):
        os.makedirs(os.path.join(self.run_dir, "assembly", "clair3", "sample-a"))
        self.assertEqual(self._value(self._argv({}), "--data-type"), "nanopore")

    def test_per_sample_model_dict_and_na_scheme_are_skipped(self):
        argv = self._argv({"clair3_model": {"sample-a": "m"}, "scheme": "NA"})
        self.assertNotIn("--clair3-model", argv)
        self.assertNotIn("--scheme", argv)

    def test_missing_scheme_file_is_left_out(self):
        argv = self._argv({"scheme": "/nowhere/primers.bed"})
        self.assertNotIn("--scheme", argv)

    def test_every_parser_flag_is_driven_by_both_callers(self):
        """The CLI and rule report must not drift from build_report.py: every
        option the script declares (bar the debugging waiver) is emitted by
        report_argv for a full config and referenced in rules/report.smk."""
        declared = set()
        for action in build_report.build_parser()._actions:
            declared.update(o for o in action.option_strings if o.startswith("--"))
        declared -= {"--help", "--no-self-check"}
        scheme = os.path.join(self._tmp.name, "primers.bed")
        open(scheme, "w").close()
        argv = self._argv(
            {
                "data": "nanopore",
                "reference": {"S": "/s.fa"},
                "scheme": scheme,
                "af_threshold": 0.5,
                "minimum_length": 50,
                "minimum_map_quality": 30,
                "clair3_model": "m",
            }
        )
        self.assertEqual(declared - set(argv), set(), "flags missing from report_argv")
        with open(os.path.join(RULES_DIR, "report.smk")) as fh:
            smk = fh.read()
        in_rule = set(re.findall(r"--[a-z][a-z0-9-]+", smk))
        self.assertEqual(declared - in_rule, set(), "flags missing from rules/report.smk")


class Test_CreateReport(_RunDir):
    def test_builds_a_page_from_a_spine_only_run(self):
        _write_summary(self.run_dir, [_ok_row("sample-a"), _ok_row("sample-b", cov=40.0)])
        self.write_config(data="illumina", minimum_depth=20, consensus_coverage_threshold=70)
        out = report.create_report(self.run_dir)
        self.assertEqual(out, os.path.join(self.run_dir, "report.html"))
        with open(out) as fh:
            page = fh.read()
        self.assertIn('"sample-a"', page)
        self.assertNotIn("/*{{DATA}}*/", page)
        self.assertIsNone(re.search(r'(src|href)="http', page))

    def test_output_and_label_are_honoured_without_config(self):
        _write_summary(self.run_dir, [_ok_row("sample-a")])
        out = os.path.join(self._tmp.name, "elsewhere.html")
        self.assertEqual(report.create_report(self.run_dir, output=out, label="Lab run"), out)
        with open(out) as fh:
            self.assertIn("Lab run", fh.read())

    def test_self_check_failure_is_a_report_error(self):
        """A pooled cov<T> FASTA that disagrees with summary.tsv fails the build."""
        _write_summary(self.run_dir, [_ok_row("sample-a"), _ok_row("sample-b", cov=40.0)])
        self.write_config(data="illumina", minimum_depth=20, consensus_coverage_threshold=70)
        os.makedirs(os.path.join(self.run_dir, "consensus"))
        with open(os.path.join(self.run_dir, "consensus", "consensus.cov70.fasta"), "w") as fh:
            fh.write(">sample-a\nACGT\n>sample-b\nACGT\n")
        with self.assertRaises(ReportError) as ctx:
            report.create_report(self.run_dir)
        self.assertIn("report build failed", str(ctx.exception))

    def test_missing_summary_is_reported(self):
        with self.assertRaises(ReportError):
            report.create_report(self.run_dir)


if __name__ == "__main__":
    unittest.main()
