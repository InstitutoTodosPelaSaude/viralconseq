"""Tests for argument validation helpers (viralconseq.validators).

Covers numeric range checks, sample-sheet integrity and identifier sanitization.
"""

import csv
import os
import tempfile
import unittest
from unittest.mock import patch

from viralconseq.constants import ViralQCDatabase
from viralconseq.exceptions import (
    ConfigurationError,
    SampleSheetError,
    ValidationError,
    ViralQCDatabaseNotFoundError,
)
from viralconseq.validators import (
    absolutise_sample_paths,
    ensure_within_base,
    get_samples_from_args,
    missing_viralqc_database_files,
    resolve_resource_budget,
    sanitize_identifier,
    validate_config_dict,
    validate_flag_strings,
    validate_numeric_parameters,
    validate_sample_sheet,
    validate_viralqc_database,
)


class TestValidateNumericParameters(unittest.TestCase):
    def test_valid_values_pass(self):
        args = {
            "threads": 4,
            "threads_total": 8,
            "af_threshold": 0.5,
            "minimum_coverage": 20,
        }
        # Should not raise.
        self.assertIsNone(validate_numeric_parameters(args))

    def test_absent_keys_pass(self):
        self.assertIsNone(validate_numeric_parameters({"data_type": "illumina"}))

    def test_zero_threads_rejected(self):
        with self.assertRaises(ValidationError):
            validate_numeric_parameters({"threads": 0})

    def test_negative_threads_total_rejected(self):
        with self.assertRaises(ValidationError):
            validate_numeric_parameters({"threads_total": -1})

    def test_af_threshold_above_one_rejected(self):
        with self.assertRaises(ValidationError):
            validate_numeric_parameters({"af_threshold": 5.0})


def _touch(path: str) -> str:
    with open(path, "w") as fh:
        fh.write("")
    return path


class Test_SampleSheetIntegrity(unittest.TestCase):
    """Guardrails against silent sample-sheet data loss / corruption."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.r1 = _touch(os.path.join(self.tmp, "s_R1.fastq.gz"))
        self.r2 = _touch(os.path.join(self.tmp, "s_R2.fastq.gz"))
        self.np = _touch(os.path.join(self.tmp, "s.fastq.gz"))

    def tearDown(self):
        self._tmp.cleanup()

    def _sheet(self, rows):
        path = os.path.join(self.tmp, "sheet.csv")
        with open(path, "w", newline="") as fh:
            csv.writer(fh).writerows(rows)
        return path

    def test_duplicate_sample_ids_rejected(self):
        """Two rows with the same id must error, not silently collapse to one."""
        sheet = self._sheet([["dup", self.r1, self.r2], ["dup", self.r1, self.r2]])
        with self.assertRaises(SampleSheetError):
            validate_sample_sheet(sheet, "illumina")

    def test_ragged_row_rejected_not_crash(self):
        """A short/ragged row raises SampleSheetError (not a TypeError from NaN)."""
        sheet = self._sheet([["s1", self.r1, self.r2], ["s2", self.r1]])
        with self.assertRaises(SampleSheetError):
            validate_sample_sheet(sheet, "illumina")

    def test_extra_column_rejected(self):
        """A nanopore sheet with an Illumina-style 3rd column is rejected."""
        sheet = self._sheet([["s1", self.np, "extra"]])
        with self.assertRaises(SampleSheetError):
            validate_sample_sheet(sheet, "nanopore")

    def test_trailing_comma_tolerated(self):
        """A stray trailing empty field is tolerated, not treated as a real column."""
        sheet = self._sheet([["s1", self.r1, self.r2, ""]])
        samples = validate_sample_sheet(sheet, "illumina")
        self.assertEqual(samples, {"s1": [self.r1, self.r2]})

    def test_sample_named_NA_kept_as_string(self):
        """A sample literally named 'NA' stays a string (pandas would coerce to NaN)."""
        sheet = self._sheet([["NA", self.np]])
        samples = validate_sample_sheet(sheet, "nanopore")
        self.assertEqual(list(samples.keys()), ["NA"])

    def test_numeric_sample_name_kept_as_string(self):
        """A numeric-looking id like '001' is not coerced to an int."""
        sheet = self._sheet([["001", self.np]])
        samples = validate_sample_sheet(sheet, "nanopore")
        self.assertEqual(list(samples.keys()), ["001"])

    def test_empty_sample_name_rejected(self):
        sheet = self._sheet([["", self.np]])
        with self.assertRaises(SampleSheetError):
            validate_sample_sheet(sheet, "nanopore")

    def test_unsafe_sample_name_rejected(self):
        """A sample id with a path separator must be rejected (injection surface)."""
        sheet = self._sheet([["../evil", self.np]])
        with self.assertRaises(SampleSheetError):
            validate_sample_sheet(sheet, "nanopore")


class Test_AbsolutiseSamplePaths(unittest.TestCase):
    """Sample FASTQ paths must be absolute in the config: Snakemake runs inside
    the run directory, so a relative path would resolve against the wrong base."""

    def test_relative_paths_resolved_against_cwd(self):
        samples = {"s1": ["reads/a_R1.fastq.gz", "reads/a_R2.fastq.gz"]}
        out = absolutise_sample_paths(samples)
        self.assertEqual(out["s1"][0], os.path.abspath("reads/a_R1.fastq.gz"))
        self.assertTrue(all(os.path.isabs(p) for p in out["s1"]))

    def test_absolute_paths_untouched_and_base_dir_honoured(self):
        samples = {"s1": ["/abs/a.fastq.gz"], "s2": ["b.fastq.gz"]}
        out = absolutise_sample_paths(samples, base_dir="/base")
        self.assertEqual(out["s1"], ["/abs/a.fastq.gz"])
        self.assertEqual(out["s2"], [os.path.abspath("/base/b.fastq.gz")])

    def test_get_samples_from_args_returns_absolute_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            old = os.getcwd()
            os.chdir(tmp)
            try:
                _touch("s.fastq.gz")
                with open("sheet.csv", "w") as fh:
                    fh.write("s1,s.fastq.gz\n")
                samples = get_samples_from_args(
                    {"sample_sheet": "sheet.csv", "data_type": "nanopore"}
                )
            finally:
                os.chdir(old)
        self.assertEqual(samples["s1"], [os.path.join(os.path.realpath(tmp), "s.fastq.gz")])


def _good_config(data="nanopore", **overrides):
    config = {
        "samples": {"sample-a": ["/r/a.fastq.gz"]},
        "data": data,
        "output": "/out/run/",
        "threads": 2,
        "reference": "/ref.fa",
        "scheme": "NA",
        "minimum_depth": 20,
        "minimum_length": 50,
        "af_threshold": 0.51,
        "run_viralqc": True,
        "viralqc_db": "/db",
        "run_viralqc_ram": 1,
        "viralqc_extra_flags": "",
        "minimap2_consensus_align_flags": "-a --sam-hit-only",
    }
    if data == "nanopore":
        config.update(
            chunk_size=10000,
            clair3_model="r941_prom_hac_g360+g422",
            variant_quality=20,
            variant_depth=10,
            minimum_map_quality=30,
            infer_consensus_sequence_ram=2,
        )
    else:
        config["samples"] = {"sample-a": ["/r/a_R1.fastq.gz", "/r/a_R2.fastq.gz"]}
        config.update(
            adapters="NA",
            trim_head=0,
            trim_tail=0,
            cut_front_mean_quality=10,
            cut_tail_mean_quality=10,
            cut_right_window_size=4,
            cut_right_mean_quality=15,
            af_isnv_threshold=0.0,
            run_isnv=False,
        )
    config.update(overrides)
    return config


class Test_ValidateConfigDict(unittest.TestCase):
    def test_generated_shapes_pass(self):
        validate_config_dict(_good_config("nanopore"))
        validate_config_dict(_good_config("illumina"))
        validate_config_dict(_good_config("illumina", reference={"S": "/S.fa", "L": "/L.fa"}))
        # legacy space-joined sample form and no viralQC
        validate_config_dict(
            _good_config("illumina", samples={"s": "/a_R1.fq /a_R2.fq"}, run_viralqc=False)
        )

    def test_rejections(self):
        cases = {
            "not a mapping": ["x"],
            "bad data": _good_config(data="pacbio"),
            "missing key": {k: v for k, v in _good_config().items() if k != "minimum_depth"},
            "missing ram of memory rule": {
                k: v for k, v in _good_config().items() if k != "infer_consensus_sequence_ram"
            },
            "wrong file count": _good_config(samples={"s": ["/a.fq", "/b.fq"]}),
            "empty samples": _good_config(samples={}),
            "bool as number": _good_config(minimum_depth=True),
            "float where int": _good_config(chunk_size=1.5),
            "af above one": _good_config(af_threshold=1.5),
            "zero threads": _good_config(threads=0),
            "cpus below one": _good_config(map_reads_cpus=0),
            "run flag not bool": _good_config(run_viralqc="yes"),
            "bad flag string": _good_config(viralqc_extra_flags="--x; rm -rf /"),
            "budget below rule": _good_config(max_memory_mb=1024),
            "empty model": _good_config(clair3_model=""),
            "segment without path": _good_config(reference={"S": ""}),
        }
        for label, config in cases.items():
            with self.subTest(case=label):
                with self.assertRaises(ConfigurationError):
                    validate_config_dict(config)

    def test_budget_at_or_above_largest_rule_accepted(self):
        validate_config_dict(_good_config(max_memory_mb=2048))
        validate_config_dict(_good_config(max_memory_mb=0))

    def test_message_names_the_key(self):
        with self.assertRaises(ConfigurationError) as ctx:
            validate_config_dict(_good_config(variant_depth=-1))
        self.assertIn("variant_depth", str(ctx.exception))


class Test_ResolveResourceBudget(unittest.TestCase):
    RULES = ["map_reads", "infer_consensus_sequence", "run_viralqc"]

    def test_explicit_budget_wins_and_is_converted_to_mb(self):
        args = {"max_memory": 8}
        with patch("viralconseq.validators.detect_memory_mb", return_value=32768):
            source = resolve_resource_budget(args, self.RULES)
        self.assertEqual(source, "--max-memory")
        self.assertEqual(args["max_memory_mb"], 8192)
        self.assertEqual(args["memory_detected_mb"], 32768)

    def test_zero_disables_the_budget(self):
        args = {"max_memory": 0}
        with patch("viralconseq.validators.detect_memory_mb", return_value=None):
            resolve_resource_budget(args, self.RULES)
        self.assertEqual(args["max_memory_mb"], 0)

    def test_explicit_budget_below_largest_rule_refused(self):
        with self.assertRaises(ValidationError):
            resolve_resource_budget({"max_memory": 1}, self.RULES)  # Clair3 wants 2 GB
        # ...unless the operator lowered that rule's figure too.
        args = {"max_memory": 1, "infer_consensus_sequence_ram": 1}
        with patch("viralconseq.validators.detect_memory_mb", return_value=None):
            resolve_resource_budget(args, self.RULES)
        self.assertEqual(args["max_memory_mb"], 1024)

    def test_detected_budget_takes_headroom(self):
        args = {}
        with patch("viralconseq.validators.detect_memory_mb", return_value=32768):
            source = resolve_resource_budget(args, self.RULES)
        self.assertEqual(args["max_memory_mb"], 29491)
        self.assertIn("detected 32.0 GB", source)

    def test_small_machine_clamps_up_to_one_job(self):
        args = {}
        with patch("viralconseq.validators.detect_memory_mb", return_value=1024):
            source = resolve_resource_budget(args, self.RULES)
        self.assertEqual(args["max_memory_mb"], 2048)
        self.assertIn("clamped", source)

    def test_undetectable_memory_runs_without_budget(self):
        args = {}
        with patch("viralconseq.validators.detect_memory_mb", return_value=None):
            source = resolve_resource_budget(args, self.RULES)
        self.assertEqual(args["max_memory_mb"], 0)
        self.assertEqual(source, "undetectable")


class Test_FlagStrings(unittest.TestCase):
    """Config-only tool flags are interpolated unquoted into shell commands."""

    def test_plain_flags_accepted(self):
        validate_flag_strings(
            {
                "viralqc_extra_flags": "--blast-pident 75 --verbose",
                "minimap2_consensus_align_flags": "-a --sam-hit-only --secondary=no --score-N=0",
            }
        )
        validate_flag_strings({"viralqc_extra_flags": ""})
        validate_flag_strings({})

    def test_unbalanced_quote_rejected(self):
        with self.assertRaises(ValidationError):
            validate_flag_strings({"viralqc_extra_flags": "--x '"})

    def test_shell_metacharacters_rejected(self):
        for bad in ("--x; rm -rf /", "$(id)", "a | b", "--x > out", "a\nb", "`id`"):
            with self.subTest(value=bad):
                with self.assertRaises(ValidationError):
                    validate_flag_strings({"minimap2_consensus_align_flags": bad})

    def test_non_string_rejected(self):
        with self.assertRaises(ValidationError):
            validate_flag_strings({"viralqc_extra_flags": ["--x"]})


class Test_Sanitization(unittest.TestCase):
    """Untrusted-input guards for run names and output paths."""

    def test_sanitize_identifier_accepts_safe(self):
        self.assertEqual(sanitize_identifier("run_2026.01-A"), "run_2026.01-A")

    def test_sanitize_identifier_rejects_path_separator(self):
        for bad in ["../etc", "a/b", "a\\b", "..", ".", "", "  ", "a b", "a;rm -rf"]:
            with self.assertRaises(ValidationError):
                sanitize_identifier(bad, field="run name")

    def test_ensure_within_base_allows_child(self):
        with tempfile.TemporaryDirectory() as base:
            target = ensure_within_base("runs/x", base)
            self.assertTrue(target.startswith(os.path.abspath(base)))

    def test_ensure_within_base_rejects_traversal(self):
        with tempfile.TemporaryDirectory() as base:
            with self.assertRaises(ValidationError):
                ensure_within_base("../../etc/passwd", base)


def _make_viralqc_db(root, omit=()):
    """Create a complete viralQC database layout under *root*, minus *omit*."""
    db = os.path.join(root, "vqc")
    os.makedirs(db, exist_ok=True)
    for name in ViralQCDatabase.REQUIRED_FILES:
        if name not in omit:
            open(os.path.join(db, name), "w").close()
    for name in ViralQCDatabase.REQUIRED_DIRS:
        if name not in omit:
            os.makedirs(os.path.join(db, name), exist_ok=True)
    if ViralQCDatabase.BLAST_INDEX_LABEL not in omit:
        open(os.path.join(db, ViralQCDatabase.BLAST_INDEX_LABEL), "w").close()
    return db


class Test_ViralQCDatabase(unittest.TestCase):
    def test_missing_files_reports_absent_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(
                missing_viralqc_database_files(os.path.join(tmp, "nope")), ["<directory>"]
            )

    def test_missing_files_empty_when_complete(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(missing_viralqc_database_files(_make_viralqc_db(tmp)), [])

    def test_missing_files_lists_each_absent_entry(self):
        with tempfile.TemporaryDirectory() as tmp:
            db = _make_viralqc_db(tmp, omit=("blast.tsv", "blast_gff"))
            self.assertEqual(missing_viralqc_database_files(db), ["blast.tsv", "blast_gff/"])

    def test_validator_ok_when_complete(self):
        with tempfile.TemporaryDirectory() as tmp:
            validate_viralqc_database({"run_viralqc": True, "viralqc_db": _make_viralqc_db(tmp)})

    def test_validator_missing_dir_raises_with_setup_hint(self):
        with tempfile.TemporaryDirectory() as tmp:
            missing = os.path.join(tmp, "absent")
            with self.assertRaises(ViralQCDatabaseNotFoundError) as ctx:
                validate_viralqc_database({"run_viralqc": True, "viralqc_db": missing})
        self.assertIn(missing, str(ctx.exception))
        self.assertIn(f"viralconseq setup --viralqc-db {missing}", str(ctx.exception))
        self.assertEqual(ctx.exception.code, "viralqc_database_not_found")

    def test_missing_blast_index_is_reported(self):
        """An unindexed blast.fasta must not pass: blastn would fail and viralQC
        degrades that to "no hits", silently reporting every unmatched sequence
        as Unclassified."""
        with tempfile.TemporaryDirectory() as tmp:
            db = _make_viralqc_db(tmp, omit=(ViralQCDatabase.BLAST_INDEX_LABEL,))
            self.assertEqual(
                missing_viralqc_database_files(db), [ViralQCDatabase.BLAST_INDEX_LABEL]
            )
            with self.assertRaises(ViralQCDatabaseNotFoundError) as ctx:
                validate_viralqc_database({"run_viralqc": True, "viralqc_db": db})
        self.assertIn(ViralQCDatabase.BLAST_INDEX_LABEL, str(ctx.exception))

    def test_volume_split_blast_index_accepted(self):
        """A large database split into volumes has blast.fasta.<NN>.nin plus a
        blast.fasta.nal alias instead of a single blast.fasta.nin."""
        with tempfile.TemporaryDirectory() as tmp:
            db = _make_viralqc_db(tmp, omit=(ViralQCDatabase.BLAST_INDEX_LABEL,))
            open(os.path.join(db, "blast.fasta.00.nin"), "w").close()
            open(os.path.join(db, "blast.fasta.nal"), "w").close()
            self.assertEqual(missing_viralqc_database_files(db), [])

    def test_validator_missing_sentinel_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            db = _make_viralqc_db(tmp, omit=(ViralQCDatabase.NEXTCLADE_SENTINEL,))
            with self.assertRaises(ViralQCDatabaseNotFoundError) as ctx:
                validate_viralqc_database({"run_viralqc": True, "viralqc_db": db})
        self.assertIn(ViralQCDatabase.NEXTCLADE_SENTINEL, str(ctx.exception))

    def test_validator_skipped_when_run_viralqc_false(self):
        validate_viralqc_database({"run_viralqc": False, "viralqc_db": "/definitely/absent"})

    def test_validator_unset_db_raises(self):
        for value in (None, "", "NA"):
            with self.subTest(value=value):
                with self.assertRaises(ViralQCDatabaseNotFoundError):
                    validate_viralqc_database({"run_viralqc": True, "viralqc_db": value})

    def test_validator_defaults_run_viralqc_true_when_key_absent(self):
        with self.assertRaises(ViralQCDatabaseNotFoundError):
            validate_viralqc_database({"viralqc_db": "/definitely/absent"})


if __name__ == "__main__":
    unittest.main()
