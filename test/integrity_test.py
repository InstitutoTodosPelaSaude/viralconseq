"""Tests for content-level input-integrity validation (viralconseq.integrity)
and its orchestration in validators.validate_consensus_input_integrity.

The FASTQ check has two layers: ``integrity.scan_fastq`` is the pure-Python
diagnoser (every issue code lives there) and ``integrity.validate_fastq`` puts
``seqkit.detect`` in front of it. Scanner-behaviour tests drive ``scan_fastq``
directly, because quality-range and whitespace problems are blind spots for
seqkit; the end-to-end tests need the real binary and skip without it.
"""

import concurrent.futures
import gzip
import os
import random
import shutil
import tempfile
import unittest
from unittest import mock

from viralconseq import integrity, seqkit
from viralconseq.exceptions import InputIntegrityError, ValidationError
from viralconseq.validators import validate_consensus_input_integrity

_SEQKIT = shutil.which("seqkit")
_NO_SEQKIT = "seqkit is not on PATH: install it (bioconda) to run the seqkit-backed tests"


def _codes(report):
    return [issue.code for issue in report.errors]


def _warn_codes(report):
    return [issue.code for issue in report.warnings]


def _require_seqkit(test):
    if _SEQKIT is None:
        test.skipTest(_NO_SEQKIT)


class _TmpBase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def _write(self, name, text, mode="w"):
        path = os.path.join(self.tmp, name)
        with open(path, mode) as fh:
            fh.write(text)
        return path

    def _write_gz(self, name, text):
        path = os.path.join(self.tmp, name)
        with gzip.open(path, "wt") as fh:
            fh.write(text)
        return path

    def _write_truncated_gz(self, name, n_reads=500):
        # A gzip stream cut in half: the first byte decompresses fine (so the
        # probe defers to seqkit) but the stream ends before its trailer.
        rng = random.Random(0)
        reads = []
        for i in range(n_reads):
            seq = "".join(rng.choice("ACGT") for _ in range(100))
            reads.append(f"@read{i}\n{seq}\n+\n{'I' * 100}\n")
        blob = gzip.compress("".join(reads).encode())
        return self._write(name, blob[: len(blob) // 2], mode="wb")


class TestScanFastq(_TmpBase):
    """The pure-Python diagnoser: every issue code, no seqkit involved."""

    def test_valid_fastq(self):
        p = self._write("ok.fastq", "@r1\nACGT\n+\nIIII\n@r2\nACGTN\n+\nIIIII\n")
        self.assertEqual(_codes(integrity.scan_fastq(p)), [])

    def test_valid_gzipped_fastq(self):
        p = self._write_gz("ok.fastq.gz", "@r1\nACGT\n+\nIIII\n")
        self.assertEqual(_codes(integrity.scan_fastq(p)), [])

    def test_truncated_record(self):
        p = self._write("trunc.fastq", "@r1\nACGT\n+\nIIII\n@r2\nACGT\n")
        self.assertIn("fastq_truncated", _codes(integrity.scan_fastq(p)))

    def test_seq_qual_length_mismatch(self):
        p = self._write("mm.fastq", "@r1\nACGT\n+\nIII\n")
        self.assertIn("fastq_length_mismatch", _codes(integrity.scan_fastq(p)))

    def test_bad_sequence_char(self):
        p = self._write("badseq.fastq", "@r1\nACGZ\n+\nIIII\n")
        self.assertIn("fastq_bad_sequence_char", _codes(integrity.scan_fastq(p)))

    def test_bad_quality_char(self):
        # A space in the quality line is outside ASCII 33-126: a seqkit blind
        # spot, so it must be caught here.
        p = self._write("badqual.fastq", "@r1\nACGT\n+\nII I\n")
        self.assertIn("fastq_bad_quality_char", _codes(integrity.scan_fastq(p)))

    def test_whitespace_in_sequence_line(self):
        # Also invisible to seqkit when the lengths still match.
        p = self._write("space.fastq", "@r1\nAC GT\n+\nIIIII\n")
        self.assertIn("fastq_bad_sequence_char", _codes(integrity.scan_fastq(p)))

    def test_bad_header_marker(self):
        # A FASTA mislabeled as FASTQ: first line starts with '>'.
        p = self._write("wrong.fastq", ">r1\nACGT\n+\nIIII\n")
        self.assertIn("fastq_bad_header", _codes(integrity.scan_fastq(p)))

    def test_bad_separator(self):
        p = self._write("badsep.fastq", "@r1\nACGT\n*\nIIII\n")
        self.assertIn("fastq_bad_separator", _codes(integrity.scan_fastq(p)))

    def test_empty_file(self):
        p = self._write("empty.fastq", "")
        self.assertIn("fastq_empty", _codes(integrity.scan_fastq(p)))

    def test_corrupt_gzip_magic_but_garbage(self):
        # gzip magic bytes but not a valid archive -> unreadable, not a crash.
        p = self._write("bad.fastq.gz", b"\x1f\x8b not a real gzip stream", mode="wb")
        self.assertIn("fastq_unreadable", _codes(integrity.scan_fastq(p)))

    def test_truncated_gzip_unreadable(self):
        p = self._write_truncated_gz("cut.fastq.gz")
        self.assertIn("fastq_unreadable", _codes(integrity.scan_fastq(p)))

    def test_issue_cap(self):
        # Many mismatched records; issue count must be bounded by the cap.
        p = self._write("many.fastq", "@r\nAC\n+\nI\n" * 100)
        report = integrity.scan_fastq(p)
        self.assertLessEqual(len(report.issues), integrity.MAX_ISSUES_PER_FILE)

    def test_blank_lines_between_records_tolerated(self):
        # A stray blank line must not desync grouping or fake a truncation (F5).
        p = self._write("blanks.fastq", "@r1\nACGT\n+\nIIII\n\n@r2\nACGT\n+\nIIII\n")
        self.assertEqual(_codes(integrity.scan_fastq(p)), [])

    def test_truncation_reported_even_at_issue_cap(self):
        # 20 bad records + a truncated tail: fastq_truncated must still surface (F14).
        p = self._write("trunc_cap.fastq", "@r\nAC\n+\nI\n" * 25 + "@last\nACGT\n")
        codes = _codes(integrity.scan_fastq(p))
        self.assertIn("fastq_truncated", codes)

    def test_binary_file_reported_unreadable(self):
        p = self._write("bin.fastq", b"\x00\x01\x02\xff\xfe non-text bytes", mode="wb")
        self.assertIn("fastq_unreadable", _codes(integrity.scan_fastq(p)))

    def test_utf8_bom_fastq_tolerated(self):
        p = self._write("bom.fastq", "\ufeff@r1\nACGT\n+\nIIII\n")
        self.assertEqual(_codes(integrity.scan_fastq(p)), [])


class TestValidateFastqVerdictMapping(_TmpBase):
    """validate_fastq with seqkit.detect patched: how each verdict becomes a
    report, without the binary."""

    def setUp(self):
        super().setUp()
        self.path = self._write("in.fastq", "@r1\nACGT\n+\nIIII\n")

    def _run(self, verdict, message=""):
        result = seqkit.SeqkitResult(verdict=verdict, returncode=0, message=message)
        with mock.patch("viralconseq.integrity.seqkit.detect", return_value=result) as detect:
            with mock.patch("viralconseq.integrity.scan_fastq", wraps=integrity.scan_fastq) as scan:
                report = integrity.validate_fastq(self.path)
        detect.assert_called_once_with(self.path)
        return report, scan

    def test_clean_gives_empty_report_without_scanning(self):
        report, scan = self._run(seqkit.VERDICT_CLEAN)
        self.assertEqual(report.issues, [])
        self.assertEqual(report.kind, "fastq")
        self.assertEqual(report.path, self.path)
        scan.assert_not_called()

    def test_empty_gives_fastq_empty_without_scanning(self):
        report, scan = self._run(seqkit.VERDICT_EMPTY)
        self.assertEqual(_codes(report), ["fastq_empty"])
        scan.assert_not_called()

    def test_unreadable_gives_fastq_unreadable_with_seqkit_message(self):
        report, scan = self._run(seqkit.VERDICT_UNREADABLE, message="[ERRO] unexpected EOF")
        self.assertEqual(_codes(report), ["fastq_unreadable"])
        self.assertIn("unexpected EOF", report.errors[0].message)
        scan.assert_not_called()

    def test_suspect_hands_the_file_to_scan_fastq(self):
        report, scan = self._run(seqkit.VERDICT_SUSPECT, message="[ERRO] fastx: bad fastq format")
        scan.assert_called_once_with(self.path)
        # The file is actually fine, so the diagnoser finds no error; seqkit's
        # verdict is kept as a warning rather than dropped.
        self.assertEqual(_codes(report), [])
        self.assertEqual(_warn_codes(report), ["fastq_seqkit_flagged"])

    def test_suspect_without_scanner_findings_is_a_warning(self):
        p = self._write("ok.fastq", "@r1\nACGT\n+\nIIII\n")
        with mock.patch.object(
            integrity.seqkit, "detect", return_value=seqkit.SeqkitResult("suspect", 255, "[ERRO] x")
        ):
            with mock.patch.object(
                integrity, "scan_fastq", return_value=integrity.IntegrityReport(p, "fastq")
            ):
                report = integrity.validate_fastq(p)
        self.assertEqual(_codes(report), [])
        self.assertEqual(_warn_codes(report), ["fastq_seqkit_flagged"])
        self.assertIn("[ERRO] x", report.warnings[0].message)

    def test_signal_death_warning_does_not_claim_a_zero_length_read(self):
        p = self._write("ok.fastq", "@r1\nACGT\n+\nIIII\n")
        with mock.patch.object(
            integrity.seqkit, "detect", return_value=seqkit.SeqkitResult("suspect", -9, "")
        ):
            report = integrity.validate_fastq(p)
        self.assertEqual(_warn_codes(report), ["fastq_seqkit_flagged"])
        self.assertIn("exit status -9: no error message", report.warnings[0].message)
        self.assertNotIn("zero-length", report.warnings[0].message)


class TestValidateFastq(_TmpBase):
    """End to end with the real seqkit: detector verdict -> old issue codes."""

    def setUp(self):
        _require_seqkit(self)
        super().setUp()

    def test_valid_fastq_is_clean_without_scanning(self):
        p = self._write("ok.fastq", "@r1\nACGT\n+\nIIII\n@r2\nACGTN\n+\nIIIII\n")
        with mock.patch("viralconseq.integrity.scan_fastq") as scan:
            report = integrity.validate_fastq(p)
        scan.assert_not_called()
        self.assertEqual(report.issues, [])

    def test_valid_gzipped_fastq_is_clean_without_scanning(self):
        p = self._write_gz("ok.fastq.gz", "@r1\nACGT\n+\nIIII\n")
        with mock.patch("viralconseq.integrity.scan_fastq") as scan:
            report = integrity.validate_fastq(p)
        scan.assert_not_called()
        self.assertEqual(report.issues, [])

    def test_empty_file(self):
        p = self._write("empty.fastq", "")
        self.assertEqual(_codes(integrity.validate_fastq(p)), ["fastq_empty"])

    def test_empty_gzip(self):
        p = self._write_gz("empty.fastq.gz", "")
        self.assertEqual(_codes(integrity.validate_fastq(p)), ["fastq_empty"])

    def test_truncated_gzip_unreadable(self):
        p = self._write_truncated_gz("cut.fastq.gz")
        self.assertEqual(_codes(integrity.validate_fastq(p)), ["fastq_unreadable"])

    def test_corrupt_gzip_magic_but_garbage_unreadable(self):
        p = self._write("bad.fastq.gz", b"\x1f\x8b not a real gzip stream", mode="wb")
        self.assertEqual(_codes(integrity.validate_fastq(p)), ["fastq_unreadable"])

    def test_wrapped_record_is_a_bad_separator(self):
        # seqkit accepts a line-wrapped record but fastp emits zero reads for
        # it; the probe routes it to the scanner, which names the old code.
        p = self._write("wrapped.fastq", "@r1\nACGTACGTAC\nGTACGTACGT\n+\nIIIIIIIIII\nIIIIIIIIII\n")
        self.assertIn("fastq_bad_separator", _codes(integrity.validate_fastq(p)))

    def test_truncated_record(self):
        p = self._write("trunc.fastq", "@r1\nACGT\n+\nIIII\n@r2\nACGT\n")
        self.assertEqual(_codes(integrity.validate_fastq(p)), ["fastq_truncated"])

    def test_truncated_record_in_gzip(self):
        p = self._write_gz("trunc.fastq.gz", "@r1\nACGT\n+\nIIII\n@r2\nACGT\n")
        self.assertEqual(_codes(integrity.validate_fastq(p)), ["fastq_truncated"])

    def test_seq_qual_length_mismatch(self):
        p = self._write("mm.fastq", "@r1\nACGT\n+\nIII\n")
        self.assertEqual(_codes(integrity.validate_fastq(p)), ["fastq_length_mismatch"])

    def test_bad_header_marker(self):
        p = self._write("wrong.fastq", ">r1\nACGT\n+\nIIII\n")
        self.assertEqual(_codes(integrity.validate_fastq(p)), ["fastq_bad_header"])

    def test_real_fasta_given_as_fastq_is_rejected(self):
        # seqkit reads a genuine FASTA (no "+" line) as FASTA and says nothing;
        # the probe routes it to the scanner, which names the header problem.
        fasta = ">r1\nACGT\n>r2\nACGTN\n"
        plain = self._write("reads.fastq", fasta)
        gz = os.path.join(self.tmp, "reads.fastq.gz")
        with gzip.open(gz, "wt") as fh:
            fh.write(fasta)
        for path in (plain, gz):
            with self.subTest(path=os.path.basename(path)):
                self.assertIn("fastq_bad_header", _codes(integrity.validate_fastq(path)))

    def test_cut_right_after_record_marker_is_truncated(self):
        # "@" alone at the end: seqkit prints a bare newline for the empty id.
        p = self._write("cut.fastq", "@r1\nACGT\n+\nIIII\n@")
        self.assertEqual(_codes(integrity.validate_fastq(p)), ["fastq_truncated"])

    def test_diagnosed_issue_carries_line_number(self):
        # The diagnoser's report is returned as is: same message shape as before.
        p = self._write("mm.fastq", "@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIII\n")
        report = integrity.validate_fastq(p)
        self.assertEqual(_codes(report), ["fastq_length_mismatch"])
        self.assertEqual(report.errors[0].line, 5)

    def test_quality_blind_spot_is_documented(self):
        # seqkit cannot see a space in the quality line, so validate_fastq
        # reports nothing for it. Pinned so a seqkit that starts catching it
        # (or a change to the detector command) is noticed.
        p = self._write("badqual.fastq", "@r1\nACGT\n+\nII I\n")
        self.assertEqual(_codes(integrity.validate_fastq(p)), [])


def _report(path, *codes):
    report = integrity.IntegrityReport(path=path, kind="fastq")
    for code in codes:
        report.issues.append(integrity.IntegrityIssue("error", code, code, path=path, kind="fastq"))
    return report


class TestValidateFastqFiles(_TmpBase):
    def _mixed_fixtures(self):
        return [
            self._write("ok.fastq", "@r1\nACGT\n+\nIIII\n"),
            self._write("trunc.fastq", "@r1\nACGT\n+\nIIII\n@r2\nACGT\n"),
            self._write("empty.fastq", ""),
            self._write_gz("ok.fastq.gz", "@r1\nACGT\n+\nIIII\n"),
            self._write("mm.fastq", "@r1\nACGT\n+\nIII\n"),
            self._write("bad.fastq.gz", b"\x1f\x8b not a real gzip stream", mode="wb"),
        ]

    def test_pool_matches_serial_in_input_order_with_duplicate(self):
        _require_seqkit(self)
        fixtures = self._mixed_fixtures()
        paths = fixtures + [fixtures[1]]  # the truncated file listed twice
        serial = integrity.validate_fastq_files(paths, workers=1)
        pooled = integrity.validate_fastq_files(paths, workers=4)
        key = lambda reports: [(r.path, _codes(r)) for r in reports]  # noqa: E731
        self.assertEqual(key(pooled), key(serial))
        self.assertEqual([r.path for r in pooled], paths)
        self.assertEqual(
            [_codes(r) for r in pooled],
            [
                [],
                ["fastq_truncated"],
                ["fastq_empty"],
                [],
                ["fastq_length_mismatch"],
                ["fastq_unreadable"],
                ["fastq_truncated"],
            ],
        )

    def test_single_file_never_builds_a_pool(self):
        p = self._write("ok.fastq", "@r1\nACGT\n+\nIIII\n")
        with mock.patch("viralconseq.integrity.validate_fastq", side_effect=_report) as vf:
            with mock.patch("viralconseq.integrity.concurrent.futures.ProcessPoolExecutor") as pool:
                reports = integrity.validate_fastq_files([p], workers=8)
        pool.assert_not_called()
        vf.assert_called_once_with(p)
        self.assertEqual([r.path for r in reports], [p])

    def test_one_worker_never_builds_a_pool(self):
        a = self._write("a.fastq", "@r1\nACGT\n+\nIIII\n")
        b = self._write("b.fastq", "@r1\nACGT\n+\nIIII\n")
        with mock.patch("viralconseq.integrity.validate_fastq", side_effect=_report) as vf:
            with mock.patch("viralconseq.integrity.concurrent.futures.ProcessPoolExecutor") as pool:
                reports = integrity.validate_fastq_files([a, b, a], workers=1)
        pool.assert_not_called()
        self.assertEqual([c.args[0] for c in vf.call_args_list], [a, b])
        self.assertEqual([r.path for r in reports], [a, b, a])

    def test_workers_capped_at_distinct_file_count(self):
        a = self._write("a.fastq", "@r1\nACGT\n+\nIIII\n")
        b = self._write("b.fastq", "@r1\nACGT\n+\nIIII\n")
        with mock.patch("viralconseq.integrity.concurrent.futures.ProcessPoolExecutor") as pool:
            pool.return_value.__enter__.return_value.map.side_effect = lambda fn, items: [
                _report(p) for p in items
            ]
            reports = integrity.validate_fastq_files([a, b, a, b], workers=16)
        pool.assert_called_once_with(max_workers=2)
        self.assertEqual([r.path for r in reports], [a, b, a, b])

    def test_broken_pool_falls_back_to_serial(self):
        a = self._write("a.fastq", "@r1\nACGT\n+\nIIII\n")
        b = self._write("b.fastq", "@r1\nACGT\n+\nIIII\n")
        with mock.patch("viralconseq.integrity.validate_fastq", side_effect=_report) as vf:
            with mock.patch("viralconseq.integrity.concurrent.futures.ProcessPoolExecutor") as pool:
                pool.return_value.__enter__.return_value.map.side_effect = (
                    concurrent.futures.process.BrokenProcessPool("worker died")
                )
                with self.assertLogs("viralconseq.integrity", "WARNING") as logs:
                    reports = integrity.validate_fastq_files([a, b], workers=2)
        pool.assert_called_once()
        self.assertEqual([c.args[0] for c in vf.call_args_list], [a, b])
        self.assertEqual([r.path for r in reports], [a, b])
        self.assertIn("serially", "\n".join(logs.output))

    def test_worker_exception_other_than_a_dead_pool_propagates(self):
        # Only BrokenProcessPool triggers the serial fallback; a real error from
        # a worker (here: no seqkit) must surface, not be retried serially.
        a = self._write("a.fastq", "@r1\nACGT\n+\nIIII\n")
        b = self._write("b.fastq", "@r1\nACGT\n+\nIIII\n")
        with mock.patch("viralconseq.integrity.validate_fastq", side_effect=_report) as vf:
            with mock.patch("viralconseq.integrity.concurrent.futures.ProcessPoolExecutor") as pool:
                pool.return_value.__enter__.return_value.map.side_effect = (
                    seqkit.SeqkitMissingError()
                )
                with self.assertRaises(seqkit.SeqkitMissingError):
                    integrity.validate_fastq_files([a, b], workers=2)
        vf.assert_not_called()


class TestValidateFasta(_TmpBase):
    def test_valid_reference(self):
        p = self._write("r.fasta", ">MN908947.3 desc\nACGTACGTNN\nACGT\n")
        report = integrity.validate_fasta(p)
        self.assertEqual(_codes(report), [])
        self.assertEqual(report.contig_ids, ["MN908947.3"])
        self.assertEqual(report.headers, ["MN908947.3 desc"])

    def test_protein_fasta_blocked(self):
        p = self._write("p.fasta", ">seq1\nMKLVEEFG\n")
        self.assertIn("fasta_non_nucleotide", _codes(integrity.validate_fasta(p)))

    def test_gaps_rejected(self):
        p = self._write("gap.fasta", ">seq1\nACGT--ACGT\n")
        self.assertIn("fasta_non_nucleotide", _codes(integrity.validate_fasta(p)))

    def test_duplicate_contig_ids(self):
        p = self._write("dup.fasta", ">seg\nACGT\n>seg\nTTTT\n")
        self.assertIn("fasta_duplicate_id", _codes(integrity.validate_fasta(p)))

    def test_empty_record(self):
        p = self._write("empty_rec.fasta", ">seg1\n>seg2\nACGT\n")
        self.assertIn("fasta_empty_record", _codes(integrity.validate_fasta(p)))

    def test_no_records(self):
        p = self._write("none.fasta", "not a fasta at all\n")
        codes = _codes(integrity.validate_fasta(p))
        self.assertIn("fasta_junk_before_header", codes)
        self.assertIn("fasta_no_records", codes)

    def test_require_nucleotide_false_allows_protein(self):
        p = self._write("p.fasta", ">seq1\nMKLVEEFG\n")
        self.assertEqual(_codes(integrity.validate_fasta(p, require_nucleotide=False)), [])

    def test_pipe_headers_keep_full_token_and_are_not_duplicates(self):
        # minimap2/samtools split on whitespace only, so 'gi|...' contig ids are
        # distinct and the full pipe token is the contig name (F1).
        p = self._write(
            "gi.fasta",
            ">gi|1|ref|NC_A.1| virus A\nACGT\n>gi|2|ref|NC_B.1| virus B\nTTTT\n",
        )
        report = integrity.validate_fasta(p)
        self.assertEqual(_codes(report), [])
        self.assertEqual(report.contig_ids, ["gi|1|ref|NC_A.1|", "gi|2|ref|NC_B.1|"])

    def test_utf8_bom_reference_tolerated(self):
        p = self._write("bom.fasta", "﻿>MN908947.3\nACGTACGT\n")
        report = integrity.validate_fasta(p)
        self.assertEqual(_codes(report), [])
        self.assertEqual(report.contig_ids, ["MN908947.3"])


class TestValidateBed(_TmpBase):
    def test_valid_bed_matching_contig(self):
        p = self._write("p.bed", "MN908947.3\t0\t20\tprimer_1\t.\t+\n")
        report = integrity.validate_bed(p, expected_contigs={"MN908947.3"})
        self.assertEqual(_codes(report), [])

    def test_too_few_columns(self):
        p = self._write("p.bed", "MN908947.3\t0\n")
        self.assertIn("bed_too_few_columns", _codes(integrity.validate_bed(p)))

    def test_bad_interval(self):
        p = self._write("p.bed", "chrom\t30\t10\tp1\n")
        self.assertIn("bed_bad_interval", _codes(integrity.validate_bed(p)))

    def test_non_integer_coord(self):
        p = self._write("p.bed", "chrom\tx\ty\tp1\n")
        self.assertIn("bed_non_integer_coord", _codes(integrity.validate_bed(p)))

    def test_chrom_mismatch_blocks_when_none_match(self):
        p = self._write("p.bed", "OTHER\t0\t20\tp1\n")
        self.assertIn(
            "bed_chrom_mismatch", _codes(integrity.validate_bed(p, expected_contigs={"MN908947.3"}))
        )

    def test_partial_chrom_match_downgrades_to_warning(self):
        # Whole-scheme BED (S + L) run against a subset reference (only S): S is
        # trimmed, L rows are harmless -> warning, not a blocking error (F4).
        p = self._write("p.bed", "S\t0\t5\tp1\nL\t0\t5\tp2\n")
        report = integrity.validate_bed(p, expected_contigs={"S"})
        self.assertEqual(_codes(report), [])
        self.assertIn("bed_chrom_unmatched", _warn_codes(report))

    def test_pipe_chrom_matches_reference(self):
        p = self._write("p.bed", "gi|1|ref|NC_A.1|\t0\t20\tp1\n")
        report = integrity.validate_bed(p, expected_contigs={"gi|1|ref|NC_A.1|"})
        self.assertEqual(_codes(report), [])

    def test_contig_named_track_not_dropped(self):
        # A real data row whose chrom starts with 'track' must not be skipped (F6).
        p = self._write("p.bed", "trackpox\t0\t20\tp1\n")
        report = integrity.validate_bed(p, expected_contigs={"trackpox"})
        self.assertEqual(_codes(report), [])
        self.assertNotIn("bed_no_features", _codes(report))

    def test_comment_and_track_lines_skipped(self):
        p = self._write("p.bed", "# a comment\ntrack name=x\nMN908947.3\t0\t20\tp1\n")
        self.assertEqual(_codes(integrity.validate_bed(p, expected_contigs={"MN908947.3"})), [])

    def test_empty_bed(self):
        p = self._write("p.bed", "# only comments\n")
        self.assertIn("bed_no_features", _codes(integrity.validate_bed(p)))

    def test_odd_strand_is_warning(self):
        p = self._write("p.bed", "chrom\t0\t20\tp1\t.\tX\n")
        report = integrity.validate_bed(p)
        self.assertIn("bed_bad_strand", _warn_codes(report))
        self.assertNotIn("bed_bad_strand", _codes(report))

    def test_nanopore_accession_vs_sanitized_fixit(self):
        # BED chrom matches the accession, but nanopore maps against the sanitized
        # full header, so trimming would silently do nothing -> block with fix-it.
        p = self._write("p.bed", "MN908947.3\t0\t20\tp1\n")
        report = integrity.validate_bed(
            p,
            expected_contigs={"MN908947.3_Severe_acute"},
            accession_map={"MN908947.3": "MN908947.3_Severe_acute"},
        )
        codes = _codes(report)
        self.assertIn("bed_chrom_sanitized_mismatch", codes)
        self.assertIn("MN908947.3_Severe_acute", report.errors[0].message)


INFLUENZA_HEADER = "Influenza A virus (A/New York/392/2004(H3N2)) segment 5, complete sequence"
INFLUENZA_CONTIG = "Influenza_A_virus__A_New_York_392_2004_H3N2___segment_5__complete_sequence"


class TestSanitizeNanoporeContig(unittest.TestCase):
    """Pins the Python copy of ``rule sanitize_reference``'s sed."""

    def test_shell_metacharacters_become_underscores(self):
        self.assertEqual(integrity.sanitize_nanopore_contig(INFLUENZA_HEADER), INFLUENZA_CONTIG)

    def test_legacy_character_set_is_unchanged(self):
        # The pre-0.3 rule rewrote only / \ | , ~ and space; such headers must keep
        # the contig names earlier runs produced.
        self.assertEqual(
            integrity.sanitize_nanopore_contig("MN908947.3 Severe acute|x~y,z\\w"),
            "MN908947.3_Severe_acute_x_y_z_w",
        )

    def test_one_underscore_per_utf8_byte(self):
        # Ñ is two bytes in UTF-8, as sed sees it under LC_ALL=C.
        self.assertEqual(integrity.sanitize_nanopore_contig("Ñ:a>b"), "___a_b")

    def test_leading_marker_and_whitespace_are_stripped(self):
        self.assertEqual(integrity.sanitize_nanopore_contig("  MN908947.3 x \n"), "MN908947.3_x")
        # CRLF and trailing spaces are trimmed, as the sed trims them; the header
        # is the text after the marker, so a second ``>`` is content (the sed
        # strips one) and a NBSP is not ASCII whitespace.
        self.assertEqual(integrity.sanitize_nanopore_contig("MN908947.3\r\n"), "MN908947.3")
        self.assertEqual(integrity.sanitize_nanopore_contig("MN908947.3 \t"), "MN908947.3")
        self.assertEqual(integrity.sanitize_nanopore_contig(">abc\tq"), "_abc_q")
        self.assertEqual(integrity.sanitize_nanopore_contig("a\u00a0"), "a__")


class TestOrchestration(_TmpBase):
    """validate_consensus_input_integrity end to end; needs the real seqkit."""

    def setUp(self):
        _require_seqkit(self)
        super().setUp()

    def _reference(self, header=">MN908947.3\n"):
        return self._write("ref.fasta", header + "ACGTACGTAC\n")

    def _fastq(self):
        return self._write("s1.fastq", "@r1\nACGT\n+\nIIII\n")

    def test_nanopore_metacharacter_header_warns_and_bed_uses_sanitized_name(self):
        args = {
            "data_type": "nanopore",
            "reference": self._reference(f">{INFLUENZA_HEADER}\n"),
            "primer_scheme": self._write("p.bed", f"{INFLUENZA_CONTIG}\t0\t5\tp1\n"),
        }
        with self.assertLogs("viralconseq.validators", "WARNING") as logs:
            validate_consensus_input_integrity(args, {"s1": [self._fastq()]})
        joined = "\n".join(logs.output)
        self.assertIn("fasta_header_sanitized", joined)
        self.assertIn(INFLUENZA_CONTIG, joined)

    def test_nanopore_metacharacter_header_with_accession_bed_blocks(self):
        args = {
            "data_type": "nanopore",
            "reference": self._reference(f">{INFLUENZA_HEADER}\n"),
            "primer_scheme": self._write("p.bed", "Influenza\t0\t5\tp1\n"),
        }
        with self.assertRaises(InputIntegrityError) as ctx:
            validate_consensus_input_integrity(args, {"s1": [self._fastq()]})
        issue = ctx.exception.issues[0]
        self.assertEqual(issue.code, "bed_chrom_sanitized_mismatch")
        self.assertIn(INFLUENZA_CONTIG, issue.message)

    def test_illumina_metacharacter_header_is_not_warned_about(self):
        args = {
            "data_type": "illumina",
            "reference": self._reference(f">{INFLUENZA_HEADER}\n"),
            "primer_scheme": self._write("p.bed", "Influenza\t0\t5\tp1\n"),
        }
        with self.assertNoLogs("viralconseq.validators", "WARNING"):
            validate_consensus_input_integrity(args, {"s1": [self._fastq()]})

    def test_nanopore_legacy_characters_do_not_warn(self):
        # The stock RefSeq header (comma, spaces) and a GISAID id (/ and |) were
        # rewritten by the old rule too: same contig name as before, no warning.
        header = "MN908947.3 Severe acute respiratory syndrome coronavirus 2, complete genome"
        args = {
            "data_type": "nanopore",
            "reference": self._reference(f">{header}\n"),
            "primer_scheme": self._write(
                "p.bed",
                "MN908947.3_Severe_acute_respiratory_syndrome_coronavirus_2__complete_genome\t0\t5\tp1\n",
            ),
        }
        with self.assertNoLogs("viralconseq.validators", "WARNING"):
            validate_consensus_input_integrity(args, {"s1": [self._fastq()]})

    def test_happy_path_illumina(self):
        args = {
            "data_type": "illumina",
            "reference": self._reference(),
            "primer_scheme": self._write("p.bed", "MN908947.3\t0\t20\tp1\n"),
        }
        # Must not raise.
        validate_consensus_input_integrity(args, {"s1": [self._fastq()]})

    def test_bad_fastq_blocks(self):
        args = {"data_type": "illumina", "reference": self._reference()}
        bad = self._write("s1.fastq", "@r1\nACGT\n+\nII\n")  # length mismatch
        with self.assertRaises(InputIntegrityError) as ctx:
            validate_consensus_input_integrity(args, {"s1": [bad]})
        self.assertEqual(ctx.exception.code, "input_integrity_error")
        self.assertTrue(ctx.exception.issues)
        self.assertIn("issues", ctx.exception.to_dict())

    def test_bed_chrom_mismatch_blocks(self):
        args = {
            "data_type": "illumina",
            "reference": self._reference(),
            "primer_scheme": self._write("p.bed", "WRONG\t0\t20\tp1\n"),
        }
        with self.assertRaises(InputIntegrityError):
            validate_consensus_input_integrity(args, {"s1": [self._fastq()]})

    def test_nanopore_sanitized_header_mismatch_blocks(self):
        args = {
            "data_type": "nanopore",
            "reference": self._reference(">MN908947.3 Severe acute\n"),
            "primer_scheme": self._write("p.bed", "MN908947.3\t0\t20\tp1\n"),
        }
        with self.assertRaises(InputIntegrityError) as ctx:
            validate_consensus_input_integrity(args, {"s1": [self._fastq()]})
        self.assertEqual(ctx.exception.issues[0].code, "bed_chrom_sanitized_mismatch")

    def test_skip_flag_bypasses_all_checks(self):
        args = {
            "data_type": "illumina",
            "reference": self._reference(),
            "skip_input_validation": True,
        }
        bad = self._write("s1.fastq", "totally not a fastq\n")
        # Must not raise despite the broken FASTQ.
        validate_consensus_input_integrity(args, {"s1": [bad]})

    def test_error_issues_carry_file_path(self):
        # Structured payload must map each issue to its file (F10).
        args = {"data_type": "illumina", "reference": self._reference()}
        bad = self._write("s1.fastq", "@r1\nACGT\n+\nII\n")
        with self.assertRaises(InputIntegrityError) as ctx:
            validate_consensus_input_integrity(args, {"s1": [bad]})
        issue = ctx.exception.issues[0]
        self.assertEqual(issue.path, bad)
        self.assertEqual(issue.kind, "fastq")
        self.assertEqual(ctx.exception.to_dict()["issues"][0]["path"], bad)

    def test_na_placeholder_in_segmented_dict_is_ignored(self):
        # A segmented reference dict with an 'NA' value must not try to open a
        # file named 'NA' (F3).
        ref_s = self._write("S.fasta", ">SEG_S\nACGTACGT\n")
        args = {
            "data_type": "illumina",
            "reference": {"S": ref_s, "L": "NA"},
            "primer_scheme": self._write("p.bed", "SEG_S\t0\t5\tp1\n"),
        }
        validate_consensus_input_integrity(args, {"s1": [self._fastq()]})

    def test_segmented_reference_union_of_contigs(self):
        ref_s = self._write("S.fasta", ">SEG_S\nACGTACGT\n")
        ref_l = self._write("L.fasta", ">SEG_L\nTTTTAAAA\n")
        args = {
            "data_type": "illumina",
            "reference": {"S": ref_s, "L": ref_l},
            "primer_scheme": self._write("p.bed", "SEG_S\t0\t5\tp1\nSEG_L\t0\t5\tp2\n"),
        }
        validate_consensus_input_integrity(args, {"s1": [self._fastq()]})

    def test_truncated_file_blocks_through_the_pool(self):
        # Three files over four workers: the pool path, with the broken one in
        # the middle. The issue names the right path and the old code.
        args = {"data_type": "illumina", "reference": self._reference(), "threads_total": 4}
        ok1 = self._write("s1_R1.fastq", "@r1\nACGT\n+\nIIII\n")
        bad = self._write("s1_R2.fastq", "@r1\nACGT\n+\nIIII\n@r2\nACGT\n")
        ok2 = self._write("s2.fastq", "@r1\nACGT\n+\nIIII\n")
        with self.assertRaises(InputIntegrityError) as ctx:
            validate_consensus_input_integrity(args, {"s1": [ok1, bad], "s2": [ok2]})
        issues = ctx.exception.issues
        self.assertEqual([(i.path, i.code) for i in issues], [(bad, "fastq_truncated")])
        self.assertIn(bad, str(ctx.exception))

    def test_bad_fastq_through_the_pool_matches_serial(self):
        args = {"data_type": "illumina", "reference": self._reference()}
        bad = self._write("s1.fastq", "@r1\nACGT\n+\nII\n")
        ok = self._write("s2.fastq", "@r1\nACGT\n+\nIIII\n")
        samples = {"s1": [bad], "s2": [ok]}
        outcomes = []
        for threads in (1, 2):
            with self.assertRaises(InputIntegrityError) as ctx:
                validate_consensus_input_integrity({**args, "threads_total": threads}, samples)
            outcomes.append([(i.path, i.code) for i in ctx.exception.issues])
        self.assertEqual(outcomes[0], outcomes[1])
        self.assertEqual(outcomes[0], [(bad, "fastq_length_mismatch")])


class TestOrchestrationPoolWiring(_TmpBase):
    """How validate_consensus_input_integrity sizes the FASTQ pool and reports
    a missing seqkit. The detector is patched, so no binary is needed."""

    def _reference(self):
        return self._write("ref.fasta", ">MN908947.3\nACGTACGTAC\n")

    def _fastq(self, name="s1.fastq"):
        return self._write(name, "@r1\nACGT\n+\nIIII\n")

    def test_single_file_never_builds_a_pool(self):
        args = {"data_type": "nanopore", "reference": self._reference(), "threads_total": 8}
        with mock.patch("viralconseq.integrity.validate_fastq", side_effect=_report):
            with mock.patch("viralconseq.integrity.concurrent.futures.ProcessPoolExecutor") as pool:
                validate_consensus_input_integrity(args, {"s1": [self._fastq()]})
        pool.assert_not_called()

    def test_threads_total_one_never_builds_a_pool(self):
        args = {"data_type": "illumina", "reference": self._reference(), "threads_total": 1}
        files = [self._fastq("s1_R1.fastq"), self._fastq("s1_R2.fastq")]
        with mock.patch("viralconseq.integrity.validate_fastq", side_effect=_report):
            with mock.patch("viralconseq.integrity.concurrent.futures.ProcessPoolExecutor") as pool:
                validate_consensus_input_integrity(args, {"s1": files})
        pool.assert_not_called()

    def test_absent_threads_total_means_serial(self):
        args = {"data_type": "illumina", "reference": self._reference()}
        files = [self._fastq("s1_R1.fastq"), self._fastq("s1_R2.fastq")]
        with mock.patch("viralconseq.integrity.validate_fastq", side_effect=_report):
            with mock.patch("viralconseq.integrity.concurrent.futures.ProcessPoolExecutor") as pool:
                validate_consensus_input_integrity(args, {"s1": files})
        pool.assert_not_called()

    def test_workers_follow_threads_total_capped_at_file_count(self):
        args = {"data_type": "illumina", "reference": self._reference(), "threads_total": 16}
        files = [self._fastq("s1_R1.fastq"), self._fastq("s1_R2.fastq"), self._fastq("s2.fastq")]
        with mock.patch("viralconseq.integrity.validate_fastq_files", return_value=[]) as vff:
            validate_consensus_input_integrity(args, {"s1": files[:2], "s2": files[2:]})
        vff.assert_called_once_with(files, workers=3)

    def test_missing_seqkit_is_a_validation_error(self):
        args = {"data_type": "illumina", "reference": self._reference()}
        with mock.patch("viralconseq.seqkit.shutil.which", return_value=None):
            with self.assertRaises(ValidationError) as ctx:
                validate_consensus_input_integrity(args, {"s1": [self._fastq()]})
        self.assertNotIsInstance(ctx.exception, InputIntegrityError)
        self.assertIn("seqkit", str(ctx.exception))
        self.assertIn("--skip-input-validation", str(ctx.exception))

    def test_missing_seqkit_through_the_pool_is_a_validation_error(self):
        # Two files and two workers: the pool path. The worker's error is the
        # same ValidationError the serial path gives, naming the dependency.
        args = {"data_type": "illumina", "reference": self._reference(), "threads_total": 2}
        files = [self._fastq("s1_R1.fastq"), self._fastq("s1_R2.fastq")]
        with mock.patch("viralconseq.integrity.concurrent.futures.ProcessPoolExecutor") as pool:
            pool.return_value.__enter__.return_value.map.side_effect = seqkit.SeqkitMissingError()
            with self.assertRaises(ValidationError) as ctx:
                validate_consensus_input_integrity(args, {"s1": files})
        pool.assert_called_once_with(max_workers=2)
        self.assertNotIsInstance(ctx.exception, InputIntegrityError)
        self.assertIn("seqkit", str(ctx.exception))

    def test_missing_seqkit_is_not_needed_under_skip_flag(self):
        args = {
            "data_type": "illumina",
            "reference": self._reference(),
            "skip_input_validation": True,
        }
        with mock.patch("viralconseq.seqkit.shutil.which", return_value=None):
            validate_consensus_input_integrity(args, {"s1": [self._fastq()]})


if __name__ == "__main__":
    unittest.main()
