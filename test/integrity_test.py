"""Tests for content-level input-integrity validation (viralconseq.integrity)
and its orchestration in validators.validate_consensus_input_integrity."""

import gzip
import os
import tempfile
import unittest

from viralconseq import integrity
from viralconseq.exceptions import InputIntegrityError
from viralconseq.validators import validate_consensus_input_integrity


def _codes(report):
    return [issue.code for issue in report.errors]


def _warn_codes(report):
    return [issue.code for issue in report.warnings]


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


class TestValidateFastq(_TmpBase):
    def test_valid_fastq(self):
        p = self._write("ok.fastq", "@r1\nACGT\n+\nIIII\n@r2\nACGTN\n+\nIIIII\n")
        self.assertEqual(_codes(integrity.validate_fastq(p)), [])

    def test_valid_gzipped_fastq(self):
        p = os.path.join(self.tmp, "ok.fastq.gz")
        with gzip.open(p, "wt") as fh:
            fh.write("@r1\nACGT\n+\nIIII\n")
        self.assertEqual(_codes(integrity.validate_fastq(p)), [])

    def test_truncated_record(self):
        p = self._write("trunc.fastq", "@r1\nACGT\n+\nIIII\n@r2\nACGT\n")
        self.assertIn("fastq_truncated", _codes(integrity.validate_fastq(p)))

    def test_seq_qual_length_mismatch(self):
        p = self._write("mm.fastq", "@r1\nACGT\n+\nIII\n")
        self.assertIn("fastq_length_mismatch", _codes(integrity.validate_fastq(p)))

    def test_bad_sequence_char(self):
        p = self._write("badseq.fastq", "@r1\nACGZ\n+\nIIII\n")
        self.assertIn("fastq_bad_sequence_char", _codes(integrity.validate_fastq(p)))

    def test_bad_header_marker(self):
        # A FASTA mislabeled as FASTQ: first line starts with '>'.
        p = self._write("wrong.fastq", ">r1\nACGT\n+\nIIII\n")
        self.assertIn("fastq_bad_header", _codes(integrity.validate_fastq(p)))

    def test_empty_file(self):
        p = self._write("empty.fastq", "")
        self.assertIn("fastq_empty", _codes(integrity.validate_fastq(p)))

    def test_corrupt_gzip_magic_but_garbage(self):
        # gzip magic bytes but not a valid archive -> unreadable, not a crash.
        p = self._write("bad.fastq.gz", b"\x1f\x8b not a real gzip stream", mode="wb")
        self.assertIn("fastq_unreadable", _codes(integrity.validate_fastq(p)))

    def test_issue_cap(self):
        # Many mismatched records; issue count must be bounded by the cap.
        p = self._write("many.fastq", "@r\nAC\n+\nI\n" * 100)
        report = integrity.validate_fastq(p)
        self.assertLessEqual(len(report.issues), integrity.MAX_ISSUES_PER_FILE)

    def test_blank_lines_between_records_tolerated(self):
        # A stray blank line must not desync grouping or fake a truncation (F5).
        p = self._write("blanks.fastq", "@r1\nACGT\n+\nIIII\n\n@r2\nACGT\n+\nIIII\n")
        self.assertEqual(_codes(integrity.validate_fastq(p)), [])

    def test_truncation_reported_even_at_issue_cap(self):
        # 20 bad records + a truncated tail: fastq_truncated must still surface (F14).
        p = self._write("trunc_cap.fastq", "@r\nAC\n+\nI\n" * 25 + "@last\nACGT\n")
        codes = _codes(integrity.validate_fastq(p))
        self.assertIn("fastq_truncated", codes)

    def test_binary_file_reported_unreadable(self):
        p = self._write("bin.fastq", b"\x00\x01\x02\xff\xfe non-text bytes", mode="wb")
        self.assertIn("fastq_unreadable", _codes(integrity.validate_fastq(p)))

    def test_utf8_bom_fastq_tolerated(self):
        p = self._write("bom.fastq", "﻿@r1\nACGT\n+\nIIII\n")
        self.assertEqual(_codes(integrity.validate_fastq(p)), [])


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


class TestOrchestration(_TmpBase):
    def _reference(self, header=">MN908947.3\n"):
        return self._write("ref.fasta", header + "ACGTACGTAC\n")

    def _fastq(self):
        return self._write("s1.fastq", "@r1\nACGT\n+\nIIII\n")

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


if __name__ == "__main__":
    unittest.main()
