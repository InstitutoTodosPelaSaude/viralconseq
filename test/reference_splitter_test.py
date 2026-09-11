"""Tests for viralconseq.reference_splitter.

The reference splitter turns a single multi-record reference FASTA into the per-segment
``{segment: path}`` mapping that the segmented consensus workflows already
consume. Splitting happens in the Python layer during validation, so the
Snakemake workflows are unchanged.
"""

import gzip
import os
import tempfile
import unittest

from viralconseq.reference_splitter import (
    count_records,
    sanitize_segment_name,
    split_multifasta,
)


class TestSanitizeSegmentName(unittest.TestCase):
    def test_first_pipe_token_is_taken(self):
        self.assertEqual(
            sanitize_segment_name("NC_007373.1|Influenza_A|H3N2|segment_1"),
            "NC_007373.1",
        )

    def test_leading_gt_is_stripped(self):
        self.assertEqual(sanitize_segment_name(">KM245534.1|Guaroa_virus"), "KM245534.1")

    def test_first_whitespace_token_is_taken(self):
        self.assertEqual(sanitize_segment_name("segS some description here"), "segS")

    def test_odd_characters_are_sanitised(self):
        # No whitespace/pipe, so the whole token is kept; / \ , ~ become _.
        self.assertEqual(sanitize_segment_name("a/b\\c,d~e"), "a_b_c_d_e")

    def test_shell_and_wildcard_hostile_chars_are_collapsed(self):
        # Anything outside [A-Za-z0-9._-] (globs, quotes, braces, ...) -> _.
        self.assertEqual(sanitize_segment_name("seg*1"), "seg_1")
        self.assertEqual(sanitize_segment_name("NC_1:2"), "NC_1_2")
        self.assertEqual(sanitize_segment_name("a b&c"), "a")  # first token only

    def test_accession_is_preserved(self):
        self.assertEqual(sanitize_segment_name("NC_007373.1"), "NC_007373.1")

    def test_empty_header_raises(self):
        with self.assertRaises(ValueError):
            sanitize_segment_name(">")

    def test_dot_and_dotdot_are_rejected(self):
        with self.assertRaises(ValueError):
            sanitize_segment_name(">.")
        with self.assertRaises(ValueError):
            sanitize_segment_name(">..")


class TestCountRecords(unittest.TestCase):
    def _write(self, text, suffix=".fasta"):
        fd, path = tempfile.mkstemp(suffix=suffix)
        os.close(fd)
        with open(path, "w") as fh:
            fh.write(text)
        self.addCleanup(os.remove, path)
        return path

    def test_single_record(self):
        path = self._write(">seq1\nACGT\nACGT\n")
        self.assertEqual(count_records(path), 1)

    def test_multi_record(self):
        path = self._write(">s1\nAC\n>s2\nGT\n>s3\nTT\n")
        self.assertEqual(count_records(path), 3)

    def test_empty_file(self):
        path = self._write("")
        self.assertEqual(count_records(path), 0)

    def test_gzip_record_count(self):
        fd, path = tempfile.mkstemp(suffix=".fasta.gz")
        os.close(fd)
        with gzip.open(path, "wt") as fh:
            fh.write(">a\nAC\n>b\nGT\n")
        self.addCleanup(os.remove, path)
        self.assertEqual(count_records(path), 2)


class TestSplitMultifasta(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(self._cleanup)

    def _cleanup(self):
        import shutil

        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write_fasta(self, text, name="ref.fasta"):
        path = os.path.join(self.tmp, name)
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def test_split_produces_one_file_per_record(self):
        src = self._write_fasta(
            ">NC_007373.1|Influenza_A|segment_1\nACGT\n"
            ">NC_007372.1|Influenza_A|segment_2\nTTTT\n"
        )
        out_dir = os.path.join(self.tmp, "out")
        mapping = split_multifasta(src, out_dir)

        self.assertEqual(list(mapping.keys()), ["NC_007373.1", "NC_007372.1"])
        for path in mapping.values():
            self.assertTrue(os.path.isabs(path))
            self.assertTrue(os.path.exists(path))

    def test_original_header_is_preserved_inside_each_file(self):
        src = self._write_fasta(">NC_007373.1|Influenza_A|segment_1\nACGT\n")
        out_dir = os.path.join(self.tmp, "out")
        mapping = split_multifasta(src, out_dir)
        with open(mapping["NC_007373.1"]) as fh:
            content = fh.read()
        self.assertTrue(content.startswith(">NC_007373.1|Influenza_A|segment_1\n"))
        self.assertIn("ACGT", content)

    def test_each_file_has_exactly_one_record(self):
        src = self._write_fasta(">s1\nAC\n>s2\nGT\n")
        out_dir = os.path.join(self.tmp, "out")
        mapping = split_multifasta(src, out_dir)
        for path in mapping.values():
            self.assertEqual(count_records(path), 1)

    def test_duplicate_keys_are_deduplicated(self):
        src = self._write_fasta(">dup|first\nAC\n>dup|second\nGT\n")
        out_dir = os.path.join(self.tmp, "out")
        mapping = split_multifasta(src, out_dir)
        self.assertEqual(list(mapping.keys()), ["dup", "dup_2"])

    def test_dedup_does_not_collide_with_a_preexisting_suffix_name(self):
        # "dup", "dup", "dup_2": the second "dup" must not clobber the real
        # "dup_2" record. All three keep distinct keys and files.
        src = self._write_fasta(">dup\nAC\n>dup\nGT\n>dup_2\nTT\n")
        out_dir = os.path.join(self.tmp, "out")
        mapping = split_multifasta(src, out_dir)
        self.assertEqual(len(mapping), 3)
        self.assertEqual(len(set(mapping.values())), 3)  # three distinct files
        for path in mapping.values():
            self.assertEqual(count_records(path), 1)

    def test_record_with_empty_sequence_raises(self):
        src = self._write_fasta(">seg1\n>seg2\nACGT\n")
        with self.assertRaises(ValueError):
            split_multifasta(src, os.path.join(self.tmp, "out"))

    def test_gzip_input_is_supported(self):
        src = os.path.join(self.tmp, "ref.fasta.gz")
        with gzip.open(src, "wt") as fh:
            fh.write(">segA\nAC\n>segB\nGT\n")
        out_dir = os.path.join(self.tmp, "out")
        mapping = split_multifasta(src, out_dir)
        self.assertEqual(set(mapping.keys()), {"segA", "segB"})

    def test_empty_fasta_raises(self):
        src = self._write_fasta("")
        with self.assertRaises(ValueError):
            split_multifasta(src, os.path.join(self.tmp, "out"))


if __name__ == "__main__":
    unittest.main()
