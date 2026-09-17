"""Tests for scripts/python/calculate_assembly_stats.py on real temp files;
only the samtools call is patched."""

import gzip
import os
import tempfile
import unittest
from unittest.mock import patch

from viralconseq.scripts.python import calculate_assembly_stats as stats

# Two records; the separator carries the read id and a one-base read has a "+"
# quality line, both spec-legal traps for a "+"-counting approach.
FASTQ = "@r1\nACGT\n+r1\n!!!!\n@r2\nA\n+\n+\n"


def _write(path, text, gz=False):
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        with open(path, "w") as fh:
            fh.write(text)
    return path


class Test_Counts(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def test_count_reads_plain_and_gz(self):
        self.assertEqual(stats.count_reads(_write(os.path.join(self.tmp, "a.fq"), FASTQ)), 2)
        self.assertEqual(
            stats.count_reads(_write(os.path.join(self.tmp, "a.fq.gz"), FASTQ, gz=True)), 2
        )
        self.assertEqual(stats.count_reads(_write(os.path.join(self.tmp, "e.fq"), "")), 0)

    def test_count_mapped_reads_uses_samtools(self):
        with patch.object(stats.subprocess, "run") as run:
            run.return_value.stdout = "10\n"
            self.assertEqual(stats.count_mapped_reads("s.bam"), 10)
        run.assert_called_once_with(
            ["samtools", "view", "-c", "-F", "260", "s.bam"],
            check=True,
            capture_output=True,
            text=True,
        )


class Test_CoverageStats(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def test_percentages_mean_and_median(self):
        table = _write(
            os.path.join(self.tmp, "cov.txt"),
            "chr\t1\t0\nchr\t2\t10\nchr\t3\t100\nchr\t4\t1000\nchr\t5\t10000\n",
        )
        out = stats.coverage_stats(table, minimum_depth=20)
        self.assertEqual(out["mean_depth"], 2222.0)
        self.assertEqual(out["median_depth"], 100.0)
        self.assertEqual(out["coverage_10x"], 80.0)
        self.assertEqual(out["coverage_100x"], 60.0)
        self.assertEqual(out["coverage_1000x"], 40.0)
        self.assertEqual(out["coverage_min_depth"], 60.0)

    def test_multi_contig_table_counts_every_position(self):
        table = _write(os.path.join(self.tmp, "cov.txt"), "A\t1\t30\nA\t2\t30\nB\t1\t0\nB\t2\t0\n")
        out = stats.coverage_stats(table, minimum_depth=20)
        self.assertEqual(out["coverage_min_depth"], 50.0)
        self.assertEqual(out["mean_depth"], 15.0)

    def test_empty_table_is_zero_breadth_and_na_depth(self):
        table = _write(os.path.join(self.tmp, "cov.txt"), "")
        out = stats.coverage_stats(table, minimum_depth=20)
        self.assertEqual(out["coverage_min_depth"], 0.0)
        self.assertEqual(out["mean_depth"], "NA")


class Test_ConsensusStats(unittest.TestCase):
    def test_length_and_n_over_all_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = _write(os.path.join(tmp, "c.fa"), ">s|A\nACGTNN\nNNAC\n>s|B\nNNNNNNNNNN\n")
            out = stats.consensus_stats(fasta)
        self.assertEqual(out, {"consensus_length": 20, "n_count": 14, "n_pct": 70.0})

    def test_empty_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = stats.consensus_stats(_write(os.path.join(tmp, "c.fa"), ""))
        self.assertEqual(out, {"consensus_length": 0, "n_count": 0, "n_pct": "NA"})


class Test_BuildRow(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.r1 = _write(os.path.join(self.tmp, "r1.fq"), FASTQ)
        self.r2 = _write(os.path.join(self.tmp, "r2.fq"), FASTQ)
        self.q1 = _write(os.path.join(self.tmp, "q1.fq"), "@r1\nACGT\n+\n!!!!\n")
        self.q2 = _write(os.path.join(self.tmp, "q2.fq"), "@r1\nACGT\n+\n!!!!\n")
        self.cov = _write(os.path.join(self.tmp, "cov.txt"), "c\t1\t25\nc\t2\t5\n")
        self.fa = _write(os.path.join(self.tmp, "c.fa"), ">sample-a\nACNN\n")

    def tearDown(self):
        self._tmp.cleanup()

    def _row(self, qc, segment=None):
        with patch.object(stats, "count_mapped_reads", return_value=3):
            return stats.build_row(
                "sample-a", segment, [self.r1, self.r2], qc, "s.bam", self.cov, self.fa, 20
            )

    def test_illumina_counts_both_mates_and_qc_passed(self):
        row = self._row([self.q1, self.q2])
        self.assertEqual(row["total_reads"], 4)
        self.assertEqual(row["qc_passed_reads"], 2)
        self.assertEqual(row["mapped_reads"], 3)
        self.assertEqual(row["pct_mapped"], 75.0)
        self.assertEqual(row["coverage_min_depth"], 50.0)
        self.assertEqual(row["min_depth"], 20)
        self.assertEqual((row["consensus_length"], row["n_count"], row["n_pct"]), (4, 2, 50.0))
        self.assertNotIn("segment", row)

    def test_nanopore_has_na_qc_passed_and_segment_column(self):
        row = self._row(None, segment="S")
        self.assertEqual(row["qc_passed_reads"], "NA")
        self.assertEqual(row["segment"], "S")

    def test_write_tsv_has_header_in_column_order(self):
        row = self._row([self.q1, self.q2])
        out = os.path.join(self.tmp, "stats.tsv")
        stats.write_tsv(row, out)
        with open(out) as fh:
            header, values = fh.read().rstrip("\n").split("\n")
        self.assertEqual(header.split("\t"), [c for c in stats.COLUMNS if c != "segment"])
        self.assertEqual(values.split("\t")[0], "sample-a")
        self.assertEqual(len(values.split("\t")), len(header.split("\t")))


if __name__ == "__main__":
    unittest.main()
