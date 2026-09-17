"""Tests for scripts/python/collect_consensus.py driven as the workflow runs it."""

import os
import tempfile
import unittest

from viralconseq.scripts.python import collect_consensus


def _fasta(path, records):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 4):
                fh.write(seq[i : i + 4] + "\n")
    return path


def _read(path):
    return collect_consensus.read_fasta(path)


def _summary(path, rows, segmented=False):
    header = ["sample_id"] + (["segment"] if segmented else []) + ["status", "coverage_min_depth"]
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(v) for v in row) + "\n")
    return path


class Test_CollectConsensus(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.out = os.path.join(self.tmp, "consensus")

    def tearDown(self):
        self._tmp.cleanup()

    def _consensus(self, sample, records):
        return _fasta(os.path.join(self.tmp, f"{sample}.consensus.renamed.fasta"), records)

    def test_per_sample_pooled_and_filtered(self):
        a = self._consensus("sample-a", [("sample-a", "ACGTACGTNN")])
        b = self._consensus("sample-b", [("sample-b", "ACGT")])
        c = self._consensus("sample-c", [("sample-c", "NNNN")])
        summary = _summary(
            os.path.join(self.tmp, "summary.tsv"),
            [
                ["sample-a", "ok", 95.0],
                ["sample-b", "ok", 69.99],
                ["sample-c", "empty_consensus", "NA"],
            ],
        )
        rc = collect_consensus.main(
            [
                "--consensus",
                a,
                b,
                c,
                "--summary",
                summary,
                "--threshold",
                "70",
                "--samples",
                "sample-a",
                "sample-b",
                "sample-c",
                "--output-dir",
                self.out,
            ]
        )
        self.assertEqual(rc, 0)
        self.assertEqual(
            sorted(os.listdir(self.out)),
            [
                "consensus.cov70.fasta",
                "consensus.fasta",
                "sample-a.fasta",
                "sample-b.fasta",
                "sample-c.fasta",
            ],
        )
        # one-line sequences, normalised header
        with open(os.path.join(self.out, "sample-a.fasta")) as fh:
            self.assertEqual(fh.read(), ">sample-a\nACGTACGTNN\n")
        self.assertEqual(
            [h for h, _ in _read(os.path.join(self.out, "consensus.fasta"))],
            ["sample-a", "sample-b", "sample-c"],
        )
        self.assertEqual(
            [h for h, _ in _read(os.path.join(self.out, "consensus.cov70.fasta"))], ["sample-a"]
        )

    def test_multi_contig_headers_and_segment_suffix(self):
        a = self._consensus("sample-a", [("sample-a_chr1", "AC"), ("sample-a_chr2", "GT")])
        summary = _summary(
            os.path.join(self.tmp, "summary.tsv"),
            [["sample-a", "L", "ok", 80.0], ["sample-a", "S", "ok", 10.0]],
            segmented=True,
        )
        collect_consensus.main(
            [
                "--consensus",
                a,
                "--summary",
                summary,
                "--threshold",
                "70",
                "--segment",
                "L",
                "--samples",
                "sample-a",
                "--output-dir",
                self.out,
            ]
        )
        self.assertEqual(
            sorted(os.listdir(self.out)),
            ["consensus.L.cov70.fasta", "consensus.L.fasta", "sample-a.L.fasta"],
        )
        self.assertEqual(
            [h for h, _ in _read(os.path.join(self.out, "sample-a.L.fasta"))],
            ["sample-a|chr1|L", "sample-a|chr2|L"],
        )
        # the L row (80%) passes; the S row is another segment and is ignored here
        self.assertEqual(len(_read(os.path.join(self.out, "consensus.L.cov70.fasta"))), 2)

    def test_threshold_tag_formatting_and_inclusive_bound(self):
        a = self._consensus("sample-a", [("sample-a", "ACGT")])
        summary = _summary(os.path.join(self.tmp, "summary.tsv"), [["sample-a", "ok", 70.5]])
        collect_consensus.main(
            [
                "--consensus",
                a,
                "--summary",
                summary,
                "--threshold",
                "70.5",
                "--samples",
                "sample-a",
                "--output-dir",
                self.out,
            ]
        )
        self.assertIn("consensus.cov70.5.fasta", os.listdir(self.out))
        self.assertEqual(len(_read(os.path.join(self.out, "consensus.cov70.5.fasta"))), 1)

    def test_normalise_header(self):
        self.assertEqual(
            collect_consensus.normalise_header("sample-a", "sample-a", None), "sample-a"
        )
        self.assertEqual(
            collect_consensus.normalise_header("sample-a_c|1 desc", "sample-a", "S"),
            "sample-a|c_1|S",
        )
        self.assertEqual(collect_consensus.normalise_header("", "sample-a", None), "sample-a")


if __name__ == "__main__":
    unittest.main()
