"""Tests for scripts/python/build_summary.py driven as the workflow runs it
(``main(argv)``) over small fake inputs."""

import os
import tempfile
import unittest

from viralconseq.scripts.python import build_summary

STATS_HEADER = [
    "sample_id",
    "total_reads",
    "qc_passed_reads",
    "mapped_reads",
    "pct_mapped",
    "mean_depth",
    "median_depth",
    "coverage_10x",
    "coverage_100x",
    "coverage_1000x",
    "coverage_min_depth",
    "min_depth",
    "consensus_length",
    "n_count",
    "n_pct",
]
VQC_HEADER = [
    "seqName",
    "virus",
    "clade",
    "lineage",
    "genomeQuality",
    "genomeQualityScore",
    "qc.overallStatus",
    "coverage",
    "dataset",
    "datasetVersion",
]


def _tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(v) for v in row) + "\n")
    return path


def _read(path, delimiter="\t"):
    with open(path) as fh:
        return [line.rstrip("\n").split(delimiter) for line in fh]


class Test_BuildSummary(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.out = os.path.join(self.tmp, "summary.tsv")

    def tearDown(self):
        self._tmp.cleanup()

    def _stats(self, sample, cov_min=95.5, length=29903, n=100, segment=None):
        header = list(STATS_HEADER)
        row = [
            sample,
            1000,
            900,
            800,
            80.0,
            150.0,
            140.0,
            99.0,
            90.0,
            10.0,
            cov_min,
            20,
            length,
            n,
            round(100 * n / length, 2) if length else "NA",
        ]
        if segment is not None:
            header.insert(1, "segment")
            row.insert(1, segment)
        name = f"{sample}{'.' + segment if segment else ''}.stats.tsv"
        return _tsv(os.path.join(self.tmp, name), header, [row])

    def _run(self, *argv):
        rc = build_summary.main(list(argv) + ["--output", self.out])
        self.assertEqual(rc, 0)
        rows = _read(self.out)
        header = rows[0]
        return header, [dict(zip(header, r)) for r in rows[1:]]

    def test_pinned_header_and_row_order(self):
        s_b = self._stats("sample-b")
        s_a = self._stats("sample-a")
        header, rows = self._run(
            "--stats",
            s_a,
            s_b,
            "--samples",
            "sample-b",
            "sample-a",
            "--min-depth",
            "20",
            "--data-type",
            "illumina",
        )
        self.assertEqual(header, [c for c in build_summary.COLUMNS if c != "segment"])
        self.assertEqual([r["sample_id"] for r in rows], ["sample-b", "sample-a"])
        self.assertEqual(rows[0]["status"], "ok")
        self.assertEqual(rows[0]["coverage_min_depth"], "95.5")
        self.assertEqual(rows[0]["clair3_model"], "NA")
        self.assertEqual(rows[0]["virus"], "NA")

    def test_every_sample_gets_a_row_even_without_stats(self):
        header, rows = self._run(
            "--stats",
            self._stats("sample-a"),
            "--samples",
            "sample-a",
            "sample-z",
            "--min-depth",
            "20",
            "--data-type",
            "nanopore",
        )
        self.assertEqual(rows[1]["sample_id"], "sample-z")
        self.assertEqual(rows[1]["status"], "missing_stats")
        self.assertEqual(rows[1]["total_reads"], "NA")

    def test_viralqc_join_best_coverage_and_status(self):
        vqc = _tsv(
            os.path.join(self.tmp, "results.tsv"),
            VQC_HEADER,
            [
                [
                    "sample-a|c1",
                    "SARS-CoV-2",
                    "B.1",
                    "",
                    "B",
                    20,
                    "good",
                    0.5,
                    "nextstrain/sars-cov-2",
                    "2026-01-06",
                ],
                [
                    "sample-a|c2",
                    "SARS-CoV-2",
                    "B.1.1.33",
                    "",
                    "A",
                    24,
                    "good",
                    0.99,
                    "nextstrain/sars-cov-2",
                    "2026-01-06",
                ],
            ],
        )
        status = os.path.join(self.tmp, "viralqc_status.txt")
        with open(status, "w") as fh:
            fh.write("status\tok\nexit_code\t0\n")
        _, rows = self._run(
            "--stats",
            self._stats("sample-a"),
            self._stats("sample-b"),
            "--samples",
            "sample-a",
            "sample-b",
            "--min-depth",
            "20",
            "--data-type",
            "illumina",
            "--viralqc-results",
            vqc,
            "--viralqc-status",
            status,
        )
        a, b = rows
        self.assertEqual(a["clade"], "B.1.1.33")  # highest coverage contig row
        self.assertEqual(a["genome_quality"], "A")
        self.assertEqual(a["genome_quality_score"], "24")
        self.assertEqual(a["lineage"], "NA")  # empty string -> NA
        self.assertEqual(a["viralqc_dataset_version"], "2026-01-06")
        self.assertEqual(a["status"], "ok")
        self.assertEqual(b["status"], "viralqc_missing")

    def test_viralqc_placeholder_shape_and_failed_status(self):
        vqc = _tsv(
            os.path.join(self.tmp, "results.tsv"),
            ["seqName", "genomeQuality", "inputSequenceStatus"],
            [["sample-a", "", "viralQC failed (exit 1)"]],
        )
        status = os.path.join(self.tmp, "viralqc_status.txt")
        with open(status, "w") as fh:
            fh.write("status\tfailed\nexit_code\t1\n")
        _, rows = self._run(
            "--stats",
            self._stats("sample-a"),
            "--samples",
            "sample-a",
            "--min-depth",
            "20",
            "--data-type",
            "illumina",
            "--viralqc-results",
            vqc,
            "--viralqc-status",
            status,
        )
        self.assertEqual(rows[0]["status"], "viralqc_failed")
        self.assertEqual(rows[0]["virus"], "NA")
        self.assertEqual(rows[0]["genome_quality"], "NA")

    def test_status_file_wins_then_empty_consensus(self):
        status_dir = os.path.join(self.tmp, "status")
        os.makedirs(status_dir)
        with open(os.path.join(status_dir, "sample-a.txt"), "w") as fh:
            fh.write("status\tno_mapped_reads\nmapped_reads\t3\n")
        with open(os.path.join(status_dir, "sample-b.txt"), "w") as fh:
            fh.write("status\tok\n")
        model_dir = os.path.join(self.tmp, "assembly", "clair3", "sample-a")
        os.makedirs(model_dir)
        with open(os.path.join(model_dir, "model.txt"), "w") as fh:
            fh.write("model\tr941_prom_hac_g360+g422\nmodel_dir\t/m\n")
        _, rows = self._run(
            "--stats",
            self._stats("sample-a", length=29903, n=29903),
            self._stats("sample-b", length=0, n=0),
            "--samples",
            "sample-a",
            "sample-b",
            "--min-depth",
            "20",
            "--data-type",
            "nanopore",
            "--status-files",
            os.path.join(status_dir, "sample-a.txt"),
            os.path.join(status_dir, "sample-b.txt"),
            "--model-files",
            os.path.join(model_dir, "model.txt"),
            "--clair3-model",
            "fallback_model",
        )
        self.assertEqual(rows[0]["status"], "no_mapped_reads")
        self.assertEqual(rows[0]["clair3_model"], "r941_prom_hac_g360+g422")
        self.assertEqual(rows[1]["status"], "empty_consensus")
        self.assertEqual(rows[1]["clair3_model"], "fallback_model")

    def test_segmented_rows_isnvs_and_legacy_csv(self):
        stats = [
            self._stats("sample-a", segment="S"),
            self._stats("sample-a", segment="L", cov_min=40.0),
        ]
        isnvs = _tsv(
            os.path.join(self.tmp, "isnvs.tsv"),
            ["sample", "segment", "number_of_isnvs"],
            [["sample-a", "S", 4], ["sample-a", "L", 1]],
        )
        legacy = os.path.join(self.tmp, "legacy.csv")
        header, rows = self._run(
            "--stats",
            *stats,
            "--samples",
            "sample-a",
            "--segments",
            "S",
            "L",
            "--min-depth",
            "20",
            "--data-type",
            "illumina",
            "--isnvs",
            isnvs,
            "--legacy-csv",
            legacy,
        )
        self.assertEqual(header[:3], ["sample_id", "segment", "status"])
        self.assertEqual([(r["segment"], r["isnv_count"]) for r in rows], [("S", "4"), ("L", "1")])
        legacy_rows = _read(legacy, ",")
        self.assertEqual(legacy_rows[0], build_summary.LEGACY_COLUMNS)
        self.assertEqual(legacy_rows[1][:2], ["sample-a", "S"])
        self.assertEqual(legacy_rows[1][-1], "0.955")  # fraction, not percent
        self.assertEqual(legacy_rows[2][-1], "0.4")
        self.assertEqual(legacy_rows[1][3], "900")  # number_of_trim_paired_reads = qc_passed_reads

    def test_legacy_csv_unsegmented_uses_total_reads_when_no_qc(self):
        header = [c for c in STATS_HEADER]
        row = ["sample-a", 10, "NA", 5, 50.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 20, 100, 100, 100.0]
        stats = _tsv(os.path.join(self.tmp, "a.stats.tsv"), header, [row])
        legacy = os.path.join(self.tmp, "legacy.csv")
        self._run(
            "--stats",
            stats,
            "--samples",
            "sample-a",
            "--min-depth",
            "20",
            "--data-type",
            "nanopore",
            "--legacy-csv",
            legacy,
        )
        legacy_rows = _read(legacy, ",")
        self.assertEqual(
            legacy_rows[0], [c for c in build_summary.LEGACY_COLUMNS if c != "segment"]
        )
        self.assertEqual(legacy_rows[1][2], "10")  # number_of_trim_paired_reads falls back to total


if __name__ == "__main__":
    unittest.main()
