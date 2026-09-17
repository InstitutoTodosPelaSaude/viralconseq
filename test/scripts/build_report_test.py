"""Tests for scripts/python/build_report.py: the data document, the optional joins and
the self-check, driven as the workflow runs it (``main(argv)``) over a fake run dir."""

import csv
import gzip
import json
import os
import tempfile
import unittest

from viralconseq.scripts.python import build_report

NA = "NA"

SUMMARY_HEADER = [
    "sample_id",
    "status",
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
    "isnv_count",
    "virus",
    "clade",
    "lineage",
    "genome_quality",
    "genome_quality_score",
    "qc_overall_status",
    "viralqc_dataset",
    "viralqc_dataset_version",
    "clair3_model",
]
SEGMENTED_HEADER = SUMMARY_HEADER[:1] + ["segment"] + SUMMARY_HEADER[1:]

VIRALQC_HEADER = [
    "seqName",
    "virus",
    "virus_species",
    "clade",
    "clade_display",
    "lineage",
    "genomeQuality",
    "genomeQualityScore",
    "qc.overallStatus",
    "coverage",
    "totalSubstitutions",
    "totalDeletions",
    "totalMissing",
    "qc.privateMutations.status",
    "qc.missingData.status",
    "qc.mixedSites.status",
    "qc.snpClusters.status",
    "qc.frameShifts.status",
    "qc.stopCodons.status",
    "dataset",
    "datasetVersion",
    "inputSequenceStatus",
]


def write_tsv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def write_text(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        handle.write(content)
    return path


def ok_row(sample, cov_min=90.0, median=100.0, mapped=20000, virus="SARS-CoV-2", segment=None):
    row = [
        sample,
        "ok",
        30000,
        28000,
        mapped,
        71.43,
        120.5,
        median,
        95.0,
        80.0,
        10.0,
        cov_min,
        20,
        29903,
        1200,
        4.01,
        3,
        virus,
        "24A",
        "KP.3",
        "B",
        19,
        "good",
        "sars-cov-2",
        "2025-01-01",
        NA,
    ]
    if segment is not None:
        row.insert(1, segment)
    return row


def status_row(sample, status="no_mapped_reads", segment=None, virus=NA):
    row = [sample, status, 30000, 28000, 12, 0.04] + [NA] * 20
    row[SUMMARY_HEADER.index("virus") - 1] = virus  # header carries no "segment" here
    if segment is not None:
        row.insert(1, segment)
    return row


def document(page):
    """The JSON the page carries, which is everything the template renders."""
    with open(page, encoding="utf-8") as handle:
        text = handle.read()
    start = text.index('type="application/json">') + len('type="application/json">')
    end = text.index("</script>", start)
    return json.loads(text[start:end].replace("<\\/", "</"))


class _RunDir(unittest.TestCase):
    """A temp directory holding one fake run plus helpers to drive ``main``."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.run_dir = os.path.join(self._tmp.name, "results_20260917")
        os.makedirs(self.run_dir)
        self.page = os.path.join(self._tmp.name, "report.html")

    def tearDown(self):
        self._tmp.cleanup()

    def spine(self, rows, header=None):
        return write_tsv(os.path.join(self.run_dir, "summary.tsv"), header or SUMMARY_HEADER, rows)

    def path(self, *parts):
        return os.path.join(self.run_dir, *parts)

    def run_report(self, *extra, data_type="illumina", min_depth=20, threshold=70.0):
        argv = [
            "--run-dir",
            self.run_dir,
            "--template",
            str(build_report.TEMPLATE_PATH),
            "--output",
            self.page,
            "--data-type",
            data_type,
            "--min-depth",
            str(min_depth),
            "--consensus-coverage-threshold",
            str(threshold),
            *extra,
        ]
        self.assertEqual(build_report.main(argv), 0)
        return self.page

    def doc(self, *extra, **kwargs):
        return document(self.run_report(*extra, **kwargs))


class Test_Spine(_RunDir):
    def test_every_sample_reaches_the_page_in_spine_order(self):
        """summary.tsv keeps a row for samples that assembled nothing; so must the report,
        or a sample lost to a gate would be invisible rather than visibly flagged."""
        self.spine([ok_row("sample-a"), status_row("sample-b"), status_row("sample-c", "weird")])
        data = self.doc()
        self.assertEqual([s["id"] for s in data["samples"]], ["sample-a", "sample-b", "sample-c"])
        self.assertEqual([s["status"] for s in data["samples"]], ["ok", "no_mapped_reads", "weird"])
        self.assertEqual([s["key"] for s in data["samples"]], ["sample-a", "sample-b", "sample-c"])
        self.assertEqual(data["samples"][1]["mapped_reads"], 12)

    def test_na_becomes_null_not_the_string(self):
        """The page tests values for null; the literal 'NA' would read as present."""
        self.spine([status_row("sample-a")])
        sample = self.doc()["samples"][0]
        for field in ("median_depth", "cov_min", "cov10", "isnv_count", "clair3_model"):
            self.assertIsNone(sample[field], field)
        self.assertIsNone(sample["consensus"])
        self.assertIsNone(sample["viralqc"])

    def test_values_are_read_by_name(self):
        self.spine([ok_row("sample-a", cov_min=88.5, median=42.0, mapped=1234)])
        sample = self.doc()["samples"][0]
        self.assertEqual(sample["total_reads"], 30000)
        self.assertEqual(sample["qc_passed_reads"], 28000)
        self.assertEqual(sample["mapped_reads"], 1234)
        self.assertEqual(sample["pct_mapped"], 71.43)
        self.assertEqual(sample["mean_depth"], 120.5)
        self.assertEqual(sample["median_depth"], 42.0)
        self.assertEqual((sample["cov10"], sample["cov100"], sample["cov1000"]), (95.0, 80.0, 10.0))
        self.assertEqual(sample["cov_min"], 88.5)
        self.assertEqual(sample["consensus"], {"length": 29903, "n": 1200, "n_pct": 4.01})
        self.assertEqual(sample["isnv_count"], 3)
        self.assertEqual(sample["viralqc"]["virus"], "SARS-CoV-2")
        self.assertEqual(sample["viralqc"]["grade"], "B")
        self.assertEqual(sample["viralqc"]["score"], 19)
        self.assertEqual(sample["viralqc"]["clade"], "24A")
        self.assertEqual(sample["viralqc"]["lineage"], "KP.3")
        self.assertEqual(sample["viralqc"]["dataset_version"], "2025-01-01")

    def test_an_unknown_summary_column_is_ignored(self):
        """summary.tsv has gained columns before and will again. Reading it by name means a
        new column is inert here rather than a crash or a silent off-by-one."""
        header = SUMMARY_HEADER[:5] + ["some_future_column"] + SUMMARY_HEADER[5:]
        row = ok_row("sample-a")
        row.insert(5, "whatever")
        self.spine([row], header=header)
        sample = self.doc()["samples"][0]
        self.assertEqual(sample["cov_min"], 90.0)
        self.assertEqual(sample["status"], "ok")
        self.assertEqual(sample["pct_mapped"], 71.43)

    def test_a_missing_optional_column_yields_null(self):
        header = [c for c in SUMMARY_HEADER if c != "isnv_count"]
        row = ok_row("sample-a")
        del row[SUMMARY_HEADER.index("isnv_count")]
        self.spine([row], header=header)
        self.assertIsNone(self.doc()["samples"][0]["isnv_count"])

    def test_only_the_spine_is_required(self):
        """A run missing consensus, VCF, depth, fastp, viralQC, versions and the manifest still
        produces a page -- with those parts absent, not faked."""
        self.spine([status_row("sample-a")])
        data = self.doc()
        self.assertEqual(data["coverage"], {})
        self.assertEqual(data["amplicons"], {})
        self.assertEqual(data["run"]["versions"], [])
        self.assertIsNone(data["run"]["version"])
        self.assertIsNone(data["run"]["manifest"])
        self.assertIsNone(data["run"]["viralqc_status"])
        self.assertEqual(data["run"]["extents"], {})
        sample = data["samples"][0]
        for block in ("consensus", "viralqc", "variants", "fastp", "gate"):
            self.assertIsNone(sample[block], block)

    def test_a_directory_without_a_summary_is_refused(self):
        with self.assertRaises(SystemExit) as ctx:
            self.run_report()
        self.assertIn("finished run directory", str(ctx.exception))

    def test_run_block_shape(self):
        self.spine([ok_row("sample-a")])
        run = self.doc(
            "--af-threshold",
            "0.75",
            "--minimum-length",
            "200",
            "--minimum-map-quality",
            "20",
            data_type="nanopore",
        )["run"]
        self.assertEqual(run["pipeline"], "viralconseq")
        self.assertEqual(run["data_type"], "nanopore")
        self.assertFalse(run["segmented"])
        self.assertEqual(run["facet_kind"], "none")
        self.assertEqual(run["facets"], [])
        self.assertEqual(
            run["params"],
            {
                "min_depth": 20,
                "consensus_coverage_threshold": 70.0,
                "primer_scheme": None,
                "af_threshold": 0.75,
                "minimum_length": 200,
                "minimum_map_quality": 20,
                "clair3_model": None,
            },
        )
        self.assertIn("generated", run)
        self.assertEqual(run["output_dir"], self.run_dir)


class Test_Label(_RunDir):
    def test_the_run_label_defaults_to_the_directory_name_and_can_be_overridden(self):
        self.spine([ok_row("sample-a")])
        self.assertEqual(self.doc()["run"]["label"], "results_20260917")
        self.assertEqual(self.doc("--label", "Run 42 / L1")["run"]["label"], "Run 42 / L1")


class Test_DepthTrace(_RunDir):
    def test_the_trace_is_binned_and_the_mask_is_derived(self):
        self.spine([ok_row("sample-a")])
        # 30 positions with depth 0..29 under --min-depth 20: bins of 1 bp (30/1000 -> 1)
        rows = "".join(f"MN908947.3\t{i + 1}\t{i}\n" for i in range(30))
        write_text(self.path("assembly", "coverage_stats", "sample-a.table_cov_basewise.txt"), rows)
        trace = self.doc()["coverage"]["sample-a"]
        self.assertEqual(trace["genome_length"], 30)
        self.assertEqual(trace["bin_bp"], 1)
        self.assertEqual(trace["bins"], [float(i) for i in range(30)])
        self.assertEqual(trace["mask"], [[1, 20]])
        self.assertEqual(trace["contigs"], [{"name": "MN908947.3", "length": 30, "offset": 0}])

    def test_bin_bp_scales_with_genome_length(self):
        self.spine([ok_row("sample-a")])
        rows = "".join(f"ref\t{i + 1}\t{50 if i % 2 else 10}\n" for i in range(2500))
        write_text(self.path("assembly", "coverage_stats", "sample-a.table_cov_basewise.txt"), rows)
        trace = self.doc(min_depth=5)["coverage"]["sample-a"]
        # ceil(2500 / 1000) = 3 bp per bin -> 834 bins, alternating 10/50 -> means round to 1 dp
        self.assertEqual(trace["bin_bp"], 3)
        self.assertEqual(len(trace["bins"]), 834)
        self.assertEqual(trace["bins"][0], round((10 + 50 + 10) / 3, 1))
        self.assertEqual(trace["mask"], [])

    def test_mask_runs_are_split_and_closed_at_the_end(self):
        self.assertEqual(
            build_report.depth_mask([0, 0, 30, 30, 5, 30, 1], 20), [[1, 2], [5, 5], [7, 7]]
        )
        self.assertEqual(build_report.depth_mask([], 20), [])

    def test_multi_contig_tables_are_concatenated_with_offsets(self):
        self.spine([ok_row("sample-a")])
        rows = "".join(f"segA\t{i + 1}\t100\n" for i in range(10))
        rows += "".join(f"segB\t{i + 1}\t{0 if i < 3 else 100}\n" for i in range(5))
        write_text(self.path("assembly", "coverage_stats", "sample-a.table_cov_basewise.txt"), rows)
        trace = self.doc()["coverage"]["sample-a"]
        self.assertEqual(trace["genome_length"], 15)
        self.assertEqual(
            trace["contigs"],
            [
                {"name": "segA", "length": 10, "offset": 0},
                {"name": "segB", "length": 5, "offset": 10},
            ],
        )
        # segB positions 1-3 are below depth: 11-13 in the concatenated coordinate
        self.assertEqual(trace["mask"], [[11, 13]])
        self.assertEqual(len(trace["bins"]), 15)

    def test_a_large_spine_halves_the_trace_resolution(self):
        self.spine([ok_row(f"sample-{i:03d}") for i in range(501)])
        rows = "".join(f"ref\t{i + 1}\t30\n" for i in range(1000))
        write_text(
            self.path("assembly", "coverage_stats", "sample-000.table_cov_basewise.txt"), rows
        )
        trace = self.doc()["coverage"]["sample-000"]
        self.assertEqual(trace["bin_bp"], 2)
        self.assertEqual(len(trace["bins"]), 500)


class Test_Segmented(_RunDir):
    def test_segments_keep_independent_keys_and_details(self):
        self.spine(
            [
                ok_row("sample-a", segment="S", cov_min=95.0),
                ok_row("sample-a", segment="L", cov_min=40.0),
            ],
            header=SEGMENTED_HEADER,
        )
        for segment, depth in (("S", 20), ("L", 80)):
            write_text(
                self.path("assembly", segment, "coverage_stats", "sample-a.table_cov_basewise.txt"),
                f"{segment}_ref\t1\t{depth}\n",
            )
            write_text(
                self.path("consensus", f"sample-a.{segment}.fasta"),
                f">sample-a|{segment}\nACGTNN\n",
            )
        data = self.doc("--segments", "S", "L")
        self.assertEqual(data["run"]["facet_kind"], "segment")
        self.assertEqual(data["run"]["facets"], ["S", "L"])
        self.assertTrue(data["run"]["segmented"])
        self.assertEqual([s["key"] for s in data["samples"]], ["sample-a|S", "sample-a|L"])
        self.assertEqual([s["segment"] for s in data["samples"]], ["S", "L"])
        self.assertEqual([s["facet"] for s in data["samples"]], ["S", "L"])
        self.assertEqual(data["coverage"]["sample-a|S"]["bins"], [20])
        self.assertEqual(data["coverage"]["sample-a|L"]["bins"], [80])

    def test_segments_are_ordered_as_given_with_unlisted_ones_appended(self):
        self.spine(
            [ok_row("sample-a", segment="L"), ok_row("sample-a", segment="M")],
            header=SEGMENTED_HEADER,
        )
        self.assertEqual(self.doc("--segments", "M", "S")["run"]["facets"], ["M", "S", "L"])

    def test_consensus_file_is_read_when_the_spine_has_no_length(self):
        row = status_row("sample-a", status="ok", segment="S")
        self.spine([row], header=SEGMENTED_HEADER)
        write_text(self.path("consensus", "sample-a.S.fasta"), ">sample-a|S\nACGTNNNNAC\n")
        sample = self.doc("--segments", "S")["samples"][0]
        self.assertEqual(sample["consensus"], {"length": 10, "n": 4, "n_pct": 40.0})


class Test_ViralQC(_RunDir):
    def viralqc_row(self, seq_name, coverage="0.98", grade="A", virus="SARS-CoV-2"):
        return [
            seq_name,
            virus,
            "Severe acute respiratory syndrome-related coronavirus",
            "24A",
            "24A (JN.1)",
            "KP.3",
            grade,
            "22",
            "good",
            coverage,
            "70",
            "2",
            "300",
            "good",
            "mediocre",
            "good",
            "good",
            NA,
            "good",
            "sars-cov-2",
            "2025-01-01",
            "ok",
        ]

    def test_results_are_joined_by_name_and_fill_the_detail(self):
        row = ok_row("sample-a")
        for column in ("virus", "clade", "lineage", "genome_quality", "genome_quality_score"):
            row[SUMMARY_HEADER.index(column)] = NA
        self.spine([row])
        write_tsv(
            self.path("qc", "viralqc", "outputs", "results.tsv"),
            VIRALQC_HEADER,
            [self.viralqc_row("sample-a")],
        )
        vq = self.doc()["samples"][0]["viralqc"]
        self.assertEqual(vq["virus"], "SARS-CoV-2")
        self.assertTrue(vq["species"].startswith("Severe acute"))
        self.assertEqual(vq["clade"], "24A")
        self.assertEqual(vq["lineage"], "KP.3")
        self.assertEqual(vq["grade"], "A")
        self.assertEqual(vq["score"], 22.0)
        self.assertEqual(vq["overall"], "good")
        self.assertEqual(vq["coverage"], 0.98)
        self.assertEqual((vq["substitutions"], vq["deletions"], vq["missing"]), (70, 2, 300))
        self.assertEqual(vq["checks"]["missingData"], "mediocre")
        # a check viralQC did not report stays absent rather than becoming "good"
        self.assertIsNone(vq["checks"]["frameShifts"])
        self.assertEqual(vq["dataset"], "sars-cov-2")
        self.assertEqual(vq["input_status"], "ok")

    def test_the_spine_wins_over_results_tsv_for_shared_fields(self):
        self.spine([ok_row("sample-a")])  # grade B in the spine
        write_tsv(
            self.path("qc", "viralqc", "outputs", "results.tsv"),
            VIRALQC_HEADER,
            [self.viralqc_row("sample-a", grade="A")],
        )
        vq = self.doc()["samples"][0]["viralqc"]
        self.assertEqual(vq["grade"], "B")
        self.assertEqual(vq["substitutions"], 70)

    def test_the_failure_placeholder_shape_is_tolerated(self):
        self.spine([status_row("sample-a", status="viralqc_failed")])
        write_tsv(
            self.path("qc", "viralqc", "outputs", "results.tsv"),
            ["seqName", "genomeQuality", "inputSequenceStatus"],
            [["sample-a", NA, "not_processed"]],
        )
        write_text(
            self.path("qc", "viralqc", "viralqc_status.txt"),
            "status\tfailed\nexit_code\t1\ninput_records\t1\nresult_rows\t0\n",
        )
        data = self.doc()
        vq = data["samples"][0]["viralqc"]
        self.assertIsNone(vq["grade"])
        self.assertIsNone(vq["virus"])
        self.assertEqual(vq["input_status"], "not_processed")
        self.assertEqual(
            data["run"]["viralqc_status"],
            {"status": "failed", "exit_code": 1, "input_records": 1, "result_rows": 0},
        )

    def test_multi_contig_rows_pick_the_best_covered_one(self):
        row = ok_row("sample-a")
        row[SUMMARY_HEADER.index("clade")] = NA
        self.spine([row])
        write_tsv(
            self.path("qc", "viralqc", "outputs", "results.tsv"),
            VIRALQC_HEADER,
            [
                self.viralqc_row("sample-a|contigA", coverage="0.10", grade="D"),
                self.viralqc_row("sample-a|contigB", coverage="0.95", grade="A"),
                self.viralqc_row("sample-a|contigC", coverage=NA, grade="C"),
            ],
        )
        vq = self.doc()["samples"][0]["viralqc"]
        self.assertEqual(vq["coverage"], 0.95)

    def test_segmented_rows_join_on_sample_pipe_segment(self):
        rows = [ok_row("sample-a", segment="S"), ok_row("sample-a", segment="L")]
        for row in rows:
            row[SEGMENTED_HEADER.index("clade")] = NA
        self.spine(rows, header=SEGMENTED_HEADER)
        write_tsv(
            self.path("qc", "viralqc", "outputs", "results.tsv"),
            VIRALQC_HEADER,
            [
                self.viralqc_row("sample-a|S", coverage="0.5"),
                self.viralqc_row("sample-a|L", coverage="0.9"),
            ],
        )
        data = self.doc("--segments", "S", "L")
        self.assertEqual([s["viralqc"]["coverage"] for s in data["samples"]], [0.5, 0.9])

    def test_a_single_virus_is_not_a_facet(self):
        self.spine([ok_row("sample-a"), ok_row("sample-b")])
        run = self.doc()["run"]
        self.assertEqual(run["facet_kind"], "none")
        self.assertEqual(run["facets"], [])
        self.assertEqual([s["facet"] for s in self.doc()["samples"]], [None, None])

    def test_two_viruses_become_facets(self):
        self.spine(
            [
                ok_row("sample-a", virus="SARS-CoV-2"),
                ok_row("sample-b", virus="Influenza A"),
                status_row("sample-c"),
            ]
        )
        data = self.doc()
        self.assertEqual(data["run"]["facet_kind"], "virus")
        self.assertEqual(data["run"]["facets"], ["SARS-CoV-2", "Influenza A"])
        self.assertEqual([s["facet"] for s in data["samples"]], ["SARS-CoV-2", "Influenza A", None])

    def test_unclassified_is_not_a_second_virus(self):
        """viralQC names an unassignable consensus "Unclassified"; an all-N
        sample must not turn a single-virus run into a two-facet one."""
        self.spine(
            [
                ok_row("sample-a", virus="SARS-CoV-2"),
                ok_row("sample-b", virus="SARS-CoV-2"),
                status_row("sample-c", virus="Unclassified"),
            ]
        )
        run = self.doc()["run"]
        self.assertEqual(run["facet_kind"], "none")
        self.assertEqual(run["facets"], [])


class Test_PerSampleArtefacts(_RunDir):
    def test_fastp_block_is_illumina_only(self):
        self.spine([ok_row("sample-a")])
        write_text(
            self.path("qc", "reports", "trim.sample-a_fastp.json"),
            json.dumps(
                {
                    "summary": {
                        "before_filtering": {"total_reads": 30000, "q30_rate": 0.91},
                        "after_filtering": {
                            "total_reads": 28000,
                            "q30_rate": 0.95,
                            "read1_mean_length": 148,
                            "read2_mean_length": 147,
                        },
                    },
                    "filtering_result": {
                        "passed_filter_reads": 28000,
                        "low_quality_reads": 1500,
                        "too_short_reads": 500,
                    },
                    "duplication": {"rate": 0.12},
                    "adapter_cutting": {"adapter_trimmed_reads": 4000},
                }
            ),
        )
        fastp = self.doc(data_type="illumina")["samples"][0]["fastp"]
        self.assertEqual(fastp["total_reads_before"], 30000)
        self.assertEqual(fastp["total_reads_after"], 28000)
        self.assertEqual(fastp["q30_before"], 0.91)
        self.assertEqual(fastp["q30_after"], 0.95)
        self.assertEqual(fastp["duplication"], 0.12)
        self.assertEqual((fastp["read1_mean_length"], fastp["read2_mean_length"]), (148, 147))
        self.assertEqual(fastp["passed_filter"], 28000)
        self.assertEqual(fastp["low_quality"], 1500)
        self.assertEqual(fastp["too_short"], 500)
        self.assertEqual(fastp["adapter_trimmed"], 4000)
        self.assertIsNone(self.doc(data_type="nanopore")["samples"][0]["fastp"])

    def test_gz_vcf_records_are_counted_per_platform_layout(self):
        self.spine([ok_row("sample-a")])
        vcf = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "ref\t100\t.\tA\tG\t50\tPASS\t.\n"
            "ref\t250\t.\tC\tT\t50\tPASS\t.\n"
            "ref\t300\t.\tG\tT\t50\tPASS\t.\n"
        )
        illumina = self.path(
            "assembly", "consensus", "final_consensus", "sample-a.consensus.vcf.gz"
        )
        os.makedirs(os.path.dirname(illumina))
        with gzip.open(illumina, "wt") as handle:
            handle.write(vcf)
        self.assertEqual(self.doc(data_type="illumina")["samples"][0]["variants"], {"records": 3})
        # the nanopore layout names the file differently, so the Illumina one is not picked up
        self.assertIsNone(self.doc(data_type="nanopore")["samples"][0]["variants"])
        with gzip.open(
            self.path("assembly", "consensus", "final_consensus", "sample-a.vcf.gz"), "wt"
        ) as h:
            h.write(vcf[: vcf.rindex("ref")])
        self.assertEqual(self.doc(data_type="nanopore")["samples"][0]["variants"], {"records": 2})

    def test_gate_and_clair3_model_files_are_joined(self):
        self.spine([status_row("sample-a")])
        write_text(
            self.path("assembly", "status", "sample-a.txt"),
            "status\tno_mapped_reads\nmapped_reads\t12\nminimum_mapped_reads\t100\n",
        )
        write_text(
            self.path("assembly", "clair3", "sample-a", "model.txt"),
            "model\tr1041_e82_400bps_sup_v500\nmodel_dir\t/models/x\n",
        )
        sample = self.doc(data_type="nanopore")["samples"][0]
        self.assertEqual(
            sample["gate"],
            {"status": "no_mapped_reads", "mapped_reads": 12, "minimum_mapped_reads": 100},
        )
        self.assertEqual(sample["clair3_model"], "r1041_e82_400bps_sup_v500")

    def test_the_spine_clair3_model_wins_over_model_txt(self):
        row = ok_row("sample-a")
        row[SUMMARY_HEADER.index("clair3_model")] = "from_spine"
        self.spine([row])
        write_text(self.path("assembly", "clair3", "sample-a", "model.txt"), "model\tfrom_file\n")
        self.assertEqual(self.doc(data_type="nanopore")["samples"][0]["clair3_model"], "from_spine")


class Test_PrimerScheme(_RunDir):
    BED = (
        "MN908947.3\t30\t54\tnCoV-2019_1_LEFT\t1\t+\n"
        "MN908947.3\t385\t410\tnCoV-2019_1_RIGHT\t1\t-\n"
        "MN908947.3\t320\t342\tnCoV-2019_2_LEFT\t2\t+\n"
        "MN908947.3\t704\t726\tnCoV-2019_2_RIGHT\t2\t-\n"
    )

    def basewise(self, length=800, depth=lambda i: 100):
        rows = "".join(f"MN908947.3\t{i + 1}\t{depth(i)}\n" for i in range(length))
        write_text(self.path("assembly", "coverage_stats", "sample-a.table_cov_basewise.txt"), rows)

    def test_extent_and_amplicon_means(self):
        self.spine([ok_row("sample-a")])
        self.basewise(depth=lambda i: 100 if i < 400 else 10)
        bed = write_text(os.path.join(self._tmp.name, "scheme.primer.bed"), self.BED)
        data = self.doc("--scheme", bed)
        self.assertEqual(
            data["run"]["extents"],
            {"*": {"start": 30, "end": 726, "genome_length": 800, "outside_bp": 30 + 74}},
        )
        self.assertEqual(data["run"]["params"]["primer_scheme"], "scheme.primer.bed")
        means = data["amplicons"]["sample-a"]
        self.assertEqual(len(means), 2)
        self.assertEqual(means[0], round((370 * 100 + 10 * 10) / 380, 1))  # 30..410
        self.assertEqual(means[1], round((80 * 100 + 326 * 10) / 406, 1))  # 320..726

    def test_amplicons_are_omitted_when_primer_names_do_not_parse(self):
        self.spine([ok_row("sample-a")])
        self.basewise()
        bed = write_text(
            os.path.join(self._tmp.name, "odd.bed"),
            "MN908947.3\t30\t54\tprimerA\n" "MN908947.3\t385\t410\tprimerB\n",
        )
        data = self.doc("--scheme", bed)
        self.assertEqual(data["run"]["extents"]["*"]["start"], 30)
        self.assertEqual(data["amplicons"], {})

    def test_nanopore_runs_sanitise_bed_chromosome_names(self):
        self.spine([ok_row("sample-a")])
        rows = "".join(f"hCoV-19_Wuhan_2019\t{i + 1}\t100\n" for i in range(100))
        rows += "".join(f"other_contig\t{i + 1}\t100\n" for i in range(50))
        write_text(self.path("assembly", "coverage_stats", "sample-a.table_cov_basewise.txt"), rows)
        bed = write_text(
            os.path.join(self._tmp.name, "scheme.bed"),
            "hCoV-19/Wuhan|2019\t10\t30\tx_1_LEFT\n" "hCoV-19/Wuhan|2019\t70\t90\tx_1_RIGHT\n",
        )
        # unsanitised, the chrom matches nothing on a two-contig reference: no extent
        self.assertEqual(self.doc("--scheme", bed, data_type="illumina")["run"]["extents"], {})
        extent = self.doc("--scheme", bed, data_type="nanopore")["run"]["extents"]["*"]
        self.assertEqual((extent["start"], extent["end"], extent["genome_length"]), (10, 90, 150))

    def test_offsets_apply_on_a_multi_contig_reference(self):
        self.spine([ok_row("sample-a")])
        rows = "".join(f"A\t{i + 1}\t100\n" for i in range(100))
        rows += "".join(f"B\t{i + 1}\t100\n" for i in range(100))
        write_text(self.path("assembly", "coverage_stats", "sample-a.table_cov_basewise.txt"), rows)
        bed = write_text(
            os.path.join(self._tmp.name, "scheme.bed"),
            "B\t10\t30\tp_1_LEFT\nB\t70\t90\tp_1_RIGHT\n",
        )
        extent = self.doc("--scheme", bed)["run"]["extents"]["*"]
        self.assertEqual((extent["start"], extent["end"]), (110, 190))

    def test_a_missing_scheme_file_is_ignored(self):
        self.spine([ok_row("sample-a")])
        self.basewise()
        data = self.doc("--scheme", os.path.join(self._tmp.name, "nope.bed"))
        self.assertEqual(data["run"]["extents"], {})
        self.assertEqual(data["run"]["params"]["primer_scheme"], "nope.bed")


class Test_SelfCheck(_RunDir):
    def test_a_mismatched_filtered_fasta_is_refused(self):
        """The filtered pooled FASTA is the pipeline's own answer to which genomes are usable.
        If the page counts a different set, one of them is wrong and neither is worth publishing."""
        self.spine([ok_row("sample-a", cov_min=90.0), ok_row("sample-b", cov_min=10.0)])
        write_text(
            self.path("consensus", "consensus.cov70.fasta"), ">sample-a\nACGT\n>sample-b\nACGT\n"
        )
        with self.assertRaises(SystemExit) as ctx:
            self.run_report()
        self.assertIn("holds 2 records", str(ctx.exception))
        self.assertIn("counts 1", str(ctx.exception))
        self.assertIn("sample-b", str(ctx.exception))

    def test_identities_not_only_counts_are_compared(self):
        self.spine([ok_row("sample-a", cov_min=90.0)])
        write_text(self.path("consensus", "consensus.cov70.fasta"), ">sample-zzz\nACGT\n")
        with self.assertRaises(SystemExit) as ctx:
            self.run_report()
        self.assertIn("missing identities", str(ctx.exception))

    def test_multi_contig_headers_reduce_to_the_sample(self):
        self.spine([ok_row("sample-a", cov_min=90.0)])
        write_text(
            self.path("consensus", "consensus.cov70.fasta"),
            ">sample-a|contigA\nACGT\n>sample-a|contigB\nACGT\n",
        )
        self.assertTrue(os.path.isfile(self.run_report()))

    def test_segmented_fastas_are_checked_per_segment(self):
        self.spine(
            [
                ok_row("sample-a", segment="S", cov_min=90.0),
                ok_row("sample-a", segment="L", cov_min=10.0),
            ],
            header=SEGMENTED_HEADER,
        )
        write_text(self.path("consensus", "consensus.S.cov70.fasta"), ">sample-a|S\nACGT\n")
        write_text(self.path("consensus", "consensus.L.cov70.fasta"), "")
        self.assertTrue(os.path.isfile(self.run_report("--segments", "S", "L")))
        write_text(self.path("consensus", "consensus.L.cov70.fasta"), ">sample-a|L\nACGT\n")
        with self.assertRaises(SystemExit) as ctx:
            self.run_report("--segments", "S", "L")
        self.assertIn("consensus.L.cov70.fasta", str(ctx.exception))

    def test_other_thresholds_are_ignored(self):
        self.spine([ok_row("sample-a", cov_min=90.0)])
        write_text(self.path("consensus", "consensus.cov80.fasta"), ">obsolete\nACGT\n")
        self.assertTrue(os.path.isfile(self.run_report()))

    def test_the_self_check_can_be_waived(self):
        self.spine([ok_row("sample-a", cov_min=90.0)])
        write_text(
            self.path("consensus", "consensus.cov70.fasta"), ">sample-a\nACGT\n>sample-b\nACGT\n"
        )
        self.assertTrue(os.path.isfile(self.run_report("--no-self-check")))


class Test_Page(_RunDir):
    def test_the_page_is_self_contained(self):
        """A run report that needs the network to render is not a record."""
        self.spine([ok_row("sample-a")])
        with open(self.run_report(), encoding="utf-8") as handle:
            page = handle.read()
        self.assertIn("<script", page)
        self.assertNotIn("/*{{DATA}}*/", page)
        for pattern in ('src="http', "src='http", 'href="http', "href='http"):
            self.assertNotIn(pattern, page)

    def test_script_closers_inside_the_payload_are_escaped(self):
        self.spine([ok_row("sample-a")])
        self.assertEqual(
            self.doc("--label", "</script><b>x</b>")["run"]["label"], "</script><b>x</b>"
        )
        with open(self.page, encoding="utf-8") as handle:
            self.assertEqual(handle.read().count("</script>"), 2)

    def test_a_template_without_the_marker_is_refused(self):
        self.spine([ok_row("sample-a")])
        template = write_text(os.path.join(self._tmp.name, "t.html"), "<html></html>")
        with self.assertRaises(SystemExit) as ctx:
            build_report.main(
                [
                    "--run-dir",
                    self.run_dir,
                    "--template",
                    template,
                    "--output",
                    self.page,
                    "--data-type",
                    "illumina",
                    "--min-depth",
                    "20",
                    "--consensus-coverage-threshold",
                    "70",
                ]
            )
        self.assertIn("marker", str(ctx.exception))

    def test_an_external_resource_in_the_template_is_refused(self):
        self.spine([ok_row("sample-a")])
        template = write_text(
            os.path.join(self._tmp.name, "t.html"),
            '<script src="https://cdn.example/x.js"></script><script>/*{{DATA}}*/</script>',
        )
        with self.assertRaises(SystemExit) as ctx:
            build_report.main(
                [
                    "--run-dir",
                    self.run_dir,
                    "--template",
                    template,
                    "--output",
                    self.page,
                    "--data-type",
                    "illumina",
                    "--min-depth",
                    "20",
                    "--consensus-coverage-threshold",
                    "70",
                ]
            )
        self.assertIn("external", str(ctx.exception))


class Test_Provenance(_RunDir):
    def test_versions_and_manifest_are_optional_and_feed_the_version(self):
        self.spine([ok_row("sample-a")])
        write_tsv(
            self.path("versions.tsv"),
            ["component", "version"],
            [["viralconseq", "0.2.0"], ["snakemake", "7.32.4"], ["viralqc_db", "/db/viralqc"]],
        )
        write_text(
            self.path("run_manifest.json"),
            json.dumps(
                {
                    "viralconseq_version": "0.1.9",
                    "created_utc": "2026-09-17T10:00:00Z",
                    "finished_utc": "2026-09-17T11:00:00Z",
                    "status": "completed",
                    "config_sha256": "abc123",
                    "data_type": "illumina",
                    "run_name": "r1",
                    "samples": {
                        "sample-a": [
                            {"path": "/data/a_R1.fastq.gz", "sha256": "f00", "size_bytes": 10},
                            {"path": "/data/a_R2.fastq.gz", "sha256": "ba7", "size_bytes": 20},
                        ]
                    },
                }
            ),
        )
        run = self.doc()["run"]
        self.assertEqual(run["version"], "0.2.0")
        self.assertEqual(
            run["versions"],
            [["viralconseq", "0.2.0"], ["snakemake", "7.32.4"], ["viralqc_db", "/db/viralqc"]],
        )
        self.assertEqual(run["manifest"]["created_utc"], "2026-09-17T10:00:00Z")
        self.assertEqual(run["manifest"]["finished_utc"], "2026-09-17T11:00:00Z")
        self.assertEqual(run["manifest"]["status"], "completed")
        self.assertEqual(run["manifest"]["config_sha256"], "abc123")
        self.assertEqual(
            run["manifest"]["inputs"],
            [
                {
                    "sample": "sample-a",
                    "path": "/data/a_R1.fastq.gz",
                    "sha256": "f00",
                    "size_bytes": 10,
                },
                {
                    "sample": "sample-a",
                    "path": "/data/a_R2.fastq.gz",
                    "sha256": "ba7",
                    "size_bytes": 20,
                },
            ],
        )

    def test_the_manifest_version_is_the_fallback(self):
        self.spine([ok_row("sample-a")])
        write_text(self.path("run_manifest.json"), json.dumps({"viralconseq_version": "0.1.9"}))
        run = self.doc()["run"]
        self.assertEqual(run["version"], "0.1.9")
        self.assertEqual(run["manifest"]["inputs"], [])
        self.assertNotIn("viralconseq_version", run["manifest"])


class Test_Parser(unittest.TestCase):
    def test_required_flags(self):
        parser = build_report.build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args(["--run-dir", "x"])
        args = parser.parse_args(
            [
                "--run-dir",
                "x",
                "--template",
                "t",
                "--output",
                "o",
                "--data-type",
                "nanopore",
                "--min-depth",
                "20",
                "--consensus-coverage-threshold",
                "70",
                "--segments",
                "S",
                "L",
                "--clair3-model",
                "r1041",
                "--no-self-check",
            ]
        )
        self.assertEqual(args.segments, ["S", "L"])
        self.assertEqual(args.clair3_model, "r1041")
        self.assertTrue(args.no_self_check)
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "--run-dir",
                    "x",
                    "--template",
                    "t",
                    "--output",
                    "o",
                    "--data-type",
                    "pacbio",
                    "--min-depth",
                    "20",
                    "--consensus-coverage-threshold",
                    "70",
                ]
            )


if __name__ == "__main__":
    unittest.main()
