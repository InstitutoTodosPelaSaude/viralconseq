"""Tests for scripts/python/filter_variants.py driven as the workflow runs it."""

import contextlib
import gzip
import io
import os
import tempfile
import unittest

from viralconseq.scripts.python import filter_variants

HEADER = [
    "##fileformat=VCFv4.2",
    "##FILTER=<ID=PASS,Description='All filters passed'>",
    "##FILTER=<ID=LowQual,Description='Low quality variant'>",
    "##FILTER=<ID=RefCall,Description='Reference call'>",
    "##INFO=<ID=P,Number=0,Type=Flag,Description='Pileup'>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>",
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description='Genotype Quality'>",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description='Read Depth'>",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description='Allelic depths'>",
    "##FORMAT=<ID=AF,Number=1,Type=Float,Description='Allele frequency'>",
    "##contig=<ID=ref,length=1000>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-a",
]

DEFAULT_ARGS = [
    "--variant-quality",
    "10",
    "--frameshift-variant-quality",
    "50",
    "--variant-depth",
    "10",
    "--af-threshold",
    "0.6",
]


def _rec(pos, ref, alt, qual, ad="2,18", filt="PASS", gt="1/1", fmt="GT:GQ:DP:AD:AF", extra=None):
    """One Clair3-style data line; ``ad`` None drops the AD field value."""
    total = sum(int(x) for x in ad.split(",")) if ad and ad != "." else 20
    sample = {"GT": gt, "GQ": "20", "DP": str(total), "AD": ad, "AF": "0.9"}
    if extra:
        sample.update(extra)
    keys = fmt.split(":")
    return "\t".join(
        [
            "ref",
            str(pos),
            ".",
            ref,
            alt,
            str(qual),
            filt,
            "P",
            fmt,
            ":".join(sample[k] for k in keys),
        ]
    )


def _write_vcf(path, lines, gz=False):
    text = "\n".join(HEADER + list(lines)) + "\n"
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        with open(path, "w") as fh:
            fh.write(text)
    return path


def _read_tsv(path):
    with open(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh]
    header, body = rows[0], rows[1:]
    return header, [dict(zip(header, row)) for row in body]


def _read_bed(path):
    with open(path) as fh:
        return [tuple(line.rstrip("\n").split("\t")) for line in fh if line.strip()]


def _read_vcf_records(path):
    with open(path) as fh:
        return [line.rstrip("\n") for line in fh if not line.startswith("#")]


class Test_FilterVariants(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.vcf_out = os.path.join(self.tmp, "s.pass.vcf")
        self.bed_out = os.path.join(self.tmp, "s.masked_variants.bed")
        self.tsv_out = os.path.join(self.tmp, "s.variant_filter.tsv")

    def tearDown(self):
        self._tmp.cleanup()

    def _run(self, lines, gz=False, args=None, name="s.norm.vcf"):
        path = _write_vcf(os.path.join(self.tmp, name), lines, gz=gz)
        argv = ["--vcf", path] + (args or DEFAULT_ARGS)
        argv += ["--output-vcf", self.vcf_out, "--output-bed", self.bed_out]
        argv += ["--output-tsv", self.tsv_out]
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            rc = filter_variants.main(argv)
        self.assertEqual(rc, 0)
        _, rows = _read_tsv(self.tsv_out)
        return rows, _read_bed(self.bed_out), _read_vcf_records(self.vcf_out), stderr.getvalue()

    def test_snv_passes_at_quality_bar_and_masks_below(self):
        rows, bed, vcf, _ = self._run([_rec(100, "A", "G", 10.0), _rec(200, "C", "T", 9.99)])
        self.assertEqual(
            [(r["verdict"], r["reason"]) for r in rows], [("pass", ""), ("mask", "quality")]
        )
        self.assertEqual(len(vcf), 1)
        self.assertTrue(vcf[0].startswith("ref\t100\t"))
        self.assertEqual(bed, [("ref", "199", "200")])

    def test_in_frame_deletion_passes(self):
        rows, bed, vcf, _ = self._run([_rec(100, "ACGT", "A", 20)])
        self.assertEqual(rows[0]["verdict"], "pass")
        self.assertEqual(rows[0]["frameshift"], "0")
        self.assertEqual(len(vcf), 1)
        self.assertEqual(bed, [])

    def test_frameshift_deletion_masked_with_full_ref_span(self):
        rows, bed, vcf, _ = self._run([_rec(100, "AC", "A", 30)])
        self.assertEqual(rows[0]["frameshift"], "1")
        self.assertEqual((rows[0]["verdict"], rows[0]["reason"]), ("mask", "frameshift_quality"))
        self.assertEqual(bed, [("ref", "99", "101")])
        self.assertEqual(vcf, [])

    def test_frameshift_deletion_passes_at_frameshift_bar(self):
        rows, bed, vcf, _ = self._run([_rec(100, "AC", "A", 50)])
        self.assertEqual(rows[0]["verdict"], "pass")
        self.assertEqual(len(vcf), 1)
        self.assertEqual(bed, [])

    def test_insertion_masks_anchor_base_only(self):
        rows, bed, _, _ = self._run([_rec(100, "A", "AT", 5)])
        self.assertEqual((rows[0]["verdict"], rows[0]["reason"]), ("mask", "frameshift_quality"))
        self.assertEqual(bed, [("ref", "99", "100")])

    def test_allele_ratio_failure_discarded(self):
        rows, bed, vcf, _ = self._run([_rec(100, "A", "G", 40, ad="10,5")])
        self.assertEqual((rows[0]["verdict"], rows[0]["reason"]), ("discard", "allele_ratio"))
        self.assertEqual(rows[0]["ratio"], "0.3333")
        self.assertEqual(bed, [])
        self.assertEqual(vcf, [])

    def test_allele_ratio_at_threshold_is_not_discarded(self):
        # ratio 0.6 >= 0.6 clears rule 3; the record then fails only the ALT floor.
        rows, _, _, _ = self._run([_rec(100, "A", "G", 40, ad="4,6")])
        self.assertNotEqual(rows[0]["reason"], "allele_ratio")
        self.assertEqual((rows[0]["verdict"], rows[0]["reason"]), ("mask", "alt_depth"))
        rows, _, _, _ = self._run([_rec(100, "A", "G", 40, ad="40,60")])
        self.assertEqual(rows[0]["verdict"], "pass")

    def test_zero_allele_reads_masked_not_discarded(self):
        rows, bed, _, _ = self._run([_rec(100, "A", "G", 40, ad="0,0")])
        self.assertEqual((rows[0]["verdict"], rows[0]["reason"]), ("mask", "no_allele_reads"))
        self.assertEqual(rows[0]["ratio"], "0.0000")
        self.assertEqual(bed, [("ref", "99", "100")])

    def test_alt_depth_failure_masked(self):
        rows, bed, _, _ = self._run([_rec(100, "A", "G", 40, ad="2,8")])
        self.assertEqual((rows[0]["verdict"], rows[0]["reason"]), ("mask", "alt_depth"))
        self.assertEqual(bed, [("ref", "99", "100")])

    def test_refcall_discarded(self):
        rows, bed, vcf, _ = self._run([_rec(100, "A", "G", 40, filt="RefCall;LowQual")])
        self.assertEqual((rows[0]["verdict"], rows[0]["reason"]), ("discard", "refcall"))
        self.assertEqual(rows[0]["filter"], "RefCall;LowQual")
        self.assertEqual(bed, [])
        self.assertEqual(vcf, [])

    def test_lowqual_stamp_ignored(self):
        rows, _, vcf, _ = self._run([_rec(100, "A", "G", 40, filt="LowQual")])
        self.assertEqual(rows[0]["verdict"], "pass")
        self.assertEqual(len(vcf), 1)

    def test_homozygous_reference_genotype_discarded(self):
        rows, bed, _, _ = self._run([_rec(100, "A", "G", 40, gt="0/0")])
        self.assertEqual((rows[0]["verdict"], rows[0]["reason"]), ("discard", "refcall"))
        self.assertEqual(bed, [])

    def test_missing_ad_masked(self):
        rows, bed, _, _ = self._run(
            [
                _rec(100, "A", "G", 40, ad="."),
                _rec(200, "A", "G", 40, fmt="GT:GQ:DP:AF"),
            ]
        )
        self.assertEqual([(r["verdict"], r["reason"]) for r in rows], [("mask", "missing_ad")] * 2)
        self.assertEqual([(r["ad_ref"], r["ad_alt"], r["ratio"]) for r in rows], [("NA",) * 3] * 2)
        self.assertEqual(bed, [("ref", "99", "100"), ("ref", "199", "200")])

    def test_missing_qual_masked(self):
        rows, bed, _, _ = self._run([_rec(100, "A", "G", ".")])
        self.assertEqual((rows[0]["verdict"], rows[0]["reason"]), ("mask", "quality"))
        self.assertEqual(rows[0]["qual"], ".")
        self.assertEqual(bed, [("ref", "99", "100")])

    def test_same_pos_pass_and_fail_suppresses_fail(self):
        rows, bed, vcf, _ = self._run([_rec(100, "A", "G", 40), _rec(100, "A", "T", 5, ad="15,5")])
        self.assertEqual(rows[0]["verdict"], "pass")
        self.assertEqual(
            (rows[1]["verdict"], rows[1]["reason"]), ("suppress", "allele_ratio;pass_at_site")
        )
        self.assertEqual(bed, [])
        self.assertEqual(len(vcf), 1)
        self.assertIn("\tG\t", vcf[0])

    def test_same_pos_pass_and_mask_suppresses_mask(self):
        rows, bed, _, _ = self._run([_rec(100, "A", "G", 40), _rec(100, "A", "T", 5, ad="2,18")])
        self.assertEqual(
            (rows[1]["verdict"], rows[1]["reason"]), ("suppress", "quality;pass_at_site")
        )
        self.assertEqual(bed, [])

    def test_two_passes_at_one_pos_keep_higher_qual(self):
        rows, _, vcf, _ = self._run([_rec(100, "A", "G", 20), _rec(100, "A", "T", 30)])
        self.assertEqual(
            (rows[0]["verdict"], rows[0]["reason"]), ("suppress", "multiple_pass_at_site")
        )
        self.assertEqual(rows[1]["verdict"], "pass")
        self.assertEqual(len(vcf), 1)
        self.assertIn("\tT\t", vcf[0])

    def test_two_passes_equal_qual_keep_higher_alt_depth(self):
        rows, _, _, _ = self._run(
            [_rec(100, "A", "G", 20, ad="2,18"), _rec(100, "A", "T", 20, ad="1,19")]
        )
        self.assertEqual([r["verdict"] for r in rows], ["suppress", "pass"])

    def test_mask_trimmed_around_passing_snv(self):
        # masked 5-bp deletion covers [99, 104); the passing SNV at 102 removes [101, 102)
        rows, bed, vcf, _ = self._run([_rec(100, "ACGTA", "A", 5), _rec(102, "G", "T", 40)])
        self.assertEqual(
            (rows[0]["verdict"], rows[0]["reason"]), ("mask", "frameshift_quality;mask_trimmed")
        )
        self.assertEqual(rows[1]["verdict"], "pass")
        self.assertEqual(bed, [("ref", "99", "101"), ("ref", "102", "104")])
        self.assertEqual(len(vcf), 1)

    def test_mask_inside_passing_deletion_span_suppressed(self):
        rows, bed, vcf, _ = self._run([_rec(100, "ACGT", "A", 40), _rec(102, "G", "T", 5)])
        self.assertEqual(rows[0]["verdict"], "pass")
        self.assertEqual(
            (rows[1]["verdict"], rows[1]["reason"]), ("suppress", "quality;inside_pass_span")
        )
        self.assertEqual(bed, [])
        self.assertEqual(len(vcf), 1)

    def test_bed_merges_overlapping_intervals_and_orders_chroms_by_appearance(self):
        lines = [
            _rec(100, "ACGT", "A", 5),  # [99, 103)
            _rec(102, "GTAC", "G", 5),  # [101, 105)
            _rec(300, "A", "G", 5),
        ]
        lines.append(lines[2].replace("ref\t300", "seg2\t50"))
        lines.insert(0, lines[2].replace("ref\t300", "seg2\t10"))
        rows, bed, _, _ = self._run(lines)
        self.assertEqual(
            bed,
            [
                ("seg2", "9", "10"),
                ("seg2", "49", "50"),
                ("ref", "99", "105"),
                ("ref", "299", "300"),
            ],
        )

    def test_header_only_input(self):
        rows, bed, vcf, err = self._run([])
        self.assertEqual(rows, [])
        self.assertEqual(bed, [])
        self.assertEqual(vcf, [])
        with open(self.vcf_out) as fh:
            self.assertEqual(fh.read(), "\n".join(HEADER) + "\n")
        with open(self.tsv_out) as fh:
            self.assertEqual(fh.read(), "\t".join(filter_variants.TSV_COLUMNS) + "\n")
        self.assertEqual(
            err.strip(), "filter_variants: 0 records, 0 pass, 0 mask, 0 discard, 0 suppress"
        )

    def test_gzip_input_identical_to_plain(self):
        lines = [
            _rec(100, "A", "G", 40),
            _rec(200, "AC", "A", 5),
            _rec(300, "A", "G", 40, ad="10,5"),
        ]
        plain = self._run(lines)
        with open(self.vcf_out) as fh:
            plain_vcf = fh.read()
        gz = self._run(lines, gz=True, name="s.norm.vcf.gz")
        with open(self.vcf_out) as fh:
            gz_vcf = fh.read()
        self.assertEqual(plain, gz)
        self.assertEqual(plain_vcf, gz_vcf)
        self.assertEqual(len(plain[2]), 1)
        self.assertEqual(plain[1], [("ref", "199", "201")])

    def test_multi_alt_record_rejected(self):
        with self.assertRaises(SystemExit) as ctx:
            self._run([_rec(100, "A", "G,T", 40)])
        self.assertIn("ref:100", str(ctx.exception))

    def test_two_sample_columns_rejected(self):
        line = _rec(100, "A", "G", 40) + "\t" + "1/1:20:20:2,18:0.9"
        with self.assertRaises(SystemExit) as ctx:
            self._run([line])
        self.assertIn("ref:100", str(ctx.exception))

    def test_ad_with_three_counts_rejected(self):
        with self.assertRaises(SystemExit) as ctx:
            self._run([_rec(100, "A", "G", 40, ad="2,18,0")])
        self.assertIn("ref:100", str(ctx.exception))

    def test_vcf_output_keeps_header_and_pass_records_verbatim(self):
        lines = [
            _rec(100, "A", "G", 40, filt="LowQual"),
            _rec(200, "A", "G", 5),
            _rec(300, "C", "T", 12.5),
        ]
        self._run(lines)
        with open(self.vcf_out) as fh:
            self.assertEqual(fh.read(), "\n".join(HEADER + [lines[0], lines[2]]) + "\n")

    def test_tsv_columns_and_stderr_summary(self):
        lines = [
            _rec(100, "A", "G", 40),  # pass
            _rec(200, "A", "G", 5),  # mask quality
            _rec(300, "A", "G", 40, ad="10,5"),  # discard ratio
            _rec(300, "A", "T", 40),  # pass -> discard above suppressed
            _rec(400, "A", "G", 40, filt="RefCall"),  # discard
        ]
        _, _, _, err = self._run(lines)
        header, rows = _read_tsv(self.tsv_out)
        self.assertEqual(
            header,
            [
                "chrom",
                "pos",
                "ref",
                "alt",
                "qual",
                "filter",
                "ad_ref",
                "ad_alt",
                "ratio",
                "frameshift",
                "verdict",
                "reason",
            ],
        )
        self.assertEqual(header, filter_variants.TSV_COLUMNS)
        self.assertEqual([r["pos"] for r in rows], ["100", "200", "300", "300", "400"])
        self.assertEqual(
            rows[0],
            {
                "chrom": "ref",
                "pos": "100",
                "ref": "A",
                "alt": "G",
                "qual": "40",
                "filter": "PASS",
                "ad_ref": "2",
                "ad_alt": "18",
                "ratio": "0.9000",
                "frameshift": "0",
                "verdict": "pass",
                "reason": "",
            },
        )
        self.assertEqual(
            [r["verdict"] for r in rows], ["pass", "mask", "suppress", "pass", "discard"]
        )
        self.assertEqual(
            err.strip(), "filter_variants: 5 records, 2 pass, 1 mask, 1 discard, 1 suppress"
        )

    def test_frameshift_modulo_table(self):
        self.assertEqual((-1) % 3, 2)
        for delta in range(-4, 5):
            ref = "A" * (1 + max(0, -delta))
            alt = "A" * (1 + max(0, delta))
            record = filter_variants.parse_record(_rec(10, ref, alt, 40))
            self.assertEqual(len(alt) - len(ref), delta)
            self.assertEqual(record.frameshift, delta % 3 != 0, msg=f"delta {delta}")
        self.assertEqual(
            [
                filter_variants.parse_record(_rec(10, "A" * (1 - d), "A", 40)).frameshift
                for d in (-4, -3, -2, -1)
            ],
            [True, False, True, True],
        )


if __name__ == "__main__":
    unittest.main()
