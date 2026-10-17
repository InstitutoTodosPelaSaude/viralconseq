"""Sort the Clair3 calls of one sample into ``pass`` / ``mask`` / ``discard``
and emit the records to apply plus a BED of positions to replace with N.

Run by ``rule infer_consensus_sequence`` (``rules/consensus_nanopore.smk``)
on the ``bcftools norm -m -`` output of Clair3, where the ``viralconseq``
package is NOT installed: standard library only, ``main(argv)``::

    python filter_variants.py --vcf <s>.norm.vcf.gz \\
        --variant-quality 10 --frameshift-variant-quality 50 \\
        --variant-depth 10 --af-threshold 0.6 \\
        --output-vcf <s>.pass.vcf --output-bed <s>.masked_variants.bed \\
        --output-tsv <s>.variant_filter.tsv

The previous ``bcftools filter`` expression silently reverted every rejected
call to the reference base. Following ARTIC's ``vcf_filter``, a call that
fails is now either *discarded* (the evidence says the reference base is
right) or *masked* (the evidence is too weak to tell either way, so the
positions of its REF span become N). One verdict per record, first match
wins; ``QUAL`` is read directly and Clair3's ``LowQual`` stamp is ignored:

1. ``RefCall`` in FILTER, or a homozygous-reference/missing GT -> ``discard``
   (``refcall``).
2. AD absent -> ``mask`` (``missing_ad``); AD ``0,0`` (every read shows a
   deletion or a third base) -> ``mask`` (``no_allele_reads``).
3. ``ad_alt / (ad_ref + ad_alt) < af_threshold`` -> ``discard``
   (``allele_ratio``): the reference allele is the majority.
4. ``ad_alt < variant_depth`` -> ``mask`` (``alt_depth``).
5. ``QUAL`` (``.`` counts as 0) below the bar -> ``mask`` (``quality``, or
   ``frameshift_quality`` under the higher bar when ``len(ALT) - len(REF)``
   is not a multiple of 3).
6. Otherwise ``pass``.

``bcftools consensus`` drops a passing record whose REF span overlaps a
``--mask`` interval, so sites are reconciled afterwards: a non-pass record
at a position that also carries a pass is *suppressed* (written nowhere),
only the best of several passes at one position survives, and mask
intervals lose every position covered by a passing record's REF span.
The TSV lists every input record with its verdict and reason.
"""

import argparse
import csv
import gzip
import sys
from collections import OrderedDict
from typing import IO, Dict, List, Optional, Sequence, Tuple

NA = "NA"
PASS = "pass"
MASK = "mask"
DISCARD = "discard"
SUPPRESS = "suppress"

GZIP_MAGIC = b"\x1f\x8b"
HOM_REF_GENOTYPES = {"0", "0/0", "0|0", ".", "./."}
TSV_COLUMNS = [
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
]

Interval = Tuple[int, int]


class Record:
    """One VCF data line with the facts the filter reads and its verdict."""

    def __init__(self, line: str, fields: Sequence[str]) -> None:
        self.line = line
        self.chrom = fields[0]
        self.pos = int(fields[1])
        self.ref = fields[3]
        self.alt = fields[4]
        self.qual_text = fields[5]
        self.qual = 0.0 if fields[5] == "." else float(fields[5])
        self.filters = fields[6].split(";")
        self.genotype: Optional[str] = None
        self.ad: Optional[Tuple[int, int]] = None
        self.verdict = PASS
        self.reason = ""
        self.intervals: List[Interval] = []

    @property
    def site(self) -> str:
        return f"{self.chrom}:{self.pos}"

    @property
    def ad_alt(self) -> int:
        return self.ad[1] if self.ad is not None else 0

    @property
    def ratio(self) -> Optional[float]:
        """ALT over REF plus ALT reads; ``None`` without AD, 0.0 for AD ``0,0``."""
        if self.ad is None:
            return None
        total = self.ad[0] + self.ad[1]
        return self.ad[1] / total if total else 0.0

    @property
    def frameshift(self) -> bool:
        return (len(self.alt) - len(self.ref)) % 3 != 0

    @property
    def span(self) -> Interval:
        """0-based half-open REF span."""
        return (self.pos - 1, self.pos - 1 + len(self.ref))


def open_text(path: str) -> IO[str]:
    """Open plain or gzip text, deciding by the magic bytes rather than the suffix."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == GZIP_MAGIC:
        return gzip.open(path, "rt")
    return open(path)


def fail(record_site: str, message: str) -> "SystemExit":
    return SystemExit(f"filter_variants: {message} at {record_site}")


def parse_record(line: str) -> Record:
    fields = line.rstrip("\r\n").split("\t")
    site = f"{fields[0]}:{fields[1]}" if len(fields) > 1 else line.strip()
    if len(fields) < 10:
        raise fail(site, "expected 8 fixed columns, FORMAT and one sample column")
    if len(fields) > 10:
        raise fail(site, f"expected a single sample column, found {len(fields) - 9}")
    if "," in fields[4]:
        raise fail(site, f"multi-allelic record (ALT {fields[4]}); run bcftools norm -m - first")
    record = Record(line, fields)
    keys = fields[8].split(":")
    values = fields[9].split(":")
    sample: Dict[str, str] = dict(zip(keys, values))
    record.genotype = sample.get("GT")
    ad_text = sample.get("AD")
    if ad_text is not None and ad_text != ".":
        parts = ad_text.split(",")
        if len(parts) != 2:
            raise fail(site, f"AD must hold two counts (ref,alt), found {len(parts)} ({ad_text})")
        try:
            record.ad = (int(parts[0]), int(parts[1]))
        except ValueError:
            raise fail(site, f"AD must hold two integers, found {ad_text}") from None
    return record


def read_vcf(path: str) -> Tuple[List[str], List[Record]]:
    headers: List[str] = []
    records: List[Record] = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#"):
                headers.append(line.rstrip("\r\n"))
            else:
                records.append(parse_record(line))
    return headers, records


def judge(
    record: Record,
    variant_quality: float,
    frameshift_quality: float,
    variant_depth: int,
    af_threshold: float,
) -> Tuple[str, str]:
    """Return ``(verdict, reason)`` for one record, first matching rule wins."""
    if "RefCall" in record.filters or record.genotype in HOM_REF_GENOTYPES:
        return DISCARD, "refcall"
    if record.ad is None:
        return MASK, "missing_ad"
    if record.ad[0] + record.ad[1] == 0:
        # Every read shows a deletion or a third base: no evidence for the
        # reference either, so a confident reference base would be a guess.
        return MASK, "no_allele_reads"
    if record.ratio < af_threshold:
        return DISCARD, "allele_ratio"
    if record.ad_alt < variant_depth:
        return MASK, "alt_depth"
    if record.frameshift:
        if record.qual < frameshift_quality:
            return MASK, "frameshift_quality"
    elif record.qual < variant_quality:
        return MASK, "quality"
    return PASS, ""


def reconcile_sites(records: List[Record]) -> None:
    """One winner per site: suppress the losers next to a pass."""
    by_site: "OrderedDict[Tuple[str, int], List[Record]]" = OrderedDict()
    for record in records:
        by_site.setdefault((record.chrom, record.pos), []).append(record)
    for group in by_site.values():
        passing = [r for r in group if r.verdict == PASS]
        if not passing:
            continue
        for record in group:
            if record.verdict != PASS:
                record.verdict = SUPPRESS
                record.reason = f"{record.reason};pass_at_site"
        if len(passing) > 1:
            best = max(passing, key=lambda r: (r.qual, r.ad_alt))
            for record in passing:
                if record is not best:
                    record.verdict = SUPPRESS
                    record.reason = "multiple_pass_at_site"


def subtract(interval: Interval, blocked: Sequence[Interval]) -> List[Interval]:
    """``interval`` minus every ``blocked`` span, as sorted disjoint pieces."""
    pieces = [interval]
    for b_start, b_end in blocked:
        kept: List[Interval] = []
        for start, end in pieces:
            if b_end <= start or b_start >= end:
                kept.append((start, end))
                continue
            if start < b_start:
                kept.append((start, b_start))
            if b_end < end:
                kept.append((b_end, end))
        pieces = kept
    return pieces


def build_mask_intervals(records: List[Record]) -> None:
    """Assign each ``mask`` record its REF span minus every passing REF span."""
    passing_spans: Dict[str, List[Interval]] = {}
    for record in records:
        if record.verdict == PASS:
            passing_spans.setdefault(record.chrom, []).append(record.span)
    for record in records:
        if record.verdict != MASK:
            continue
        pieces = subtract(record.span, passing_spans.get(record.chrom, []))
        if not pieces:
            record.verdict = SUPPRESS
            record.reason = f"{record.reason};inside_pass_span"
            continue
        if pieces != [record.span]:
            record.reason = f"{record.reason};mask_trimmed"
        record.intervals = pieces


def merged_bed(records: List[Record]) -> List[Tuple[str, int, int]]:
    """Sorted (chrom in order of first appearance, start) merged intervals."""
    by_chrom: "OrderedDict[str, List[Interval]]" = OrderedDict()
    for record in records:
        by_chrom.setdefault(record.chrom, [])
        if record.verdict == MASK:
            by_chrom[record.chrom].extend(record.intervals)
    out: List[Tuple[str, int, int]] = []
    for chrom, intervals in by_chrom.items():
        current: Optional[Interval] = None
        for start, end in sorted(intervals):
            if current is not None and start <= current[1]:
                current = (current[0], max(current[1], end))
                continue
            if current is not None:
                out.append((chrom, current[0], current[1]))
            current = (start, end)
        if current is not None:
            out.append((chrom, current[0], current[1]))
    return out


def tsv_row(record: Record) -> List[str]:
    ratio = record.ratio
    return [
        record.chrom,
        str(record.pos),
        record.ref,
        record.alt,
        record.qual_text,
        ";".join(record.filters),
        str(record.ad[0]) if record.ad is not None else NA,
        str(record.ad[1]) if record.ad is not None else NA,
        f"{ratio:.4f}" if ratio is not None else NA,
        "1" if record.frameshift else "0",
        record.verdict,
        record.reason,
    ]


def write_outputs(
    headers: Sequence[str],
    records: List[Record],
    output_vcf: str,
    output_bed: str,
    output_tsv: str,
) -> None:
    with open(output_vcf, "w") as vcf:
        for line in headers:
            vcf.write(line + "\n")
        for record in records:
            if record.verdict == PASS:
                vcf.write(record.line.rstrip("\r\n") + "\n")
    with open(output_bed, "w") as bed:
        for chrom, start, end in merged_bed(records):
            bed.write(f"{chrom}\t{start}\t{end}\n")
    with open(output_tsv, "w", newline="") as tsv:
        writer = csv.writer(tsv, delimiter="\t", lineterminator="\n")
        writer.writerow(TSV_COLUMNS)
        for record in records:
            writer.writerow(tsv_row(record))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--vcf", required=True, help="single-sample, bcftools norm -m - VCF")
    parser.add_argument("--variant-quality", type=float, required=True, help="QUAL bar")
    parser.add_argument(
        "--frameshift-variant-quality", type=float, required=True, help="QUAL bar for frameshifts"
    )
    parser.add_argument("--variant-depth", type=int, required=True, help="minimum ALT reads")
    parser.add_argument(
        "--af-threshold", type=float, required=True, help="minimum alt/(ref+alt) ratio"
    )
    parser.add_argument("--output-vcf", required=True, help="header plus the pass records")
    parser.add_argument("--output-bed", required=True, help="positions to mask with N")
    parser.add_argument("--output-tsv", required=True, help="verdict per input record")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    headers, records = read_vcf(args.vcf)
    for record in records:
        record.verdict, record.reason = judge(
            record,
            args.variant_quality,
            args.frameshift_variant_quality,
            args.variant_depth,
            args.af_threshold,
        )
    reconcile_sites(records)
    build_mask_intervals(records)
    write_outputs(headers, records, args.output_vcf, args.output_bed, args.output_tsv)
    counts = {verdict: 0 for verdict in (PASS, MASK, DISCARD, SUPPRESS)}
    for record in records:
        counts[record.verdict] += 1
    print(
        f"filter_variants: {len(records)} records, {counts[PASS]} pass, {counts[MASK]} mask, "
        f"{counts[DISCARD]} discard, {counts[SUPPRESS]} suppress",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
