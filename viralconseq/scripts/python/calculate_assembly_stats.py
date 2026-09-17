#!/usr/bin/env python
"""Per-sample assembly statistics for the viralconseq consensus workflows.

Executed by Snakemake via a ``script:`` directive in
``rules/consensus_*_common.smk``; Snakemake injects a ``snakemake`` global
(``snakemake.input`` with named members ``raw``, ``qc_passed`` (Illumina
only), ``bam``, ``table_cov``, ``consensus``; ``snakemake.output[0]``;
``snakemake.params.minimum_depth``; ``snakemake.wildcards``). ruff/mypy
cannot see the global and are silenced for this directory in ``pyproject.toml``.
Standard library only.

Writes a one-row, headered TSV: read counts (every sequenced read, so both
mates on Illumina), mapped reads, mean and median depth, breadth of coverage
at 10x / 100x / 1000x / ``minimum_depth`` as percentages (0-100, 2 dp), and
the consensus length and N content. ``summary.tsv`` (``build_summary.py``)
joins these rows with viralQC and the rest, so the coverage arithmetic lives
here only.
"""

import gzip
import statistics
import subprocess
from typing import Dict, List, Optional, Sequence, Union

Value = Union[int, float, str]

COLUMNS = [
    "sample_id",
    "segment",
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
NA = "NA"


def count_reads(fastq: str) -> int:
    """Records in a FASTQ (plain or gzip) as ``lines // 4``.

    A line count, not a count of ``+`` separators: the separator may carry
    the read id and a one-base read's quality line can itself be ``+``.
    """
    opener = gzip.open if fastq.endswith(".gz") else open
    n_lines = 0
    with opener(fastq, "rt") as handle:
        for _ in handle:
            n_lines += 1
    return n_lines // 4


def count_mapped_reads(bam: str) -> int:
    """Primary mapped records (``samtools view -c -F 260``): each Illumina
    mate counts once, matching how ``total_reads`` counts both mates."""
    result = subprocess.run(
        ["samtools", "view", "-c", "-F", "260", bam],
        check=True,
        capture_output=True,
        text=True,
    )
    return int(result.stdout.strip())


def _pct(numerator: float, denominator: float) -> Value:
    return round(100.0 * numerator / denominator, 2) if denominator else NA


def coverage_stats(table_cov: str, minimum_depth: int) -> Dict[str, Value]:
    """Depth and breadth from a ``bedtools genomecov -d`` table (chrom, pos,
    depth; every reference position, all contigs). An empty table (a BAM
    without reads) yields zeros for the breadths and NA for the depths."""
    depths: List[int] = []
    with open(table_cov) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 3:
                depths.append(int(fields[2]))
    n = len(depths)
    if n == 0:
        return {
            "mean_depth": NA,
            "median_depth": NA,
            "coverage_10x": 0.0,
            "coverage_100x": 0.0,
            "coverage_1000x": 0.0,
            "coverage_min_depth": 0.0,
        }
    return {
        "mean_depth": round(sum(depths) / n, 2),
        "median_depth": round(statistics.median(depths), 1),
        "coverage_10x": _pct(sum(1 for d in depths if d >= 10), n),
        "coverage_100x": _pct(sum(1 for d in depths if d >= 100), n),
        "coverage_1000x": _pct(sum(1 for d in depths if d >= 1000), n),
        "coverage_min_depth": _pct(sum(1 for d in depths if d >= minimum_depth), n),
    }


def consensus_stats(fasta: str) -> Dict[str, Value]:
    """Total sequence length, N count and N % over every record of a FASTA."""
    length = 0
    n_count = 0
    with open(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            seq = line.strip()
            length += len(seq)
            n_count += seq.count("N") + seq.count("n")
    return {
        "consensus_length": length,
        "n_count": n_count,
        "n_pct": _pct(n_count, length),
    }


def build_row(
    sample_id: str,
    segment: Optional[str],
    raw_fastqs: Sequence[str],
    qc_fastqs: Optional[Sequence[str]],
    bam: str,
    table_cov: str,
    consensus: str,
    minimum_depth: int,
) -> Dict[str, Value]:
    """The per-sample row (``segment`` key present only when given)."""
    total_reads = sum(count_reads(path) for path in raw_fastqs)
    qc_passed: Value = sum(count_reads(path) for path in qc_fastqs) if qc_fastqs else NA
    mapped = count_mapped_reads(bam)
    row: Dict[str, Value] = {"sample_id": sample_id}
    if segment is not None:
        row["segment"] = segment
    row.update(
        {
            "total_reads": total_reads,
            "qc_passed_reads": qc_passed,
            "mapped_reads": mapped,
            "pct_mapped": _pct(mapped, total_reads),
        }
    )
    row.update(coverage_stats(table_cov, minimum_depth))
    row["min_depth"] = minimum_depth
    row.update(consensus_stats(consensus))
    return row


def write_tsv(row: Dict[str, Value], output: str) -> None:
    columns = [c for c in COLUMNS if c in row]
    with open(output, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        handle.write("\t".join(str(row[c]) for c in columns) + "\n")


def main_from_snakemake(smk) -> None:
    """Entry point given the injected ``snakemake`` object."""
    raw = list(smk.input.raw)
    qc_passed = list(smk.input.qc_passed) if "qc_passed" in smk.input.keys() else None
    segment = getattr(smk.wildcards, "segment", None)
    row = build_row(
        sample_id=smk.wildcards.sample,
        segment=segment,
        raw_fastqs=raw,
        qc_fastqs=qc_passed,
        bam=str(smk.input.bam),
        table_cov=str(smk.input.table_cov),
        consensus=str(smk.input.consensus),
        minimum_depth=int(smk.params.minimum_depth),
    )
    write_tsv(row, smk.output[0])


if __name__ == "__main__":
    import sys

    # Snakemake does not redirect a script:'s output; do it here so the rule log
    # holds whatever this script prints.
    if snakemake.log:
        sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    main_from_snakemake(snakemake)
