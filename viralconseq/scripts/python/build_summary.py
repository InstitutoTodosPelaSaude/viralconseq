"""Build ``summary.tsv``: one row per sample (per sample and segment in
segmented runs) joining the per-sample assembly statistics with viralQC,
iSNV counts, the nanopore status file and the Clair3 model used.

Run by ``rule summary`` (``rules/collect.smk``) inside ``envs/utils.yaml``,
where the ``viralconseq`` package is NOT installed: standard library only,
every fact arrives as a flag, ``main(argv)`` is what the tests drive::

    python build_summary.py --stats <per-sample .stats.tsv ...> \\
        --samples sample-a sample-b [--segments S L] --min-depth 20 \\
        --data-type illumina|nanopore \\
        [--viralqc-results results.tsv --viralqc-status viralqc_status.txt] \\
        [--isnvs isnvs_summary.tsv] [--status-files ...] [--model-files ...] \\
        --output summary.tsv [--legacy-csv assembly_stats_summary.csv]

The header is pinned (``COLUMNS``): downstream readers and the report rely
on it. Percentages are 0-100 (2 dp), missing values are ``NA``. Every
configured sample gets a row even when its statistics are missing
(``status missing_stats``), so the table never silently drops a sample.
"""

import argparse
import csv
import os
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

NA = "NA"

COLUMNS = [
    "sample_id",
    "segment",  # segmented runs only
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

# viralQC column -> summary column
VIRALQC_COLUMNS = {
    "virus": "virus",
    "clade": "clade",
    "lineage": "lineage",
    "genomeQuality": "genome_quality",
    "genomeQualityScore": "genome_quality_score",
    "qc.overallStatus": "qc_overall_status",
    "dataset": "viralqc_dataset",
    "datasetVersion": "viralqc_dataset_version",
}

LEGACY_COLUMNS = [
    "sample_name",
    "segment",
    "number_of_reads",
    "number_of_trim_paired_reads",
    "number_of_mapped_reads",
    "average_depth",
    "percentage_above_10x",
    "percentage_above_100x",
    "percentage_above_1000x",
    "horizontal_coverage",
]

Key = Tuple[str, Optional[str]]


def read_tsv(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def read_key_values(path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with open(path) as handle:
        for line in handle:
            key, sep, value = line.rstrip("\n").partition("\t")
            if sep:
                out[key] = value
    return out


def load_stats(paths: Iterable[str]) -> Dict[Key, Dict[str, str]]:
    stats: Dict[Key, Dict[str, str]] = {}
    for path in paths:
        for row in read_tsv(path):
            stats[(row["sample_id"], row.get("segment"))] = row
    return stats


def load_isnvs(path: Optional[str]) -> Dict[Key, str]:
    if not path:
        return {}
    out: Dict[Key, str] = {}
    for row in read_tsv(path):
        out[(row["sample"], row.get("segment"))] = row.get("number_of_isnvs", NA)
    return out


def load_status_files(paths: Iterable[str], segmented: bool) -> Dict[Key, str]:
    """``assembly/status/<sample>[.<segment>].txt`` -> its ``status`` value."""
    out: Dict[Key, str] = {}
    for path in paths:
        base = os.path.basename(path)[: -len(".txt")]
        if segmented:
            sample, _, segment = base.rpartition(".")
            key: Key = (sample, segment)
        else:
            key = (base, None)
        out[key] = read_key_values(path).get("status", NA)
    return out


def load_model_files(paths: Iterable[str], segmented: bool) -> Dict[Key, str]:
    """``assembly/[<segment>/]clair3/<sample>/model.txt`` -> its ``model`` value."""
    out: Dict[Key, str] = {}
    for path in paths:
        parts = os.path.normpath(path).split(os.sep)
        sample = parts[-2]
        segment = parts[-4] if segmented and len(parts) >= 4 else None
        out[(sample, segment)] = read_key_values(path).get("model", NA)
    return out


def load_viralqc(path: Optional[str]) -> List[Dict[str, str]]:
    return read_tsv(path) if path else []


def best_viralqc_row(
    rows: List[Dict[str, str]], sample: str, segment: Optional[str]
) -> Optional[Dict[str, str]]:
    """The viralQC row for a sample: ``seqName`` equals the key
    (``sample`` or ``sample|segment``) or, for a multi-contig reference,
    starts with ``key|``; several candidates -> the highest ``coverage``."""
    key = sample if segment is None else f"{sample}|{segment}"
    candidates = [
        row
        for row in rows
        if row.get("seqName") == key or (row.get("seqName") or "").startswith(key + "|")
    ]
    if not candidates:
        return None

    def coverage(row: Dict[str, str]) -> float:
        try:
            return float(row.get("coverage") or 0)
        except ValueError:
            return 0.0

    return max(candidates, key=coverage)


def decide_status(
    stats_row: Optional[Dict[str, str]],
    status_file_value: Optional[str],
    viralqc_status: Optional[str],
    viralqc_row: Optional[Dict[str, str]],
    run_viralqc: bool,
) -> str:
    if stats_row is None:
        return "missing_stats"
    if status_file_value and status_file_value not in ("ok", NA):
        return status_file_value
    length = stats_row.get("consensus_length", NA)
    n_pct = stats_row.get("n_pct", NA)
    if length in ("0", NA) or n_pct == "100.0":
        return "empty_consensus"
    if run_viralqc:
        if viralqc_status and viralqc_status != "ok":
            return f"viralqc_{viralqc_status}"
        if viralqc_row is None:
            return "viralqc_missing"
    return "ok"


def build_rows(
    samples: Sequence[str],
    segments: Sequence[str],
    stats: Dict[Key, Dict[str, str]],
    isnvs: Dict[Key, str],
    statuses: Dict[Key, str],
    models: Dict[Key, str],
    viralqc_rows: List[Dict[str, str]],
    viralqc_status: Optional[str],
    run_viralqc: bool,
    min_depth: int,
    data_type: str,
    clair3_model: Optional[str],
) -> List[Dict[str, str]]:
    segmented = bool(segments)
    keys: List[Key] = (
        [(s, seg) for s in samples for seg in segments]
        if segmented
        else [(s, None) for s in samples]
    )
    rows: List[Dict[str, str]] = []
    for sample, segment in keys:
        stats_row = stats.get((sample, segment))
        qc_row = best_viralqc_row(viralqc_rows, sample, segment) if run_viralqc else None
        row: Dict[str, str] = {c: NA for c in COLUMNS}
        row["sample_id"] = sample
        if segmented:
            row["segment"] = segment or NA
        else:
            row.pop("segment")
        row["min_depth"] = str(min_depth)
        if stats_row:
            for column in COLUMNS:
                if column in stats_row and column not in ("sample_id", "segment"):
                    row[column] = stats_row[column]
        row["isnv_count"] = isnvs.get((sample, segment), NA)
        if qc_row:
            for source, target in VIRALQC_COLUMNS.items():
                value = qc_row.get(source)
                row[target] = value if value not in (None, "") else NA
        if data_type == "nanopore":
            row["clair3_model"] = models.get((sample, segment)) or clair3_model or NA
        row["status"] = decide_status(
            stats_row, statuses.get((sample, segment)), viralqc_status, qc_row, run_viralqc
        )
        rows.append(row)
    return rows


def write_summary(rows: List[Dict[str, str]], output: str, segmented: bool) -> None:
    columns = [c for c in COLUMNS if segmented or c != "segment"]
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, NA) for c in columns})


def _fraction(pct: str) -> str:
    try:
        return str(round(float(pct) / 100.0, 4))
    except ValueError:
        return NA


def write_legacy_csv(rows: List[Dict[str, str]], output: str, segmented: bool) -> None:
    """The pre-0.2.0 ``assembly_stats_summary.csv`` shape, derived from the
    summary rows (fractions 0-1, mean depth), kept for one release."""
    columns = [c for c in LEGACY_COLUMNS if segmented or c != "segment"]
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            qc_passed = row.get("qc_passed_reads", NA)
            values = {
                "sample_name": row["sample_id"],
                "segment": row.get("segment", NA),
                "number_of_reads": row.get("total_reads", NA),
                "number_of_trim_paired_reads": (
                    qc_passed if qc_passed != NA else row.get("total_reads", NA)
                ),
                "number_of_mapped_reads": row.get("mapped_reads", NA),
                "average_depth": row.get("mean_depth", NA),
                "percentage_above_10x": _fraction(row.get("coverage_10x", NA)),
                "percentage_above_100x": _fraction(row.get("coverage_100x", NA)),
                "percentage_above_1000x": _fraction(row.get("coverage_1000x", NA)),
                "horizontal_coverage": _fraction(row.get("coverage_min_depth", NA)),
            }
            writer.writerow([values[c] for c in columns])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--stats", nargs="*", default=[], help="per-sample .stats.tsv files")
    parser.add_argument("--samples", nargs="+", required=True, help="sample ids in sheet order")
    parser.add_argument("--segments", nargs="*", default=[], help="segment names (segmented runs)")
    parser.add_argument("--min-depth", type=int, required=True)
    parser.add_argument("--data-type", choices=["illumina", "nanopore"], required=True)
    parser.add_argument("--viralqc-results", default=None)
    parser.add_argument("--viralqc-status", default=None)
    parser.add_argument("--isnvs", default=None)
    parser.add_argument("--status-files", nargs="*", default=[])
    parser.add_argument("--model-files", nargs="*", default=[])
    parser.add_argument("--clair3-model", default=None, help="run-wide model name (fallback)")
    parser.add_argument("--output", required=True)
    parser.add_argument("--legacy-csv", default=None)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    segmented = bool(args.segments)
    run_viralqc = args.viralqc_results is not None
    viralqc_status = None
    if args.viralqc_status and os.path.isfile(args.viralqc_status):
        viralqc_status = read_key_values(args.viralqc_status).get("status")
    rows = build_rows(
        samples=args.samples,
        segments=args.segments,
        stats=load_stats(args.stats),
        isnvs=load_isnvs(args.isnvs),
        statuses=load_status_files(args.status_files, segmented),
        models=load_model_files(args.model_files, segmented),
        viralqc_rows=load_viralqc(args.viralqc_results),
        viralqc_status=viralqc_status,
        run_viralqc=run_viralqc,
        min_depth=args.min_depth,
        data_type=args.data_type,
        clair3_model=args.clair3_model,
    )
    write_summary(rows, args.output, segmented)
    if args.legacy_csv:
        write_legacy_csv(rows, args.legacy_csv, segmented)
    flagged = [r["sample_id"] for r in rows if r["status"] != "ok"]
    print(
        f"summary: {len(rows)} row(s); {len(flagged)} not ok: {', '.join(flagged) or '-'}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
