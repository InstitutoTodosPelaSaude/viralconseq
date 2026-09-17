"""Aggregate Snakemake ``*.benchmark.txt`` files into one ``benchmark.tsv``.

Run by the ``collect_benchmarks`` rule (the last rule of every workflow) inside
``envs/utils.yaml``, where the ``viralconseq`` package is NOT installed: this
script is standard-library only and takes every fact as a flag::

    python collect_benchmarks.py --logs-dir <run>/logs --samples sample-a sample-b \\
        [--segmented] --cpus map_reads=4 --exclude <own benchmark> --output benchmark.tsv

Layout contract (``LOG()``/``BENCH()`` in ``rules/common.smk``): every benchmark
lives at ``logs/<rule>/<target>.benchmark.txt`` or, in segmented workflows,
``logs/<rule>/<segment>/<target>.benchmark.txt``. ``target`` is the sample id
(``sample-<id>``) for per-sample rules and the rule name for run-level rules,
which the table reports as sample ``All``. Snakemake writes the measured
columns itself; they are taken from each file's header so a Snakemake upgrade
that adds a column widens the table instead of corrupting it. Header-only files
(a job that was killed before its first measurement) are skipped and listed on
stderr.
"""

import argparse
import csv
import os
import sys
from typing import Dict, List, Optional, Sequence

FIXED_COLUMNS = ["sample", "segment", "rule", "target", "threads"]
RUN_LEVEL_SAMPLE = "All"
NO_SEGMENT = "-"


def parse_cpus(entries: Sequence[str]) -> Dict[str, str]:
    cpus: Dict[str, str] = {}
    for entry in entries:
        rule, sep, value = entry.partition("=")
        if not sep or not rule:
            raise SystemExit(f"collect_benchmarks.py: --cpus expects RULE=N, got {entry!r}")
        cpus[rule] = value
    return cpus


def find_benchmarks(logs_dir: str, exclude: Optional[str]) -> List[str]:
    found: List[str] = []
    excluded = os.path.abspath(exclude) if exclude else None
    for root, _dirs, files in os.walk(logs_dir):
        for name in files:
            if not name.endswith(".benchmark.txt"):
                continue
            path = os.path.join(root, name)
            if excluded and os.path.abspath(path) == excluded:
                continue
            found.append(path)
    return sorted(found)


def locate(path: str, logs_dir: str, segmented: bool) -> Optional[Dict[str, str]]:
    """Return ``{rule, segment, target}`` for a benchmark path, or ``None`` when
    the path does not follow the layout contract."""
    rel = os.path.relpath(path, logs_dir)
    parts = rel.split(os.sep)
    target = parts[-1][: -len(".benchmark.txt")]
    if len(parts) == 2:
        return {"rule": parts[0], "segment": NO_SEGMENT, "target": target}
    if len(parts) == 3 and segmented:
        return {"rule": parts[0], "segment": parts[1], "target": target}
    return None


def read_rows(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader if any((v or "").strip() for v in row.values())]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--logs-dir", required=True, help="<run>/logs directory to walk")
    parser.add_argument("--output", required=True, help="benchmark.tsv to write")
    parser.add_argument(
        "--samples", nargs="*", default=[], help="configured sample ids (sample-<id>)"
    )
    parser.add_argument(
        "--segmented", action="store_true", help="expect logs/<rule>/<segment>/<target>"
    )
    parser.add_argument(
        "--cpus", nargs="*", default=[], metavar="RULE=N", help="threads given to each rule"
    )
    parser.add_argument("--exclude", default=None, help="a benchmark file to ignore (own)")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    samples = set(args.samples)
    cpus = parse_cpus(args.cpus)

    records: List[Dict[str, str]] = []
    measured_columns: List[str] = []
    skipped: List[str] = []
    for path in find_benchmarks(args.logs_dir, args.exclude):
        where = locate(path, args.logs_dir, args.segmented)
        if where is None:
            skipped.append(f"{path} (unexpected layout)")
            continue
        rows = read_rows(path)
        if not rows:
            skipped.append(f"{path} (no measurement)")
            continue
        for column in rows[0]:
            if column and column not in measured_columns:
                measured_columns.append(column)
        sample = where["target"] if where["target"] in samples else RUN_LEVEL_SAMPLE
        for row in rows:
            record = {
                "sample": sample,
                "segment": where["segment"],
                "rule": where["rule"],
                "target": where["target"],
                "threads": cpus.get(where["rule"], ""),
            }
            record.update({k: (v or "") for k, v in row.items() if k})
            records.append(record)

    # Per-sample rows first (in sample order), run-level rows last; then by
    # segment and rule so the table reads as a per-sample story.
    sample_order = {s: i for i, s in enumerate(args.samples)}

    def key(record):
        return (
            record["sample"] == RUN_LEVEL_SAMPLE,
            sample_order.get(record["sample"], len(sample_order)),
            record["segment"],
            record["rule"],
            record["target"],
        )

    records.sort(key=key)

    columns = [c for c in FIXED_COLUMNS if args.segmented or c != "segment"] + measured_columns
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        for record in records:
            writer.writerow(record)

    if skipped:
        print("collect_benchmarks.py: skipped " + ", ".join(skipped), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
