"""Build the flat ``consensus/`` directory: one self-describing FASTA per
sample, a pooled multi-FASTA without the reference, and a coverage-filtered
pooled multi-FASTA with the threshold in its name.

Run by ``rule collect_consensus`` (``rules/collect.smk``) inside
``envs/utils.yaml``; standard library only, ``main(argv)``::

    python collect_consensus.py --consensus <renamed FASTAs ...> --summary summary.tsv \\
        --threshold 70 [--segment S] --samples sample-a sample-b --output-dir <run>/consensus

Headers are normalised to ``sample-<id>[|<contig>][|<segment>]`` (the same
form viralQC reports in ``seqName``) and sequences are written on one line.
The filtered file holds the records of the samples whose ``coverage_min_depth``
in ``summary.tsv`` is at or above the threshold; a sample with ``NA`` coverage
is excluded and reported.
"""

import argparse
import csv
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

NA = "NA"
SUFFIX = ".consensus.renamed.fasta"


def read_fasta(path: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    chunks: List[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip()
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def normalise_header(header: str, sample: str, segment: Optional[str]) -> str:
    """``rename_sequences`` writes ``sample`` or ``sample_<contig>``; return
    ``sample[|contig][|segment]`` with any ``|`` inside names collapsed."""
    token = header.split()[0] if header else sample
    contig = ""
    if token != sample and token.startswith(sample + "_"):
        contig = token[len(sample) + 1 :]
    parts = [sample]
    if contig:
        parts.append(contig.replace("|", "_"))
    if segment:
        parts.append(segment.replace("|", "_"))
    return "|".join(parts)


def sample_of(path: str) -> str:
    base = os.path.basename(path)
    return base[: -len(SUFFIX)] if base.endswith(SUFFIX) else os.path.splitext(base)[0]


def coverage_by_sample(summary: str, segment: Optional[str]) -> Dict[str, str]:
    """``sample_id -> coverage_min_depth`` for the rows of this segment."""
    out: Dict[str, str] = {}
    with open(summary, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if segment is not None and row.get("segment") != segment:
                continue
            out[row["sample_id"]] = row.get("coverage_min_depth", NA)
    return out


def passes(coverage: str, threshold: float) -> bool:
    try:
        return float(coverage) >= threshold
    except (TypeError, ValueError):
        return False


def write_records(path: str, records: List[Tuple[str, str]]) -> None:
    with open(path, "w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n{seq}\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--consensus", nargs="*", default=[], help="renamed per-sample FASTAs")
    parser.add_argument("--summary", required=True, help="summary.tsv")
    parser.add_argument("--threshold", type=float, required=True, help="coverage_min_depth cutoff")
    parser.add_argument("--segment", default=None)
    parser.add_argument("--samples", nargs="+", required=True, help="sample ids in sheet order")
    parser.add_argument("--output-dir", required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    os.makedirs(args.output_dir, exist_ok=True)
    tag = f"cov{args.threshold:g}"
    suffix = f".{args.segment}" if args.segment else ""
    by_sample = {sample_of(path): path for path in args.consensus}
    coverage = coverage_by_sample(args.summary, args.segment)

    pooled: List[Tuple[str, str]] = []
    passing: List[Tuple[str, str]] = []
    excluded: List[str] = []
    for sample in args.samples:
        path = by_sample.get(sample)
        if path is None:
            print(f"collect_consensus.py: no consensus FASTA for {sample}", file=sys.stderr)
            continue
        records = [
            (normalise_header(header, sample, args.segment), seq)
            for header, seq in read_fasta(path)
        ]
        write_records(os.path.join(args.output_dir, f"{sample}{suffix}.fasta"), records)
        pooled.extend(records)
        cov = coverage.get(sample, NA)
        if passes(cov, args.threshold):
            passing.extend(records)
        elif cov == NA:
            excluded.append(f"{sample} (coverage NA)")
    write_records(os.path.join(args.output_dir, f"consensus{suffix}.fasta"), pooled)
    write_records(os.path.join(args.output_dir, f"consensus{suffix}.{tag}.fasta"), passing)
    n_pass = len({h.split("|")[0] for h, _ in passing})
    print(
        f"consensus{suffix}: {len(args.samples)} sample(s), {n_pass} at or above {args.threshold:g}% "
        f"({tag}); excluded: {', '.join(excluded) or '-'}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
