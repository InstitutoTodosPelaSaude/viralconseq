#!/usr/bin/env python3
"""Assemble the viralconseq run report (one self-contained HTML page).

Reads a finished run directory and writes a single HTML file: run-level figures, a
sortable table of every sample (or sample-segment) row, and a per-sample detail panel
with the depth trace. The page's shell -- CSS, JS, layout -- is the template passed as
``--template``; this script only produces the JSON document and splices it in at the
template's ``/*{{DATA}}*/`` marker.

``summary.tsv`` is the spine. Every other artefact is joined onto it and is optional,
so a run that is missing one still produces a page rather than a traceback:

  * ``versions.tsv``                          -- component/version pairs for the footer
  * ``run_manifest.json``                     -- provenance written by the orchestrator
  * ``qc/viralqc/outputs/results.tsv``        -- viralQC per-genome table
  * ``qc/viralqc/viralqc_status.txt``         -- how the viralQC step ended
  * ``assembly/[<seg>/]coverage_stats/<sample>.table_cov_basewise.txt`` -- depth trace
  * ``consensus/<sample>[.<seg>].fasta``      -- consensus length and N count
  * ``consensus/consensus[.<seg>].cov<T>.fasta`` -- the pooled filtered FASTA (self-check only)
  * ``assembly/[<seg>/]consensus/final_consensus/<sample>[.consensus].vcf.gz`` -- variants
  * ``qc/reports/trim.<sample>_fastp.json``   -- Illumina read QC
  * ``assembly/status/<sample>[.<seg>].txt``  -- the mapped-read gate (nanopore)
  * ``assembly/[<seg>/]clair3/<sample>/model.txt`` -- the Clair3 model that was used
  * a primer scheme BED (``--scheme``)        -- primer extent and amplicon depths

Every path is constructed from the run layout; nothing is discovered by walking the
directory. ``summary.tsv`` is read by column name, never by position, so a new column is
inert here. And the page never reports a figure the pipeline did not write: unless
``--no-self-check`` is given, the headline numbers are re-derived from ``summary.tsv``
and the build fails if they disagree.

This script runs inside a per-rule conda environment where the ``viralconseq`` package is
not installed, so it must stay on the standard library.
"""

import argparse
import csv
import gzip
import json
import math
import re
import statistics
import sys
from datetime import datetime
from pathlib import Path

NA = "NA"
PIPELINE = "viralconseq"

#: The template shipped next to this script; the workflow passes it explicitly.
TEMPLATE_PATH = Path(__file__).resolve().parents[1] / "templates" / "report.html"

#: Target number of points in a depth trace, and the smaller target used when the spine
#: is large enough that a full-resolution trace per row would bloat the page.
TRACE_POINTS = 1000
TRACE_POINTS_LARGE_RUN = 500
LARGE_RUN_ROWS = 500

#: viralQC's per-check verdicts, shown individually beside the overall status.
QC_CHECKS = [
    "privateMutations",
    "missingData",
    "mixedSites",
    "snpClusters",
    "frameShifts",
    "stopCodons",
]

#: The nanopore workflow rewrites these characters in reference contig names to ``_``
#: (Clair3 makes per-contig directories), so BED chromosome names must follow suit.
NANOPORE_UNSAFE = re.compile(r"[/\\|,~ ]")

#: ``<prefix>_<n>_(LEFT|RIGHT)`` with an optional ``_alt``-style suffix.
PRIMER_NAME = re.compile(r"^(?P<prefix>.+)_(?P<n>\d+)_(?P<side>LEFT|RIGHT)(?:_.*)?$")

#: viralQC's names for a sequence it could not assign; never a virus facet.
UNINFORMATIVE_VIRUS = {"unclassified", "unknown", "unassigned", "na", "none", ""}

STATUS_VALUES = {
    "ok",
    "no_mapped_reads",
    "empty_consensus",
    "viralqc_failed",
    "viralqc_partial",
    "viralqc_skipped",
    "viralqc_missing",
    "missing_stats",
}


# ------------------------------------------------------------------------ small readers


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def num(value, cast=float):
    """A table cell as a number, or None for the literal ``NA`` / an empty cell."""
    if value is None or value.strip() in ("", NA):
        return None
    try:
        return cast(value)
    except ValueError:
        return None


def text(value):
    """A table cell as a string, or None for the literal ``NA`` / an empty cell."""
    if value is None or value.strip() in ("", NA):
        return None
    return value.strip()


def first(*values):
    """The first value that is not None."""
    for value in values:
        if value is not None:
            return value
    return None


def read_kv(path):
    """{key: value} from a ``key<TAB>value`` file, integers cast where they parse."""
    result = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        key, _, value = line.partition("\t")
        value = value.strip()
        try:
            result[key.strip()] = int(value)
        except ValueError:
            result[key.strip()] = value
    return result


def sanitize_nanopore_contig(name):
    return NANOPORE_UNSAFE.sub("_", name)


# ---------------------------------------------------------------- per-sample artefacts


def consensus_stats(path):
    """{length, n, n_pct} for a consensus FASTA (all records concatenated)."""
    seq = "".join(
        line.strip() for line in path.read_text().splitlines() if not line.startswith(">")
    ).upper()
    if not seq:
        return None
    ns = seq.count("N")
    return {"length": len(seq), "n": ns, "n_pct": round(100.0 * ns / len(seq), 2)}


def vcf_records(path):
    """{records} for a gzip-compressed VCF: the number of non-header lines."""
    records = 0
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                records += 1
    return {"records": records}


def read_basewise(path):
    """(depths, contigs) from ``bedtools genomecov -d`` output.

    Contigs are concatenated in file order; ``contigs`` carries each one's name, length
    and offset so the page can draw the boundaries and map positions back.
    """
    depths = []
    contigs = []
    current = None
    with open(path) as handle:
        for lineno, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or not fields[0]:
                continue
            chrom, depth = fields[0], fields[2]
            if chrom != current:
                contigs.append({"name": chrom, "length": 0, "offset": len(depths)})
                current = chrom
            try:
                depths.append(int(float(depth)))
            except ValueError:
                raise SystemExit(f"{path}:{lineno}: non-numeric depth {depth!r}")
            contigs[-1]["length"] += 1
    return depths, contigs


def bin_depths(depths, target_points):
    """([mean depth per window], bin_bp): the trace, binned to about ``target_points``."""
    if not depths:
        return [], 1
    bin_bp = max(1, math.ceil(len(depths) / target_points))
    bins = [
        round(statistics.fmean(depths[i : i + bin_bp]), 1) for i in range(0, len(depths), bin_bp)
    ]
    return bins, bin_bp


def depth_mask(depths, min_depth):
    """[[start, end], ...] (1-based, inclusive) runs of positions below ``min_depth``."""
    ranges = []
    start = None
    for index, depth in enumerate(depths, start=1):
        if depth < min_depth:
            if start is None:
                start = index
        elif start is not None:
            ranges.append([start, index - 1])
            start = None
    if start is not None:
        ranges.append([start, len(depths)])
    return ranges


def fastp_block(path):
    """The read-QC figures the page shows for an Illumina sample."""
    report = json.loads(path.read_text())
    summary = report.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    filtering = report.get("filtering_result", {})
    return {
        "total_reads_before": before.get("total_reads"),
        "total_reads_after": after.get("total_reads"),
        "q30_before": before.get("q30_rate"),
        "q30_after": after.get("q30_rate"),
        "duplication": report.get("duplication", {}).get("rate"),
        "read1_mean_length": after.get("read1_mean_length"),
        "read2_mean_length": after.get("read2_mean_length"),
        "passed_filter": filtering.get("passed_filter_reads"),
        "low_quality": filtering.get("low_quality_reads"),
        "too_short": filtering.get("too_short_reads"),
        "adapter_trimmed": report.get("adapter_cutting", {}).get("adapter_trimmed_reads"),
    }


# ------------------------------------------------------------------------------ viralQC


def load_viralqc(path):
    """{sample id: [rows]} from viralQC's results table, grouped by the first ``|`` token."""
    grouped = {}
    for row in read_tsv(path):
        name = (row.get("seqName") or "").strip()
        if not name:
            continue
        grouped.setdefault(name.split("|")[0], []).append(row)
    return grouped


def pick_viralqc_row(grouped, sample, segment):
    """The row for a sample (and segment); the best-covered one when several match."""
    candidates = grouped.get(sample, [])
    if segment is not None:
        candidates = [
            row
            for row in candidates
            if row["seqName"] == f"{sample}|{segment}" or row["seqName"].endswith(f"|{segment}")
        ]
    if not candidates:
        return None

    def rank(row):
        coverage = num(row.get("coverage"))
        return -1.0 if coverage is None else coverage

    return max(candidates, key=rank)


def viralqc_block(row, spine):
    """The viralQC facts for one sample: spine columns first, results.tsv for the rest."""
    r = row or {}
    block = {
        "virus": first(text(spine.get("virus")), text(r.get("virus"))),
        "species": text(r.get("virus_species")),
        "clade": first(
            text(spine.get("clade")), text(r.get("clade")), text(r.get("clade_display"))
        ),
        "lineage": first(text(spine.get("lineage")), text(r.get("lineage"))),
        "grade": first(text(spine.get("genome_quality")), text(r.get("genomeQuality"))),
        "score": first(num(spine.get("genome_quality_score")), num(r.get("genomeQualityScore"))),
        "overall": first(text(spine.get("qc_overall_status")), text(r.get("qc.overallStatus"))),
        "coverage": num(r.get("coverage")),
        "substitutions": num(r.get("totalSubstitutions"), int),
        "deletions": num(r.get("totalDeletions"), int),
        "missing": num(r.get("totalMissing"), int),
        "checks": {check: text(r.get(f"qc.{check}.status")) for check in QC_CHECKS},
        "dataset": first(text(spine.get("viralqc_dataset")), text(r.get("dataset"))),
        "dataset_version": first(
            text(spine.get("viralqc_dataset_version")), text(r.get("datasetVersion"))
        ),
        "input_status": text(r.get("inputSequenceStatus")),
    }
    informative = {k: v for k, v in block.items() if k != "checks"}
    if row is None and all(v is None for v in informative.values()):
        return None
    return block


# -------------------------------------------------------------------------- primer scheme


def read_scheme(path, sanitize):
    """{chrom: [(start, end, name)]} from a primer BED, chrom names optionally sanitised."""
    rows = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        chrom = sanitize(fields[0]) if sanitize else fields[0]
        try:
            start, end = int(fields[1]), int(fields[2])
        except ValueError:
            continue
        name = fields[3] if len(fields) > 3 else ""
        rows.setdefault(chrom, []).append((start, end, name))
    return rows


def scheme_geometry(scheme, contigs):
    """The scheme's rows in a sample's concatenated coordinates: [(start, end, name)].

    Rows are matched to contigs by name. When nothing matches and the reference has a
    single contig the names are assumed to differ only cosmetically and every row is used.
    """
    placed = []
    for contig in contigs:
        for start, end, name in scheme.get(contig["name"], []):
            placed.append((start + contig["offset"], end + contig["offset"], name))
    if not placed and len(contigs) == 1:
        placed = [row for rows in scheme.values() for row in rows]
    return placed


def primer_extent(rows):
    """(first primer start, last primer end) over placed scheme rows, or (None, None)."""
    lo = hi = None
    for start, end, _ in rows:
        lo = start if lo is None else min(lo, start)
        hi = end if hi is None else max(hi, end)
    return lo, hi


def amplicon_spans(rows):
    """[(start, end), ...] in amplicon-number order, or None when names do not parse."""
    if not rows:
        return None
    left, right = {}, {}
    for start, end, name in rows:
        match = PRIMER_NAME.match(name)
        if not match:
            return None
        n = int(match.group("n"))
        if match.group("side") == "LEFT":
            left[n] = min(start, left.get(n, start))
        else:
            right[n] = max(end, right.get(n, end))
    spans = [(left[n], right[n]) for n in sorted(left) if n in right and right[n] > left[n]]
    return spans or None


def amplicon_means(depths, spans):
    means = []
    for start, end in spans:
        window = depths[max(0, start) : min(len(depths), end)]
        means.append(round(statistics.fmean(window), 1) if window else 0.0)
    return means


# --------------------------------------------------------------------------- run-level


def read_versions(path):
    """[[component, version], ...] from versions.tsv."""
    pairs = []
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for index, fields in enumerate(reader):
            if len(fields) < 2:
                continue
            if index == 0 and fields[0].strip().lower() == "component":
                continue
            pairs.append([fields[0].strip(), fields[1].strip()])
    return pairs


def read_manifest(path):
    raw = json.loads(path.read_text())
    inputs = []
    for sample, entries in sorted((raw.get("samples") or {}).items()):
        for entry in entries or []:
            inputs.append(
                {
                    "sample": sample,
                    "path": entry.get("path"),
                    "sha256": entry.get("sha256"),
                    "size_bytes": entry.get("size_bytes"),
                }
            )
    return {
        "viralconseq_version": raw.get("viralconseq_version"),
        "created_utc": raw.get("created_utc"),
        "finished_utc": raw.get("finished_utc"),
        "status": raw.get("status"),
        "config_sha256": raw.get("config_sha256"),
        "run_name": raw.get("run_name"),
        "inputs": inputs,
    }


# ------------------------------------------------------------------------- the document


def sample_paths(run_dir, sample, segment, data_type):
    """Every per-row artefact path, built from the run layout."""
    assembly = run_dir / "assembly" / segment if segment else run_dir / "assembly"
    stem = f"{sample}.{segment}" if segment else sample
    vcf_name = f"{sample}.consensus.vcf.gz" if data_type == "illumina" else f"{sample}.vcf.gz"
    return {
        "basewise": assembly / "coverage_stats" / f"{sample}.table_cov_basewise.txt",
        "consensus": run_dir / "consensus" / f"{stem}.fasta",
        "vcf": assembly / "consensus" / "final_consensus" / vcf_name,
        "fastp": run_dir / "qc" / "reports" / f"trim.{sample}_fastp.json",
        "gate": run_dir / "assembly" / "status" / f"{stem}.txt",
        "clair3": assembly / "clair3" / sample / "model.txt",
    }


def build(run_dir, rows, data_type, segments, params, label, versions, manifest, viralqc_status):
    """The JSON document the template renders."""
    segmented = bool(segments) or any("segment" in row for row in rows[:1])
    min_depth = params["min_depth"]
    target_points = TRACE_POINTS_LARGE_RUN if len(rows) > LARGE_RUN_ROWS else TRACE_POINTS

    viralqc_path = run_dir / "qc" / "viralqc" / "outputs" / "results.tsv"
    viralqc_rows = load_viralqc(viralqc_path) if viralqc_path.is_file() else {}

    samples = []
    coverage = {}
    raw_depths = {}

    for row in rows:
        sample = row["sample_id"].strip()
        segment = text(row.get("segment")) if segmented else None
        key = f"{sample}|{segment}" if segment else sample
        paths = sample_paths(run_dir, sample, segment, data_type)

        record = {
            "id": sample,
            "key": key,
            "segment": segment,
            "facet": None,
            "status": text(row.get("status")) or NA,
            "total_reads": num(row.get("total_reads"), int),
            "qc_passed_reads": num(row.get("qc_passed_reads"), int),
            "mapped_reads": num(row.get("mapped_reads"), int),
            "pct_mapped": num(row.get("pct_mapped")),
            "mean_depth": num(row.get("mean_depth")),
            "median_depth": num(row.get("median_depth")),
            "cov10": num(row.get("coverage_10x")),
            "cov100": num(row.get("coverage_100x")),
            "cov1000": num(row.get("coverage_1000x")),
            "cov_min": num(row.get("coverage_min_depth")),
            "consensus": None,
            "isnv_count": num(row.get("isnv_count"), int),
            "clair3_model": text(row.get("clair3_model")),
            "viralqc": None,
            "variants": None,
            "fastp": None,
            "gate": None,
        }

        spine_consensus = num(row.get("consensus_length"), int)
        if spine_consensus is not None:
            record["consensus"] = {
                "length": spine_consensus,
                "n": num(row.get("n_count"), int),
                "n_pct": num(row.get("n_pct")),
            }
        if paths["consensus"].is_file():
            stats = consensus_stats(paths["consensus"])
            if stats is not None and record["consensus"] is None:
                record["consensus"] = stats

        if paths["vcf"].is_file():
            record["variants"] = vcf_records(paths["vcf"])

        if data_type == "illumina" and paths["fastp"].is_file():
            record["fastp"] = fastp_block(paths["fastp"])

        if paths["gate"].is_file():
            gate = read_kv(paths["gate"])
            record["gate"] = {
                "status": gate.get("status"),
                "mapped_reads": gate.get("mapped_reads"),
                "minimum_mapped_reads": gate.get("minimum_mapped_reads"),
            }

        if record["clair3_model"] is None and paths["clair3"].is_file():
            model = read_kv(paths["clair3"]).get("model")
            record["clair3_model"] = str(model) if model not in (None, "") else None

        record["viralqc"] = viralqc_block(
            pick_viralqc_row(viralqc_rows, sample, segment) if viralqc_rows else None, row
        )

        if paths["basewise"].is_file():
            depths, contigs = read_basewise(paths["basewise"])
            if depths:
                bins, bin_bp = bin_depths(depths, target_points)
                coverage[key] = {
                    "bins": bins,
                    "bin_bp": bin_bp,
                    "genome_length": len(depths),
                    "mask": depth_mask(depths, min_depth),
                    "contigs": contigs,
                }
                raw_depths[key] = depths

        samples.append(record)

    # Facets: segments when the run is segmented; otherwise the virus viralQC named, but
    # only when there are at least two of them -- one virus is not a way to split a run.
    if segmented:
        facet_kind = "segment"
        facets = list(segments or [])
        for record in samples:
            record["facet"] = record["segment"]
            if record["segment"] and record["segment"] not in facets:
                facets.append(record["segment"])
    else:
        viruses = []
        for record in samples:
            virus = (record["viralqc"] or {}).get("virus")
            if virus and virus.lower() not in UNINFORMATIVE_VIRUS and virus not in viruses:
                viruses.append(virus)
        if len(viruses) >= 2:
            facet_kind = "virus"
            facets = viruses
            for record in samples:
                virus = (record["viralqc"] or {}).get("virus")
                record["facet"] = virus if virus in viruses else None
        else:
            facet_kind = "none"
            facets = []

    extents = {}
    amplicons = {}
    scheme = params.get("_scheme_rows")
    if scheme:
        groups = {}
        for record in samples:
            if record["key"] in coverage:
                groups.setdefault(record["facet"] or "*", []).append(record["key"])
        for facet, keys in groups.items():
            longest = max(keys, key=lambda k: coverage[k]["genome_length"])
            placed = scheme_geometry(scheme, coverage[longest]["contigs"])
            start, end = primer_extent(placed)
            if start is None:
                continue
            genome_length = coverage[longest]["genome_length"]
            extents[facet] = {
                "start": start,
                "end": end,
                "genome_length": genome_length,
                "outside_bp": start + max(0, genome_length - end),
            }
            spans = amplicon_spans(placed)
            if spans:
                for key in keys:
                    amplicons[key] = amplicon_means(raw_depths[key], spans)

    version = first(
        next((v for c, v in versions if c == PIPELINE), None),
        (manifest or {}).get("viralconseq_version"),
    )
    public_params = {k: v for k, v in params.items() if not k.startswith("_")}
    public_manifest = None
    if manifest:
        public_manifest = {k: v for k, v in manifest.items() if k != "viralconseq_version"}

    return {
        "run": {
            "label": label,
            "generated": datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z"),
            "output_dir": str(run_dir),
            "pipeline": PIPELINE,
            "version": version,
            "data_type": data_type,
            "segmented": segmented,
            "facet_kind": facet_kind,
            "facets": facets,
            "params": public_params,
            "versions": versions,
            "viralqc_status": viralqc_status,
            "manifest": public_manifest,
            "extents": extents,
        },
        "samples": samples,
        "coverage": coverage,
        "amplicons": amplicons,
    }


# ----------------------------------------------------------------------------- self-check


def self_check(document, run_dir, rows, threshold):
    """Refuse to write a page whose headline figures contradict summary.tsv."""
    samples = document["samples"]
    segmented = document["run"]["segmented"]
    problems = []

    if len(samples) != len(rows):
        problems.append(f"{len(samples)} rows on the page, {len(rows)} in summary.tsv")

    for record, row in zip(samples, rows):
        if record["id"] != row["sample_id"].strip():
            problems.append(f"row order diverged at {record['id']!r}")
            break

    for record, row in zip(samples, rows):
        for field, column in (
            ("cov_min", "coverage_min_depth"),
            ("median_depth", "median_depth"),
        ):
            if record[field] != num(row.get(column)):
                problems.append(f"{record['key']}: {column} does not match the page")
        if record["mapped_reads"] != num(row.get("mapped_reads"), int):
            problems.append(f"{record['key']}: mapped_reads does not match the page")

    # The filtered pooled FASTA is the pipeline's own answer to "which genomes are usable";
    # the page must arrive at the same set from summary.tsv alone.
    facets = document["run"]["facets"] if segmented else [None]
    for segment in facets:
        name = (
            f"consensus.{segment}.cov{threshold:g}.fasta"
            if segment
            else (f"consensus.cov{threshold:g}.fasta")
        )
        fasta = run_dir / "consensus" / name
        if not fasta.is_file():
            continue
        headers = [
            line[1:].split()[0]
            for line in fasta.read_text().splitlines()
            if line.startswith(">") and line[1:].strip()
        ]
        actual = {header.split("|")[0] for header in headers}
        expected = {
            s["id"]
            for s in samples
            if s["segment"] == segment and s["cov_min"] is not None and s["cov_min"] >= threshold
        }
        if actual != expected:
            problems.append(
                f"{fasta.name}: holds {len(headers)} records, the page counts "
                f"{len(expected)} at >={threshold:g}%; missing identities: "
                f"{sorted(expected - actual)}, unexpected identities: {sorted(actual - expected)}"
            )

    if problems:
        raise SystemExit("report would contradict summary.tsv:\n  " + "\n  ".join(problems))


# ----------------------------------------------------------------------------------- CLI


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--run-dir", required=True, help="A finished viralconseq run directory.")
    parser.add_argument("--template", required=True, help="The HTML shell to fill in.")
    parser.add_argument("--output", required=True, help="HTML file to write.")
    parser.add_argument("--data-type", required=True, choices=["illumina", "nanopore"])
    parser.add_argument(
        "--min-depth",
        type=int,
        required=True,
        help="Depth below which a position is masked to N (the coverage_min_depth threshold).",
    )
    parser.add_argument(
        "--consensus-coverage-threshold",
        type=float,
        required=True,
        help="Percent of the genome covered at --min-depth for a consensus to be pooled.",
    )
    parser.add_argument(
        "--segments",
        nargs="+",
        default=None,
        help="Segment names, in display order. Presence selects the segmented run layout.",
    )
    parser.add_argument("--scheme", default=None, help="Primer scheme BED, if one was used.")
    parser.add_argument("--label", default=None, help="Run name for the page header.")
    parser.add_argument("--af-threshold", type=float, default=None, help="Display only.")
    parser.add_argument("--minimum-length", type=int, default=None, help="Display only.")
    parser.add_argument("--minimum-map-quality", type=int, default=None, help="Display only.")
    parser.add_argument("--clair3-model", default=None, help="Display only.")
    parser.add_argument(
        "--no-self-check",
        action="store_true",
        help="Write the page even if it disagrees with summary.tsv. For debugging a run "
        "whose tables are already inconsistent; never for a report anyone reads.",
    )
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    run_dir = Path(args.run_dir)
    summary = run_dir / "summary.tsv"
    if not summary.is_file():
        raise SystemExit(f"{summary}: not found -- is this a finished run directory?")
    rows = read_tsv(summary)
    if rows and "sample_id" not in rows[0]:
        raise SystemExit(f"{summary}: no 'sample_id' column")

    template_path = Path(args.template)
    template = template_path.read_text(encoding="utf-8")
    marker = "/*{{DATA}}*/"
    if marker not in template:
        raise SystemExit(f"{template_path}: no {marker} marker")

    versions_path = run_dir / "versions.tsv"
    versions = read_versions(versions_path) if versions_path.is_file() else []

    manifest_path = run_dir / "run_manifest.json"
    manifest = read_manifest(manifest_path) if manifest_path.is_file() else None

    status_path = run_dir / "qc" / "viralqc" / "viralqc_status.txt"
    viralqc_status = read_kv(status_path) if status_path.is_file() else None

    scheme_rows = None
    scheme_path = Path(args.scheme) if args.scheme else None
    if scheme_path is not None and scheme_path.is_file():
        sanitize = sanitize_nanopore_contig if args.data_type == "nanopore" else None
        scheme_rows = read_scheme(scheme_path, sanitize)

    params = {
        "min_depth": args.min_depth,
        "consensus_coverage_threshold": args.consensus_coverage_threshold,
        "primer_scheme": scheme_path.name if scheme_path is not None else None,
        "af_threshold": args.af_threshold,
        "minimum_length": args.minimum_length,
        "minimum_map_quality": args.minimum_map_quality,
        "clair3_model": args.clair3_model,
        "_scheme_rows": scheme_rows,
    }

    document = build(
        run_dir=run_dir,
        rows=rows,
        data_type=args.data_type,
        segments=args.segments,
        params=params,
        label=args.label or run_dir.name,
        versions=versions,
        manifest=manifest,
        viralqc_status=viralqc_status,
    )

    if not args.no_self_check:
        self_check(document, run_dir, rows, args.consensus_coverage_threshold)

    payload = json.dumps(document, separators=(",", ":"))
    # A literal </script> inside a JSON string would close the host element early.
    page = template.replace(marker, payload.replace("</", "<\\/"))

    external = re.findall(r"(?:src|href)\s*=\s*[\"']https?:", page)
    if external:
        raise SystemExit(
            f"{template_path}: the page would load {len(external)} external "
            "resource(s); the report must render with no network"
        )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page, encoding="utf-8")

    ids = {s["id"] for s in document["samples"]}
    complete = sum(
        1
        for s in document["samples"]
        if s["cov_min"] is not None and s["cov_min"] >= args.consensus_coverage_threshold
    )
    print(
        f"{len(ids)} samples ({len(document['samples'])} rows, {complete} at "
        f">={args.consensus_coverage_threshold:g}% covered), "
        f"{len(document['coverage'])} depth traces -> {out} [{len(page) / 1024:.0f} kB]"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
