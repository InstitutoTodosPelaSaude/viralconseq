# Collection step: the run-level summary.tsv and the flat consensus/ directory.
#
# Included by the four entry-point workflows after rules/viralqc.smk and before
# organize_files. Expects: config (output, data, samples, minimum_depth,
# run_viralqc, run_isnv), SEGMENT_WILDCARD / SEGMENTS, PY / LOG / BENCH from
# rules/common.smk, rules.calculate_assembly_statistics, rules.run_viralqc and,
# on nanopore, rules.check_mapped_reads / rules.infer_consensus_sequence.
#
# summary.tsv has one row per sample (per sample and segment when segmented),
# a pinned header (scripts/python/build_summary.py COLUMNS), percentages 0-100
# and NA for missing values, plus a status column. The pre-0.2.0
# assembly/assembly_stats_summary.csv is still written from the same rows for
# one release (deprecated).

import os

COLLECT_DATA = str(config.get("data", "illumina"))
COLLECT_RUN_VIRALQC = bool(config.get("run_viralqc", True))
COLLECT_RUN_ISNV = COLLECT_DATA == "illumina" and bool(config.get("run_isnv", False))
SEGMENT_KEYS = list(SEGMENTS.keys()) if SEGMENT_WILDCARD else []
COLLECT_COV_T = float(config["consensus_coverage_threshold"])
COLLECT_COV_TAG = f"cov{COLLECT_COV_T:g}"
# ".{segment}" in per-file names of segmented runs, "" otherwise.
COLLECT_SEG_SUFFIX = ".{segment}" if SEGMENT_WILDCARD else ""


def _per_sample(pattern):
    """Expand a per-sample (and, when segmented, per-segment) output pattern."""
    if SEGMENT_WILDCARD:
        return expand(pattern, sample=config["samples"], segment=SEGMENT_KEYS)
    return expand(pattern, sample=config["samples"])


rule summary:
    conda:
        "../envs/utils.yaml"
    input:
        script = os.path.join(PY, "build_summary.py"),
        stats = _per_sample(rules.calculate_assembly_statistics.output.stats),
        viralqc = rules.run_viralqc.output.results if COLLECT_RUN_VIRALQC else [],
        viralqc_status = rules.run_viralqc.output.status if COLLECT_RUN_VIRALQC else [],
        # Literal path: summarize_isnvs lives in the entry file, after this include.
        isnvs = config["output"] + "isnvs/isnvs_summary.tsv" if COLLECT_RUN_ISNV else [],
        status_files = (
            _per_sample(rules.check_mapped_reads.output.status) if COLLECT_DATA == "nanopore" else []
        ),
        model_files = (
            _per_sample(rules.infer_consensus_sequence.output.model_txt)
            if COLLECT_DATA == "nanopore"
            else []
        ),
    output:
        summary = config["output"] + "summary.tsv",
        # Deprecated alias of the old shape; removed in 0.3.0.
        legacy = config["output"] + "assembly/assembly_stats_summary.csv",
    params:
        samples = list(config["samples"].keys()),
        segments = (["--segments"] + SEGMENT_KEYS) if SEGMENT_KEYS else [],
        min_depth = config["minimum_depth"],
        data_type = COLLECT_DATA,
        viralqc = (
            (lambda wildcards, input: ["--viralqc-results", str(input.viralqc), "--viralqc-status", str(input.viralqc_status)])
            if COLLECT_RUN_VIRALQC
            else []
        ),
        isnvs = (lambda wildcards, input: ["--isnvs", str(input.isnvs)]) if COLLECT_RUN_ISNV else [],
        status_files = (
            (lambda wildcards, input: ["--status-files"] + list(input.status_files))
            if COLLECT_DATA == "nanopore"
            else []
        ),
        model_files = (
            (lambda wildcards, input: ["--model-files"] + list(input.model_files))
            if COLLECT_DATA == "nanopore"
            else []
        ),
    log:
        LOG("summary", target="summary", per_segment=False)
    benchmark:
        BENCH("summary", target="summary", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        python {input.script:q} --stats {input.stats:q} --samples {params.samples:q} {params.segments:q} \
            --min-depth {params.min_depth} --data-type {params.data_type} \
            {params.viralqc:q} {params.isnvs:q} {params.status_files:q} {params.model_files:q} \
            --output {output.summary:q} --legacy-csv {output.legacy:q}
        """


rule collect_consensus:
    # consensus/<sample>[.<segment>].fasta (one-line sequences, headers
    # sample-<id>[|<contig>][|<segment>]), consensus[.<segment>].fasta pooled
    # without the reference, and consensus[.<segment>].cov<T>.fasta holding
    # the samples at or above the coverage threshold.
    conda:
        "../envs/utils.yaml"
    input:
        script = os.path.join(PY, "collect_consensus.py"),
        consensus = expand(
            rules.rename_sequences.output.consensus_renamed,
            sample=config["samples"],
            allow_missing=True,
        ),
        summary = rules.summary.output.summary,
    output:
        per_sample = expand(
            config["output"] + "consensus/{sample}" + COLLECT_SEG_SUFFIX + ".fasta",
            sample=config["samples"],
            allow_missing=True,
        ),
        pooled = config["output"] + "consensus/consensus" + COLLECT_SEG_SUFFIX + ".fasta",
        filtered = config["output"] + "consensus/consensus" + COLLECT_SEG_SUFFIX + "." + COLLECT_COV_TAG + ".fasta",
    params:
        outdir = config["output"] + "consensus",
        threshold = COLLECT_COV_T,
        samples = list(config["samples"].keys()),
        segment = (lambda wildcards: ["--segment", wildcards.segment]) if SEGMENT_WILDCARD else [],
    log:
        LOG("collect_consensus", target="collect_consensus", per_segment=True)
    benchmark:
        BENCH("collect_consensus", target="collect_consensus", per_segment=True)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        python {input.script:q} --consensus {input.consensus:q} --summary {input.summary:q} \
            --threshold {params.threshold} {params.segment:q} --samples {params.samples:q} \
            --output-dir {params.outdir:q}
        """
