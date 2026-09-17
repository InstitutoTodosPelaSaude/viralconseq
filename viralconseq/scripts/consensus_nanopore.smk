# Shared preamble first: exact-sample wildcard constraints and run.log hooks.
include: "rules/common.smk"

SEGMENT_WILDCARD = ""
rule sanitize_reference:
    conda:
        "envs/clair3.yaml"
    input: config["reference"]
    output: 
        fasta = config["output"] + "reference/reference.sanitized.fasta",
        fai = config["output"] + "reference/reference.sanitized.fasta.fai"
    log:
        LOG("sanitize_reference", target="reference", per_segment=False)
    benchmark:
        BENCH("sanitize_reference", target="reference", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        mkdir -p $(dirname {output.fasta})
        sed '/^>/s/[\\/|,~ ]/_/g' {input} > {output.fasta}
        samtools faidx {output.fasta}
        """

REFERENCE = rules.sanitize_reference.output.fasta

# The analysis products; everything else (versions.tsv, config.yml,
# benchmark.tsv) is provenance about them. rules/provenance.smk's
# collect_benchmarks depends on this list so it runs last.
TERMINAL_INPUTS = [
    config['output'] + "assembly/consensus/final_consensus/samples_alignment.fasta",
    config['output'] + "qc/viralqc/outputs/results.tsv" if config.get("run_viralqc", True) else [],
    config['output'] + "summary.tsv",
]

rule all:
    default_target: True
    input:
        TERMINAL_INPUTS,
        config['output'] + "versions.tsv",
        config['output'] + "config.yml",
        config['output'] + "benchmark.tsv"

def get_map_input_fastqs(wildcards):
    reads = config["samples"][wildcards.sample]
    if isinstance(reads, str):
        reads = reads.split()
    return reads

include: "rules/alignment_nanopore.smk"
include: "rules/consensus_nanopore.smk"
include: "rules/stats.smk"
include: "rules/consensus_nanopore_common.smk"
include: "rules/viralqc.smk"
include: "rules/collect.smk"

# ``calculate_assembly_statistics`` and ``align_consensus_to_reference_genome``
# are defined in the included ``consensus_nanopore_common.smk``. The
# ``calculate_assembly_stats.py`` helper expects three fastq inputs
# (raw_r1, raw_r2, trimmed) — Nanopore passes the same fastq for all three.

rule organize_files:
    conda:
        "envs/utils.yaml"
    input:
        vcf_files = expand(rules.infer_consensus_sequence.output.vcf, sample=config["samples"]),
        vcf_raw_files = expand(rules.infer_consensus_sequence.output.vcf_raw, sample=config["samples"]),
        model_files = expand(rules.infer_consensus_sequence.output.model_txt, sample=config["samples"]),
        status_files = expand(rules.check_mapped_reads.output.status, sample=config["samples"]),
        stats_files = expand(rules.calculate_assembly_statistics.output.stats, sample=config["samples"]),
        table_cov = expand(rules.calculate_coverage_basewise.output.table_cov, sample=config["samples"]),
        consensus_files = expand(rules.rename_sequences.output.consensus_renamed, sample=config["samples"]),
        raw_mapped_reads = expand(rules.map_reads.output.bam, sample=config["samples"]),
        trimmed_mapped_reads = expand(rules.trim_primer_sequences.output.bam, sample=config["samples"]),
        viralqc_files = expand(rules.split_viralqc_results.output.tsv, sample=config["samples"]) if config.get("run_viralqc", True) else [],
    output:
        # Sentinel: the symlink tree has no single file to declare. benchmark.tsv
        # is produced by collect_benchmarks (rules/provenance.smk).
        sentinel = touch(config['output'] + "samples/.organized")
    params:
        outdir = config['output'],
        samples = " ".join(config["samples"].keys())
    log:
        LOG("organize_files", target="organize_files", per_segment=False)
    benchmark:
        BENCH("organize_files", target="organize_files", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        mkdir -p {params.outdir}samples/
        for sample in {params.samples}; do
            mkdir -p {params.outdir}samples/$sample;
        done
        for _file in {input.vcf_files}; do
            sample=$(basename $_file .vcf.gz);
            ln -sf $_file {params.outdir}samples/$sample/consensus.vcf.gz;
            ln -sf $_file.tbi {params.outdir}samples/$sample/consensus.vcf.gz.tbi;
        done
        for _file in {input.vcf_raw_files}; do
            sample=$(basename $_file .raw.vcf.gz);
            ln -sf $_file {params.outdir}samples/$sample/raw.vcf.gz;
            ln -sf $_file.tbi {params.outdir}samples/$sample/raw.vcf.gz.tbi;
        done
        for _file in {input.model_files}; do
            sample=$(basename $(dirname $_file));
            ln -sf $_file {params.outdir}samples/$sample/clair3_model.txt;
        done
        for _file in {input.status_files}; do
            sample=$(basename $_file .txt);
            ln -sf $_file {params.outdir}samples/$sample/status.txt;
        done
        for _file in {input.stats_files}; do
            sample=$(basename $_file .stats.tsv);
            ln -sf $_file {params.outdir}samples/$sample/stats.tsv;
        done
        for _file in {input.table_cov}; do
            sample=$(basename $_file .table_cov_basewise.txt);
            ln -sf $_file {params.outdir}samples/$sample/table_cov_basewise.txt;
        done
        for _file in {input.consensus_files}; do
            sample=$(basename $_file .consensus.renamed.fasta);
            ln -sf $_file {params.outdir}samples/$sample/consensus.fasta;
        done
        for _file in {input.raw_mapped_reads}; do
            sample=$(basename $_file .sorted.bam);
            ln -sf $_file {params.outdir}samples/$sample/raw_mapped_reads.bam;
            ln -sf $_file.bai {params.outdir}samples/$sample/raw_mapped_reads.bam.bai;
        done
        for _file in {input.trimmed_mapped_reads}; do
            sample=$(basename $_file .sorted.bam);
            ln -sf $_file {params.outdir}samples/$sample/trimmed_mapped_reads.bam;
            ln -sf $_file.bai {params.outdir}samples/$sample/trimmed_mapped_reads.bam.bai;
        done

        for _file in {input.viralqc_files} ""; do
            if [ -z "$_file" ]; then continue; fi
            sample=$(basename $_file .viralqc.tsv);
            ln -sf $_file {params.outdir}samples/$sample/viralqc.tsv;
        done

        """

# Provenance and the benchmark collector last: collect_benchmarks depends on
# TERMINAL_INPUTS and on rules.organize_files, both defined above.
include: "rules/provenance.smk"
