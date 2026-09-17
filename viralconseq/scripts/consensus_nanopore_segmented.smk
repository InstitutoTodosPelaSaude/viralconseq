# Shared preamble first: exact-sample wildcard constraints and run.log hooks.
include: "rules/common.smk"

SEGMENTS = config["reference"]  # dict: {"S": "/path/S.fa", "L": "/path/L.fa", ...}

# The analysis products; everything else (versions.tsv, config.yml,
# benchmark.tsv) is provenance about them. rules/provenance.smk's
# collect_benchmarks depends on this list so it runs last.
TERMINAL_INPUTS = [
    expand(
        config['output'] + "assembly/{segment}/consensus/final_consensus/samples_alignment.fasta",
        segment=SEGMENTS.keys()
    ),
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

rule sanitize_reference:
    conda:
        "envs/clair3.yaml"
    input: lambda wildcards: SEGMENTS[wildcards.segment]
    output:
        fasta = config["output"] + "reference/{segment}.sanitized.fasta",
        fai = config["output"] + "reference/{segment}.sanitized.fasta.fai"
    log:
        LOG("sanitize_reference", target="{segment}", per_segment=False)
    benchmark:
        BENCH("sanitize_reference", target="{segment}", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        mkdir -p $(dirname {output.fasta})
        sed '/^>/s/[\\/|,~ ]/_/g' {input} > {output.fasta}
        samtools faidx {output.fasta}
        """

def get_map_input_fastqs(wildcards):
    reads = config["samples"][wildcards.sample]
    if isinstance(reads, str):
        reads = reads.split()
    return reads

REFERENCE = rules.sanitize_reference.output.fasta
SEGMENT_WILDCARD = "{segment}/"

include: "rules/alignment_nanopore.smk"
include: "rules/consensus_nanopore.smk"
include: "rules/stats.smk"
include: "rules/consensus_nanopore_common.smk"
include: "rules/viralqc.smk"
include: "rules/collect.smk"

rule organize_files:
    conda:
        "envs/utils.yaml"
    input:
        vcf_files = expand(
            rules.infer_consensus_sequence.output.vcf,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        vcf_raw_files = expand(
            rules.infer_consensus_sequence.output.vcf_raw,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        model_files = expand(
            rules.infer_consensus_sequence.output.model_txt,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        status_files = expand(
            rules.check_mapped_reads.output.status,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        stats_files = expand(
            rules.calculate_assembly_statistics.output.stats,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        table_cov = expand(
            rules.calculate_coverage_basewise.output.table_cov,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        consensus_files = expand(
            rules.rename_sequences.output.consensus_renamed,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        raw_mapped_reads = expand(
            rules.map_reads.output.bam,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        trimmed_mapped_reads = expand(
            rules.trim_primer_sequences.output.bam,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        viralqc_files = expand(
            rules.split_viralqc_results.output.tsv, sample=config["samples"]
        ) if config.get("run_viralqc", True) else [],
    output:
        # Sentinel: the symlink tree has no single file to declare. benchmark.tsv
        # is produced by collect_benchmarks (rules/provenance.smk).
        sentinel = touch(config['output'] + "samples/.organized")
    params:
        outdir = config['output'],
        samples = " ".join(config["samples"].keys()),
        segments = " ".join(SEGMENTS.keys())
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
            for segment in {params.segments}; do
                mkdir -p {params.outdir}samples/$sample/$segment;
            done
        done
        for _file in {input.vcf_files}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .vcf.gz);
            ln -sf $_file {params.outdir}samples/$sample/$segment/consensus.vcf.gz;
            ln -sf $_file.tbi {params.outdir}samples/$sample/$segment/consensus.vcf.gz.tbi;
        done
        for _file in {input.vcf_raw_files}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo "$rel" | cut -d'/' -f1);
            sample=$(basename $_file .raw.vcf.gz);
            ln -sf $_file {params.outdir}samples/$sample/$segment/raw.vcf.gz;
            ln -sf $_file.tbi {params.outdir}samples/$sample/$segment/raw.vcf.gz.tbi;
        done
        for _file in {input.model_files}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo "$rel" | cut -d'/' -f1);
            sample=$(basename $(dirname $_file));
            ln -sf $_file {params.outdir}samples/$sample/$segment/clair3_model.txt;
        done
        for _file in {input.status_files}; do
            # assembly/status/<sample>.<segment>.txt
            base=$(basename $_file .txt); segment=${{base##*.}}; sample=${{base%.*}};
            ln -sf $_file {params.outdir}samples/$sample/$segment/status.txt;
        done
        for _file in {input.stats_files}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo "$rel" | cut -d'/' -f1);
            sample=$(basename $_file .stats.tsv);
            ln -sf $_file {params.outdir}samples/$sample/$segment/stats.tsv;
        done
        for _file in {input.table_cov}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .table_cov_basewise.txt);
            ln -sf $_file {params.outdir}samples/$sample/$segment/table_cov_basewise.txt;
        done
        for _file in {input.consensus_files}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .consensus.renamed.fasta);
            ln -sf $_file {params.outdir}samples/$sample/$segment/consensus.fasta;
        done
        for _file in {input.raw_mapped_reads}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .sorted.bam);
            ln -sf $_file {params.outdir}samples/$sample/$segment/raw_mapped_reads.bam;
            ln -sf $_file.bai {params.outdir}samples/$sample/$segment/raw_mapped_reads.bam.bai;
        done
        for _file in {input.trimmed_mapped_reads}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .sorted.bam);
            ln -sf $_file {params.outdir}samples/$sample/$segment/trimmed_mapped_reads.bam;
            ln -sf $_file.bai {params.outdir}samples/$sample/$segment/trimmed_mapped_reads.bam.bai;
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
