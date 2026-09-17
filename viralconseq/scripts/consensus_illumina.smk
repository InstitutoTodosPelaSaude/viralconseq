# Shared preamble first: exact-sample wildcard constraints and run.log hooks.
include: "rules/common.smk"

SEGMENT_WILDCARD = ""
REFERENCE = config["reference"]

# The analysis products; everything else (versions.tsv, config.yml,
# benchmark.tsv) is provenance about them. rules/provenance.smk's
# collect_benchmarks depends on this list so it runs last.
TERMINAL_INPUTS = [
    config['output'] + "assembly/consensus/final_consensus/samples_alignment.fasta",
    config['output'] + "isnvs/isnvs_summary.tsv" if config.get("run_isnv", False) else [],
    config['output'] + "qc/viralqc/outputs/results.tsv" if config.get("run_viralqc", True) else [],
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

include: "rules/qc_illumina.smk"
include: "rules/alignment_illumina.smk"
include: "rules/consensus_illumina.smk"
include: "rules/stats.smk"
include: "rules/consensus_illumina_common.smk"
include: "rules/viralqc.smk"

rule unify_assembly_statistics_reports:
    conda:
        "envs/utils.yaml"
    input:
        reports = expand(rules.calculate_assembly_statistics.output.stats_summary, sample=config["samples"])
    output:
        unified_stats_summary = config['output'] + "assembly/assembly_stats_summary.csv"
    log:
        LOG("unify_assembly_statistics_reports", target="unify_assembly_statistics_reports", per_segment=False)
    benchmark:
        BENCH("unify_assembly_statistics_reports", target="unify_assembly_statistics_reports", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        echo \"sample_name,number_of_reads,number_of_trim_paired_reads,number_of_mapped_reads,average_depth,percentage_above_10x,percentage_above_100x,percentage_above_1000x,horizontal_coverage\" > {output.unified_stats_summary} ;
        cat {input.reports} >> {output.unified_stats_summary}
        """

rule summarize_isnvs:
    conda:
        "envs/utils.yaml"
    input:
        vcf_files = expand(rules.detect_isnv.output.vcf, sample=config["samples"])
    output:
        isnvs_summary = config['output'] +  "isnvs/isnvs_summary.tsv"
    log:
        LOG("summarize_isnvs", target="summarize_isnvs", per_segment=False)
    benchmark:
        BENCH("summarize_isnvs", target="summarize_isnvs", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        echo -e "sample\\tnumber_of_isnvs" > {output.isnvs_summary};
        for _file in {input.vcf_files}; do
            sample=$(basename $_file .isnvs.vcf.gz);
            isnv_count=$(bcftools view -H $_file | wc -l);
            echo -e "$sample\\t$isnv_count" >> {output.isnvs_summary};
        done
        """

rule organize_files:
    conda:
        "envs/utils.yaml"
    input:
        fastp_reports = expand(rules.perform_qc.output.html, sample=config["samples"]),
        vcf_files = expand(rules.generate_vcf_consensus.output.vcf, sample=config["samples"]),
        isn_vcf_files = expand(rules.detect_isnv.output.vcf, sample=config["samples"]) if config.get("run_isnv", False) else [],
        stats_summary = expand(rules.calculate_assembly_statistics.output.stats_summary, sample=config["samples"]),
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
        for _file in {input.fastp_reports}; do
            sample=$(basename $_file _fastp.html | sed 's/^trim.//');
            ln -sf $_file {params.outdir}samples/$sample/fastp.html;
        done
        for _file in {input.vcf_files}; do
            sample=$(basename $_file .consensus.vcf.gz);
            ln -sf $_file {params.outdir}samples/$sample/consensus.vcf.gz;
            ln -sf $_file.tbi {params.outdir}samples/$sample/consensus.vcf.gz.tbi;
        done
        for _file in {input.isn_vcf_files} ""; do
            if [ -z "$_file" ]; then continue; fi
            sample=$(basename $_file .isnvs.vcf.gz);
            ln -sf $_file {params.outdir}samples/$sample/isnvs.vcf.gz;
            ln -sf $_file.tbi {params.outdir}samples/$sample/isnvs.vcf.gz.tbi;
        done
        for _file in {input.stats_summary}; do
            sample=$(basename $_file .stats_summary.csv);
            ln -sf $_file {params.outdir}samples/$sample/stats_summary.csv;
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
