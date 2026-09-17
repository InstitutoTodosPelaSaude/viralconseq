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
        config["output"] + "logs/consensus_nanopore/sanitize_reference/reference.log"
    benchmark:
        config["output"] + "logs/consensus_nanopore/sanitize_reference/reference.benchmark.txt"
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        mkdir -p $(dirname {output.fasta})
        sed '/^>/s/[\\/|,~ ]/_/g' {input} > {output.fasta}
        samtools faidx {output.fasta}
        """

REFERENCE = rules.sanitize_reference.output.fasta

rule all:
    default_target: True
    input:
        config['output'] + "assembly/consensus/final_consensus/samples_alignment.fasta",
        config['output'] + "qc/viralqc/outputs/results.tsv" if config.get("run_viralqc", True) else [],
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

# ``calculate_assembly_statistics`` and ``align_consensus_to_reference_genome``
# are defined in the included ``consensus_nanopore_common.smk``. The
# ``calculate_assembly_stats.py`` helper expects three fastq inputs
# (raw_r1, raw_r2, trimmed) — Nanopore passes the same fastq for all three.

rule unify_assembly_statistics_reports:
    conda:
        "envs/utils.yaml"
    input:
        reports = expand(rules.calculate_assembly_statistics.output.stats_summary, sample=config["samples"])
    output:
        unified_stats_summary = config['output'] + "assembly/assembly_stats_summary.csv"
    log:
        config['output'] + "logs/consensus_nanopore/unify_assembly_statistics_reports/unify_assembly_statistics_reports.log"
    benchmark:
        config['output'] + "logs/consensus_nanopore/unify_assembly_statistics_reports/unify_assembly_statistics_reports.benchmark.txt"
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        echo \"sample_name,number_of_reads,number_of_trim_paired_reads,number_of_mapped_reads,average_depth,percentage_above_10x,percentage_above_100x,percentage_above_1000x,horizontal_coverage\" > {output.unified_stats_summary} ;
        cat {input.reports} >> {output.unified_stats_summary}
        """


rule organize_files:
    conda:
        "envs/utils.yaml"
    input:
        vcf_files = expand(rules.infer_consensus_sequence.output.vcf, sample=config["samples"]),
        vcf_raw_files = expand(rules.infer_consensus_sequence.output.vcf_raw, sample=config["samples"]),
        table_cov = expand(rules.calculate_coverage_basewise.output.table_cov, sample=config["samples"]),
        consensus_files = expand(rules.rename_sequences.output.consensus_renamed, sample=config["samples"]),
        raw_mapped_reads = expand(rules.map_reads.output.bam, sample=config["samples"]),
        trimmed_mapped_reads = expand(rules.trim_primer_sequences.output.bam, sample=config["samples"]),
        viralqc_files = expand(rules.split_viralqc_results.output.tsv, sample=config["samples"]) if config.get("run_viralqc", True) else [],
    output:
        config['output'] + "benchmark.tsv"
    params:
        outdir = config['output'],
        samples = " ".join(config["samples"].keys()),
        own_benchmark = config['output'] + "logs/consensus_nanopore/organize_files/organize_files.benchmark.txt"
    log:
        config['output'] + "logs/consensus_nanopore/organize_files/organize_files.log"
    benchmark:
        config['output'] + "logs/consensus_nanopore/organize_files/organize_files.benchmark.txt"
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

        # Benchmark aggregation. Drop this rule's own benchmark from a previous
        # run first: Snakemake clears log: but not benchmark: before a job, so a
        # rerun would otherwise ingest a stale row for organize_files itself.
        rm -f {params.own_benchmark}
        echo -e "sample\\ttask\\tseconds\\th:m:s\\tmax_rss\\tmax_vms\\tmax_uss\\tmax_pss\\tio_in\\tio_out\\tmean_load\\tcpu_time" > {output}
        find {params.outdir} -name "*.benchmark.txt" | while read -r file; do
            task=$(basename $(dirname $file))
            sample=$(basename $file .benchmark.txt)

            matched=false
            for s in {params.samples}; do
                if [[ "$sample" == "$s" ]]; then
                    matched=true
                    break
                fi
            done

            if [[ "$matched" == "false" ]]; then
                sample="All"
            fi

            tail -n +2 $file | awk -v sample=$sample -v task=$task '{{print sample"\\t"task"\\t"$0}}' >> {output}
        done
        """
