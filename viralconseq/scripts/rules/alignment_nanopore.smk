# Nanopore alignment rules: map_reads, trim_primer_sequences
# These rules expect the following variables to be defined in the entry-point workflow:
# - REFERENCE: path to reference genome (str)
# - config: standard Snakemake config dict


rule map_reads:
    conda:
        "../envs/alignment.yaml"
    input:
        reference = REFERENCE,
        fastq = get_map_input_fastqs
    params:
        minimum_map_quality = config["minimum_map_quality"]
    output:
        bam = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/raw/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/raw/{sample}.sorted.bam.bai",
    log:
        LOG("map_reads")
    benchmark:
        BENCH("map_reads")
    threads: cpus("map_reads")
    shell:
        """
        set -euo pipefail
        exec 2> {log}
        minimap2 -a -t {threads} -x map-ont {input.reference} {input.fastq} |
        samtools view --min-MQ {params.minimum_map_quality} -bS -F 4 - |
        samtools sort -o {output.bam} -
        samtools index {output.bam} {output.bam_index}
        """

rule trim_primer_sequences:
    conda:
        "../envs/alignment.yaml"
    input:
        bam = rules.map_reads.output.bam,
        bam_index = rules.map_reads.output.bam_index
    output:
        bam = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/trimmed/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/trimmed/{sample}.sorted.bam.bai",
        trimmed_info = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/trimmed/{sample}.trimmed.txt"
    params:
        bed = config.get("scheme", "NA"),
        minimum_length = config["minimum_length"],
        path = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/raw/"
    log:
        LOG("trim_primer_sequences")
    benchmark:
        BENCH("trim_primer_sequences")
    threads: cpus("trim_primer_sequences")
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        if [ "{params.bed}" = "NA" ]; then
            cp {input.bam} {output.bam};
            cp {input.bam_index} {output.bam_index};
            touch {output.trimmed_info};
            echo "No primer scheme detected, assuming sequence data came from untargeted sequencing approach (bam files copied to trimmed dir for convenience)." > {params.path}notes.txt
        else
            samtools ampliconclip \
                --both-ends \
                --hard-clip \
                --filter-len {params.minimum_length} \
                -b {params.bed} \
                -f {output.trimmed_info} \
                {input.bam} > {output.bam}
            samtools index {output.bam} {output.bam_index}
        fi
        """
