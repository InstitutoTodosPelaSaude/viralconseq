# Illumina consensus rules: infer_consensus_sequence, detect_isnv, generate_vcf_consensus
# These rules expect the following variables to be defined in the entry-point workflow:
# - REFERENCE: path to reference genome (str)
# - config: standard Snakemake config dict


rule detect_isnv:
    conda:
        "../envs/consensus.yaml"
    input:
        reference = REFERENCE,
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index
    output:
        bam = temp(config['output'] + "assembly/" + SEGMENT_WILDCARD + "isnvs/{sample}.lofreq.sorted.bam"),
        bam_index = temp(config['output'] + "assembly/" + SEGMENT_WILDCARD + "isnvs/{sample}.lofreq.sorted.bam.bai"),
        vcf_tmp = temp(config['output'] + "assembly/" + SEGMENT_WILDCARD + "isnvs/{sample}.lofreq.tmp.vcf"),
        vcf = config['output'] + "assembly/" + SEGMENT_WILDCARD + "isnvs/{sample}.isnvs.vcf.gz",
        vcf_index = config['output'] + "assembly/" + SEGMENT_WILDCARD + "isnvs/{sample}.isnvs.vcf.gz.tbi"
    params:
        af_min_threshold = config.get("af_isnv_threshold", 0.05)
    log:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/lofreq/{sample}.log"
    benchmark:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/lofreq/{sample}.benchmark.txt"
    threads: config.get("detect_isnv_cpus", 2)
    resources:
        mem_mb = config.get("detect_isnv_ram", 4) * 1024
    shell:
        """
        set -euo pipefail
        lofreq indelqual \
            -f {input.reference} \
            -o {output.bam} \
            --dindel \
            {input.bam}
        samtools index {output.bam} {output.bam_index}

        lofreq call-parallel \
            --pp-threads {threads} \
            --call-indels \
            -f {input.reference} \
            -o {output.vcf_tmp} \
            {output.bam}

        bcftools view -i 'INFO/AF<0.5 & INFO/AF>={params.af_min_threshold}' {output.vcf_tmp} -Oz -o {output.vcf}
        tabix {output.vcf}
        """

rule infer_consensus_sequence:
    conda:
        "../envs/alignment.yaml"
    input:
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index
    output:
        consensus = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.consensus.fasta"
    params:
        minimum_depth = config.get("minimum_depth", 10),
        af_threshold = config.get("af_threshold", 0.5)
    benchmark:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/samtools/consensus/{sample}.benchmark.txt"
    log:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/samtools/consensus/{sample}.log"
    shell:
        "samtools consensus -a -d {params.minimum_depth} -m simple -q -c {params.af_threshold} --show-ins yes {input.bam} -o {output.consensus}"

rule generate_vcf_consensus:
    conda:
        "../envs/alignment.yaml"
    input:
        reference = REFERENCE,
        consensus = rules.infer_consensus_sequence.output.consensus
    output:
        vcf = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.consensus.vcf.gz",
        vcf_index = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.consensus.vcf.gz.tbi"
    benchmark:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/gsaalign/{sample}.benchmark.txt"
    log:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/gsaalign/{sample}.log"
    shell:
        """
        set -euo pipefail
        out_prefix=$(echo {output.vcf} | sed 's/.vcf.gz//')

        write_mock_vcf() {{
            echo "##fileformat=VCFv4.2" > $out_prefix.vcf
            printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" >> $out_prefix.vcf
            bgzip -f $out_prefix.vcf
            tabix -p vcf {output.vcf} || touch {output.vcf_index}
        }}

        # Only align when the consensus has real bases. "N" is not sequence
        # content, so an all-N / near-empty consensus (a zero-coverage sample
        # against a divergent reference) must NOT be sent to GSAlign -- it would
        # produce no VCF and fail the whole run under `set -euo pipefail`.
        if grep -v "^>" {input.consensus} | grep -qi "[ACGT]"; then
            # Clear any stale VCF from a previous run so the emptiness check below
            # reflects only this GSAlign invocation (the raw .vcf is not a tracked
            # Snakemake output, so it is not cleaned automatically).
            rm -f "$out_prefix.vcf"
            # GSAlign may legitimately emit no VCF on a degenerate query; capture
            # its exit status so a genuine tool failure is logged distinctly, but
            # neither case aborts the run -- both fall back to a mock VCF.
            gsalign_rc=0
            GSAlign \
                -r {input.reference} \
                -q {input.consensus} \
                -o $out_prefix \
                -fmt 1 \
                -sen || gsalign_rc=$?
            rm -f $out_prefix.maf
            if [ "$gsalign_rc" -ne 0 ]; then
                echo "Warning: GSAlign exited with status $gsalign_rc for {wildcards.sample}; falling back to a mock VCF." >&2
            fi
            if [ -s "$out_prefix.vcf" ]; then
                bgzip -f $out_prefix.vcf
                tabix -p vcf {output.vcf} || touch {output.vcf_index}
            else
                echo "Warning: GSAlign produced no VCF for {wildcards.sample}. Creating a mock VCF." >&2
                write_mock_vcf
            fi
        else
            echo "Warning: Consensus sequence for {wildcards.sample} is empty. Creating a mock VCF." >&2
            write_mock_vcf
        fi
        """
