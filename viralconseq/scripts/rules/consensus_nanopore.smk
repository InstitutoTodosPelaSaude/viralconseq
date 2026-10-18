# Nanopore consensus rules: infer_consensus_sequence (clair3 + bcftools)
# These rules expect the following variables to be defined in the entry-point workflow:
# - REFERENCE: path to reference genome (str)
# - config: standard Snakemake config dict (clair3_model: a name or a
#   {sample: name} mapping; clair3_model_dir)
# - LOG / BENCH / cpus / ram_mb / clair3_model_for / clair3_model_files from
#   rules/common.smk

# Clair3 2.x: models are two PyTorch checkpoints (pileup.pt, full_alignment.pt)
# in <model_dir>/<model>/, declared as inputs below. run_clair3.sh checks only
# that the directory exists and a stale (TensorFlow-era) directory dies inside
# torch.load with a bare exit 2, hence the exit-code decoding.
import os

CLAIR3_MODEL_DIR = str(config["clair3_model_dir"])

# Per-sample status written before variant calling: ``status\tok`` or
# ``status\tno_mapped_reads`` (fewer primary mapped reads than
# minimum_mapped_reads). Read by summary.tsv; a no_mapped_reads sample gets an
# all-N consensus instead of a Clair3 crash on an (almost) empty BAM.
STATUS_FILE = (
    config['output'] + "assembly/status/{sample}" + (".{segment}" if SEGMENT_WILDCARD else "") + ".txt"
)


rule check_mapped_reads:
    conda:
        "../envs/alignment.yaml"
    input:
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index
    output:
        status = STATUS_FILE
    params:
        minimum = config["minimum_mapped_reads"]
    log:
        LOG("check_mapped_reads")
    benchmark:
        BENCH("check_mapped_reads")
    shell:
        """
        set -euo pipefail
        exec 2> {log}
        n=$(samtools view -c -F 260 {input.bam})
        status=ok
        if [ "$n" -lt {params.minimum} ]; then status=no_mapped_reads; fi
        mkdir -p "$(dirname {output.status})"
        printf 'status\t%s\nmapped_reads\t%s\nminimum_mapped_reads\t%s\n' "$status" "$n" {params.minimum} > {output.status}
        """



rule infer_consensus_sequence:
    conda:
        "../envs/clair3.yaml"
    input:
        # The depth-normalised BAM (rules/normalise.smk); check_mapped_reads
        # and the coverage statistics stay on the trimmed BAM.
        bam = rules.normalise_depth.output.bam,
        bam_index = rules.normalise_depth.output.bam_index,
        reference = REFERENCE,
        fai = REFERENCE + ".fai",
        status = rules.check_mapped_reads.output.status,
        model_files = clair3_model_files,
        script = os.path.join(PY, "filter_variants.py")
    output:
        vcf_raw = config['output'] + "assembly/" + SEGMENT_WILDCARD + "clair3/{sample}/{sample}.raw.vcf.gz",
        vcf_raw_index = config['output'] + "assembly/" + SEGMENT_WILDCARD + "clair3/{sample}/{sample}.raw.vcf.gz.tbi",
        vcf_norm = temp(config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.norm.vcf.gz"),
        vcf = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.vcf.gz",
        vcf_index = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.vcf.gz.tbi",
        vcf_pass = temp(config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.pass.vcf"),
        filter_tsv = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.variant_filter.tsv",
        low_cov_bed = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.low_cov.bed",
        variant_mask_bed = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.masked_variants.bed",
        mask_bed = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.mask.bed",
        consensus = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.consensus.fasta",
        model_txt = config['output'] + "assembly/" + SEGMENT_WILDCARD + "clair3/{sample}/model.txt"
    params:
        output_prefix_dir = config['output'] + "assembly/" + SEGMENT_WILDCARD + "clair3/{sample}",
        model_dir = CLAIR3_MODEL_DIR,
        minimum_depth = config["minimum_depth"],
        af_threshold = config["af_threshold"],
        chunk_size = config["chunk_size"],
        clair3_model = lambda wildcards: clair3_model_for(wildcards.sample),
        variant_quality = config["variant_quality"],
        frameshift_variant_quality = config["frameshift_variant_quality"],
        minimum_map_quality = config["minimum_map_quality"],
        variant_depth = config["variant_depth"]
    benchmark:
        BENCH("infer_consensus_sequence")
    log:
        LOG("infer_consensus_sequence")
    threads: cpus("infer_consensus_sequence")
    resources:
        mem_mb = ram_mb("infer_consensus_sequence")
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        model_path="{params.model_dir}/{params.clair3_model}"
        # One definition for both branches: the empty-sample branch runs it on
        # a header-only VCF so every sample gets the same verdict table.
        run_filter() {{
            python {input.script:q} --vcf {output.vcf_norm:q} \
                --variant-quality {params.variant_quality} \
                --frameshift-variant-quality {params.frameshift_variant_quality} \
                --variant-depth {params.variant_depth} --af-threshold {params.af_threshold} \
                --output-vcf {output.vcf_pass:q} --output-bed {output.variant_mask_bed:q} \
                --output-tsv {output.filter_tsv:q}
        }}

        if [ "$(awk -F'\t' '$1=="status"{{print $2}}' {input.status})" = no_mapped_reads ]; then
            n_mapped=$(awk -F'\t' '$1=="mapped_reads"{{print $2}}' {input.status})
            echo "WARNING: {wildcards.sample}: $n_mapped primary mapped reads (< $(awk -F'\t' '$1=="minimum_mapped_reads"{{print $2}}' {input.status})); writing an all-N consensus and empty variant calls instead of running Clair3" >&2
            mkdir -p {params.output_prefix_dir} "$(dirname {output.consensus})"
            # One all-N record per reference contig, lengths from the .fai.
            : > {output.consensus}
            while IFS=$'\t' read -r name len _rest; do
                printf '>%s\n' "$name" >> {output.consensus}
                head -c "$len" /dev/zero | tr '\\0' 'N' >> {output.consensus}
                printf '\n' >> {output.consensus}
            done < {input.fai}
            write_empty_vcf() {{
                {{
                    echo "##fileformat=VCFv4.2"
                    awk -F'\t' '{{printf "##contig=<ID=%s,length=%s>\\n", $1, $2}}' {input.fai}
                    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                }} | bgzip -c > "$1"
                tabix -f -p vcf "$1"
            }}
            write_empty_vcf {output.vcf_raw}
            write_empty_vcf {output.vcf_norm}
            write_empty_vcf {output.vcf}
            # The filter on a header-only VCF writes header-only outputs, so the
            # verdict table has the same columns as every other sample's.
            run_filter
            : > {output.low_cov_bed}
            : > {output.mask_bed}
            printf 'model\t%s\nmodel_dir\t%s\nstatus\tno_mapped_reads\n' {params.clair3_model:q} "{params.model_dir}" > {output.model_txt}
            exit 0
        fi

        # Value-taking flags first, bare boolean flags last with nothing after
        # them: Clair3 2.x parses booleans as nargs="?" and would swallow a
        # following token as their value.
        rc=0
        run_clair3.sh \
            --bam_fn={input.bam} \
            --ref_fn={input.reference} \
            --model_path="$model_path" \
            --output={params.output_prefix_dir} \
            --platform=ont \
            --threads={threads} \
            --chunk_size={params.chunk_size} \
            --qual={params.variant_quality} \
            --min_mq={params.minimum_map_quality} \
            --var_pct_full=1.0 \
            --enable_long_indel \
            --haploid_sensitive \
            --no_phasing_for_fa \
            --include_all_ctgs || rc=$?
        if [ "$rc" -ne 0 ]; then
            echo "ERROR: run_clair3.sh exited $rc for {wildcards.sample} (model $model_path)" >&2
            case "$rc" in
                1) echo "  exit 1: argument or environment error - see {params.output_prefix_dir}/run_clair3.log" >&2 ;;
                2) echo "  exit 2: the model checkpoints failed to load (torch.load); expected $model_path/pileup.pt and full_alignment.pt from Clair3 2.x. Fetch them with 'viralconseq setup --clair3-models {params.clair3_model}'." >&2 ;;
                127) echo "  exit 127: run_clair3.sh is not on PATH - the clair3 env is broken; rebuild it with 'viralconseq setup'." >&2 ;;
            esac
            exit "$rc"
        fi
        printf 'model\t%s\nmodel_dir\t%s\n' {params.clair3_model:q} "{params.model_dir}" > {output.model_txt}

        cp {params.output_prefix_dir}/merge_output.vcf.gz {output.vcf_raw}
        cp {params.output_prefix_dir}/merge_output.vcf.gz.tbi {output.vcf_raw_index}

        bcftools norm -m - -f {input.reference} {output.vcf_raw} > {output.vcf_norm}
        # Pass / mask / discard (scripts/python/filter_variants.py), first match wins:
        #   discard  RefCall, or allele ratio below 0.5: the reference allele is
        #            the majority, keep it;
        #   mask     allele ratio between 0.5 and af_threshold (neither allele
        #            clears the bar), ALT reads below variant_depth, or QUAL below
        #            variant_quality (frameshift_variant_quality for an indel whose
        #            length is not a multiple of three): Clair3 saw something and
        #            the filter is not sure, so the position becomes N instead of
        #            a confident reference base;
        #   pass     written into the consensus.
        # The allele ratio is ALT / (REF + ALT) from FORMAT/AD, not Clair3's
        # FORMAT/AF (ALT / DP): DP counts reads that carry a deletion at the site,
        # and at a variant that creates a homopolymer half the ONT reads do, which
        # pushed true variants under the threshold and reverted them to the
        # reference. The denominator drops every read showing neither allele, so
        # the ratio is deliberately not relative to depth and variant_depth is the
        # only absolute floor left on the ALT read count. QUAL is read directly;
        # Clair3's LowQual stamp (--qual) only matters in raw.vcf.gz.
        # bcftools consensus drops any passing record whose REF span touches a
        # --mask interval, so the script never masks inside a passing span.
        run_filter
        bgzip -c {output.vcf_pass:q} > {output.vcf:q}
        tabix {output.vcf:q}

        samtools depth -J -a {input.bam} | \
            awk '$3 < int({params.minimum_depth}) {{print $1 "\t" $2-1 "\t" $2}}' > {output.low_cov_bed}
        # mask.bed is what bcftools consensus applies: positions below the depth
        # gate plus the rejected calls; the two parts are kept so each N can be
        # attributed.
        cat {output.low_cov_bed:q} {output.variant_mask_bed:q} | sort -k1,1 -k2,2n -k3,3n -u > {output.mask_bed:q}

        bcftools consensus -f {input.reference:q} --mask {output.mask_bed:q} {output.vcf:q} > {output.consensus:q}

        find {params.output_prefix_dir} -mindepth 1 ! -name '*.raw.vcf.gz' ! -name '*.raw.vcf.gz.tbi' ! -name 'model.txt' -exec rm -rf {{}} +
        """
