# Read-depth normalisation, shared by all four workflows. Included after the
# alignment module (it reads rules.trim_primer_sequences) and before the
# consensus module (infer_consensus_sequence reads rules.normalise_depth).
#
# Deep amplicon runs hand every read to the consensus caller for no gain in
# accuracy; Clair3 rescales its pileup tensor above 216x. normalise_depth.py
# keeps, at every reference position, at least min(normalise_depth, original
# depth) records, preserving the strand ratio and (with a scheme) the amplicon
# mix, and re-emits them in input order so the BAM stays coordinate-sorted and
# byte-reproducible. Only infer_consensus_sequence reads the normalised BAM;
# depth statistics, the nanopore mapped-read gate and the Illumina iSNV step
# stay on the trimmed BAM. `normalise_depth: 0` switches the step off: the
# trimmed BAM is copied and header-only tables are written, so the DAG is
# fixed either way (the trim_primer_sequences `scheme: NA` pattern).
#
# Expects from the including workflow: SEGMENT_WILDCARD, LOG/BENCH/cpus/PY
# (rules/common.smk) and rules.trim_primer_sequences.
import os

NORMALISED_DIR = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/normalised/"


rule normalise_depth:
    conda:
        "../envs/alignment.yaml"
    input:
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index,
        script = os.path.join(PY, "normalise_depth.py")
    output:
        bam = NORMALISED_DIR + "{sample}.sorted.bam",
        bam_index = NORMALISED_DIR + "{sample}.sorted.bam.bai",
        depth = NORMALISED_DIR + "{sample}.depth.tsv",
        amplicons = NORMALISED_DIR + "{sample}.amplicons.tsv",
        summary = NORMALISED_DIR + "{sample}.normalise.txt"
    params:
        target = config["normalise_depth"],
        scheme = config.get("scheme", "NA")
    log:
        LOG("normalise_depth")
    benchmark:
        BENCH("normalise_depth")
    threads: cpus("normalise_depth")
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        run_normalise() {{
            python {input.script:q} --bam "$1" --bai {input.bam_index:q} \
                --output-bam {output.bam:q} --output-bai {output.bam_index:q} \
                --target {params.target} --scheme {params.scheme:q} \
                --depth-tsv {output.depth:q} --amplicons-tsv {output.amplicons:q} \
                --summary {output.summary:q} --threads {threads}
        }}
        rc=0
        run_normalise {input.bam:q} || rc=$?
        case $rc in
            0) ;;
            3)
                # Defensive: trim_primer_sequences indexes its output and
                # samtools index refuses an unsorted BAM, so a trimmed BAM out
                # of coordinate order should never reach this rule. If one does
                # (a hand-made input), sort a copy and retry. --no-PG keeps the
                # output header, and so the BAM, reproducible.
                echo "normalise_depth: input is not coordinate-sorted; sorting a copy and retrying"
                scratch={output.bam:q}.sorted_input.bam
                samtools sort --no-PG -@ {threads} -o "$scratch" {input.bam:q}
                rc=0
                run_normalise "$scratch" || rc=$?
                rm -f "$scratch"
                if [ "$rc" -ne 0 ]; then echo "normalise_depth.py exited with status $rc on the sorted copy"; exit "$rc"; fi
                ;;
            *) echo "normalise_depth.py exited with status $rc"; exit $rc ;;
        esac
        if [ "{params.target}" -ne 0 ]; then
            samtools index {output.bam:q} {output.bam_index:q}
        fi
        """
