# viralQC post-consensus genome quality control.
#
# Expects the following to be defined in the entry-point workflow *before*
# this file is included:
# - config: output, samples, data, viralqc_db (required when run_viralqc is
#   on), optional run_viralqc (default True), viralqc_extra_flags,
#   run_viralqc_cpus, run_viralqc_ram
# - SEGMENT_WILDCARD: "" (single reference) or "{segment}/" (segmented)
# - SEGMENTS: segmented workflows only, dict segment -> reference path
# - rules.rename_sequences (rules/stats.smk)
#
# One viralQC invocation per run. Every sample's renamed consensus FASTA (all
# segments in segmented runs) is merged into ONE input FASTA. Headers are
# normalised to ``sample-<id>[|<contig>][|<segment>]``: rename_sequences writes
# ``sample-<id>`` (or ``sample-<id>_<contig>`` for a multi-contig single
# reference); the contig suffix is re-joined with ``|`` and segmented runs
# append ``|<segment>``, so ``|`` is the only separator and a sample id can be
# matched unambiguously. viralQC itemizes headers internally and restores them
# in the ``seqName`` column, so the results table has exactly one row per input
# record. The reference is NOT included.
#
# Lenient failure: if ``vqc run`` fails, run_viralqc still succeeds. It writes a
# placeholder results.tsv (one row per input record, ``inputSequenceStatus``
# explaining the failure), records the outcome in viralqc_status.txt and the
# rule log, and prints a WARNING. Consensus outputs are never affected.

import os

from snakemake.exceptions import WorkflowError

RUN_VIRALQC = config.get("run_viralqc", True)

# Snakemake-facing paths are built as plain strings on config["output"] (which
# may be relative in the dryrun configs) so that ``rule all`` targets and rule
# outputs match textually. Only the values handed to ``vqc`` are made absolute.
VIRALQC_DIR = config["output"] + "qc/viralqc/"
VIRALQC_DB = str(config.get("viralqc_db", "") or "")
VIRALQC_LOG_PREFIX = config["output"] + "logs/consensus_" + str(config.get("data", "run")) + "/"

if RUN_VIRALQC and VIRALQC_DB in ("", "NA"):
    raise WorkflowError(
        "run_viralqc is enabled but config key 'viralqc_db' is not set. "
        "Run 'viralconseq setup' to download the viralQC databases and point "
        "viralqc_db at them, or set run_viralqc: false."
    )

if SEGMENT_WILDCARD:
    VIRALQC_CONSENSUS_INPUTS = expand(
        rules.rename_sequences.output.consensus_renamed,
        sample=config["samples"],
        segment=SEGMENTS.keys(),
    )
else:
    VIRALQC_CONSENSUS_INPUTS = expand(
        rules.rename_sequences.output.consensus_renamed,
        sample=config["samples"],
    )


rule prepare_viralqc_input:
    conda:
        "../envs/utils.yaml"
    input:
        consensus_files = VIRALQC_CONSENSUS_INPUTS
    output:
        fasta = VIRALQC_DIR + "input.fasta"
    params:
        segmented = "1" if SEGMENT_WILDCARD else ""
    log:
        VIRALQC_LOG_PREFIX + "prepare_viralqc_input/prepare_viralqc_input.log"
    shell:
        """
        set -euo pipefail
        exec 2> {log}
        : > {output.fasta}
        for _file in {input.consensus_files}; do
            # .../assembly/[<segment>/]consensus/final_consensus/<sample>.consensus.renamed.fasta
            sample=$(basename $_file .consensus.renamed.fasta)
            segment=""
            if [ -n "{params.segmented}" ]; then
                segment=$(basename $(dirname $(dirname $(dirname $_file))))
            fi
            # Header -> sample[|contig][|segment]; contig = what rename_sequences
            # appended after "<sample>_" for multi-contig references.
            awk -v s="$sample" -v seg="$segment" '
                /^>/ {{
                    h = substr($1, 2); out = s
                    if (h != s && substr(h, 1, length(s) + 1) == s "_") out = out "|" substr(h, length(s) + 2)
                    else if (h != s) out = out "|" h
                    if (seg != "") out = out "|" seg
                    print ">" out; next
                }}
                {{print}}' $_file >> {output.fasta}
        done
        n_records=$(grep -c '^>' {output.fasta} || true)
        n_unique=$( (grep '^>' {output.fasta} || true) | sort -u | wc -l)
        echo "prepared $n_records record(s) for viralQC" >&2
        if [ "$n_records" -eq 0 ]; then
            echo "ERROR: no consensus records to hand to viralQC" >&2
            exit 1
        fi
        if [ "$n_records" -ne "$n_unique" ]; then
            echo "ERROR: duplicate FASTA headers in {output.fasta}" >&2
            exit 1
        fi
        """


rule run_viralqc:
    conda:
        "../envs/viralqc.yaml"
    input:
        fasta = rules.prepare_viralqc_input.output.fasta,
        blast_fasta = VIRALQC_DB + "/blast.fasta",
        blast_tsv = VIRALQC_DB + "/blast.tsv",
        nextclade_ok = VIRALQC_DB + "/.nextclade_datasets_ok",
    output:
        results = VIRALQC_DIR + "outputs/results.tsv",
        status = VIRALQC_DIR + "viralqc_status.txt",
    params:
        # vqc passes --output-dir verbatim to ``snakemake --directory`` and
        # runs from the job cwd, so hand it absolute paths.
        workdir = os.path.abspath(VIRALQC_DIR),
        db = os.path.abspath(VIRALQC_DB) if VIRALQC_DB not in ("", "NA") else VIRALQC_DB,
        extra_flags = config.get("viralqc_extra_flags", "") or "",
    threads: config.get("run_viralqc_cpus", 2)
    resources:
        mem_mb = int(config.get("run_viralqc_ram", 4)) * 1024
    log:
        VIRALQC_LOG_PREFIX + "run_viralqc/run_viralqc.log"
    benchmark:
        VIRALQC_LOG_PREFIX + "run_viralqc/run_viralqc.benchmark.txt"
    shell:
        """
        set -euo pipefail
        # Fresh nested-Snakemake working directory. Snakemake only removes the
        # declared outputs before a rerun; a stale outputs/ or .snakemake/ left
        # by an earlier attempt (e.g. the placeholder written below) would let
        # the inner Snakemake report "Nothing to be done" and exit 0.
        rm -rf {params.workdir}/outputs {params.workdir}/.snakemake
        mkdir -p {params.workdir}/outputs

        exit_code=0
        env -u SNAKEMAKE_PROFILE vqc run \
            --input {input.fasta} \
            --output-dir {params.workdir} \
            --output-file results.tsv \
            --datasets-dir {params.db} \
            --blast-database {input.blast_fasta} \
            --blast-database-metadata {input.blast_tsv} \
            --cores {threads} \
            --verbose {params.extra_flags} > {log} 2>&1 || exit_code=$?

        n_records=$(grep -c '^>' {input.fasta} || true)
        if [ -s {output.results} ]; then
            n_rows=$(($(wc -l < {output.results}) - 1))
            if [ "$exit_code" -eq 0 ]; then status=ok; else status=partial; fi
        else
            status=failed
            n_rows=0
            [ "$exit_code" -eq 0 ] && exit_code=99
            mkdir -p {params.workdir}/outputs
            {{
                printf 'seqName\\tgenomeQuality\\tinputSequenceStatus\\n'
                grep '^>' {input.fasta} | sed 's/^>//' | \
                    awk -v code="$exit_code" 'BEGIN{{OFS="\\t"}} {{print $0, "", "viralQC failed (exit " code ")"}}'
            }} > {output.results}
        fi

        printf 'status\\t%s\\nexit_code\\t%s\\ninput_records\\t%s\\nresult_rows\\t%s\\nresults\\t%s\\nlog\\t%s\\n' \
            "$status" "$exit_code" "$n_records" "$n_rows" {output.results} {log} > {output.status}

        if [ "$status" != ok ]; then
            msg="WARNING: viralQC $status (exit $exit_code). Consensus outputs are unaffected; see {output.status} and {log}. To retry after fixing the cause, delete {params.workdir} and rerun."
            echo "$msg" >> {log}
            echo "viralconseq: $msg" >&2
        fi
        """


rule split_viralqc_results:
    # Per-sample slice of results.tsv for the samples/<sample>/ browse tree.
    # Rows match on seqName == sample or seqName starting with "sample|"
    # (contig and/or segment suffixes). Works on the placeholder table too, so
    # organize_files never blocks on a viralQC failure.
    conda:
        "../envs/utils.yaml"
    input:
        results = rules.run_viralqc.output.results
    output:
        tsv = VIRALQC_DIR + "per_sample/{sample}.viralqc.tsv"
    shell:
        """
        set -euo pipefail
        awk -F'\\t' -v s="{wildcards.sample}" \
            'NR==1 || $1==s || index($1, s "|")==1' \
            {input.results} > {output.tsv}
        """
