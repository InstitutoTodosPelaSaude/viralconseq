# Downloads the viralQC databases consumed by rules/viralqc.smk. Run only by
# ``viralconseq setup`` via the Snakemake Python API with
# ``config={"viralqc_db": <dir>}``; never part of a consensus run.
#
# Both downloads target the same directory: the nextclade datasets (one
# sub-directory per virus) and the BLAST reference set (blast.fasta + index,
# blast.tsv, blast_gff/). viralQC has no completion marker for the nextclade
# datasets, so this workflow writes ``.nextclade_datasets_ok`` itself; the
# Python layer (constants.ViralQCDatabase) treats that file, blast.fasta,
# blast.tsv and blast_gff/ as the database's required entries.

import os

VIRALQC_DB = os.path.abspath(config["viralqc_db"])


rule all:
    input:
        VIRALQC_DB + "/blast.fasta",
        VIRALQC_DB + "/blast.tsv",
        VIRALQC_DB + "/blast_gff",
        VIRALQC_DB + "/.nextclade_datasets_ok",


rule get_nextclade_datasets:
    conda:
        "envs/viralqc.yaml"
    output:
        touch(VIRALQC_DB + "/.nextclade_datasets_ok")
    params:
        db = VIRALQC_DB
    threads: workflow.cores
    # Dataset fetches from data.clades.nextstrain.org occasionally time out;
    # a retry resumes where the previous attempt stopped (datasets already on
    # disk are skipped by viralQC's own workflow).
    retries: 2
    log:
        VIRALQC_DB + "/logs/get_nextclade_datasets.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.db}
        # vqc runs a nested Snakemake in the *current* directory; work in a
        # scratch dir so no .snakemake/ is left behind in the user's cwd.
        scratch=$(mktemp -d)
        trap 'rm -rf "$scratch"' EXIT
        cd "$scratch"
        env -u SNAKEMAKE_PROFILE vqc get-nextclade-datasets \
            --datasets-dir {params.db} --cores {threads} --verbose > {log} 2>&1
        """


rule get_blast_database:
    conda:
        "envs/viralqc.yaml"
    output:
        fasta = VIRALQC_DB + "/blast.fasta",
        tsv = VIRALQC_DB + "/blast.tsv",
        gff_dir = directory(VIRALQC_DB + "/blast_gff"),
    params:
        db = VIRALQC_DB
    threads: workflow.cores
    # The NCBI download is large; retry transient failures (viralQC >= 1.2.1
    # adds its own retry loop as well).
    retries: 2
    log:
        VIRALQC_DB + "/logs/get_blast_database.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.db}
        scratch=$(mktemp -d)
        trap 'rm -rf "$scratch"' EXIT
        cd "$scratch"
        env -u SNAKEMAKE_PROFILE vqc get-blast-database \
            --output-dir {params.db} --cores {threads} --verbose > {log} 2>&1
        """
