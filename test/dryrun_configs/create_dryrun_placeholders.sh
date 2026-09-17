#!/bin/bash
# Create the empty placeholder inputs referenced by test/dryrun_configs/*.yaml
# so `snakemake -n` can resolve every path. Run from the repository root.

mkdir -p data/reads
mkdir -p data/references

# Dummy FASTQ files
touch data/reads/SAMPLE1_R1.fastq.gz
touch data/reads/SAMPLE1_R2.fastq.gz
touch data/reads/SAMPLE2_R1.fastq.gz
touch data/reads/SAMPLE2_R2.fastq.gz
touch data/reads/SAMPLE_NP.fastq.gz
touch data/reads/SAMPLE_NP_SEG.fastq.gz

# Dummy reference files (single + per-segment)
touch data/references/reference.fasta
touch data/references/segment1.fasta
touch data/references/segment2.fasta
touch data/references/segment3.fasta

# Dummy primer scheme
touch data/references/primers.bed

# Dummy viralQC database (layout per constants.ViralQCDatabase)
mkdir -p data/viralqc_db/blast_gff
touch data/viralqc_db/blast.fasta data/viralqc_db/blast.tsv data/viralqc_db/.nextclade_datasets_ok

echo "Placeholder files created successfully in data/."
echo "You can now run 'snakemake -n' for any of the dry-run configs."

# Clair3 model checkpoints declared as inputs by rules/consensus_nanopore.smk
mkdir -p data/clair3_models/r1041_e82_400bps_sup_v500
touch data/clair3_models/r1041_e82_400bps_sup_v500/pileup.pt data/clair3_models/r1041_e82_400bps_sup_v500/full_alignment.pt
