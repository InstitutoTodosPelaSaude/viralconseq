# Tutorial

This tutorial walks you from a fresh viralconseq install to finished consensus genomes on real example data. It assumes no prior experience with the tool. By the end, you will know how to point the pipeline at your own viral high-throughput sequencing data and interpret the output.

viralconseq orchestrates a Snakemake workflow that turns raw viral sequencing reads into **a polished consensus genome per sample**, given a reference genome for the virus you sequenced.

## Pick a workflow

| If you have…                                                                  | …use                                                        |
|-------------------------------------------------------------------------------|-------------------------------------------------------------|
| Illumina paired-end reads (amplicon or shotgun) and a reference genome        | `viralconseq consensus illumina`                            |
| Nanopore reads and a reference genome                                          | `viralconseq consensus nanopore`                            |
| A multi-segment virus (influenza, bunyaviruses, …)                             | either subcommand with `--reference multi.fasta` (auto-split) |
| A fragmented reference (several contigs) that should be treated as one genome  | either subcommand with `--reference contigs.fasta --single-reference` |

## How to read this tutorial

1. [Setup](setup.md) — install viralconseq, build the per-rule environments, download the example data and generate sample sheets. **Do this once.**
2. [Consensus pipeline](consensus.md) — walk through reference-guided consensus assembly on Illumina and Nanopore SARS-CoV-2 data, plus how to handle segmented viruses.

```{toctree}
:maxdepth: 1
:hidden:

setup
consensus
```

## About the example data

The worked examples use a small SARS-CoV-2 dataset:

- `my_test_data/illumina_data/` — two paired-end Illumina amplicon samples: `itps-0001_R{1,2}.fastq.gz` and `itps-0002_R{1,2}.fastq.gz`.
- `my_test_data/nanopore_data/` — two Nanopore samples: `barcode05.itps-0003.fastq.gz` and `barcode09.itps-0004.fastq.gz`.

```{note}
[Setup §3](setup.md#3-download-the-example-data) walks you through downloading the tarball and the ARTIC reference / primer scheme. If you prefer to use your own SARS-CoV-2 FASTQs, only the paths and sample IDs in the commands need to change; a suitable reference genome is NCBI [`MN908947.3`](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3), available as `nCoV-2019.reference.fasta` in the public [ARTIC](https://github.com/artic-network/primer-schemes) primer-scheme repository.
```

## Further reading

The tutorial covers the common path through viralconseq. Once you are productive, the reference pages have the exhaustive detail:

- [Commands reference](../commands.md) — every CLI option, with defaults and units.
- [Output layout](../output.md) — the full directory tree the pipeline produces.
- [Notes](../notes.md) — segmented viruses, header sanitization and other advanced topics.
