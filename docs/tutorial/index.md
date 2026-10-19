# Tutorial

This tutorial takes you from a fresh viralconseq install to finished consensus genomes on real example data. It assumes no prior experience with the tool. By the end, you will know how to point the pipeline at your own viral sequencing data and how to read the output.

viralconseq turns raw viral sequencing reads into **one consensus genome per sample**, given a reference genome for the virus you sequenced.

## Pick a workflow

| If you have…                                                                  | …use                                                        |
|-------------------------------------------------------------------------------|-------------------------------------------------------------|
| Illumina paired-end reads (amplicon or shotgun) and a reference genome        | `viralconseq consensus illumina`                            |
| Nanopore reads and a reference genome                                          | `viralconseq consensus nanopore`                            |
| A multi-segment virus (influenza, bunyaviruses, …)                             | either subcommand with `--reference multi.fasta` (auto-split) |
| A fragmented reference (several contigs) that should be treated as one genome  | either subcommand with `--reference contigs.fasta --single-reference` |

## How to read this tutorial

1. [Setup](setup.md) — install viralconseq, download the example data and generate sample sheets. **Do this once.**
2. [Consensus pipeline](consensus.md) — run the Illumina and Nanopore workflows on SARS-CoV-2 data, read the outputs, and learn the few settings worth changing.

```{toctree}
:maxdepth: 1
:hidden:

setup
consensus
```

## About the example data

The worked examples use a small SARS-CoV-2 dataset with two Illumina amplicon samples (`my_test_data/illumina_data/`) and two Nanopore samples (`my_test_data/nanopore_data/`), plus the ARTIC nCoV-2019 reference genome and primer scheme. [Setup](setup.md#3-download-the-example-data) shows how to download them. If you prefer your own SARS-CoV-2 FASTQs, only the paths and sample IDs in the commands change.

## Further reading

Once you are productive, the reference pages have the full detail:

- [Commands reference](../commands.md) — every CLI option, with defaults.
- [Output layout](../output.md) — every file the pipeline produces and the columns of `summary.tsv`.
- [Notes](../notes.md) — segmented viruses, header sanitisation, amplicon handling and the consensus QC step.
