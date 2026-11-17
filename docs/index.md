# viralconseq Documentation

**viralconseq** infers consensus genome sequences from viral high-throughput sequencing data. Give it a sample sheet, a reference genome and, for amplicon data, a primer scheme; it validates the inputs and runs a reference-guided Snakemake workflow for Illumina paired-end or Nanopore reads, against a single or a segmented reference. It runs on Linux and macOS and processes a whole sequencing run on a regular computer.

## Main features

- **Reference-guided consensus assembly** for Illumina paired-end and Nanopore reads
- **Amplicon support** — primer clipping from a BED scheme, with primer/reference consistency checks
- **Segmented viruses** — a multi-record reference FASTA is split per segment and every segment is assembled in one run
- **Consensus QC** — every consensus sequence is assigned a virus and clade and graded A–D by [viralQC](https://github.com/InstitutoTodosPelaSaude/viralQC)
- **Intra-host variant calling** (Illumina, optional) with LoFreq
- **Input validation before the run starts** — truncated FASTQs, protein FASTAs or mismatched primer schemes are rejected up front
- **Nanopore models resolved for you** — the Clair3 model is read from the basecalling tag in the reads; barcodes without mapped reads get an all-N consensus instead of failing the run
- **One place to look** — `summary.tsv` with a status per sample, a share-ready `consensus/` directory and a self-contained interactive `report.html`
- **Reproducible** — every run records tool versions, its configuration, input checksums and per-rule benchmarks, and can be replayed with `viralconseq rerun`

## Where to start

New to viralconseq? Install it, then follow the [tutorial](tutorial/index.md): it takes you from a fresh install to finished consensus genomes on example data. The reference pages describe every command, option and output file. The developer pages describe how the package and the workflows are put together.

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
tutorial/index
```

```{toctree}
:maxdepth: 2
:caption: Reference

usage
commands
output
notes
```

```{toctree}
:maxdepth: 1
:caption: For developers

architecture
workflow
embedding
```

```{toctree}
:maxdepth: 1
:caption: About

citation
```

## Quick links

- [GitHub repository](https://github.com/InstitutoTodosPelaSaude/viralconseq)
- [Issues / bugs](https://github.com/InstitutoTodosPelaSaude/viralconseq/issues)
- **License:** MIT
