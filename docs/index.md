# viralconseq Documentation

**viralconseq** infers consensus genome sequences from viral high-throughput sequencing data. It is a Python package that validates inputs, writes a Snakemake configuration and launches one of four reference-guided consensus workflows (Illumina or Nanopore reads, single or segmented reference). All bioinformatics steps run inside per-rule conda environments managed by Snakemake.

viralconseq runs on *nix systems and processes entire sequencing runs in minimal time on a regular computer.

## Main features

- **Reference-guided consensus assembly** for Illumina paired-end and Nanopore long reads
- **Amplicon support** — primer clipping from a BED scheme, with primer/reference consistency checks
- **Segmented viruses** — a multi-record reference FASTA is split per segment and every segment is assembled in one run
- **Intra-host variant calling** (Illumina, optional) with LoFreq
- **Content-level input validation** — truncated FASTQs, protein FASTAs or mismatched primer schemes are rejected before the run starts
- **Reproducibility** — every run writes a `run_manifest.json` with input checksums and the tool version; per-rule conda environments are pinned

## Documentation contents

```{toctree}
:maxdepth: 2

installation
tutorial/index
usage
commands
architecture
embedding
output
notes
citation
```

## Quick links

- [GitHub repository](https://github.com/filiperomero2/viralconseq)
- [Issues / bugs](https://github.com/filiperomero2/viralconseq/issues)
- **License:** MIT
