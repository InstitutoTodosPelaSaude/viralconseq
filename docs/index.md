# viralconseq Documentation

**viralconseq** infers consensus genome sequences from viral high-throughput sequencing data. It is a Python package that validates inputs, writes a Snakemake configuration and launches one of four reference-guided consensus workflows (Illumina or Nanopore reads, single or segmented reference). All bioinformatics steps run inside per-rule conda environments managed by Snakemake.

viralconseq runs on *nix systems and processes entire sequencing runs in minimal time on a regular computer.

## Main features

- **Reference-guided consensus assembly** for Illumina paired-end and Nanopore long reads
- **Amplicon support** — primer clipping from a BED scheme, with primer/reference consistency checks
- **Segmented viruses** — a multi-record reference FASTA is split per segment and every segment is assembled in one run
- **Intra-host variant calling** (Illumina, optional) with LoFreq
- **Consensus QC** — every consensus sequence is assigned a virus and clade and scored A–D by [viralQC](https://github.com/InstitutoTodosPelaSaude/viralQC) (Nextclade + BLAST); opt out with `--no-run-viralqc`
- **Content-level input validation** — truncated FASTQs, protein FASTAs or mismatched primer schemes are rejected before the run starts
- **One place to look** — `summary.tsv` with a status per sample, a share-ready `consensus/` directory and a self-contained interactive `report.html`
- **Nanopore models resolved for you** — `--clair3-model auto` reads the basecalling model from the reads; models are cached by `viralconseq setup`; barcodes without mapped reads degrade to an all-N consensus instead of failing the run
- **Reproducibility** — every run writes `versions.tsv`, a copy of its config, `benchmark.tsv` and a `run_manifest.json` with input checksums; `viralconseq rerun <run>/config.yml` replays it; per-rule conda environments are pinned

## Documentation contents

```{toctree}
:maxdepth: 2

installation
tutorial/index
usage
commands
architecture
workflow
embedding
output
notes
citation
```

## Quick links

- [GitHub repository](https://github.com/filiperomero2/viralconseq)
- [Issues / bugs](https://github.com/filiperomero2/viralconseq/issues)
- **License:** MIT
