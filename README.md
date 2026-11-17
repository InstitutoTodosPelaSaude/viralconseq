# viralconseq

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CI](https://github.com/InstitutoTodosPelaSaude/viralconseq/actions/workflows/ci.yaml/badge.svg)](https://github.com/InstitutoTodosPelaSaude/viralconseq/actions/workflows/ci.yaml)
[![PyPI](https://img.shields.io/pypi/v/viralconseq.svg)](https://pypi.org/project/viralconseq/)
[![Documentation](https://readthedocs.org/projects/viralconseq/badge/?version=latest)](https://viralconseq.readthedocs.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/InstitutoTodosPelaSaude/viralconseq/blob/main/LICENSE)

viralconseq infers consensus genome sequences from viral high-throughput sequencing data. Give it a sample sheet, a reference genome and, for amplicon data, a primer scheme; it validates the inputs and runs a reference-guided Snakemake workflow for Illumina paired-end or Nanopore reads, against a single or a segmented reference. It runs on Linux and macOS and processes a whole sequencing run on a regular computer.

> **Full documentation:** <https://viralconseq.readthedocs.io/en/latest/>

## Features

- Reference-guided consensus assembly for **Illumina** paired-end and **Nanopore** reads
- **Amplicon** primer clipping from a BED scheme, with primer/reference consistency checks
- **Segmented viruses**: a multi-record reference FASTA is split per segment and every segment is assembled in one run
- **Consensus QC with [viralQC](https://github.com/InstitutoTodosPelaSaude/viralQC)**: virus and clade assignment plus a genome-quality grade (A–D) for every consensus sequence
- Optional **intra-host variant calling** (LoFreq) for Illumina data
- **Input validation before the run starts**: truncated FASTQs, protein FASTAs or mismatched primer schemes are rejected up front
- **Nanopore models resolved for you**: the Clair3 model is read from the basecalling tag in the reads; barcodes without mapped reads get an all-N consensus instead of failing the run
- **One place to look**: `summary.tsv` with a status per sample, a share-ready `consensus/` directory and a self-contained interactive `report.html`
- **Reproducible**: every run records tool versions, its configuration, input checksums and per-rule benchmarks, and can be replayed with `viralconseq rerun`

## Installation

viralconseq needs Python 3.10 or 3.11, conda (or mamba/micromamba) and seqkit. The project's `environment.yml` provides all of them:

```bash
conda env create -n viralconseq -f https://raw.githubusercontent.com/InstitutoTodosPelaSaude/viralconseq/main/environment.yml
conda activate viralconseq
pip install viralconseq
viralconseq setup --pipelines all      # once: tool environments, viralQC databases, Clair3 models
```

See the [installation page](https://viralconseq.readthedocs.io/en/latest/installation.html) for installing from source (`bash setup.sh`), where the caches live and troubleshooting.

## Quick start

```bash
viralconseq create-samplesheet --input <run-dir> --output samples.csv
viralconseq consensus illumina --sample-sheet samples.csv --reference ref.fasta \
    --primer-scheme primers.bed --run-name run1 --config-file run1.yml --output results/
viralconseq consensus nanopore --sample-sheet samples.csv --reference ref.fasta \
    --run-name run1 --config-file run1.yml --output results/
viralconseq rerun results/run1/config.yml --dry-run                # replay a run from its config
viralconseq create-report results/run1                             # rebuild report.html
```

Every subcommand has a `--help`. The [tutorial](https://viralconseq.readthedocs.io/en/latest/tutorial/index.html) walks through a complete run on example data, and the [commands reference](https://viralconseq.readthedocs.io/en/latest/commands.html) lists every option.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for the development setup, the test suites and the release flow.

## Citation

A scientific publication describing viralconseq is being prepared. Meanwhile, please cite this repository (see `CITATION.cff`). References for the upstream tools (fastp, MultiQC, minimap2, SAMtools/BCFtools, BEDtools, LoFreq, Clair3, GSAlign, gofasta, viralQC, Nextclade, BLAST+) are listed in the [documentation](https://viralconseq.readthedocs.io/en/latest/citation.html).

## License

MIT — see [`LICENSE`](https://github.com/InstitutoTodosPelaSaude/viralconseq/blob/main/LICENSE).
