# viralconseq

[![CI](https://github.com/filiperomero2/viralconseq/actions/workflows/ci.yaml/badge.svg)](https://github.com/filiperomero2/viralconseq/actions/workflows/ci.yaml)
[![PyPI](https://img.shields.io/pypi/v/viralconseq.svg)](https://pypi.org/project/viralconseq/)
[![Documentation](https://readthedocs.org/projects/viralconseq/badge/?version=latest)](https://viralconseq.readthedocs.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/filiperomero2/viralconseq/blob/main/LICENSE)

viralconseq infers consensus genome sequences from viral high-throughput sequencing data. It is a Python package that validates inputs, writes a Snakemake configuration and launches one of four reference-guided consensus workflows: Illumina paired-end or Nanopore reads, against a single or a segmented reference. viralconseq runs on *nix systems and processes entire sequencing runs in minimal time on a regular computer.

> **Full documentation:** <https://viralconseq.readthedocs.io/en/latest/>

## Features

- Reference-guided consensus assembly for **Illumina** (fastp → minimap2 → samtools/GSAlign) and **Nanopore** (minimap2 → Clair3 → bcftools) reads
- **Amplicon** primer clipping from a BED scheme, with primer/reference consistency checks
- **Segmented viruses**: a multi-record reference FASTA is split per segment and every segment is assembled in one run
- Optional **intra-host variant calling** (LoFreq) for Illumina data
- **Content-level input validation**: truncated FASTQs, protein FASTAs or mismatched primer schemes are rejected before the run starts
- Per-sample coverage statistics, a multi-sample alignment ready for phylogenetics, and a `run_manifest.json` with input checksums for reproducibility

## Installation

Install the CLI from PyPI:

```bash
pip install viralconseq
```

> **conda/mamba is required at runtime.** viralconseq orchestrates Snakemake, which builds
> the per-rule tool environments (aligner, variant callers, QC tools) via `--use-conda` on
> first run. Make sure conda or mamba is installed and on your `PATH`; pre-build those
> environments up front with `viralconseq setup --pipelines all`.

To install from source for development instead:

```bash
git clone https://github.com/filiperomero2/viralconseq.git
cd viralconseq
conda env create -n viralconseq -f environment.yml
conda activate viralconseq
pip install -e ".[dev]"
```

Per-rule conda environments under `viralconseq/scripts/envs/` are managed automatically by Snakemake; the top-level `environment.yml` only installs viralconseq itself and its core runtime dependencies.

## Quick start

```bash
viralconseq setup --pipelines all                                  # pre-build per-rule conda envs (once)
viralconseq create-samplesheet --input <run-dir> --output samples.csv
viralconseq consensus illumina --sample-sheet samples.csv --reference ref.fasta \
    --primer-scheme primers.bed --run-name run1 --config-file run1.yml --output results/
viralconseq consensus nanopore --sample-sheet samples.csv --reference ref.fasta \
    --run-name run1 --config-file run1.yml --output results/
```

Global options: `--log-level {DEBUG,INFO,WARNING,ERROR}` and `--json-logs`
(e.g. `viralconseq --log-level DEBUG consensus illumina ...`).

Each subcommand has its own `--help`; the same information is exhaustively documented in the `docs/` Sphinx site (rendered on ReadTheDocs at the link above), including a step-by-step [tutorial](https://viralconseq.readthedocs.io/en/latest/tutorial/index.html).

## Tests

```bash
make test          # unit suite (unittest)
make test-dryrun   # snakemake -n against every workflow (pytest)
make lint          # black + ruff
make typecheck     # mypy
```

## Citation

A scientific publication describing viralconseq is being prepared. Meanwhile, please cite this repository (see `CITATION.cff`). Primary references for the upstream tools (fastp, MultiQC, minimap2, SAMtools/BCFtools, BEDtools, LoFreq, Clair3, GSAlign, gofasta) are listed in the [documentation](https://viralconseq.readthedocs.io/en/latest/citation.html).

## Acknowledgements

viralconseq is the consensus-inference half of [ViralUnity](https://github.com/InstitutoTodosPelaSaude/ViralUnity), developed at the Instituto Todos pela Saúde together with Felippe Nacif. It was split into a standalone package (from ViralUnity v1.5.0) so that the consensus workflows can be validated, benchmarked and published independently of the metagenomics pipeline.

## License

MIT — see [`LICENSE`](https://github.com/filiperomero2/viralconseq/blob/main/LICENSE).
