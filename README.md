# viralconseq

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![CI](https://github.com/InstitutoTodosPelaSaude/viralconseq/actions/workflows/ci.yaml/badge.svg)](https://github.com/InstitutoTodosPelaSaude/viralconseq/actions/workflows/ci.yaml)
[![PyPI](https://img.shields.io/pypi/v/viralconseq.svg)](https://pypi.org/project/viralconseq/)
[![Documentation](https://readthedocs.org/projects/viralconseq/badge/?version=latest)](https://viralconseq.readthedocs.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/InstitutoTodosPelaSaude/viralconseq/blob/main/LICENSE)

viralconseq infers consensus genome sequences from viral high-throughput sequencing data. It is a Python package that validates inputs, writes a Snakemake configuration and launches one of four reference-guided consensus workflows: Illumina paired-end or Nanopore reads, against a single or a segmented reference. viralconseq runs on *nix systems and processes entire sequencing runs in minimal time on a regular computer.

> **Full documentation:** <https://viralconseq.readthedocs.io/en/latest/>

## Features

- Reference-guided consensus assembly for **Illumina** (fastp → minimap2 → samtools/GSAlign) and **Nanopore** (minimap2 → Clair3 → bcftools) reads
- **Amplicon** primer clipping from a BED scheme, with primer/reference consistency checks
- **Segmented viruses**: a multi-record reference FASTA is split per segment and every segment is assembled in one run
- Optional **intra-host variant calling** (LoFreq) for Illumina data
- **Consensus QC with [viralQC](https://github.com/InstitutoTodosPelaSaude/viralQC)**: virus and clade assignment plus a genome-quality score (A–D) for every consensus sequence (Nextclade + BLAST), on by default
- **Content-level input validation**: truncated FASTQs, protein FASTAs or mismatched primer schemes are rejected before the run starts
- **Nanopore models resolved for you**: `--clair3-model auto` reads the basecalling model from the reads; models are cached by `viralconseq setup`; barcodes without mapped reads degrade to an all-N consensus instead of failing the run
- **Resources that bind**: detected cores and a memory budget are passed to Snakemake by default (`--threads-total`, `--max-memory`)
- **One place to look**: `summary.tsv` with a status column per sample, a share-ready `consensus/` directory and a self-contained interactive `report.html`
- Per-sample coverage statistics, a multi-sample alignment ready for phylogenetics, `versions.tsv`, and a `run_manifest.json` with input checksums for reproducibility

## Installation

viralconseq is published on [PyPI](https://pypi.org/project/viralconseq/). It needs Python 3.10
or 3.11 and conda at runtime (Snakemake builds the per-rule tool environments with
`--use-conda`), so install it into a conda environment created from the project's
`environment.yml`:

```bash
conda env create -n viralconseq -f https://raw.githubusercontent.com/InstitutoTodosPelaSaude/viralconseq/main/environment.yml
conda activate viralconseq
pip install viralconseq
```

Then run `viralconseq setup --pipelines all` once, as in the Quick start below: it pre-builds those
environments and downloads the viralQC databases (about 1 GB).

To install from source instead, `setup.sh` does the whole thing in one command: it
finds conda/mamba/micromamba (even when they are not on `PATH`), creates or updates the
`viralconseq` environment from `environment.yml`, installs the package into it and then runs
`viralconseq setup --pipelines all` (per-rule envs, viralQC databases, Clair3 models):

```bash
git clone https://github.com/InstitutoTodosPelaSaude/viralconseq.git
cd viralconseq
bash setup.sh                        # everything; re-running is a no-op for what is already there
bash setup.sh --no-setup             # environment + package only
bash setup.sh --skip-clair3-models   # Illumina-only site: no Clair3 models
```

`bash setup.sh --help` lists the pass-through options (`--conda-prefix`, `--viralqc-db`,
`--clair3-model-dir`, `--clair3-models`, `--skip-viralqc-db`, `--cores`, `--env-name`). Or do
the same steps by hand:

```bash
conda env create -n viralconseq -f environment.yml
conda activate viralconseq
pip install -e ".[dev]"
```

Per-rule conda environments under `viralconseq/scripts/envs/` are managed automatically by Snakemake; the top-level `environment.yml` only holds the core runtime dependencies (Python, Snakemake, conda).

## Quick start

```bash
viralconseq setup --pipelines all                                  # pre-build per-rule conda envs + download viralQC databases (once)
viralconseq create-samplesheet --input <run-dir> --output samples.csv
viralconseq consensus illumina --sample-sheet samples.csv --reference ref.fasta \
    --primer-scheme primers.bed --run-name run1 --config-file run1.yml --output results/
viralconseq consensus nanopore --sample-sheet samples.csv --reference ref.fasta \
    --run-name run1 --config-file run1.yml --output results/
viralconseq rerun results/run1/config.yml --dry-run                # replay a run from its config
viralconseq create-report results/run1                             # rebuild report.html
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

A scientific publication describing viralconseq is being prepared. Meanwhile, please cite this repository (see `CITATION.cff`). Primary references for the upstream tools (fastp, MultiQC, minimap2, SAMtools/BCFtools, BEDtools, LoFreq, Clair3, GSAlign, gofasta, viralQC, Nextclade, BLAST+) are listed in the [documentation](https://viralconseq.readthedocs.io/en/latest/citation.html).

## Acknowledgements

viralconseq is the consensus-inference half of [ViralUnity](https://github.com/InstitutoTodosPelaSaude/ViralUnity), developed at the Instituto Todos pela Saúde together with Felippe Nacif. It was split into a standalone package (from ViralUnity v1.5.0) so that the consensus workflows can be validated, benchmarked and published independently of the metagenomics pipeline.

## License

MIT — see [`LICENSE`](https://github.com/InstitutoTodosPelaSaude/viralconseq/blob/main/LICENSE).
