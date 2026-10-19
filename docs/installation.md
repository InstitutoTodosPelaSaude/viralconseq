# Installation

## Requirements

- Linux or macOS (on Apple Silicon see the note under [Troubleshooting](#troubleshooting))
- Python 3.10 or 3.11
- conda, mamba or micromamba: the bioinformatics tools run in per-step conda environments that viralconseq builds for you
- [seqkit](https://bioinf.shenwei.me/seqkit/): used to validate every FASTQ before a run starts

The project's `environment.yml` provides Python, Snakemake, conda and seqkit in one environment, so the two routes below both start from it.

## Install from PyPI

```bash
conda env create -n viralconseq -f https://raw.githubusercontent.com/InstitutoTodosPelaSaude/viralconseq/main/environment.yml
conda activate viralconseq
pip install viralconseq
```

If you install with pip into an environment of your own instead, make sure conda and seqkit are on your `PATH` (`conda install -c bioconda seqkit`, or the static seqkit binary from its GitHub releases). Without seqkit, only `--skip-input-validation` runs.

## Install from source

`setup.sh` does everything in one command: it finds conda, mamba or micromamba, creates or updates the `viralconseq` environment from `environment.yml`, installs the package into it and runs the first-time setup described below. Re-running it is safe; what is already in place is skipped.

```bash
git clone https://github.com/InstitutoTodosPelaSaude/viralconseq.git
cd viralconseq
bash setup.sh                        # environment, package, tool environments, databases, Clair3 models
bash setup.sh --no-setup             # environment and package only
bash setup.sh --skip-clair3-models   # Illumina-only site
```

`bash setup.sh --help` lists the remaining options (`--conda-prefix`, `--viralqc-db`, `--clair3-model-dir`, `--clair3-models`, `--skip-viralqc-db`, `--cores`, `--env-name`).

The same steps by hand:

```bash
conda env create -n viralconseq -f environment.yml
conda activate viralconseq
pip install -e .
```

## First-time setup

After installing, run once:

```bash
viralconseq setup --pipelines all
```

This builds the per-step tool environments, downloads the viralQC databases (Nextclade datasets and the NCBI RefSeq viral BLAST set) and downloads the Clair3 models that `viralconseq consensus nanopore` needs. Expect about 1 GB on disk for the databases and 15–60 minutes in total, depending on your connection. Every later run reuses the three caches:

| Cache | Default location | Override |
|---|---|---|
| Tool environments | `~/.cache/viralconseq/conda-envs/` | `--conda-prefix`, `$VIRALCONSEQ_CONDA_PREFIX` |
| viralQC databases | `~/.cache/viralconseq/viralqc-db/` | `--viralqc-db`, `$VIRALCONSEQ_VIRALQC_DB` |
| Clair3 models | `~/.cache/viralconseq/clair3-models/` | `--clair3-model-dir`, `$VIRALCONSEQ_CLAIR3_MODELS` |

The environment variables are honoured by `viralconseq setup` and by every `viralconseq consensus` run, which is the simplest way to share the caches on a cluster. `viralconseq setup --pipelines consensus-illumina` builds only what the Illumina workflow needs; `--skip-viralqc-db` and `--skip-clair3-models` skip a download (run the pipeline with `--no-run-viralqc` if you skipped the databases). The [commands reference](commands.md#viralconseq-setup) lists every option.

On a node without internet access, run the setup elsewhere and copy the three cache directories over. The consensus workflows then run offline, but the viralQC step still contacts the Nextclade server for a small index on every run, so pass `--no-run-viralqc` there.

## Verify

```bash
viralconseq --version
viralconseq --help
```

## Troubleshooting

**Environment creation fails with a `repodata_shards.msgpack.zst` 404.** conda 26 queries a sharded repodata file that bioconda does not publish and mishandles the resulting 404 (`CreateCondaEnvironmentException` with a chained `HTTPError: 404` and `JSONDecodeError`). `environment.yml` pins `conda<26` and `conda-libmamba-solver<26` to avoid this, so rebuild the `viralconseq` environment from the pinned file:

```bash
conda deactivate
conda env remove -n viralconseq
conda env create -n viralconseq -f environment.yml
conda activate viralconseq
pip install viralconseq      # or pip install -e . from a checkout
```

If you cannot rebuild the environment (a shared cluster install, for example), fall back to the classic solver with `conda config --set solver classic`.

**macOS with Apple Silicon.** The Nanopore variant caller, Clair3 2.x, pulls in a CPU build of PyTorch (1–2 GB) and is published by bioconda for linux-64 and osx-64 only. Create the environment with `micromamba env create -f environment.yml --platform osx-64` and build the tool environments once with `viralconseq setup`; the Illumina workflow is unaffected.

**WSL.** Keep the caches on the Linux filesystem (the defaults under `~/.cache` already are), not under `/mnt/c`.

## Development install

To work on viralconseq itself, install the `dev` extras (linters, type checker and the test suites) and read [CONTRIBUTING.md](https://github.com/InstitutoTodosPelaSaude/viralconseq/blob/main/CONTRIBUTING.md):

```bash
pip install -e ".[dev]"
```
