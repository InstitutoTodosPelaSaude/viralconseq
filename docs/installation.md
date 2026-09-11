# Installation

viralconseq is a Python package that launches Snakemake workflows. The core runtime dependencies are listed in the conda environment file `environment.yml`; the per-rule bioinformatics tools (minimap2, samtools, fastp, Clair3, …) live in `viralconseq/scripts/envs/` and are built by Snakemake on demand.

## Install from PyPI (recommended)

The viralconseq CLI is published on PyPI, so the quickest install is:

```bash
pip install viralconseq
```

```{note}
conda (or mamba/micromamba) is still required **at runtime**. viralconseq orchestrates
Snakemake, which builds the per-rule tool environments (aligner, variant callers, QC tools)
via `--use-conda` on first run — so conda/mamba must be installed and on your `PATH`. Pre-build
those environments once with `viralconseq setup --pipelines all` (see [First-time environment
setup](#first-time-environment-setup) below).
```

## Clone and create the environment

To work from a source checkout instead (required for development):

```bash
git clone https://github.com/filiperomero2/viralconseq.git
cd viralconseq/
conda env create -f environment.yml
conda activate viralconseq
```

Or with **micromamba** (recommended on macOS with Apple Silicon):

```bash
micromamba env create -f environment.yml --platform osx-64
micromamba activate viralconseq
```

```{warning}
On macOS with Apple Silicon (M1 or later), the `viralconseq/scripts/envs/clair3.yaml` environment may not install correctly due to compatibility constraints in the clair3 dependencies.
```

## Troubleshooting

### `CreateCondaEnvironmentException` with a `repodata_shards.msgpack.zst` 404

If the first pipeline run aborts on `Creating conda environment .../qc.yaml ...` with a chained `HTTPError: 404` for `repodata_shards.msgpack.zst` followed by a `JSONDecodeError`, you are hitting a known incompatibility between conda 26.x's libmamba shards path and bioconda (which does not publish shards). `environment.yml` already pins `conda<26` and `conda-libmamba-solver<26` to avoid this, so the fix is to **rebuild the `viralconseq` env from the pinned file**:

```bash
conda deactivate
conda env remove -n viralconseq
conda env create -n viralconseq -f environment.yml
conda activate viralconseq
pip install -e .
```

If you cannot rebuild the env (for example, on a shared cluster install), the manual escape hatch is to fall back to the classic solver:

```bash
conda config --set solver classic
```

This is why `environment.yml` pins `conda<26` / `conda-libmamba-solver<26`: newer
solvers query `repodata_shards.msgpack.zst` from bioconda (not published) and
mishandle the resulting 404. If you hit environment-creation errors, set the
classic solver as shown above.

## Verify the installation

```bash
viralconseq --version
viralconseq --help
```

## First-time environment setup

After installing viralconseq, build the per-rule conda environments once with:

```bash
viralconseq setup --pipelines all
```

This downloads and resolves every workflow dependency into `~/.cache/viralconseq/conda-envs/` (override with `--conda-prefix PATH`, or set `$VIRALCONSEQ_CONDA_PREFIX`). Future `viralconseq consensus` runs reuse the cached envs and do not need to re-create them per working directory, which both speeds up first runs and isolates env-creation failures from real pipeline runs.

Run `viralconseq setup --pipelines consensus-illumina --dry-run` first to inspect what would be built.

## Development install

To work on viralconseq itself, install the optional `dev` extras (linters and tests):

```bash
pip install -e ".[dev]"
```

See [CONTRIBUTING.md](https://github.com/filiperomero2/viralconseq/blob/main/CONTRIBUTING.md) for the full development workflow.
