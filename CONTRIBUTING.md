# Contributing to viralconseq

Thanks for your interest in improving viralconseq. This document covers
the basics of setting up a development environment, running the test
suite, keeping the code lint-clean, and getting changes merged.

## Setup

viralconseq targets Python 3.10 – 3.11. The top-level conda environment
provides Snakemake and the core runtime dependencies; the bioinformatics
tools themselves (minimap2, samtools, fastp, Clair3, …) live in per-rule
conda environments under `viralconseq/scripts/envs/` that Snakemake builds
on demand.

```bash
# Clone and create the conda environment
git clone https://github.com/filiperomero2/viralconseq.git
cd viralconseq
conda env create -f environment.yml
conda activate viralconseq

# Install viralconseq in editable mode with dev dependencies
pip install -e ".[dev]"

viralconseq --version
```

The `dev` extras install `black`, `ruff`, `mypy`, `pytest` and `pytest-cov`
(see `pyproject.toml`'s `[project.optional-dependencies]`).

## Running tests

```bash
# Full unit suite (unittest discover under the hood)
make test

# Snakemake dry-run integration tests. Four dryrun configs cover
# consensus_{illumina,nanopore}{,_segmented}.
make test-dryrun

# Opt-in end-to-end suite on real SARS-CoV-2 data (downloads data, runs the
# full pipelines; needs the per-rule envs built via `viralconseq setup`).
make test-empirical
```

If `pytest` fails inside the dryrun tests with
`AttributeError: module 'pulp' has no attribute 'list_solvers'`, the
local `pulp` is too new for Snakemake 7.32; pin it with
`pip install 'pulp<2.8'`.

## Linting, formatting and type checks

```bash
make lint        # black --check + ruff check
make typecheck   # mypy viralconseq/
make format      # black + ruff --fix
```

CI runs all three and will fail the PR on diffs or type errors. The
canonical configurations live in `pyproject.toml` (`[tool.black]`,
`[tool.ruff]`, `[tool.ruff.lint]`, `[tool.mypy]`).

## Changing the workflows

The Snakemake files under `viralconseq/scripts/` are the substantive
code. When you change rule wiring, run `make test-dryrun`: it catches
missing inputs, broken `expand` patterns and circular dependencies that
the Python suite cannot see. New per-sample outputs need an `ln -sf`
block in `organize_files` to appear under `samples/<sample>/`.

## Submitting changes

1. Branch from `main`. Keep each commit focused on one logical change so
   the history stays bisectable.
2. Open a pull request against `main`. If your change addresses an
   existing issue, link it in the PR description.
3. Make sure CI is green: lint, type check, unit tests and dryrun tests
   must all pass.
4. A maintainer will review and merge.

## Releasing

Cutting a versioned release (bumping `__version__`, tagging, publishing
to PyPI) is documented in [`RELEASING.md`](RELEASING.md).

## Code of conduct

Participation in this project is governed by our
[Code of Conduct](CODE_OF_CONDUCT.md) (Contributor Covenant 2.1).
