# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

viralconseq is a Python package whose only job at runtime is to validate inputs, write a YAML config, and launch one of four Snakemake workflows for reference-guided viral consensus inference. All of the actual bioinformatics — QC, alignment, primer clipping, variant and consensus calling, coverage statistics — lives in the Snakemake files under `viralconseq/scripts/`. Treat the Python layer as a CLI + orchestration shim; treat the `.smk` rule files as the substantive code.

viralconseq was extracted from the consensus half of [ViralUnity](https://github.com/InstitutoTodosPelaSaude/ViralUnity) v1.5.0. The retained rule modules, env YAMLs and helper scripts were copied byte-for-byte; keep it that way unless a change is deliberately about pipeline behaviour.

## Common commands

Setup is conda-based (the per-rule tool stacks live in `viralconseq/scripts/envs/*.yaml` and are pulled in by Snakemake's `--use-conda`):

```bash
conda env create -n viralconseq -f environment.yml
conda activate viralconseq
pip install -e ".[dev]"
```

Development loop:

```bash
make test            # python -m unittest discover ./test -p *test.py
make test-dryrun     # pytest test/dryrun_test.py -v   (snakemake -n on every workflow)
make lint            # black --check viralconseq/ test/  +  ruff check viralconseq/ test/
make typecheck       # mypy viralconseq/   (gating in CI)
make format          # black + ruff --fix
```

Run a single test file or test:

```bash
python -m unittest test.consensus_test -v
python -m unittest test.consensus_test.Test_ValidateArgs.test_validate_args_success
pytest test/dryrun_test.py -v -k consensus_illumina
```

Pipeline invocations (the user runs these against real data; for editing, prefer `--create-config-only` plus the dryrun tests to avoid long runs):

```bash
viralconseq create-samplesheet --input <runs-dir> --output samples.csv
viralconseq setup --pipelines all
viralconseq consensus illumina --sample-sheet ... --reference ... --config-file ... --output ...
viralconseq consensus nanopore --sample-sheet ... --reference ... --config-file ... --output ...
```

If `pytest` on the dryrun suite fails with `AttributeError: module 'pulp' has no attribute 'list_solvers'`, pin `pip install 'pulp<2.8'` — Snakemake 7.32 is incompatible with newer pulp.

## Architecture

### CLI → orchestrator → Snakemake

The single console entry point `viralconseq` (declared in `pyproject.toml`'s `[project.scripts]`) lives in `viralconseq/cli.py` and is a Click group with three subcommands: `consensus` (group with `illumina` / `nanopore`, options in `consensus_cli.py`), `setup` (`setup_cli.py`) and `create-samplesheet` (`create_samplesheet.py`). `consensus_cli.py` hands a plain `args` dict to `main()` in `viralconseq/consensus.py`.

`consensus.main()` is a thin wrapper around `_orchestrator.run_pipeline(...)` in `viralconseq/_orchestrator.py`. The orchestrator runs four callbacks in order: `resolve_paths` → `validate` → `generate_config` → `run_workflow_fn`, with shared error handling, a provenance manifest, and a `--create-config-only` short-circuit. `consensus.py` still owns `validate_args`, `generate_config_file`, and `run_snakemake_workflow` so existing test patches at those module-level names keep working — do not move those names into `_orchestrator` without updating the tests.

`run_snakemake_workflow` picks the workflow file by formatting a path: `viralconseq/scripts/consensus_{illumina,nanopore}[_segmented].smk`. Segmentation is selected by detecting that `args["reference"]` is a dict rather than a string. The dict is built either from repeatable `--segmented-reference SEGMENT=PATH` options, or by `viralconseq/reference_splitter.py` auto-splitting a single multi-record `--reference` FASTA during validation (segment keys = sanitized first token of each header; `--single-reference` opts out). All splitting happens in the Python layer, so the `.smk` files are unchanged — they always receive the per-segment dict.

### Config file is the contract

`viralconseq/config_generator.py` (`ConfigGenerator`) writes the YAML that Snakemake reads. The keys it emits are the same strings hard-coded in the `.smk` files (e.g. `config["samples"]`, `config["output"]`, `config["run_isnv"]`). If you add a new pipeline option, you touch four places: the Click option in `consensus_cli.py`, the `validators.py` checks, a `ConfigGenerator.add_*` setter, and the rule(s) that read it. Constants for these key names live in `viralconseq/constants.py` (`ConfigKeys`) but the `.smk` files reference the raw strings directly, so renaming a key means grepping `viralconseq/scripts/` too.

`ConfigGenerator.add_resource_settings(args, rule_list)` emits one `{rule}_cpus` / `{rule}_ram` pair per rule. The rule lists are declared as class attributes on `ResourceDefaults` in `constants.py` (`CONSENSUS_ILLUMINA_RULES`, `CONSENSUS_NANOPORE_RULES`). When you add a computationally heavy rule, add it to the right list so its resources land in the generated config.

One tool-level flag is config-only (no CLI), defaulted to the historical value for backwards compatibility: `minimap2_consensus_align_flags`. Users edit it by re-running with `--create-config-only`, editing the YAML, and rerunning Snakemake directly. Don't promote it to a CLI flag without discussion.

### Input validation

`validators.py` does existence and cross-dependency checks; `integrity.py` streams FASTQ / FASTA / BED content (pure stdlib) and collects `IntegrityIssue`s, which `validators.validate_consensus_input_integrity` aggregates into one `InputIntegrityError`. Integrity checks run even under `--create-config-only`; `--skip-input-validation` bypasses them.

### Snakemake workflows

Each top-level `.smk` in `viralconseq/scripts/` is small — it sets up wildcard helpers, `rule all`, and `include:`s rule modules from `viralconseq/scripts/rules/` (`qc_illumina`, `alignment_{illumina,nanopore}`, `consensus_{illumina,nanopore}`, `consensus_{illumina,nanopore}_common`, `stats`).

Cross-cutting conventions to know before editing rules:

- **Per-rule conda envs.** Every rule has `conda: "envs/<name>.yaml"` (relative to the workflow file; `../envs/` from `rules/`). Adding a new tool means either reusing an env or adding a YAML there. Snakemake's `--use-conda` is enabled in `_orchestrator.run_workflow`. `setup_cli.py --dry-run` scans the `include:`/`conda:` directives to list the envs a workflow needs; the real build lets Snakemake create them via `conda_create_envs_only`.
- **Optional outputs are computed conditionally.** `rule all` appends `isnvs/isnvs_summary.tsv` only when `run_isnv` is set, and `organize_files` mirrors that with `expand(... if <flag> else [])`. Whenever you add an optional step, edit both places.
- **`organize_files` is the symlink terminus.** It is the last rule before `benchmark.tsv` and creates the per-sample `samples/<sample>/...` symlinks that users actually browse. New per-sample outputs need a `ln -sf` block there to be discoverable.
- **Reference sanitization (nanopore).** The nanopore workflow sanitizes reference FASTA headers (replacing `/ \ | , ~` and spaces with `_`) before use because Clair3 makes per-contig directories from the seq IDs. Don't bypass this; `integrity.sanitize_nanopore_contig` must stay in sync with the `sed` in `sanitize_reference`.
- **Multi-contig single reference.** `align_consensus_to_reference_genome` builds the cross-sample alignment one reference contig at a time (gofasta aborts on multi-contig references) and concatenates them into `samples_alignment.fasta`.

### Sample sheets

CSV with no header. Illumina has 3 columns (`sample_id,R1,R2`), Nanopore has 2 (`sample_id,fastq`). `create-samplesheet` builds them by scanning a run directory; the parser (`viralconseq/validators.py:validate_sample_sheet`) keys off the data type given on the CLI and rejects rows with the wrong column count.

Sample names are prefixed with `sample-` inside the generated YAML by `ConfigGenerator.add_samples`, and the `.smk` files refer to `sample-<id>` everywhere. The prefix is visible in the outputs: `samples/sample-<id>/`, `assembly/coverage_stats/sample-<id>.table_cov_basewise.txt`, the `sample_name` column of `assembly_stats_summary.csv` and the consensus FASTA headers. Only `benchmark.tsv` strips it. Tests, docs and config inspection should expect the prefixed form.

## Tests

The unittest suite under `test/` covers Python-layer behaviour (CLI parsing, validators, integrity checks, config generation, path resolution, provenance). The Snakemake dryrun suite (`test/dryrun_test.py`) runs `snakemake -n` against every workflow + a YAML in `test/dryrun_configs/`; the placeholder fixture `create_dryrun_placeholders.sh` writes empty input files so paths resolve. `test/empirical_test.py` (opt-in, `make test-empirical`) runs the real SARS-CoV-2 scenario in `test/empirical/scenarios/`. CI runs lint + mypy, the unit suite + dryruns (one `test` job per Python version), a Sphinx docs build, a Docker build smoke test and a conda-env smoke test.

When changing rule wiring, run `make test-dryrun` — it catches missing inputs, broken `expand` patterns, and circular dependencies that the Python suite cannot see.

## Editing notes

- Don't rename `validate_args`, `generate_config_file`, or `run_snakemake_workflow` in `viralconseq/consensus.py` — multiple tests patch those exact module-level names.
- The `viralconseq/scripts/python/*.py` scripts run inside Snakemake's `script:` directive, which injects a `snakemake` global. Ruff would otherwise flag F821; `pyproject.toml` already silences this for that path. Don't add `from snakemake import snakemake` to those files, and don't make them import from the `viralconseq` package (they run inside per-rule conda envs where it is not installed).
- `viralconseq` is a published console script (`pyproject.toml: [project.scripts]`). After a fresh checkout, `pip install -e .` is required before `viralconseq --help` works — `make install-dev` does this.
- `viralconseq/__init__.py` must stay import-free: `docs/conf.py` imports `__version__` from it on ReadTheDocs, where the runtime dependencies are not installed.
- Version is single-sourced from `viralconseq/__init__.py:__version__`. Bumping it for release also requires editing the `Dockerfile` LABEL and `CITATION.cff` — see `RELEASING.md`.
- `my_results/`, `my_test_data/`, `data/` and `results/` are local working directories, not test fixtures; they are gitignored. Use `test/dryrun_configs/` for additions to the regression suite.

## Versioning policy

Follow Semantic Versioning (`MAJOR.MINOR.PATCH`) with these conventions:

- **PATCH** (`x.y.Z`) — small code adjustments: bug fixes, minor tweaks, docs/CI/packaging changes that don't add user-facing features.
- **MINOR** (`x.Y.0`) — larger blocks of changes: a group of significant features or a substantial enhancement. Reset PATCH to 0.
- **MAJOR** (`X.0.0`) — reserved for exceptional circumstances and bumped **manually by the maintainer only**. Do not bump MAJOR on your own; if a change seems to warrant it, propose it and let the user decide.

When cutting any release, follow the full `RELEASING.md` flow (bump `__version__` + `Dockerfile` LABEL + `CITATION.cff`, update `CHANGELOG.md`, tag `vX.Y.Z`) — pushing the tag triggers the PyPI publish workflow.
