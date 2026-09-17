# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

viralconseq is a Python package whose only job at runtime is to validate inputs, write a YAML config, and launch one of four Snakemake workflows for reference-guided viral consensus inference. All of the actual bioinformatics — QC, alignment, primer clipping, variant and consensus calling, coverage statistics, and the final consensus QC with viralQC — lives in the Snakemake files under `viralconseq/scripts/`. Treat the Python layer as a CLI + orchestration shim; treat the `.smk` rule files as the substantive code.

viralconseq was extracted from the consensus half of [ViralUnity](https://github.com/InstitutoTodosPelaSaude/ViralUnity) v1.5.0. The rule modules, env YAMLs and helper scripts started as byte-for-byte copies; since 0.1.1 they diverge freely, and every divergence that changes behaviour or outputs gets a `CHANGELOG.md` line.

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
bash setup.sh                              # conda env + pip install -e + viralconseq setup, idempotent
viralconseq create-samplesheet --input <runs-dir> --output samples.csv [--prefix RUN1]
viralconseq setup --pipelines all          # per-rule envs, viralQC databases, default Clair3 models
viralconseq consensus illumina --sample-sheet ... --reference ... --config-file ... --output ...
viralconseq consensus nanopore --sample-sheet ... --reference ... --config-file ... --output ...
viralconseq rerun <run>/config.yml --dry-run          # replay a run from its config (--unlock, --set KEY=VALUE)
viralconseq create-report <run>                       # rebuild report.html from a finished run directory
```

Real-data fixtures for the two workflows live in `my_test_data/` (untagged 2022 R9.4.1 Guppy nanopore reads, so `--clair3-model auto` fails there by design; pass `--clair3-model r941_prom_hac_g360+g422`). The scratch run directories go under `results/` (gitignored).

If `pytest` on the dryrun suite fails with `AttributeError: module 'pulp' has no attribute 'list_solvers'`, pin `pip install 'pulp<2.8'` — Snakemake 7.32 is incompatible with newer pulp.

## Architecture

### CLI → orchestrator → Snakemake

The single console entry point `viralconseq` (declared in `pyproject.toml`'s `[project.scripts]`) lives in `viralconseq/cli.py` and is a Click group with five subcommands: `consensus` (group with `illumina` / `nanopore`, options in `consensus_cli.py`), `rerun` (`rerun_cli.py`), `create-report` (`report_cli.py` over `report.py`), `setup` (`setup_cli.py`) and `create-samplesheet` (`create_samplesheet.py`). `consensus_cli.py` hands a plain `args` dict to `main()` in `viralconseq/consensus.py`.

`consensus.main()` is a thin wrapper around `_orchestrator.run_pipeline(...)` in `viralconseq/_orchestrator.py`. The orchestrator runs four callbacks in order: `resolve_paths` → `validate` → `generate_config` → `run_workflow_fn`, with shared error handling, a provenance manifest, and a `--create-config-only` short-circuit. `consensus.py` still owns `validate_args`, `generate_config_file`, and `run_snakemake_workflow` so existing test patches at those module-level names keep working — do not move those names into `_orchestrator` without updating the tests.

`run_snakemake_workflow` picks the workflow file by formatting a path: `viralconseq/scripts/consensus_{illumina,nanopore}[_segmented].smk`. Segmentation is selected by detecting that `args["reference"]` is a dict rather than a string. The dict is built either from repeatable `--segmented-reference SEGMENT=PATH` options, or by `viralconseq/reference_splitter.py` auto-splitting a single multi-record `--reference` FASTA during validation (segment keys = sanitized first token of each header; `--single-reference` opts out). All splitting happens in the Python layer, so the `.smk` files are unchanged — they always receive the per-segment dict.

### Config file is the contract

`viralconseq/config_generator.py` (`ConfigGenerator`) writes the YAML that Snakemake reads. The keys it emits are the same strings hard-coded in the `.smk` files (e.g. `config["samples"]`, `config["output"]`, `config["run_isnv"]`). If you add a new pipeline option, you touch four places: the Click option in `consensus_cli.py`, the `validators.py` checks, a `ConfigGenerator.add_*` setter (pick the right `SECTION_*` so the key lands under the right commented section), and the rule(s) that read it. Constants for these key names live in `viralconseq/constants.py` (`ConfigKeys`) but the `.smk` files reference the raw strings directly, so renaming a key means grepping `viralconseq/scripts/` too.

Analysis parameters are **required** keys: rules read `config["key"]` with no fallback (the CLI default is the only default), and `rules/common.smk` holds the `_REQUIRED_KEYS` list that refuses a YAML missing one at parse time. A new required key goes into that list and into the mirror list in `test/rule_inventory_test.py`. Only sentinel/optional keys (`scheme`, `adapters`, `run_isnv`, `run_viralqc`, the two `*_flags` strings, `*_cpus`/`*_ram`, `viralconseq_version`) are read with `config.get`.

Snakemake runs with `<output>/<run_name>/` as its working directory (`_orchestrator.run_workflow` passes `workdir=`), so every path in the YAML must be absolute: `resolve_path_args` handles the CLI path options and `validators.absolutise_sample_paths` the sample FASTQs. `.snakemake/` therefore lives inside the run directory, next to `logs/run.log` (written by the `onstart`/`onsuccess`/`onerror` hooks in `rules/common.smk`) and `logs/snakemake.log` (copied by `provenance.copy_snakemake_log`).

**Resources bind by default.** `ResourceDefaults` in `constants.py` declares the threaded rules (`CONSENSUS_ILLUMINA_RULES`, `CONSENSUS_NANOPORE_RULES`, each with a `--<rule>-cpus` option written as `<rule>_cpus` only when given) and `MEMORY_RULES` (`infer_consensus_sequence` 2 GB on nanopore, `run_viralqc` 1 GB), the only rules with a `--<rule>-ram` option and a `mem_mb` declaration. Rules read them through `cpus("<own name>")` (`--<rule>-cpus` > `--threads` > 2) and `ram_mb("<own name>")` in `rules/common.smk`. `--threads-total` defaults to the cores available to the process minus one (`detect_cores`), and a memory budget (`MemTotal` capped by the cgroup limit, minus 10 %; `detect_memory_mb`/`memory_budget_mb`) is passed as `resources={"mem_mb": ...}`; `--max-memory 0` disables it, a positive value below the largest per-rule figure is refused, and `describe_resources` prints the resolved budget at start. Never add a `mem_mb` figure that was not read out of a run's `benchmark.tsv` (`max_rss`); an unmeasured declaration only throttles the run.

**Clair3 models live outside the conda env.** `~/.cache/viralconseq/clair3-models` (`$VIRALCONSEQ_CLAIR3_MODELS`, `--clair3-model-dir`) holds one `<name>/{pileup.pt,full_alignment.pt}` pair per model (Clair3 2.0.2, PyTorch checkpoints). Only `viralconseq setup --clair3-models` downloads them (`scripts/clair3_setup.smk` + `scripts/python/fetch_clair3_model.py`, ZIP-validated, atomic publish); analysis rules never fetch. `--clair3-model auto` (the default) reads the basecall tag from each sample's FASTQ headers (`clair3_models.py`, a port of artic's `choose_model` over `constants.Clair3Models.MANIFEST`) and resolves per sample; untagged reads fail pre-flight with a hint naming the flag, and the rule declares the checkpoints as inputs so a missing model is a `MissingInputException`, not a Clair3 crash. `_with_mv` (move-table) models are rejected by name.

Two tool-level settings are config-only (no CLI): `minimap2_consensus_align_flags` (historical default kept for backwards compatibility) and `viralqc_extra_flags` (extra `vqc run` arguments, default empty). Users edit them by re-running with `--create-config-only`, editing the YAML, and rerunning Snakemake directly. Don't promote them to CLI flags without discussion.

### Input validation

`validators.py` does existence and cross-dependency checks; `integrity.py` streams FASTQ / FASTA / BED content (pure stdlib) and collects `IntegrityIssue`s, which `validators.validate_consensus_input_integrity` aggregates into one `InputIntegrityError`. Integrity checks run even under `--create-config-only`; `--skip-input-validation` bypasses them. The viralQC database check (`validators.validate_viralqc_database`, error code `viralqc_database_not_found`) is an existence check and is bypassed only by `--no-run-viralqc`; the required directory layout is single-sourced in `constants.ViralQCDatabase` (including the `makeblastdb` index next to `blast.fasta`) and must stay in sync with `scripts/viralqc_setup.smk`'s `rule all`. That check is the authority: `run_viralqc` declares only the plain files as rule inputs (a `directory()` input would break `setup`, whose placeholder helper can only create files), so `blast_gff/` and the BLAST index are enforced in Python, not by Snakemake.

### Consensus QC (viralQC)

`rules/viralqc.smk` (included by all four workflows) merges every sample's renamed consensus FASTA into `qc/viralqc/input.fasta` (segmented runs suffix headers with `|<segment>`), runs one `vqc run` inside `envs/viralqc.yaml` (viralQC pinned from PyPI plus nextclade/BLAST/seqtk; it nests its own Snakemake, hence `snakemake-minimal` and `pulp<2.8` in that env), and slices `results.tsv` per sample. The `run_viralqc` rule is lenient by design: a `vqc` failure writes a placeholder table and `viralqc_status.txt` and the run still succeeds. `viralconseq setup` downloads the databases through `scripts/viralqc_setup.smk` using the same env cache.

### Snakemake workflows

Each top-level `.smk` in `viralconseq/scripts/` is small — it includes `rules/common.smk` first (config guard, wildcard constraints, `run.log` hooks, `LOG`/`BENCH`/`cpus`/`ram_mb`/`clair3_model_*` helpers), sets up `TERMINAL_INPUTS` and `rule all` (`default_target: True`), then `include:`s the rule modules from `viralconseq/scripts/rules/` in DAG order (`qc_illumina`, `alignment_{illumina,nanopore}`, `consensus_{illumina,nanopore}`, `stats`, `consensus_{illumina,nanopore}_common`, `viralqc`, `collect`, `report`), defines `organize_files`, and includes `rules/provenance.smk` last. `docs/workflow.md` documents every rule's command and flags.

Cross-cutting conventions to know before editing rules:

- **Per-rule conda envs.** Every rule has `conda: "envs/<name>.yaml"` (relative to the workflow file; `../envs/` from `rules/`). Adding a new tool means either reusing an env or adding a YAML there. Snakemake's `--use-conda` is enabled in `_orchestrator.run_workflow`. `setup_cli.py --dry-run` scans the `include:`/`conda:` directives to list the envs a workflow needs; the real build lets Snakemake create them via `conda_create_envs_only`.
- **Optional outputs are computed conditionally.** `TERMINAL_INPUTS` appends `isnvs/isnvs_summary.tsv` only when `run_isnv` is set, `qc/viralqc/outputs/results.tsv` only when `run_viralqc` is set (default on) and `report.html` only when `run_report` is set (default on); `organize_files` and the collection rules mirror the flags with `... if <flag> else []`. Whenever you add an optional step, edit every place that lists it, and keep it ON in `ConfigGenerator.write_skeleton` so `setup` materialises its env.
- **`organize_files` builds the browse tree; `collect_benchmarks` is the terminus.** `organize_files` creates the per-sample `samples/<sample>/...` symlinks and touches `samples/.organized`; new per-sample outputs need a `ln -sf` block there to be discoverable. Each entry file defines `TERMINAL_INPUTS` (the analysis products) before `rule all`; `rules/provenance.smk`, included at the very end, defines `collect_benchmarks`, which depends on `TERMINAL_INPUTS`, the sentinel, `versions.tsv` and `config.yml` and writes `benchmark.tsv`. A new terminal artifact goes into `TERMINAL_INPUTS`, not into `rule all` directly.
- **Logs and benchmarks come from `LOG()`/`BENCH()`** (`rules/common.smk`): `log: LOG("<own rule name>")` for per-sample rules (per segment in segmented workflows), `LOG("<rule>", target="<rule>", per_segment=False)` for run-level rules. Every log lands under `logs/<rule>/[<segment>/]<target>`, which is the layout `scripts/python/collect_benchmarks.py` relies on; the inventory test rejects any other form.
- **Two kinds of helper scripts.** `scripts/python/*.py` run either via `script:` (the two legacy ones, injected `snakemake` global) or via `shell:` as stdlib argparse programs exposing `main(argv)` (`tool_versions.py`, `collect_benchmarks.py`, and every new one). Prefer the second form: tests can drive them exactly as the workflow does. Neither kind may import the `viralconseq` package.
- **Every new tool needs a version probe.** `scripts/python/tool_versions.py` holds the probe table; the `versions_<env>` rule for its environment in `rules/provenance.smk` lists the tools to probe.
- **Reference sanitization (nanopore).** The nanopore workflow sanitizes reference FASTA headers (replacing `/ \ | , ~` and spaces with `_`) before use because Clair3 makes per-contig directories from the seq IDs. Don't bypass this; `integrity.sanitize_nanopore_contig` must stay in sync with the `sed` in `sanitize_reference`.
- **Multi-contig single reference.** `align_consensus_to_reference_genome` builds the cross-sample alignment one reference contig at a time (gofasta aborts on multi-contig references) and concatenates them into `samples_alignment.fasta`.
- **Empty nanopore samples degrade, they do not fail.** `check_mapped_reads` writes `assembly/status/<sample>[.<segment>].txt` (`status ok|no_mapped_reads`, key–value TSV); below `--minimum-mapped-reads` (default 10) `infer_consensus_sequence` skips Clair3 and writes an all-N consensus, header-only VCFs and a `model.txt` carrying the status, so every downstream rule sees a well-formed input and `summary.tsv` shows `no_mapped_reads`. Nanopore `fastq_empty` is therefore only a warning at validation; Illumina keeps it as an error.
- **The report is a port of dezi_pipe's template.** `scripts/templates/report.html` is vanilla JS + hand-drawn SVG with a `/*{{DATA}}*/` marker; `scripts/python/build_report.py` reads `summary.tsv` by column name, builds every other path from the run layout (no globbing), and refuses to write a page whose figures disagree with `summary.tsv` (`self_check`). `viralconseq create-report` (`report.py`) drives the same script; `test/report_test.py` pins its flags against the script's parser and against `rules/report.smk`, so add a flag in all three places. No external `src=`/`href=` is allowed on the page. There is no JS runtime in CI: after editing the template, open a real run's `report.html` in a browser (headless Chrome with `--dump-dom` works) and check the console.

### Shell conventions in rules

A Snakemake `shell:` body is a Python format string first and a bash script second. The traps below have each cost a debugging session; `test/rule_inventory_test.py` enforces the mechanical ones.

- `set -euo pipefail` first, then `exec > {log} 2>&1` (or `exec 2> {log}` when stdout is the output). Snakemake deletes a job's `log:` before it runs but **not** its `benchmark:`; `{benchmark}` is not a shell placeholder. Where several processes share one log they must all **append** (`>>`); a truncating `2>` erases what the others wrote. Keep a one-line console notice where a user needs it (`exec 3>&2` before the redirect, then `echo ... >&3`).
- **`set -e` exempts the left side of `&&` and of `||`.** `a && b` swallows a failure of `a`; write separate statements. Under `pipefail`, `grep -q` closing a pipe early yields status 141: prefer a single `awk` that reads everything (the GSAlign emptiness test in `rules/consensus_illumina.smk` is the reference).
- Braces are doubled (`{{ }}`) for awk/bash, and **backslash escapes are resolved by Python before bash sees them**, including inside `#` comments: write `\0`, `\n`, `\t` (so `tr '\0' 'N'`, `printf "...\n"`), never ` `/`
`, or the shell receives a NUL byte or a literal line break. Heredocs containing `{}` are eaten by the formatter; put such code in `scripts/python/` instead.
- External tools that exit non-zero for a known reason get `|| rc=$?` followed by a `case` that names the cause (the `run_clair3.sh` block: 1 = arguments/env, 2 = checkpoint failed to load, 127 = env broken). Clair3 2.x parses bare boolean flags with `nargs="?"`, so value flags go first and booleans last with nothing positional after them.
- Rule names are load-bearing: `LOG("<own name>")`, `cpus("<own name>")`, the `versions_<env>` probes and `benchmark.tsv`'s `rule` column all use them; renaming a rule means grepping `rules/`, `test/rule_inventory_test.py`, `ResourceDefaults` and `docs/workflow.md`.
- Snakemake shell strings are interpolated **unquoted** unless you write `{x:q}`; every path and every list param goes through `:q`. The two `*_flags` config strings are the one deliberate exception, which is why `validate_flag_strings` and the parse-time shlex guard exist.

### Sample sheets

CSV with no header. Illumina has 3 columns (`sample_id,R1,R2`), Nanopore has 2 (`sample_id,fastq`). `create-samplesheet` builds them by scanning a run directory; the parser (`viralconseq/validators.py:validate_sample_sheet`) keys off the data type given on the CLI and rejects rows with the wrong column count.

`summary.tsv` (built by `scripts/python/build_summary.py` in `rules/collect.smk`) is the run's "start here" table: one row per sample with a `status` column, the per-sample statistics (`scripts/python/calculate_assembly_stats.py`, written as `coverage_stats/<sample>.stats.tsv`), iSNV counts, viralQC columns and the Clair3 model. Its header is pinned in `build_summary.COLUMNS`; percentages are 0–100 and missing values `NA`. New per-sample facts go into the stats script or into `build_summary.py` as a new column at the end of `COLUMNS`, never into a side table. `assembly/assembly_stats_summary.csv` is a deprecated derivation kept until 0.3.0.

Sample names are prefixed with `sample-` inside the generated YAML by `ConfigGenerator.add_samples`, and the `.smk` files refer to `sample-<id>` everywhere. The prefix is visible in every output: `samples/sample-<id>/`, `assembly/coverage_stats/sample-<id>.table_cov_basewise.txt`, the `sample_id` column of `summary.tsv`, the `sample` column of `benchmark.tsv`, the `samples` block of `run_manifest.json` and the consensus FASTA headers. Tests, docs and config inspection should expect the prefixed form; the bare id exists only in the sample sheet and the `args` dict.

## Tests

The unittest suite under `test/` covers Python-layer behaviour (CLI parsing, validators, integrity checks, config generation, path resolution, provenance). The Snakemake dryrun suite (`test/dryrun_test.py`) runs `snakemake -n` against every workflow + a YAML in `test/dryrun_configs/`; the placeholder fixture `create_dryrun_placeholders.sh` writes empty input files so paths resolve. `test/empirical_test.py` (opt-in, `make test-empirical`) runs the real SARS-CoV-2 scenario in `test/empirical/scenarios/`. CI runs lint + mypy, the unit suite + dryruns (one `test` job per Python version), a Sphinx docs build, a Docker build smoke test and a conda-env smoke test.

When changing rule wiring, run `make test-dryrun` — it catches missing inputs, broken `expand` patterns, and circular dependencies that the Python suite cannot see.

`test/rule_inventory_test.py` is a text-level check over the four entry workflows and their includes: every rule must declare `log:` and `benchmark:` through `LOG("<own name>", ...)`/`BENCH(...)`, every `shell:` body must `set -euo pipefail` and redirect into `{log}`, analysis keys may not carry a literal fallback, a rule may read only its own `<rule>_cpus`/`<rule>_ram` key, `ResourceDefaults.CONSENSUS_*_RULES` must name exactly the rules that read one, and every dry-run config must hold the required keys. Add a new rule with those directives from the start or the suite fails.

## Editing notes

- Don't rename `validate_args`, `generate_config_file`, or `run_snakemake_workflow` in `viralconseq/consensus.py` — multiple tests patch those exact module-level names.
- Two `viralconseq/scripts/python/*.py` scripts (`calculate_assembly_stats.py`, `rename_sequences.py`) still run inside Snakemake's `script:` directive, which injects a `snakemake` global. Ruff would otherwise flag F821; `pyproject.toml` already silences this for that path. Don't add `from snakemake import snakemake` to those files, and don't make any script there import from the `viralconseq` package (they run inside per-rule conda envs where it is not installed). The Python layer may import them the other way round (`report.py` drives `build_report.main`), which is why they must stay stdlib-only.
- Package data is explicit: `pyproject.toml` `[tool.setuptools.package-data]` and `MANIFEST.in` list `*.smk`, `envs/*.yaml`, `templates/*.html` and `python/*.py`; a new kind of shipped file needs both, and the CI `wheel-smoke` job asserts they arrive in an installed wheel.
- `viralconseq rerun --set KEY=VALUE` rewrites the YAML (previous copy in `.bak`) rather than passing overrides in memory: the config on disk is the run's provenance. Keep it that way.
- Snakemake reserves some names for rule keywords' attributes: `copy`, `count`, `index`, `keys`, `items` (and other `list`/`dict` method names) cannot be used as `input:`/`output:`/`params:`/`log:` keys — the parser raises `AttributeError: ... is reserved for internal use` at the rule line. Use `config_copy`, `n_reads`, and so on.
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
