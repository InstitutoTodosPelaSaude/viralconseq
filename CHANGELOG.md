# Changelog

All notable changes to this project will be documented in this file.

The format is loosely based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The release process is documented in [RELEASING.md](RELEASING.md).

## [Unreleased]

### Added

- `<run>/report.html`: a self-contained interactive report of the run (rule
  `report`, `scripts/python/build_report.py`, shell
  `scripts/templates/report.html`; vanilla JS and hand-drawn SVG, no external
  resources, light/dark theme, print CSS). Headline figures, a sortable and
  filterable table of every `summary.tsv` row, a per-sample panel with the
  depth trace and masked regions, variant counts, fastp figures (Illumina) or
  the Clair3 model and mapped-read gate (nanopore), primer extents and
  amplicon depths when a scheme was given, the viralQC verdict with its six
  checks, and tool versions and run parameters in the footer. Segmented runs
  facet by segment; unsegmented runs with several viruses facet by virus. The
  page self-checks against `summary.tsv` before it is written. Toggle with
  `--report/--no-report` (config key `run_report`, section `report`).
- `<run>/consensus/`: a flat, share-ready copy of the consensus sequences
  (rule `collect_consensus`, `scripts/python/collect_consensus.py`):
  `sample-<id>[.<segment>].fasta` per sample with headers normalised to
  `sample-<id>[|<contig>][|<segment>]` and one-line sequences,
  `consensus[.<segment>].fasta` pooling every sample without the reference,
  and `consensus[.<segment>].cov<T>.fasta` with the samples whose
  `coverage_min_depth` reaches the new `--consensus-coverage-threshold`
  (percent, default 70, config key `consensus_coverage_threshold`).

- `<run>/versions.tsv`: the versions of every tool that ran, probed inside the
  per-rule conda environments at run time (`component<TAB>version`), with the
  viralconseq and Snakemake versions on top and, when viralQC ran, the database
  directory and the date its Nextclade datasets were downloaded. Produced by
  the new `scripts/rules/provenance.smk` (one `versions_<env>` rule per
  environment plus `versions`) and `scripts/python/tool_versions.py`. The
  config gains a `viralconseq_version` key (section `provenance`).
- `<run>/config.yml`: a copy of the resolved config inside the run directory
  (rule `run_config`), so a run is self-describing wherever `--config-file`
  pointed.

- Resources that bind. `--threads-total` defaults to the cores available to
  the process minus one (CPU affinity, so container and cgroup limits count).
  A memory budget is detected from the machine (`MemTotal` capped by the
  cgroup limit, minus 10 %) and passed to Snakemake as `--resources mem_mb=`,
  so the `mem_mb` declared by Clair3 (`infer_consensus_sequence`, 2 GB) and
  viralQC (`run_viralqc`, 1 GB) now limits how many of those jobs run at once.
  `--max-memory GB` overrides the budget, `--max-memory 0` disables it, and a
  positive value below the largest per-rule figure is refused. A one-line
  summary of cores and budget is printed at start. The config records
  `threads_total`, `max_memory_mb` and `memory_detected_mb`.

- Automatic Clair3 model selection. `--clair3-model` now defaults to `auto`:
  the basecall model tag Dorado/MinKNOW write into every read header
  (`basecall_model_version_id=` or `RG:Z:`) is read from the first reads of
  each sample and mapped to a Clair3 model the way ARTIC's `choose_model`
  does (same manifest, same tie-breaks, same refusal to fall back to a
  Guppy-era model for Dorado reads). Samples basecalled with different models
  get a per-sample mapping in the config; a sample whose reads disagree is an
  error. Reads without a tag (Guppy-era, re-headered) stop the run before any
  work with the model to pass by hand.
- `--clair3-model-dir` (default `~/.cache/viralconseq/clair3-models`, or
  `$VIRALCONSEQ_CLAIR3_MODELS`): models live in viralconseq's own cache, not
  in the conda env, and the chosen model's two checkpoints are validated
  (zip-level) before the run starts and declared as rule inputs. A missing or
  incomplete model fails with the `viralconseq setup --clair3-models NAME`
  command to run. New error codes `clair3_model_unresolved`,
  `clair3_model_mixed_within_sample`, `clair3_model_not_found`.
- `<run>/summary.tsv`, the run's "start here" table: one row per sample (per
  sample and segment when segmented) with a `status` column (`ok`,
  `no_mapped_reads`, `empty_consensus`, `viralqc_failed` / `_partial` /
  `_skipped` / `_missing`, `missing_stats`), read counts, mean and median
  depth, breadth of coverage in percent, consensus length and N content, iSNV
  count, viralQC virus / clade / lineage / grade / score / dataset, and the
  Clair3 model. Built by the new `scripts/python/build_summary.py` (rule
  `summary`, `rules/collect.smk`); pinned header. Every sample of the sheet
  has a row.
- Nanopore samples without mapped reads no longer abort the run. A new rule,
  `check_mapped_reads`, writes `assembly/status/<sample>.txt` (`status`,
  `mapped_reads`, `minimum_mapped_reads`); below `--minimum-mapped-reads`
  (default 10, `0` disables) the sample gets an all-N consensus, header-only
  VCFs and an empty low-coverage BED instead of a Clair3 crash, viralQC grades
  it `D`, and the run continues. A nanopore FASTQ with zero reads is now an
  input-integrity warning rather than an error (Illumina keeps the error).
  `samples/<sample>/status.txt` links the status file.
- `viralconseq setup --clair3-models NAME...` (default: `r1041_e82_400bps_sup_v500`,
  `r1041_e82_400bps_hac_v500`, `r941_prom_hac_g360+g422`; `all` for the whole
  29-model manifest) downloads Clair3 models into `--clair3-model-dir` from the
  Clair3 authors' server with the ARTIC mirror as fallback, validating each
  checkpoint before publishing it; `--skip-clair3-models` opts out and
  `--dry-run` lists what would be fetched. `scripts/python/fetch_clair3_model.py`
  is the stdlib downloader and can be run on its own.
- `viralconseq rerun CONFIG.yml`: run a workflow again from a saved config,
  to resume an interrupted run, `--dry-run` a change, `--unlock` a stale lock,
  `--keep-going`, or `--set KEY=VALUE` (parsed as YAML, written back with the
  previous file kept as `.bak`). The config is validated first with the new
  `validators.validate_config_dict` (required keys, types and bounds, flag
  strings, memory budget), which the generated config is also checked against
  before it is written and which the workflows apply at parse time to
  hand-edited YAMLs.

### Changed

- Per-sample assembly statistics are a headered TSV
  (`assembly/coverage_stats/<sample>.stats.tsv`, linked as
  `samples/<sample>/stats.tsv` on both platforms) with the new columns
  (median depth, consensus length and N content, `pct_mapped`), percentages
  0–100, and read counts that include both Illumina mates (0.1.x counted R1
  reads against mapped mates, so `pct_mapped` could exceed 100).
  `assembly/assembly_stats_summary.csv` is still written with its old header
  and fractions, derived from the same rows, and is **deprecated**: it will be
  removed in 0.3.0. The `unify_assembly_statistics_reports` rule is gone.
- `scripts/envs/utils.yaml` no longer carries pandas or biopython (nothing
  used them); `calculate_assembly_stats.py` is standard-library only, and the
  `pandas` dev extra is gone.
- **Breaking for nanopore runs:** `--clair3-model` no longer defaults to
  `r1041_e82_400bps_sup_v500`, and the model is read from `--clair3-model-dir`
  instead of the conda environment. Runs on untagged reads (Guppy-era data)
  must now pass `--clair3-model NAME`, and every model must be fetched once
  with `viralconseq setup --clair3-models NAME`.
- Nanopore variant calling uses Clair3 2.0.2 (was 1.2). Clair3 2.x ships
  PyTorch models (`pileup.pt` + `full_alignment.pt`); the `run_clair3.sh`
  invocation passes value flags before the bare boolean flags (2.x parses
  booleans as optional-value arguments), decodes Clair3's exit codes into a
  readable message (a stale TensorFlow-era model directory is a bare exit 2)
  and records the model used in `assembly/clair3/<sample>/model.txt`, linked
  as `samples/<sample>/clair3_model.txt`.
- Per-rule threads follow `--<rule>-cpus` > `--threads`: a `--<rule>-cpus`
  option now defaults to unset and the config carries `<rule>_cpus` only when
  you gave it; every rule otherwise uses the `--threads` baseline (rules used
  to default to 2 regardless of `--threads`).
- `mem_mb` is declared only by rules with a measured peak. The
  `--perform-qc-ram`, `--map-reads-ram`, `--trim-primer-sequences-ram` and
  `--detect-isnv-ram` options are gone (they never bound: no budget was ever
  passed). `--infer-consensus-sequence-ram` and `--run-viralqc-ram` remain and
  default to the measured 2 GB and 1 GB instead of 4 GB.
- All per-rule logs and benchmarks now live under one tree,
  `<run>/logs/<rule>/[<segment>/]<target>.{log,benchmark.txt}` (`target` is
  `sample-<id>` for per-sample rules and the rule name for run-level ones).
  `assembly/logs/`, `logs/fastp/` and `logs/consensus_<data>/` are gone.
- `benchmark.tsv` is produced by a new terminal rule, `collect_benchmarks`
  (stdlib `scripts/python/collect_benchmarks.py`), instead of a bash loop in
  `organize_files`. Columns: `sample`, `[segment,]`, `rule`, `target`,
  `threads`, then Snakemake's measured columns read from each file's header.
  `rule` holds the rule name (`ampliconclip`, `gsaalign`, `consensus` used to
  be tool names), run-level rows are `All`, and every benchmarked rule appears
  (two were missed before). `organize_files` now only builds the `samples/`
  symlink tree and marks completion with `samples/.organized`.

## [0.1.1] - 2026-09-17

### Changed

- Snakemake now runs with `<output>/<run_name>/` as its working directory, so
  `.snakemake/` (locks, metadata, transcripts) lives next to the results
  instead of in whatever directory the CLI was launched from. The newest
  transcript is also copied to `<run>/logs/snakemake.log` and recorded in
  `run_manifest.json` (`snakemake_log`). Sample FASTQ paths are written to the
  config as absolute paths (they used to be copied verbatim from the sheet).
  A run interrupted mid-job resumes on the next invocation instead of failing
  with `IncompleteFilesException`.
- `run_manifest.json` keys its `samples` block by `sample-<id>`, like every
  other output.
- New shared preamble `scripts/rules/common.smk`, included first by all four
  workflows: the `sample` wildcard now matches exactly the configured sample
  ids (no ambiguity between ids that are prefixes of one another), `rule all`
  is the default target (a bare `snakemake -s <workflow>` plans the whole
  DAG), and `onstart`/`onsuccess`/`onerror` append to `<run>/logs/run.log`.
- Analysis parameters (`minimum_depth`, `af_threshold`, `af_isnv_threshold`,
  `chunk_size`, `clair3_model`, `variant_quality`, `variant_depth`,
  `minimum_map_quality`, `minimum_length`) are required config keys: the rules
  read them without a fallback, so the CLI is the single source of defaults.
  The old rule-side fallbacks had drifted from the CLI (min depth 10 vs 20,
  AF 0.5/0.7 vs 0.51, chunk 50000 vs 10000, variant depth 5 vs 10, MAPQ 20 vs
  30). A hand-edited YAML missing a key fails at parse time with the key named.
- `--threads`, `--threads-total`, every `--<rule>-cpus`/`--<rule>-ram` and
  `setup --threads` reject values below 1 at the command line.
- The generated config YAML is written atomically (temporary sibling file,
  then rename) and is organised into commented sections (`run`, `consensus`,
  `read_qc`, `isnv`, `clair3`, `viralqc`, `resources`) that explain each group
  of keys, so the file doubles as the run's documentation. Keys and values are
  unchanged; only comments and ordering differ.
- The config-only flag strings `viralqc_extra_flags` and
  `minimap2_consensus_align_flags` are validated with `shlex` and refused if
  they contain a shell metacharacter, both by the Python layer and at
  Snakefile parse time for hand-edited YAMLs (the rules splice them unquoted
  into shell commands).
- Every rule now declares `log:` and `benchmark:` and redirects its shell
  output into the log (`sanitize_reference`, nanopore `trim_primer_sequences`,
  `calculate_coverage_basewise`, `rename_sequences`, nanopore
  `calculate_assembly_statistics`, `split_viralqc_results`, `organize_files`
  and the nanopore `unify_assembly_statistics_reports` had no log or benchmark;
  `perform_qc`, Illumina `trim_primer_sequences`, `detect_isnv`, Illumina
  `infer_consensus_sequence`, `generate_multiqc_report` and `summarize_isnvs`
  declared a log and never wrote to it). The nanopore `trim_primer_sequences`
  rule now also honours `--trim-primer-sequences-cpus/-ram`. Tool output that
  used to scroll past on the console is in `logs/` and `assembly/logs/`.
- `benchmark.tsv` keeps the `sample-<id>` prefix in its `sample` column, like
  every other output; it was the one table that stripped it. `organize_files`
  also drops its own stale benchmark row before aggregating on a rerun.
- The rule modules copied from ViralUnity v1.5.0 are no longer kept
  byte-identical; divergences are listed here as they happen.

### Fixed

- Illumina runs linked `samples/<sample>/consensus.fasta` to the consensus whose
  header still carried the reference name (`>MN908947.3`); it now points at the
  renamed FASTA (`>sample-<id>`), as the nanopore workflow always did.
- `generate_vcf_consensus` (GSAlign) declared a log but never wrote to it, so a
  fallback to the header-only mock VCF left an empty log and only a stderr
  warning. GSAlign's output, exit status and the fallback reason now land in
  `assembly/logs/gsaalign/<sample>.log`; a one-line `WARNING` still reaches the
  console.
- Illumina consensus VCFs were always the header-only mock file, even for
  consensuses with real substitutions. The "does the consensus contain any
  A/C/G/T" test was a `grep -v | grep -q` pipeline under `set -o pipefail`: the
  early exit of `grep -q` killed the upstream `grep` with SIGPIPE, the pipeline
  returned 141 and every sample took the mock branch. The test is now a single
  `awk` invocation, and GSAlign runs (and its VCF is kept) for every consensus
  with called bases.

### Added

- **Consensus QC with viralQC, on by default.** Every run now ends with one
  [viralQC](https://github.com/InstitutoTodosPelaSaude/viralQC) 1.2.0 invocation
  over all final consensus sequences (segmented runs merge every segment with
  `sample-<id>|<segment>` headers). Outputs: `qc/viralqc/outputs/results.tsv`
  (virus, clade, genome-quality grade A–D per sequence), `qc/viralqc/viralqc_status.txt`,
  and a per-sample slice symlinked at `samples/sample-<id>/viralqc.tsv`.
  New rule module `scripts/rules/viralqc.smk` and env `scripts/envs/viralqc.yaml`.
- `--run-viralqc/--no-run-viralqc` (default on) and `--viralqc-db PATH`
  (default `$VIRALCONSEQ_VIRALQC_DB` or `~/.cache/viralconseq/viralqc-db`) on both
  `consensus` subcommands; `--run-viralqc-cpus/--run-viralqc-ram` resource flags;
  config-only key `viralqc_extra_flags`.
- `viralconseq setup` now also downloads the viralQC databases (Nextclade
  datasets + NCBI RefSeq viral BLAST set, ~1 GB on disk) into `--viralqc-db`
  via `scripts/viralqc_setup.smk`; `--skip-viralqc-db` opts out; rerunning is a
  no-op once the databases exist; `--dry-run` reports the database step.
- `ViralQCDatabaseNotFoundError` (code `viralqc_database_not_found`): a run with
  viralQC enabled aborts before any work starts if the database directory is
  missing or incomplete, naming the directory and the `viralconseq setup` command.
  Completeness includes the `makeblastdb` index next to `blast.fasta`, because
  viralQC degrades a `blastn` failure to "no hits" and would otherwise report
  every unmatched sequence as `Unclassified` without any error.
- `run_manifest.json` records whether viralQC ran and the database directory
  used (`viralqc` entry).

### Changed

- A consensus run now requires `viralconseq setup` to have been run once (or
  `--no-run-viralqc`). Existing per-rule environments are reused; only the new
  `viralqc.yaml` env is built.
- If viralQC fails at run time (a tool error), the run still completes: a
  placeholder `results.tsv` and `viralqc_status.txt` are written and a
  `WARNING` is printed. QC grades never affect the exit code. Delete
  `qc/viralqc/` and rerun to retry the step.
- When a run produces no consensus sequences at all (for example no reads map),
  viralQC is skipped with `status skipped` and an empty results table instead of
  failing the run; samples with an empty consensus are warned about and are
  absent from the table.
- The QC step needs outbound HTTPS even with the databases in place
  (`nextclade sort` fetches its reference minimizer index from the Nextclade
  server on every run); use `--no-run-viralqc` on air-gapped nodes.

## [0.1.0] - 2026-09-11

First release of viralconseq as a standalone package. viralconseq is the
consensus-inference half of [ViralUnity](https://github.com/InstitutoTodosPelaSaude/ViralUnity),
extracted from ViralUnity v1.5.0 so that it can be validated, benchmarked and
published independently of the metagenomics pipeline.

### Pipeline behaviour

The rule modules (`scripts/rules/*.smk`), per-rule conda environment files
(`scripts/envs/*.yaml`) and helper scripts (`scripts/python/*.py`) are
**byte-identical** to ViralUnity v1.5.0. The four entry-point workflows
(`consensus_{illumina,nanopore}{,_segmented}.smk`) differ only by the removal of
the HTML-report / annotation-staging rules and the `viralconseq:` message
prefix. Alignment, primer clipping, variant calling, consensus calling,
coverage statistics and the multi-sample alignment are unchanged, and the
generated Snakemake config is identical key for key.

### Removed (relative to ViralUnity v1.5.0)

- The metagenomics pipeline (`meta`), database download helpers (`get-databases`)
  and `build-deacon-index`.
- The interactive HTML report (`report.html`), the `report` subcommand and the
  `--generate-html-report/--no-generate-html-report` flag. A redesigned report will
  return in a future release.
- `--gene-annotation` / `--segmented-gene-annotation` (GFF3 inputs), including their
  auto-splitting and warn-only validation. Their only consumer was the report's
  annotation track. Primer-scheme (BED) handling is untouched.
- `plotly`, `jinja2`, `biopython` and `pandas` are no longer runtime dependencies of
  the Python layer (`pandas` is still pinned inside `envs/utils.yaml`, where the
  statistics script runs, and is a `dev` extra for its unit test).

### Changed (relative to ViralUnity v1.5.0)

- Package, console script and import name: `viralconseq` (`viralconseq consensus
  illumina|nanopore`, `viralconseq setup`, `viralconseq create-samplesheet`).
- `setup --pipelines` accepts `consensus-illumina`, `consensus-nanopore` and `all`.
- Environment variable `VIRALCONSEQ_CONDA_PREFIX` and default cache directory
  `~/.cache/viralconseq/conda-envs` (previously `VIRALUNITY_CONDA_PREFIX` /
  `~/.cache/viralunity/conda-envs`). Point the new variable at an existing cache to
  reuse already-built environments: the per-rule environment files are identical, so
  their content hashes match.
- `run_manifest.json` records the tool version under `viralconseq_version`.
- Base exception classes are `ViralConseqError` / `ViralConseqFileNotFoundError`
  (error code `viralconseq_error`).
- Snakemake `onsuccess`/`onerror` messages are prefixed with `viralconseq:`.
- The Snakemake dry-run test suite targets `rule all` explicitly (as the CLI does),
  so the non-segmented nanopore dry-run now exercises the whole DAG instead of only
  `sanitize_reference`, and two extra dry-run configs cover primer clipping and
  iSNV calling (`consensus_illumina__isnv_primers.yaml`, `consensus_nanopore__primers.yaml`).

### Lineage

Capabilities inherited from ViralUnity's consensus pipeline (see ViralUnity's
[CHANGELOG](https://github.com/InstitutoTodosPelaSaude/ViralUnity/blob/main/CHANGELOG.md)
up to v1.5.0 for their history):

- Illumina and Nanopore reference-guided consensus workflows with amplicon primer
  clipping and per-sample coverage statistics.
- Segmented-virus workflows driven either by repeatable `--segmented-reference
  SEGMENT=PATH` or by auto-splitting a multi-record `--reference` FASTA
  (`--single-reference` opts out).
- Multi-contig single references, with per-contig cross-sample alignments.
- Content-level input-integrity validation of FASTQ, reference FASTA and primer BED
  (`--skip-input-validation` bypasses it).
- Optional intra-host SNV calling with LoFreq (`--run-isnv`).
- `setup` / `--conda-prefix` shared per-rule environment cache.
- Provenance manifest (`run_manifest.json`), structured machine-readable errors,
  JSON logging, and in-process embedding via `viralconseq.consensus.main`.
- Configurable per-rule CPU/RAM (`--<rule>-cpus`, `--<rule>-ram`) and the
  config-only `minimap2_consensus_align_flags` override.
