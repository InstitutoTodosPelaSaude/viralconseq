# Changelog

All notable changes to this project will be documented in this file.

The format is loosely based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The release process is documented in [RELEASING.md](RELEASING.md).

## [Unreleased]

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
