# Architecture

viralconseq is a thin Python orchestration layer over Snakemake workflows. The
Python code validates inputs, writes a YAML config, and launches one of four
Snakemake workflows; all the bioinformatics lives in the `.smk` rule files.

## Data flow

```
CLI (click)                       viralconseq/cli.py  →  consensus_cli.py
   │  parses options, configures logging (--log-level / --json-logs)
   ▼
main(args)                        viralconseq/consensus.py
   │  thin wrapper around…
   ▼
_orchestrator.run_pipeline(...)   viralconseq/_orchestrator.py
   │  resolve_paths → validate → generate_config → [run manifest] → run_workflow
   ▼
ConfigGenerator → YAML            viralconseq/config_generator.py
   │  emits the exact keys the .smk files read (config["samples"], …)
   ▼
snakemake(workflow.smk, config)   viralconseq/scripts/consensus_<datatype>[_segmented].smk
   │  rule all → include: rules/*.smk
   ▼
per-rule conda envs               viralconseq/scripts/envs/*.yaml  (--use-conda)
```

## Module map

| Module | Responsibility |
|---|---|
| `cli.py` | Top-level click group; wires `--log-level` / `--json-logs`. |
| `consensus_cli.py` | Click options for `consensus illumina` / `consensus nanopore` → plain `args` dict. |
| `consensus.py` | Owns `validate_args`, `generate_config_file`, `run_snakemake_workflow`; calls the orchestrator. |
| `create_samplesheet.py` | `create-samplesheet` subcommand. |
| `setup_cli.py` | `setup` subcommand: pre-builds per-rule conda envs into a shared cache. |
| `_orchestrator.py` | Shared `run_pipeline` skeleton (resolve → validate → config → manifest → run) with structured error handling. |
| `validators.py` | File existence, sample-sheet parsing, reference/primer checks, input sanitization, and content-level input-integrity orchestration (`validate_consensus_input_integrity`). |
| `integrity.py` | Streaming, pure-stdlib content validators for FASTQ / FASTA / BED; collect `IntegrityIssue`s rather than raising. |
| `reference_splitter.py` | Splits a multi-record `--reference` FASTA into the per-segment `{segment: path}` dict the segmented workflows consume. |
| `config_generator.py` | Writes the YAML config (the contract with the `.smk` files). |
| `constants.py` | `ConfigKeys`, `DataType`, `ResourceDefaults` (per-workflow rule lists). |
| `exceptions.py` | Typed error hierarchy with machine-readable `code`s. |
| `logging_config.py` | Central logging (run id, text/JSON). |
| `provenance.py` | `run_manifest.json` (version, config, input checksums). |
| `scripts/*.smk`, `scripts/rules/*.smk` | The actual workflows. |
| `scripts/python/*.py` | Helpers run via Snakemake's `script:` directive. |

## The config is the contract

`ConfigGenerator` emits keys that are hard-coded as strings in the `.smk` files
(e.g. `config["samples"]`, `config["output"]`, `config["run_isnv"]`).
Adding a pipeline option means touching four places: the click option in
`consensus_cli.py`, the `validators.py` check, a `ConfigGenerator.add_*` setter,
and the rule(s) that read it. `config["samples"]` maps `sample-<id>` to a **list**
of FASTQ paths (readers tolerate the legacy space-joined string form too).

## Workflow selection

`run_snakemake_workflow` formats a path:
`scripts/consensus_{illumina,nanopore}[_segmented].smk`.
Segmentation is chosen when `args["reference"]` is a dict — built either from
repeatable `--segmented-reference SEGMENT=PATH`, or by `reference_splitter.py`
auto-splitting a single multi-record `--reference` FASTA during validation
(unless `--single-reference`).
