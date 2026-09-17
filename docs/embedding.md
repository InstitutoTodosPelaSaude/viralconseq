# Embedding viralconseq in a service

viralconseq is designed to be driven by a web service as its bioinformatics core.
This page documents the integration seam: how to call it in-process, the input
contract you must enforce, and what it emits for observability and provenance.
It stops short of prescribing an HTTP API or job queue — those are the caller's
responsibility.

## Calling the pipeline

`viralconseq.consensus.main(args)` is the entry point; it delegates to
`viralconseq._orchestrator.run_pipeline`, which runs
`resolve_paths → validate → generate_config → write_run_manifest → run_workflow`
and returns an exit code (`0` success, `1` failure).

```python
from viralconseq.consensus import main as consensus_main
from viralconseq.logging_config import configure_logging

configure_logging(level="INFO", json_logs=True, run_id=job_id)
exit_code = consensus_main({
    "data_type": "illumina",
    "sample_sheet": "/jobs/<job>/samples.csv",
    "config_file": "/jobs/<job>/config.yaml",
    "output": "/jobs/<job>/out",
    "run_name": job_id,          # must be a safe identifier — see below
    "reference": "/refs/MN908947.3.fasta",
    "primer_scheme": "/refs/nCoV-2019.scheme.bed",
    "adapters": None,            # illumina: required key, None = fastp auto-detection
    "threads": 8, "threads_total": 8,
    "conda_prefix": "/srv/viralconseq/conda-envs",   # shared per-rule env cache
    "viralqc_db": "/srv/viralconseq/viralqc-db",     # populated by `viralconseq setup`
})
```

Any option accepted by the CLI can be passed as a key (CLI `--foo-bar` becomes
`foo_bar`). The `args` dict is the *post-Click* argument set, so a few keys are
required rather than defaulted: `data_type`, `sample_sheet`, `config_file`, `output`,
`run_name`, `threads`, `threads_total`, `conda_prefix`, a reference (`reference` or
`segmented_reference`), `viralqc_db` (unless `run_viralqc` is `False`) and, for
Illumina, `adapters` (`None` is fine). Tuning
parameters (`minimum_coverage`, `af_threshold`, fastp/Clair3 options, per-rule
`*_cpus`/`*_ram`) fall back to the CLI defaults when omitted. Always pass
`conda_prefix`: without it Snakemake builds the per-rule envs under the job's
working directory, i.e. once per job. Build a fresh dict per call — validation
normalises some values in place (e.g. an absent `primer_scheme` becomes `"NA"`).
The option lists in `viralconseq/consensus_cli.py` are the reference for the full
key set and defaults.

## Per-job isolation

Snakemake runs with the run directory `<output>/<run_name>/` as its working
directory, so its `.snakemake/` lock/state directory and per-run transcripts
live next to the results (a copy of the newest transcript is written to
`<output>/<run_name>/logs/snakemake.log`). Two jobs therefore only collide when
they share an `output` + `run_name`: give every job its own run directory. The
in-process `snakemake()` call is still not safe to run concurrently *within one
process*; use one process per job. Interrupted runs are resumed on the next
invocation (`force_incomplete`), so a killed job can simply be relaunched with
the same arguments.

## Input contract (enforce on untrusted uploads)

- **Sample ids** and **`run_name`** become shell tokens and filesystem path
  components. The sample-sheet parser rejects unsafe ids; for service-supplied
  values use `validators.sanitize_identifier(value, field=...)`.
- Constrain user-supplied output/config paths to a per-job base with
  `validators.ensure_within_base(path, base)` to block path traversal.
- Sample sheets are validated up-front: duplicate ids, ragged rows, wrong column
  counts, and missing files all raise `SampleSheetError` / `ViralConseqFileNotFoundError`.
- **Inputs are content-validated, not just existence-checked**, in
  `validate_args` before any workflow runs (and under `create_config_only`, so a
  bad upload is caught in a config-only preflight). Every sample FASTQ is streamed
  (record structure, seq/quality lengths, gzip integrity, truncation); the
  reference FASTA must be nucleotide-only with unique contig ids; a primer BED's
  chrom names must match the reference contigs. Blocking problems raise
  `InputIntegrityError`. Pass `skip_input_validation=True` in `args` to bypass.
- **Not yet covered:** the FASTQ scan streams each file fully with no size or
  time budget — enforce your own upload-size limit upstream before validation.

## Observability

Call `configure_logging(level, json_logs=True, run_id=<job id>)` once per process
(or per job in a worker). Every log record is stamped with the run id; JSON mode
emits one object per line for ingestion.

## Structured errors

Expected failures raise subclasses of `viralconseq.exceptions.ViralConseqError`,
each with a machine-readable `code` and a `to_dict()` payload:

```python
from viralconseq.exceptions import ViralConseqError
try:
    ...
except ViralConseqError as e:
    return {"status": "error", **e.to_dict()}   # {"error","code","message"}
```

`run_pipeline` already logs these as `[<code>] <message>` and returns `1`;
catch them earlier if you need the structured payload.

`ViralQCDatabaseNotFoundError` (code `viralqc_database_not_found`) is raised at
validation time when viralQC is enabled but its database directory is missing or
incomplete; the message names the directory and the `viralconseq setup` command.

`ConfigurationError` (code `configuration_error`) is raised by
`validators.validate_config_dict` when a config mapping (the one about to be
written, a saved YAML given to `viralconseq rerun`, or one hand-edited and run
with `snakemake -s`) misses a required key or holds a value the rules cannot use;
the message names the key.

Nanopore runs resolve their Clair3 model before any work starts and can raise
`Clair3ModelUnresolvedError` (`clair3_model_unresolved`: the reads carry no
basecall tag or it matches no model; pass `clair3_model` explicitly),
`Clair3ModelMixedError` (`clair3_model_mixed_within_sample`: one sample's reads
were basecalled with different models) or `Clair3ModelNotFoundError`
(`clair3_model_not_found`: the model is not in `clair3_model_dir`; the message
holds the `viralconseq setup --clair3-models` command). After validation
`args["clair3_model"]` is a single name or a `{sample id: name}` mapping.

`InputIntegrityError` (code `input_integrity_error`) additionally carries an
`issues` list in its `to_dict()`; each issue has `path`, `kind`
(`fastq`/`fasta`/`bed`), `code`, `severity`, `message`, and an optional
`line` — surface these to tell the user exactly which file (and line) failed.

## Provenance

Each run writes `<output>/<run_name>/run_manifest.json` with the viralconseq
version (`viralconseq_version`), a UTC timestamp, the resolved config path, whether
viralQC ran and against which database directory (`viralqc`), a `sha256`+size
for every input FASTQ (keyed `sample-<id>`), the outcome (`status`,
`finished_utc`) and the path of the copied Snakemake transcript
(`snakemake_log`) — enough to reproduce a result by record.
Persist it alongside job outputs. Per-rule tool versions are pinned in the
per-rule conda environment files under `viralconseq/scripts/envs/`.
