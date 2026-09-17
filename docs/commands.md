# Commands Reference

## `viralconseq setup`

Pre-builds the per-rule conda environments declared in the Snakemake workflows into a shared cache directory and downloads the viralQC databases into a second cache directory. Run once after installing viralconseq (or after upgrading); subsequent `viralconseq consensus` runs reuse both caches. Rerunning is a no-op once everything is in place.

```bash
viralconseq setup --pipelines all
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--conda-prefix` | `$VIRALCONSEQ_CONDA_PREFIX` or `~/.cache/viralconseq/conda-envs` | Directory where per-rule envs are cached. Reused by every later pipeline run that points at the same prefix. |
| `--pipelines` | `all` | Repeatable. Pick from `consensus-illumina`, `consensus-nanopore`, or `all`. Segmented variants share envs with their non-segmented counterparts. |
| `--threads` | `4` | Cores given to Snakemake while materializing envs and downloading databases. |
| `--viralqc-db` | `$VIRALCONSEQ_VIRALQC_DB` or `~/.cache/viralconseq/viralqc-db` | Directory where the viralQC databases (Nextclade datasets + BLAST reference set) are downloaded. `viralconseq consensus` reads the same location by default. |
| `--skip-viralqc-db` | off | Only build conda envs; do not download the viralQC databases. |
| `--clair3-models` | `r1041_e82_400bps_sup_v500 r1041_e82_400bps_hac_v500 r941_prom_hac_g360+g422` | Clair3 models to download (repeatable or comma-separated; `all` fetches the whole manifest, 29 models). About 20 MB each, from the Clair3 authors' server with the ARTIC mirror as fallback; each model is validated before it is put in place, and a model already present is not fetched again. |
| `--clair3-model-dir` | `$VIRALCONSEQ_CLAIR3_MODELS` or `~/.cache/viralconseq/clair3-models` | Directory the models go into; `viralconseq consensus nanopore` reads the same location by default. |
| `--skip-clair3-models` | off | Do not download Clair3 models. |
| `--dry-run` | off | Print the envs, databases and models that would be created and exit without invoking conda or downloading anything. |

### Examples

**Build everything once after install:**

```bash
viralconseq setup --pipelines all
```

**Build only what you need:**

```bash
viralconseq setup --pipelines consensus-illumina
```

**Inspect first:**

```bash
viralconseq setup --pipelines consensus-illumina --dry-run
```

**Envs only (e.g. a node that must not download databases):**

```bash
viralconseq setup --skip-viralqc-db --skip-clair3-models
```

Run the pipeline with `--no-run-viralqc` afterwards, or copy a populated `--viralqc-db` directory over from another machine; copy the `--clair3-model-dir` directory the same way (or `cp -r $CONDA_PREFIX/bin/models/<name>` from a Clair3 2.x environment, which ships the same files).

**Fetch one more Clair3 model** (the run told you which with a `clair3_model_not_found` error):

```bash
viralconseq setup --skip-viralqc-db --clair3-models r1041_e82_400bps_hac_v420
```

**Use a shared cache on a cluster:**

```bash
export VIRALCONSEQ_CONDA_PREFIX=/shared/viralconseq-envs
viralconseq setup --pipelines all
```

Every subsequent `viralconseq consensus` invocation will pick up `$VIRALCONSEQ_CONDA_PREFIX` automatically; no per-run flag needed.

---

## `viralconseq create-samplesheet`

Before running the pipeline you need a CSV sample sheet. The `create-samplesheet` command generates it automatically from a sequencing run directory.

```bash
viralconseq create-samplesheet --input /path/to/run/ --output samples.csv
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--input` | *(required)* | Run directory: contains FASTQ files directly (level 0) or one subdirectory per sample (level 1). |
| `--output` | *(required)* | Output CSV file path. |
| `--level` | `1` | `0` = files in `--input`, `1` = files in subdirectories (one per sample). |
| `--separator` | `-` | Character used to split the file/directory name and extract the sample ID (`-`, `_`, or `.`). |
| `--pattern` | `R1` | Pattern identifying the first read file at level 0 (`R1` for Illumina, `barcode` for Nanopore). |

### Output format

The generated CSV has no header row:

```
# Illumina (3 columns)
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz

# Nanopore (2 columns)
barcode01,/path/to/barcode01.fastq
```

### Examples

**Illumina** — FASTQs organised in one subdirectory per sample (level 1, default):

```bash
viralconseq create-samplesheet \
    --input /path/to/illumina_run/ \
    --output /path/to/samples.csv
```

**Illumina** — all FASTQs directly in the run directory (level 0):

```bash
viralconseq create-samplesheet \
    --input /path/to/illumina_run/ \
    --output /path/to/samples.csv \
    --level 0 \
    --separator _ \
    --pattern R1
```

**Nanopore** — one FASTQ per barcode subdirectory:

```bash
viralconseq create-samplesheet \
    --input /path/to/nanopore_run/ \
    --output /path/to/samples_nano.csv \
    --separator _ \
    --pattern barcode
```

---

## `viralconseq consensus`

The consensus pipeline takes raw reads to processed consensus genome sequences with a single command. Select the data type as a subcommand (`illumina` or `nanopore`).

### Options — shared (both data types)

| Option | Default | Description |
|--------|---------|-------------|
| `--sample-sheet` | *(required)* | CSV file with sample IDs and file paths. |
| `--config-file` | *(required)* | Path for the YAML config file to be created. |
| `--output` | *(required)* | Base output directory. |
| `--reference` | — | Reference genome FASTA (mutually exclusive with `--segmented-reference`). A multi-record FASTA is treated as a segmented virus and split into one reference per record, with segment names taken from the headers (first `\|`/whitespace token, sanitized). |
| `--single-reference` | off | Treat a multi-record `--reference` as a single reference (align all records together in one pass) instead of auto-splitting it into segments. |
| `--segmented-reference` | — | Per-segment reference: `SEGMENT=PATH` (repeatable). |
| `--primer-scheme` | — | Primer scheme BED file (amplicon sequencing only). |
| `--minimum-coverage` | `20` | Minimum depth for consensus base inclusion. |
| `--minimum-read-length` | `50` | Minimum read length threshold. |
| `--af-threshold` | `0.51` | Min allele frequency to call variant into consensus. |
| `--run-name` | `undefined` | Name for the sequencing run. |
| `--threads` | `1` | Threads per individual task (at least 1); the baseline every rule uses unless a `--<rule>-cpus` override is given. |
| `--threads-total` | cores available minus one | Total cores Snakemake may use at once across all running jobs (at least 1). Detected from the CPUs available to the process (affinity, so container and cgroup limits are honoured), leaving one core free. |
| `--max-memory` | detected | Memory budget in GB for the rules that declare memory (Clair3 and viralQC): they run concurrently only while their declared RAM fits it. Detected as `MemTotal` capped by the cgroup limit, minus 10 % headroom; `0` disables the budget. Refused when positive but below the largest per-rule figure. |
| `--create-config-only` | off | Only generate the config file; do not run the workflow. |
| `--skip-input-validation` | off | Skip content-level integrity checks of the input files (FASTQ/FASTA/BED). Existence checks still run. |
| `--conda-prefix` | `~/.cache/viralconseq/conda-envs` | Cache directory for per-rule conda envs. Picked up from `$VIRALCONSEQ_CONDA_PREFIX` if set. Pre-warm with `viralconseq setup`. |
| `--run-viralqc` / `--no-run-viralqc` | on | Run viralQC on the final consensus sequences (virus and clade assignment plus genome-quality scoring via Nextclade + BLAST). |
| `--viralqc-db` | `~/.cache/viralconseq/viralqc-db` | Directory with the viralQC databases. Picked up from `$VIRALCONSEQ_VIRALQC_DB` if set. Populate once with `viralconseq setup`. Ignored with `--no-run-viralqc`. |

### Input integrity validation

Before writing the config, the consensus pipeline validates the *content* (not
just the existence) of its critical inputs and refuses to start on a file that
would break the run or silently produce wrong results:

- **FASTQ** (sample reads) — streamed end to end: 4-line record structure, the
  `@`/`+` marker lines, matching sequence/quality lengths, valid sequence and
  quality characters, and gzip integrity. Catches truncated uploads and corrupt
  archives.
- **Reference FASTA** — at least one record, unique contig ids, non-empty
  sequences, and a nucleotide-only (IUPAC) alphabet — so a protein FASTA or a
  mislabeled file is rejected up front.
- **Primer BED** (when provided) — ≥3 columns, valid integer intervals, and
  chrom names that match the reference contigs. If *no* chrom matches any
  reference contig it is a hard error (the whole scheme is wrong for this
  reference, and `samtools ampliconclip` would silently trim nothing); if some
  chroms match, the unmatched rows are downgraded to warnings (e.g. a
  whole-scheme BED reused on a subset of segments). On Nanopore the reference
  headers are sanitized before mapping, so a BED chrom that matches only the bare
  accession is reported with a fix-it hint.

All problems in a file are reported together. Pass `--skip-input-validation` to
bypass these checks (existence checks still run).

### Consensus QC (viralQC)

Unless `--no-run-viralqc` is given, the last step of every run merges all final
consensus sequences into one FASTA and runs [viralQC](https://github.com/InstitutoTodosPelaSaude/viralQC)
on it. viralQC identifies the virus of each sequence (Nextclade datasets first,
BLAST against the NCBI RefSeq viral set for anything else), assigns a clade
where a dataset exists, and scores each genome A–D from coverage, private
mutations, frameshifts and stop codons. Results land in `qc/viralqc/outputs/results.tsv`
(one row per consensus sequence) with a per-sample slice symlinked at
`samples/sample-<id>/viralqc.tsv`; see [Output layout](output.md#viralqc-results).

The step needs the databases downloaded by `viralconseq setup`. If the
`--viralqc-db` directory is missing or incomplete, the run aborts **before any
work starts** with error code `viralqc_database_not_found` naming the directory
and the exact `viralconseq setup --viralqc-db ...` command to run. This check is
an existence check, so it is not bypassed by `--skip-input-validation` and it
also runs under `--create-config-only`; `--no-run-viralqc` is the only opt-out.

If viralQC itself fails at run time (a tool error, not a poor QC verdict), the
run still completes: a placeholder `results.tsv` is written (one row per
sequence, `inputSequenceStatus` explaining the failure), `qc/viralqc/viralqc_status.txt`
records the outcome, and a `WARNING` is printed. QC verdicts themselves are data
in the `genomeQuality` column and never affect the exit code. Snakemake then
considers the step complete: to retry viralQC after fixing the cause, delete
`<output>/<run_name>/qc/viralqc/` and rerun the same command.

Note that the step needs outbound HTTPS even with the databases in place:
`nextclade sort` downloads its reference minimizer index from
`data.clades.nextstrain.org` on every run. On an air-gapped node use
`--no-run-viralqc`.

### Options — Illumina only

| Option | Default | Description |
|--------|---------|-------------|
| `--adapters` | — | Adapter sequences FASTA (fastp QC). |
| `--trim-head` | `0` | Bases to trim from 5′ end. |
| `--trim-tail` | `0` | Bases to trim from 3′ end. |
| `--cut-front-mean-quality` | `10` | fastp cut_front quality threshold. |
| `--cut-tail-mean-quality` | `10` | fastp cut_tail quality threshold. |
| `--cut-right-window-size` | `4` | fastp cut_right window size. |
| `--cut-right-mean-quality` | `15` | fastp cut_right quality threshold. |
| `--af-isnv-threshold` | `0.0` | Min allele frequency for iSNV analysis. |
| `--run-isnv` | off | Run intra-host SNV analysis (LoFreq). |

### Options — Nanopore only

| Option | Default | Description |
|--------|---------|-------------|
| `--chunk-size` | `10000` | Chunk size for clair3 processing. |
| `--clair3-model` | `auto` | Clair3 model for variant calling. `auto` reads the basecall model tag Dorado/MinKNOW write into every read header (`basecall_model_version_id=` or `RG:Z:`) from the first reads of each sample and picks the matching model per sample, the way ARTIC's `choose_model` does (a run whose samples were basecalled with different models gets a per-sample mapping in the config). Reads without a tag (Guppy-era, re-headered, SRA dumps) fail before any work starts with a message naming this option; pass the model yourself, e.g. `r941_prom_hac_g360+g422` for R9.4.1 Guppy hac reads. Move-table (`*_with_mv`) models are refused. |
| `--clair3-model-dir` | `~/.cache/viralconseq/clair3-models` | Directory with one sub-directory per model (`pileup.pt` + `full_alignment.pt`). Picked up from `$VIRALCONSEQ_CLAIR3_MODELS` if set. Every model the run needs must be present and complete or the run stops before any work with the exact `viralconseq setup --clair3-models NAME` command to fetch it. |
| `--variant-quality` | `20` | Minimum variant quality (clair3). |
| `--variant-depth` | `10` | Minimum alt allele depth (clair3). |
| `--minimum-map-quality` | `30` | Minimum mapping quality (clair3). |
| `--minimum-mapped-reads` | `10` | Samples with fewer primary mapped reads than this are not sent to Clair3: they get an all-N consensus, header-only VCFs and `status no_mapped_reads` in `assembly/status/<sample>.txt`, and the run continues (a barcode with zero reads is a warning, not an error, on nanopore). A crash guard for (near-)empty barcodes, not a QC threshold; `0` disables it. |

### Examples

**Illumina — single reference:**

```bash
viralconseq consensus illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --reference /path/to/references/viral_genome.fasta \
    --primer-scheme /path/to/primers.bed \
    --adapters /path/to/adapters.fa \
    --threads 2 \
    --threads-total 4
```

**Illumina — segmented genome, single multi-FASTA** (simplest; one file with all
segments, segment names taken from the FASTA headers):

```bash
viralconseq consensus illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example_segmented.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --reference /path/to/references/all_segments.fasta \
    --threads 2 \
    --threads-total 4
```

**Illumina — segmented genome, one file per segment** (explicit segment names via
repeatable `--segmented-reference`):

```bash
viralconseq consensus illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example_segmented.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --segmented-reference L=/path/to/L_segment.fasta \
    --segmented-reference S=/path/to/S_segment.fasta \
    --threads 2 \
    --threads-total 4
```

**Illumina — single reference with multiple contigs (fragmented genome)**: when a
genome is only available as several contigs/scaffolds but should be treated as **one**
reference (all contigs aligned together in a single pass, e.g. to report one
whole-assembly horizontal coverage), pass the multi-record FASTA with
`--single-reference` so it is *not* auto-split into segments:

```bash
viralconseq consensus illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --reference /path/to/fragmented_assembly.fasta \
    --single-reference \
    --threads 2 \
    --threads-total 4
```

Behavior notes for a multi-contig single reference:

- **Coverage/depth statistics** (`summary.tsv`) are aggregated across all contigs:
  `coverage_min_depth` is the percentage of *all* reference positions at or above
  the threshold, and `mean_depth` / `median_depth` are whole-assembly figures; the
  viralQC columns report the contig with the highest coverage.
- **Consensus** keeps one record per contig (contigs are never fused).
- The cross-sample alignment is written **per contig** under
  `consensus/final_consensus/per_contig_alignments/<contig>.fasta` (one MSA each).
  A concatenation of all of them is also kept as `samples_alignment.fasta` (one
  padded block per contig, not a single rectangular alignment) for backward
  compatibility.

**Nanopore:**

```bash
viralconseq consensus nanopore \
    --sample-sheet /path/to/samplesheet_nano.csv \
    --config-file /path/to/example_nano.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --reference /path/to/reference.fasta \
    --threads 4 \
    --threads-total 4
```

**Config only** (generates config without running the workflow):

```bash
viralconseq consensus illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example.yml \
    --output /path/to/example_output \
    --reference /path/to/reference.fasta \
    --create-config-only
```

---

## `viralconseq rerun`

Run a workflow again from the config YAML a previous `viralconseq consensus`
wrote, without retyping the command:

```bash
viralconseq rerun <output>/config.yml                       # resume / finish an interrupted run
viralconseq rerun <output>/config.yml --dry-run             # show what would run
viralconseq rerun <output>/config.yml --unlock              # release a stale Snakemake lock
viralconseq rerun <output>/config.yml --set minimum_depth=30 --set af_threshold=0.6
```

### Options

| Option | Description |
|--------|-------------|
| `CONFIG_FILE` | The YAML written by `--config-file` (or its copy `<run>/config.yml`). |
| `--dry-run` | Plan only (`snakemake -n`). |
| `--unlock` | Release the lock left by a run that was killed, then exit. |
| `--keep-going` | Keep running independent jobs after one fails. |
| `--set KEY=VALUE` | Change a config key (repeatable). The value is parsed as YAML (`30` is an integer, `true` a boolean, `'x'` a string). Only keys already in the file may be set. The file is rewritten, with the previous copy kept as `CONFIG_FILE.bak`, so the config on disk always matches the run. |
| `--conda-prefix` | Directory of the cached per-rule conda envs (`$VIRALCONSEQ_CONDA_PREFIX` or `~/.cache/viralconseq/conda-envs`). |

The config is validated before Snakemake starts (required keys, types and
bounds, flag strings, memory budget against the per-rule figures); an invalid
file or override fails with `[configuration_error] ...` and nothing runs. The
workflow (Illumina or nanopore, single or segmented reference) is chosen from
the `data` and `reference` keys. `run_manifest.json` in the run directory, when
present, gets its `status` updated.

---

## Configuration file overrides

A few tool-level parameters are tunable only through the YAML config file produced by `--config-file` (they are not exposed as CLI flags because they rarely need to change). The defaults preserve the historical behaviour, so most users can ignore this section.

Open the generated YAML config file after running with `--create-config-only` and edit the corresponding key (`minimap2_consensus_align_flags` sits in the `# --- consensus ---` section, `viralqc_extra_flags` in `# --- viralqc ---`). Both values are spliced unquoted into a shell command, so they must be plain tool flags: an unbalanced quote or a shell metacharacter (`; & | < > ` $ \\`) is refused when the workflow is parsed.

| Config key | Default | Effect |
|------------|---------|--------|
| `minimap2_consensus_align_flags` | `-a --sam-hit-only --secondary=no --score-N=0` | Flags passed to `minimap2` when re-aligning the per-sample consensus back to the reference for the final multiple-sequence alignment. |
| `viralqc_extra_flags` | `""` | Extra flags appended to the `vqc run` command line (e.g. `--blast-task dc-megablast --blast-pident 75`); see `vqc run --help` inside `envs/viralqc.yaml`. |

Example: edit the YAML config to

```yaml
minimap2_consensus_align_flags: "-a --sam-hit-only --secondary=no"
```

then run it with `viralconseq rerun example.yml` (or `viralconseq rerun example.yml --set
minimap2_consensus_align_flags="-a --sam-hit-only --secondary=no"` to skip the manual edit).
Driving Snakemake directly still works; point `--directory` at the run directory (the
`output` key of the YAML) so Snakemake's `.snakemake/` state lands where `viralconseq
consensus` would have put it:

```bash
snakemake -s "$(python -c 'import viralconseq, os; print(os.path.dirname(viralconseq.__file__))')/scripts/consensus_illumina.smk" \
    --configfile example.yml --directory <output>/<run_name> \
    --use-conda --conda-prefix ~/.cache/viralconseq/conda-envs -j 4
```

The analysis keys the CLI writes (`minimum_depth`, `af_threshold`, `chunk_size`,
`clair3_model`, …) are required: the rules read them without a fallback, and a
YAML missing one is refused at parse time with the key named. Regenerate the
file with `--create-config-only` rather than writing it from scratch.

---

## Per-rule CPU and RAM overrides

Three levels bound a run:

1. `--threads-total` — the cores Snakemake may use at once (default: available cores minus one).
2. `--threads` — the threads each rule gets (default 1).
3. `--<rule>-cpus` — an override for one computationally significant rule. Precedence is `--<rule>-cpus` > `--threads`; the config carries `<rule>_cpus` only when you set it.

The rules with a `--<rule>-cpus` option (flags replace `_` with `-`):

- **Illumina** — `perform_qc`, `map_reads`, `trim_primer_sequences`, `detect_isnv`, `run_viralqc`.
- **Nanopore** — `map_reads`, `trim_primer_sequences`, `infer_consensus_sequence`, `run_viralqc`.

Memory is handled separately, and only where it was measured. Two rules declare memory and have a `--<rule>-ram` option (GB): `--infer-consensus-sequence-ram` (Clair3, default **2**) and `--run-viralqc-ram` (viralQC, default **1**). The declarations bind: Snakemake is given the `--max-memory` budget and runs those jobs concurrently only while their declared RAM fits it (with the default budget on a 16 GB machine, Clair3 runs at most 7 samples at once). Every other rule peaked well under 200 MB on the test data and declares nothing. The values land in the `# --- resources ---` section of the config as `<rule>_ram`, `threads_total`, `max_memory_mb` and `memory_detected_mb`.

```bash
viralconseq consensus nanopore ... \
    --threads 2 --threads-total 16 --max-memory 24 \
    --infer-consensus-sequence-cpus 4 --infer-consensus-sequence-ram 3
```

The run summary line printed at start (`resources: 16 core(s) for Snakemake | memory budget 24.0 GB ...`) shows what was resolved.
