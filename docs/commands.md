# Commands Reference

## Global options

These go before the subcommand:

| Option | Description |
|--------|-------------|
| `--version` | Print the viralconseq version and exit. |
| `--log-level` | `DEBUG`, `INFO` (default), `WARNING` or `ERROR`. |
| `--json-logs` | Write one JSON object per log line instead of text. |

```bash
viralconseq --log-level DEBUG consensus illumina ...
```

---

## `viralconseq setup`

Builds the per-step conda environments the workflows run in, downloads the viralQC databases and downloads the Clair3 models, each into its own cache directory. Run it once after installing or upgrading viralconseq; every later `viralconseq consensus` run reuses the caches, and rerunning it is a no-op once everything is in place. See [Installation](installation.md#first-time-setup) for the cache locations.

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
| `--clair3-models` | `r1041_e82_400bps_sup_v500 r1041_e82_400bps_hac_v500 r941_prom_hac_g360+g422` | Clair3 models to download (repeatable or comma-separated; `all` fetches every model in the manifest). About 20 MB each; a model already present is not fetched again. |
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
| `--prefix` | — | Prepend `PREFIX_` to every sample ID (a run code, so `barcode05` from two runs stays distinct: `RUN1_barcode05`). Must be a plain identifier. |

Finder and Windows droppings next to the reads (`.DS_Store`, `._*`, `__MACOSX/`,
`Thumbs.db`, `desktop.ini`, `*Zone.Identifier`, `*.temp`, `*.tmp`, `*.part`,
`*.crdownload`) are ignored: they neither become samples nor count as a
sample's files.

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
| `--minimum-coverage` | `10` | Minimum depth for consensus base inclusion. |
| `--normalise-depth` | `100` | Call the consensus from at most this many reads per reference position (after primer clipping; the statistics and iSNV calls still see every read). With a primer scheme every amplicon keeps a share of the positions it covers, so always pass the scheme for amplicon data. At a capped position the allele fraction is measured on that sample of reads, so a site near `--af-threshold` can land on either side; raise the cap or set `0` (off) when such sites matter. Must be `0` or at least `--minimum-coverage` and, on nanopore, `--variant-depth`. See [Workflow: Depth normalisation](workflow.md#depth-normalisation). |
| `--minimum-read-length` | `50` | Minimum read length threshold. |
| `--af-threshold` | `0.6` | Minimum allele fraction for a base to be written into the consensus; a position where no base reaches it becomes `N`. On Illumina it is the fraction of base counts at the position. On Nanopore it is ALT reads over REF plus ALT reads for a Clair3 call: a call below 0.5 is discarded (the reference is kept), one between 0.5 and the bar is masked as `N`. |
| `--consensus-coverage-threshold` | `70` | Percent of reference positions at or above `--minimum-coverage` (`coverage_min_depth` in `summary.tsv`) a sample needs to enter `consensus/consensus.cov<T>.fasta`; 0–100. Every sample is always in `consensus/consensus.fasta`. |
| `--run-name` | `undefined` | Name for the sequencing run. |
| `--threads` | `1` | Threads per individual task (at least 1); the baseline every rule uses unless a `--<rule>-cpus` override is given. |
| `--threads-total` | cores available minus one | Total cores Snakemake may use at once across all running jobs (at least 1). Container and cgroup limits are honoured. |
| `--max-memory` | detected | Memory budget in GB for the steps that declare memory (Clair3 and viralQC): they run concurrently only while their declared RAM fits it. Detected from the machine (cgroup limit honoured) minus 10 % headroom; `0` disables the budget. A positive value below the largest per-step figure is refused. |
| `--create-config-only` | off | Only generate the config file; do not run the workflow. |
| `--skip-input-validation` | off | Skip the content checks of the input files (FASTQ, FASTA, BED). Existence checks still run. Also the only way to run without `seqkit` on `PATH`. |
| `--conda-prefix` | `~/.cache/viralconseq/conda-envs` | Cache directory for per-rule conda envs. Picked up from `$VIRALCONSEQ_CONDA_PREFIX` if set. Pre-warm with `viralconseq setup`. |
| `--run-viralqc` / `--no-run-viralqc` | on | Run viralQC on the final consensus sequences (virus and clade assignment plus genome-quality scoring via Nextclade + BLAST). |
| `--report` / `--no-report` | on | Write `<run>/report.html`, a self-contained interactive page built from `summary.tsv`: run-level figures, a sortable sample table, per-sample depth traces with the masked regions, viralQC verdicts and tool versions. Opens offline in any browser. |
| `--viralqc-db` | `~/.cache/viralconseq/viralqc-db` | Directory with the viralQC databases. Picked up from `$VIRALCONSEQ_VIRALQC_DB` if set. Populate once with `viralconseq setup`. Ignored with `--no-run-viralqc`. |

### Input integrity validation

Before writing the config, the pipeline reads its inputs and refuses to start
on a file that would break the run or silently produce wrong results:

- **FASTQ** (sample reads) — record structure, matching sequence and quality
  lengths, the DNA alphabet and gzip integrity, checked by `seqkit` with up to
  `--threads-total` files at once. Truncated uploads and corrupt archives are
  caught here. A line-wrapped (multi-line) FASTQ is rejected too, because
  fastp cannot read it; unwrap it with `seqkit seq -w 0`.
- **Reference FASTA** — at least one record, unique contig ids, non-empty
  sequences and a nucleotide-only alphabet, so a protein FASTA or a mislabelled
  file is rejected.
- **Primer BED** (when given) — at least three columns, valid integer intervals,
  and chrom names that match the reference contigs. If no chrom matches any
  contig the run stops (the scheme is wrong for this reference); if some match,
  the unmatched rows are reported as warnings, which is what you see when a
  whole-scheme BED is reused on a subset of segments.

All problems in a file are reported together, with the file and line named.
Pass `--skip-input-validation` to bypass these checks (existence checks still
run).

### Consensus QC (viralQC)

Unless `--no-run-viralqc` is given, the last step of every run hands all
consensus sequences to [viralQC](https://github.com/InstitutoTodosPelaSaude/viralQC),
which identifies the virus of each sequence (Nextclade datasets first, BLAST
against the NCBI RefSeq viral set for anything else), assigns a clade where a
dataset exists, and grades each genome A–D from coverage, private mutations,
frameshifts and stop codons. The headline columns go into `summary.tsv`; the
full table is `qc/viralqc/outputs/results.tsv`, with a per-sample slice at
`samples/sample-<id>/viralqc.tsv` (see [Output layout](output.md#viralqc-results)).

The step needs the databases downloaded by `viralconseq setup`. If the
`--viralqc-db` directory is missing or incomplete the run stops before any
work starts with error code `viralqc_database_not_found` and the exact
`viralconseq setup` command to run; `--no-run-viralqc` is the only way around
the check. It also needs outbound HTTPS even with the databases in place
(`nextclade sort` fetches a small index on every run), so use
`--no-run-viralqc` on an air-gapped node.

If viralQC itself fails at run time (a tool error, not a poor QC verdict), the
run still completes: `results.tsv` is a placeholder, `qc/viralqc/viralqc_status.txt`
records what happened and a warning is printed. QC verdicts are data in the
`genome_quality` column and never affect the exit code. To retry viralQC after
fixing the cause, delete `<output>/<run_name>/qc/viralqc/` and rerun the same
command.

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
| `--clair3-model` | `auto` | Clair3 model for variant calling. `auto` reads the basecalling model tag that Dorado and MinKNOW write into every read header and picks the matching model per sample. Reads without a tag (Guppy-era, re-headered, SRA dumps) stop the run before any work starts with a message naming this option; pass the model yourself then, e.g. `r941_prom_hac_g360+g422` for R9.4.1 Guppy hac reads. Move-table (`*_with_mv`) models are refused. |
| `--clair3-model-dir` | `~/.cache/viralconseq/clair3-models` | Directory holding the downloaded models, one sub-directory each. Picked up from `$VIRALCONSEQ_CLAIR3_MODELS` if set. A model the run needs but that is not there stops the run before any work starts, with the `viralconseq setup --clair3-models NAME` command to fetch it. |
| `--variant-quality` | `10` | Minimum Clair3 QUAL for a SNV or in-frame indel to enter the consensus. A call below it is masked (`N`), not reverted to the reference base. |
| `--frameshift-variant-quality` | `50` | The same bar for an indel that shifts the reading frame. Such calls are rare in viruses and usually errors, so the bar is higher. A frameshifting indel below it is masked (`N`). |
| `--variant-depth` | `10` | Minimum number of reads carrying the variant (counted in the depth-capped BAM) for a Clair3 call to enter the consensus; a call below it is masked (`N`). `--normalise-depth` must be `0` or at least this value. |
| `--minimum-map-quality` | `30` | Minimum mapping quality (clair3). |
| `--minimum-mapped-reads` | `10` | A sample with fewer primary mapped reads than this is not sent to Clair3: it gets an all-N consensus and `status no_mapped_reads` in `summary.tsv`, and the run continues. A crash guard for (near-)empty barcodes, not a QC threshold; `0` disables it. |

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
  `assembly/consensus/final_consensus/per_contig_alignments/<contig>.fasta` (one MSA each).
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
viralconseq rerun <output>/config.yml --set normalise_depth=0   # call from every read; diff the consensus
```

### Options

| Option | Description |
|--------|-------------|
| `CONFIG_FILE` | The YAML written by `--config-file` (or its copy `<run>/config.yml`). |
| `--dry-run` | Plan only (`snakemake -n`). |
| `--unlock` | Release the lock left by a run that was killed, then exit. |
| `--keep-going` | Keep running independent jobs after one fails. |
| `--set KEY=VALUE` | Change a config key (repeatable). The value is parsed as YAML (`30` is an integer, `true` a boolean, `'x'` a string). Only keys already in the file, or keys the current release requires for the data type (so a 0.2.x config can be completed with `--set normalise_depth=100` and, on nanopore, `--set frameshift_variant_quality=50`), may be set. The file is rewritten, with the previous copy kept as `CONFIG_FILE.bak`, so the config on disk always matches the run. |
| `--conda-prefix` | Directory of the cached per-rule conda envs (`$VIRALCONSEQ_CONDA_PREFIX` or `~/.cache/viralconseq/conda-envs`). |

The config is validated before Snakemake starts (required keys, types and
bounds, flag strings, memory budget against the per-rule figures); an invalid
file or override fails with `[configuration_error] ...` and nothing runs. The
workflow (Illumina or nanopore, single or segmented reference) is chosen from
the `data` and `reference` keys. `run_manifest.json` in the run directory, when
present, gets its `status` updated.

---

## `viralconseq create-report`

Rebuild `report.html` for a finished run directory (`<output>/<run_name>`),
for example after `--no-report`, for a run made by an older version, or to
give the page another header label:

```bash
viralconseq create-report results/run1
viralconseq create-report results/run1 -o /tmp/run1_report.html --label "Run 1, re-analysed"
```

### Options

| Option | Description |
|--------|-------------|
| `RUN_DIR` | The run directory; must contain `summary.tsv`. |
| `-o`, `--output` | Where to write the page (default `RUN_DIR/report.html`). |
| `--label` | Run name shown in the page header (default: the directory's name). |

The parameters (data type, minimum depth, coverage threshold, segments, primer
scheme, display-only thresholds) are read from `RUN_DIR/config.yml`, the copy
the workflow leaves behind; without it the data type is inferred from the
layout and the CLI defaults are used, with a warning. The page is built by the
same script as the workflow's `report` rule and fails with `[report_error]` if
its figures would contradict `summary.tsv`.

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
2. `--threads` — the threads each step gets (default 1).
3. `--<rule>-cpus` — an override for one computationally significant step (flags replace `_` with `-`). Precedence is `--<rule>-cpus` > `--threads`. Illumina: `perform_qc`, `map_reads`, `trim_primer_sequences`, `normalise_depth`, `detect_isnv`, `run_viralqc`. Nanopore: `map_reads`, `trim_primer_sequences`, `normalise_depth`, `infer_consensus_sequence`, `run_viralqc`.

Memory is declared only where it was measured: `--infer-consensus-sequence-ram` (Clair3, default **2** GB) and `--run-viralqc-ram` (viralQC, default **1** GB). Those two steps run concurrently only while their declared RAM fits the `--max-memory` budget; every other step peaked well under 200 MB on the test data and declares nothing.

```bash
viralconseq consensus nanopore ... \
    --threads 2 --threads-total 16 --max-memory 24 \
    --infer-consensus-sequence-cpus 4 --infer-consensus-sequence-ram 3
```

The summary line printed at start (`resources: 16 core(s) for Snakemake | memory budget 24.0 GB ...`) shows what was resolved. The [workflow reference](workflow.md#resources) explains how the values reach Snakemake.
