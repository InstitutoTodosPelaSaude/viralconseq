# Commands Reference

## `viralconseq setup`

Pre-builds the per-rule conda environments declared in the Snakemake workflows into a shared cache directory. Run once after installing viralconseq (or after upgrading); subsequent `viralconseq consensus` runs reuse the cached envs and skip env creation.

```bash
viralconseq setup --pipelines all
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--conda-prefix` | `$VIRALCONSEQ_CONDA_PREFIX` or `~/.cache/viralconseq/conda-envs` | Directory where per-rule envs are cached. Reused by every later pipeline run that points at the same prefix. |
| `--pipelines` | `all` | Repeatable. Pick from `consensus-illumina`, `consensus-nanopore`, or `all`. Segmented variants share envs with their non-segmented counterparts. |
| `--threads` | `4` | Cores given to Snakemake while materializing envs. |
| `--dry-run` | off | Print the envs that would be created and exit without invoking conda. |

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
| `--threads` | `1` | Threads per individual task. |
| `--threads-total` | `1` | Total threads for the workflow. |
| `--create-config-only` | off | Only generate the config file; do not run the workflow. |
| `--skip-input-validation` | off | Skip content-level integrity checks of the input files (FASTQ/FASTA/BED). Existence checks still run. |
| `--conda-prefix` | `~/.cache/viralconseq/conda-envs` | Cache directory for per-rule conda envs. Picked up from `$VIRALCONSEQ_CONDA_PREFIX` if set. Pre-warm with `viralconseq setup`. |

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
| `--clair3-model` | `r1041_e82_400bps_sup_v500` | Clair3 model for variant calling. |
| `--variant-quality` | `20` | Minimum variant quality (clair3). |
| `--variant-depth` | `10` | Minimum alt allele depth (clair3). |
| `--minimum-map-quality` | `30` | Minimum mapping quality (clair3). |

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

- **Coverage/depth statistics** (`assembly_stats_summary.csv`) are aggregated across
  all contigs: `horizontal_coverage` is the fraction of *all* reference positions at
  or above the threshold, and `average_depth` is the whole-assembly mean.
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

## Configuration file overrides

A few tool-level parameters are tunable only through the YAML config file produced by `--config-file` (they are not exposed as CLI flags because they rarely need to change). The defaults preserve the historical behaviour, so most users can ignore this section.

Open the generated YAML config file after running with `--create-config-only` and edit the corresponding key under the `# parameters` section:

| Config key | Default | Effect |
|------------|---------|--------|
| `minimap2_consensus_align_flags` | `-a --sam-hit-only --secondary=no --score-N=0` | Flags passed to `minimap2` when re-aligning the per-sample consensus back to the reference for the final multiple-sequence alignment. |

Example: edit the YAML config to

```yaml
minimap2_consensus_align_flags: "-a --sam-hit-only --secondary=no"
```

then run Snakemake directly against the edited config:

```bash
snakemake -s viralconseq/scripts/consensus_illumina.smk --configfile example.yml --use-conda -j 4 all
```

---

## Per-rule CPU and RAM overrides

`viralconseq consensus` auto-generates a `--<rule>-cpus` and a `--<rule>-ram` option for each computationally significant Snakemake rule, so you can size individual steps without touching the global `--threads` / `--threads-total`. Every such option defaults to **2** CPUs / **4** GB and is written into the generated YAML as `{rule}_cpus` / `{rule}_ram`.

The exact set depends on the data type — run the relevant subcommand's `--help` to list them all (flags replace `_` with `-` and append `-cpus` / `-ram`):

- **Illumina** — `perform_qc`, `map_reads`, `trim_primer_sequences`, `detect_isnv`.
- **Nanopore** — `map_reads`, `trim_primer_sequences`, `infer_consensus_sequence`.

```bash
viralconseq consensus illumina ... \
    --map-reads-cpus 8 --map-reads-ram 16 \
    --perform-qc-cpus 4
```
