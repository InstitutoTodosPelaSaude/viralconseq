# Output Layout

After a successful run, the output directory (`<output>/<run_name>/`) is organised as follows:

```
{output}/{run_name}/
├── run_manifest.json                 # version, config path, input checksums, outcome
├── summary.tsv                       # START HERE: one row per sample, status + stats + viralQC
├── versions.tsv                      # tool versions probed at run time (component, version)
├── config.yml                        # copy of the resolved config this run used
├── report.html                       # self-contained interactive report (unless --no-report)
├── consensus/                        # share-ready FASTAs (one-line sequences)
│   ├── sample-{id}[.{segment}].fasta #   per sample, header sample-{id}[|contig][|segment]
│   ├── consensus[.{segment}].fasta   #   every sample pooled, reference not included
│   └── consensus[.{segment}].cov{T}.fasta  # samples with coverage_min_depth >= T (--consensus-coverage-threshold)
├── logs/
│   ├── run.log                       # start/end lines written by the workflow itself
│   ├── snakemake.log                 # copy of the newest Snakemake transcript
│   └── <rule>/[<segment>/]<target>.{log,benchmark.txt}   # one directory per rule
├── other/versions/                   # per-environment fragments behind versions.tsv
├── .snakemake/                       # Snakemake's own state (locks, metadata, transcripts)
├── assembly/
│   ├── assembly_stats_summary.csv    # deprecated pre-0.2.0 shape of summary.tsv (removed in 0.3.0)
│   ├── coverage_stats/
│   │   └── {sample}.table_cov_basewise.txt
│   ├── consensus/                    # non-segmented
│   │   └── final_consensus/
│   │       ├── samples_alignment.fasta
│   │       └── per_contig_alignments/   # only for a multi-contig --single-reference
│   └── {segment}/                    # segmented (one per segment)
│       └── consensus/
│           └── final_consensus/
│               └── samples_alignment.fasta
├── input_references/                 # segmented: per-segment reference FASTAs split
│                                     #   from a multi-record --reference
├── reference/                        # nanopore only: reference with sanitised headers
│   ├── reference.sanitized.fasta     #   single reference
│   └── {segment}.sanitized.fasta     #   segmented run (one per segment)
├── qc/reports/multiqc_report.html    # illumina only
├── qc/viralqc/                       # unless --no-run-viralqc
│   ├── input.fasta                   # merged consensus sequences handed to viralQC
│   ├── outputs/results.tsv           # one row per consensus sequence
│   ├── outputs/...                   # viralQC's own sub-directories (nextclade_results/, blast_results/, gff_files/, ...)
│   ├── per_sample/                   # per-sample slices of results.tsv
│   └── viralqc_status.txt            # ok / partial / failed and why
├── isnvs/isnvs_summary.tsv           # illumina + --run-isnv only
├── samples/
│   └── sample-{sample_id}/           # symlinks to the per-sample results (note the prefix)
│       ├── consensus.fasta
│       ├── consensus.vcf.gz
│       ├── raw.vcf.gz                # nanopore
│       ├── clair3_model.txt          # nanopore: the Clair3 model used
│       ├── status.txt                # nanopore: ok / no_mapped_reads (see below)
│       ├── isnvs.vcf.gz              # illumina + --run-isnv
│       ├── fastp.html                # illumina
│       ├── stats.tsv                 # this sample's row of the assembly statistics
│       ├── table_cov_basewise.txt    # nanopore
│       ├── viralqc.tsv               # unless --no-run-viralqc
│       ├── raw_mapped_reads.bam
│       └── trimmed_mapped_reads.bam
└── benchmark.tsv                     # per-task runtime
```

Every sample id from the sample sheet appears with a `sample-` prefix everywhere:
in the output tree (`samples/sample-<id>/`), in the `sample_id` column of
`summary.tsv`, in the `sample` column of `benchmark.tsv`, in the
`samples` block of `run_manifest.json`, in the coverage-table file names and in the
consensus FASTA headers.

In segmented runs the per-sample symlinks are nested one level deeper, under
`samples/sample-{sample_id}/{segment}/`, and the per-segment results live under
`assembly/{segment}/`.

On nanopore runs `assembly/status/<sample>[.<segment>].txt` records, per sample,
`status` (`ok`, or `no_mapped_reads` when fewer than `--minimum-mapped-reads`
primary reads mapped), `mapped_reads` and the threshold. A `no_mapped_reads`
sample is not sent to Clair3: its consensus is all `N` (one record per reference
contig), its VCFs are header-only, it is absent from `samples_alignment.fasta`
(nothing aligns) and viralQC grades it `D` with `no usable sequence`; the run
completes and the other samples are unaffected. `assembly/clair3/<sample>/model.txt`
records the Clair3 model and model directory used (and the status, when skipped).

Intermediate files are kept as well (useful for debugging, safe to delete):
`assembly/mapped_reads/{raw,trimmed}/` (BAMs and the ampliconclip `*.trimmed.txt`
reports), `assembly/consensus/final_consensus/` (per-sample consensus FASTAs and VCFs
before symlinking, `aln.consensus.sam`, the indel-masked alignment),
`assembly/clair3/{sample}/` (nanopore variant calls), `assembly/isnvs/` (LoFreq VCFs
with `--run-isnv`), `qc/data/` and `qc/reports/` (fastp-trimmed reads and per-sample
fastp reports, illumina), and `logs/<rule>/` (per-rule logs and
`*.benchmark.txt`). Every rule writes its tool output into its own log, so a
failure is diagnosed from `logs/<rule>/<sample>.log` (per-sample rules) or
`logs/<rule>/<rule>.log` (run-level rules), with a `<segment>/` level in between
for segmented runs; nothing is only on the console. Snakemake runs with the run
directory as its working directory, so `.snakemake/` is here too and can be
deleted once the run is finished.

## Key files

| File | Description |
|------|-------------|
| `summary.tsv` | One row per sample (per sample and segment when segmented): status, read counts, depth, breadth of coverage, consensus length and N content, iSNV count, viralQC virus / clade / grade, Clair3 model. The place to start; columns below |
| `assembly/assembly_stats_summary.csv` | Deprecated: the pre-0.2.0 table (fractions 0–1, mean depth) derived from the same rows; removed in 0.3.0 |
| `report.html` | Interactive, self-contained report of the run (no external resources): headline figures, a sortable/filterable table of every `summary.tsv` row, a per-sample panel with the depth trace, masked regions, variant counts, read QC (Illumina) or Clair3 model and mapped-read gate (nanopore), and the viralQC verdict; footer with tool versions and run parameters. The page refuses to render numbers that disagree with `summary.tsv`. Rebuild it with `viralconseq create-report <run>` |
| `consensus/consensus.fasta` | Every consensus pooled into one multi-FASTA (one-line sequences, no reference record); `consensus/consensus.cov<T>.fasta` keeps only the samples whose `coverage_min_depth` reaches `--consensus-coverage-threshold` (default 70 %), ready to share or submit. Segmented runs write one pair per segment |
| `samples/sample-{id}/consensus.fasta` | Final consensus sequence |
| `samples/sample-{id}/consensus.vcf.gz` | Variants relative to the reference |
| `assembly/coverage_stats/sample-{id}.table_cov_basewise.txt` | Per-base coverage table (`RNAME`, `POS`, `DEPTH`) |
| `samples/sample-{id}/raw_mapped_reads.bam` | Reads mapped to the reference, before primer clipping |
| `samples/sample-{id}/trimmed_mapped_reads.bam` | Primer-clipped BAM used for consensus calling |
| `assembly/consensus/final_consensus/samples_alignment.fasta` | All consensus sequences aligned to the reference (MSA-ready) |
| `qc/viralqc/outputs/results.tsv` | viralQC table: virus, clade and genome-quality score per consensus sequence (see below) |
| `samples/sample-{id}/viralqc.tsv` | The rows of `results.tsv` belonging to one sample |
| `run_manifest.json` | Provenance: version, timestamp, config path and hash, input SHA-256 checksums, outcome, paths of `logs/snakemake.log`, `versions.tsv` and `config.yml` |
| `versions.tsv` | `component<TAB>version`: viralconseq, Snakemake, every tool probed inside its conda environment at run time (minimap2, samtools, bedtools, GSAlign, gofasta, fastp, MultiQC, LoFreq, bcftools, Clair3, Python, seqtk, viralQC, Nextclade, BLAST as applicable) and, when viralQC ran, the database directory and dataset download date |
| `config.yml` | The resolved configuration this run used (a copy of the file `--config-file` pointed at) |
| `logs/run.log` | One line per run start and end, written by the workflow (survives a `.snakemake/` clean-up) |
| `logs/snakemake.log` | Copy of the newest Snakemake transcript for this run directory |
| `benchmark.tsv` | Runtime and resource usage per rule execution (see below) |

## Benchmark columns

`benchmark.tsv` has one row per rule execution, per-sample rows first (in sample
sheet order) and run-level rows (`sample` = `All`) last:

| Column | Meaning |
|---|---|
| `sample` | `sample-<id>`, or `All` for run-level rules |
| `segment` | Segmented runs only: the segment, or `-` for run-level rules |
| `rule` | The Snakemake rule name (`map_reads`, `run_viralqc`, …) |
| `target` | What the rule ran on: the sample id, or the rule name for run-level rules |
| `threads` | The CPUs the rule was given (`--<rule>-cpus`, else `--threads`); empty for rules without a `--<rule>-cpus` option, which run single-threaded |
| `s`, `h:m:s`, `max_rss`, `max_vms`, `max_uss`, `max_pss`, `io_in`, `io_out`, `mean_load`, `cpu_time` | Snakemake's own measurements (seconds, memory in MB, I/O in MB) |

## summary.tsv columns

`summary.tsv` is tab-separated with a pinned header, one row per sample in sample
sheet order (one row per sample and segment in segmented runs, which add a
`segment` column after `sample_id`). Percentages are 0–100 with two decimals;
a value that does not exist for a row is `NA`. Every sample of the sheet has a
row, even when nothing was assembled for it.

| Column | Meaning |
|---|---|
| `sample_id` | `sample-<id>` |
| `status` | `ok`, or why the row is not a normal genome (table below) |
| `total_reads` | Sequenced reads in the input FASTQ(s); both mates on Illumina |
| `qc_passed_reads` | Reads kept by fastp (Illumina); `NA` on nanopore |
| `mapped_reads` | Primary mapped reads in the primer-clipped BAM (each Illumina mate counts once, so `pct_mapped` is bounded by 100) |
| `pct_mapped` | `100 × mapped_reads / total_reads` |
| `mean_depth`, `median_depth` | Depth over all reference positions (all contigs), zeros included |
| `coverage_10x`, `coverage_100x`, `coverage_1000x` | Percent of reference positions at or above each depth |
| `coverage_min_depth` | Percent of reference positions at or above `--minimum-coverage` (the value in `min_depth`); the completeness figure the `consensus/*.cov<T>.fasta` filter and the report use |
| `min_depth` | The `--minimum-coverage` threshold behind `coverage_min_depth` |
| `consensus_length`, `n_count`, `n_pct` | Length of the consensus (all records), number and percent of `N` |
| `isnv_count` | Intra-host variants called (Illumina with `--run-isnv`), else `NA` |
| `virus`, `clade`, `lineage` | viralQC identification (a multi-contig reference reports the contig with the highest coverage) |
| `genome_quality`, `genome_quality_score` | viralQC grade A–D and its 0–24 score |
| `qc_overall_status` | Nextclade's overall verdict (`good` / `mediocre` / `bad`) as reported by viralQC |
| `viralqc_dataset`, `viralqc_dataset_version` | Nextclade dataset and version behind the clade call |
| `clair3_model` | Clair3 model used for this sample (nanopore), else `NA` |

`status` values, in the order they are tested:

| Status | Meaning |
|---|---|
| `missing_stats` | No statistics were produced for the sample (should not happen in a completed run) |
| `no_mapped_reads` | Nanopore: fewer than `--minimum-mapped-reads` reads mapped; the consensus is all `N` and Clair3 was skipped |
| `empty_consensus` | The consensus has no called bases (length 0 or 100 % `N`) |
| `viralqc_failed`, `viralqc_partial`, `viralqc_skipped` | viralQC did not deliver a verdict for the run (see `qc/viralqc/viralqc_status.txt`); consensus outputs are unaffected |
| `viralqc_missing` | viralQC ran but has no row for this sample |
| `ok` | A genome with statistics and, when viralQC ran, a QC verdict |

`samples/sample-<id>/stats.tsv` is the single-row per-sample statistics file the
summary is built from (same statistics columns, no viralQC). The deprecated
`assembly/assembly_stats_summary.csv` keeps the pre-0.2.0 header
(`sample_name`, `number_of_reads`, `number_of_trim_paired_reads`,
`number_of_mapped_reads`, `average_depth`, `percentage_above_{10,100,1000}x`,
`horizontal_coverage`) with fractions 0–1; note that `number_of_reads` and
`number_of_mapped_reads` now count both Illumina mates, where 0.1.x counted R1
reads against mapped mates.

## viralQC results

`qc/viralqc/outputs/results.tsv` is written by [viralQC](https://viralqc.readthedocs.io/)
and has one row per consensus sequence (the reference is not included). The most
useful columns (read them by name, not position):

| Column | Meaning |
|---|---|
| `seqName` | The consensus header, normalised to `sample-<id>[\|<contig>][\|<segment>]`: `sample-<id>` for a single reference, `sample-<id>\|<segment>` in segmented runs, `sample-<id>\|<contig>` for a multi-contig single reference |
| `virus`, `virus_species`, `segment` | Virus identified by Nextclade sort or BLAST (`Unclassified` when nothing matched) |
| `clade` (+ `lineage`, `subclade`, …) | Clade/lineage from the matching Nextclade dataset, when one exists |
| `genomeQuality` | Overall grade A–D (empty when the virus could not be identified) |
| `genomeQualityScore` | The numeric score behind the grade (0–24) |
| `coverage`, `cdsCoverage`, `targetRegionsCoverage` | Fraction of the genome / CDS / target regions covered by called bases |
| `privateMutationsQuality`, `frameShiftsQuality`, `stopCodonsQuality`, … | The individual QC components |
| `dataset`, `datasetVersion` | Nextclade dataset (and version) used |
| `inputSequenceStatus` | Empty when the sequence was analysed; explains why it was not (e.g. no usable bases) |

The full column reference is in viralQC's [output documentation](https://viralqc.readthedocs.io/en/latest/output.html).
`samples/sample-<id>/viralqc.tsv` holds the header plus that sample's rows (all
segments in a segmented run).

`qc/viralqc/viralqc_status.txt` records `status` (`ok`, `partial` or `failed`),
the `vqc` exit code, the number of input records and result rows, and the log
path. If viralQC fails (a tool error, not a bad QC verdict) the run still
completes and `results.tsv` is a placeholder with one row per sequence whose
`inputSequenceStatus` reads `viralQC failed (exit N)`; check
`logs/run_viralqc/run_viralqc.log`.
