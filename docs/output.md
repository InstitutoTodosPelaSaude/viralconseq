# Output Layout

After a successful run, the output directory (`<output>/<run_name>/`) is organised as follows:

```
{output}/{run_name}/
├── run_manifest.json                 # version, config path, input checksums, outcome
├── versions.tsv                      # tool versions probed at run time (component, version)
├── config.yml                        # copy of the resolved config this run used
├── logs/
│   ├── run.log                       # start/end lines written by the workflow itself
│   ├── snakemake.log                 # copy of the newest Snakemake transcript
│   └── <rule>/[<segment>/]<target>.{log,benchmark.txt}   # one directory per rule
├── other/versions/                   # per-environment fragments behind versions.tsv
├── .snakemake/                       # Snakemake's own state (locks, metadata, transcripts)
├── assembly/
│   ├── assembly_stats_summary.csv    # per-sample QC metrics
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
│       ├── isnvs.vcf.gz              # illumina + --run-isnv
│       ├── fastp.html                # illumina
│       ├── stats_summary.csv         # illumina
│       ├── table_cov_basewise.txt    # nanopore
│       ├── viralqc.tsv               # unless --no-run-viralqc
│       ├── raw_mapped_reads.bam
│       └── trimmed_mapped_reads.bam
└── benchmark.tsv                     # per-task runtime
```

Every sample id from the sample sheet appears with a `sample-` prefix everywhere:
in the output tree (`samples/sample-<id>/`), in the `sample_name` column of
`assembly_stats_summary.csv`, in the `sample` column of `benchmark.tsv`, in the
`samples` block of `run_manifest.json`, in the coverage-table file names and in the
consensus FASTA headers.

In segmented runs the per-sample symlinks are nested one level deeper, under
`samples/sample-{sample_id}/{segment}/`, and the per-segment results live under
`assembly/{segment}/`.

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
| `assembly/assembly_stats_summary.csv` | Read counts, mapped reads, average depth, breadth of coverage per sample (and per segment) |
| `samples/sample-{id}/consensus.fasta` | Final consensus sequence |
| `samples/sample-{id}/consensus.vcf.gz` | Variants relative to the reference |
| `assembly/coverage_stats/sample-{id}.table_cov_basewise.txt` | Per-base coverage table (`RNAME`, `POS`, `DEPTH`) |
| `samples/sample-{id}/raw_mapped_reads.bam` | Reads mapped to the reference, before primer clipping |
| `samples/sample-{id}/trimmed_mapped_reads.bam` | Primer-clipped BAM used for consensus calling |
| `assembly/consensus/final_consensus/samples_alignment.fasta` | All consensus sequences aligned to the reference (MSA-ready) |
| `qc/viralqc/outputs/results.tsv` | viralQC table: virus, clade and genome-quality score per consensus sequence (see below) |
| `samples/sample-{id}/viralqc.tsv` | The rows of `results.tsv` belonging to one sample |
| `run_manifest.json` | Provenance: version, timestamp, config path and hash, input SHA-256 checksums, outcome, paths of `logs/snakemake.log`, `versions.tsv` and `config.yml` |
| `versions.tsv` | `component<TAB>version`: viralconseq, Snakemake, every tool probed inside its conda environment at run time (minimap2, samtools, bedtools, GSAlign, gofasta, fastp, MultiQC, LoFreq, bcftools, Clair3, Python, pandas, viralQC, Nextclade, BLAST as applicable) and, when viralQC ran, the database directory and dataset download date |
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
| `threads` | The CPUs the rule was given (`--<rule>-cpus`); empty for rules without a resource option |
| `s`, `h:m:s`, `max_rss`, `max_vms`, `max_uss`, `max_pss`, `io_in`, `io_out`, `mean_load`, `cpu_time` | Snakemake's own measurements (seconds, memory in MB, I/O in MB) |

## Assembly statistics columns

`assembly_stats_summary.csv` has one row per sample (one row per sample and segment in
segmented runs, which add a `segment` column after `sample_name`):

| Column | Meaning |
|---|---|
| `sample_name` | Sample id, with the `sample-` prefix |
| `number_of_reads` | Raw reads in the input FASTQ(s) |
| `number_of_trim_paired_reads` | Reads retained after QC (Illumina); equals the raw count on Nanopore |
| `number_of_mapped_reads` | Reads mapped to the reference |
| `average_depth` | Mean depth over all reference positions |
| `percentage_above_10x` / `_100x` / `_1000x` | Fraction of reference positions at or above each depth |
| `horizontal_coverage` | Fraction of reference positions at or above `--minimum-coverage` |

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
