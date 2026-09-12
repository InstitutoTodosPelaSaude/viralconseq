# Output Layout

After a successful run, the output directory (`<output>/<run_name>/`) is organised as follows:

```
{output}/{run_name}/
├── run_manifest.json                 # version, config path, input checksums
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
│       ├── raw_mapped_reads.bam
│       └── trimmed_mapped_reads.bam
└── benchmark.tsv                     # per-task runtime
```

Every sample id from the sample sheet appears with a `sample-` prefix in the output
tree (`samples/sample-<id>/`), in the `sample_name` column of
`assembly_stats_summary.csv`, in the coverage-table file names and in the consensus
FASTA headers. Only `benchmark.tsv` reports the bare id.

In segmented runs the per-sample symlinks are nested one level deeper, under
`samples/sample-{sample_id}/{segment}/`, and the per-segment results live under
`assembly/{segment}/`.

Intermediate files are kept as well (useful for debugging, safe to delete):
`assembly/mapped_reads/{raw,trimmed}/` (BAMs and the ampliconclip `*.trimmed.txt`
reports), `assembly/consensus/final_consensus/` (per-sample consensus FASTAs and VCFs
before symlinking, `aln.consensus.sam`, the indel-masked alignment),
`assembly/clair3/{sample}/` (nanopore variant calls), `assembly/isnvs/` (LoFreq VCFs
with `--run-isnv`), `qc/data/` and `qc/reports/` (fastp-trimmed reads and per-sample
fastp reports, illumina), and `logs/` / `assembly/logs/` (per-rule logs and
`*.benchmark.txt`).

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
| `run_manifest.json` | Provenance: version, timestamp, config path, input SHA-256 checksums |
| `benchmark.tsv` | Runtime and resource usage per task |

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
