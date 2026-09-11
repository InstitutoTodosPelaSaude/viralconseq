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
├── reference/
│   └── reference.sanitized.fasta     # reference with sanitised headers (nanopore)
├── qc/reports/multiqc_report.html    # illumina only
├── isnvs/isnvs_summary.tsv           # illumina + --run-isnv only
├── samples/
│   └── {sample_name}/                # symlinks to the per-sample results
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

In segmented runs the per-sample symlinks are nested one level deeper, under
`samples/{sample_name}/{segment}/`.

## Key files

| File | Description |
|------|-------------|
| `assembly/assembly_stats_summary.csv` | Read counts, mapped reads, average depth, breadth of coverage per sample (and per segment) |
| `samples/{sample}/consensus.fasta` | Final consensus sequence |
| `samples/{sample}/consensus.vcf.gz` | Variants relative to the reference |
| `assembly/coverage_stats/{sample}.table_cov_basewise.txt` | Per-base coverage table (`RNAME`, `POS`, `DEPTH`) |
| `samples/{sample}/raw_mapped_reads.bam` | Reads mapped to the reference, before primer clipping |
| `samples/{sample}/trimmed_mapped_reads.bam` | Primer-clipped BAM used for consensus calling |
| `assembly/consensus/final_consensus/samples_alignment.fasta` | All consensus sequences aligned to the reference (MSA-ready) |
| `run_manifest.json` | Provenance: version, timestamp, config path, input SHA-256 checksums |
| `benchmark.tsv` | Runtime and resource usage per task |

## Assembly statistics columns

`assembly_stats_summary.csv` has one row per sample (per segment in segmented runs):

| Column | Meaning |
|---|---|
| `sample_name` | Sample id (and segment, where applicable) |
| `number_of_reads` | Raw reads in the input FASTQ(s) |
| `number_of_trim_paired_reads` | Reads retained after QC (Illumina); equals the raw count on Nanopore |
| `number_of_mapped_reads` | Reads mapped to the reference |
| `average_depth` | Mean depth over all reference positions |
| `percentage_above_10x` / `_100x` / `_1000x` | Fraction of reference positions at or above each depth |
| `horizontal_coverage` | Fraction of reference positions at or above `--minimum-coverage` |
