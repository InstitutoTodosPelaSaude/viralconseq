# Notes

## Segmented viruses

viralconseq natively supports the assembly of segmented viral genomes. Instead of running the pipeline multiple times, pass multiple segment references with `--segmented-reference`:

```bash
viralconseq consensus illumina \
    --sample-sheet samples.csv \
    --config-file config_segmented.yml \
    --output /path/to/output \
    --segmented-reference S=/path/to/s_segment.fasta \
    --segmented-reference L=/path/to/l_segment.fasta
```

Alternatively, pass a single multi-record FASTA to `--reference`: when it holds
more than one sequence viralconseq treats it as segmented and splits it into one
reference per record automatically, deriving segment names from the FASTA headers
(first `|`/whitespace token, sanitised to `[A-Za-z0-9._-]`). Use `--single-reference`
to opt out and align all records together as one reference. The split per-segment
references are written under `<output>/<run_name>/input_references/`.

When either form is detected, the pipeline dispatches a specialised modular workflow that processes each segment independently in parallel. Results are organised under `assembly/{segment}/` and `samples/{sample_name}/{segment_name}/`.

## Multi-contig single reference

A genome that is only available as several contigs but should be treated as **one**
reference (e.g. to report one whole-assembly horizontal coverage) can be passed as a
multi-record FASTA together with `--single-reference`. Coverage statistics are then
aggregated across all contigs, the consensus keeps one record per contig, and the
cross-sample alignment is written per contig under
`consensus/final_consensus/per_contig_alignments/<contig>.fasta` (with a concatenation
kept as `samples_alignment.fasta`).

## Reference header sanitization

The nanopore workflow automatically sanitizes reference FASTA headers before use. Special characters (`/`, `\`, `|`, `,`, `~`, and spaces) in sequence identifiers are replaced with underscores (`_`). This prevents issues with downstream tools such as Clair3 that use the sequence ID to create output directories. The sanitized copy is written to `reference/reference.sanitized.fasta`, and the input-integrity check reports primer-BED chrom names against the sanitized contig names.

## Amplicon vs. shotgun data

`--primer-scheme` is optional. When a BED file is given, `samtools ampliconclip` removes primer sequences from the mapped reads before consensus calling and both the raw and the clipped BAM are kept. Without a scheme the clipping step is a pass-through and the two BAMs are identical.

## Per-rule conda environments

Every rule runs inside a pinned conda environment from `viralconseq/scripts/envs/`. Snakemake builds them on first use; `viralconseq setup` pre-builds them into a shared cache (`~/.cache/viralconseq/conda-envs/` or `$VIRALCONSEQ_CONDA_PREFIX`) so that real runs never pay the environment-creation cost.

## Running tests

```bash
make test          # unit suite
make test-dryrun   # snakemake -n against every workflow
```
