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

The nanopore workflow automatically sanitizes reference FASTA headers before use. Every character outside `[A-Za-z0-9._-]` (one `_` per byte, under `LC_ALL=C`) in a header line is replaced with an underscore (`_`). Clair3 uses the sequence ID to create output directories, passes it to `/bin/sh` unquoted and builds `contig:start-end` region strings from it, so parentheses, colons, quotes and other shell metacharacters (an unedited influenza RefSeq header, for instance) would otherwise break the run. A header carrying a character the old `/ \ | , ~` + space rule did not rewrite (so its contig name is not what earlier releases produced) is reported as a `fasta_header_sanitized` warning naming the new contig; the stock RefSeq header, which only has spaces and a comma, stays silent. The sanitized copy is written to `reference/reference.sanitized.fasta`, and the input-integrity check reports primer-BED chrom names against the sanitized contig names.

## Amplicon vs. shotgun data

`--primer-scheme` is optional. When a BED file is given, `samtools ampliconclip` removes primer sequences from the mapped reads before consensus calling and both the raw and the clipped BAM are kept. Without a scheme the clipping step is a pass-through and the two BAMs are identical. In both cases the clipped BAM is then depth-normalised to `--normalise-depth` records per position (default 100) before consensus calling; with a scheme the cap is split per strand *and* amplicon, so every amplicon keeps a share of the positions it covers (an overlap becomes a roughly even mix of the amplicons covering it, not their original proportions) and reads that span two amplicons (skip-amplicon products) keep their original share rather than being dropped or kept wholesale. Without a scheme the amplicon that starts first owns each overlap, so pass the BED whenever the data are amplicons. Allele fractions at capped positions are then a sample of 100 records; the statistics, the mapped-read gate and the iSNV calls still see every read, so an iSNV frequency and the consensus call at the same site come from different read sets (see [Workflow: Depth normalisation](workflow.md#depth-normalisation) for the measured strand and amplicon shifts and the sampling probabilities).

## Per-rule conda environments

Every rule runs inside a pinned conda environment from `viralconseq/scripts/envs/`. Snakemake builds them on first use; `viralconseq setup` pre-builds them into a shared cache (`~/.cache/viralconseq/conda-envs/` or `$VIRALCONSEQ_CONDA_PREFIX`) so that real runs never pay the environment-creation cost.

## Consensus QC (viralQC)

The final stage of every run is [viralQC](https://github.com/InstitutoTodosPelaSaude/viralQC)
(`vqc`, pinned to 1.2.0 in `viralconseq/scripts/envs/viralqc.yaml`). All final
consensus sequences of the run are merged into `qc/viralqc/input.fasta`
(reference excluded; in segmented runs headers become `sample-<id>|<segment>`
so every segment can go into the same table) and analysed in one `vqc run`.
The databases are local, but `nextclade sort` still fetches its reference
minimizer index from `data.clades.nextstrain.org` on each run, so the QC step
needs outbound HTTPS; use `--no-run-viralqc` on air-gapped nodes.

**Databases.** `viralconseq setup` downloads the Nextclade datasets and the
NCBI RefSeq viral BLAST set into one directory: `--viralqc-db`,
`$VIRALCONSEQ_VIRALQC_DB`, or `~/.cache/viralconseq/viralqc-db` by default.
The directory is considered complete when it contains `blast.fasta` together
with its `makeblastdb` index (`blast.fasta.nin`, or per-volume `.nin` files plus
a `blast.fasta.nal` alias for a split database), `blast.tsv`, `blast_gff/` and
the marker `.nextclade_datasets_ok` (written by `setup` once the Nextclade
download finished). The index is checked separately from the FASTA on purpose:
without it `blastn` fails, and viralQC treats a `blastn` failure as "no hits",
so the run would otherwise succeed with every sequence that no Nextclade
dataset matched silently reported as `Unclassified`. `viralconseq consensus` checks
this during argument validation, before any work starts, and aborts with
`viralqc_database_not_found` otherwise; `--no-run-viralqc` skips both the check and the step. To refresh
the databases (new Nextclade datasets, new RefSeq release) delete the directory
and rerun `viralconseq setup`. The download also leaves NCBI's taxonomy dump
under `$HOME/.taxonkit` (used by taxonkit during the BLAST-set build).

**Failure policy.** A viralQC tool error never fails the run: the rule writes
a placeholder `results.tsv`, records the outcome in `qc/viralqc/viralqc_status.txt`
and prints a `WARNING`. QC grades are data, not exit codes. To retry after
fixing the cause, delete `qc/viralqc/` in the run directory and rerun.

**Tuning.** `run_viralqc` gets its own `--run-viralqc-cpus/--run-viralqc-ram`;
other `vqc run` options (BLAST task, identity, e-value, Nextclade sort
thresholds) can be appended through the config-only key `viralqc_extra_flags`.

## Running tests

```bash
make test          # unit suite
make test-dryrun   # snakemake -n against every workflow
```

The unit suite includes `test/rule_inventory_test.py`, a text-level audit of the
workflow files: every rule must declare `log:` and `benchmark:`, every `shell:`
body must fail fast and redirect into its log, a rule may read only its own
`_cpus`/`_ram` key, analysis keys may not carry a literal fallback, and every
dry-run config must hold the required keys. A new rule that breaks a convention
fails the suite rather than hiding its output.
