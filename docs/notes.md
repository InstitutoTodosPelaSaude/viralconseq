# Notes

## Segmented viruses

viralconseq assembles segmented viral genomes in one run. Pass a single multi-record FASTA to `--reference`: when it holds more than one sequence, viralconseq treats it as segmented and splits it into one reference per record, naming the segments after the FASTA headers (the first token of each header, sanitised to letters, digits, `.`, `_` and `-`). The split references are written under `<output>/<run_name>/input_references/`.

To name the segments yourself, or supply each from its own file, pass `--segmented-reference SEGMENT=PATH` once per segment instead:

```bash
viralconseq consensus illumina \
    --sample-sheet samples.csv \
    --config-file config_segmented.yml \
    --output /path/to/output \
    --run-name run1 \
    --segmented-reference S=/path/to/s_segment.fasta \
    --segmented-reference L=/path/to/l_segment.fasta
```

Either way the pipeline processes every segment in parallel. Results are organised under `assembly/<segment>/` and `samples/sample-<id>/<segment>/`, `consensus/` holds one set of files per segment, and `summary.tsv` has one row per sample and segment.

## Multi-contig single reference

A genome that is only available as several contigs but should be treated as **one** reference (for example to report one whole-assembly coverage figure) can be passed as a multi-record FASTA together with `--single-reference`. Coverage statistics are then aggregated across all contigs, the consensus keeps one record per contig, and the cross-sample alignment is written per contig under `assembly/consensus/final_consensus/per_contig_alignments/<contig>.fasta`, with a concatenation kept as `samples_alignment.fasta`.

## Reference header sanitization

The Nanopore workflow rewrites the reference FASTA headers before use: every character other than letters, digits, `.`, `_` and `-` becomes `_`. Clair3 turns the sequence ID into a directory name and a shell argument, so parentheses, colons, quotes, spaces and similar characters (an unedited influenza RefSeq header, for instance) would otherwise break the run. The sanitised copy is written to `reference/reference.sanitized.fasta` (`reference/<segment>.sanitized.fasta` in segmented runs) and its contig names are the ones that appear in the BAM and VCF files. The input check compares primer-BED chrom names against the sanitised names, and warns (`fasta_header_sanitized`) when a header contains a character that releases before 0.3.0 would have left in place, because the contig name has changed since then.

## Amplicon vs. shotgun data

`--primer-scheme` is optional. With a BED file, `samtools ampliconclip` removes the primer sequences from the mapped reads before consensus calling, and both the raw and the clipped BAM are kept. Without a scheme the clipping step is a pass-through.

In both cases the consensus is called from every clipped read unless `--normalise-depth` sets a cap (off by default; `100` is a reasonable value for deep amplicon runs). With a scheme the cap is split per strand and per amplicon, so every amplicon keeps a share of the positions it covers; without one, the amplicon that starts first would own each overlap, which is why a capped run should always have its scheme. The statistics, the Nanopore mapped-read gate and the iSNV calls always see every read. The [workflow reference](workflow.md#depth-normalisation) describes what the cap does and does not preserve.

## Consensus QC (viralQC)

The final stage of every run is [viralQC](https://github.com/InstitutoTodosPelaSaude/viralQC). All consensus sequences of the run are merged into `qc/viralqc/input.fasta` (in segmented runs the headers become `sample-<id>|<segment>`) and analysed in one `vqc run`.

**Databases.** `viralconseq setup` downloads the Nextclade datasets and the NCBI RefSeq viral BLAST set into one directory (`--viralqc-db`, `$VIRALCONSEQ_VIRALQC_DB`, or `~/.cache/viralconseq/viralqc-db` by default). `viralconseq consensus` checks that the directory is complete before any work starts and stops with `viralqc_database_not_found` otherwise; `--no-run-viralqc` skips both the check and the step. To refresh the databases (new Nextclade datasets, a new RefSeq release) delete the directory and rerun `viralconseq setup`. The databases are local, but `nextclade sort` still fetches a small index from the Nextclade server on every run, so the step needs outbound HTTPS; use `--no-run-viralqc` on air-gapped nodes.

**Failure policy.** A viralQC tool error never fails the run: the step writes a placeholder `results.tsv`, records the outcome in `qc/viralqc/viralqc_status.txt` and prints a warning. QC grades are data, not exit codes. To retry after fixing the cause, delete `qc/viralqc/` in the run directory and rerun.

**Tuning.** `--run-viralqc-cpus` and `--run-viralqc-ram` size the step; other `vqc run` options (BLAST task, identity, e-value, Nextclade sort thresholds) can be appended through the config-only key `viralqc_extra_flags` (see [Configuration file overrides](commands.md#configuration-file-overrides)).
