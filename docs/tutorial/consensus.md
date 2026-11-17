# Consensus pipeline tutorial

`viralconseq consensus` turns raw viral reads into one consensus genome per sample, aligned to a reference you supply. By the end of this page you will have run it on Illumina and Nanopore SARS-CoV-2 data and you will know how to read each output file.

```{note}
This page assumes you have completed [Setup](setup.md): viralconseq is installed, `viralconseq setup` has run, and `samples_illumina.csv` and `samples_nanopore.csv` point at `my_test_data/`.
```

## When to use it

Use `viralconseq consensus` when you know which virus is in the sample and you have a reference genome for it. The pipeline maps the reads, calls variants, and writes one consensus FASTA per sample plus a multi-sample alignment ready for tree building. It does not classify reads taxonomically: if you do not know what is in a sample, profile it with a metagenomic classifier first.

## What it does

Both data types follow the same idea — *map → clip primers → call → consensus → align all samples to the reference* — with a different variant caller:

```text
Illumina:
  raw FASTQ ──► fastp ──► minimap2 ──► primer clipping ──► [optional depth cap]
            ──► samtools consensus ──► coverage statistics ──► multi-sample alignment
            └► [optional] LoFreq intra-host variants

Nanopore:
  raw FASTQ ──► minimap2 ──► primer clipping ──► [optional depth cap]
            ──► Clair3 + variant filter + bcftools consensus
            ──► coverage statistics ──► multi-sample alignment
```

The tools are wired in. What you choose is:

- the **reference** (`--reference`; a multi-record FASTA is treated as a segmented virus and split per segment),
- whether the data is **amplicon** (`--primer-scheme primers.bed`) or shotgun,
- the **depth and allele-fraction thresholds** that decide which positions become `N` (`--minimum-coverage`, `--af-threshold`).

## Worked example 1 — Illumina (SARS-CoV-2 amplicon)

### Run

```bash
viralconseq consensus illumina \
    --sample-sheet  samples_illumina.csv \
    --config-file   results/consensus_illumina/config.yml \
    --output        results/consensus_illumina/ \
    --run-name      sarscov2 \
    --reference     databases/refs/nCoV-2019.reference.fasta \
    --primer-scheme databases/refs/nCoV-2019.bed \
    --threads 2 \
    --threads-total 4
```

- `--sample-sheet` — the CSV from Setup.
- `--config-file` — where to write the run's configuration file; keep it with the results.
- `--output` and `--run-name` — everything lands under `results/consensus_illumina/sarscov2/`.
- `--reference` — the reference FASTA.
- `--primer-scheme` — the primer BED. Omit it for shotgun data.
- `--threads` and `--threads-total` — threads per step and cores for the whole run. Leave `--threads-total` out to use every core but one.

Snakemake prints one line per step as it runs; each tool's own output goes to a log under `logs/<step>/`, so a failure is diagnosed there. The two example samples finish in a few minutes on a laptop.

```{tip}
`--create-config-only` writes the configuration file and stops, so you can inspect it. `viralconseq rerun <config>` then runs it; `viralconseq rerun` is also how you resume an interrupted run.
```

### Tour the outputs

The files worth knowing under `results/consensus_illumina/sarscov2/`:

```text
summary.tsv                              # START HERE: one row per sample, status + stats + QC
report.html                              # interactive report, opens in any browser
consensus/
├── sample-itps-0001.fasta               # one consensus per sample, ready to share
├── consensus.fasta                      # all samples pooled
└── consensus.cov70.fasta                # only samples covering >= 70 % of the genome
samples/sample-<id>/                     # everything about one sample (note the sample- prefix)
├── consensus.fasta
├── consensus.vcf.gz                     # differences from the reference
├── fastp.html                           # read QC report
├── raw_mapped_reads.bam                 # as mapped
├── trimmed_mapped_reads.bam             # primer-clipped; the source of the depth statistics
├── normalised_mapped_reads.bam          # what the consensus was called from (same as trimmed unless --normalise-depth)
├── stats.tsv                            # this sample's statistics
└── viralqc.tsv                          # this sample's consensus QC
assembly/coverage_stats/sample-<id>.table_cov_basewise.txt   # depth at every position
assembly/consensus/final_consensus/samples_alignment.fasta  # all samples + reference, aligned
qc/reports/multiqc_report.html           # combined read QC
qc/viralqc/outputs/results.tsv           # consensus QC for every sample
benchmark.tsv, versions.tsv, config.yml  # runtime per step, tool versions, the configuration used
logs/                                    # one directory per step
```

**`summary.tsv`** is the place to start: one row per sample with a `status` (`ok`, or why not), read counts, mean and median depth, the percentage of the genome covered at the depth gate (`coverage_min_depth`), the N content of the consensus (`n_pct`), and the viralQC `virus`, `clade` and `genome_quality`. A low `coverage_min_depth` or a high `n_pct` tells you a sample is poorly covered without opening a BAM. Every column is described in [Output layout](../output.md#summarytsv-columns).

**`report.html`** shows the same table interactively, with a depth trace per sample. It is a single file that works offline, so it can be mailed or archived with the run.

**`consensus/`** holds the sequences you will share: one FASTA per sample, all samples pooled, and a pooled file restricted to the samples whose `coverage_min_depth` reaches `--consensus-coverage-threshold` (default 70 %).

**`samples/sample-<id>/consensus.fasta`** is the sample's genome. Runs of `N` mark positions with fewer than `--minimum-coverage` reads, or where no base reached `--af-threshold`. The three BAMs show the reads at each stage: as mapped, after primer clipping (what the statistics count), and what the consensus was called from (the same reads, unless you set a depth cap with `--normalise-depth`).

**`samples/sample-<id>/consensus.vcf.gz`** lists how the consensus differs from the reference:

```bash
bcftools view results/consensus_illumina/sarscov2/samples/sample-itps-0001/consensus.vcf.gz | head -30
```

**`assembly/coverage_stats/sample-<id>.table_cov_basewise.txt`** has the depth at every reference position (`RNAME`, `POS`, `DEPTH`), handy for spotting drop-outs:

```bash
awk '$3 < 10' results/consensus_illumina/sarscov2/assembly/coverage_stats/sample-itps-0001.table_cov_basewise.txt | head
```

**`assembly/consensus/final_consensus/samples_alignment.fasta`** is every consensus aligned to the reference, ready for a tree builder such as IQ-TREE.

**`versions.tsv`** records the exact version of every tool that touched the data; paste it into a methods section. **`benchmark.tsv`** gives the runtime and memory of every step, useful when planning a larger run.

### Consensus QC (viralQC)

The last step hands every consensus to [viralQC](https://viralqc.readthedocs.io/), which identifies the virus, assigns a clade where a Nextclade dataset exists, and grades the genome. The headline columns are already in `summary.tsv`; the full table is `qc/viralqc/outputs/results.tsv`:

```bash
cut -f 1-4 results/consensus_illumina/sarscov2/qc/viralqc/outputs/results.tsv | column -t -s$'\t'
```

- **`virus` / `clade`** — a sanity check that the reference matched the sample. For the example data expect *Severe acute respiratory syndrome coronavirus 2* and a Pango lineage.
- **`genomeQuality`** — A (complete, clean) to D (fragmentary or suspicious), from coverage, private mutations, frameshifts and premature stop codons. B is fine for most surveillance use, C means "look at the sample", D is not suitable for phylogenetics.

Pass `--no-run-viralqc` to skip the step, for instance on a machine without internet access (viralQC fetches a small Nextclade index on every run). If viralQC itself fails, the run still completes and `qc/viralqc/viralqc_status.txt` says why.

### Optional: intra-host variants

To look for sub-consensus mutations within a sample, add `--run-isnv`:

```bash
viralconseq consensus illumina \
    --sample-sheet  samples_illumina.csv \
    --config-file   results/consensus_illumina_isnv/config.yml \
    --output        results/consensus_illumina_isnv/ \
    --run-name      sarscov2 \
    --reference     databases/refs/nCoV-2019.reference.fasta \
    --primer-scheme databases/refs/nCoV-2019.bed \
    --run-isnv \
    --af-isnv-threshold 0.05 \
    --threads 2 --threads-total 4
```

This runs LoFreq on the primer-clipped reads of every sample. The variants land in `samples/sample-<id>/isnvs.vcf.gz`, the per-sample counts in `isnvs/isnvs_summary.tsv` and the `isnv_count` column of `summary.tsv`. Only sub-consensus alleles (frequency below 0.5) at or above `--af-isnv-threshold` are kept.

## Worked example 2 — Nanopore (SARS-CoV-2)

Same virus, different chemistry. Nanopore's higher per-base error rate calls for a dedicated variant caller, Clair3, instead of the pileup-based consensus used for Illumina.

### Run

```bash
viralconseq consensus nanopore \
    --sample-sheet  samples_nanopore.csv \
    --config-file   results/consensus_nanopore/config.yml \
    --output        results/consensus_nanopore/ \
    --run-name      sarscov2 \
    --reference     databases/refs/nCoV-2019.reference.fasta \
    --primer-scheme databases/refs/nCoV-2019.bed \
    --clair3-model  r941_prom_hac_g360+g422 \
    --threads 4 --threads-total 4
```

What differs from the Illumina run:

- **No read QC step.** Nanopore reads go straight to the mapper; the variant caller is expected to absorb noisy bases.
- **`--clair3-model`** names the Clair3 model matching your basecaller. Leave it at its default, `auto`, for recent data: the pipeline reads the basecalling model from the read headers and picks the model per sample. The tutorial reads were basecalled in 2022 with Guppy and carry no such tag, so `auto` would stop and ask for a model; `r941_prom_hac_g360+g422` is the right one for them. `viralconseq setup` downloads the common models; `viralconseq setup --clair3-models NAME` fetches any other.
- **Reference contig names are made shell-safe.** Clair3 uses the contig name as a directory name, so every character other than letters, digits, `.`, `_` and `-` becomes `_`. If a downstream tool needs the exact contig names used in the BAM and VCF, read them from `reference/reference.sanitized.fasta` in the run directory. See [Notes](../notes.md#reference-header-sanitization).

### What is different in the output

The layout matches the Illumina run, with a few additions under `samples/sample-<id>/`:

- `raw.vcf.gz` — every call Clair3 made; `consensus.vcf.gz` — the calls that passed the filter and were written into the consensus.
- `variant_filter.tsv` — every call with its verdict: written, masked as `N`, or discarded with the reference kept. `summary.tsv` counts the masked ones in `masked_variant_sites`.
- `mask.bed` — the intervals written as `N`, from low depth or masked calls.
- `clair3_model.txt` and `status.txt` — the model used, and whether the sample had enough mapped reads to be called at all. A barcode with fewer than `--minimum-mapped-reads` (default 10) mapped reads gets an all-N consensus and `status no_mapped_reads` in `summary.tsv`; the run continues.

There is no `fastp.html` or `multiqc_report.html`. `consensus.fasta` and `samples_alignment.fasta` behave exactly as on Illumina.

## Segmented viruses

For a multi-segment genome (influenza, the bunyaviruses, …), pass a single multi-record FASTA to `--reference`: viralconseq treats a reference with more than one record as segmented, splits it into one reference per record, and names the segments after the FASTA headers.

```bash
viralconseq consensus illumina \
    --sample-sheet samples_flu.csv \
    --config-file  results/flu/config.yml \
    --output       results/flu/ \
    --run-name     flu \
    --reference    refs/flu/influenza_all_segments.fasta \
    --threads 4 --threads-total 8
```

To name the segments yourself, or supply each from its own file, repeat `--segmented-reference SEGMENT=PATH` instead (it replaces `--reference`):

```bash
viralconseq consensus illumina \
    --sample-sheet samples_flu.csv \
    --config-file  results/flu/config.yml \
    --output       results/flu/ \
    --run-name     flu \
    --segmented-reference HA=refs/flu/HA.fasta \
    --segmented-reference NA=refs/flu/NA.fasta \
    --segmented-reference PB1=refs/flu/PB1.fasta \
    --threads 4 --threads-total 8
```

Every segment is assembled in the same run. `summary.tsv` gets one row per sample and segment, `consensus/` one set of files per segment, and the per-sample files sit under `samples/sample-<id>/<segment>/`. Pass `--single-reference` when a multi-record FASTA is one genome in several contigs rather than several segments. See [Notes](../notes.md#segmented-viruses).

## Tuning consensus quality

The settings you will reach for most often:

- **`--minimum-coverage`** (default `10`) — a position with fewer reads than this becomes `N`. Raise it (for example to `20`) for high-confidence assemblies; lower it only when you would rather see a tentative base than a hole.
- **`--af-threshold`** (default `0.6`) — the fraction of reads a base needs to be written into the consensus. A position where no base reaches it becomes `N`. Raise it (for example to `0.7`) for noisier data.
- **`--normalise-depth`** (default `0`, off) — cap the reads the consensus is called from at this many per position, for example `100`, to keep a very deep amplicon run fast. At a capped position the allele fraction is measured on that sample of reads, so a site close to `--af-threshold` can land on either side; pass `--primer-scheme` with a cap so every amplicon keeps its share.
- **`--minimum-read-length`** (default `50`) — reads shorter than this are dropped after trimming and primer clipping.

Nanopore has three more:

- **`--variant-quality`** (default `10`) — the Clair3 quality a call needs to be written. A call below it is masked as `N`, not silently reverted to the reference.
- **`--frameshift-variant-quality`** (default `50`) — the same bar for an indel that shifts the reading frame. Such indels are rare in viral genomes and usually sequencing errors, so the bar is higher.
- **`--variant-depth`** (default `10`) — the number of reads that must carry the variant. A call with fewer is masked.

Every option is in the [commands reference](../commands.md#viralconseq-consensus); the reasoning behind the filters is in the [workflow reference](../workflow.md).

## Common pitfalls

- **Sample sheet columns must match the subcommand.** `consensus illumina` needs three columns and `consensus nanopore` two; a mismatch is rejected with the offending row named.
- **Primer scheme contig names must match the reference.** viralconseq checks this before running and stops if no primer matches any reference contig, instead of silently clipping nothing.
- **Inputs are checked before the run.** Every FASTQ, the reference and the primer BED are read up front; a truncated or mislabelled file stops the run immediately rather than deep inside the workflow. `--skip-input-validation` bypasses the check.
- **`--reference` and `--segmented-reference` are mutually exclusive.** Pass one or the other.
- **`viralqc_database_not_found` at start.** The consensus QC databases are missing: run `viralconseq setup` once, or pass `--no-run-viralqc`.
- **Reading contig names after a Nanopore run.** Use the names in `reference/reference.sanitized.fasta`, not those in your input FASTA.

## Reference

- [Commands reference — `viralconseq consensus`](../commands.md#viralconseq-consensus): every option.
- [Output layout](../output.md): every output file and column.
- [Notes](../notes.md): segmented viruses, amplicon handling, header sanitisation, consensus QC.
