# Consensus pipeline tutorial

The consensus pipeline (`viralconseq consensus`) turns raw viral reads into a polished consensus genome per sample, aligned to a reference you supply. By the end of this page you will have run it end-to-end on Illumina and Nanopore SARS-CoV-2 data and you will know how to read each output file.

```{note}
This page assumes you have completed [Setup](setup.md): viralconseq is installed, the per-rule environments are built, and you have generated `samples_illumina.csv` and `samples_nanopore.csv` from `my_test_data/`.
```

## When to use it

Use `viralconseq consensus` when you already know what virus is in the sample and you have a reference genome for it. The pipeline maps reads, calls variants, and produces one consensus FASTA per sample plus a multi-sample alignment ready for tree-building. viralconseq does not perform taxonomic classification: if you do not know what is in the sample, profile it with a metagenomic classifier first.

## What it does

The two data types share the same overall idea — *map → trim primers → call → consensus → align all samples to the reference* — but the variant caller and quality-control steps differ:

```text
Illumina:
  raw FASTQ ──► fastp ──► minimap2 -x sr ──► samtools ampliconclip (BED primers)
            ──► samtools consensus ──► per-base coverage ──► multi-sample alignment
            └► [optional] LoFreq (iSNVs)

Nanopore:
  raw FASTQ ──► (sanitize reference headers) ──► minimap2 -x map-ont
            ──► samtools ampliconclip ──► Clair3 + bcftools (variants + consensus)
            ──► per-base coverage ──► multi-sample alignment
```

You do not pick these tools — they are wired in. You pick:

- the **reference** (`--reference` — a multi-record FASTA is auto-split per segment for influenza-style genomes — or `--segmented-reference SEG=PATH` once per segment to name them yourself),
- whether the data is **amplicon** (`--primer-scheme primers.bed` clips amplicon primers) or shotgun,
- the **coverage and allele-frequency thresholds** that decide which bases become `N` (`--minimum-coverage`, `--af-threshold`).

## Worked example 1 — Illumina (SARS-CoV-2 amplicon)

We will assemble the two paired-end samples in `my_test_data/illumina_data/`. The reference and primer scheme below assume the [ARTIC nCoV-2019](https://github.com/artic-network/primer-schemes) bundle, which is the source of `nCoV-2019.reference.fasta` (NCBI `MN908947.3`) and the matching `nCoV-2019.bed`. Substitute paths for whichever copies you have.

### Sample sheet

Already generated in [Setup](setup.md):

```text
itps-0001,my_test_data/illumina_data/itps-0001_R1.fastq.gz,my_test_data/illumina_data/itps-0001_R2.fastq.gz
itps-0002,my_test_data/illumina_data/itps-0002_R1.fastq.gz,my_test_data/illumina_data/itps-0002_R2.fastq.gz
```

### Run

```bash
viralconseq consensus illumina \
    --sample-sheet  samples_illumina.csv \
    --config-file   results/consensus_illumina/config.yml \
    --output        results/consensus_illumina/ \
    --run-name      sarscov2 \
    --reference     databases/refs/nCoV-2019.reference.fasta \
    --primer-scheme databases/refs/nCoV-2019.bed \
    --minimum-coverage 20 \
    --threads 2 \
    --threads-total 4
```

What each flag does:

- `--sample-sheet` — the CSV from above.
- `--config-file` — where to write the generated Snakemake YAML config (handy to keep alongside the results for reproducibility).
- `--output` — base directory for all results.
- `--run-name` — appended to `--output`; everything lands under `results/consensus_illumina/sarscov2/`.
- `--reference` — the reference FASTA. A single-record FASTA is used as-is; a multi-record FASTA is auto-split into one reference per segment (unless `--single-reference`). See [Segmented viruses](#segmented-viruses) below.
- `--primer-scheme` — BED file of primer coordinates. Omit it for shotgun data; the primer-clipping step then becomes a pass-through.
- `--minimum-coverage 20` — positions with fewer than 20 reads after primer clipping become `N` in the consensus.
- `--threads 2 --threads-total 4` — 2 threads per task, 4 cores total across the workflow. Leave `--threads-total` out to use every core but one; a memory budget is detected from the machine too (`--max-memory` overrides it) so the memory-hungry steps never oversubscribe RAM.

```{tip}
Add `--create-config-only` to write the YAML config and stop. The file is organised in commented sections, so you can inspect or edit it before running `snakemake --configfile <config> -s "$(python -c 'import viralconseq, os; print(os.path.dirname(viralconseq.__file__))')/scripts/consensus_illumina.smk" --directory <output>/<run_name> --use-conda --conda-prefix ~/.cache/viralconseq/conda-envs -j 4` yourself.
```

Snakemake prints a job table as each rule runs (`perform_qc`, `map_reads`, `trim_primer_sequences`, `infer_consensus_sequence`, `calculate_assembly_statistics`, `generate_multiqc_report`, `align_consensus_to_reference_genome`, …); the tools' own output goes to per-rule logs under `logs/<rule>/`, so the console stays readable and nothing is lost. On a typical laptop the bundled SARS-CoV-2 samples finish in a few minutes; consult `results/consensus_illumina/sarscov2/benchmark.tsv` after the run for per-rule timing.

### Tour the outputs

After the run, the relevant paths under `results/consensus_illumina/sarscov2/` are:

```text
samples/sample-<id>/           # note the sample- prefix on every sample directory
├── consensus.fasta            # final consensus sequence
├── consensus.vcf.gz           # variants relative to the reference
├── fastp.html                 # fastp QC report
├── raw_mapped_reads.bam       # post-mapping, before primer clipping
├── trimmed_mapped_reads.bam   # primer-clipped — the BAM used for consensus
└── stats.tsv                  # this sample's row of the assembly statistics
summary.tsv                              # START HERE: one row per sample, status + stats + viralQC
assembly/
├── coverage_stats/sample-<id>.table_cov_basewise.txt  # per-base depth
└── consensus/final_consensus/
    └── samples_alignment.fasta          # all samples + reference, MSA-ready
qc/reports/multiqc_report.html           # combined fastp/QC report
qc/viralqc/outputs/results.tsv           # viralQC: virus, clade, genome-quality grade per consensus
samples/<sample>/viralqc.tsv             # that sample's rows of the viralQC table
benchmark.tsv                            # wall time + memory per rule per sample
versions.tsv                             # every tool's version, probed in its own env
config.yml                               # the configuration this run used
logs/run.log, logs/snakemake.log         # start/end lines and the Snakemake transcript
```

How to read each one:

**`samples/sample-<id>/consensus.fasta`** — your finished genome. Long runs of `N` indicate stretches with coverage below `--minimum-coverage` (or where every read disagreed with the reference but no allele exceeded `--af-threshold`). A first sanity check is the proportion of non-N bases — `summary.tsv` reports it as `coverage_min_depth` (and `n_pct`).

**`samples/sample-<id>/consensus.vcf.gz`** — the differences between your sample and the reference, called from the consensus FASTA via GSAlign. Inspect with:

```bash
bcftools view results/consensus_illumina/sarscov2/samples/sample-itps-0001/consensus.vcf.gz | head -30
```

**`assembly/coverage_stats/sample-<id>.table_cov_basewise.txt`** — three columns: `RNAME`, `POS`, `DEPTH`. Useful for spotting drop-outs:

```bash
awk '$3 < 20' results/consensus_illumina/sarscov2/assembly/coverage_stats/sample-itps-0001.table_cov_basewise.txt | head
```

**`samples/sample-<id>/{raw,trimmed}_mapped_reads.bam`** — both exist deliberately. `raw_mapped_reads.bam` is what minimap2 produced; `trimmed_mapped_reads.bam` is the same BAM after `samtools ampliconclip` removed primer sequences (only different when you passed `--primer-scheme`). The consensus is called from the trimmed BAM.

**`summary.tsv`** — one row per sample (`sample_id` carries the `sample-` prefix) with a `status`, read counts (`total_reads`, `qc_passed_reads`, `mapped_reads`, `pct_mapped`), depth (`mean_depth`, `median_depth`), breadth (`coverage_10x` … `coverage_min_depth`, in percent), the consensus `n_pct`, and the viralQC `virus`, `clade` and `genome_quality`. The quickest way to spot low-coverage or poorly-mapping samples without opening a single BAM; a `status` other than `ok` says what went wrong. The pre-0.2.0 `assembly/assembly_stats_summary.csv` is still written (deprecated).

**`assembly/consensus/final_consensus/samples_alignment.fasta`** — all per-sample consensuses plus the reference, aligned (built by `minimap2` followed by `gofasta sam toMultiAlign`). Drop this straight into a tree-builder such as IQ-TREE for a quick phylogeny.

**`benchmark.tsv`** — every rule execution's runtime, memory, and CPU (`sample` is `sample-<id>`, or `All` for run-level rules; `rule` names the Snakemake rule). Useful when you scale up to a real run.

**`versions.tsv`** — the exact version of every tool that touched your data, read from inside its conda environment when it ran. Paste it into a methods section, or diff it against another run's when results differ.

**`logs/`** — `run.log` records when the run started and how it ended; `snakemake.log` is the transcript Snakemake printed; the per-rule sub-directories hold each tool's output. When a rule fails, its log is the first thing to read.

### Consensus QC (viralQC)

The last rule of every run hands all consensus sequences to [viralQC](https://viralqc.readthedocs.io/), which identifies the virus, assigns a clade where a Nextclade dataset exists, and grades each genome. The headline table is `qc/viralqc/outputs/results.tsv`:

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} {print $c["seqName"],$c["virus"],$c["clade"],$c["genomeQuality"],$c["coverage"],$c["inputSequenceStatus"]}' OFS='\t' \
    results/consensus_illumina/sarscov2/qc/viralqc/outputs/results.tsv | column -t -s$'\t'
```

Columns worth checking first (names are stable; positions are not, hence the name-based `awk`):

- **`virus` / `clade`** — a sanity check that the reference matched the sample. For the example data expect *Severe acute respiratory syndrome coronavirus 2* and a Pango lineage / Nextstrain clade.
- **`genomeQuality`** — A (complete, clean) to D (fragmentary or suspicious), derived from `genomeQualityScore`, which combines coverage, private-mutation counts, frameshifts and premature stop codons reported by Nextclade. Treat B as fine for most surveillance use, C as "look at the sample" and D as not suitable for phylogenetics.
- **`coverage`** — should agree with `coverage_min_depth` in `summary.tsv` (there in percent).
- **`inputSequenceStatus`** — empty for analysed sequences; set when viralQC could not analyse a record (e.g. all `N`).

`samples/<sample>/viralqc.tsv` contains the header plus that sample's rows. Pass `--no-run-viralqc` to skip the step (for instance on a node without the databases or without internet access: `nextclade sort` fetches a small index from the Nextclade server on every run); if viralQC itself fails, the run still completes and `qc/viralqc/viralqc_status.txt` says why. Because the step is then considered done, delete `qc/viralqc/` and rerun the same command to retry it after fixing the cause.

### Optional: intra-host SNVs

If you care about within-sample variation (e.g., looking for sub-consensus mutations), add `--run-isnv`:

```bash
viralconseq consensus illumina \
    --sample-sheet samples_illumina.csv \
    --config-file  results/consensus_illumina_isnv/config.yml \
    --output       results/consensus_illumina_isnv/ \
    --reference    databases/refs/nCoV-2019.reference.fasta \
    --primer-scheme databases/refs/nCoV-2019.bed \
    --run-isnv \
    --af-isnv-threshold 0.05 \
    --threads 2 --threads-total 4
```

This adds a LoFreq pass on the primer-clipped BAM. The output VCFs land under `assembly/isnvs/<sample>.isnvs.vcf.gz`, and a per-sample tally is written to `isnvs/isnvs_summary.tsv`. `--af-isnv-threshold` sets the minimum allele frequency reported (LoFreq filters consensus-level variants out so this captures the sub-consensus range).

## Worked example 2 — Nanopore (SARS-CoV-2)

Same virus, different chemistry. The reads are in `my_test_data/nanopore_data/`. Nanopore's higher per-base error rate calls for a deep-learning variant caller (Clair3) instead of the pileup-based consensus used for Illumina.

### Sample sheet (2 columns)

```text
barcode05,my_test_data/nanopore_data/barcode05.itps-0003.fastq.gz
barcode09,my_test_data/nanopore_data/barcode09.itps-0004.fastq.gz
```

### Run

```bash
viralconseq consensus nanopore \
    --sample-sheet samples_nanopore.csv \
    --config-file  results/consensus_nanopore/config.yml \
    --output       results/consensus_nanopore/ \
    --run-name     sarscov2 \
    --reference    databases/refs/nCoV-2019.reference.fasta \
    --clair3-model r941_prom_hac_g360+g422 \
    --minimum-coverage 20 \
    --minimum-map-quality 30 \
    --threads 4 --threads-total 4
```

The interesting differences from the Illumina invocation:

- No `fastp` step — the workflow skips QC entirely for Nanopore. The variant caller is expected to absorb noisy bases.
- `--clair3-model` names the Clair3 model matching your basecaller. Leave it at its default, `auto`, for recent data: the pipeline reads the basecall model tag Dorado writes into every read header and picks the model for you, per sample. The tutorial reads are 2022 Guppy R9.4.1 reads that carry no such tag, so `auto` would stop with a message asking for a model, and `r941_prom_hac_g360+g422` is the right one for them. Whatever the model, it must be present in `~/.cache/viralconseq/clair3-models` (or `--clair3-model-dir`); `viralconseq setup` fetches the common R10.4.1 hac/sup models and this R9 one by default, and `viralconseq setup --clair3-models NAME` fetches any other name from the [Clair3 model zoo](https://github.com/HKU-BAL/Clair3#pre-trained-models).
- `--minimum-map-quality 30` filters reads with MAPQ below 30 before variant calling. Tighten this for very noisy runs.
- The pipeline silently sanitizes the reference's FASTA headers (replacing `/`, `|`, `,`, `~`, and spaces with `_`) before use — clair3 turns the seq ID into a directory name, so the sanitization avoids cryptic filesystem errors. The sanitized copy lands at `results/consensus_nanopore/sarscov2/reference/reference.sanitized.fasta`.

### What is different in the output

Most of the output layout matches the Illumina run, with a couple of additions and substitutions:

- `reference/reference.sanitized.fasta` — the sanitized reference used by every downstream rule. If you need the exact contig names that appear in the BAM and VCF, look here, not at your input FASTA.
- `samples/<sample>/consensus.vcf.gz` — produced by Clair3 + `bcftools` rather than from the consensus directly. You will see `QUAL` values reflective of clair3's confidence.
- No `fastp.html` per sample, no `multiqc_report.html`.

The headline `consensus.fasta` and `samples_alignment.fasta` files behave identically to the Illumina case.

```{tip}
`--af-threshold` means something slightly different on each platform. On Illumina it is the fraction of base counts at the position (`samtools consensus -c`) and defaults to `0.51`. On Nanopore it is ALT reads over REF plus ALT reads from Clair3's `AD` and defaults to `0.6`; reads showing a deletion at the site do not count, so a variant that creates a homopolymer is not lost to the reads that collapse it. Clair3 additionally filters calls by its own quality score (`--variant-quality`, `--variant-depth`). If you see too many heterozygous-looking sites in your consensus, raise the threshold: `--af-threshold 0.7` or higher.
```

## Segmented viruses

For multi-segment genomes (influenza, the bunyaviruses, etc.) the simplest way is
to pass a single multi-record FASTA to `--reference`: viralconseq treats a
reference with more than one record as segmented and splits it into one reference
per record, taking segment names from the FASTA headers. Pass `--single-reference`
to opt out (align all records together as one reference — e.g. a fragmented genome
that should be treated as a single reference).

```bash
viralconseq consensus illumina \
    --sample-sheet samples_flu.csv \
    --config-file  results/flu/config.yml \
    --output       results/flu/ \
    --reference    refs/flu/influenza_all_segments.fasta \
    --threads 4 --threads-total 8
```

To name segments yourself, or supply each from a separate file, repeat
`--segmented-reference SEG=PATH` once per segment instead (mutually exclusive with
`--reference`):

```bash
viralconseq consensus illumina \
    --sample-sheet samples_flu.csv \
    --config-file  results/flu/config.yml \
    --output       results/flu/ \
    --segmented-reference HA=refs/flu/HA.fasta \
    --segmented-reference NA=refs/flu/NA.fasta \
    --segmented-reference PB1=refs/flu/PB1.fasta \
    # … one per segment …
    --threads 4 --threads-total 8
```

The workflow dispatches a specialised segmented variant that processes each segment in parallel. Outputs are keyed by segment under `assembly/<segment>/consensus/final_consensus/` and per-sample symlinks land under `samples/<sample>/<segment>/`.

See [Notes — Segmented viruses](../notes.md#segmented-viruses) for the full discussion.

## Tuning consensus quality

The three knobs you will reach for most often:

- **`--minimum-coverage`** (default `20`) — any reference position with fewer than this many reads after primer clipping becomes `N`. Lower it (e.g. `10`) for low-depth samples where you would rather see a tentative base than a hole; raise it for high-confidence assemblies.
- **`--af-threshold`** (default `0.51` on Illumina, `0.6` on Nanopore) — minimum allele fraction to call a variant into the consensus. On Illumina it is a fraction of base counts, so `0.51` is the majority allele; on Nanopore it is ALT over REF plus ALT reads. Raise it (e.g. `0.7`) for noisier data; lower it to capture ambiguity codes.
- **`--minimum-read-length`** (default `50`) — on Illumina, reads below this length are dropped by fastp at QC time and by `samtools ampliconclip --filter-len` after primer clipping; on Nanopore it is applied only by `samtools ampliconclip --filter-len`, i.e. only when a primer scheme is given.

Nanopore has a few extra knobs:

- **`--variant-quality`** (default `15`) — Clair3 minimum QUAL. Calls below it are flagged `LowQual` and never reach the consensus.
- **`--variant-depth`** (default `10`) — minimum alt-allele read support.
- **`--chunk-size`** (default `10000`) — Clair3 chunk size; tune only if Clair3 runs out of memory on huge references.

The full set is in the [Commands reference](../commands.md#viralconseq-consensus).

## Common pitfalls

- **Sample sheet columns must match the subcommand.** `consensus illumina` requires 3 columns and `consensus nanopore` 2; a mismatch is rejected with a `SampleSheetError` that names the offending row.
- **Reference header sanitization is silent.** On Nanopore, if your downstream tooling expects the exact contig names used in the BAM/VCF, read them from `reference/reference.sanitized.fasta` (or `reference/<segment>.sanitized.fasta` in segmented runs), not your input FASTA.
- **Primer scheme contig names must match the reference.** The BED file must use the same chromosome/contig names as the reference FASTA. viralconseq checks this before running and aborts with an `InputIntegrityError` if the primer BED chrom matches no reference contig (so `samtools ampliconclip` would clip nothing) — no more silently un-clipped runs.
- **Inputs are content-validated before the run.** FASTQ, reference FASTA, and primer BED are streamed and checked (record structure, gzip integrity, nucleotide alphabet, chrom matching); a broken, truncated, or mismatched file stops the run up front rather than failing deep inside Snakemake. Pass `--skip-input-validation` to bypass.
- **`--reference` and `--segmented-reference` are mutually exclusive.** Pass one or the other.
- **`viralqc_database_not_found` at start-up.** viralQC is on by default and its databases live in `~/.cache/viralconseq/viralqc-db` (or `$VIRALCONSEQ_VIRALQC_DB`). Run `viralconseq setup` once (see [Setup §2](setup.md#2-build-per-rule-environments-download-the-viralqc-databases-and-the-clair3-models)) or pass `--no-run-viralqc`.

## Reference

- [Commands reference — `viralconseq consensus`](../commands.md#viralconseq-consensus): the complete flag table.
- [Output layout](../output.md): the full output tree.
- [Notes — Segmented viruses](../notes.md#segmented-viruses): more on multi-segment runs and the segmented output layout.
