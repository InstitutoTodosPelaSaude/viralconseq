# Workflow Reference

What the Snakemake workflows actually run. For every rule: the conda
environment it runs in, the commands, each flag and why it is set, the threads
and memory the rule declares, and the peak memory measured on the test fixture.
Read this when you need to know which tool and settings produced an output
file. The files themselves are described in [Output Layout](output.md); the
Python layer that writes the config and launches Snakemake is described in
[Architecture](architecture.md).

Everything below is read from `viralconseq/scripts/*.smk`,
`viralconseq/scripts/rules/*.smk`, `viralconseq/scripts/envs/*.yaml` and
`viralconseq/scripts/python/*.py`. Where a flag's purpose is not stated in the
rule file and is taken from the tool's documentation or from domain knowledge,
the reason is marked *(inferred)*.

## The four workflows

| Entry file | Reads | Reference |
|---|---|---|
| `consensus_illumina.smk` | Illumina paired-end | one FASTA (one or more contigs) |
| `consensus_illumina_segmented.smk` | Illumina paired-end | one FASTA per segment |
| `consensus_nanopore.smk` | Nanopore | one FASTA (one or more contigs) |
| `consensus_nanopore_segmented.smk` | Nanopore | one FASTA per segment |

Each entry file is short and has the same shape:

1. `include: "rules/common.smk"` **first**. Global `wildcard_constraints` only
   apply to rules parsed after them, so the shared preamble must precede every
   rule. It provides the required-key guard (`_REQUIRED_KEYS`; a hand-edited
   YAML missing an analysis key fails at parse time), the full config check
   through `viralconseq.validators.validate_config_dict` when the package is
   importable, the shell-metacharacter check on the two `*_flags` strings, the
   `sample` wildcard constraint (a sample wildcard matches only the configured
   ids, so ids that are prefixes of one another never match ambiguously), the
   `LOG()` / `BENCH()` / `cpus()` / `ram_mb()` helpers, the Clair3 model
   helpers, and the `onstart` / `onsuccess` / `onerror` hooks that write
   `logs/run.log`.
2. `TERMINAL_INPUTS`, the list of analysis products, followed by `rule all`,
   which is the **default target** (`default_target: True`). `rule all` asks
   for `TERMINAL_INPUTS` plus `versions.tsv`, `config.yml` and `benchmark.tsv`.
3. The rule modules, included in DAG order (table below), then the
   entry-file rules `organize_files` (all four) and `summarize_isnvs`
   (Illumina), and finally `include: "rules/provenance.smk"`.
4. `collect_benchmarks` (in `rules/provenance.smk`) is the **terminal rule**:
   it depends on `TERMINAL_INPUTS`, the `organize_files` sentinel,
   `versions.tsv` and `config.yml`, and writes `benchmark.tsv`. Nothing runs
   after it, which is why its own benchmark is the only one missing from the
   table it writes.

`TERMINAL_INPUTS` and the conditions under which each entry is requested:

| Product | Condition |
|---|---|
| `assembly/[<segment>/]consensus/final_consensus/samples_alignment.fasta` | always (one per segment) |
| `isnvs/isnvs_summary.tsv` | Illumina and `run_isnv` |
| `qc/viralqc/outputs/results.tsv` | `run_viralqc` (default on) |
| `summary.tsv` | always |
| `consensus/consensus[.<segment>].cov<T>.fasta` | always; `T` is `consensus_coverage_threshold` |
| `report.html` | `run_report` (default on) |

Rule modules in include order:

| Module | Rules | Included by |
|---|---|---|
| `rules/common.smk` | (helpers, guards, hooks; no rules) | all four |
| `rules/qc_illumina.smk` | `perform_qc` | Illumina |
| `rules/alignment_illumina.smk` / `alignment_nanopore.smk` | `map_reads`, `trim_primer_sequences` | per platform |
| `rules/consensus_illumina.smk` | `detect_isnv`, `infer_consensus_sequence`, `generate_vcf_consensus` | Illumina |
| `rules/consensus_nanopore.smk` | `check_mapped_reads`, `infer_consensus_sequence` | Nanopore |
| `rules/stats.smk` | `calculate_coverage_basewise`, `rename_sequences` | all four |
| `rules/consensus_illumina_common.smk` | `calculate_assembly_statistics`, `generate_multiqc_report`, `align_consensus_to_reference_genome` | Illumina |
| `rules/consensus_nanopore_common.smk` | `calculate_assembly_statistics`, `align_consensus_to_reference_genome` | Nanopore |
| `rules/viralqc.smk` | `prepare_viralqc_input`, `run_viralqc`, `split_viralqc_results` | all four |
| `rules/collect.smk` | `summary`, `collect_consensus` | all four |
| `rules/report.smk` | `report` | all four |
| (entry file) | `organize_files`; Illumina also `summarize_isnvs`; Nanopore also `sanitize_reference` | all four |
| `rules/provenance.smk` | `versions_<env>`, `versions`, `run_config`, `collect_benchmarks` | all four (last) |

**Segmented variants.** The segmented entry files set
`SEGMENT_WILDCARD = "{segment}/"`, which every per-sample rule splices into its
output path (`assembly/{segment}/...`), and make `REFERENCE` a per-segment
value (a function of the `segment` wildcard on Illumina; the per-segment
sanitised reference on Nanopore). Per-sample rules therefore run once per
sample **and segment**; run-level rules run once. The single-reference entry
files set `SEGMENT_WILDCARD = ""` and the rule modules are unchanged. The
per-segment FASTAs come from the Python layer (`reference_splitter.py`, written
to `input_references/`), never from the workflow.

**How Snakemake is launched.** `_orchestrator.run_workflow` calls
`snakemake()` with the config file, `cores=threads_total`, `use_conda=True`,
`workdir=<output>/<run_name>` (so every path in the YAML is absolute and only
`.snakemake/` moves), `force_incomplete=True` (an interrupted run resumes),
`targets=["all"]`, and `resources={"mem_mb": <budget>}` when the memory budget
is positive (see [Resources](#resources)).

## Reading the tables

- **Env** is `envs/<name>.yaml` under `viralconseq/scripts/`. Tool versions are
  pinned there and probed at run time into `versions.tsv`:

  | Env | Pinned tools |
  |---|---|
  | `alignment` | minimap2 2.30, samtools 1.23, bedtools 2.31, GSAlign 1.0, gofasta 1.2, deacon 0.15, Python 3.11 |
  | `qc` | fastp 1.1, MultiQC 1.21, Python 3.11 |
  | `consensus` | LoFreq 2.1, bcftools 1.23, samtools 1.23, racon 1.5, gofasta 1.2, Python 3.11 |
  | `clair3` | Clair3 2.0.2, samtools 1.23, bcftools 1.23, Python 3.11 |
  | `utils` | Python 3.11, samtools 1.23, bcftools 1.23, seqtk 1.5 |
  | `viralqc` | viralQC 1.2.0 (PyPI), Nextclade 3.15, BLAST 2.16, seqtk 1.5, snakemake-minimal 7.32.4, pulp < 2.8, pandas 2.2, ncbi-datasets-cli 18.9, TaxonKit 0.20 |

  `bgzip` and `tabix` come with the htslib that samtools and bcftools pull in.
  `run_config` has no env: it is a `run:` block executed by the Snakemake
  process itself.
- **Command(s)** summarise the shell body; the exact text is in the rule file.
  Parameters in braces are config keys (`{minimum_depth}` is
  `config["minimum_depth"]`, set by `--minimum-coverage`).
- **Threads / mem_mb**: `cpus(<rule>)` means the rule has a `--<rule>-cpus`
  option and reads `<rule>_cpus`, else `threads`; "1" means the rule declares no
  `threads:` and Snakemake gives it one core. Only two rules declare `mem_mb`.
- **Measured max_rss** is Snakemake's `max_rss` column of `benchmark.tsv`,
  rounded to whole MB, taken as the maximum over the rows of that rule in two
  real runs of the SARS-CoV-2 fixture: 2 Illumina samples of about 200k read
  pairs each, and 2 Nanopore samples of 8–17k reads each (plus two near-empty
  control samples), both with `--threads 2`. Cells read *Illumina / Nanopore*;
  a single figure means the rule exists on one platform only. These are
  small-input figures: Clair3 and viralQC in particular scale with input, which
  is why their declared memory carries headroom.

## Read QC (Illumina)

| Rule | Env | Command(s) | Flags and why | Threads / mem_mb | Measured max_rss |
|---|---|---|---|---|---|
| `perform_qc` | `qc` | `fastp -i R1 -I R2 -o trim.p.R1 -O trim.p.R2 --unpaired1 trim.u.R1 --unpaired2 trim.u.R2 ... --json --html --thread {threads} <adapter arg>` | `--length_required {minimum_length}`: drop reads shorter than `--minimum-read-length` after trimming. `--trim_front1/2 {trim_head}`, `--trim_tail1/2 {trim_tail}`: fixed head/tail trimming of both mates. `--cut_front` + `--cut_front_mean_quality`, `--cut_tail` + `--cut_tail_mean_quality`: sliding-window quality trimming from the 5' and 3' ends (window size fastp default). `--cut_right` + `--cut_right_window_size` + `--cut_right_mean_quality`: trim from the first low-quality window to the read end. `--adapter_fasta <adapters>` when `--adapters` was given, else `--detect_adapter_for_pe` (auto-detect from mate overlap). | `cpus(perform_qc)` | 88 |
| `generate_multiqc_report` | `qc` | `multiqc -f --cl-config "extra_fn_clean_exts: ['_R1']" -o qc/reports/ qc/reports/` | `-f`: overwrite a report from an earlier attempt. `--cl-config extra_fn_clean_exts`: strip `_R1` so the sample name matches the sheet *(inferred)*. Inputs are the per-sample fastp JSONs so the report is rebuilt when one changes. | 1 | 109 |

Only the paired survivors (`qc/data/trim.p.<sample>_R{1,2}.fastq.gz`) go on to
`map_reads`; the unpaired files are kept but unused. `generate_multiqc_report`
is a run-level rule and is in the DAG only because the Illumina
`align_consensus_to_reference_genome` lists its report as an input; nothing in
`rule all` asks for it directly. Nanopore reads are not trimmed or filtered
before mapping.

## Reference sanitisation (Nanopore)

| Rule | Env | Command(s) | Flags and why | Threads / mem_mb | Measured max_rss |
|---|---|---|---|---|---|
| `sanitize_reference` | `clair3` | `sed '/^>/s/[\/|,~ ]/_/g' <reference> > reference/reference.sanitized.fasta; samtools faidx` | The `sed` rewrites `/ \ \| , ~` and spaces in header lines to `_`: Clair3 creates one directory per contig named after the sequence id. `faidx`: the `.fai` is a declared input of `infer_consensus_sequence` (contig lengths for the all-N branch). | 1 | 2 |

Segmented runs run it once per segment (`reference/<segment>.sanitized.fasta`,
log target `<segment>`). Every downstream Nanopore rule uses the sanitised
copy as `REFERENCE`; the Illumina workflows use the reference as given. The
Python layer keeps `integrity.sanitize_nanopore_contig` in step with this
`sed`, so a primer BED whose chromosome names differ only by these characters
still validates.

## Alignment and primer clipping

| Rule | Env | Command(s) | Flags and why | Threads / mem_mb | Measured max_rss |
|---|---|---|---|---|---|
| `map_reads` (Illumina) | `alignment` | `minimap2 -a -t {threads} -x sr <ref> R1 R2`, piped into `samtools view -bS -F 4 -`, into `samtools sort -o raw/<sample>.sorted.bam -`; then `samtools index` | `-a`: SAM output. `-x sr`: short-read preset; two FASTQs make minimap2 treat them as pairs. `-F 4`: drop unmapped reads. `sort`: coordinate order for indexing and every downstream tool. | `cpus(map_reads)` | 196 |
| `map_reads` (Nanopore) | `alignment` | as above with `-x map-ont` and `samtools view --min-MQ {minimum_map_quality} -bS -F 4 -` | `-x map-ont`: ONT preset. `--min-MQ`: drop alignments below `--minimum-map-quality` (default 30) before counting or calling. | `cpus(map_reads)` | 55 |
| `trim_primer_sequences` (both) | `alignment` | `scheme` set: `samtools ampliconclip --both-ends --hard-clip --filter-len {minimum_length} -b <scheme.bed> -f trimmed/<sample>.trimmed.txt raw.bam > trimmed/<sample>.sorted.bam; samtools index`. `scheme` is `NA`: `cp` BAM and index to `trimmed/`, `touch` an empty `.trimmed.txt`, write `raw/notes.txt` saying the data are treated as untargeted. | `--both-ends`: clip primers at both ends of a read, not only the 5' end. `--hard-clip`: remove the primer bases instead of soft-clipping, so no caller can count them *(inferred)*. `--filter-len`: discard reads shorter than `--minimum-read-length` after clipping. `-b`: the primer BED (`--primer-scheme`). `-f`: per-read clipping report. | `cpus(trim_primer_sequences)` | 7 / 5 |

`trim_primer_sequences` declares threads through `cpus()` (it has a
`--trim-primer-sequences-cpus` option) but the shell does not pass `{threads}`
to `samtools ampliconclip`; the value only limits how many of these jobs run at
once. Both BAMs are kept: `raw/` (`raw_mapped_reads.bam` in the browse tree)
and `trimmed/` (`trimmed_mapped_reads.bam`), and every rule after this one
reads the trimmed BAM, whether or not a scheme was given.

## Mapped-read gate (Nanopore)

| Rule | Env | Command(s) | Flags and why | Threads / mem_mb | Measured max_rss |
|---|---|---|---|---|---|
| `check_mapped_reads` | `alignment` | `n=$(samtools view -c -F 260 trimmed.bam)`; writes `assembly/status/<sample>[.<segment>].txt` with `status` (`ok` or `no_mapped_reads`), `mapped_reads` and `minimum_mapped_reads` | `-c`: count only. `-F 260`: exclude unmapped (4) and secondary (256) records, i.e. count primary mapped reads. `status` is `no_mapped_reads` when `n < minimum_mapped_reads` (`--minimum-mapped-reads`, default 10). | 1 | 5 |

The status file is a declared input of the Nanopore `infer_consensus_sequence`,
which branches on it, and of `summary` and `report`. There is no equivalent
gate on Illumina: `samtools consensus` tolerates an empty BAM, and
`generate_vcf_consensus` handles the resulting all-N consensus itself.

## Variant calling and consensus

### Illumina

| Rule | Env | Command(s) | Flags and why | Threads / mem_mb | Measured max_rss |
|---|---|---|---|---|---|
| `infer_consensus_sequence` | `alignment` | `samtools consensus -a -d {minimum_depth} -m simple -q -c {af_threshold} --show-ins yes trimmed.bam -o final_consensus/<sample>.consensus.fasta` | `-a`: emit every reference position, uncovered ones as `N`. `-d`: positions below `--minimum-coverage` become `N`. `-m simple`: frequency-based calling. `-q`: weight the base counts by base quality *(inferred)*. `-c`: an allele is called when it holds at least `--af-threshold` (default 0.51) of the counts. `--show-ins yes`: keep insertions relative to the reference. | 1 | 8 |
| `generate_vcf_consensus` | `alignment` | If the consensus has any A/C/G/T (an `awk` test): `GSAlign -r <ref> -q <consensus> -o <prefix> -fmt 1 -sen`, then `bgzip` + `tabix -p vcf` on `<prefix>.vcf` and `rm <prefix>.maf`. Otherwise, or when GSAlign exits non-zero or writes an empty VCF: a header-only "mock" VCF, and a one-line `WARNING` on the console. | `-fmt 1`: MAF alignment output, deleted afterwards (only the VCF is wanted). `-sen`: sensitive mode, for a consensus divergent from the reference *(inferred)*. The gate exists because GSAlign produces no VCF for an all-N query and the run would otherwise abort under `set -euo pipefail`. | 1 | 2 |
| `detect_isnv` (`run_isnv` only) | `consensus` | `lofreq indelqual --dindel -f <ref> -o isnvs/<sample>.lofreq.sorted.bam trimmed.bam; samtools index`; `lofreq call-parallel --pp-threads {threads} --call-indels -f <ref> -o tmp.vcf lofreq.bam`; `bcftools view -i 'INFO/AF<0.5 & INFO/AF>={af_isnv_threshold}' tmp.vcf -Oz -o isnvs/<sample>.isnvs.vcf.gz; tabix` | `indelqual --dindel`: insert indel base qualities with the Dindel model, required before `--call-indels`. `--pp-threads`: parallel calling over regions. `--call-indels`: report indels as well as SNVs. The `bcftools view` filter keeps only sub-consensus alleles (`AF < 0.5`) at or above `--af-isnv-threshold` (default 0); LoFreq's own default filters still apply. The re-qualified BAM and the raw VCF are `temp()`. | `cpus(detect_isnv)` | 74 |
| `summarize_isnvs` (entry file; `run_isnv` only) | `utils` | for each sample: `bcftools view -H <isnvs.vcf.gz> \| wc -l`, appended to `isnvs/isnvs_summary.tsv` (`sample`, `[segment,]`, `number_of_isnvs`) | `-H`: records only, no header, so the line count is the variant count. | 1 | 1 |

The Illumina `consensus.vcf.gz` is therefore derived from the *consensus
sequence* (GSAlign aligns the consensus back to the reference); it describes
what the consensus differs in, not read-level evidence. Read-level minority
variants are the separate iSNV table.

### Nanopore

| Rule | Env | Command(s) | Flags and why | Threads / mem_mb | Measured max_rss |
|---|---|---|---|---|---|
| `infer_consensus_sequence` | `clair3` | **Gate**: if the status file says `no_mapped_reads`, write one all-`N` record per reference contig (lengths from the `.fai`), header-only `raw`, `norm` and final VCFs (with `##contig` lines), an empty `low_cov.bed`, a `model.txt` with `status no_mapped_reads`, print a `WARNING`, `exit 0`. **Otherwise**: `run_clair3.sh --bam_fn=trimmed.bam --ref_fn=<ref> --model_path=<clair3_model_dir>/<model> --output=assembly/[<segment>/]clair3/<sample> --platform=ont --threads={threads} --chunk_size={chunk_size} --qual={variant_quality} --min_mq={minimum_map_quality} --enable_long_indel --haploid_sensitive --no_phasing_for_fa --include_all_ctgs`; decode exit codes 1 / 2 / 127; `cp merge_output.vcf.gz(.tbi)` to `<sample>.raw.vcf.gz`; `bcftools norm -m - -f <ref> raw.vcf.gz > <sample>.norm.vcf.gz`; `bcftools filter -i 'FILTER="PASS" && FORMAT/AF >= {af_threshold} && FORMAT/AD[0:1] >= {variant_depth}' norm -O z -o <sample>.vcf.gz; tabix`; `samtools depth -J -a trimmed.bam \| awk '$3 <= {minimum_depth}'` into `<sample>.low_cov.bed`; `bcftools consensus -f <ref> --mask low_cov.bed <sample>.vcf.gz > <sample>.consensus.fasta`; `find <clair3 dir> ... -exec rm -rf` everything but the raw VCF, its index and `model.txt`. | `--platform=ont`: ONT models and error profile. `--threads`: Clair3's parallelism. `--chunk_size` (`--chunk-size`, default 10000 bp): size of the regions called in parallel. `--qual` (`--variant-quality`, default 20): calls below this QUAL are flagged `LowQual`, then dropped by `FILTER="PASS"`. `--min_mq` (`--minimum-map-quality`): ignore alignments below this mapping quality. `--enable_long_indel`: call indels longer than 50 bp *(inferred from Clair3 docs)*. `--haploid_sensitive`: haploid calling in which both 0/1 and 1/1 count as variants, so mixed sites reach the AF filter instead of being dropped *(inferred)*. `--no_phasing_for_fa`: skip WhatsHap phasing before full-alignment calling; there is nothing to phase in a haploid genome *(inferred)*. `--include_all_ctgs`: call every contig (Clair3 otherwise restricts itself to human chromosome names). Boolean flags come last because Clair3 2.x parses them as `nargs="?"` and would swallow the next token. `bcftools norm -m -`: split multi-allelic records; `-f`: left-align and normalise indels against the reference. The `filter` keeps PASS calls whose allele fraction is at least `--af-threshold` and whose ALT read count (`AD[0:1]`) is at least `--variant-depth` (default 10). `samtools depth -a`: every position; `-J`: reads with a deletion at the position still count, so a real deletion is not masked *(inferred)*; the `awk` masks positions at or below `--minimum-coverage`. `bcftools consensus --mask`: apply the filtered variants and write `N` over the masked intervals. The model checkpoints `pileup.pt` and `full_alignment.pt` are declared rule inputs, so a missing model is a `MissingInputException` at DAG time rather than a Clair3 crash. | `cpus(infer_consensus_sequence)`; `mem_mb = ram_mb(infer_consensus_sequence)` (default 2 GB) | 416 |

`raw.vcf.gz` is Clair3's `merge_output.vcf.gz` untouched; `consensus.vcf.gz`
is what survived normalisation and the filter and what `bcftools consensus`
applied. `clair3/<sample>/model.txt` records the model name and directory. The
model is `config["clair3_model"]`, either one name or a `{sample: name}`
mapping when `--clair3-model auto` detected different basecallers per sample
(`clair3_model_for()` in `rules/common.smk`).

## Coverage statistics and renaming

| Rule | Env | Command(s) | Flags and why | Threads / mem_mb | Measured max_rss |
|---|---|---|---|---|---|
| `calculate_coverage_basewise` (both) | `alignment` | `bedtools genomecov -d -ibam trimmed.bam > coverage_stats/<sample>.table_cov_basewise.txt` | `-d`: depth at every reference position (1-based), zero-depth positions included. `-ibam`: BAM input. | 1 | 8 / 8 |
| `rename_sequences` (both) | `utils` | `script:` `python/rename_sequences.py` — rewrites the header(s) of `<sample>.consensus.fasta` into `<sample>.consensus.renamed.fasta` | One record: header becomes `sample-<id>` (from the file name). Several records (a multi-contig reference): one record per contig headed `sample-<id>_<contig>`, so contigs are never fused into one chimeric sequence. Sequence bodies are copied verbatim. | 1 | 1 / 23 |
| `calculate_assembly_statistics` (both) | `utils` | `script:` `python/calculate_assembly_stats.py` — one headered row in `coverage_stats/<sample>.stats.tsv` | Reads: `total_reads` = records in the input FASTQ(s) (both mates on Illumina); `qc_passed_reads` from the fastp paired outputs (Illumina only, `NA` on Nanopore); `mapped_reads` = `samtools view -c -F 260` on the trimmed BAM. Depth and breadth from the basewise table: mean and median depth, percent of positions at or above 10x / 100x / 1000x / `minimum_depth` (`coverage_min_depth`). Consensus length, N count and N percent from the renamed FASTA. | 1 | 33 / 33 |

`calculate_assembly_statistics` and the two Python `script:` rules run in
`envs/utils.yaml`, where the `viralconseq` package is not installed; the
scripts are standard-library only and redirect their own stdout/stderr into
the rule log. The `coverage_min_depth` figure counts positions with depth
**at or above** `minimum_depth`.

## Multi-sample alignment

| Rule | Env | Command(s) | Flags and why | Threads / mem_mb | Measured max_rss |
|---|---|---|---|---|---|
| `align_consensus_to_reference_genome` (both; run-level, once per segment) | `alignment` | `cat <ref> <every renamed consensus> > final_consensus/consensus.fasta`; `minimap2 {minimap2_consensus_align_flags} <ref> consensus.fasta -o aln.consensus.sam`; one reference contig: `gofasta sam toMultiAlign --pad -s aln.consensus.sam -o samples_alignment.fasta`; several contigs: split the SAM per contig (`@HD`, that contig's `@SQ`, `@PG`, records with that `RNAME`), skip contigs with no records, run `gofasta sam toMultiAlign --pad` per contig into `per_contig_alignments/<contig>.fasta` and concatenate into `samples_alignment.fasta`; finally `sed '/^>/ ! s/-/N/g'` into `aln.consensus.indelsMasked.fasta` | `minimap2_consensus_align_flags` is a config-only key, default `-a --sam-hit-only --secondary=no --score-N=0`: `-a` SAM; `--sam-hit-only` no records for unmapped queries (an all-N consensus); `--secondary=no` one alignment per consensus; `--score-N=0` runs of `N` in a consensus are not scored as mismatches *(inferred)*. `gofasta sam toMultiAlign --pad`: a reference-coordinate alignment (insertions relative to the reference are dropped, deletions become `-`), padded with `N` to the reference length. gofasta aborts on a multi-contig reference and on a header-only SAM, hence the per-contig loop and the skip. The `sed` turns deletion gaps into `N` for tools that cannot read `-`. | 1 | 1 / 3 |

The combined input starts with the reference, so the reference is the first
record of `samples_alignment.fasta`. On Illumina this rule also lists the
MultiQC report as an input, which is what pulls `generate_multiqc_report` into
the DAG.

## Consensus QC (viralQC; `run_viralqc`, default on)

| Rule | Env | Command(s) | Flags and why | Threads / mem_mb | Measured max_rss |
|---|---|---|---|---|---|
| `prepare_viralqc_input` (run-level) | `utils` | for every renamed consensus of every sample (and segment): `awk` rewrites each header to `sample-<id>[\|<contig>][\|<segment>]` and appends the record to `qc/viralqc/input.fasta`; warns on a file with no records; fails on duplicate headers | Headers are normalised so `\|` is the only separator and a sample id can be matched unambiguously in viralQC's `seqName`. A `\|` inside a contig or segment name is replaced by `_`. The reference is **not** included. Zero records in total is a warning, not an error: `run_viralqc` then writes an empty table. | 1 | 2 / 1 |
| `run_viralqc` (run-level) | `viralqc` | `rm -rf qc/viralqc/outputs qc/viralqc/.snakemake` (a stale nested-Snakemake state would make `vqc` report "Nothing to be done"); if `input.fasta` has no records: header-only `results.tsv`, `viralqc_status.txt` with `status skipped`, `exit 0`; otherwise `env -u SNAKEMAKE_PROFILE vqc run --input input.fasta --output-dir <abs qc/viralqc> --output-file results.tsv --datasets-dir <abs viralqc_db> --blast-database <db>/blast.fasta --blast-database-metadata <db>/blast.tsv --cores {threads} --verbose {viralqc_extra_flags}`; then classify the outcome | `env -u SNAKEMAKE_PROFILE`: the nested Snakemake must not pick up the operator's profile. `--output-dir` and `--datasets-dir` are made absolute because `vqc` hands them to `snakemake --directory`. `--cores {threads}`: viralQC's own parallelism. `viralqc_extra_flags` is a config-only key, empty by default. Declared inputs are `blast.fasta`, `blast.tsv` and `.nextclade_datasets_ok` of the database directory (`--viralqc-db`); the rest of the layout is checked in Python before the run. **Lenient**: `status` is `ok` (exit 0, table written), `partial` (non-zero exit but a table) or `failed` (no table: a placeholder with one row per input record and `inputSequenceStatus` = `viralQC failed (exit N)` is written); the run continues in every case and a `WARNING` goes to the console. | `cpus(run_viralqc)`; `mem_mb = ram_mb(run_viralqc)` (default 1 GB) | 349 / 273 |
| `split_viralqc_results` (per sample) | `utils` | `awk -F'\t' 'NR==1 \|\| $1==s \|\| index($1, s "\|")==1' results.tsv > qc/viralqc/per_sample/<sample>.viralqc.tsv` | Header plus the rows whose `seqName` is the sample or starts with `sample\|` (contig and segment suffixes). Works on the placeholder table too, so `organize_files` never blocks on a viralQC failure. | 1 | 1 / 2 |

`rules/viralqc.smk` refuses to parse when `run_viralqc` is on and `viralqc_db`
is unset or `NA`. With `--no-run-viralqc` none of the three rules is in the
DAG, `versions_viralqc` is not probed, and `summary` / `report` receive no
viralQC inputs.

## Collection

| Rule | Env | Command(s) | Flags and why | Threads / mem_mb | Measured max_rss |
|---|---|---|---|---|---|
| `summary` (run-level) | `utils` | `python python/build_summary.py --stats <every .stats.tsv> --samples <ids> [--segments <keys>] --min-depth {minimum_depth} --data-type <illumina\|nanopore> [--viralqc-results results.tsv --viralqc-status viralqc_status.txt] [--isnvs isnvs/isnvs_summary.tsv] [--status-files ...] [--model-files ...] --output summary.tsv --legacy-csv assembly/assembly_stats_summary.csv` | Joins the per-sample statistics with viralQC (`run_viralqc`), iSNV counts (Illumina and `run_isnv`), the mapped-read status and Clair3 model (Nanopore) into the pinned `build_summary.COLUMNS` header, one row per sample (per sample and segment when segmented), and derives the `status` column. Every configured sample gets a row (`missing_stats` when its statistics are absent). The legacy CSV is the deprecated pre-0.2.0 shape. | 1 | 10 / 10 |
| `collect_consensus` (run-level, once per segment) | `utils` | `python python/collect_consensus.py --consensus <renamed FASTAs> --summary summary.tsv --threshold {consensus_coverage_threshold} [--segment <key>] --samples <ids> --output-dir consensus` | Writes `consensus/<sample>[.<segment>].fasta` (one-line sequences, headers `sample-<id>[\|<contig>][\|<segment>]`), `consensus/consensus[.<segment>].fasta` (all samples pooled, no reference) and `consensus/consensus[.<segment>].cov<T>.fasta` with the samples whose `coverage_min_depth` in `summary.tsv` is at or above `T` (`--consensus-coverage-threshold`, default 70); a sample with `NA` coverage is excluded and named in the log. | 1 | 9 / 11 |

## Report (`run_report`, default on)

| Rule | Env | Command(s) | Flags and why | Threads / mem_mb | Measured max_rss |
|---|---|---|---|---|---|
| `report` (run-level) | `utils` | `python python/build_report.py --run-dir <run> --template templates/report.html --output report.html --data-type <illumina\|nanopore> --min-depth {minimum_depth} --consensus-coverage-threshold {consensus_coverage_threshold} --label <run dir basename> [--segments <keys>] [--scheme <bed>] [--af-threshold N] [--minimum-length N] [--minimum-map-quality N] [--clair3-model <name>]` | Builds the self-contained HTML page from `summary.tsv` and the artefacts around it. `--scheme` is passed only when `scheme` is not `NA` (primer extent and amplicon depths); the last four flags are display-only footer values, `--clair3-model` only when the config holds a single model name. The rule declares as inputs everything the script reads, so it runs after them: `summary.tsv`, `versions.tsv`, `config.yml`, the basewise tables, the per-sample and filtered `consensus/` FASTAs, the per-sample VCFs (`generate_vcf_consensus` on Illumina, `infer_consensus_sequence` on Nanopore), viralQC results and status (`run_viralqc`), fastp JSONs (Illumina), `isnvs_summary.tsv` (`run_isnv`), status and model files (Nanopore). `run_manifest.json` is read when present but not declared (it is written by the CLI, not the workflow). The script re-derives the headline figures from `summary.tsv` and fails if the page would disagree. | 1 | 17 / 16 |

`viralconseq create-report <run>` runs the same script outside Snakemake.

## Provenance and terminus

| Rule | Env | Command(s) | Flags and why | Threads / mem_mb | Measured max_rss |
|---|---|---|---|---|---|
| `versions_alignment` | `alignment` | `python python/tool_versions.py --tools minimap2 samtools bedtools GSAlign gofasta --output other/versions/alignment.tsv` | Probes each tool's version command inside the env that will run it. Always in the DAG. | 1 | 1 / 3 |
| `versions_qc` | `qc` | `... --tools fastp multiqc` | Illumina only. | 1 | 73 |
| `versions_consensus` | `consensus` | `... --tools lofreq bcftools` | Illumina and `run_isnv` only (no env is built just to read a version). | 1 | 1 |
| `versions_clair3` | `clair3` | `... --tools clair3 samtools bcftools` | Nanopore only. | 1 | 1 |
| `versions_utils` | `utils` | `... --tools python samtools bcftools seqtk` | Always. | 1 | 1 / 19 |
| `versions_viralqc` | `viralqc` | `... --tools nextclade blastn --python-dists viralQC` | `run_viralqc` only; viralQC's version is read from Python package metadata. | 1 | 1 / 3 |
| `versions` | `utils` | `printf` the `component`/`version` header, `viralconseq` (`viralconseq_version` from the config) and `snakemake`; `tail -n +2` every fragment through `awk '!seen[$0]++'`; when viralQC ran, `viralqc_nextclade_datasets_built` (second line of `.nextclade_datasets_ok`, `unknown` if empty) and `viralqc_db` | The `awk` lists a tool pinned in two envs (samtools, bcftools) once when the versions agree and twice when they differ. | 1 | 1 / 1 |
| `run_config` | none (`run:` block) | `shutil.copyfile(<last --configfile>, config.yml)`; when the config file *is* `<run>/config.yml` (Snakemake would already have deleted it), `yaml.safe_dump` of the in-memory config instead | Copying preserves the section comments `ConfigGenerator` writes. Runs in the Snakemake process, which is why its `max_rss` is the interpreter's footprint, not a tool's. | 1 | 97 / 97 |
| `organize_files` (entry file) | `utils` | `mkdir -p samples/<sample>[/<segment>]` and one `ln -sf` per artefact; `touch samples/.organized` | Per-platform link sets (see [Output Layout](output.md)): Illumina links `fastp.html`, `consensus.vcf.gz(.tbi)`, `isnvs.vcf.gz(.tbi)` (`run_isnv`), `stats.tsv`, `consensus.fasta`, `raw_mapped_reads.bam(.bai)`, `trimmed_mapped_reads.bam(.bai)`, `viralqc.tsv` (`run_viralqc`); Nanopore links `consensus.vcf.gz`, `raw.vcf.gz`, `clair3_model.txt`, `status.txt`, `stats.tsv`, `table_cov_basewise.txt`, `consensus.fasta`, both BAMs and `viralqc.tsv`. The sentinel exists because a symlink tree has no single file to declare. | 1 | 1 / 5 |
| `collect_benchmarks` | `utils` | `python python/collect_benchmarks.py --logs-dir logs/ --samples <ids> [--segmented] --cpus <rule>=<threads> ... --exclude <its own benchmark> --output benchmark.tsv` | Walks `logs/` for `*.benchmark.txt`, derives `rule`, `segment` and `target` from the path, labels run-level rows `All`, fills the `threads` column from `--cpus` (one entry per rule in `THREADED_RULES`), and orders per-sample rows first. Depends on `TERMINAL_INPUTS`, the `organize_files` sentinel, `versions.tsv` and `config.yml`, so it is the last job of every run. | 1 | — (excludes itself) |

## Resources

Two helpers in `rules/common.smk` turn config keys into Snakemake directives:

- `cpus(rule)` returns `int(config.get("<rule>_cpus", config["threads"]))`. It
  is called in the `threads:` directive of exactly the rules that have a
  `--<rule>-cpus` option (`ResourceDefaults.CONSENSUS_ILLUMINA_RULES` /
  `CONSENSUS_NANOPORE_RULES` in `constants.py`, mirrored by `THREADED_RULES`
  in `common.smk` and checked by the rule-inventory test):

  | Workflow | Rules calling `cpus()` |
  |---|---|
  | Illumina | `perform_qc`, `map_reads`, `trim_primer_sequences`, `detect_isnv`, `run_viralqc` |
  | Nanopore | `map_reads`, `trim_primer_sequences`, `infer_consensus_sequence`, `run_viralqc` |

  Every other rule has no `threads:` directive and receives one core.
- `ram_mb(rule)` returns `int(config["<rule>_ram"]) * 1024` and is used in the
  `resources: mem_mb=` directive of the two rules in
  `ResourceDefaults.MEMORY_RULES`: `infer_consensus_sequence` on Nanopore
  (default **2 GB**; Clair3 peaked at 0.53 GB on the fixture) and
  `run_viralqc` (default **1 GB**; Nextclade + BLAST peaked at 0.38 GB). Their
  `<rule>_ram` keys are **required** (`common.smk` adds them to
  `_REQUIRED_KEYS`; `ConfigGenerator` always writes them). No other rule
  declares memory: the figures are only ever taken from a measured run.

How the CLI options combine (`consensus_cli.py`, `validators.resolve_resource_budget`,
`_orchestrator.run_workflow`):

| Option | Config key | Effect |
|---|---|---|
| `--threads` (default 1) | `threads` | Baseline thread count of every rule that calls `cpus()`. |
| `--<rule>-cpus` | `<rule>_cpus` | Written only when given; overrides `--threads` for that rule. Snakemake caps a rule's threads at `--threads-total`. |
| `--threads-total` (default: cores available to the process, minus one) | `threads_total` (record only) | `snakemake(cores=...)`: jobs run concurrently while the sum of their threads fits. Single-threaded rules count one each. |
| `--<rule>-ram` (GB; `infer_consensus_sequence`, `run_viralqc` only) | `<rule>_ram` | The `mem_mb` the rule declares; counts against the budget. |
| `--max-memory` (GB) | `max_memory_mb`, `memory_detected_mb` (records) | Passed as `snakemake(resources={"mem_mb": budget})` when positive, so the two declarations bind: those jobs run concurrently only while their declared sum fits. Default: detected memory (`MemTotal` capped by the cgroup limit) less 10 % headroom, raised to the largest per-rule figure when the machine is smaller than one such job. An explicit value that is positive but below the largest `<rule>_ram` is refused. `0` disables the budget (declarations are then informational). When memory cannot be detected the run proceeds without a budget and says so. |

The rules that do not declare memory are bounded by `--threads-total` alone.
`threads_total`, `max_memory_mb` and `memory_detected_mb` are written to the
`# --- resources ---` section of the config as a record of what the run was
given; the workflows do not read them, so a bare `snakemake -s ... --configfile`
rerun must pass `--cores` and `--resources mem_mb=` itself.

## Logs and benchmarks

Every rule declares `log:` and `benchmark:` through the `LOG()` / `BENCH()`
helpers, which build

```
logs/<rule>/[<segment>/]<target>.log
logs/<rule>/[<segment>/]<target>.benchmark.txt
```

`<target>` is `sample-<id>` for per-sample rules and the rule's own name for
run-level rules (`benchmark.tsv` reports those as sample `All`). The
`<segment>/` level appears only in segmented workflows and only for rules that
run once per segment: per-sample rules that are not per segment (`perform_qc`,
`split_viralqc_results`) and all run-level rules omit it, while the two
run-level rules that do run per segment (`align_consensus_to_reference_genome`,
`collect_consensus`) keep it (`logs/collect_consensus/<segment>/collect_consensus.log`).
`sanitize_reference` is the one exception in naming: its target is `reference`
(single) or `<segment>` (segmented, without a segment directory).

Every `shell:` body starts with `set -euo pipefail` and redirects into its log
(`exec > {log} 2>&1`, or `exec 2> {log}` where stdout is data, as in
`map_reads`, `check_mapped_reads`, `prepare_viralqc_input`, `split_viralqc_results`,
`versions` and `collect_benchmarks`); the `script:` rules redirect from Python.
Nothing is only on the console except the one-line `WARNING`s of
`generate_vcf_consensus`, the Nanopore `infer_consensus_sequence` gate and
`run_viralqc`. `collect_benchmarks.py` relies on this layout to attribute each
`*.benchmark.txt` to a rule, segment and sample, and the rule-inventory test
rejects any other form. `logs/run.log` (start and end lines from the hooks) and
`logs/snakemake.log` (a copy of the Snakemake transcript, made by the CLI) sit
next to the rule directories. The columns of `benchmark.tsv` are described in
[Output Layout](output.md#benchmark-columns).
