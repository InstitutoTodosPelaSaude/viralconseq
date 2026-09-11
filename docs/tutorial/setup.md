# Setup

This page covers everything you need before running the pipeline: installing viralconseq, building the per-rule environments, downloading the example data and generating sample sheets. Most of it is one-time work; once you finish, you can move straight to the [consensus walkthrough](consensus.md).

## 1. Install viralconseq

viralconseq is a Python package; the per-rule bioinformatics tools (minimap2, samtools, fastp, Clair3, …) live in conda environments under `viralconseq/scripts/envs/` that Snakemake builds and picks up automatically.

```bash
pip install viralconseq
```

Or, from a source checkout (needed for development):

```bash
git clone https://github.com/filiperomero2/viralconseq.git
cd viralconseq
conda env create -n viralconseq -f environment.yml
conda activate viralconseq
pip install -e .
```

Either way, conda (or mamba/micromamba) must be available on your `PATH` at runtime.

Verify:

```bash
viralconseq --version
viralconseq --help
```

See [Installation](../installation.md) for the macOS Apple Silicon caveat (the `clair3` environment can be sensitive there), for development installs (`pip install -e ".[dev]"`), and for the [Troubleshooting](../installation.md#troubleshooting) note if the first pipeline run fails on conda env creation.

## 2. Build per-rule environments

Each pipeline rule runs in its own conda environment. viralconseq caches these envs across runs so they only need to be built once; pre-warm the cache now:

```bash
viralconseq setup --pipelines all
```

This materializes every per-rule env into `~/.cache/viralconseq/conda-envs/` (override with `--conda-prefix PATH` or `$VIRALCONSEQ_CONDA_PREFIX`). The subsequent pipeline runs in this tutorial — and every future run on this machine — reuse the cache and skip env creation. Restricting to a single workflow (e.g. `--pipelines consensus-illumina`) is faster if you only plan to run one flavour; add the other later.

If `setup` fails with `CreateCondaEnvironmentException`, see [Troubleshooting](../installation.md#troubleshooting).

## 3. Download the example data

The example dataset used throughout this tutorial is hosted as a single tarball. Download it once into your working directory:

```bash
curl -L -o my_test_data.tar.gz https://itps-nimbus.nyc3.cdn.digitaloceanspaces.com/my_test_data.tar.gz
tar -xzf my_test_data.tar.gz
rm my_test_data.tar.gz
```

You should now have:

```text
my_test_data/
├── illumina_data/
│   ├── itps-0001_R1.fastq.gz
│   ├── itps-0001_R2.fastq.gz
│   ├── itps-0002_R1.fastq.gz
│   └── itps-0002_R2.fastq.gz
└── nanopore_data/
    ├── barcode05.itps-0003.fastq.gz
    └── barcode09.itps-0004.fastq.gz
```

The tutorial additionally needs the ARTIC nCoV-2019 V3 reference FASTA and primer BED:

```bash
mkdir -p databases/refs
curl -L -o databases/refs/nCoV-2019.reference.fasta \
  https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/V3/nCoV-2019.reference.fasta
curl -L -o databases/refs/nCoV-2019.bed \
  https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/V3/nCoV-2019.scheme.bed
```

The rest of this tutorial points `--reference` / `--primer-scheme` at `databases/refs/`.

## 4. Generate sample sheets

A *sample sheet* is a no-header CSV that tells viralconseq which FASTQ files belong to which sample:

- **Illumina** — 3 columns: `sample_id,R1_path,R2_path`
- **Nanopore** — 2 columns: `sample_id,fastq_path`

You can write the CSV by hand, but `viralconseq create-samplesheet` builds it for you by scanning a run directory.

The example data has all FASTQs directly inside the run directory (rather than one subdirectory per sample), so we pass `--level 0` and a small naming convention:

```bash
# Illumina: FASTQs in my_test_data/illumina_data/, sample ID = portion before first '_'
viralconseq create-samplesheet \
    --input  my_test_data/illumina_data/ \
    --output samples_illumina.csv \
    --level 0 \
    --separator _ \
    --pattern R1
```

This writes:

```text
itps-0001,my_test_data/illumina_data/itps-0001_R1.fastq.gz,my_test_data/illumina_data/itps-0001_R2.fastq.gz
itps-0002,my_test_data/illumina_data/itps-0002_R1.fastq.gz,my_test_data/illumina_data/itps-0002_R2.fastq.gz
```

```bash
# Nanopore: FASTQs in my_test_data/nanopore_data/, sample ID = portion before first '.'
viralconseq create-samplesheet \
    --input  my_test_data/nanopore_data/ \
    --output samples_nanopore.csv \
    --level 0 \
    --separator . \
    --pattern barcode
```

Which writes:

```text
barcode05,my_test_data/nanopore_data/barcode05.itps-0003.fastq.gz
barcode09,my_test_data/nanopore_data/barcode09.itps-0004.fastq.gz
```

For your own data, the same command works with `--level 1` (one subdirectory per sample, the default) if your run is organised that way; see the [Commands reference](../commands.md#viralconseq-create-samplesheet) for the full table of options.

## 5. Verify the setup

You should now have:

```text
.
├── databases/refs/
│   ├── nCoV-2019.reference.fasta
│   └── nCoV-2019.bed
├── my_test_data/
├── samples_illumina.csv
└── samples_nanopore.csv
```

A quick sanity check that the workflow environments are in place:

```bash
viralconseq setup --pipelines all --dry-run
```

You are now ready to [run the consensus pipeline](consensus.md).
