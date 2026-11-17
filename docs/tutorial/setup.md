# Setup

This page covers everything you need before running the pipeline: installing viralconseq, downloading the example data and generating sample sheets. It is one-time work; once you finish, move on to the [consensus walkthrough](consensus.md).

## 1. Install viralconseq

viralconseq is a Python package. It needs conda and seqkit at runtime, and the project's `environment.yml` provides both:

```bash
conda env create -n viralconseq -f https://raw.githubusercontent.com/InstitutoTodosPelaSaude/viralconseq/main/environment.yml
conda activate viralconseq
pip install viralconseq
viralconseq --version
```

See [Installation](../installation.md) for installing from source, the macOS note and troubleshooting.

## 2. Run the first-time setup

The bioinformatics tools run in their own conda environments, the consensus QC step needs its databases and the Nanopore workflow needs Clair3 models. One command fetches all of them; it takes a while the first time and is skipped on every later run:

```bash
viralconseq setup --pipelines all
```

Where the downloads go and how to skip parts of them is described in [Installation](../installation.md#first-time-setup).

## 3. Download the example data

The example reads are a single tarball. Download it into your working directory:

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

The tutorial also needs the ARTIC nCoV-2019 V3 reference FASTA and primer BED:

```bash
mkdir -p databases/refs
curl -L -o databases/refs/nCoV-2019.reference.fasta \
  https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/V3/nCoV-2019.reference.fasta
curl -L -o databases/refs/nCoV-2019.bed \
  https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/V3/nCoV-2019.scheme.bed
```

## 4. Generate sample sheets

A *sample sheet* is a CSV without a header that tells viralconseq which FASTQ files belong to which sample: `sample_id,R1,R2` for Illumina, `sample_id,fastq` for Nanopore. You can write it by hand, or let `viralconseq create-samplesheet` build it by scanning a run directory.

The example FASTQs sit directly in the run directory rather than in one subdirectory per sample, so we pass `--level 0` and tell the command how to derive the sample ID from the file name:

```bash
# Illumina: sample ID = part of the file name before the first '_'
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
# Nanopore: sample ID = part of the file name before the first '.'
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

For your own data, the default `--level 1` expects one subdirectory per sample, which is how most sequencers deliver a run. The [commands reference](../commands.md#viralconseq-create-samplesheet) lists the options.

## 5. Check what you have

```text
.
├── databases/refs/
│   ├── nCoV-2019.reference.fasta
│   └── nCoV-2019.bed
├── my_test_data/
├── samples_illumina.csv
└── samples_nanopore.csv
```

You are ready to [run the consensus pipeline](consensus.md).
