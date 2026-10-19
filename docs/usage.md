# Usage Overview

viralconseq is one command with five subcommands:

```
viralconseq [--version] [--log-level LEVEL] [--json-logs]
├── setup                Build the tool environments and download the viralQC databases and Clair3 models (once)
├── create-samplesheet   Generate a sample-sheet CSV from a run directory
├── consensus
│   ├── illumina         Reference-guided consensus assembly for Illumina paired-end data
│   └── nanopore         Reference-guided consensus assembly for Nanopore data
├── rerun                Run a workflow again from the config file a previous run wrote
└── create-report        Rebuild report.html for a finished run directory
```

Use `--help` at any level for the full option list, for example `viralconseq consensus illumina --help`. The [commands reference](commands.md) documents every option.

Global options go before the subcommand: `--log-level {DEBUG,INFO,WARNING,ERROR}` and `--json-logs` (one JSON object per log line), e.g. `viralconseq --log-level DEBUG consensus illumina ...`.

## General workflow

1. **Once per machine**: build the tool environments and download the viralQC databases and Clair3 models with `viralconseq setup --pipelines all` (see [Installation](installation.md)).
2. **Per run**: generate a sample sheet with `viralconseq create-samplesheet`, or write the CSV by hand.
3. Run `viralconseq consensus illumina` or `viralconseq consensus nanopore` with the sample sheet, the reference and, for amplicon data, the primer scheme.
4. Read `summary.tsv` and `report.html` in the run directory; the consensus sequences are in `consensus/` (see [Output layout](output.md)).

```{tip}
Use absolute paths for the sample sheet, reference, primer scheme and output directory to avoid mistakes.
```
