# Usage Overview

viralconseq exposes three subcommands:

```
viralconseq [--version] [--log-level LEVEL] [--json-logs]
├── create-samplesheet   Generate a sample-sheet CSV from a run directory
├── setup                Pre-build the per-rule conda environments into a shared cache
└── consensus
    ├── illumina         Reference-guided consensus assembly for Illumina paired-end data
    └── nanopore         Reference-guided consensus assembly for Nanopore data
```

Use `--help` at any level for the full option list:

```bash
viralconseq --help
viralconseq consensus --help
viralconseq consensus illumina --help
viralconseq consensus nanopore --help
```

Global options go before the subcommand: `--log-level {DEBUG,INFO,WARNING,ERROR}` and
`--json-logs` (one JSON object per log line), e.g. `viralconseq --log-level DEBUG consensus illumina ...`.

## General workflow

1. **Build the per-rule environments and download the viralQC databases once** with `viralconseq setup --pipelines all`
2. **Generate a sample sheet** with `viralconseq create-samplesheet` (or write the CSV by hand)
3. **Run the pipeline** with `viralconseq consensus illumina` or `viralconseq consensus nanopore`

```{tip}
Always use absolute paths to avoid mistakes when specifying file locations.
```
