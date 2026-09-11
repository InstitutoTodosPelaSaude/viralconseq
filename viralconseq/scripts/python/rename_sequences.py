#!/usr/bin/env python
"""Rename a per-sample consensus FASTA's header to the sample name.

This script is executed by Snakemake via a ``script:`` directive (see
``rules/consensus_*.smk::rename_sequences``). Snakemake injects a
``snakemake`` global at runtime (``snakemake.input``, ``snakemake.output``)
which ruff/mypy cannot see; F821/ignore_errors are configured for this
directory in ``pyproject.toml``.

The script reads a single-sample consensus FASTA (where the original header
may be anything), strips it, and writes a single record whose header is
``>{sample_name}`` derived from the input filename (``.consensus.fasta``
suffix stripped).
"""

import re


def rename_sequences(input_path: str, output: str) -> None:
    """Rewrite the FASTA at ``input_path`` so its headers carry the sample name.

    The sample name is the basename of ``input_path`` minus the
    ``.consensus.fasta`` suffix. A single-record FASTA is written with the
    header ``>{sample_name}``. A multi-record FASTA (e.g. a multi-contig
    reference) keeps one record *per contig*, headed ``>{sample_name}_{origid}``,
    so distinct contigs are never silently fused into one chimeric sequence.
    The sequence body of each record is preserved verbatim.
    """
    sample_name = re.sub(".+/", "", input_path)
    sample_name = re.sub(".consensus.fasta", "", sample_name)
    with open(input_path) as f:
        lines = f.readlines()

    # Split the FASTA into (original_header, body_lines) records.
    records = []
    current_header = None
    current_body: list = []
    for line in lines:
        if line.startswith(">"):
            if current_header is not None or current_body:
                records.append((current_header, current_body))
            current_header = line[1:].strip()
            current_body = []
        else:
            current_body.append(line)
    if current_header is not None or current_body:
        records.append((current_header, current_body))

    multi_record = len(records) > 1
    chunks = []
    for index, (header, body) in enumerate(records):
        if multi_record:
            original_id = header.split()[0] if header else f"contig{index + 1}"
            name = f"{sample_name}_{original_id}"
        else:
            name = sample_name
        chunks.append(">" + name + "\n" + "".join(body))

    with open(output, "w") as f:
        f.write("".join(chunks))


if __name__ == "__main__":
    print("Running rename_sequences.py")
    rename_sequences(snakemake.input[0], snakemake.output[0])
    print("Finished rename_sequences.py")
    exit()
