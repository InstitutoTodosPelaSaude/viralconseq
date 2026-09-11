"""Rules shared between the segmented and non-segmented Nanopore consensus
workflows.

The including snakefile must define:
* ``SEGMENT_WILDCARD``       -- ``"{segment}/"`` (segmented) or ``""`` (single)
* ``REFERENCE``              -- path string / callable taking ``wildcards`` and
                                returning a path (segmented) -- forwarded to
                                ``params.reference``
"""


rule calculate_assembly_statistics:
    conda:
        "../envs/utils.yaml"
    input:
        get_map_input_fastqs,
        get_map_input_fastqs,
        get_map_input_fastqs,
        rules.trim_primer_sequences.output.bam,
        rules.calculate_coverage_basewise.output.table_cov,
        rules.rename_sequences.output.consensus_renamed
    output:
        stats_summary = temp(config['output'] + "assembly/" + SEGMENT_WILDCARD + "coverage_stats/{sample}.stats_summary.csv")
    params:
        minimum_depth = config["minimum_depth"]
    script:
        "../python/calculate_assembly_stats.py"


rule align_consensus_to_reference_genome:
    conda:
        "../envs/alignment.yaml"
    input:
        # Literal path (instead of ``rules.unify_assembly_statistics_reports.``)
        # to avoid an include-order cycle: ``unify_assembly_statistics_reports``
        # lives in the top-level snakefile and itself references
        # ``rules.calculate_assembly_statistics`` defined here.
        stats = config['output'] + "assembly/assembly_stats_summary.csv",
        consensus_files = expand(
            rules.rename_sequences.output.consensus_renamed,
            sample=config["samples"],
            allow_missing=True
        )
    output:
        aln_consensus = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/samples_alignment.fasta",
        combined = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/consensus.fasta",
        sam = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/aln.consensus.sam",
        masked = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/aln.consensus.indelsMasked.fasta"
    params:
        reference = REFERENCE,
        minimap2_flags = config.get("minimap2_consensus_align_flags", "-a --sam-hit-only --secondary=no --score-N=0")
    log:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/align_consensus_to_reference_genome.log"
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        cat {params.reference} {input.consensus_files} > {output.combined}
        minimap2 {params.minimap2_flags} {params.reference} {output.combined} -o {output.sam}
        # gofasta sam toMultiAlign builds an MSA against a single reference and
        # aborts on a multi-contig reference (e.g. a fragmented genome supplied
        # with --single-reference). Align one reference contig at a time and
        # concatenate; a single-contig reference (the common case, and every
        # segmented per-segment run) takes the same one-pass path as before.
        grep '^>' {params.reference} | sed 's/^>//' | cut -d' ' -f1 > {output.sam}.contigs
        if [ "$(wc -l < {output.sam}.contigs)" -le 1 ]; then
            gofasta sam toMultiAlign --pad -s {output.sam} -o {output.aln_consensus}
        else
            : > {output.aln_consensus}
            per_contig_dir=$(dirname {output.aln_consensus})/per_contig_alignments
            mkdir -p "$per_contig_dir"
            while read -r contig; do
                awk '$1 == "@HD"' {output.sam} > {output.sam}.one
                awk -F'\t' -v c="$contig" '$1 == "@SQ" && $2 == "SN:"c' {output.sam} >> {output.sam}.one
                awk '$1 == "@PG"' {output.sam} >> {output.sam}.one
                awk -F'\t' -v c="$contig" '$0 !~ /^@/ && $3 == c' {output.sam} >> {output.sam}.one
                # Skip contigs with no mapped records: minimap2 --sam-hit-only
                # still emits their @SQ line, and gofasta aborts on a header-only
                # SAM (e.g. an all-N consensus for an uncovered contig).
                if grep -qv '^@' {output.sam}.one; then
                    # Write each contig's alignment as its own MSA under
                    # per_contig_alignments/ (a multi-contig reference is really
                    # one alignment per contig), then concatenate into the declared
                    # {output.aln_consensus} so the Snakemake DAG output is unchanged.
                    safe=$(printf '%s' "$contig" | sed 's/[^A-Za-z0-9._-]/_/g')
                    gofasta sam toMultiAlign --pad -s {output.sam}.one -o "$per_contig_dir/$safe.fasta"
                    cat "$per_contig_dir/$safe.fasta" >> {output.aln_consensus}
                fi
            done < {output.sam}.contigs
            rm -f {output.sam}.one
        fi
        rm -f {output.sam}.contigs
        sed '/^>/ ! s/-/N/g' {output.aln_consensus} > {output.masked}
        """
