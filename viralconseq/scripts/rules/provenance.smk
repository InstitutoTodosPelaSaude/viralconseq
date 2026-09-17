# Workflow-side provenance: which tool versions produced this run.
#
# Included LAST by the four entry-point workflows (after rules/viralqc.smk).
# Expects: config (output, data, run_isnv, run_viralqc, viralqc_db, optional
# viralconseq_version), PY / LOG / BENCH from rules/common.smk.
#
# One ``versions_<env>`` rule per conda environment probes the tools installed
# in that environment (``scripts/python/tool_versions.py``, stdlib only) and
# writes a fragment under other/versions/. Only the environments the DAG uses
# are probed, so no environment is built just to read a version. ``versions``
# concatenates the fragments into <run>/versions.tsv (component<TAB>version)
# with the viralconseq and Snakemake versions on top and, when viralQC ran, the
# database directory and the date its Nextclade datasets were downloaded.

import os

from snakemake import __version__ as SNAKEMAKE_VERSION

PROVENANCE_DATA = str(config.get("data", "illumina"))
PROVENANCE_RUN_VIRALQC = bool(config.get("run_viralqc", True))
PROVENANCE_VIRALQC_DB = str(config.get("viralqc_db", "") or "")

VERSION_ENVS = ["alignment", "utils"]
VERSION_ENVS.append("qc" if PROVENANCE_DATA == "illumina" else "clair3")
if PROVENANCE_DATA == "illumina" and config.get("run_isnv", False):
    VERSION_ENVS.append("consensus")
if PROVENANCE_RUN_VIRALQC:
    VERSION_ENVS.append("viralqc")

VERSIONS_DIR = config["output"] + "other/versions/"
TOOL_VERSIONS_PY = os.path.join(PY, "tool_versions.py")


rule versions_alignment:
    conda:
        "../envs/alignment.yaml"
    input:
        script = TOOL_VERSIONS_PY
    output:
        fragment = VERSIONS_DIR + "alignment.tsv"
    log:
        LOG("versions_alignment", target="versions_alignment", per_segment=False)
    benchmark:
        BENCH("versions_alignment", target="versions_alignment", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        python {input.script:q} --tools minimap2 samtools bedtools GSAlign gofasta \
            --output {output.fragment:q}
        cat {output.fragment:q}
        """


rule versions_qc:
    conda:
        "../envs/qc.yaml"
    input:
        script = TOOL_VERSIONS_PY
    output:
        fragment = VERSIONS_DIR + "qc.tsv"
    log:
        LOG("versions_qc", target="versions_qc", per_segment=False)
    benchmark:
        BENCH("versions_qc", target="versions_qc", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        python {input.script:q} --tools fastp multiqc --output {output.fragment:q}
        cat {output.fragment:q}
        """


rule versions_consensus:
    conda:
        "../envs/consensus.yaml"
    input:
        script = TOOL_VERSIONS_PY
    output:
        fragment = VERSIONS_DIR + "consensus.tsv"
    log:
        LOG("versions_consensus", target="versions_consensus", per_segment=False)
    benchmark:
        BENCH("versions_consensus", target="versions_consensus", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        python {input.script:q} --tools lofreq bcftools --output {output.fragment:q}
        cat {output.fragment:q}
        """


rule versions_clair3:
    conda:
        "../envs/clair3.yaml"
    input:
        script = TOOL_VERSIONS_PY
    output:
        fragment = VERSIONS_DIR + "clair3.tsv"
    log:
        LOG("versions_clair3", target="versions_clair3", per_segment=False)
    benchmark:
        BENCH("versions_clair3", target="versions_clair3", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        python {input.script:q} --tools clair3 samtools bcftools --output {output.fragment:q}
        cat {output.fragment:q}
        """


rule versions_utils:
    conda:
        "../envs/utils.yaml"
    input:
        script = TOOL_VERSIONS_PY
    output:
        fragment = VERSIONS_DIR + "utils.tsv"
    log:
        LOG("versions_utils", target="versions_utils", per_segment=False)
    benchmark:
        BENCH("versions_utils", target="versions_utils", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        python {input.script:q} --tools python --python-dists pandas \
            --output {output.fragment:q}
        cat {output.fragment:q}
        """


rule versions_viralqc:
    conda:
        "../envs/viralqc.yaml"
    input:
        script = TOOL_VERSIONS_PY
    output:
        fragment = VERSIONS_DIR + "viralqc.tsv"
    log:
        LOG("versions_viralqc", target="versions_viralqc", per_segment=False)
    benchmark:
        BENCH("versions_viralqc", target="versions_viralqc", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        python {input.script:q} --tools nextclade blastn --python-dists viralQC \
            --output {output.fragment:q}
        cat {output.fragment:q}
        """


rule versions:
    conda:
        "../envs/utils.yaml"
    input:
        fragments = expand(VERSIONS_DIR + "{env}.tsv", env=VERSION_ENVS),
        nextclade_marker = (
            [PROVENANCE_VIRALQC_DB + "/.nextclade_datasets_ok"]
            if PROVENANCE_RUN_VIRALQC and PROVENANCE_VIRALQC_DB not in ("", "NA")
            else []
        ),
    output:
        versions = config["output"] + "versions.tsv"
    params:
        viralconseq = str(config.get("viralconseq_version", "unknown")),
        snakemake = SNAKEMAKE_VERSION,
        viralqc_db = (
            os.path.abspath(PROVENANCE_VIRALQC_DB)
            if PROVENANCE_RUN_VIRALQC and PROVENANCE_VIRALQC_DB not in ("", "NA")
            else ""
        ),
    log:
        LOG("versions", target="versions", per_segment=False)
    benchmark:
        BENCH("versions", target="versions", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec 2> {log}
        {{
            printf 'component\\tversion\\n'
            printf 'viralconseq\\t%s\\n' {params.viralconseq:q}
            printf 'snakemake\\t%s\\n' {params.snakemake:q}
            # A tool pinned in two environments (samtools, bcftools) is listed
            # once when the versions agree, twice when they differ.
            for fragment in {input.fragments}; do
                tail -n +2 "$fragment"
            done | awk '!seen[$0]++'
            if [ -n {params.viralqc_db:q} ]; then
                # The marker's second line is the UTC date the datasets were
                # fetched (older setups wrote an empty marker: report unknown).
                built=$(sed -n 2p {input.nextclade_marker:q} 2>/dev/null || true)
                printf 'viralqc_nextclade_datasets_built\\t%s\\n' "${{built:-unknown}}"
                printf 'viralqc_db\\t%s\\n' {params.viralqc_db:q}
            fi
        }} > {output.versions}
        """

