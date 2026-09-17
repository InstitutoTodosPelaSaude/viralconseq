# Run report: one self-contained HTML page (<run>/report.html) built from
# summary.tsv and the artefacts around it by scripts/python/build_report.py
# with the shell in scripts/templates/report.html. No external resources.
#
# Included by the four entry-point workflows after rules/collect.smk (it uses
# its helpers and rules.summary / rules.collect_consensus) and before
# rules/provenance.smk, so versions.tsv and config.yml are referenced by their
# literal paths. Gated by config["run_report"] (default true) in the entry
# files' TERMINAL_INPUTS; run_manifest.json is read when present but not
# declared (absent under a bare snakemake invocation).

import os

REPORT_DATA = str(config.get("data", "illumina"))
REPORT_RUN_VIRALQC = bool(config.get("run_viralqc", True))
REPORT_RUN_ISNV = REPORT_DATA == "illumina" and bool(config.get("run_isnv", False))
REPORT_SCHEME = str(config.get("scheme", "NA"))
REPORT_TEMPLATE = os.path.join(workflow.basedir, "templates", "report.html")
REPORT_LABEL = os.path.basename(os.path.normpath(config["output"]))

if REPORT_DATA == "illumina":
    _REPORT_VCF = rules.generate_vcf_consensus.output.vcf
else:
    _REPORT_VCF = rules.infer_consensus_sequence.output.vcf


def _report_display_flags():
    """Parameters shown in the page footer (display only, never analysed)."""
    flags = []
    for key, flag in (
        ("af_threshold", "--af-threshold"),
        ("minimum_length", "--minimum-length"),
        ("minimum_map_quality", "--minimum-map-quality"),
    ):
        if key in config:
            flags += [flag, str(config[key])]
    if isinstance(config.get("clair3_model"), str):
        flags += ["--clair3-model", config["clair3_model"]]
    return flags


rule report:
    conda:
        "../envs/utils.yaml"
    input:
        script = os.path.join(PY, "build_report.py"),
        template = REPORT_TEMPLATE,
        summary = rules.summary.output.summary,
        # Literal paths: rules/provenance.smk is included after this module.
        versions = config["output"] + "versions.tsv",
        config_copy = config["output"] + "config.yml",
        basewise = _per_sample(rules.calculate_coverage_basewise.output.table_cov),
        consensus = _per_segment(rules.collect_consensus.output.per_sample),
        filtered = _per_segment(rules.collect_consensus.output.filtered),
        vcfs = _per_sample(_REPORT_VCF),
        viralqc = rules.run_viralqc.output.results if REPORT_RUN_VIRALQC else [],
        viralqc_status = rules.run_viralqc.output.status if REPORT_RUN_VIRALQC else [],
        fastp = (
            expand(rules.perform_qc.output.json, sample=config["samples"])
            if REPORT_DATA == "illumina"
            else []
        ),
        isnvs = config["output"] + "isnvs/isnvs_summary.tsv" if REPORT_RUN_ISNV else [],
        status_files = (
            _per_sample(rules.check_mapped_reads.output.status) if REPORT_DATA == "nanopore" else []
        ),
        model_files = (
            _per_sample(rules.infer_consensus_sequence.output.model_txt)
            if REPORT_DATA == "nanopore"
            else []
        ),
        scheme = REPORT_SCHEME if REPORT_SCHEME != "NA" else [],
    output:
        report = config["output"] + "report.html",
    params:
        run_dir = config["output"],
        data_type = REPORT_DATA,
        min_depth = config["minimum_depth"],
        threshold = COLLECT_COV_T,
        label = REPORT_LABEL,
        segments = (["--segments"] + SEGMENT_KEYS) if SEGMENT_KEYS else [],
        scheme = ["--scheme", REPORT_SCHEME] if REPORT_SCHEME != "NA" else [],
        display = _report_display_flags(),
    log:
        LOG("report", target="report", per_segment=False)
    benchmark:
        BENCH("report", target="report", per_segment=False)
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        python {input.script:q} --run-dir {params.run_dir:q} --template {input.template:q} \
            --output {output.report:q} --data-type {params.data_type} \
            --min-depth {params.min_depth} --consensus-coverage-threshold {params.threshold} \
            --label {params.label:q} {params.segments:q} {params.scheme:q} {params.display:q}
        """
