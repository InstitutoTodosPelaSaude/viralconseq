# Shared preamble for the four consensus workflows. Included FIRST by every
# entry-point .smk, before ``rule all`` and before any other include, because
# global wildcard constraints only apply to rules parsed after them.
#
# Provides:
# - exact-sample wildcard constraints (a sample wildcard matches only the ids
#   in config["samples"], so ids that are prefixes of one another can never
#   produce an ambiguous match)
# - the run log ``<output>/logs/run.log`` written by onstart/onsuccess/onerror
#
# Expects ``config["samples"]`` and ``config["output"]`` (with its trailing
# slash, as ConfigGenerator writes it).

import datetime
import os
import re
import shlex

from snakemake.exceptions import WorkflowError

# Analysis parameters are REQUIRED keys: the rules read config["key"] with no
# fallback, so a default can only live in one place (the CLI). A hand-edited
# YAML missing one fails here with the full list instead of a KeyError deep in
# a rule. Optional/sentinel keys (scheme, adapters, run_isnv, run_viralqc, the
# two *_flags strings, per-rule *_cpus/*_ram, viralconseq_version) stay
# optional on purpose.
_REQUIRED_KEYS = [
    "samples",
    "data",
    "output",
    "threads",
    "reference",
    "scheme",
    "minimum_depth",
    "minimum_length",
    "af_threshold",
]
if config.get("data") == "illumina":
    _REQUIRED_KEYS += [
        "adapters",
        "trim_head",
        "trim_tail",
        "cut_front_mean_quality",
        "cut_tail_mean_quality",
        "cut_right_window_size",
        "cut_right_mean_quality",
        "af_isnv_threshold",
    ]
else:
    _REQUIRED_KEYS += [
        "chunk_size",
        "clair3_model",
        "variant_quality",
        "variant_depth",
        "minimum_map_quality",
    ]
_missing = [k for k in _REQUIRED_KEYS if k not in config]
if _missing:
    raise WorkflowError(
        "config is missing required key(s): "
        + ", ".join(_missing)
        + " - regenerate it with 'viralconseq consensus ... --create-config-only'"
    )

# viralqc_extra_flags and minimap2_consensus_align_flags are interpolated
# unquoted into shell commands. Refuse anything that is not a plain list of
# tokens (same rule as validators.validate_flag_strings, for hand-edited YAMLs).
for _key in ("viralqc_extra_flags", "minimap2_consensus_align_flags"):
    _value = config.get(_key)
    if _value is None:
        continue
    try:
        shlex.split(str(_value))
    except ValueError as _exc:
        raise WorkflowError(f"config key {_key} is not a valid flag string ({_exc}): {_value!r}")
    if re.search(r"[;&|<>`$\\\n]", str(_value)):
        raise WorkflowError(
            f"config key {_key} contains a shell metacharacter: {_value!r}. "
            "Pass plain tool flags only."
        )

wildcard_constraints:
    sample="|".join(re.escape(s) for s in config["samples"]),
    segment=r"[^/]+",
    ref_key=r"[^/]+",


RUN_LOG = config["output"] + "logs/run.log"


def append_run_log(message):
    """Append one timestamped line to the run log. Never raises: this is
    provenance, and a full disk or read-only directory must not turn a
    finished analysis into a failure."""
    try:
        os.makedirs(os.path.dirname(RUN_LOG), exist_ok=True)
        stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        with open(RUN_LOG, "a") as handle:
            handle.write(f"{stamp}\t{message}\n")
    except OSError:
        pass


onstart:
    append_run_log(
        "start\tviralconseq "
        + str(config.get("viralconseq_version", "unknown"))
        + "\tdata=" + str(config.get("data", "unknown"))
    )
    append_run_log(f"samples\t{len(config['samples'])}")
    append_run_log(f"output\t{config['output']}")
    append_run_log(f"snakemake_log\t{log}")


onsuccess:
    append_run_log("end\tsuccess")
    print("viralconseq: workflow completed successfully.")


onerror:
    append_run_log("end\tfailed")
    print(f"viralconseq: workflow FAILED - see the Snakemake log: {log} and {RUN_LOG}")
