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
    "consensus_coverage_threshold",
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
        "clair3_model_dir",
        "variant_quality",
        "minimum_mapped_reads",
        "variant_depth",
        "minimum_map_quality",
    ]
# Memory-declaring rules read <rule>_ram as a required key (ConfigGenerator
# always writes it for them; see ResourceDefaults.MEMORY_RULES).
if config.get("data") != "illumina":
    _REQUIRED_KEYS.append("infer_consensus_sequence_ram")
if config.get("run_viralqc", True):
    _REQUIRED_KEYS.append("run_viralqc_ram")
_missing = [k for k in _REQUIRED_KEYS if k not in config]
if _missing:
    raise WorkflowError(
        "config is missing required key(s): "
        + ", ".join(_missing)
        + " - regenerate it with 'viralconseq consensus ... --create-config-only'"
    )

# The full contract check (types, bounds, flag strings, memory budget) lives in
# viralconseq.validators.validate_config_dict; the workflows run in the CLI
# environment where the package is installed, so it is applied here too, for
# hand-edited YAMLs driven with `snakemake -s`. A bare checkout without the
# package keeps the key-presence guard above and skips the rest.
try:
    from viralconseq.validators import validate_config_dict as _validate_config_dict
except ImportError:  # pragma: no cover - only without `pip install -e .`
    import sys as _sys

    print(
        "viralconseq: package not importable; skipping the full config check",
        file=_sys.stderr,
    )
else:
    try:
        _validate_config_dict(config)
    except Exception as _exc:
        raise WorkflowError(str(_exc)) from _exc

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
LOGS = config["output"] + "logs/"
# Helper scripts (stdlib, argparse) run via shell: live next to the workflows.
PY = os.path.join(workflow.basedir, "python")


def LOG(rule, target="{sample}", per_segment=True):
    """Log path ``<output>/logs/<rule>/[<segment>/]<target>.log``.

    ``per_segment`` inserts the segment wildcard in segmented workflows so a
    per-sample rule that runs once per segment does not overwrite its own log.
    Run-level rules pass ``target=<rule name>, per_segment=False``.
    """
    segment = globals().get("SEGMENT_WILDCARD", "") if per_segment else ""
    return LOGS + rule + "/" + segment + target + ".log"


# The rules with a --<rule>-cpus option (mirror of ResourceDefaults.CONSENSUS_*_RULES;
# the inventory test keeps them in step). Used to report each rule's threads.
THREADED_RULES = (
    ["perform_qc", "map_reads", "trim_primer_sequences", "detect_isnv", "run_viralqc"]
    if config.get("data") == "illumina"
    else ["map_reads", "trim_primer_sequences", "infer_consensus_sequence", "run_viralqc"]
)


def clair3_model_for(sample):
    """The Clair3 model of one sample: config["clair3_model"] is a name or a
    {sample: name} mapping (samples basecalled with different models)."""
    value = config.get("clair3_model")
    return value[sample] if isinstance(value, dict) else value


def clair3_model_files(wildcards):
    """The two checkpoint files of a sample's model, declared as rule inputs so
    a missing model is a MissingInputException at DAG time, not a Clair3 crash."""
    directory = os.path.join(str(config["clair3_model_dir"]), str(clair3_model_for(wildcards.sample)))
    return [os.path.join(directory, "pileup.pt"), os.path.join(directory, "full_alignment.pt")]


def cpus(rule):
    """Threads for ``rule``: ``--<rule>-cpus`` if given, else the ``--threads`` baseline."""
    return int(config.get(f"{rule}_cpus", config["threads"]))


def ram_mb(rule):
    """``mem_mb`` for a memory-declaring rule, from its ``<rule>_ram`` (GB) key.

    The key is required exactly when the rule can run, so a rule whose whole
    step is optional (``run_viralqc``) must call this from a ``resources:``
    lambda, which Snakemake evaluates per job instead of at parse time.
    """
    return int(config[f"{rule}_ram"]) * 1024


def BENCH(rule, target="{sample}", per_segment=True):
    """Benchmark path next to :func:`LOG`: ``<...>/<target>.benchmark.txt``."""
    return LOG(rule, target=target, per_segment=per_segment)[: -len(".log")] + ".benchmark.txt"


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
