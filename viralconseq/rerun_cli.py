"""``viralconseq rerun``: run a workflow again from a saved config.

The YAML written by ``viralconseq consensus`` is the contract the workflows
read, so it is also the way to replay a run: after an interruption, to plan a
change (``--dry-run``), to release a stale lock (``--unlock``), or with a few
keys changed (``--set KEY=VALUE``). Overrides are written back into the file
(the previous copy is kept as ``<config>.bak``) so the config on disk always
matches the run it describes; nothing is passed to Snakemake in memory only.
"""

import logging
import os
from typing import Any, Dict, Sequence

import click
import yaml

from viralconseq import _orchestrator
from viralconseq.config_generator import ConfigGenerator
from viralconseq.consensus import workflow_path_for
from viralconseq.consensus_cli import _default_conda_prefix
from viralconseq.constants import ResourceDefaults
from viralconseq.exceptions import ConfigurationError, ViralConseqError
from viralconseq.provenance import MANIFEST_FILENAME, record_run_completion
from viralconseq.validators import validate_config_dict

logger = logging.getLogger(__name__)


def load_config(path: str) -> Dict[str, Any]:
    """Read a config YAML; raise ``ConfigurationError`` unless it is a mapping."""
    try:
        with open(path) as handle:
            loaded = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as e:
        raise ConfigurationError(f"Cannot read config {path}: {e}") from e
    if not isinstance(loaded, dict):
        raise ConfigurationError(f"{path} does not hold a mapping at the top level")
    return loaded


def apply_overrides(config: Dict[str, Any], overrides: Sequence[str]) -> Dict[str, Any]:
    """Apply ``KEY=VALUE`` pairs; values are parsed as YAML (``4`` -> int,
    ``true`` -> bool, ``'x'`` -> str). Only keys already in the config may be
    set, which catches typos before they become silently ignored keys."""
    for entry in overrides:
        key, sep, raw = entry.partition("=")
        if not sep or not key.strip():
            raise click.BadParameter(f"expected KEY=VALUE, got {entry!r}", param_hint="--set")
        key = key.strip()
        if key not in config:
            raise click.BadParameter(
                f"{key!r} is not a key of this config (keys: {', '.join(sorted(config))})",
                param_hint="--set",
            )
        try:
            config[key] = yaml.safe_load(raw) if raw.strip() else ""
        except yaml.YAMLError as e:
            raise click.BadParameter(f"{key}: {e}", param_hint="--set") from e
    return config


def args_from_config(config: Dict[str, Any], config_file: str, conda_prefix: str) -> Dict[str, Any]:
    """The subset of the pipeline ``args`` dict ``run_workflow`` needs."""
    output = str(config["output"]).rstrip("/")
    args = {
        "config_file": os.path.abspath(config_file),
        "data_type": config["data"],
        "threads_total": int(config.get("threads_total") or 1),
        "max_memory_mb": int(config.get("max_memory_mb") or 0),
        "memory_detected_mb": int(config.get("memory_detected_mb") or 0),
        "memory_budget_source": "saved config",
        "conda_prefix": conda_prefix,
        # output/run_name only feed run_dir_for: <output>/<run_name> must be
        # the run directory, which is exactly config["output"].
        "output": os.path.dirname(output),
        "run_name": os.path.basename(output),
    }
    # describe_resources reads each memory-declaring rule's allowance from
    # args; without these the line printed at start would quote the built-in
    # default while the rule itself uses the saved <rule>_ram value.
    for rule in ResourceDefaults.MEMORY_RULES:
        if f"{rule}_ram" in config:
            args[f"{rule}_ram"] = config[f"{rule}_ram"]
    return args


@click.command(name="rerun")
@click.argument("config_file", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--dry-run", is_flag=True, default=False, help="Plan the run and stop (snakemake -n)."
)
@click.option(
    "--unlock",
    is_flag=True,
    default=False,
    help="Release the lock left by an interrupted run, then exit.",
)
@click.option(
    "--keep-going",
    is_flag=True,
    default=False,
    help="Keep running independent jobs after one fails.",
)
@click.option(
    "--set",
    "overrides",
    multiple=True,
    metavar="KEY=VALUE",
    help="Change a config key before running (repeatable). The value is parsed as "
    "YAML; the file is rewritten and the previous copy kept as CONFIG_FILE.bak.",
)
@click.option(
    "--conda-prefix",
    default=_default_conda_prefix,
    show_default="$VIRALCONSEQ_CONDA_PREFIX or ~/.cache/viralconseq/conda-envs",
    help="Directory where per-rule conda envs are cached.",
)
def rerun(
    config_file: str,
    dry_run: bool,
    unlock: bool,
    keep_going: bool,
    overrides: Sequence[str],
    conda_prefix: str,
) -> None:
    """Run a workflow again from the config CONFIG_FILE written by a previous run.

    Use it to resume an interrupted run, to preview what would run (--dry-run),
    to release a stale lock (--unlock), or to change a setting (--set) without
    retyping the whole command. The config is validated before Snakemake starts.
    """
    try:
        config = load_config(config_file)
        apply_overrides(config, overrides)
        validate_config_dict(config)
        if overrides:
            ConfigGenerator.from_dict(config_file, config).save(backup=True)
            click.echo(f"config updated; previous copy kept as {config_file}.bak", err=True)
    except ViralConseqError as e:
        raise click.ClickException(f"[{e.code}] {e}") from e

    args = args_from_config(config, config_file, conda_prefix)
    workflow = workflow_path_for(config["data"], config.get("reference"))
    rule_list = _orchestrator.rule_list_for_args(args)
    logger.info(_orchestrator.describe_resources(args, rule_list))

    try:
        successful = _orchestrator.run_workflow(
            workflow, args, dryrun=dry_run, unlock=unlock, keepgoing=keep_going
        )
    except ViralConseqError as e:
        raise click.ClickException(f"[{e.code}] {e}") from e

    if not (dry_run or unlock):
        manifest = os.path.join(str(config["output"]), MANIFEST_FILENAME)
        if os.path.isfile(manifest):
            try:
                record_run_completion(manifest, status="success" if successful else "failed")
            except Exception as e:  # provenance is best-effort
                logger.warning(f"Could not update {manifest}: {e}")
    raise SystemExit(0 if successful else 1)
