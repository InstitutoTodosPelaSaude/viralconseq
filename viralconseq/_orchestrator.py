"""Shared pipeline-orchestration helpers used by ``viralconseq.consensus``.

The pipeline module still owns its own ``validate_args``, ``generate_config_file``,
``run_snakemake_workflow``, and ``main`` so that existing test patches at
those module-level names continue to work; the helpers in this module are
called from inside those functions and keep the orchestration flow
(resolve paths -> validate -> write config -> run Snakemake) in one place.
"""

import logging
import os
from typing import Any, Callable, Dict, Optional

from snakemake import snakemake

from viralconseq.config_generator import ConfigGenerator
from viralconseq.exceptions import ValidationError, ViralConseqError
from viralconseq.provenance import record_run_completion, write_run_manifest

logger = logging.getLogger(__name__)


def start_config(args: Dict[str, Any], samples: Dict[str, list]) -> ConfigGenerator:
    """Instantiate a ``ConfigGenerator`` and add the keys every workflow
    always needs: samples, output, threads.

    The caller is responsible for adding pipeline-specific keys, per-rule
    resource settings, and calling ``generator.save()``.
    """
    logger.info("Generating configuration file")
    generator = ConfigGenerator(args["config_file"])
    generator.add_samples(samples, args["data_type"])
    generator.add_output(args["output"], args["run_name"])
    generator.add_threads(args["threads"])
    return generator


def run_workflow(workflow_path: str, args: Dict[str, Any]) -> bool:
    """Run a Snakemake workflow with the kwargs both pipelines share.

    Args:
        workflow_path: Absolute path to the ``.smk`` file to execute.
        args: Pipeline argument dict. Must contain ``config_file`` and
            ``threads_total``.

    Raises:
        ValidationError: If ``workflow_path`` does not exist.
    """
    logger.info("Starting Snakemake workflow")

    if not os.path.isfile(workflow_path):
        raise ValidationError(f"Workflow file not found: {workflow_path}")

    successful = snakemake(
        workflow_path,
        configfiles=[args["config_file"]],
        cores=args["threads_total"],
        use_conda=True,
        conda_prefix=args.get("conda_prefix"),
        targets=["all"],
    )

    if successful:
        logger.info("Snakemake workflow completed successfully")
    else:
        logger.error("Snakemake workflow failed")

    return successful


def run_pipeline(
    args: Dict[str, Any],
    *,
    resolve_paths: Callable[[Dict[str, Any]], object],
    validate: Callable[[Dict[str, Any]], Optional[Dict[str, list]]],
    generate_config: Callable[[Dict[str, list], Dict[str, Any]], None],
    run_workflow_fn: Callable[[Dict[str, Any]], bool],
    skip_when_no_samples: bool = False,
) -> int:
    """The try/except ``main`` skeleton shared by both pipelines.

    Args:
        args: Pipeline argument dict.
        resolve_paths: Called first to convert relative path args to
            absolute paths.
        validate: Called to validate args and return the samples dict.
        generate_config: Called to write the Snakemake config file.
        run_workflow_fn: Called to execute the Snakemake workflow. The
            pipeline-specific module passes its own
            ``run_snakemake_workflow`` so test patches at that name
            keep working.
        skip_when_no_samples: If ``True``, return 0 (success) when the
            validated samples dict is empty. The consensus pipeline does not
            opt in (an empty sample set is a validation error there).

    Returns:
        Exit code (0 for success, 1 for failure).
    """
    try:
        resolve_paths(args)
        samples = validate(args)

        if skip_when_no_samples and (samples is None or len(samples) == 0):
            logger.warning("No samples were provided.")
            return 0

        # ``validate`` is typed as returning Optional; narrow for the type
        # checker (the consensus pipeline always returns a dict).
        if samples is None:
            samples = {}

        generate_config(samples, args)

        if args.get("create_config_only", False):
            logger.info("Config file created. Exiting without running workflow.")
            return 0

        # Provenance: record version, config hash, and input checksums for the
        # run. Best-effort — a manifest failure must never abort an analysis.
        manifest_path: Optional[str] = None
        try:
            manifest_path = write_run_manifest(args, samples)
            logger.info(f"Wrote run manifest: {manifest_path}")
        except Exception as e:
            logger.warning(f"Could not write run manifest: {e}")

        successful = run_workflow_fn(args)

        # Patch the manifest with the actual outcome (it was written pre-run).
        if manifest_path is not None:
            try:
                record_run_completion(manifest_path, status="success" if successful else "failed")
            except Exception as e:
                logger.warning(f"Could not update run manifest with completion status: {e}")

        return 0 if successful else 1

    except ViralConseqError as e:
        # Expected, user-facing failures (bad inputs, missing files/DBs, config
        # errors). Log a clean message with the machine-readable error code so a
        # caller/service can key off it; no stack trace for these.
        logger.error(f"[{e.code}] {e}")
        return 1
    except Exception as e:
        # Genuinely unexpected: keep the traceback for debugging.
        logger.exception(f"Unexpected error: {e}")
        return 1
