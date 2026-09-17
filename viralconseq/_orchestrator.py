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
from viralconseq.provenance import (
    SNAKEMAKE_LOG_COPY,
    copy_snakemake_log,
    record_run_completion,
    write_run_manifest,
)

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


def run_dir_for(args: Dict[str, Any]) -> Optional[str]:
    """``<output>/<run_name>`` (absolute), or ``None`` when args lack either.

    This is Snakemake's working directory for the run, so ``.snakemake/``
    (locks, metadata, per-run logs) lives next to the results instead of in
    whatever directory the CLI was launched from.
    """
    output = args.get("output")
    run_name = args.get("run_name")
    if not output or not run_name:
        return None
    return os.path.abspath(os.path.join(output, run_name))


def run_workflow(workflow_path: str, args: Dict[str, Any]) -> bool:
    """Run a Snakemake workflow with the kwargs every workflow shares.

    Args:
        workflow_path: Absolute path to the ``.smk`` file to execute.
        args: Pipeline argument dict. Must contain ``config_file`` and
            ``threads_total``; ``output`` and ``run_name`` select the run
            directory used as Snakemake's working directory.

    Raises:
        ValidationError: If ``workflow_path`` does not exist.
    """
    logger.info("Starting Snakemake workflow")

    if not os.path.isfile(workflow_path):
        raise ValidationError(f"Workflow file not found: {workflow_path}")

    run_dir = run_dir_for(args)
    conda_prefix = args.get("conda_prefix")
    if conda_prefix:
        # Snakemake resolves a relative prefix against its working directory,
        # which is now the run directory rather than the caller's cwd.
        conda_prefix = os.path.abspath(os.path.expanduser(str(conda_prefix)))

    successful = snakemake(
        workflow_path,
        configfiles=[os.path.abspath(args["config_file"])],
        cores=args["threads_total"],
        use_conda=True,
        conda_prefix=conda_prefix,
        # Run inside the run directory: every path in the config is absolute
        # (resolve_path_args + absolutise_sample_paths), so only .snakemake/
        # moves. force_incomplete resumes a run that was interrupted mid-job
        # instead of stopping with IncompleteFilesException.
        workdir=run_dir,
        force_incomplete=True,
        targets=["all"],
    )

    if run_dir is not None:
        copied = copy_snakemake_log(run_dir)
        if copied:
            logger.info(f"Snakemake log copied to {copied}")

    if successful:
        logger.info("Snakemake workflow completed successfully")
    else:
        logger.error("Snakemake workflow failed")

    return successful


def run_pipeline(
    args: Dict[str, Any],
    *,
    resolve_paths: Callable[[Dict[str, Any]], object],
    validate: Callable[[Dict[str, Any]], Dict[str, list]],
    generate_config: Callable[[Dict[str, list], Dict[str, Any]], None],
    run_workflow_fn: Callable[[Dict[str, Any]], bool],
) -> int:
    """The try/except ``main`` skeleton of the pipeline entry point.

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

    Returns:
        Exit code (0 for success, 1 for failure).
    """
    try:
        resolve_paths(args)
        samples = validate(args)

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
                extra: Dict[str, Any] = {}
                run_dir = run_dir_for(args)
                if run_dir is not None:
                    log_copy = os.path.join(run_dir, SNAKEMAKE_LOG_COPY)
                    if os.path.isfile(log_copy):
                        extra["snakemake_log"] = log_copy
                record_run_completion(
                    manifest_path, status="success" if successful else "failed", extra=extra
                )
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
