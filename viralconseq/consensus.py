#!/usr/bin/env python

"""Consensus pipeline entry point.

Validates the CLI arguments, writes the Snakemake config YAML and launches one
of the four consensus workflows (illumina/nanopore, single/segmented reference).
"""

import logging
import os
from typing import Any, Dict

from viralconseq import _orchestrator
from viralconseq.constants import Clair3Models, DataType, ResourceDefaults
from viralconseq.validators import (
    CONSENSUS_PATH_ARG_KEYS,
    get_samples_from_args,
    resolve_path_args,
    resolve_resource_budget,
    sanitize_identifier,
    validate_clair3_model,
    validate_config_dict,
    validate_consensus_input_integrity,
    validate_consensus_requirements,
    validate_flag_strings,
    validate_illumina_requirements,
    validate_numeric_parameters,
    validate_viralqc_database,
)

# Set up logging
logger = logging.getLogger(__name__)


def rule_list_for(args: Dict[str, Any]) -> list:
    """The ``ResourceDefaults`` rule list for the run's data type."""
    if args.get("data_type") == DataType.NANOPORE:
        return ResourceDefaults.CONSENSUS_NANOPORE_RULES
    return ResourceDefaults.CONSENSUS_ILLUMINA_RULES


def validate_args(args: Dict[str, Any]) -> Dict[str, list]:
    """Validate all pipeline arguments.

    Args:
        args: Dictionary of pipeline arguments

    Returns:
        Dictionary mapping sample names to file paths

    Raises:
        ValidationError: If validation fails
    """
    logger.info("Validating pipeline arguments")

    # Reject unsafe run names before they become output-path / DAG components
    # (run_name is spliced into config["output"] and every rule's output path).
    if args.get("run_name"):
        sanitize_identifier(args["run_name"], "run_name")

    # Range-check numeric parameters (threads, thresholds, coverage, ...).
    validate_numeric_parameters(args)
    # Extra tool flags are spliced unquoted into shell commands.
    validate_flag_strings(args)
    # Memory budget: explicit --max-memory, else detected from the machine.
    args["memory_budget_source"] = resolve_resource_budget(args, rule_list_for(args))

    # Get and validate samples
    samples = get_samples_from_args(args)

    logger.info(f"Found {len(samples)} samples")

    # Validate consensus-specific requirements
    validate_consensus_requirements(args)

    # Handle primer scheme - set to "NA" if not provided
    if not args.get("primer_scheme"):
        logger.info("A primer scheme was not provided (untargeted sequencing).")
        args["primer_scheme"] = "NA"
    else:
        logger.info("A primer scheme was provided (Amplicon sequencing)...")

    # Validate data-type specific requirements
    validate_illumina_requirements(args)

    # viralQC databases must exist before we write a config that points at
    # them. Existence-level check: not gated by --skip-input-validation and
    # also run under --create-config-only; --no-run-viralqc skips it.
    validate_viralqc_database(args)

    # ``validate_consensus_requirements`` parses ``--segmented-reference``
    # (``L=/path/L.fasta``) into a dict stored under ``reference``. Now that
    # the dict exists, resolve its values to absolute paths.
    if isinstance(args.get("reference"), dict):
        resolve_path_args(args, ("reference",))

    # Content-level integrity of the actual input files (FASTQ/FASTA/BED),
    # now that references are split and all paths are absolute. Gated by
    # --skip-input-validation. Kept as a module-level name so it can be patched
    # in tests alongside validate_consensus_requirements.
    validate_consensus_input_integrity(args, samples)

    # Nanopore: choose (or check) the Clair3 model per sample and make sure it
    # is on disk. Last, so the reads have already been validated.
    validate_clair3_model(args, samples)

    logger.info("All arguments validated successfully")
    return samples


def generate_config_file(samples: Dict[str, list], args: Dict[str, Any]) -> None:
    """Generate configuration file for the pipeline.

    Args:
        samples: Dictionary mapping sample names to file paths
        args: Dictionary of pipeline arguments

    Raises:
        ValidationError: If config generation fails
    """
    data_type = args["data_type"]
    generator = _orchestrator.start_config(args, samples)

    # Add consensus-specific settings
    generator.add_consensus_settings(
        reference=args["reference"],
        primer_scheme=args.get("primer_scheme", "NA"),
        minimum_coverage=args.get("minimum_coverage", 20),
        minimap2_consensus_align_flags=args.get(
            "minimap2_consensus_align_flags",
            "-a --sam-hit-only --secondary=no --score-N=0",
        ),
    )

    # Add workflow_path (consensus-specific). Derive from this module's location
    # rather than sys.path[0], which is the entry-point script's directory and is
    # wrong when viralconseq is imported and called in-process (e.g. by a service).
    generator.add_workflow_path(os.path.join(os.path.abspath(os.path.dirname(__file__)), "scripts"))

    # Add Nanopore-specific settings if needed
    if data_type == DataType.NANOPORE:
        generator.add_consensus_nanopore_settings(
            minimum_read_length=args.get("minimum_read_length", 50),
            af_threshold=args.get("af_threshold", 0.51),
            chunk_size=args.get("chunk_size", 10000),
            clair3_model=args.get("clair3_model", Clair3Models.AUTO),
            clair3_model_dir=args.get("clair3_model_dir") or Clair3Models.default_dir(),
            variant_quality=args.get("variant_quality", 20),
            variant_depth=args.get("variant_depth", 10),
            minimum_map_quality=args.get("minimum_map_quality", 30),
        )

    # Add Illumina-specific settings if needed
    if data_type == DataType.ILLUMINA:
        generator.add_illumina_settings(
            adapters=args["adapters"],
            minimum_read_length=args.get("minimum_read_length", 50),
            trim_head=args.get("trim_head", 0),
            trim_tail=args.get("trim_tail", 0),
            cut_front_mean_quality=args.get("cut_front_mean_quality", 10),
            cut_tail_mean_quality=args.get("cut_tail_mean_quality", 10),
            cut_right_window_size=args.get("cut_right_window_size", 4),
            cut_right_mean_quality=args.get("cut_right_mean_quality", 15),
            af_threshold=args.get("af_threshold", 0.51),
            af_isnv_threshold=args.get("af_isnv_threshold", 0),
            run_isnv=args.get("run_isnv", False),
        )

    # Consensus QC (viralQC) - shared by both data types
    generator.add_viralqc_settings(
        run_viralqc=args.get("run_viralqc", True),
        viralqc_db=args.get("viralqc_db") or "NA",
        viralqc_extra_flags=args.get("viralqc_extra_flags", ""),
    )

    # Add resource settings
    if data_type == DataType.ILLUMINA:
        generator.add_resource_settings(args, ResourceDefaults.CONSENSUS_ILLUMINA_RULES)
    else:
        generator.add_resource_settings(args, ResourceDefaults.CONSENSUS_NANOPORE_RULES)

    # Save config file
    # Self-check the contract before writing: the same check ``viralconseq
    # rerun`` and the Snakefile apply to a saved YAML.
    validate_config_dict(generator.config)
    generator.save()

    logger.info(f"Configuration file generated: {args['config_file']}")


def workflow_path_for(data_type: str, reference: Any) -> str:
    """Absolute path of the workflow file for a data type and reference shape.

    Segmentation is selected by the reference being a mapping (segment ->
    FASTA) rather than a single path.
    """
    thisdir = os.path.abspath(os.path.dirname(__file__))
    segmented_suffix = "_segmented" if isinstance(reference, dict) else ""
    return os.path.join(thisdir, "scripts", f"consensus_{data_type}{segmented_suffix}.smk")


def run_snakemake_workflow(args: Dict[str, Any]) -> bool:
    """Run the Snakemake workflow for the consensus pipeline."""
    workflow_path = workflow_path_for(args["data_type"], args.get("reference"))
    return _orchestrator.run_workflow(workflow_path, args)


def main(args: Dict[str, Any]) -> int:
    """Main entry point for the consensus pipeline.

    Returns:
        Exit code (0 for success, 1 for failure).
    """
    return _orchestrator.run_pipeline(
        args,
        resolve_paths=lambda a: resolve_path_args(a, CONSENSUS_PATH_ARG_KEYS),
        validate=validate_args,
        generate_config=generate_config_file,
        run_workflow_fn=run_snakemake_workflow,
    )
