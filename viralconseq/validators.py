"""Validation functions for viralconseq pipeline arguments and data."""

import csv
import glob
import logging
import os
import re
import shlex
from collections import Counter
from typing import Any, Dict, List, Optional, cast

from viralconseq import clair3_models, integrity
from viralconseq.constants import (
    Clair3Models,
    DataType,
    ResourceDefaults,
    ViralQCDatabase,
    detect_memory_mb,
    memory_budget_mb,
)
from viralconseq.exceptions import (
    AdaptersNotFoundError,
    Clair3ModelNotFoundError,
    Clair3ModelUnresolvedError,
    ConfigurationError,
    InputIntegrityError,
    PrimerSchemeNotFoundError,
    ReferenceNotFoundError,
    SampleConfigurationNotFoundError,
    SampleSheetError,
    ValidationError,
    ViralConseqFileNotFoundError,
    ViralQCDatabaseNotFoundError,
)
from viralconseq.reference_splitter import count_records, split_multifasta

logger = logging.getLogger(__name__)


def validate_file_exists(file_path: str, description: str = "File") -> None:
    """Validate that a file exists.

    Args:
        file_path: Path to the file
        description: Description of the file for error messages

    Raises:
        ViralConseqFileNotFoundError: If the file does not exist
    """
    if not os.path.isfile(file_path):
        raise ViralConseqFileNotFoundError(f"{description} does not exist: {file_path}")


# ---------------------------------------------------------------------------
# Untrusted-input sanitization
# ---------------------------------------------------------------------------
#
# Sample ids and run names become shell tokens and filesystem path components
# inside the Snakemake rules (``sample-<id>``, ``<output>/<run_name>/...``).
# When viralconseq is embedded in a service that ingests uploaded data, these
# must be constrained so they cannot inject path traversal or shell
# metacharacters. The helpers here are opt-in checks the CLI/service layer can
# call; they intentionally allow only a conservative identifier charset.

_SAFE_IDENTIFIER_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def sanitize_identifier(value: Any, field: str = "identifier") -> str:
    """Validate that *value* is a safe filesystem/shell identifier.

    Allows letters, digits, ``.``, ``_`` and ``-`` (and must start with a
    letter/digit). Rejects empty values, path separators, ``..``, NUL, and any
    other metacharacter.

    Returns:
        The stripped value.

    Raises:
        ValidationError: If the value is not a safe identifier.
    """
    if value is None or not str(value).strip():
        raise ValidationError(f"{field} must be a non-empty string.")
    stripped = str(value).strip()
    if stripped in (".", "..") or not _SAFE_IDENTIFIER_RE.match(stripped):
        raise ValidationError(
            f"{field} may only contain letters, digits, '.', '_', '-' and must not "
            f"contain path separators or spaces: {value!r}"
        )
    return stripped


def ensure_within_base(path: str, base: str) -> str:
    """Resolve *path* against *base* and ensure it does not escape *base*.

    Absolute paths are resolved as-is; relative paths are joined to *base*.

    Returns:
        The absolute, normalized target path.

    Raises:
        ValidationError: If the resolved path is outside *base*.
    """
    base_abs = os.path.abspath(base)
    target = os.path.abspath(path if os.path.isabs(path) else os.path.join(base_abs, path))
    if os.path.commonpath([base_abs, target]) != base_abs:
        raise ValidationError(f"Path escapes base directory {base!r}: {path!r}")
    return target


# Objective validity bounds for numeric parameters. These are correctness
# constraints (a thread count must be >= 1; an allele frequency is a fraction in
# [0, 1]), NOT scientific tuning choices — analysis knobs such as
# z_score_threshold / log10_ratio_threshold are intentionally left unbounded.
_NUMERIC_BOUNDS = {
    "threads": (1, None),
    "threads_total": (1, None),
    "max_memory": (0, None),
    "minimum_mapped_reads": (0, None),
    "minimum_coverage": (1, None),
    "minimum_depth": (1, None),
    "minimum_length": (0, None),
    "minimum_read_length": (0, None),
    "af_threshold": (0.0, 1.0),
    "af_isnv_threshold": (0.0, 1.0),
}


# Config keys whose value is spliced unquoted into a rule's shell command line
# as extra tool flags. They are config-only (no CLI option), so the check runs
# both here (embedding API) and at Snakefile parse time (hand-edited YAML).
FLAG_STRING_KEYS = ("viralqc_extra_flags", "minimap2_consensus_align_flags")
_SHELL_META_RE = re.compile(r"[;&|<>`$\\\n]")


def validate_flag_strings(args: Dict[str, Any]) -> None:
    """Reject flag strings that would not survive, or would escape, the shell.

    A flag string must split cleanly with ``shlex`` and contain no shell
    metacharacters (``; & | < > ` $ \\`` or a newline): the rules interpolate
    these values unquoted, so anything beyond plain tokens is a command
    injection, not a flag.

    Raises:
        ValidationError: On an unbalanced quote or a metacharacter.
    """
    for key in FLAG_STRING_KEYS:
        value = args.get(key)
        if value is None:
            continue
        if not isinstance(value, str):
            raise ValidationError(f"{key} must be a string, got {value!r}.")
        try:
            shlex.split(value)
        except ValueError as e:
            raise ValidationError(f"{key} is not a valid flag string ({e}): {value!r}") from e
        if _SHELL_META_RE.search(value):
            raise ValidationError(
                f"{key} contains a shell metacharacter (one of ; & | < > ` $ \\ or a newline): "
                f"{value!r}. Pass plain tool flags only."
            )


# Keys every generated config carries, by data type. Mirrored by
# ``_REQUIRED_KEYS`` in rules/common.smk (the Snakefile fallback when this
# package is not importable) and by test/rule_inventory_test.py.
CONFIG_REQUIRED_COMMON = (
    "samples",
    "data",
    "output",
    "threads",
    "reference",
    "scheme",
    "minimum_depth",
    "minimum_length",
    "af_threshold",
)
CONFIG_REQUIRED_ILLUMINA = (
    "adapters",
    "trim_head",
    "trim_tail",
    "cut_front_mean_quality",
    "cut_tail_mean_quality",
    "cut_right_window_size",
    "cut_right_mean_quality",
    "af_isnv_threshold",
)
CONFIG_REQUIRED_NANOPORE = (
    "chunk_size",
    "clair3_model",
    "variant_quality",
    "variant_depth",
    "minimum_map_quality",
    "minimum_mapped_reads",
)
# key -> (lower bound, upper bound, integer required)
_CONFIG_NUMERIC = {
    "threads": (1, None, True),
    "threads_total": (1, None, True),
    "minimum_depth": (1, None, True),
    "minimum_length": (0, None, True),
    "af_threshold": (0.0, 1.0, False),
    "af_isnv_threshold": (0.0, 1.0, False),
    "trim_head": (0, None, True),
    "trim_tail": (0, None, True),
    "cut_front_mean_quality": (0, None, True),
    "cut_tail_mean_quality": (0, None, True),
    "cut_right_window_size": (1, None, True),
    "cut_right_mean_quality": (0, None, True),
    "chunk_size": (1, None, True),
    "variant_quality": (0, None, True),
    "variant_depth": (0, None, True),
    "minimum_map_quality": (0, None, True),
    "minimum_mapped_reads": (0, None, True),
    "max_memory_mb": (0, None, True),
    "memory_detected_mb": (0, None, True),
}


def validate_config_dict(config: Any) -> None:
    """Check a Snakemake config mapping the way the rules will read it.

    Used on the dict ``ConfigGenerator`` is about to write, by ``viralconseq
    rerun`` on a saved (possibly hand-edited) YAML, and at Snakefile parse time
    through ``rules/common.smk``. Standard library only; no coercion: a value
    with the wrong type is an error, not something to repair.

    Raises:
        ConfigurationError: Naming the first offending key.
    """

    def fail(message: str) -> None:
        raise ConfigurationError(f"invalid config: {message}")

    if not isinstance(config, dict):
        fail("expected a mapping at the top level")
    data = config.get("data")
    if data not in (DataType.ILLUMINA, DataType.NANOPORE):
        fail(f"data must be 'illumina' or 'nanopore', got {data!r}")

    required = list(CONFIG_REQUIRED_COMMON)
    required += CONFIG_REQUIRED_ILLUMINA if data == DataType.ILLUMINA else CONFIG_REQUIRED_NANOPORE
    if data == DataType.NANOPORE:
        required.extend(["infer_consensus_sequence_ram", "clair3_model_dir"])
    if config.get("run_viralqc", True):
        required.extend(["viralqc_db", "run_viralqc_ram"])
    missing = [key for key in required if key not in config]
    if missing:
        fail(
            "missing required key(s): "
            + ", ".join(missing)
            + " - regenerate the file with 'viralconseq consensus ... --create-config-only'"
        )

    samples = config["samples"]
    if not isinstance(samples, dict) or not samples:
        fail("samples must be a non-empty mapping of sample id -> FASTQ path(s)")
    expected_files = 2 if data == DataType.ILLUMINA else 1
    for sample, paths in samples.items():
        if not isinstance(sample, str) or not sample.strip():
            fail(f"sample id {sample!r} is not a non-empty string")
        if isinstance(paths, str):
            paths = paths.split()
        if not isinstance(paths, list) or len(paths) != expected_files:
            fail(
                f"sample {sample!r} must list exactly {expected_files} FASTQ path(s), got {paths!r}"
            )
        if not all(isinstance(p, str) and p.strip() for p in paths):
            fail(f"sample {sample!r} has an empty FASTQ path")

    reference = config["reference"]
    if isinstance(reference, dict):
        if not reference or not all(isinstance(v, str) and v.strip() for v in reference.values()):
            fail("reference must map every segment to a FASTA path")
    elif not isinstance(reference, str) or not reference.strip():
        fail("reference must be a FASTA path or a mapping of segment -> FASTA path")

    for key in ("output", "scheme"):
        if not isinstance(config[key], str) or not config[key].strip():
            fail(f"{key} must be a non-empty string")
    if data == DataType.ILLUMINA and not isinstance(config["adapters"], str):
        fail("adapters must be a path or 'NA'")
    if data == DataType.NANOPORE:
        model = config["clair3_model"]
        names = list(model.values()) if isinstance(model, dict) else [model]
        if isinstance(model, dict) and set(model) - set(samples):
            fail("clair3_model names samples that are not in samples")
        for name in names:
            if not isinstance(name, str) or not Clair3Models.NAME_RE.match(name):
                fail(f"clair3_model {name!r} is not a valid model name")
            if name.endswith(Clair3Models.DWELL_TIME_SUFFIX):
                fail(f"clair3_model {name!r} is a move-table model viralconseq cannot run")
        if (
            not isinstance(config["clair3_model_dir"], str)
            or not config["clair3_model_dir"].strip()
        ):
            fail("clair3_model_dir must be a directory path")

    for key in ("run_isnv", "run_viralqc"):
        if key in config and not isinstance(config[key], bool):
            fail(f"{key} must be true or false, got {config[key]!r}")

    def check_number(key: str, lo: Any, hi: Any, integer: bool) -> None:
        value = config[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            fail(f"{key} must be a number, got {value!r}")
        if integer and not isinstance(value, int):
            fail(f"{key} must be an integer, got {value!r}")
        if lo is not None and value < lo:
            fail(f"{key} must be >= {lo}, got {value}")
        if hi is not None and value > hi:
            fail(f"{key} must be <= {hi}, got {value}")

    for key, (lo, hi, integer) in _CONFIG_NUMERIC.items():
        if key in config:
            check_number(key, lo, hi, integer)
    for key in config:
        if isinstance(key, str) and (key.endswith("_cpus") or key.endswith("_ram")):
            check_number(key, 1, None, True)

    try:
        validate_flag_strings(config)
    except ValidationError as exc:
        fail(str(exc))

    budget = config.get("max_memory_mb", 0)
    if isinstance(budget, int) and budget > 0:
        for rule in ResourceDefaults.MEMORY_RULES:
            ram_key = f"{rule}_ram"
            if ram_key in config and budget < int(config[ram_key]) * 1024:
                fail(
                    f"max_memory_mb ({budget}) is below {ram_key} ({config[ram_key]} GB); "
                    "raise the budget or set max_memory_mb: 0"
                )


def resolve_resource_budget(args: Dict[str, Any], rule_list: List[str]) -> str:
    """Fill ``args["max_memory_mb"]`` and ``args["memory_detected_mb"]``.

    An explicit ``--max-memory`` (0 included) wins; it is refused when it is
    positive but below the largest ``<rule>_ram`` in the run, because Snakemake
    would then clamp every memory-declaring job down to the budget and
    serialise the run without saying so. Otherwise the budget is the detected
    memory less ``ResourceDefaults.MEMORY_HEADROOM``, raised to the largest
    per-rule figure when the machine is smaller than one such job.

    Returns:
        A short label saying where the budget came from (for the run summary).

    Raises:
        ValidationError: On an explicit budget below the largest rule figure.
    """
    memory_rules = ResourceDefaults.memory_rules_for(rule_list)
    largest_gb = max(
        (int(args.get(f"{rule}_ram") or ResourceDefaults.ram_for(rule)) for rule in memory_rules),
        default=0,
    )
    explicit = args.get("max_memory")
    if explicit is not None:
        budget = int(explicit) * 1024
        if 0 < budget < largest_gb * 1024:
            raise ValidationError(
                f"--max-memory {explicit} GB is below the largest per-rule figure "
                f"({largest_gb} GB); raise it, or pass 0 to run without a budget."
            )
        args["memory_detected_mb"] = detect_memory_mb() or 0
        args["max_memory_mb"] = budget
        return "--max-memory"
    detected = detect_memory_mb()
    args["memory_detected_mb"] = detected or 0
    if detected is None:
        args["max_memory_mb"] = 0
        logger.warning(
            "Could not detect the memory available to this process; running without a "
            "memory budget (as --max-memory 0 would). Pass --max-memory to set one."
        )
        return "undetectable"
    budget, clamped = memory_budget_mb(detected, largest_gb)
    args["max_memory_mb"] = budget
    if clamped:
        logger.warning(
            f"Detected {detected / 1024:.1f} GB usable, below the largest per-rule figure "
            f"({largest_gb} GB); budget set to {largest_gb} GB so memory-declaring rules run "
            "one at a time. Pass --max-memory to override."
        )
        return f"detected {detected / 1024:.1f} GB, clamped up to one job"
    return f"detected {detected / 1024:.1f} GB, {ResourceDefaults.MEMORY_HEADROOM:.0%} headroom"


def validate_numeric_parameters(args: Dict[str, Any]) -> None:
    """Range-check numeric parameters that have objective validity bounds.

    Only keys present in ``args`` are checked. Raises on the first out-of-range
    value.

    Raises:
        ValidationError: If a present numeric parameter is non-numeric or falls
            outside its valid range.
    """

    def _as_number(key: str) -> Optional[float]:
        if key not in args or args[key] is None:
            return None
        val = args[key]
        # bool is a subclass of int; reject it explicitly.
        if isinstance(val, bool) or not isinstance(val, (int, float)):
            raise ValidationError(f"{key} must be a number, got {val!r}.")
        return val

    for key, (lo, hi) in _NUMERIC_BOUNDS.items():
        val = _as_number(key)
        if val is None:
            continue
        if lo is not None and val < lo:
            raise ValidationError(f"{key} must be >= {lo}, got {val}.")
        if hi is not None and val > hi:
            raise ValidationError(f"{key} must be <= {hi}, got {val}.")


def validate_sample_sheet(sample_sheet_path: str, data_type: str) -> Dict[str, List[str]]:
    """Validate and parse sample sheet file.

    Args:
        sample_sheet_path: Path to the sample sheet CSV file
        data_type: Type of sequencing data (illumina or nanopore)

    Returns:
        Dictionary mapping sample names to file paths

    Raises:
        SampleSheetError: If the sample sheet is invalid
        ViralConseqFileNotFoundError: If the sample sheet file does not exist
    """
    validate_file_exists(sample_sheet_path, "Sample sheet file")

    expected_columns = 3 if data_type == DataType.ILLUMINA else 2

    try:
        with open(sample_sheet_path, newline="", encoding="utf-8") as handle:
            rows = list(csv.reader(handle))
    except (OSError, UnicodeDecodeError) as e:
        raise SampleSheetError(f"Failed to read sample sheet: {e}") from e

    samples: Dict[str, List[str]] = {}

    for line_number, row in enumerate(rows, start=1):
        # Skip blank lines (rows with no cells, or only empty cells).
        if not row or all(not cell.strip() for cell in row):
            continue

        # Tolerate trailing empty fields (e.g. a stray trailing comma) but
        # reject rows whose real column count does not match the data type,
        # instead of silently NaN-padding them the way pandas did.
        trimmed = list(row)
        while trimmed and not trimmed[-1].strip():
            trimmed.pop()

        if len(trimmed) != expected_columns:
            raise SampleSheetError(
                f"{str(data_type).capitalize()} sample sheet requires exactly "
                f"{expected_columns} columns; row {line_number} has {len(trimmed)}: {row!r}"
            )

        sample_name = trimmed[0].strip()
        if not sample_name:
            raise SampleSheetError(f"Empty sample name on row {line_number}.")

        # Sample ids become shell tokens and path components in the workflow
        # (``sample-<id>``); reject anything that could inject a path or
        # metacharacter when the sheet comes from untrusted input.
        if sample_name in (".", "..") or not _SAFE_IDENTIFIER_RE.match(sample_name):
            raise SampleSheetError(
                f"Unsafe sample id '{sample_name}' on row {line_number}: use only letters, "
                f"digits, '.', '_', '-' (no path separators or spaces)."
            )

        if sample_name in samples:
            raise SampleSheetError(
                f"Duplicate sample id '{sample_name}' on row {line_number}. "
                f"Sample ids must be unique so no sample is silently dropped."
            )

        file_paths = [cell.strip() for cell in trimmed[1:]]
        for offset, file_path in enumerate(file_paths, start=2):
            if not file_path:
                raise SampleSheetError(
                    f"Missing file path in column {offset} for sample "
                    f"'{sample_name}' on row {line_number}."
                )

        if data_type == DataType.ILLUMINA:
            validate_file_exists(file_paths[0], f"R1 file for sample {sample_name}")
            validate_file_exists(file_paths[1], f"R2 file for sample {sample_name}")
        else:  # nanopore
            validate_file_exists(file_paths[0], f"File for sample {sample_name}")

        samples[sample_name] = file_paths

    if not samples:
        raise SampleSheetError("No valid samples found in sample sheet")

    return samples


def validate_illumina_requirements(args: Dict[str, Any]) -> None:
    """Validate Illumina-specific requirements (fastp: adapters optional)."""
    if args.get("data_type") != DataType.ILLUMINA:
        return

    adapters = args.get("adapters")
    if adapters and str(adapters).strip() and str(adapters).strip() != "NA":
        try:
            validate_file_exists(adapters, "Illumina adapter sequences file")
        except ViralConseqFileNotFoundError as e:
            raise AdaptersNotFoundError(
                f"Illumina adapter sequences file does not exist: {e}"
            ) from e


def _segment_input_dir(args: Dict[str, Any]) -> str:
    """Directory where auto-split per-segment inputs are written.

    Mirrors ``ConfigGenerator.add_output``'s ``output/run_name/`` layout so the
    generated per-segment references live inside the run's directory.
    """
    return os.path.join(args["output"], args.get("run_name") or "", "input_references")


def _maybe_split_multifasta_reference(args: Dict[str, Any]) -> bool:
    """Auto-split a single multi-record ``--reference`` into per-segment files.

    A ``--reference`` FASTA with more than one record is treated as a segmented
    virus: it is split into one file per record and ``args["reference"]`` is
    replaced with the ``{segment: path}`` mapping the segmented workflows
    consume. ``--single-reference`` forces the historical single-reference
    behaviour (all records aligned together in one pass).

    Returns:
        True if the reference was split (segmented mode engaged), else False.
    """
    reference = args.get("reference")
    if args.get("single_reference"):
        return False
    if not isinstance(reference, str):
        return False
    if count_records(reference) <= 1:
        return False

    try:
        segment_map = split_multifasta(reference, _segment_input_dir(args))
    except ValueError as e:
        raise ValidationError(str(e)) from e
    logger.info(
        "Reference '%s' contains %d records; running in segmented mode with " "segments: %s",
        reference,
        len(segment_map),
        ", ".join(segment_map),
    )
    args["reference"] = segment_map
    return True


def validate_consensus_requirements(args: Dict[str, Any]) -> None:
    """Validate consensus pipeline requirements.

    Exactly one of --reference or --segmented-reference must be provided.
    When --segmented-reference is used, the values are parsed from
    SEGMENT=PATH format and stored as a dict in args["reference"].

    Args:
        args: Dictionary of pipeline arguments

    Raises:
        ValidationError: If consensus requirements are not met
    """
    reference = args.get("reference")
    segmented_reference = args.get("segmented_reference")

    if reference and segmented_reference:
        raise ValidationError(
            "--reference and --segmented-reference are mutually exclusive. "
            "Please provide only one."
        )

    if not reference and not segmented_reference:
        raise ValidationError(
            "A reference is required. Provide --reference for a single reference "
            "or --segmented-reference for segmented viruses."
        )

    if segmented_reference:
        if isinstance(segmented_reference, dict):
            parsed_segments = segmented_reference
        else:
            parsed_segments = {}
            for entry in segmented_reference:
                if "=" not in entry:
                    raise ValidationError(
                        f"Invalid segmented reference format: '{entry}'. "
                        f"Expected SEGMENT=PATH (e.g. S=/path/to/S.fasta)"
                    )
                segment_name, segment_path = entry.split("=", 1)
                segment_name = segment_name.strip()
                segment_path = segment_path.strip()
                if not segment_name or not segment_path:
                    raise ValidationError(
                        f"Invalid segmented reference format: '{entry}'. "
                        f"Both segment name and path are required."
                    )
                parsed_segments[segment_name] = segment_path

        args["reference"] = parsed_segments
        args["segmented_reference"] = None
        reference = parsed_segments

    if isinstance(reference, dict):
        for segment_name, segment_path in reference.items():
            try:
                validate_file_exists(segment_path, f"Reference file for segment '{segment_name}'")
            except ViralConseqFileNotFoundError as e:
                raise ReferenceNotFoundError(str(e)) from e
    else:
        try:
            validate_file_exists(cast(str, reference), "Reference sequence file")
        except ViralConseqFileNotFoundError as e:
            raise ReferenceNotFoundError(f"Reference sequence file does not exist: {e}") from e
        _maybe_split_multifasta_reference(args)

    primer_scheme = args.get("primer_scheme")
    if primer_scheme:
        try:
            validate_file_exists(primer_scheme, "Primer scheme file")
        except ViralConseqFileNotFoundError as e:
            raise PrimerSchemeNotFoundError(f"Primer scheme file does not exist: {e}") from e


# ---------------------------------------------------------------------------
# Content-level input integrity (consensus)
# ---------------------------------------------------------------------------


def _is_real_path(value: Any) -> bool:
    """True if *value* is a usable path string (not None/empty/the ``NA`` sentinel)."""
    return isinstance(value, str) and value.strip() not in ("", "NA")


def _reference_paths(reference: Any) -> List[str]:
    """Return the reference FASTA path(s), whether single or per-segment dict."""
    if isinstance(reference, dict):
        return [p for p in reference.values() if _is_real_path(p)]
    if _is_real_path(reference):
        return [reference]
    return []


def validate_consensus_input_integrity(args: Dict[str, Any], samples: Dict[str, List[str]]) -> None:
    """Run content-level integrity checks on the consensus pipeline's inputs.

    Streams every sample FASTQ, the reference FASTA(s) and the primer BED (if
    any), aggregating problems across all files. Blocking (``error``) issues are
    raised together as a single :class:`InputIntegrityError`; ``warning`` issues
    (e.g. an odd BED strand column) are logged. Cross-file BED chrom matching is
    data-type aware: the nanopore pipeline sanitizes reference headers before
    mapping, so the expected BAM contig names differ from Illumina's raw ids.

    Skipped entirely when ``args['skip_input_validation']`` is set. Assumes file
    existence has already been checked (so failures here are about *content*).

    Raises:
        InputIntegrityError: If any blocking issue is found.
    """
    if args.get("skip_input_validation"):
        logger.info("Input integrity validation skipped (--skip-input-validation).")
        return

    logger.info("Running input integrity checks")
    reports: List[integrity.IntegrityReport] = []

    # FASTQ: every read file of every sample.
    for file_paths in samples.values():
        for file_path in file_paths:
            reports.append(integrity.validate_fastq(file_path))

    # Reference FASTA(s): collect the headers so we can derive the exact contig
    # names the pipeline will map against.
    all_headers: List[str] = []
    for ref_path in _reference_paths(args.get("reference")):
        report = integrity.validate_fasta(ref_path)
        reports.append(report)
        all_headers.extend(report.headers)

    data_type = args.get("data_type")
    if data_type == DataType.NANOPORE:
        expected_contigs = {integrity.sanitize_nanopore_contig(h) for h in all_headers}
        accession_map = {
            integrity.header_token(h): integrity.sanitize_nanopore_contig(h) for h in all_headers
        }
    else:
        expected_contigs = {integrity.header_token(h) for h in all_headers}
        accession_map = None

    # Primer BED (critical when provided): chrom names must match the reference.
    primer_scheme = args.get("primer_scheme")
    if isinstance(primer_scheme, str) and primer_scheme.strip() not in ("", "NA"):
        reports.append(
            integrity.validate_bed(
                primer_scheme,
                expected_contigs=expected_contigs or None,
                accession_map=accession_map or None,
            )
        )

    errors = [(report.path, issue) for report in reports for issue in report.errors]
    warnings = [(report.path, issue) for report in reports for issue in report.warnings]
    if args.get("data_type") == DataType.NANOPORE:
        # A barcode with zero reads is the everyday case the nanopore workflow
        # degrades per sample (all-N consensus, status no_mapped_reads); it must
        # not stop a 96-sample run. Illumina keeps the error: fastp/minimap2
        # have no such path.
        demoted = [(p, i) for p, i in errors if i.code == "fastq_empty"]
        errors = [(p, i) for p, i in errors if i.code != "fastq_empty"]
        warnings.extend(demoted)

    for path, issue in warnings:
        loc = f" (line {issue.line})" if issue.line else ""
        logger.warning(
            "Input integrity warning in %s [%s]%s: %s",
            os.path.basename(path),
            issue.code,
            loc,
            issue.message,
        )

    if errors:
        grouped: Dict[str, List[integrity.IntegrityIssue]] = {}
        for path, issue in errors:
            grouped.setdefault(path, []).append(issue)
        lines = [f"Input integrity check failed with {len(errors)} problem(s):"]
        for path, issues in grouped.items():
            lines.append(f"  {path}:")
            for issue in issues:
                loc = f" (line {issue.line})" if issue.line else ""
                lines.append(f"    - [{issue.code}]{loc} {issue.message}")
        raise InputIntegrityError("\n".join(lines), issues=[issue for _, issue in errors])

    logger.info("Input integrity checks passed")


# ---------------------------------------------------------------------------
# viralQC database
# ---------------------------------------------------------------------------


def missing_viralqc_database_files(db_dir: str) -> List[str]:
    """Return the required viralQC database entries missing from *db_dir*.

    An empty list means the database is complete; ``["<directory>"]`` means the
    directory itself does not exist. Shared with ``viralconseq setup`` (dry-run
    report, "already present" short-circuit, post-download check).
    """
    if not os.path.isdir(db_dir):
        return ["<directory>"]
    missing = [
        name
        for name in ViralQCDatabase.REQUIRED_FILES
        if not os.path.isfile(os.path.join(db_dir, name))
    ]
    missing += [
        f"{name}/"
        for name in ViralQCDatabase.REQUIRED_DIRS
        if not os.path.isdir(os.path.join(db_dir, name))
    ]
    if not any(glob.glob(os.path.join(db_dir, p)) for p in ViralQCDatabase.BLAST_INDEX_GLOBS):
        missing.append(ViralQCDatabase.BLAST_INDEX_LABEL)
    return missing


def viralqc_setup_hint(db_dir: str) -> str:
    """The exact command that populates *db_dir* (used in errors and docs)."""
    return f"viralconseq setup --viralqc-db {db_dir}"


def missing_clair3_model_files(model_dir: str, name: str) -> List[str]:
    """Checkpoints missing or invalid for ``<model_dir>/<name>``; ``[]`` when usable."""
    return clair3_models.missing_checkpoints(model_dir, name)


def clair3_setup_hint(model_dir: str, name: str) -> str:
    """The command that fetches ``name`` into ``model_dir``; single source for messages."""
    return f"viralconseq setup --clair3-models {name} --clair3-model-dir {model_dir}"


def validate_clair3_model(args: Dict[str, Any], samples: Dict[str, List[str]]) -> None:
    """Resolve and check the Clair3 model(s) of a nanopore run, before any work.

    ``args["clair3_model"]`` is ``"auto"`` (default), a model name, or a mapping
    ``{sample id: name or "auto"}`` (embedding API). ``auto`` reads the basecall
    tag of each sample's reads; a sample whose FASTQ holds no reads takes the
    run's majority model (Clair3 never runs for it anyway). Every resolved
    model must be a valid name, not a move-table (``_with_mv``) model, and
    present and complete in ``args["clair3_model_dir"]``.

    On success ``args["clair3_model"]`` is a single name when all samples
    agree, else ``{sample id: name}``; ``args["clair3_model_dir"]`` is set.

    Raises:
        Clair3ModelUnresolvedError, Clair3ModelMixedError,
        Clair3ModelNotFoundError, ValidationError.
    """
    if args.get("data_type") != DataType.NANOPORE:
        return
    model_dir = args.get("clair3_model_dir") or Clair3Models.default_dir()
    requested = args.get("clair3_model") or Clair3Models.AUTO
    if isinstance(requested, dict):
        unknown = sorted(set(requested) - set(samples))
        if unknown:
            raise ValidationError(
                f"clair3_model names samples not in the sample sheet: {', '.join(unknown)}"
            )
    per_sample = clair3_models.models_by_sample(samples, requested)
    resolved = [name for name in per_sample.values() if name]
    if not resolved:
        raise Clair3ModelUnresolvedError(
            "no sample yielded a Clair3 model (every FASTQ is empty?); " f"{clair3_models.HINT}."
        )
    majority = Counter(resolved).most_common(1)[0][0]
    final = {sample: (name or majority) for sample, name in per_sample.items()}
    for name in sorted(set(final.values())):
        if not Clair3Models.NAME_RE.match(name):
            raise ValidationError(f"invalid Clair3 model name {name!r}")
        if name.endswith(Clair3Models.DWELL_TIME_SUFFIX):
            raise Clair3ModelUnresolvedError(
                f"{name} is a move-table model that needs Clair3's --enable_dwell_time and mv:B "
                f"tags in the BAM, which viralconseq does not provide; use "
                f"{name[: -len(Clair3Models.DWELL_TIME_SUFFIX)]}."
            )
        missing = missing_clair3_model_files(model_dir, name)
        if missing:
            if name in Clair3Models.MANIFEST:
                fix = f"Fetch it with:\n    {clair3_setup_hint(model_dir, name)}"
            else:
                fix = (
                    f"{name} is not in the model manifest ({Clair3Models.MANIFEST_SOURCE}); "
                    f"place its pileup.pt and full_alignment.pt under {model_dir}/{name}/ yourself."
                )
            raise Clair3ModelNotFoundError(
                f"Clair3 model {name} not found or incomplete at {model_dir}/{name} "
                f"(missing: {', '.join(missing)}).\n{fix}"
            )
    names = set(final.values())
    args["clair3_model"] = final if len(names) > 1 else names.pop()
    args["clair3_model_dir"] = model_dir
    if isinstance(args["clair3_model"], dict):
        logger.info("Clair3 models: %s", ", ".join(f"{s}={m}" for s, m in final.items()))
    else:
        logger.info("Clair3 model: %s (%s)", args["clair3_model"], model_dir)


def validate_viralqc_database(args: Dict[str, Any]) -> None:
    """Fail fast when viralQC is enabled but its databases are not on disk.

    Existence-level check: it runs regardless of ``skip_input_validation`` and
    also under ``create_config_only`` (the config would otherwise point at a
    missing database). Skipped when ``args["run_viralqc"]`` is False.

    Raises:
        ViralQCDatabaseNotFoundError: If ``viralqc_db`` is unset, is not a
            directory, or lacks any required entry.
    """
    if not args.get("run_viralqc", True):
        logger.info("viralQC disabled (--no-run-viralqc); skipping database check.")
        return
    db_dir = args.get("viralqc_db")
    if _is_path_sentinel(db_dir):
        raise ViralQCDatabaseNotFoundError(
            "viralQC is enabled but no database directory is set (viralqc_db). "
            "Pass --viralqc-db PATH, set $VIRALCONSEQ_VIRALQC_DB, or use --no-run-viralqc."
        )
    db_dir = cast(str, db_dir)
    missing = missing_viralqc_database_files(db_dir)
    if missing:
        raise ViralQCDatabaseNotFoundError(
            f"viralQC database not found or incomplete at {db_dir} "
            f"(missing: {', '.join(missing)}).\n"
            f"Download it with:\n    {viralqc_setup_hint(db_dir)}\n"
            f"or pass --no-run-viralqc to skip consensus QC."
        )
    logger.info("viralQC database found: %s", db_dir)


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------
#
# CLI arguments that point to a filesystem location. The list is kept here
# (next to the validators) so it can be reused both at validation time and
# inside ``consensus.main`` to make the paths absolute before anything else
# sees them. Resolving paths up-front means a user who runs
# ``viralconseq consensus illumina ... --reference refs/ref.fasta --config-file
# scratch/run.yml`` does NOT silently get the reference looked up under
# ``scratch/refs/...`` just because the generated config happens to live in
# ``scratch/``.

CONSENSUS_PATH_ARG_KEYS = (
    "sample_sheet",
    "config_file",
    "output",
    "reference",
    "primer_scheme",
    "adapters",
    "viralqc_db",
    "clair3_model_dir",
)


def _is_path_sentinel(value: Any) -> bool:
    """Return True if a path-typed argument value should be left untouched.

    Sentinels are: ``None``, non-string scalars (ints, bools), the empty
    string, and the literal placeholder ``"NA"`` (after stripping
    whitespace). This matches the conventions used by the existing
    validators in this module.
    """
    if value is None:
        return True
    if not isinstance(value, str):
        return True
    stripped = value.strip()
    if not stripped:
        return True
    if stripped == "NA":
        return True
    return False


def resolve_path_args(
    args: Dict[str, Any],
    keys,
    base_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """Rewrite path-typed argument values to absolute paths in place.

    Each value listed in ``keys`` that is a non-sentinel relative path string
    is replaced with ``os.path.abspath(os.path.join(base_dir, value))``.
    Absolute paths, sentinels (``None``, ``""``, ``"NA"``), and non-string
    values are left unchanged. Missing keys are ignored.

    The ``reference`` argument of the consensus pipeline can be a dict
    (segmented reference, ``segment -> path``); each value of the dict is
    resolved while preserving the keys.

    Args:
        args: Mutable dict of CLI arguments.
        keys: Iterable of argument keys whose values are filesystem paths.
        base_dir: Base directory to resolve relative paths against.
            Defaults to the current working directory at call time.

    Returns:
        The same ``args`` dict (modified in place). Returned for chaining.
    """
    base_dir = base_dir or os.getcwd()

    for key in keys:
        if key not in args:
            continue
        value = args[key]

        if isinstance(value, dict):
            args[key] = {
                seg: (
                    v
                    if _is_path_sentinel(v) or os.path.isabs(v)
                    else os.path.abspath(os.path.join(base_dir, v))
                )
                for seg, v in value.items()
            }
            continue

        if _is_path_sentinel(value):
            continue
        if os.path.isabs(value):
            continue

        args[key] = os.path.abspath(os.path.join(base_dir, value))

    return args


def absolutise_sample_paths(
    samples: Dict[str, List[str]], base_dir: Optional[str] = None
) -> Dict[str, List[str]]:
    """Return ``samples`` with every relative FASTQ path made absolute.

    Snakemake runs with the run directory as its working directory (so that
    ``.snakemake/`` lives next to the results), which means every path in the
    generated config must be absolute. Sample-sheet paths are resolved against
    the current working directory, the same base ``validate_file_exists``
    checked them against.
    """
    base_dir = base_dir or os.getcwd()
    return {
        sample: [
            path if os.path.isabs(path) else os.path.abspath(os.path.join(base_dir, path))
            for path in paths
        ]
        for sample, paths in samples.items()
    }


def get_samples_from_args(args: Dict[str, Any]) -> Dict[str, List[str]]:
    """Extract and validate samples from arguments.

    Args:
        args: Dictionary of pipeline arguments

    Returns:
        Dictionary mapping sample names to file paths

    Raises:
        ValidationError: If samples cannot be determined from arguments
    """
    sample_sheet = args.get("sample_sheet")
    samples = args.get("samples")
    data_type = args.get("data_type")

    if sample_sheet:
        # A path was provided: it must exist. Do not silently fall back to
        # `samples`, which would report a misleading "nothing provided" error
        # when the real problem is a mistyped or missing sample-sheet path.
        validate_file_exists(sample_sheet, "Sample sheet file")
        return absolutise_sample_paths(validate_sample_sheet(sample_sheet, cast(str, data_type)))
    if samples:
        return absolutise_sample_paths(samples)
    raise SampleConfigurationNotFoundError(
        "Either 'sample_sheet' or 'samples' must be provided. " f"Sample sheet path: {sample_sheet}"
    )
