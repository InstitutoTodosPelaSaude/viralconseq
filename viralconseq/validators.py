"""Validation functions for viralconseq pipeline arguments and data."""

import csv
import logging
import os
import re
from typing import Any, Dict, List, Optional, cast

from viralconseq import integrity
from viralconseq.constants import DataType
from viralconseq.exceptions import (
    AdaptersNotFoundError,
    InputIntegrityError,
    PrimerSchemeNotFoundError,
    ReferenceNotFoundError,
    SampleConfigurationNotFoundError,
    SampleSheetError,
    ValidationError,
    ViralConseqFileNotFoundError,
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
    "minimum_coverage": (1, None),
    "minimum_depth": (1, None),
    "minimum_length": (0, None),
    "minimum_read_length": (0, None),
    "af_threshold": (0.0, 1.0),
    "af_isnv_threshold": (0.0, 1.0),
}


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
        return validate_sample_sheet(sample_sheet, cast(str, data_type))
    if samples:
        return samples
    raise SampleConfigurationNotFoundError(
        "Either 'sample_sheet' or 'samples' must be provided. " f"Sample sheet path: {sample_sheet}"
    )
