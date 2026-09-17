"""Custom exceptions for viralconseq pipeline.

Each exception carries a stable, machine-readable ``code`` so an embedding
service can surface a structured error payload (``to_dict()``) to callers
instead of only a log line.
"""

from typing import Optional


class ViralConseqError(Exception):
    """Base exception for all viralconseq errors."""

    #: Stable, machine-readable identifier for this error class. Overridden by
    #: subclasses; may also be overridden per-instance via the ``code`` kwarg.
    code = "viralconseq_error"

    def __init__(self, message: str = "", *, code: Optional[str] = None):
        super().__init__(message)
        self.message = message
        if code is not None:
            self.code = code

    def to_dict(self) -> dict:
        """Return a structured, JSON-serializable representation of the error."""
        return {
            "error": type(self).__name__,
            "code": self.code,
            "message": str(self),
        }


class ValidationError(ViralConseqError):
    """Raised when validation of arguments or data fails."""

    code = "validation_error"


class InputIntegrityError(ValidationError):
    """Raised when a content-level input-integrity check fails.

    Carries the list of blocking :class:`~viralconseq.integrity.IntegrityIssue`
    objects that caused the failure so an embedding service can surface a
    structured, per-file payload rather than only the aggregated message. It
    subclasses :class:`ValidationError` so existing ``except ViralConseqError``
    handling (clean exit, no traceback) applies unchanged.
    """

    code = "input_integrity_error"

    def __init__(self, message: str = "", *, issues=None, code: Optional[str] = None):
        super().__init__(message, code=code)
        self.issues = list(issues) if issues else []

    def to_dict(self) -> dict:
        payload = super().to_dict()
        payload["issues"] = [
            issue.as_dict() if hasattr(issue, "as_dict") else issue for issue in self.issues
        ]
        return payload


class ViralConseqFileNotFoundError(ViralConseqError):
    """Raised when a required file or directory is not found."""

    code = "file_not_found"


class ConfigurationError(ViralConseqError):
    """Raised when there's an error in configuration."""

    code = "configuration_error"


class SampleSheetError(ViralConseqError):
    """Raised when there's an error processing the sample sheet."""

    code = "sample_sheet_error"


class SampleConfigurationNotFoundError(ViralConseqError):
    """Raised when the sample configuration is not found."""

    code = "sample_configuration_not_found"


class ReferenceNotFoundError(ViralConseqError):
    """Raised when the reference sequence file is not found."""

    code = "reference_not_found"


class PrimerSchemeNotFoundError(ViralConseqError):
    """Raised when the primer scheme file is not found."""

    code = "primer_scheme_not_found"


class AdaptersNotFoundError(ViralConseqError):
    """Raised when the Illumina adapter sequences file is not found or not provided."""

    code = "adapters_not_found"


class Clair3ModelError(ViralConseqError):
    """Base for Clair3 model resolution problems (nanopore runs)."""

    code = "clair3_model_error"


class Clair3ModelUnresolvedError(Clair3ModelError):
    """No Clair3 model could be chosen from the reads' basecall tags.

    Raised when the reads carry no ``basecall_model_version_id``/``RG:Z`` tag
    (Guppy-era or re-headered reads), when the tag matches no manifest model,
    or when only Guppy-era models match a Dorado tag. The message says which
    ``--clair3-model NAME`` to pass instead.
    """

    code = "clair3_model_unresolved"


class Clair3ModelMixedError(Clair3ModelError):
    """The reads of one sample resolve to different Clair3 models."""

    code = "clair3_model_mixed_within_sample"


class Clair3ModelNotFoundError(Clair3ModelError):
    """The chosen model is not in the model directory (or is incomplete).

    Raised at validation time, before any work starts; the message names the
    directory and the ``viralconseq setup --clair3-models`` command.
    """

    code = "clair3_model_not_found"


class ReportError(ViralConseqError):
    """Raised by ``viralconseq create-report`` when the directory is not a
    finished run or the page cannot be built (for example when it fails its
    self-check against ``summary.tsv``)."""

    code = "report_error"


class ViralQCDatabaseNotFoundError(ViralConseqError):
    """Raised when viralQC is enabled but its database directory is missing or
    incomplete. The message names the directory and the ``viralconseq setup``
    command that populates it."""

    code = "viralqc_database_not_found"
