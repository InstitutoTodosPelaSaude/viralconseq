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
