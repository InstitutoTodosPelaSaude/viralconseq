"""Content-level integrity validation for consensus-pipeline input files.

viralconseq's existing validators check only that inputs *exist*. When viralconseq
is the analysis engine behind a web service, existence is not enough: a truncated
upload, a corrupt gzip, a protein FASTA, or a primer BED whose chrom names don't
match the reference all pass an existence check and then fail deep inside
Snakemake -- or, worse, silently produce wrong results (``samtools ampliconclip``
clips nothing when no BED chrom matches a mapped contig). This module adds
streaming, dependency-free validators for the three consensus input formats
(FASTQ, FASTA, BED) so bad inputs are rejected up front with clear, structured
errors.

Each validator returns an :class:`IntegrityReport` that *collects*
:class:`IntegrityIssue` objects rather than raising, so the orchestration layer
(``validators.validate_consensus_input_integrity``) can aggregate problems across
every input and decide which severities block a run. The formats here are simple
enough that hand-rolled streaming parsers are both easy and give full control over
error codes and messages.
"""

import gzip
import re
import zlib
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Set

# Cap the number of issues recorded per file. A file with more than this many
# problems is clearly broken; reporting the first N keeps error payloads bounded
# while still letting a user fix several issues in one round-trip.
MAX_ISSUES_PER_FILE = 20

_GZIP_MAGIC = b"\x1f\x8b"

# IUPAC nucleotide alphabet (including U for RNA and the ambiguity codes), plus N.
# Checked case-insensitively. Gaps ('-', '.') are intentionally excluded: a
# consensus reference must not contain gaps.
_IUPAC_NT: Set[str] = set("ACGTURYSWKMBDHVN")

# Mirrors the nanopore ``rule sanitize_reference`` sed
# (``sed '/^>/s/[\\/|,~ ]/_/g'`` in consensus_nanopore*.smk): backslash, '/', '|',
# ',', '~' and space are replaced with '_' across the whole header. Because spaces
# are removed, the sanitized header is a single token, so it becomes the BAM @SQ
# contig name in the nanopore pipeline.
_NANOPORE_SANITIZE_RE = re.compile(r"[/\\|,~ ]")

# A FASTA/mapper contig id is the first WHITESPACE-delimited token of the header.
# This must match how minimap2/samtools name the BAM @SQ contig and how the
# consensus workflows list contigs (`grep '^>' | sed 's/^>//' | cut -d' ' -f1`):
# both split on whitespace only and keep a literal '|'. (reference_splitter's
# sanitize_segment_name additionally splits on '|' to derive internal segment
# *directory* keys -- a different concern from the BAM contig identity used for
# primer-BED chrom matching, so it is intentionally not reused here.)


@dataclass
class IntegrityIssue:
    """A single problem found in an input file.

    ``severity`` is ``"error"`` (blocking) or ``"warning"`` (non-blocking).
    ``code`` is a stable, machine-readable slug for the service to key off.
    """

    severity: str
    code: str
    message: str
    line: Optional[int] = None
    path: Optional[str] = None
    kind: Optional[str] = None

    def as_dict(self) -> dict:
        """JSON-serializable representation for structured error payloads."""
        payload: Dict[str, Any] = {
            "severity": self.severity,
            "code": self.code,
            "message": self.message,
        }
        if self.line is not None:
            payload["line"] = self.line
        if self.path is not None:
            payload["path"] = self.path
        if self.kind is not None:
            payload["kind"] = self.kind
        return payload


@dataclass
class IntegrityReport:
    """The collected issues for one validated file.

    ``contig_ids``/``headers`` are populated by :func:`validate_fasta` only and
    are used by the orchestration layer for cross-file (BED) name matching.
    """

    path: str
    kind: str
    issues: List[IntegrityIssue] = field(default_factory=list)
    contig_ids: List[str] = field(default_factory=list)
    headers: List[str] = field(default_factory=list)

    @property
    def errors(self) -> List[IntegrityIssue]:
        return [i for i in self.issues if i.severity == "error"]

    @property
    def warnings(self) -> List[IntegrityIssue]:
        return [i for i in self.issues if i.severity == "warning"]


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _is_gzip(path: str) -> bool:
    """Return True if *path* begins with the gzip magic bytes.

    We sniff the content rather than trust the ``.gz`` suffix: a mislabeled file
    (a plain FASTQ named ``.gz``, or vice versa) is exactly the kind of breakage
    these checks exist to catch.
    """
    with open(path, "rb") as handle:
        return handle.read(2) == _GZIP_MAGIC


def open_maybe_gzip(path: str):
    """Open *path* as text, transparently decompressing if it is gzipped.

    Uses ``utf-8-sig`` so a leading UTF-8 BOM (common in files exported from
    Windows/Excel -- exactly the untrusted-upload case) is stripped rather than
    glued to the first line, which would otherwise make a valid FASTA/FASTQ look
    malformed. Invalid non-UTF-8 bytes surface as ``UnicodeDecodeError``, which
    the validators catch and report as an unreadable/corrupt file.
    """
    if _is_gzip(path):
        return gzip.open(path, "rt", encoding="utf-8-sig")
    return open(path, "rt", encoding="utf-8-sig")


def sample_fastq_headers(path: str, n: int = 20) -> List[str]:
    """Return the first ``n`` FASTQ header lines (without the leading ``@``).

    Used to read the basecall model tag Dorado/MinKNOW write into every read
    header. Tolerates blank lines like :func:`validate_fastq`; on a read error
    (corrupt gzip, truncated file) returns what was read so far, so the caller
    decides how to report an empty result.
    """
    headers: List[str] = []
    try:
        with open_maybe_gzip(path) as handle:
            record: List[str] = []
            for raw in handle:
                line = _strip_eol(raw)
                if not record and not line.strip():
                    continue
                record.append(line)
                if len(record) == 4:
                    if record[0].startswith("@"):
                        headers.append(record[0][1:])
                    record = []
                    if len(headers) >= n:
                        break
    except _READ_ERRORS:
        pass
    return headers


def header_token(header: str) -> str:
    """Return the contig id: first whitespace-delimited token of a FASTA header."""
    parts = header.lstrip(">").strip().split()
    return parts[0] if parts else ""


def sanitize_nanopore_contig(header: str) -> str:
    """Return the BAM contig name the nanopore pipeline derives from *header*.

    The nanopore workflow sanitizes reference headers (``/ \\ | , ~`` and spaces
    become ``_``) before mapping, so the mapped contig name is the whole
    sanitized header rather than just the accession.
    """
    return _NANOPORE_SANITIZE_RE.sub("_", header.lstrip(">").strip())


def _strip_eol(line: str) -> str:
    """Strip a trailing newline (handling CRLF) without touching other chars."""
    return line.rstrip("\r\n")


def _is_directive_line(line: str) -> bool:
    """Return True for a BED comment or UCSC ``track``/``browser`` directive.

    ``track``/``browser`` must be a standalone leading word, so a real data row
    whose chrom merely starts with those letters (e.g. a contig ``trackpox``) is
    not silently dropped.
    """
    if line.startswith("#"):
        return True
    first = line.split(None, 1)
    return bool(first) and first[0] in ("track", "browser")


# Type alias for the per-report issue-appending closure.
_Adder = Callable[..., None]


def _make_adder(report: IntegrityReport, default_severity: str = "error") -> _Adder:
    """Return a closure that appends issues to *report* up to the per-file cap.

    Each issue is stamped with the report's ``path`` and ``kind`` so a service
    consuming the structured payload can map every issue back to its file.
    """

    def add(
        code: str, message: str, line: Optional[int] = None, severity: Optional[str] = None
    ) -> None:
        if len(report.issues) >= MAX_ISSUES_PER_FILE:
            return
        report.issues.append(
            IntegrityIssue(
                severity or default_severity,
                code,
                message,
                line,
                path=report.path,
                kind=report.kind,
            )
        )

    return add


# Errors raised while reading/decoding a possibly-corrupt or mislabeled file.
# UnicodeDecodeError covers a binary/garbage upload opened as text.
_READ_ERRORS = (OSError, EOFError, gzip.BadGzipFile, zlib.error, UnicodeDecodeError)


# ---------------------------------------------------------------------------
# FASTQ
# ---------------------------------------------------------------------------


def _check_fastq_record(record: List[str], start_line: int, add: _Adder) -> None:
    """Validate one 4-line FASTQ record beginning at ``start_line``.

    All issues are anchored to the record's header line: blank lines between
    records are tolerated (see :func:`validate_fastq`), so the four physical
    lines are not necessarily contiguous and per-line offsets would be
    misleading.
    """
    header, seq, plus, qual = record
    if not header.startswith("@"):
        add(
            "fastq_bad_header",
            f"Record at line {start_line}: header does not start with '@'.",
            start_line,
        )
    if not plus.startswith("+"):
        add(
            "fastq_bad_separator",
            f"Record at line {start_line}: separator line does not start with '+'.",
            start_line,
        )
    if len(seq) != len(qual):
        add(
            "fastq_length_mismatch",
            f"Record at line {start_line}: sequence length ({len(seq)}) != quality length ({len(qual)}).",
            start_line,
        )
    bad_seq = sorted(set(seq.upper()) - _IUPAC_NT)
    if bad_seq:
        add(
            "fastq_bad_sequence_char",
            f"Record at line {start_line}: invalid sequence character(s) {bad_seq}.",
            start_line,
        )
    if any(not (33 <= ord(c) <= 126) for c in qual):
        add(
            "fastq_bad_quality_char",
            f"Record at line {start_line}: quality contains out-of-range character(s) (expected ASCII 33-126).",
            start_line,
        )


def validate_fastq(path: str) -> IntegrityReport:
    """Full streaming integrity scan of a (optionally gzipped) FASTQ file.

    Reads the whole file once, verifying 4-line record structure, the ``@``/``+``
    marker lines, matching sequence/quality lengths, and valid sequence and
    quality characters. Blank lines (a common, benign artifact) are tolerated so
    they do not desynchronize record grouping. A trailing partial record is
    reported as truncation; a read/decode error (truncated or corrupt gzip,
    non-text bytes) is reported as unreadable.
    """
    report = IntegrityReport(path=path, kind="fastq")
    add = _make_adder(report)

    content_seen = False
    record: List[str] = []
    record_start = 0
    try:
        with open_maybe_gzip(path) as handle:
            # Scan to EOF even past the issue cap: truncation is only knowable at
            # the end, and a full scan is the same cost as the happy path. add()
            # still bounds how many issues are recorded.
            for line_no, raw in enumerate(handle, 1):
                line = _strip_eol(raw)
                if not line.strip():
                    continue  # tolerate blank lines between/after records
                content_seen = True
                if not record:
                    record_start = line_no
                record.append(line)
                if len(record) == 4:
                    _check_fastq_record(record, record_start, add)
                    record = []
    except _READ_ERRORS as exc:
        add(
            "fastq_unreadable",
            f"Could not read FASTQ (possibly a truncated or corrupt gzip): {exc}",
        )
        return report

    if not content_seen:
        add("fastq_empty", "FASTQ file is empty.")
        return report
    # A leftover partial record at EOF means truncation -- a critical signal, so
    # record it even if the issue cap was already reached.
    if record:
        report.issues.append(
            IntegrityIssue(
                "error",
                "fastq_truncated",
                f"File ends with an incomplete FASTQ record ({len(record)} of 4 lines); "
                "the file may be truncated.",
                record_start,
                path=report.path,
                kind=report.kind,
            )
        )
    return report


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------


def validate_fasta(path: str, require_nucleotide: bool = True) -> IntegrityReport:
    """Streaming integrity scan of a (optionally gzipped) reference FASTA.

    Verifies at least one record, no content before the first ``>`` header, that
    every record has a header and a non-empty sequence, that contig ids (first
    header token) are unique, and -- when ``require_nucleotide`` -- that
    sequences use only the IUPAC nucleotide alphabet. Populates ``contig_ids``
    and ``headers`` for downstream cross-file matching.
    """
    report = IntegrityReport(path=path, kind="fasta")
    add = _make_adder(report)

    seen_ids: Set[str] = set()
    current_header: Optional[str] = None
    current_header_line = 0
    current_has_seq = False
    n_records = 0
    saw_junk = False

    def _close_record() -> None:
        if current_header is not None and not current_has_seq:
            add(
                "fasta_empty_record",
                f"Record '{current_header}' (line {current_header_line}) has no sequence.",
                current_header_line,
            )

    try:
        with open_maybe_gzip(path) as handle:
            for line_no, raw in enumerate(handle, 1):
                line = _strip_eol(raw)
                if line.startswith(">"):
                    _close_record()
                    n_records += 1
                    current_header = line[1:].strip()
                    current_header_line = line_no
                    current_has_seq = False
                    token = header_token(line)
                    if not token:
                        add(
                            "fasta_bad_header",
                            f"Line {line_no}: FASTA header has no usable id.",
                            line_no,
                        )
                    elif token in seen_ids:
                        add(
                            "fasta_duplicate_id",
                            f"Duplicate contig id '{token}' (line {line_no}); "
                            "minimap2/samtools require unique contig names.",
                            line_no,
                        )
                    else:
                        seen_ids.add(token)
                        report.contig_ids.append(token)
                        report.headers.append(current_header)
                elif current_header is None:
                    if line.strip() and not saw_junk:
                        saw_junk = True
                        add(
                            "fasta_junk_before_header",
                            f"Line {line_no}: content before the first '>' header; not a valid FASTA.",
                            line_no,
                        )
                else:
                    seq = line.strip()
                    if seq:
                        current_has_seq = True
                        if require_nucleotide:
                            bad = sorted(set(seq.upper()) - _IUPAC_NT)
                            if bad:
                                add(
                                    "fasta_non_nucleotide",
                                    f"Record '{current_header}' (line {line_no}): non-nucleotide "
                                    f"character(s) {bad}; consensus requires a nucleotide reference.",
                                    line_no,
                                )
                if len(report.issues) >= MAX_ISSUES_PER_FILE:
                    break
            else:
                _close_record()
    except _READ_ERRORS as exc:
        add("fasta_unreadable", f"Could not read FASTA (possibly a corrupt gzip): {exc}")
        return report

    if n_records == 0:
        add("fasta_no_records", "FASTA file contains no sequences.")
    return report


# ---------------------------------------------------------------------------
# BED (primer scheme)
# ---------------------------------------------------------------------------


def _format_contig_hint(expected_contigs: Set[str], limit: int = 5) -> str:
    """Render a short, deterministic sample of expected contig names."""
    names = sorted(expected_contigs)
    shown = ", ".join(names[:limit])
    if len(names) > limit:
        shown += f", ... (+{len(names) - limit} more)"
    return shown


def validate_bed(
    path: str,
    expected_contigs: Optional[Set[str]] = None,
    accession_map: Optional[dict] = None,
) -> IntegrityReport:
    """Streaming integrity scan of a primer-scheme BED file.

    Verifies at least one interval, >=3 tab columns, integer ``start``/``end``
    with ``0 <= start < end``, and (when ``expected_contigs`` is given) that
    chrom names match the reference contigs. Chrom matching is partial-aware: a
    chrom matching no reference contig blocks the run only when NO chrom in the
    file matches (a wrong scheme -- ``samtools ampliconclip`` would clip
    nothing); if some chroms match, unmatched rows are downgraded to warnings (a
    whole-scheme BED reused on a subset of segments). ``accession_map`` (nanopore
    only: raw accession -> sanitized full header) always yields a blocking fix-it
    error when a chrom matches the accession but not the sanitized header the
    pipeline will use. An out-of-range strand column is a warning.
    """
    report = IntegrityReport(path=path, kind="bed")
    add = _make_adder(report)
    n_features = 0
    matched_any = False
    # Deferred plain-mismatch rows: their severity depends on whether ANY chrom
    # matched (see below), which is only known after the whole file is read.
    unmatched: List[tuple] = []

    try:
        with open_maybe_gzip(path) as handle:
            for line_no, raw in enumerate(handle, 1):
                line = _strip_eol(raw)
                if not line.strip() or _is_directive_line(line):
                    continue
                cols = line.split("\t")
                if len(cols) < 3:
                    add(
                        "bed_too_few_columns",
                        f"Line {line_no}: BED requires at least 3 tab-separated columns, found {len(cols)}.",
                        line_no,
                    )
                    continue

                n_features += 1
                chrom, start_s, end_s = cols[0], cols[1], cols[2]

                try:
                    start, end = int(start_s), int(end_s)
                    if start < 0 or end < 0:
                        add(
                            "bed_negative_coord",
                            f"Line {line_no}: coordinates must be non-negative.",
                            line_no,
                        )
                    elif start >= end:
                        add(
                            "bed_bad_interval",
                            f"Line {line_no}: start ({start}) must be < end ({end}).",
                            line_no,
                        )
                except ValueError:
                    add(
                        "bed_non_integer_coord",
                        f"Line {line_no}: start/end must be integers, found {start_s!r}/{end_s!r}.",
                        line_no,
                    )

                if len(cols) >= 6 and cols[5] not in ("+", "-", "."):
                    add(
                        "bed_bad_strand",
                        f"Line {line_no}: strand column is {cols[5]!r}, expected '+', '-' or '.'.",
                        line_no,
                        severity="warning",
                    )

                if expected_contigs:
                    if chrom in expected_contigs:
                        matched_any = True
                    elif accession_map and chrom in accession_map:
                        # Nanopore accession-vs-sanitized: always a blocking error
                        # (ampliconclip would clip nothing on that contig).
                        add(
                            "bed_chrom_sanitized_mismatch",
                            f"Line {line_no}: primer BED chrom '{chrom}' matches the reference accession but "
                            f"nanopore sanitizes the header to '{accession_map[chrom]}'; samtools ampliconclip "
                            f"would clip nothing. Simplify the reference header to the bare accession, or set "
                            f"the BED chrom to '{accession_map[chrom]}'.",
                            line_no,
                        )
                    else:
                        unmatched.append((line_no, chrom))
    except _READ_ERRORS as exc:
        add("bed_unreadable", f"Could not read BED (possibly a corrupt gzip): {exc}")
        return report

    if n_features == 0:
        add("bed_no_features", "Primer BED file contains no intervals.")

    # A chrom that matches NO reference contig means those primers are silently
    # ignored by ampliconclip. If NONE of the BED's chroms match, the scheme is
    # wrong for this reference -> block. If some matched (e.g. a whole-scheme BED
    # reused on a subset of segments), the unmatched rows are harmless extra
    # primers -> warn.
    hint = _format_contig_hint(expected_contigs) if expected_contigs else ""
    for line_no, chrom in unmatched:
        if matched_any:
            add(
                "bed_chrom_unmatched",
                f"Line {line_no}: primer BED chrom '{chrom}' does not match any reference contig "
                f"({hint}); its primers will not be trimmed (harmless if not part of this run).",
                line_no,
                severity="warning",
            )
        else:
            add(
                "bed_chrom_mismatch",
                f"Line {line_no}: primer BED chrom '{chrom}' does not match any reference contig "
                f"({hint}); primer trimming would silently do nothing.",
                line_no,
            )

    return report
