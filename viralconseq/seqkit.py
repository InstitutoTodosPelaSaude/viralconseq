"""seqkit-backed fast detector for the FASTQ integrity check.

``integrity.scan_fastq`` is a pure-Python line scanner that names the exact
problem in a FASTQ file, but it costs about 1.3 s per 8 MB of gzip. seqkit reads
the same file about four times faster on one core, so it is used here as a *detector*
in front of the scanner: when seqkit says a file is clean the scanner is skipped,
and when seqkit flags anything the scanner runs to produce the precise issue
code. This module is the seqkit side only; wiring it into ``integrity`` is the
caller's job.

The detector command is ``seqkit seq --validate-seq --seq-type dna --max-len 0
--name --only-id --threads 1 -- <path>`` (:data:`SEQKIT_ARGS`). Measured against
seqkit v2.14.0:

* exit 0 and empty stdout for a clean file -- plain or gzip, multi-member gzip,
  a suffix that lies about compression either way, CRLF line ends, a UTF-8 BOM,
  blank lines between records;
* exit 255 and an ANSI-coloured ``[ERRO] ...`` line on stderr for a last record
  cut after the sequence or ``+`` line (``fastx: unequal sequence and
  quality``), a truncated gzip (``unexpected EOF``), gzip magic followed by
  garbage (``fastx: gzip: invalid header``), a mid-file length mismatch, bad
  header or bad separator (``fastx: bad fastq format``), a non-DNA letter
  (``seq: invalid DNAredundant letter: Z``), a directory, a permission error or a
  missing file;
* exit 0 but a read id on stdout for a file cut after a header line or inside a
  sequence line without a trailing newline: the cut yields a zero-length read,
  which ``--max-len 0`` prints. "exit 0 with non-empty stdout" therefore means
  "suspicious, run the Python scanner";
* exit 0 and empty stdout for an empty file (0 bytes, or a gzip of an empty
  stream), which is why :func:`probe_fastq` checks emptiness *before* seqkit;
  it also reads a FASTA given as FASTQ without complaint (``>`` records are
  FASTA to it), accepts a line-wrapped (multi-line) FASTQ that fastp turns
  into zero reads, and decompresses bzip2/xz/zstd transparently, so the probe
  reports the first two as ``suspect`` and the last as ``unreadable``.

Known blind spots -- seqkit does **not** catch these, and nothing here claims
to: quality characters outside ASCII 33-126 (a space and byte 0x80 both pass),
and whitespace or gap characters (``-``, ``.``, space, tab) inside sequence
lines (``--seq-type dna`` accepts them as part of the redundant alphabet). Only
the Python scanner sees those, so a ``clean`` verdict from this module is "no
structural damage", not "every byte is valid".

seqkit colours stderr even without a tty; :func:`strip_ansi` removes the escape
sequences before matching, and ``[WARN]``/``[INFO]`` lines (for example ``[WARN]
you may switch on flag -g/--remove-gaps to remove spaces``) are ignored.

The module is stdlib-only and importable without seqkit installed: only
:func:`find_seqkit` and :func:`detect` touch the binary.
"""

import gzip
import os
import re
import shutil
import subprocess
import zlib
from dataclasses import dataclass, field
from typing import Any, List, Optional

#: Arguments to ``seqkit`` that make up the detector command (the binary and the
#: input path are added around them by :func:`detect`).
SEQKIT_ARGS: List[str] = [
    "seq",
    "--validate-seq",
    "--seq-type",
    "dna",
    "--max-len",
    "0",
    "--name",
    "--only-id",
    "--threads",
    "1",
]

#: The four verdicts a :class:`SeqkitResult` can carry.
VERDICT_CLEAN = "clean"
VERDICT_EMPTY = "empty"
VERDICT_UNREADABLE = "unreadable"
VERDICT_SUSPECT = "suspect"

_GZIP_MAGIC = b"\x1f\x8b"
_BOM = b"\xef\xbb\xbf"
_WHITESPACE = b" \t\r\n\x0b\x0c"
# Compression the pipeline's tools (minimap2, fastp, Clair3) cannot read; seqkit
# would decompress them transparently and report the FASTQ as clean.
_OTHER_COMPRESSION = ((b"BZh", "bzip2"), (b"\xfd7zXZ\x00", "xz"), (b"\x28\xb5\x2f\xfd", "zstd"))
_PROBE_CHUNK = 4096
_PROBE_LIMIT = 1 << 20

_ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")

# seqkit prefixes its error lines with "[ERRO]" (sic). Matched after ANSI
# stripping, so the regex sees the bare prefix.
_ERRO_PREFIX = "[ERRO]"

# Fragments of seqkit error messages that mean the file could not be read at
# all, as opposed to being a readable but malformed FASTQ. Compared lowercase.
_UNREADABLE_FRAGMENTS = (
    "unexpected eof",
    "gzip: invalid header",
    "input is a directory",
    "permission denied",
    "no such file",
)


class SeqkitMissingError(RuntimeError):
    """Raised when the seqkit binary cannot be found or executed."""

    def __init__(self, message: Optional[str] = None) -> None:
        super().__init__(
            message
            or (
                "seqkit is required for the FASTQ integrity check: install it in the "
                "viralconseq environment or from bioconda, or pass "
                "--skip-input-validation."
            )
        )


@dataclass
class SeqkitResult:
    """Outcome of :func:`detect` (or of :func:`probe_fastq` alone).

    ``verdict`` is one of ``"clean"`` (seqkit saw no structural problem; the
    Python scanner may be skipped), ``"empty"`` (no record at all),
    ``"unreadable"`` (the file could not be opened or decompressed; the caller
    can report that without rescanning) or ``"suspect"`` (seqkit flagged
    something, or printed a zero-length read id; run the Python scanner to name
    the issue).

    ``returncode`` and ``command`` describe the seqkit invocation so the caller
    can log it. When the verdict came from :func:`probe_fastq`, seqkit was never
    run: ``returncode`` is 0 and ``command`` is empty. ``message`` is the first
    ANSI-stripped ``[ERRO]`` line from seqkit, the ``OSError`` text from the
    probe, or an empty string.
    """

    verdict: str
    returncode: int = 0
    message: str = ""
    command: List[str] = field(default_factory=list)

    @property
    def ran_seqkit(self) -> bool:
        """True when the verdict came from a seqkit run rather than the probe."""
        return bool(self.command)


def strip_ansi(text: str) -> str:
    """Remove ANSI SGR colour sequences (``ESC [ ... m``) from ``text``."""
    return _ANSI_RE.sub("", text)


def find_seqkit() -> str:
    """Return the absolute path of the ``seqkit`` binary on ``PATH``.

    Raises :class:`SeqkitMissingError` when it is not found.
    """
    found = shutil.which("seqkit")
    if not found:
        raise SeqkitMissingError()
    return os.path.abspath(found)


def _first_content(path: str) -> bytes:
    """The first bytes of ``path`` (decompressed when gzip) after any UTF-8 BOM
    and leading ASCII whitespace; ``b""`` when there are none within
    :data:`_PROBE_LIMIT` bytes. Raises ``OSError``/``EOFError``/``zlib.error``
    when the file or its gzip header cannot be read."""
    with open(path, "rb") as fh:
        head = fh.read(len(_GZIP_MAGIC))
        if head != _GZIP_MAGIC:
            fh.seek(0)
            return _skip_whitespace(fh)
    with gzip.open(path, "rb") as gz:
        return _skip_whitespace(gz)


def _skip_whitespace(handle: Any) -> bytes:
    seen = 0
    first = True
    while seen < _PROBE_LIMIT:
        chunk = handle.read(_PROBE_CHUNK)
        if not chunk:
            return b""
        if first and chunk.startswith(_BOM):
            chunk = chunk[len(_BOM) :]
        first = False
        content = chunk.lstrip(_WHITESPACE)
        if content:
            return content
        seen += len(chunk)
    return b"\0"  # whitespace beyond the limit: not empty, let seqkit decide


def probe_fastq(path: str) -> Optional[SeqkitResult]:
    """Cheap Python prelude run before seqkit.

    Returns an ``unreadable`` result when the path is missing, a directory,
    unreadable (permissions), compressed with something other than gzip
    (bzip2, xz, zstd: the pipeline's tools read plain or gzip FASTQ only) or
    starts with the gzip magic but has a corrupt gzip header; an ``empty``
    result when the file holds nothing but an optional UTF-8 BOM and
    whitespace (0-byte file, a gzip of an empty stream, an "Excel export" of
    nothing); a ``suspect`` result when the first content byte is ``>`` (a
    FASTA mislabelled as FASTQ, which seqkit would read happily and the
    Python scanner names) or when the first record's third line is not a
    ``+`` separator (a line-wrapped FASTQ: seqkit reads it, fastp emits no
    reads, the scanner reports ``fastq_bad_separator``); and ``None``
    otherwise, meaning "ask seqkit".
    """
    try:
        with open(path, "rb") as fh:
            head = fh.read(8)
    except OSError as exc:
        return SeqkitResult(verdict=VERDICT_UNREADABLE, message=str(exc))
    if not head:
        return SeqkitResult(verdict=VERDICT_EMPTY)
    for magic, name in _OTHER_COMPRESSION:
        if head.startswith(magic):
            return SeqkitResult(
                verdict=VERDICT_UNREADABLE,
                message=f"{name}-compressed input is not supported (plain or gzip FASTQ only)",
            )
    try:
        content = _first_content(path)
    except (OSError, EOFError, zlib.error) as exc:
        # gzip.BadGzipFile is an OSError; EOFError/zlib.error cover a stream
        # that ends or breaks inside the header. A gzip truncated later in the
        # stream usually yields its first bytes fine and is left to seqkit.
        return SeqkitResult(verdict=VERDICT_UNREADABLE, message=str(exc))
    if not content:
        return SeqkitResult(verdict=VERDICT_EMPTY)
    if content.startswith(b">"):
        return SeqkitResult(
            verdict=VERDICT_SUSPECT, message="file starts with '>' (a FASTA given as FASTQ?)"
        )
    if _first_record_is_wrapped(path):
        return SeqkitResult(
            verdict=VERDICT_SUSPECT,
            message="the first record's third line is not a '+' separator (a line-wrapped FASTQ?)",
        )
    return None


def _first_record_is_wrapped(path: str) -> bool:
    """True when the first record's third line does not start with ``+``: a
    multi-line (line-wrapped) FASTQ, which seqkit reads happily but fastp
    turns into zero reads without a word. The file was readable when this is
    called. Fewer than three complete lines within :data:`_PROBE_LIMIT` bytes
    (one giant read, or a file cut short) is left to seqkit."""
    try:
        with open(path, "rb") as fh:
            raw = fh.read(_PROBE_LIMIT)
        if raw.startswith(_GZIP_MAGIC):
            with gzip.open(path, "rb") as gz:
                raw = gz.read(_PROBE_LIMIT)
    except (OSError, EOFError, zlib.error):
        return False
    if raw.startswith(_BOM):
        raw = raw[len(_BOM) :]
    lines = raw.lstrip(_WHITESPACE).split(b"\n")
    if len(lines) < 4:  # the third line must be complete, so a fourth has begun
        return False
    return not lines[2].lstrip(b" \t\r").startswith(b"+")


def _first_error_line(stderr: str) -> str:
    """Return the first ``[ERRO]`` line of seqkit's stderr, ANSI-stripped."""
    for raw in stderr.splitlines():
        line = strip_ansi(raw).strip()
        if line.startswith(_ERRO_PREFIX):
            return line
    return ""


def _is_unreadable_message(message: str) -> bool:
    lowered = message.lower()
    return any(fragment in lowered for fragment in _UNREADABLE_FRAGMENTS)


def detect(path: str, seqkit: Optional[str] = None) -> SeqkitResult:
    """Run the probe and then, if needed, the seqkit detector on ``path``.

    ``seqkit`` is the binary to run; ``None`` means :func:`find_seqkit`. The
    verdict mapping is: probe result when it is not ``None``; exit 0 with empty
    stdout gives ``clean``; exit 0 with non-empty stdout (zero-length read ids)
    gives ``suspect``; a non-zero exit gives ``suspect``, except that an
    ``[ERRO]`` message naming an I/O failure (``unexpected EOF``, ``gzip:
    invalid header``, ``input is a directory``, ``permission denied``, ``no such
    file``) gives ``unreadable``.
    """
    probed = probe_fastq(path)
    if probed is not None:
        return probed

    binary = seqkit if seqkit is not None else find_seqkit()
    command = [binary, *SEQKIT_ARGS, "--", path]
    try:
        completed = subprocess.run(
            command,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            errors="replace",
            check=False,
        )
    except FileNotFoundError:
        raise SeqkitMissingError(
            f"seqkit binary not found at {binary!r}: install seqkit in the viralconseq "
            "environment or from bioconda, or pass --skip-input-validation."
        ) from None

    message = _first_error_line(completed.stderr or "")
    stdout = completed.stdout or ""

    if completed.returncode == 0:
        # seqkit prints nothing at all for a clean file; a bare newline is the
        # id of a zero-length read with an empty name (a file cut after "@").
        verdict = VERDICT_CLEAN if not stdout else VERDICT_SUSPECT
    elif _is_unreadable_message(message):
        verdict = VERDICT_UNREADABLE
    else:
        verdict = VERDICT_SUSPECT

    return SeqkitResult(
        verdict=verdict,
        returncode=completed.returncode,
        message=message,
        command=command,
    )
