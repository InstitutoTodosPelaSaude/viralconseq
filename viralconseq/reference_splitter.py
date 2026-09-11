"""Split a single multi-record reference into per-segment files.

Segmented consensus runs historically required one ``--segmented-reference
SEGMENT=PATH`` per segment. To make the common case easier, a user can instead
pass a single multi-FASTA to ``--reference``; when it holds more than one
record the pipeline treats it as segmented. This module turns such a file into
the exact ``{segment: path}`` mapping the segmented workflows already consume,
so nothing under ``viralconseq/scripts/`` has to change.

Segment names are derived from the FASTA headers: the first whitespace/pipe
token (usually the accession), reduced to a filesystem/shell-safe allowlist
(``[A-Za-z0-9._-]``; anything else becomes ``_``). Each split file keeps its
record's original header verbatim, so the coverage-table contig id is
unaffected.
"""

import gzip
import logging
import os
import re
from typing import Dict, List, Set, Tuple

logger = logging.getLogger(__name__)

# Segment keys become directory names and Snakemake wildcard values, so keep
# them to a conservative filesystem/shell-safe allowlist. This is a superset of
# the nanopore header sanitisation (which maps / \ | , ~ and space to _): any
# character outside the allowlist is collapsed to _.
_SEGMENT_NAME_RE = re.compile(r"[^A-Za-z0-9._-]")


def _open_text(path: str):
    """Open a plain or gzipped text file for reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def sanitize_segment_name(header: str) -> str:
    """Derive a filesystem-safe segment key from a FASTA header.

    Strips a leading ``>``, takes the first whitespace/pipe-delimited token,
    then maps any character outside ``[A-Za-z0-9._-]`` to ``_`` so the key is a
    safe directory name and Snakemake wildcard value.

    Raises:
        ValueError: if the header yields an empty or unusable token
            (empty, ``.`` or ``..``).
    """
    token = header.lstrip(">").strip()
    token = re.split(r"[\s|]", token)[0] if token else ""
    key = _SEGMENT_NAME_RE.sub("_", token)
    if not key or key in (".", ".."):
        raise ValueError(f"Could not derive a valid segment name from header: {header!r}")
    return key


def parse_fasta(path: str) -> List[Tuple[str, List[str]]]:
    """Parse a (optionally gzipped) FASTA into ``(header, body_lines)`` records.

    ``header`` is the header line without the leading ``>`` and without the
    trailing newline. ``body_lines`` are the raw sequence lines (newlines kept)
    so each record can be re-emitted verbatim.
    """
    records: List[Tuple[str, List[str]]] = []
    current_header = None
    current_body: List[str] = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, current_body))
                current_header = line[1:].strip()
                current_body = []
            elif current_header is not None:
                current_body.append(line)
    if current_header is not None:
        records.append((current_header, current_body))
    return records


def count_records(path: str) -> int:
    """Count the number of records (``>`` headers) in a FASTA file."""
    count = 0
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count


def _dedup_key(key: str, seen: Set[str]) -> str:
    """Return ``key`` (or ``key_2``, ``key_3`` ...) so every key is unique.

    Both the original and any generated suffix key are registered in ``seen``,
    so a header that already sanitises to ``foo_2`` cannot silently collide with
    the de-duplicated form of a duplicate ``foo``.
    """
    if key not in seen:
        seen.add(key)
        return key
    i = 2
    new_key = f"{key}_{i}"
    while new_key in seen:
        i += 1
        new_key = f"{key}_{i}"
    seen.add(new_key)
    logger.warning(
        "Duplicate segment name %r derived from reference headers; using %r instead.",
        key,
        new_key,
    )
    return new_key


def split_multifasta(path: str, out_dir: str) -> Dict[str, str]:
    """Split a multi-record FASTA into one single-record file per segment.

    Each record's original header is preserved verbatim inside its file; the
    segment key (dict key, and later the wildcard/directory name) is the
    sanitised header token. Returns ``{segment_key: absolute_path}`` in the
    original record order.

    Raises:
        ValueError: if ``path`` contains no records.
    """
    records = parse_fasta(path)
    if not records:
        raise ValueError(f"Reference FASTA contains no sequences: {path}")

    os.makedirs(out_dir, exist_ok=True)
    mapping: Dict[str, str] = {}
    seen: Set[str] = set()
    for header, body in records:
        if not "".join(body).strip():
            raise ValueError(
                f"Reference record {header!r} in {path} has no sequence; "
                "cannot use it as a segment."
            )
        key = _dedup_key(sanitize_segment_name(header), seen)
        seg_path = os.path.abspath(os.path.join(out_dir, f"{key}.fasta"))
        text = ">" + header + "\n" + "".join(body)
        if not text.endswith("\n"):
            text += "\n"
        with open(seg_path, "w") as fh:
            fh.write(text)
        mapping[key] = seg_path
    return mapping
