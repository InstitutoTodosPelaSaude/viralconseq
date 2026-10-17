#!/usr/bin/env python
"""Cap the read depth of a coordinate-sorted BAM before consensus calling.

Run by ``rule normalise_depth`` (``rules/normalise.smk``) inside
``envs/alignment.yaml`` (samtools + python); standard library only, no
``viralconseq`` import. Deep amplicon runs (dengue at 30,000x, SARS-CoV-2 at
3,000x) hand every read to ``samtools consensus`` or Clair3 for no gain in
accuracy and a large cost in time and memory, while Clair3 rescales its
pileup tensor above 216x. This script keeps, for every reference position,
at least ``min(target, original depth)`` records and drops the rest, chosen
so that the strand ratio (and, with a primer scheme, the amplicon mix) of
the original data is preserved at every position.

Algorithm (one left-to-right sweep per contig):

* A **unit** is one record, or both mates of a proper pair (flags 0x1 and
  0x2, mate mapped, same contig), matched **by QNAME** -- ``samtools
  ampliconclip --hard-clip`` moves POS without rewriting the mate's PNEXT, so
  PNEXT is unreliable. A unit is kept or dropped whole. Secondary records
  (0x100) and unmapped records are dropped and counted; supplementary records
  (0x800) are single-record units on their own strand.
* A record covers the reference positions its ``M``/``D``/``=``/``X`` CIGAR
  operations span (``N`` consumes but does not cover, like ``samtools depth
  -J``); ``I``/``S``/``H``/``P`` consume nothing.
* **Strata.** Without a scheme the strata are the two strands. With a scheme
  they are (strand x amplicon): a record is assigned to the amplicon whose
  outer span contains it (10 bp tolerance; several -> the nearest primer
  ends, tie -> smaller amplicon number). A record that fits no single
  amplicon and contains a whole primer interval is a skip-amplicon product
  (``interior_primer``: it runs from one amplicon into the next, across the
  primers between them); one that fits no amplicon and contains no primer is
  ``unassigned``. Containing a neighbour's primer while fitting one amplicon
  is normal in a tiling scheme and is not a signature of anything.
  Both classes form one extra ``other`` stratum per strand, capped
  proportionally like any amplicon (so their share of the kept depth equals
  their share of the original depth) and counted in the reports. They used
  to be kept unconditionally; in LoCost libraries they cluster at primer
  edges carrying artefacts, and once the normal reads were capped they
  dominated the pileup and Clair3 called the artefact (benchmark v4).
* At each position ``p`` the target ``min(T, original depth)`` is split
  among the strata in the original proportions (largest-remainder rounding,
  never above a stratum's own depth). Units starting at ``p`` are pushed onto their
  stratum's pending set, ordered by furthest reach first (start plus aligned
  query bases, so a deletion or insertion does not change a read's rank) and
  then input order; while a stratum is below its target the furthest-reaching
  pending unit is popped and kept
  (both mates); a deficit a stratum cannot fill is filled from the other
  strata, so the total lower bound always holds.

Invariants, self-checked from the kept flags before any output is written:
final depth >= ``min(T, original)`` at every position (so every read is kept
wherever the original depth is below the target and ``--minimum-coverage``
masks are unchanged there), and per-strand depth <= ``2T`` for single-record
units (exactly ``T`` for staggered fixed-length reads; amplicon overlaps may
reach ``2T``, as in ARTIC's normalisation).

The target counts **records**, as ``samtools consensus -d`` and ``samtools
depth`` do: an overlapping Illumina pair contributes two. Kept records are
re-emitted in input order, so the output is coordinate-sorted and two runs
over the same input are byte-identical (the ``@PG`` line carries no paths or
timestamps). Input that is not in coordinate order exits with status
:data:`EXIT_UNSORTED` so the caller can sort and retry.

``--target 0`` switches normalisation off: the input BAM and index are
copied and header-only tables are written, without opening the BAM.
"""

import argparse
import array
import bisect
import heapq
import os
import re
import shutil
import subprocess
import sys
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, TextIO, Tuple

EXIT_UNSORTED = 3
EXIT_SELF_CHECK = 4

FLAG_PAIRED = 0x1
FLAG_PROPER = 0x2
FLAG_UNMAPPED = 0x4
FLAG_MATE_UNMAPPED = 0x8
FLAG_REVERSE = 0x10
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
# Same pattern as build_report.PRIMER_NAME (both stdlib scripts; kept in step
# by test/rule_inventory_test.py).
PRIMER_NAME = re.compile(r"^(?P<prefix>.+)_(?P<n>\d+)_(?P<side>LEFT|RIGHT)(?:_.*)?$")
#: Slack, in bp, allowed when testing whether an amplicon's outer span contains a record.
AMPLICON_TOLERANCE = 10

CLASS_NORMAL = 0
CLASS_INTERIOR_PRIMER = 1
CLASS_UNASSIGNED = 2
CLASS_NAMES = {
    CLASS_NORMAL: "normal",
    CLASS_INTERIOR_PRIMER: "interior_primer",
    CLASS_UNASSIGNED: "unassigned",
}

MODE_OFF = "off"
MODE_POSITION = "position"
MODE_AMPLICON = "amplicon"

DEPTH_COLUMNS = [
    "contig",
    "pos",
    "before",
    "after",
    "before_fwd",
    "before_rev",
    "after_fwd",
    "after_rev",
]
AMPLICON_COLUMNS = [
    "contig",
    "amplicon",
    "pool",
    "start",
    "end",
    "mean_before",
    "mean_after",
    "reads_before",
    "reads_after",
    "interior_primer_before",
    "interior_primer_after",
]
SUMMARY_KEYS = [
    "mode",
    "target",
    "records_in",
    "records_kept",
    "records_dropped",
    "secondary_dropped",
    "unmapped_dropped",
    "orphan_mates",
    "unassigned_before",
    "unassigned_after",
    "interior_primer_before",
    "interior_primer_after",
    "mean_span_before",
    "mean_span_after",
    "positions_capped",
    "pct_positions_capped",
    "self_check",
]
NA = "NA"


class NormaliseError(Exception):
    """A problem with the input that is not a bug: message for the log."""


class UnsortedInputError(NormaliseError):
    """The records are not in coordinate order."""


class SelfCheckError(Exception):
    """The kept set violates an invariant: a bug, never written out."""


# ---------------------------------------------------------------------------
# SAM parsing
# ---------------------------------------------------------------------------


def cigar_intervals(start: int, cigar: str) -> Tuple[List[Tuple[int, int]], int]:
    """Half-open reference intervals ``[a, b)`` a record covers, from its
    0-based ``start`` and CIGAR, and the number of aligned query bases
    (``M``/``I``/``=``/``X``). ``M``/``D``/``=``/``X`` cover, ``N`` splits,
    the rest consume no reference. ``*`` gives no interval."""
    intervals: List[Tuple[int, int]] = []
    pos = start
    run_start = start
    covering = False
    query = 0
    for length, op in CIGAR_RE.findall(cigar):
        n = int(length)
        if op in "MDX=":
            if not covering:
                run_start = pos
                covering = True
            pos += n
            if op != "D":
                query += n
        elif op == "I":
            query += n
        elif op == "N":
            if covering:
                intervals.append((run_start, pos))
                covering = False
            pos += n
    if covering:
        intervals.append((run_start, pos))
    return intervals, query


class Amplicon:
    """One amplicon of a tiling scheme: outer span (primer start to primer
    end) and insert (between the primers), with alternate primers merged."""

    __slots__ = (
        "contig",
        "name",
        "n",
        "pool",
        "outer_start",
        "outer_end",
        "insert_start",
        "insert_end",
    )

    def __init__(self, contig: str, name: str, n: int, pool: str) -> None:
        self.contig = contig
        self.name = name
        self.n = n
        self.pool = pool
        self.outer_start = -1
        self.outer_end = -1
        self.insert_start = -1
        self.insert_end = -1

    def add_left(self, start: int, end: int) -> None:
        self.outer_start = start if self.outer_start < 0 else min(self.outer_start, start)
        self.insert_start = max(self.insert_start, end)

    def add_right(self, start: int, end: int) -> None:
        self.outer_end = max(self.outer_end, end)
        self.insert_end = start if self.insert_end < 0 else min(self.insert_end, start)

    @property
    def complete(self) -> bool:
        return self.outer_start >= 0 and self.outer_end >= 0 and self.insert_start < self.insert_end


class Scheme:
    """Amplicons and primer intervals of a scheme BED, per contig."""

    def __init__(self) -> None:
        self.amplicons: Dict[str, List[Amplicon]] = {}
        self.primers: Dict[str, List[Tuple[int, int]]] = {}
        self.notes: List[str] = []

    def contigs(self) -> List[str]:
        return sorted(self.amplicons)


def _read_bed(path: str) -> List[Tuple[str, int, int, str, str, str]]:
    rows = []
    with open(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise NormaliseError(f"scheme BED line has fewer than 4 columns: {line!r}")
            pool = fields[4] if len(fields) > 4 and fields[4].strip() else NA
            strand = fields[5].strip() if len(fields) > 5 else ""
            rows.append((fields[0], int(fields[1]), int(fields[2]), fields[3], pool, strand))
    return rows


def parse_scheme(path: str) -> Scheme:
    """Parse a primer BED into amplicons. Primer names are paired by the
    ``<prefix>_<n>_(LEFT|RIGHT)[_alt]`` convention; when a name does not parse,
    the whole scheme falls back to pairing the i-th ``+`` primer with the i-th
    ``-`` primer along each contig (strand column; alternate primers cannot be
    recognised then), and a note says so. Neither working is an error."""
    scheme = Scheme()
    rows = _read_bed(path)
    if not rows:
        raise NormaliseError(f"scheme BED {path} holds no primer")
    parsed = [PRIMER_NAME.match(name) for _, _, _, name, _, _ in rows]
    by_key: Dict[Tuple[str, str, int], Amplicon] = {}
    if all(parsed):
        for (contig, start, end, name, pool, _strand), match in zip(rows, parsed):
            assert match is not None
            n = int(match.group("n"))
            key = (contig, match.group("prefix"), n)
            amp = by_key.get(key)
            if amp is None:
                amp = by_key[key] = Amplicon(contig, f"{match.group('prefix')}_{n}", n, pool)
            if match.group("side") == "LEFT":
                amp.add_left(start, end)
            else:
                amp.add_right(start, end)
    else:
        bad = [name for (_, _, _, name, _, _), m in zip(rows, parsed) if m is None]
        scheme.notes.append(
            f"{len(bad)} primer name(s) do not follow <prefix>_<n>_LEFT/RIGHT (e.g. {bad[0]!r}); "
            "pairing consecutive LEFT/RIGHT primers by the strand column instead"
        )
        if any(strand not in ("+", "-") for _, _, _, _, _, strand in rows):
            raise NormaliseError(
                "scheme BED primer names do not parse and the strand column is missing or "
                "not +/-: cannot pair primers into amplicons"
            )
        per_contig: Dict[str, List[Tuple[str, int, int, str, str, str]]] = {}
        for row in rows:
            per_contig.setdefault(row[0], []).append(row)
        for contig, contig_rows in per_contig.items():
            # In a tiling scheme LEFT n+1 lies before RIGHT n, so primers are
            # paired index-wise along the genome, not by alternation.
            lefts = sorted((r for r in contig_rows if r[5] == "+"), key=lambda r: (r[1], r[2]))
            rights = sorted((r for r in contig_rows if r[5] == "-"), key=lambda r: (r[1], r[2]))
            if len(lefts) != len(rights):
                raise NormaliseError(
                    f"scheme BED: {contig} has {len(lefts)} + and {len(rights)} - primers; cannot "
                    "pair them without <prefix>_<n>_LEFT/RIGHT names"
                )
            for n, (left, right) in enumerate(zip(lefts, rights), start=1):
                amp = by_key[(contig, "amplicon", n)] = Amplicon(
                    contig, f"amplicon_{n}", n, left[4]
                )
                amp.add_left(left[1], left[2])
                amp.add_right(right[1], right[2])
    for amp in by_key.values():
        if not amp.complete:
            scheme.notes.append(
                f"amplicon {amp.name} on {amp.contig} lacks a LEFT or RIGHT primer; ignored"
            )
            continue
        scheme.amplicons.setdefault(amp.contig, []).append(amp)
    for contig, amps in scheme.amplicons.items():
        amps.sort(key=lambda a: (a.outer_start, a.outer_end, a.n))
    for contig, start, end, _name, _pool, _strand in rows:
        scheme.primers.setdefault(contig, []).append((start, end))
    for primers in scheme.primers.values():
        primers.sort()
    if not scheme.amplicons:
        raise NormaliseError(f"scheme BED {path} yields no complete amplicon")
    return scheme


class _ContigAssigner:
    """Amplicon assignment and interior-primer detection for one contig."""

    def __init__(self, amplicons: Sequence[Amplicon], primers: Sequence[Tuple[int, int]]) -> None:
        self.amplicons = list(amplicons)
        self.starts = [a.outer_start for a in self.amplicons]
        self.max_len = max((a.outer_end - a.outer_start for a in self.amplicons), default=0)
        self.primers = list(primers)
        self.primer_starts = [p[0] for p in self.primers]

    def assign(self, start: int, end: int) -> Tuple[int, int]:
        """``(amplicon index or -1, class)`` for a record spanning ``[start, end)``."""
        tol = AMPLICON_TOLERANCE
        lo = bisect.bisect_left(self.starts, start - self.max_len - tol)
        hi = bisect.bisect_right(self.starts, start + tol)
        best = -1
        best_key: Optional[Tuple[int, int]] = None
        for j in range(lo, hi):
            amp = self.amplicons[j]
            if amp.outer_start - tol <= start and end <= amp.outer_end + tol:
                key = (abs(start - amp.insert_start) + abs(end - amp.insert_end), amp.n)
                if best_key is None or key < best_key:
                    best, best_key = j, key
        if best >= 0:
            # Fits one amplicon. A tiling scheme puts the neighbours' primers
            # inside this amplicon's insert, so containing a whole primer is
            # normal here and not a skip-amplicon signature.
            return best, CLASS_NORMAL
        k = bisect.bisect_left(self.primer_starts, start)
        while k < len(self.primers) and self.primers[k][0] < end:
            if self.primers[k][1] <= end:
                return -1, CLASS_INTERIOR_PRIMER
            k += 1
        return -1, CLASS_UNASSIGNED


# ---------------------------------------------------------------------------
# Record store
# ---------------------------------------------------------------------------


class Records:
    """Per-record fields of every alignment line, in input order, in compact
    arrays (about 30 bytes per record: a 5.8 M-record BAM measured 277 MB for the
    whole run, where Python ints would cost ~650 MB for the fields alone)."""

    def __init__(self) -> None:
        self.start = array.array("i")
        self.end = array.array("i")
        self.rank = array.array("i")  # start + aligned query bases: the selection key
        self.contig = array.array("i")  # index into header order; -1 = not a unit
        self.strand = array.array("b")
        self.cls = array.array("b")
        self.amplicon = array.array("i")
        self.mate = array.array("i")
        self.paired = bytearray()  # proper-pair record whose mate was looked for
        self.multi: Dict[int, List[Tuple[int, int]]] = {}  # records with N operations
        self.kept = bytearray()
        self.n = 0
        self.records_in = 0
        self.secondary_dropped = 0
        self.unmapped_dropped = 0
        self.orphan_mates = 0

    def append(
        self,
        start: int,
        end: int,
        rank: int,
        contig: int,
        strand: int,
        cls: int,
        amplicon: int,
        intervals: Optional[List[Tuple[int, int]]],
    ) -> int:
        self.start.append(start)
        self.end.append(end)
        self.rank.append(rank)
        self.contig.append(contig)
        self.strand.append(strand)
        self.cls.append(cls)
        self.amplicon.append(amplicon)
        self.mate.append(-1)
        self.paired.append(0)
        self.kept.append(0)
        if intervals is not None:
            self.multi[self.n] = intervals
        self.n += 1
        return self.n - 1

    def intervals(self, i: int) -> List[Tuple[int, int]]:
        multi = self.multi.get(i)
        return multi if multi is not None else [(self.start[i], self.end[i])]


def parse_header(header_lines: Sequence[str]) -> Tuple[List[str], List[int]]:
    names: List[str] = []
    lengths: List[int] = []
    for line in header_lines:
        if line.startswith("@SQ"):
            fields = dict(f.split(":", 1) for f in line.rstrip("\n").split("\t")[1:] if ":" in f)
            names.append(fields["SN"])
            lengths.append(int(fields["LN"]))
    if not names:
        raise NormaliseError("BAM header has no @SQ line")
    return names, lengths


def load_records(
    lines: Iterable[str], contig_names: Sequence[str], scheme: Optional[Scheme]
) -> Records:
    """Pass 1: parse every alignment line into :class:`Records`, pairing
    mates by QNAME and refusing input that is not coordinate-sorted. ``lines``
    must not include header lines."""
    recs = Records()
    contig_index = {name: i for i, name in enumerate(contig_names)}
    assigners: Dict[int, Optional[_ContigAssigner]] = {}
    open_pairs: Dict[str, int] = {}
    last_contig = -1
    last_pos = -1
    seen_contigs = set()
    for line in lines:
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t", 9)
        if len(fields) < 9:
            raise NormaliseError(f"malformed SAM line: {line[:80]!r}")
        recs.records_in += 1
        qname = fields[0]
        flag = int(fields[1])
        rname = fields[2]
        if flag & FLAG_UNMAPPED or rname == "*":
            recs.unmapped_dropped += 1
            recs.append(0, 0, 0, -1, 0, CLASS_NORMAL, -1, None)
            continue
        if flag & FLAG_SECONDARY:
            recs.secondary_dropped += 1
            recs.append(0, 0, 0, -1, 0, CLASS_NORMAL, -1, None)
            continue
        try:
            contig = contig_index[rname]
        except KeyError:
            raise NormaliseError(f"record on {rname!r}, which has no @SQ line") from None
        pos = int(fields[3]) - 1
        if contig != last_contig:
            if contig in seen_contigs or contig < last_contig:
                raise UnsortedInputError(
                    f"records of {rname} are not contiguous: input is not coordinate-sorted"
                )
            seen_contigs.add(contig)
            last_contig = contig
            last_pos = -1
            open_pairs.clear()
        elif pos < last_pos:
            raise UnsortedInputError(
                f"record {qname} at {rname}:{pos + 1} follows position {last_pos + 1}: input is not coordinate-sorted"
            )
        last_pos = pos
        intervals, query_bases = cigar_intervals(pos, fields[5])
        if not intervals:
            recs.unmapped_dropped += 1
            recs.append(0, 0, 0, -1, 0, CLASS_NORMAL, -1, None)
            continue
        # Ordering and pushing follow POS (the sort key); coverage follows the
        # intervals, which start after POS only when the CIGAR opens with N.
        start, end = pos, intervals[-1][1]
        multi = len(intervals) > 1 or intervals[0][0] != pos
        strand = 1 if flag & FLAG_REVERSE else 0
        cls, amplicon = CLASS_NORMAL, -1
        if scheme is not None:
            if contig not in assigners:
                amps = scheme.amplicons.get(rname)
                assigners[contig] = (
                    _ContigAssigner(amps, scheme.primers.get(rname, [])) if amps else None
                )
            assigner = assigners[contig]
            if assigner is None:
                cls = CLASS_UNASSIGNED
            else:
                amplicon, cls = assigner.assign(start, end)
        i = recs.append(
            start,
            end,
            pos + query_bases,
            contig,
            strand,
            cls,
            amplicon,
            intervals if multi else None,
        )
        pair_eligible = (
            flag & FLAG_PAIRED
            and flag & FLAG_PROPER
            and not flag & FLAG_MATE_UNMAPPED
            and not flag & FLAG_SUPPLEMENTARY
            and fields[6] in ("=", rname)
        )
        if pair_eligible:
            recs.paired[i] = 1
            mate = open_pairs.pop(qname, None)
            if mate is None:
                open_pairs[qname] = i
            else:
                recs.mate[i] = mate
                recs.mate[mate] = i
    recs.orphan_mates = sum(1 for i in range(recs.n) if recs.paired[i] and recs.mate[i] < 0)
    return recs


# ---------------------------------------------------------------------------
# The sweep
# ---------------------------------------------------------------------------


def split_target(total: int, parts: Sequence[int]) -> List[int]:
    """Split ``total`` among ``parts`` in proportion, each share capped at its
    part, by largest remainder (ties to the lower index). ``sum(parts) <=
    total`` returns the parts themselves."""
    s = sum(parts)
    if total >= s:
        return list(parts)
    if total <= 0:
        return [0] * len(parts)
    quotas = [total * part / s for part in parts]
    shares = [int(q) for q in quotas]
    remainder = total - sum(shares)
    order = sorted(range(len(parts)), key=lambda k: (-(quotas[k] - shares[k]), k))
    for k in order[:remainder]:
        shares[k] += 1
    return shares


class ContigResult:
    """Depth tables of one contig after the sweep."""

    def __init__(self, name: str, length: int) -> None:
        self.name = name
        self.length = length
        self.before_fwd = array.array("i", bytes(4 * length))
        self.before_rev = array.array("i", bytes(4 * length))
        self.after_fwd = array.array("i", bytes(4 * length))
        self.after_rev = array.array("i", bytes(4 * length))


class Normaliser:
    """Chooses the kept records of a :class:`Records` store."""

    def __init__(
        self,
        recs: Records,
        contig_names: Sequence[str],
        contig_lengths: Sequence[int],
        target: int,
        scheme: Optional[Scheme],
    ) -> None:
        self.recs = recs
        self.contig_names = list(contig_names)
        self.contig_lengths = list(contig_lengths)
        self.target = target
        self.scheme = scheme
        self.n_amplicons = (
            {c: len(scheme.amplicons.get(c, [])) for c in contig_names} if scheme else {}
        )
        self.results: List[ContigResult] = []
        self.positions_capped = 0

    # stratum ids: amplicon mode -> 2 * amplicon + strand, "other" = 2 * n_amplicons + strand;
    # position mode -> strand.
    def _stratum(self, i: int, n_amp: int) -> int:
        strand = self.recs.strand[i]
        if self.scheme is None:
            return strand
        if self.recs.cls[i] != CLASS_NORMAL or self.recs.amplicon[i] < 0:
            return 2 * n_amp + strand
        return 2 * self.recs.amplicon[i] + strand

    def run(self) -> None:
        recs = self.recs
        n = recs.n
        i = 0
        for c, (name, length) in enumerate(zip(self.contig_names, self.contig_lengths)):
            # The contig's records form one contiguous block (load_records refused anything else).
            while i < n and recs.contig[i] < 0:
                i += 1
            lo = i
            while i < n and (recs.contig[i] == c or recs.contig[i] < 0):
                i += 1
            hi = i
            ordinals = array.array("i", (k for k in range(lo, hi) if recs.contig[k] == c))
            self.results.append(self._sweep_contig(c, name, length, ordinals))

    def _sweep_contig(
        self, c: int, name: str, length: int, ordinals: Sequence[int]
    ) -> ContigResult:
        recs = self.recs
        target = self.target
        n_amp = self.n_amplicons.get(name, 0)
        result = ContigResult(name, length)
        strata: Dict[int, Tuple[array.array, array.array]] = {}

        def diffs(s: int) -> Tuple[array.array, array.array]:
            pair = strata.get(s)
            if pair is None:
                pair = strata[s] = (
                    array.array("i", bytes(4 * (length + 1))),
                    array.array("i", bytes(4 * (length + 1))),
                )
            return pair

        stratum_of = array.array("i", bytes(4 * recs.n)) if ordinals else array.array("i")
        for k in ordinals:
            s = self._stratum(k, n_amp)
            stratum_of[k] = s
            orig, _kept = diffs(s)
            for a, b in recs.intervals(k):
                if b > length:
                    raise NormaliseError(
                        f"a record on {name} reaches position {b}, past the @SQ length {length}: "
                        "the BAM header does not match the reference"
                    )
                orig[a] += 1
                orig[b] -= 1
        ids = sorted(strata)
        orig_cur = {s: 0 for s in ids}
        kept_cur = {s: 0 for s in ids}
        # Pending units per stratum: for every rank (start + aligned query
        # bases, a reach that an indel does not change: ranking by the true
        # reference end made every deletion-bearing read reach further than
        # its peers and systematically enriched deletions, benchmark v4/v5),
        # the ordinals of the units with that rank in input order (a compact
        # array, 4 bytes per unit), and a heap of the live ranks. Selection
        # order is furthest rank first, then input order; expiry uses the true
        # coverage end per record. At 30,000x every read of an amplicon starts
        # at the same clipped position, so millions of units are pending at
        # once; tuples per unit would cost ~140 bytes each.
        ends_heap: Dict[int, List[int]] = {s: [] for s in ids}
        groups: Dict[int, Dict[int, array.array]] = {s: {} for s in ids}
        heads: Dict[int, Dict[int, int]] = {s: {} for s in ids}
        kept = recs.kept
        end = recs.end
        mate = recs.mate

        def credit(k: int, p: int) -> None:
            """Mark record ``k`` kept and add its coverage: positions before
            ``p`` go straight into the finished rows (a first mate kept through
            its mate after it ended), ``p`` into the running sum, later
            positions into the difference array."""
            kept[k] = 1
            s = stratum_of[k]
            _orig, kdiff = diffs(s)
            after = result.after_rev if s & 1 else result.after_fwd
            for a, b in recs.intervals(k):
                past_end = min(b, p)
                if a < past_end:
                    for q in range(a, past_end):
                        after[q] += 1
                if b <= p:
                    continue
                if a <= p:
                    kept_cur[s] += 1
                else:
                    kdiff[a] += 1
                kdiff[b] -= 1

        def keep_unit(k: int, p: int) -> None:
            credit(k, p)
            m = mate[k]
            if m >= 0 and not kept[m]:
                credit(m, p)

        rank = recs.rank

        def push(k: int) -> None:
            s = stratum_of[k]
            r = rank[k]
            group = groups[s].get(r)
            if group is None:
                group = groups[s][r] = array.array("i")
                heads[s][r] = 0
                heapq.heappush(ends_heap[s], -r)
            group.append(k)

        def pop_live(s: int, p: int) -> int:
            """The furthest-ranked pending unit of stratum ``s`` whose coverage
            still reaches past ``p``, or -1. Units already kept through their
            mate and expired units (dropped unless a mate is still to come)
            are skipped; an exhausted rank group is discarded."""
            heap = ends_heap[s]
            while heap:
                r = -heap[0]
                group = groups[s][r]
                i = heads[s][r]
                while i < len(group) and (kept[group[i]] or end[group[i]] <= p):
                    i += 1
                if i >= len(group):
                    heapq.heappop(heap)
                    del groups[s][r]
                    del heads[s][r]
                    continue
                heads[s][r] = i + 1
                return group[i]
            return -1

        pending = 0  # next ordinal (in `ordinals`) to push
        n_ord = len(ordinals)
        capped = 0
        for p in range(length):
            for s in ids:
                orig, kdiff = strata[s]
                orig_cur[s] += orig[p]
                kept_cur[s] += kdiff[p]
            while pending < n_ord and recs.start[ordinals[pending]] == p:
                k = ordinals[pending]
                pending += 1
                if kept[k]:
                    continue  # pre-credited through its mate
                m = mate[k]
                if m >= 0 and kept[m]:
                    credit(k, p)
                    continue
                push(k)
            if pending < n_ord and recs.start[ordinals[pending]] < p:
                raise UnsortedInputError("record order violated during the sweep")
            orig_total = sum(orig_cur.values())
            if orig_total == 0:
                continue
            t_total = min(target, orig_total)
            if orig_total > target:
                capped += 1
            active = [s for s in ids if orig_cur[s] > 0]
            shares = split_target(t_total, [orig_cur[s] for s in active])
            for s, share in zip(active, shares):
                while kept_cur[s] < share:
                    k = pop_live(s, p)
                    if k < 0:
                        break
                    keep_unit(k, p)
            total_kept = sum(kept_cur.values())
            if total_kept < t_total:
                for s in active:
                    while total_kept < t_total:
                        k = pop_live(s, p)
                        if k < 0:
                            break
                        keep_unit(k, p)
                        total_kept = sum(kept_cur.values())
                    if total_kept >= t_total:
                        break
            # per-strand tallies for the depth table
            bf = br = af = ar = 0
            for s in ids:
                o, kk = orig_cur[s], kept_cur[s]
                if s & 1:
                    br += o
                    ar += kk
                else:
                    bf += o
                    af += kk
            result.before_fwd[p] = bf
            result.before_rev[p] = br
            result.after_fwd[p] = af
            result.after_rev[p] = ar
        self.positions_capped += capped
        return result


# ---------------------------------------------------------------------------
# Self-check and tables
# ---------------------------------------------------------------------------


def self_check(
    recs: Records, results: Sequence[ContigResult], contig_index: Dict[str, int], target: int
) -> None:
    """Recompute the kept and original depth from the record flags and refuse
    a kept set that breaks ``after >= min(target, before)`` anywhere, or whose
    depth differs from the sweep's own accounting."""
    per_contig: Dict[int, Tuple[array.array, array.array, array.array, array.array]] = {}
    for res in results:
        c = contig_index[res.name]
        per_contig[c] = tuple(array.array("i", bytes(4 * (res.length + 1))) for _ in range(4))  # type: ignore[assignment]
    for i in range(recs.n):
        c = recs.contig[i]
        if c < 0:
            continue
        bf, br, af, ar = per_contig[c]
        before = br if recs.strand[i] else bf
        after = ar if recs.strand[i] else af
        for a, b in recs.intervals(i):
            before[a] += 1
            before[b] -= 1
            if recs.kept[i]:
                after[a] += 1
                after[b] -= 1
    for res in results:
        bf, br, af, ar = per_contig[contig_index[res.name]]
        cb_f = cb_r = ca_f = ca_r = 0
        for p in range(res.length):
            cb_f += bf[p]
            cb_r += br[p]
            ca_f += af[p]
            ca_r += ar[p]
            if (cb_f, cb_r, ca_f, ca_r) != (
                res.before_fwd[p],
                res.before_rev[p],
                res.after_fwd[p],
                res.after_rev[p],
            ):
                raise SelfCheckError(
                    f"{res.name}:{p + 1}: sweep accounting ({res.before_fwd[p]},{res.before_rev[p]} -> "
                    f"{res.after_fwd[p]},{res.after_rev[p]}) disagrees with the kept flags "
                    f"({cb_f},{cb_r} -> {ca_f},{ca_r})"
                )
            if ca_f + ca_r < min(target, cb_f + cb_r):
                raise SelfCheckError(
                    f"{res.name}:{p + 1}: kept depth {ca_f + ca_r} is below min(target, original) = "
                    f"{min(target, cb_f + cb_r)}"
                )


def write_depth_table(results: Sequence[ContigResult], handle: TextIO) -> None:
    handle.write("\t".join(DEPTH_COLUMNS) + "\n")
    for res in results:
        for p in range(res.length):
            bf, br, af, ar = (
                res.before_fwd[p],
                res.before_rev[p],
                res.after_fwd[p],
                res.after_rev[p],
            )
            handle.write(f"{res.name}\t{p + 1}\t{bf + br}\t{af + ar}\t{bf}\t{br}\t{af}\t{ar}\n")


def _mean(total: int, length: int) -> str:
    return f"{total / length:.2f}" if length > 0 else NA


def amplicon_rows(
    recs: Records,
    results: Sequence[ContigResult],
    contig_index: Dict[str, int],
    scheme: Optional[Scheme],
) -> List[List[str]]:
    """Rows of the amplicon table: per amplicon (with a scheme) or per contig,
    plus an ``unassigned`` and an ``interior_primer`` row per contig in scheme
    mode. Mean depths are covered bases over the span (the insert for an
    amplicon, the contig otherwise), before and after."""
    # (contig, amplicon idx or -1, class) -> [records, covered bases, kept records, kept covered bases]
    acc: Dict[Tuple[int, int, int], List[int]] = {}
    for i in range(recs.n):
        c = recs.contig[i]
        if c < 0:
            continue
        key = (c, recs.amplicon[i], recs.cls[i])
        cell = acc.setdefault(key, [0, 0, 0, 0])
        span = sum(b - a for a, b in recs.intervals(i))
        cell[0] += 1
        cell[1] += span
        if recs.kept[i]:
            cell[2] += 1
            cell[3] += span

    def total(c: int, amp: Optional[int], cls: Optional[int]) -> List[int]:
        out = [0, 0, 0, 0]
        for (cc, a, k), cell in acc.items():
            if cc == c and (amp is None or a == amp) and (cls is None or k == cls):
                for j in range(4):
                    out[j] += cell[j]
        return out

    rows: List[List[str]] = []
    for res in results:
        c = contig_index[res.name]
        if scheme is None:
            n, bases, n_kept, bases_kept = total(c, None, None)
            rows.append(
                [
                    res.name,
                    "contig",
                    NA,
                    "0",
                    str(res.length),
                    _mean(bases, res.length),
                    _mean(bases_kept, res.length),
                    str(n),
                    str(n_kept),
                    "0",
                    "0",
                ]
            )
            continue
        for j, amp in enumerate(scheme.amplicons.get(res.name, [])):
            n, bases, n_kept, bases_kept = total(c, j, CLASS_NORMAL)
            span = max(1, amp.insert_end - amp.insert_start)
            rows.append(
                [
                    res.name,
                    amp.name,
                    amp.pool,
                    str(amp.outer_start),
                    str(amp.outer_end),
                    _mean(bases, span),
                    _mean(bases_kept, span),
                    str(n),
                    str(n_kept),
                    "0",
                    "0",
                ]
            )
        for cls, label in (
            (CLASS_UNASSIGNED, "unassigned"),
            (CLASS_INTERIOR_PRIMER, "interior_primer"),
        ):
            n, bases, n_kept, bases_kept = total(c, None, cls)
            interior = (str(n), str(n_kept)) if cls == CLASS_INTERIOR_PRIMER else ("0", "0")
            rows.append(
                [
                    res.name,
                    label,
                    NA,
                    "0",
                    str(res.length),
                    _mean(bases, res.length),
                    _mean(bases_kept, res.length),
                    str(n),
                    str(n_kept),
                    *interior,
                ]
            )
    return rows


def write_amplicon_table(rows: Sequence[Sequence[str]], handle: TextIO) -> None:
    handle.write("\t".join(AMPLICON_COLUMNS) + "\n")
    for row in rows:
        handle.write("\t".join(row) + "\n")


def summary_values(
    recs: Records, results: Sequence[ContigResult], target: int, mode: str, positions_capped: int
) -> Dict[str, str]:
    n_units = n_kept = 0
    span_before = span_after = 0
    unassigned = [0, 0]
    interior = [0, 0]
    for i in range(recs.n):
        if recs.contig[i] < 0:
            continue
        kept = recs.kept[i]
        span = sum(b - a for a, b in recs.intervals(i))
        n_units += 1
        span_before += span
        if kept:
            n_kept += 1
            span_after += span
        if recs.cls[i] == CLASS_UNASSIGNED:
            unassigned[0] += 1
            unassigned[1] += kept
        elif recs.cls[i] == CLASS_INTERIOR_PRIMER:
            interior[0] += 1
            interior[1] += kept
    positions = sum(r.length for r in results)
    return {
        "mode": mode,
        "target": str(target),
        "records_in": str(recs.records_in),
        "records_kept": str(n_kept),
        "records_dropped": str(recs.records_in - n_kept),
        "secondary_dropped": str(recs.secondary_dropped),
        "unmapped_dropped": str(recs.unmapped_dropped),
        "orphan_mates": str(recs.orphan_mates),
        "unassigned_before": str(unassigned[0]),
        "unassigned_after": str(unassigned[1]),
        "interior_primer_before": str(interior[0]),
        "interior_primer_after": str(interior[1]),
        "mean_span_before": _mean(span_before, n_units),
        "mean_span_after": _mean(span_after, n_kept),
        "positions_capped": str(positions_capped),
        "pct_positions_capped": f"{100.0 * positions_capped / positions:.2f}" if positions else NA,
        "self_check": "ok",
    }


def write_summary(values: Dict[str, str], handle: TextIO) -> None:
    for key in SUMMARY_KEYS:
        handle.write(f"{key}\t{values.get(key, NA)}\n")


def write_placeholders(depth_path: str, amplicons_path: str, summary_path: str) -> None:
    """Header-only tables and an ``off`` summary for ``--target 0``."""
    with open(depth_path, "w") as handle:
        handle.write("\t".join(DEPTH_COLUMNS) + "\n")
    with open(amplicons_path, "w") as handle:
        handle.write("\t".join(AMPLICON_COLUMNS) + "\n")
    with open(summary_path, "w") as handle:
        write_summary({"mode": MODE_OFF, "target": "0", "self_check": NA}, handle)


# ---------------------------------------------------------------------------
# Emission
# ---------------------------------------------------------------------------


def program_line(header_lines: Sequence[str], target: int, mode: str) -> str:
    """A ``@PG`` line without paths or timestamps, chained to the last ``@PG``."""
    previous = None
    for line in header_lines:
        if line.startswith("@PG"):
            for field in line.rstrip("\n").split("\t")[1:]:
                if field.startswith("ID:"):
                    previous = field[3:]
    fields = ["@PG", "ID:normalise_depth", "PN:normalise_depth"]
    if previous:
        fields.append(f"PP:{previous}")
    fields.append(f"CL:normalise_depth.py --target {target} --mode {mode}")
    return "\t".join(fields) + "\n"


def emit(lines: Iterable[str], recs: Records, out: TextIO, target: int, mode: str) -> int:
    """Pass 2: copy the header (``@HD`` stamped ``SO:coordinate``, our ``@PG``
    appended) and the kept alignment lines, in input order. Returns the
    number of records written."""
    header: List[str] = []
    written = 0
    ordinal = 0
    header_done = False
    for line in lines:
        if not header_done and line.startswith("@"):
            header.append(line)
            continue
        if not header_done:
            header_done = True
            _write_header(header, out, target, mode)
        if recs.kept[ordinal]:
            out.write(line)
            written += 1
        ordinal += 1
    if not header_done:
        _write_header(header, out, target, mode)
    if ordinal != recs.n:
        raise NormaliseError(
            f"second pass saw {ordinal} records, first pass {recs.n}: input changed underneath"
        )
    return written


def _write_header(header: List[str], out: TextIO, target: int, mode: str) -> None:
    hd_seen = False
    for line in header:
        if line.startswith("@HD"):
            hd_seen = True
            fields = [f for f in line.rstrip("\n").split("\t") if not f.startswith("SO:")]
            fields.append("SO:coordinate")
            out.write("\t".join(fields) + "\n")
        else:
            out.write(line)
    if not hd_seen:
        out.write("@HD\tVN:1.6\tSO:coordinate\n")
    out.write(program_line(header, target, mode))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def analyse(
    header: Sequence[str], record_lines: Iterable[str], target: int, scheme: Optional[Scheme]
) -> Tuple[Records, List[ContigResult], Dict[str, str], List[List[str]]]:
    """Pass 1 and the sweep over a stream of alignment lines: returns the
    records (with kept flags), per-contig depth tables, summary values and
    amplicon rows. Raises :class:`SelfCheckError` rather than returning a
    kept set that breaks an invariant. ``record_lines`` is consumed once and
    never stored, so the footprint is the compact :class:`Records` arrays."""
    names, lengths = parse_header(header)
    mode = MODE_AMPLICON if scheme is not None else MODE_POSITION
    recs = load_records(record_lines, names, scheme)
    normaliser = Normaliser(recs, names, lengths, target, scheme)
    normaliser.run()
    contig_index = {name: i for i, name in enumerate(names)}
    self_check(recs, normaliser.results, contig_index, target)
    values = summary_values(recs, normaliser.results, target, mode, normaliser.positions_capped)
    rows = amplicon_rows(recs, normaliser.results, contig_index, scheme)
    return recs, normaliser.results, values, rows


def normalise(
    sam_lines: Sequence[str], target: int, scheme: Optional[Scheme]
) -> Tuple[Records, List[ContigResult], Dict[str, str], List[List[str]], List[str]]:
    """:func:`analyse` over SAM text held in memory (tests); also returns the header lines."""
    header = [line for line in sam_lines if line.startswith("@")]
    recs, results, values, rows = analyse(
        header, (line for line in sam_lines if not line.startswith("@")), target, scheme
    )
    return recs, results, values, rows, header


def _split_header(lines: Iterator[str]) -> Tuple[List[str], Iterator[str]]:
    """Consume the header lines of a SAM stream; return them and an iterator
    over the remaining (alignment) lines."""
    header: List[str] = []
    first: Optional[str] = None
    for line in lines:
        if line.startswith("@"):
            header.append(line)
        else:
            first = line
            break

    def rest() -> Iterator[str]:
        if first is not None:
            yield first
        yield from lines

    return header, rest()


def _samtools_view(samtools: str, bam: str) -> subprocess.Popen:
    return subprocess.Popen(
        [samtools, "view", "-h", "--no-PG", bam],
        stdout=subprocess.PIPE,
        text=True,
        errors="surrogateescape",
    )


def _iter_lines(proc: subprocess.Popen) -> Iterator[str]:
    assert proc.stdout is not None
    for line in proc.stdout:
        yield line
    proc.stdout.close()
    if proc.wait() != 0:
        raise NormaliseError(f"samtools view exited with status {proc.returncode}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--bam", required=True, help="coordinate-sorted input BAM (the primer-clipped BAM)"
    )
    parser.add_argument(
        "--bai", default=None, help="its index; copied with --target 0 (default: <bam>.bai)"
    )
    parser.add_argument("--output-bam", required=True)
    parser.add_argument(
        "--output-bai",
        default=None,
        help="written by --target 0 only (default: <output-bam>.bai); otherwise index the output yourself",
    )
    parser.add_argument(
        "--target",
        type=int,
        required=True,
        help="records to keep per position; 0 switches normalisation off",
    )
    parser.add_argument(
        "--scheme",
        default=None,
        help="primer scheme BED; 'NA' or absent = position mode (strands only)",
    )
    parser.add_argument("--depth-tsv", required=True)
    parser.add_argument("--amplicons-tsv", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument(
        "--threads", type=int, default=1, help="threads for samtools view when writing the BAM"
    )
    parser.add_argument("--samtools", default="samtools")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.target < 0:
        print("--target must be >= 0", file=sys.stderr)
        return 2
    bai = args.bai or args.bam + ".bai"
    out_bai = args.output_bai or args.output_bam + ".bai"
    for path in (args.output_bam, args.depth_tsv, args.amplicons_tsv, args.summary):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    if args.target == 0:
        shutil.copyfile(args.bam, args.output_bam)
        shutil.copyfile(bai, out_bai)
        write_placeholders(args.depth_tsv, args.amplicons_tsv, args.summary)
        print(
            f"normalise_depth: off (target 0); copied {os.path.basename(args.bam)}", file=sys.stderr
        )
        return 0

    scheme: Optional[Scheme] = None
    if args.scheme and args.scheme != "NA":
        try:
            scheme = parse_scheme(args.scheme)
        except NormaliseError as exc:
            print(f"normalise_depth: {exc}", file=sys.stderr)
            return 1
        for note in scheme.notes:
            print(f"normalise_depth: scheme: {note}", file=sys.stderr)

    try:
        header, record_lines = _split_header(_iter_lines(_samtools_view(args.samtools, args.bam)))
        if scheme is not None:
            names, _lengths = parse_header(header)
            if not any(name in scheme.amplicons for name in names):
                print(
                    f"normalise_depth: WARNING: no contig of the BAM ({', '.join(names)}) appears in the scheme "
                    f"({', '.join(scheme.contigs())}); every read is 'unassigned' and kept",
                    file=sys.stderr,
                )
        recs, results, values, rows = analyse(header, record_lines, args.target, scheme)
    except UnsortedInputError as exc:
        print(f"normalise_depth: {exc}", file=sys.stderr)
        return EXIT_UNSORTED
    except NormaliseError as exc:
        print(f"normalise_depth: {exc}", file=sys.stderr)
        return 1
    except SelfCheckError as exc:
        print(f"normalise_depth: SELF-CHECK FAILED, nothing written: {exc}", file=sys.stderr)
        return EXIT_SELF_CHECK

    mode = values["mode"]
    writer = subprocess.Popen(
        [
            args.samtools,
            "view",
            "-b",
            "--no-PG",
            "-@",
            str(max(1, args.threads)),
            "-o",
            args.output_bam,
            "-",
        ],
        stdin=subprocess.PIPE,
        text=True,
        errors="surrogateescape",
    )
    assert writer.stdin is not None
    try:
        written = emit(
            _iter_lines(_samtools_view(args.samtools, args.bam)),
            recs,
            writer.stdin,
            args.target,
            mode,
        )
    except NormaliseError as exc:
        writer.stdin.close()
        writer.wait()
        print(f"normalise_depth: {exc}", file=sys.stderr)
        return 1
    writer.stdin.close()
    if writer.wait() != 0:
        print(
            f"normalise_depth: samtools view -b exited with status {writer.returncode}",
            file=sys.stderr,
        )
        return 1
    with open(args.depth_tsv, "w") as handle:
        write_depth_table(results, handle)
    with open(args.amplicons_tsv, "w") as handle:
        write_amplicon_table(rows, handle)
    with open(args.summary, "w") as handle:
        write_summary(values, handle)
    print(
        f"normalise_depth: {mode} mode, target {args.target}: kept {written} of {recs.records_in} records; "
        f"{values['positions_capped']} position(s) capped ({values['pct_positions_capped']} %); self-check ok",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
