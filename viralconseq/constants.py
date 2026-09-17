"""Constants used throughout the viralconseq pipeline."""

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class DataType:
    """Sequencing data types supported by viralconseq."""

    ILLUMINA = "illumina"
    NANOPORE = "nanopore"

    @classmethod
    def is_valid(cls, value: str) -> bool:
        """Check if a value is a valid data type."""
        return value in [cls.ILLUMINA, cls.NANOPORE]


class ConfigKeys:
    """Keys used in configuration files."""

    SAMPLES = "samples"
    DATA = "data"
    OUTPUT = "output"
    THREADS = "threads"
    REFERENCE = "reference"
    SCHEME = "scheme"
    ADAPTERS = "adapters"
    MINIMUM_LENGTH = "minimum_length"
    MINIMUM_DEPTH = "minimum_depth"
    TRIM_HEAD = "trim_head"
    TRIM_TAIL = "trim_tail"
    CUT_FRONT_MEAN_QUALITY = "cut_front_mean_quality"
    CUT_TAIL_MEAN_QUALITY = "cut_tail_mean_quality"
    CUT_RIGHT_WINDOW_SIZE = "cut_right_window_size"
    CUT_RIGHT_MEAN_QUALITY = "cut_right_mean_quality"
    MINIMAP2_CONSENSUS_ALIGN_FLAGS = "minimap2_consensus_align_flags"

    AF_THRESHOLD = "af_threshold"
    AF_ISNV_THRESHOLD = "af_isnv_threshold"
    CHUNK_SIZE = "chunk_size"
    CLAIR3_MODEL = "clair3_model"
    VARIANT_QUALITY = "variant_quality"
    VARIANT_DEPTH = "variant_depth"
    MINIMUM_MAP_QUALITY = "minimum_map_quality"
    RUN_ISNV = "run_isnv"
    RUN_VIRALQC = "run_viralqc"
    VIRALQC_DB = "viralqc_db"
    VIRALQC_EXTRA_FLAGS = "viralqc_extra_flags"
    VIRALCONSEQ_VERSION = "viralconseq_version"
    THREADS_TOTAL = "threads_total"
    MAX_MEMORY_MB = "max_memory_mb"
    MEMORY_DETECTED_MB = "memory_detected_mb"


class ViralQCDatabase:
    """Layout contract of the viralQC database directory.

    Shared by ``validators.validate_viralqc_database``, ``viralconseq setup``
    and ``scripts/viralqc_setup.smk`` (whose ``rule all`` produces exactly these
    entries). ``rules/viralqc.smk`` declares the files as rule inputs.
    """

    ENV_VAR = "VIRALCONSEQ_VIRALQC_DB"
    NEXTCLADE_SENTINEL = ".nextclade_datasets_ok"
    REQUIRED_FILES = ("blast.fasta", "blast.tsv", NEXTCLADE_SENTINEL)
    REQUIRED_DIRS = ("blast_gff",)

    # ``makeblastdb`` writes ``blast.fasta.nin`` (and siblings) next to the
    # FASTA, or one ``blast.fasta.<NN>.nin`` per volume plus a
    # ``blast.fasta.nal`` alias when the database is split. The index must be
    # validated separately from ``blast.fasta`` itself: without it ``blastn``
    # fails, and viralQC deliberately degrades a blastn failure to "no hits",
    # so the run would succeed with every unmatched sequence silently reported
    # as Unclassified.
    BLAST_INDEX_GLOBS = ("blast.fasta*.nin", "blast.fasta.nal")
    BLAST_INDEX_LABEL = "blast.fasta.nin"

    @classmethod
    def default_dir(cls) -> str:
        """Default database directory, resolved at call time.

        ``$VIRALCONSEQ_VIRALQC_DB`` first, then ``~/.cache/viralconseq/viralqc-db``.
        Used as the click default of ``--viralqc-db`` on ``consensus`` and
        ``setup`` so both commands agree on the location.
        """
        # ``or`` (not the get() default) so an exported-but-empty
        # VIRALCONSEQ_VIRALQC_DB falls back instead of resolving to "".
        return os.environ.get(cls.ENV_VAR, "") or str(
            Path.home() / ".cache" / "viralconseq" / "viralqc-db"
        )


class SampleSheetPattern:
    """Patterns for identifying sample files."""

    R1 = "R1"
    BARCODE = "barcode"


class SampleSheetSeparator:
    """Separators for parsing sample names."""

    UNDERSCORE = "_"
    HYPHEN = "-"
    DOT = "."


class ResourceDefaults:
    """Resource contract for the computational Snakemake rules.

    Threads: every rule in ``CONSENSUS_*_RULES`` has a ``--<rule>-cpus`` option;
    the rule reads ``<rule>_cpus`` when the operator set it and falls back to
    the ``--threads`` baseline (``cpus()`` in ``rules/common.smk``).

    Memory: only rules with a *measured* peak declare ``mem_mb``, and they are
    the only ones with a ``--<rule>-ram`` option. Their figures (GB) come from
    ``benchmark.tsv`` ``max_rss`` on the SARS-CoV-2 test data, rounded up with
    headroom; never add or change one that was not read from a run. A budget
    (``--max-memory``, detected from the machine by default) is handed to
    Snakemake as ``--resources mem_mb=``, so these declarations bind: rules
    declaring memory run concurrently only while their sum fits the budget.
    """

    DEFAULT_CPUS = 2
    #: Fraction of detected memory left to the OS and unmeasured rules.
    MEMORY_HEADROOM = 0.10
    #: rule -> GB. Clair3 peaked at 0.53 GB, viralQC (nextclade + BLAST) at 0.38 GB
    #: on two samples; both scale with input, hence the headroom.
    MEMORY_RULES: Dict[str, int] = {
        "infer_consensus_sequence": 2,
        "run_viralqc": 1,
    }

    # Illumina computational rules
    CONSENSUS_ILLUMINA_RULES = [
        "perform_qc",
        "map_reads",
        "trim_primer_sequences",
        "detect_isnv",
        "run_viralqc",
    ]

    # Nanopore computational rules
    CONSENSUS_NANOPORE_RULES = [
        "map_reads",
        "trim_primer_sequences",
        "infer_consensus_sequence",
        "run_viralqc",
    ]

    @classmethod
    def ram_for(cls, rule: str) -> int:
        """Default GB for a memory-declaring rule (KeyError for any other)."""
        return cls.MEMORY_RULES[rule]

    @classmethod
    def memory_rules_for(cls, rule_list: List[str]) -> List[str]:
        """The rules of ``rule_list`` that declare memory, in list order."""
        return [rule for rule in rule_list if rule in cls.MEMORY_RULES]


# ----------------------------------------------------------- machine detection
#
# Both helpers answer "what may this run use", not "what does this machine
# have", so they consult the limits imposed on the process (CPU affinity, the
# process's own cgroup and its ancestors) before the hardware totals.


def detect_cores() -> int:
    """Cores to give Snakemake by default: those available to this process, less one.

    ``sched_getaffinity`` rather than ``cpu_count``: the latter reports the
    host's CPUs even when cgroup or container limits confine the process to
    fewer (the WSL and cluster case). One core is left free so a long run does
    not take the whole machine with it; ``--threads-total`` overrides this.
    """
    try:
        available = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        # AttributeError: not Linux. OSError: a sandbox that blocks the call.
        available = os.cpu_count() or 1
    return max(1, available - 1)


#: cgroup v1 writes 2^63 rounded down to a page for "unlimited"; v2 writes the
#: literal ``max``. Anything this large is not a limit.
_NO_LIMIT_BYTES = 1 << 62


def _read_text(path: str) -> Optional[str]:
    """One file read; the seam the tests replace. None when absent or unreadable."""
    try:
        with open(path) as handle:
            return handle.read()
    except OSError:
        return None


def _cgroup_memory_limit_bytes() -> Optional[int]:
    """The tightest memory limit on this process's cgroup or any of its ancestors.

    ``/proc/self/cgroup`` names the cgroup the process is in. cgroup v2: one
    ``0::/path`` line, limit at ``/sys/fs/cgroup/<path>/memory.max``. v1: the
    line whose controller list contains ``memory``, limit at
    ``/sys/fs/cgroup/memory/<path>/memory.limit_in_bytes``. A limit may sit on
    any ancestor (a Slurm job slice, a Docker container), so the path is walked
    upwards and the minimum kept. ``max``, unparsable text and values at or
    above 2^62 mean no limit.
    """
    text = _read_text("/proc/self/cgroup")
    if not text:
        return None
    candidates = []
    for line in text.splitlines():
        parts = line.split(":", 2)
        if len(parts) != 3:
            continue
        _, controllers, path = parts
        if controllers == "":
            root, leaf = "/sys/fs/cgroup", "memory.max"
        elif "memory" in controllers.split(","):
            root, leaf = "/sys/fs/cgroup/memory", "memory.limit_in_bytes"
        else:
            continue
        segments = [segment for segment in path.split("/") if segment]
        while True:
            candidates.append("/".join([root, *segments, leaf]))
            if not segments:
                break
            segments.pop()
    limits = []
    for candidate in candidates:
        value = _read_text(candidate)
        if value is None:
            continue
        value = value.strip()
        if not value.isdigit():
            continue
        number = int(value)
        if 0 < number < _NO_LIMIT_BYTES:
            limits.append(number)
    return min(limits) if limits else None


def detect_memory_mb() -> Optional[int]:
    """Memory this process may use, in MB, or None when it cannot be known.

    Linux: ``MemTotal`` from ``/proc/meminfo``, capped by the cgroup limit when
    one is set. Elsewhere: physical pages times page size via ``os.sysconf``.
    None is a real answer (the CLI then runs without a budget and says so)
    rather than a guess.
    """
    total_mb = None
    meminfo = _read_text("/proc/meminfo")
    if meminfo:
        for line in meminfo.splitlines():
            if line.startswith("MemTotal:"):
                fields = line.split()
                if len(fields) >= 2 and fields[1].isdigit():
                    total_mb = int(fields[1]) // 1024
                break
    if total_mb is None:
        try:
            total_mb = (os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")) // (1024 * 1024)
        except (AttributeError, OSError, ValueError):
            return None
    limit = _cgroup_memory_limit_bytes()
    if limit is not None:
        total_mb = min(total_mb, limit // (1024 * 1024))
    return total_mb if total_mb > 0 else None


def memory_budget_mb(
    detected_mb: int,
    largest_rule_gb: int,
    headroom: float = ResourceDefaults.MEMORY_HEADROOM,
) -> Tuple[int, bool]:
    """The budget to hand Snakemake, and whether it had to be raised to fit one job.

    ``detected_mb`` less ``headroom``; if that is below the largest per-rule
    figure the budget becomes exactly that figure, so memory-declaring rules
    run one at a time instead of the run being refused. Pure arithmetic so the
    clamp is testable without touching ``/proc``.
    """
    budget = int(detected_mb * (1 - headroom))
    floor = int(largest_rule_gb) * 1024
    if budget < floor:
        return floor, True
    return budget, False
