"""Constants used throughout the viralconseq pipeline."""

import os
from pathlib import Path


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
    """Default resource allocation for computational Snakemake rules.

    Each rule gets <rule>_cpus and <rule>_ram (GB) in the config YAML.
    """

    DEFAULT_CPUS = 2
    DEFAULT_RAM = 4  # GB

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
