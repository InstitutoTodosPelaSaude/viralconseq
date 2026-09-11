"""Constants used throughout the viralconseq pipeline."""


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
    ]

    # Nanopore computational rules
    CONSENSUS_NANOPORE_RULES = [
        "map_reads",
        "trim_primer_sequences",
        "infer_consensus_sequence",
    ]
