"""Configuration file generation for the viralconseq workflows."""

import contextlib
import os
import tempfile
import textwrap
from typing import Any, Dict, List, Optional, Tuple, Union

import yaml

from viralconseq.constants import ConfigKeys, DataType, ResourceDefaults, ViralQCDatabase
from viralconseq.exceptions import ConfigurationError

# The viralQC database inputs ``rules/viralqc.smk`` declares, as paths relative
# to the placeholder directory. Derived from the database contract so the two
# cannot drift apart.
_VIRALQC_DB_PLACEHOLDERS = [f"viralqc_db/{name}" for name in ViralQCDatabase.REQUIRED_FILES]


class ConfigGenerator:
    """Generates YAML configuration files for Snakemake workflows."""

    # Sections of the written YAML, in output order, each with the comment
    # block printed above it. The YAML doubles as the run's documentation, so
    # a reader can tell what a key does and whether it is safe to edit.
    SECTIONS: List[Tuple[str, str]] = [
        (
            "run",
            "Samples (sample-<id> -> absolute FASTQ paths), data type, run directory "
            "and the thread baseline. Written from the sample sheet and the CLI; "
            "regenerate with --create-config-only rather than editing by hand.",
        ),
        (
            "consensus",
            "Reference-guided consensus calling: reference FASTA (or segment -> FASTA), "
            "primer scheme BED or NA, minimum depth for a called base, allele-frequency "
            "threshold. minimap2_consensus_align_flags is config-only (no CLI flag) and "
            "must hold plain tool flags.",
        ),
        (
            "read_qc",
            "Read filtering. Illumina: fastp adapters and trimming. minimum_length is "
            "fastp's --length_required on Illumina and the ampliconclip read-length "
            "filter on nanopore.",
        ),
        (
            "isnv",
            "Intra-host variant calling with LoFreq (Illumina only, run_isnv).",
        ),
        (
            "clair3",
            "Clair3 variant calling (nanopore): model, chunk size, variant quality and "
            "depth filters, minimum mapping quality.",
        ),
        (
            "viralqc",
            "Consensus QC with viralQC (virus, clade, genome-quality grade). "
            "viralqc_extra_flags is config-only and must hold plain tool flags.",
        ),
        (
            "provenance",
            "Written by viralconseq for the record; not read as a parameter.",
        ),
        (
            "resources",
            "Per-rule CPUs and RAM in GB; override on the command line with "
            "--<rule>-cpus / --<rule>-ram.",
        ),
    ]
    SECTION_RUN = "run"
    SECTION_CONSENSUS = "consensus"
    SECTION_READ_QC = "read_qc"
    SECTION_ISNV = "isnv"
    SECTION_CLAIR3 = "clair3"
    SECTION_VIRALQC = "viralqc"
    SECTION_PROVENANCE = "provenance"
    SECTION_RESOURCES = "resources"
    # Historical alias: keys with no better home land in the first section.
    SECTION_PARAMETERS = SECTION_RUN

    def __init__(self, config_path: str):
        """Initialize config generator.

        Args:
            config_path: Path where the config file will be written
        """
        self.config_path = config_path
        self.config: Dict[str, Any] = {}
        # Track which section each key belongs to
        self._sections: Dict[str, str] = {}

    def _set(self, key: str, value: Any, section: str) -> None:
        """Set a config key and tag it to a section.

        Args:
            key: Configuration key name
            value: Configuration value
            section: One of the ``SECTION_*`` names
        """
        self.config[key] = value
        self._sections[key] = section

    def add_samples(self, samples: Dict[str, List[str]], data_type: str) -> None:
        """Add samples to configuration.

        Args:
            samples: Dictionary mapping sample names to file paths
            data_type: Type of sequencing data (illumina or nanopore)
        """
        formatted_samples = {}
        for sample_name, file_paths in samples.items():
            key = f"sample-{sample_name}"
            if data_type == DataType.ILLUMINA:
                if len(file_paths) != 2:
                    raise ConfigurationError(
                        f"Illumina sample {sample_name} must have 2 files, "
                        f"found {len(file_paths)}"
                    )
                # Store R1/R2 as a list rather than a space-joined string so a
                # file path containing a space is not silently split by the
                # workflow's ``.split()`` on the sample value.
                formatted_samples[key] = [file_paths[0], file_paths[1]]
            else:
                if len(file_paths) != 1:
                    raise ConfigurationError(
                        f"Nanopore sample {sample_name} must have 1 file, "
                        f"found {len(file_paths)}"
                    )
                formatted_samples[key] = [file_paths[0]]

        self._set(ConfigKeys.SAMPLES, formatted_samples, self.SECTION_PARAMETERS)
        self._set(ConfigKeys.DATA, data_type, self.SECTION_PARAMETERS)

    def add_output(self, output_dir: str, run_name: str) -> None:
        """Add output directory to configuration.

        Args:
            output_dir: Base output directory
            run_name: Name of the run
        """
        self._set(
            ConfigKeys.OUTPUT,
            os.path.join(output_dir, run_name, ""),
            self.SECTION_PARAMETERS,
        )

    def add_threads(self, threads: int) -> None:
        """Add thread count to configuration.

        Args:
            threads: Number of threads
        """
        self._set(ConfigKeys.THREADS, threads, self.SECTION_PARAMETERS)

    def add_consensus_nanopore_settings(
        self,
        minimum_read_length: int,
        af_threshold: float,
        chunk_size: int,
        clair3_model: str,
        variant_quality: int,
        variant_depth: int,
        minimum_map_quality: int,
    ) -> None:
        """Add Nanopore consensus-specific settings to configuration.

        Args:
            minimum_read_length: Minimum read length threshold
            af_threshold: Allele frequency threshold to call a variant into consensus
            chunk_size: Size of chunks to process [clair3]
            clair3_model: Model to use for variant calling [clair3]
            variant_quality: Minimum variant quality to call a variant into consensus [clair3]
            variant_depth: Minimum alt allele depth to call a variant into consensus [clair3]
            minimum_map_quality: Minimum map quality to call a variant into consensus [clair3]
        """
        self._set(ConfigKeys.MINIMUM_LENGTH, minimum_read_length, self.SECTION_READ_QC)
        self._set(ConfigKeys.AF_THRESHOLD, af_threshold, self.SECTION_CONSENSUS)
        C = self.SECTION_CLAIR3
        self._set(ConfigKeys.CHUNK_SIZE, chunk_size, C)
        self._set(ConfigKeys.CLAIR3_MODEL, clair3_model, C)
        self._set(ConfigKeys.VARIANT_QUALITY, variant_quality, C)
        self._set(ConfigKeys.VARIANT_DEPTH, variant_depth, C)
        self._set(ConfigKeys.MINIMUM_MAP_QUALITY, minimum_map_quality, C)

    def add_illumina_settings(
        self,
        adapters: str,
        minimum_read_length: int,
        trim_head: Optional[int] = None,
        trim_tail: Optional[int] = None,
        cut_front_mean_quality: int = 20,
        cut_tail_mean_quality: int = 20,
        cut_right_window_size: int = 4,
        cut_right_mean_quality: int = 20,
        af_threshold: float = 0.5,
        af_isnv_threshold: float = 0.05,
        run_isnv: bool = False,
    ) -> None:
        """Add Illumina-specific settings to configuration (fastp QC).

        Args:
            adapters: Path to adapters file or "NA" for auto-detection
            minimum_read_length: Minimum read length threshold
            trim_head: Bases to trim from 5'
            trim_tail: Bases to trim from 3'
            cut_front_mean_quality: fastp cut_front mean quality threshold
            cut_tail_mean_quality: fastp cut_tail mean quality threshold
            cut_right_window_size: fastp cut_right window size
            cut_right_mean_quality: fastp cut_right mean quality threshold
            af_threshold: Allele frequency threshold to call a variant into consensus
            af_isnv_threshold: Minimum allele frequency threshold to call a variant into iSNV analysis
            run_isnv: Whether to run iSNV analysis
        """
        Q = self.SECTION_READ_QC
        self._set(ConfigKeys.ADAPTERS, adapters, Q)
        self._set(ConfigKeys.MINIMUM_LENGTH, minimum_read_length, Q)
        self._set(ConfigKeys.TRIM_HEAD, trim_head if trim_head is not None else 0, Q)
        self._set(ConfigKeys.TRIM_TAIL, trim_tail if trim_tail is not None else 0, Q)
        self._set(ConfigKeys.CUT_FRONT_MEAN_QUALITY, cut_front_mean_quality, Q)
        self._set(ConfigKeys.CUT_TAIL_MEAN_QUALITY, cut_tail_mean_quality, Q)
        self._set(ConfigKeys.CUT_RIGHT_WINDOW_SIZE, cut_right_window_size, Q)
        self._set(ConfigKeys.CUT_RIGHT_MEAN_QUALITY, cut_right_mean_quality, Q)
        self._set(ConfigKeys.AF_THRESHOLD, af_threshold, self.SECTION_CONSENSUS)
        self._set(ConfigKeys.AF_ISNV_THRESHOLD, af_isnv_threshold, self.SECTION_ISNV)
        self._set(ConfigKeys.RUN_ISNV, run_isnv, self.SECTION_ISNV)

    def add_consensus_settings(
        self,
        reference: Union[str, Dict[str, str]],
        primer_scheme: str,
        minimum_coverage: int,
        minimap2_consensus_align_flags: str = "-a --sam-hit-only --secondary=no --score-N=0",
    ) -> None:
        """Add consensus-specific settings to configuration.

        Args:
            reference: Path to reference genome (str) or dict mapping
            segment names to paths for segmented viruses
            primer_scheme: Path to primer scheme file or "NA"
            minimum_coverage: Minimum coverage for consensus
            minimap2_consensus_align_flags: Flags passed to minimap2 when
                re-aligning the per-sample consensus FASTA back to the
                reference for the final multiple-sequence alignment. The
                default keeps the historical behaviour.
        """
        P = self.SECTION_CONSENSUS
        self._set(ConfigKeys.REFERENCE, reference, P)
        self._set(ConfigKeys.SCHEME, primer_scheme, P)
        self._set(ConfigKeys.MINIMUM_DEPTH, minimum_coverage, P)
        self._set(
            ConfigKeys.MINIMAP2_CONSENSUS_ALIGN_FLAGS,
            minimap2_consensus_align_flags,
            P,
        )

    def add_viralqc_settings(
        self,
        run_viralqc: bool = True,
        viralqc_db: str = "NA",
        viralqc_extra_flags: str = "",
    ) -> None:
        """Add consensus-QC (viralQC) settings to configuration.

        Args:
            run_viralqc: Run viralQC on the final consensus sequences.
            viralqc_db: Absolute path to the viralQC database directory. Written
                even when ``run_viralqc`` is False so the YAML can be switched on
                later; ``"NA"`` only when no directory is known.
            viralqc_extra_flags: Extra flags appended to the ``vqc run`` command
                line. Config-only (no CLI flag), like
                ``minimap2_consensus_align_flags``.
        """
        P = self.SECTION_VIRALQC
        self._set(ConfigKeys.RUN_VIRALQC, run_viralqc, P)
        self._set(ConfigKeys.VIRALQC_DB, viralqc_db, P)
        self._set(ConfigKeys.VIRALQC_EXTRA_FLAGS, viralqc_extra_flags, P)

    def add_workflow_path(self, workflow_path: str) -> None:
        """Add workflow path to configuration.

        Args:
            workflow_path: Path to the workflow directory
        """
        self._set("workflow_path", workflow_path, self.SECTION_PARAMETERS)

    def add_resource_settings(self, args: Dict[str, Any], rule_names: list) -> None:
        """Add per-rule resource settings (CPUs and RAM) to configuration.

        For each rule name in rule_names, writes ``<rule>_cpus`` and
        ``<rule>_ram`` keys to the config dict.  Values are taken from
        *args* if present; otherwise the defaults from
        ``ResourceDefaults`` are used.

        Args:
            args: Dictionary of pipeline arguments (from the CLI).
            rule_names: List of Snakemake rule name strings that should
                receive resource entries.
        """
        R = self.SECTION_RESOURCES
        for rule in rule_names:
            cpus_key = f"{rule}_cpus"
            ram_key = f"{rule}_ram"
            self._set(cpus_key, args.get(cpus_key, ResourceDefaults.DEFAULT_CPUS), R)
            self._set(ram_key, args.get(ram_key, ResourceDefaults.DEFAULT_RAM), R)

    @classmethod
    def write_skeleton(
        cls,
        pipeline: str,
        data_type: str,
        config_path: str,
        placeholder_dir: str,
    ) -> str:
        """Write a placeholder YAML config that lets Snakemake parse a workflow.

        Used by ``viralconseq setup`` to materialize per-rule conda envs via
        ``snakemake(..., conda_create_envs_only=True)``. Even in
        envs-only mode Snakemake walks the DAG and verifies that rule
        inputs exist, so the caller is expected to populate
        ``placeholder_dir`` with empty files for the FASTQs and references
        emitted here. ``setup_cli`` does
        this for the package; tests use the same convention.

        Args:
            pipeline: ``"consensus"`` (the only supported value; kept as an
                explicit argument so the skeleton API stays stable).
            data_type: ``"illumina"`` or ``"nanopore"``.
            config_path: Where to write the YAML.
            placeholder_dir: Directory that contains the placeholder input
                files referenced by the generated config.

        Returns:
            ``config_path``.
        """
        gen = cls(config_path)
        root = placeholder_dir.rstrip("/")

        if data_type == DataType.ILLUMINA:
            placeholder_samples = {
                "skel": [
                    f"{root}/reads/skel_R1.fastq.gz",
                    f"{root}/reads/skel_R2.fastq.gz",
                ],
            }
        else:
            placeholder_samples = {"skel": [f"{root}/reads/skel.fastq.gz"]}

        gen.add_samples(placeholder_samples, data_type)
        gen.add_output(f"{root}/output", "skeleton_run")
        gen.add_threads(1)

        # Every optional feature flag is enabled below so Snakemake's DAG
        # walk reaches every rule and ``conda_create_envs_only=True``
        # materializes every per-rule env that pipeline could possibly need.
        # If a flag is left off here, its env is silently skipped by
        # ``viralconseq setup`` and falls back to dynamic env creation at
        # the user's first real run — exactly the failure mode setup
        # exists to prevent.
        if pipeline == "consensus":
            gen.add_consensus_settings(
                reference=f"{root}/references/skel.reference.fasta",
                primer_scheme="NA",
                minimum_coverage=20,
            )
            if data_type == DataType.ILLUMINA:
                # run_isnv=True pulls the detect_isnv (LoFreq) rule into
                # the DAG, which uses envs/consensus.yaml.
                gen.add_illumina_settings(
                    adapters="NA",
                    minimum_read_length=50,
                    run_isnv=True,
                )
            else:
                gen.add_consensus_nanopore_settings(
                    minimum_read_length=50,
                    af_threshold=0.51,
                    chunk_size=10000,
                    clair3_model="r1041_e82_400bps_sup_v500",
                    variant_quality=20,
                    variant_depth=10,
                    minimum_map_quality=30,
                )
            # run_viralqc=True pulls run_viralqc (envs/viralqc.yaml) into the
            # DAG; the placeholder DB files are listed in SKELETON_PLACEHOLDERS.
            gen.add_viralqc_settings(run_viralqc=True, viralqc_db=f"{root}/viralqc_db")
            gen.add_workflow_path(".")
        else:
            raise ValueError(f"Unknown pipeline: {pipeline!r} (expected 'consensus')")

        gen.save()
        return config_path

    # Empty files (relative to ``placeholder_dir``) that
    # ``ConfigGenerator.write_skeleton`` expects to exist on disk for
    # ``snakemake --conda-create-envs-only`` to succeed. Kept here so
    # the skeleton and its placeholder-file expectations stay in lock-step.
    SKELETON_PLACEHOLDERS: Dict[str, Dict[str, List[str]]] = {
        "consensus": {
            "illumina": [
                "reads/skel_R1.fastq.gz",
                "reads/skel_R2.fastq.gz",
                "references/skel.reference.fasta",
                *_VIRALQC_DB_PLACEHOLDERS,
            ],
            "nanopore": [
                "reads/skel.fastq.gz",
                "references/skel.reference.fasta",
                *_VIRALQC_DB_PLACEHOLDERS,
            ],
        },
    }

    def save(self) -> None:
        """Write the YAML atomically, one commented section at a time.

        Keys are grouped by the section tags assigned in ``_set()`` and written
        in ``SECTIONS`` order; empty sections are skipped and insertion order is
        kept within a section. The file is written to a temporary sibling and
        renamed into place, so a crash mid-write never leaves a half config.

        Raises:
            ConfigurationError: If the file or its directory cannot be written.
        """
        config_dir = os.path.dirname(self.config_path)
        if config_dir:  # Only create directory if path contains a directory component
            os.makedirs(config_dir, exist_ok=True)

        grouped: Dict[str, Dict[str, Any]] = {name: {} for name, _ in self.SECTIONS}
        for key, value in self.config.items():
            section = self._sections.get(key, self.SECTION_RUN)
            grouped.setdefault(section, {})[key] = value

        tmp_path: Optional[str] = None
        try:
            fd, tmp_path = tempfile.mkstemp(
                dir=config_dir or ".", prefix=".config-", suffix=".tmp", text=True
            )
            with os.fdopen(fd, "w") as f:
                f.write("# viralconseq run configuration. Generated by the CLI; the Snakemake\n")
                f.write("# workflows read exactly these keys.\n")
                for section, comment in self.SECTIONS:
                    items = grouped.get(section)
                    if not items:
                        continue
                    f.write(f"\n# --- {section} ---\n")
                    for line in textwrap.wrap(comment, width=76):
                        f.write(f"# {line}\n")
                    yaml.dump(items, f, default_flow_style=False, sort_keys=False)
            os.replace(tmp_path, self.config_path)
        except (OSError, IOError) as e:
            if tmp_path is not None:
                with contextlib.suppress(OSError):
                    os.unlink(tmp_path)
            raise ConfigurationError(
                f"Failed to write config file to {self.config_path}: {e}"
            ) from e
