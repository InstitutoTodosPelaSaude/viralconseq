"""Configuration file generation for the viralconseq workflows."""

import contextlib
import os
import shutil
import tempfile
import textwrap
from typing import Any, Dict, List, Optional, Tuple, Union

import yaml

from viralconseq import __version__
from viralconseq.constants import (
    Clair3Models,
    ConfigKeys,
    DataType,
    ResourceDefaults,
    ViralQCDatabase,
)
from viralconseq.exceptions import ConfigurationError

# The viralQC database inputs ``rules/viralqc.smk`` declares, as paths relative
# to the placeholder directory. Derived from the database contract so the two
# cannot drift apart.
_VIRALQC_DB_PLACEHOLDERS = [f"viralqc_db/{name}" for name in ViralQCDatabase.REQUIRED_FILES]
# The Clair3 model inputs ``rules/consensus_nanopore.smk`` declares (empty files
# suffice: the DAG needs them to exist, validation is Python pre-flight only).
_SKELETON_CLAIR3_MODEL = "r1041_e82_400bps_sup_v500"
_CLAIR3_MODEL_PLACEHOLDERS = [
    f"clair3_models/{_SKELETON_CLAIR3_MODEL}/{checkpoint}"
    for checkpoint in Clair3Models.CHECKPOINTS
]


# Section of each known key, for configs loaded back from disk (``from_dict``).
_KEY_SECTIONS: Dict[str, str] = {
    ConfigKeys.SAMPLES: "run",
    ConfigKeys.DATA: "run",
    ConfigKeys.OUTPUT: "run",
    ConfigKeys.THREADS: "run",
    "workflow_path": "run",
    ConfigKeys.REFERENCE: "consensus",
    ConfigKeys.SCHEME: "consensus",
    ConfigKeys.MINIMUM_DEPTH: "consensus",
    ConfigKeys.MINIMAP2_CONSENSUS_ALIGN_FLAGS: "consensus",
    ConfigKeys.AF_THRESHOLD: "consensus",
    ConfigKeys.CONSENSUS_COVERAGE_THRESHOLD: "consensus",
    ConfigKeys.ADAPTERS: "read_qc",
    ConfigKeys.MINIMUM_LENGTH: "read_qc",
    ConfigKeys.TRIM_HEAD: "read_qc",
    ConfigKeys.TRIM_TAIL: "read_qc",
    ConfigKeys.CUT_FRONT_MEAN_QUALITY: "read_qc",
    ConfigKeys.CUT_TAIL_MEAN_QUALITY: "read_qc",
    ConfigKeys.CUT_RIGHT_WINDOW_SIZE: "read_qc",
    ConfigKeys.CUT_RIGHT_MEAN_QUALITY: "read_qc",
    ConfigKeys.AF_ISNV_THRESHOLD: "isnv",
    ConfigKeys.RUN_ISNV: "isnv",
    ConfigKeys.CHUNK_SIZE: "clair3",
    ConfigKeys.CLAIR3_MODEL: "clair3",
    ConfigKeys.CLAIR3_MODEL_DIR: "clair3",
    ConfigKeys.MINIMUM_MAPPED_READS: "clair3",
    ConfigKeys.VARIANT_QUALITY: "clair3",
    ConfigKeys.VARIANT_DEPTH: "clair3",
    ConfigKeys.MINIMUM_MAP_QUALITY: "clair3",
    ConfigKeys.RUN_VIRALQC: "viralqc",
    ConfigKeys.VIRALQC_DB: "viralqc",
    ConfigKeys.VIRALQC_EXTRA_FLAGS: "viralqc",
    ConfigKeys.RUN_REPORT: "report",
    ConfigKeys.VIRALCONSEQ_VERSION: "provenance",
    ConfigKeys.THREADS_TOTAL: "resources",
    ConfigKeys.MAX_MEMORY_MB: "resources",
    ConfigKeys.MEMORY_DETECTED_MB: "resources",
}


def sample_key(sample_name: str) -> str:
    """The config/output key of a sample: ``sample-<id>``."""
    return f"sample-{sample_name}"


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
            "report",
            "The self-contained run report (report.html): run_report toggles it "
            "(--report / --no-report).",
        ),
        (
            "provenance",
            "Written by viralconseq for the record; not read as a parameter.",
        ),
        (
            "resources",
            "Run-wide budgets given to Snakemake (threads_total cores; max_memory_mb, "
            "0 = no memory budget) and per-rule overrides: <rule>_cpus (absent = the "
            "threads baseline) and <rule>_ram in GB for the rules that declare memory. "
            "Override with --<rule>-cpus / --<rule>-ram / --threads-total / --max-memory.",
        ),
    ]
    SECTION_RUN = "run"
    SECTION_CONSENSUS = "consensus"
    SECTION_READ_QC = "read_qc"
    SECTION_ISNV = "isnv"
    SECTION_CLAIR3 = "clair3"
    SECTION_VIRALQC = "viralqc"
    SECTION_REPORT = "report"
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

    @classmethod
    def section_for(cls, key: str) -> str:
        """The section a config key belongs to (for configs loaded from disk)."""
        if key.endswith("_cpus") or key.endswith("_ram"):
            return cls.SECTION_RESOURCES
        return _KEY_SECTIONS.get(key, cls.SECTION_RUN)

    @classmethod
    def from_dict(cls, config_path: str, config: Dict[str, Any]) -> "ConfigGenerator":
        """A generator holding an existing config (e.g. one loaded from YAML),
        with every key placed in its section, ready to ``save()``."""
        generator = cls(config_path)
        for key, value in config.items():
            generator._set(key, value, cls.section_for(key))
        return generator

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
            key = sample_key(sample_name)
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
        clair3_model: Union[str, Dict[str, str]],
        variant_quality: int,
        variant_depth: int,
        minimum_map_quality: int,
        clair3_model_dir: str = "",
        minimum_mapped_reads: int = 10,
    ) -> None:
        """Add Nanopore consensus-specific settings to configuration.

        Args:
            minimum_read_length: Minimum read length threshold
            af_threshold: Allele fraction threshold to call a variant into consensus,
                ALT reads over REF plus ALT reads at the site
            chunk_size: Size of chunks to process [clair3]
            clair3_model: Model for variant calling [clair3]: one name for every
                sample, or ``{sample id: name}`` (re-keyed ``sample-<id>`` like
                ``samples``) when samples were basecalled with different models
            clair3_model_dir: Directory holding ``<model>/pileup.pt`` and
                ``full_alignment.pt`` (``constants.Clair3Models``)
            minimum_mapped_reads: Below this many primary mapped reads a sample
                is not sent to Clair3 but gets an all-N consensus and a
                ``no_mapped_reads`` status (0 disables the floor)
            variant_quality: Minimum variant quality to call a variant into consensus [clair3]
            variant_depth: Minimum alt allele depth to call a variant into consensus [clair3]
            minimum_map_quality: Minimum map quality to call a variant into consensus [clair3]
        """
        self._set(ConfigKeys.MINIMUM_LENGTH, minimum_read_length, self.SECTION_READ_QC)
        self._set(ConfigKeys.AF_THRESHOLD, af_threshold, self.SECTION_CONSENSUS)
        C = self.SECTION_CLAIR3
        self._set(ConfigKeys.CHUNK_SIZE, chunk_size, C)
        if isinstance(clair3_model, dict):
            clair3_model = {sample_key(sample): name for sample, name in clair3_model.items()}
        self._set(ConfigKeys.CLAIR3_MODEL, clair3_model, C)
        if clair3_model_dir:
            self._set(ConfigKeys.CLAIR3_MODEL_DIR, clair3_model_dir, C)
        self._set(ConfigKeys.MINIMUM_MAPPED_READS, int(minimum_mapped_reads), C)
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

    def add_collect_settings(self, consensus_coverage_threshold: float = 70.0) -> None:
        """Settings of the collection step (``rules/collect.smk``).

        Args:
            consensus_coverage_threshold: ``coverage_min_depth`` (percent) a
                sample needs to enter ``consensus/consensus.cov<T>.fasta``.
        """
        self._set(
            ConfigKeys.CONSENSUS_COVERAGE_THRESHOLD,
            float(consensus_coverage_threshold),
            self.SECTION_CONSENSUS,
        )

    def add_report_settings(self, run_report: bool = True) -> None:
        """Toggle the self-contained run report (``rules/report.smk``)."""
        self._set(ConfigKeys.RUN_REPORT, bool(run_report), self.SECTION_REPORT)

    def add_provenance(self, viralconseq_version: str) -> None:
        """Record the viralconseq version that wrote this config.

        The workflow copies it into ``versions.tsv``; nothing reads it as a
        parameter, so a hand-written config may omit it.
        """
        self._set(ConfigKeys.VIRALCONSEQ_VERSION, viralconseq_version, self.SECTION_PROVENANCE)

    def add_workflow_path(self, workflow_path: str) -> None:
        """Add workflow path to configuration.

        Args:
            workflow_path: Path to the workflow directory
        """
        self._set("workflow_path", workflow_path, self.SECTION_PARAMETERS)

    def add_resource_settings(self, args: Dict[str, Any], rule_names: list) -> None:
        """Add per-rule resource settings to the configuration.

        ``<rule>_cpus`` is written only when the operator set ``--<rule>-cpus``
        (the rule otherwise falls back to the ``threads`` baseline, see
        ``cpus()`` in ``rules/common.smk``). ``<rule>_ram`` (GB) is written for
        every rule in ``ResourceDefaults.MEMORY_RULES`` that appears in
        ``rule_names``, from ``--<rule>-ram`` or the measured default: those
        rules read it as a required key.

        Args:
            args: Dictionary of pipeline arguments (from the CLI).
            rule_names: Snakemake rule names that have resource options.
        """
        R = self.SECTION_RESOURCES
        for rule in rule_names:
            cpus = args.get(f"{rule}_cpus")
            if cpus is not None:
                self._set(f"{rule}_cpus", int(cpus), R)
            if rule in ResourceDefaults.MEMORY_RULES:
                ram = args.get(f"{rule}_ram")
                self._set(
                    f"{rule}_ram",
                    int(ram) if ram is not None else ResourceDefaults.ram_for(rule),
                    R,
                )

    def add_run_resources(
        self, threads_total: int, max_memory_mb: int, memory_detected_mb: int
    ) -> None:
        """Record the run-wide budgets Snakemake was given.

        The rules do not read these; ``viralconseq rerun`` does, and a reader of
        the config sees what bounded the run.
        """
        R = self.SECTION_RESOURCES
        self._set(ConfigKeys.THREADS_TOTAL, int(threads_total), R)
        self._set(ConfigKeys.MAX_MEMORY_MB, int(max_memory_mb), R)
        self._set(ConfigKeys.MEMORY_DETECTED_MB, int(memory_detected_mb), R)

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
                    af_threshold=0.6,
                    chunk_size=10000,
                    clair3_model=_SKELETON_CLAIR3_MODEL,
                    variant_quality=15,
                    variant_depth=10,
                    minimum_map_quality=30,
                    clair3_model_dir=f"{root}/clair3_models",
                )
            # run_viralqc=True pulls run_viralqc (envs/viralqc.yaml) into the
            # DAG; the placeholder DB files are listed in SKELETON_PLACEHOLDERS.
            gen.add_viralqc_settings(run_viralqc=True, viralqc_db=f"{root}/viralqc_db")
            gen.add_collect_settings()
            gen.add_report_settings(run_report=True)
            gen.add_workflow_path(".")
            gen.add_provenance(__version__)
            rules = (
                ResourceDefaults.CONSENSUS_ILLUMINA_RULES
                if data_type == DataType.ILLUMINA
                else ResourceDefaults.CONSENSUS_NANOPORE_RULES
            )
            gen.add_resource_settings({}, rules)
            gen.add_run_resources(threads_total=1, max_memory_mb=0, memory_detected_mb=0)
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
                *_CLAIR3_MODEL_PLACEHOLDERS,
            ],
        },
    }

    def save(self, backup: bool = False) -> None:
        """Write the YAML atomically, one commented section at a time.

        Keys are grouped by the section tags assigned in ``_set()`` and written
        in ``SECTIONS`` order; empty sections are skipped and insertion order is
        kept within a section. The file is written to a temporary sibling and
        renamed into place, so a crash mid-write never leaves a half config.

        Args:
            backup: Keep the existing file as ``<path>.bak`` before replacing it
                (``viralconseq rerun --set`` uses this so an edit is reversible).

        Raises:
            ConfigurationError: If the file or its directory cannot be written.
        """
        config_dir = os.path.dirname(self.config_path)
        if config_dir:  # Only create directory if path contains a directory component
            os.makedirs(config_dir, exist_ok=True)
        if backup and os.path.isfile(self.config_path):
            try:
                shutil.copy2(self.config_path, self.config_path + ".bak")
            except OSError as e:
                raise ConfigurationError(
                    f"Failed to back up {self.config_path} to {self.config_path}.bak: {e}"
                ) from e

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
            # mkstemp creates the staging file 0600 and os.replace preserves
            # that mode. The config is the run's contract (rerun, snakemake -s)
            # and lives in a results tree that is often group-readable, so
            # restore the mode a plain open() would have produced.
            umask = os.umask(0o077)
            os.umask(umask)
            os.chmod(tmp_path, 0o666 & ~umask)
            os.replace(tmp_path, self.config_path)
        except (OSError, IOError) as e:
            if tmp_path is not None:
                with contextlib.suppress(OSError):
                    os.unlink(tmp_path)
            raise ConfigurationError(
                f"Failed to write config file to {self.config_path}: {e}"
            ) from e
