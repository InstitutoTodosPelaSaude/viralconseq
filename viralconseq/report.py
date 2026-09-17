"""Rebuild ``report.html`` from a finished run directory (``viralconseq create-report``).

The workflow's ``rule report`` (``scripts/rules/report.smk``) and this module
drive the same stdlib script, ``scripts/python/build_report.py``, so a page
regenerated after the fact carries exactly what the run would have produced.
The parameters come from ``<run>/config.yml`` when it is present; a run
directory without one still yields a page from ``summary.tsv`` and the
artefacts next to it, with defaults from :class:`constants.ReportDefaults`
and the data type inferred from the layout.
"""

import logging
import os
from typing import Any, Dict, List, Optional, Sequence

import yaml

from viralconseq import scripts
from viralconseq.constants import ConfigKeys, DataType, ReportDefaults
from viralconseq.exceptions import ConfigurationError, ReportError
from viralconseq.scripts.python import build_report

logger = logging.getLogger(__name__)

#: The page shell shipped with the package, the same file ``rule report`` reads.
TEMPLATE_PATH = os.path.join(os.path.dirname(scripts.__file__), "templates", "report.html")
SUMMARY_FILENAME = "summary.tsv"
CONFIG_FILENAME = "config.yml"
VERSIONS_FILENAME = "versions.tsv"
REPORT_FILENAME = "report.html"


def load_run_config(run_dir: str) -> Dict[str, Any]:
    """The mapping in ``<run>/config.yml``, ``{}`` when the file is absent.

    Raises:
        ConfigurationError: If the file exists but is unreadable or not a mapping.
    """
    path = os.path.join(run_dir, CONFIG_FILENAME)
    if not os.path.isfile(path):
        return {}
    try:
        with open(path) as handle:
            loaded = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as e:
        raise ConfigurationError(f"Cannot read {path}: {e}") from e
    if loaded is None:
        return {}
    if not isinstance(loaded, dict):
        raise ConfigurationError(f"{path} does not hold a mapping at the top level")
    return loaded


def validate_run_dir(run_dir: str) -> List[str]:
    """Check that ``run_dir`` is a finished run; return the warnings worth printing.

    Raises:
        ReportError: If the directory or its ``summary.tsv`` is missing.
    """
    if not os.path.isdir(run_dir):
        raise ReportError(f"{run_dir} is not a directory")
    if not os.path.isfile(os.path.join(run_dir, SUMMARY_FILENAME)):
        raise ReportError(
            f"{run_dir} has no {SUMMARY_FILENAME}: is this a finished viralconseq run directory "
            "(the <output>/<run_name> folder)?"
        )
    warnings = []
    for name, consequence in (
        (CONFIG_FILENAME, "parameters fall back to defaults and the data type is inferred"),
        (VERSIONS_FILENAME, "the page carries no tool versions"),
    ):
        if not os.path.isfile(os.path.join(run_dir, name)):
            warnings.append(f"{run_dir}/{name} not found: {consequence}")
    return warnings


def infer_data_type(run_dir: str) -> str:
    """Illumina or nanopore from the run layout: Clair3 output means nanopore."""
    if os.path.isdir(os.path.join(run_dir, "assembly", "clair3")):
        return DataType.NANOPORE
    if os.path.isdir(os.path.join(run_dir, "qc", "reports")):
        return DataType.ILLUMINA
    for entry in (
        os.listdir(os.path.join(run_dir, "assembly"))
        if os.path.isdir(os.path.join(run_dir, "assembly"))
        else []
    ):
        if os.path.isdir(os.path.join(run_dir, "assembly", entry, "clair3")):
            return DataType.NANOPORE
    return ReportDefaults.DATA_TYPE


def report_argv(
    run_dir: str,
    config: Dict[str, Any],
    output: str,
    label: Optional[str] = None,
    template: str = TEMPLATE_PATH,
) -> List[str]:
    """The ``build_report.py`` argument list for a run, mirroring ``rule report``."""
    data_type = str(config.get(ConfigKeys.DATA) or infer_data_type(run_dir))
    argv = [
        "--run-dir",
        run_dir,
        "--template",
        template,
        "--output",
        output,
        "--data-type",
        data_type,
        "--min-depth",
        str(int(config.get(ConfigKeys.MINIMUM_DEPTH, ReportDefaults.MIN_DEPTH))),
        "--consensus-coverage-threshold",
        str(
            float(
                config.get(
                    ConfigKeys.CONSENSUS_COVERAGE_THRESHOLD,
                    ReportDefaults.CONSENSUS_COVERAGE_THRESHOLD,
                )
            )
        ),
        "--label",
        label or os.path.basename(os.path.normpath(run_dir)),
    ]
    reference = config.get(ConfigKeys.REFERENCE)
    if isinstance(reference, dict) and reference:
        argv += ["--segments", *[str(segment) for segment in reference]]
    scheme = config.get(ConfigKeys.SCHEME)
    if isinstance(scheme, str) and scheme != "NA":
        if os.path.isfile(scheme):
            argv += ["--scheme", scheme]
        else:
            logger.warning("primer scheme %s not found; primer extents are left out", scheme)
    for key, flag in (
        (ConfigKeys.AF_THRESHOLD, "--af-threshold"),
        (ConfigKeys.MINIMUM_LENGTH, "--minimum-length"),
        (ConfigKeys.MINIMUM_MAP_QUALITY, "--minimum-map-quality"),
    ):
        if key in config and config[key] is not None:
            argv += [flag, str(config[key])]
    model = config.get(ConfigKeys.CLAIR3_MODEL)
    if isinstance(model, str):
        argv += ["--clair3-model", model]
    return argv


def create_report(run_dir: str, output: Optional[str] = None, label: Optional[str] = None) -> str:
    """Build the report for ``run_dir``; return the path written.

    Raises:
        ReportError: If the directory is not a finished run or the build fails
            (including the page failing its self-check against ``summary.tsv``).
        ConfigurationError: If ``<run>/config.yml`` is present but unusable.
    """
    run_dir = os.path.abspath(run_dir)
    for warning in validate_run_dir(run_dir):
        logger.warning(warning)
    config = load_run_config(run_dir)
    output = os.path.abspath(output or os.path.join(run_dir, REPORT_FILENAME))
    argv: Sequence[str] = report_argv(run_dir, config, output, label)
    try:
        code = build_report.main(list(argv))
    except SystemExit as e:  # build_report reports its failures this way
        raise ReportError(f"report build failed: {e.code}") from e
    if code not in (0, None):
        raise ReportError(f"report build failed with exit status {code}")
    return output
