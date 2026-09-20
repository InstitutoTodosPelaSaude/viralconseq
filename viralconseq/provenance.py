"""Run provenance manifest.

For public-health reporting, results must be reproducible-by-record: given an
output, you must be able to recover *which* pipeline version, config, and exact
input files produced it. This module writes a ``run_manifest.json`` into the run
output directory capturing the viralconseq version, a timestamp, the resolved
config path, and a checksum/size for every input FASTQ.

Tool and database versions are captured at rule-execution time, inside the
per-rule conda envs, by ``rules/provenance.smk`` into ``<run>/versions.tsv``;
the manifest records that file's path once the run has finished. This module
covers the orchestration-level provenance the Python layer can record reliably.
"""

from __future__ import annotations

import datetime
import glob
import hashlib
import json
import os
import shutil
from typing import Any, Dict, Optional

from viralconseq import __version__
from viralconseq.config_generator import sample_key
from viralconseq.validators import _is_path_sentinel

MANIFEST_FILENAME = "run_manifest.json"
SNAKEMAKE_LOG_COPY = os.path.join("logs", "snakemake.log")


def _sha256(path: str, chunk_size: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def _describe_input(path: str) -> Dict[str, Any]:
    """Return a provenance record for a single input file."""
    record: Dict[str, Any] = {"path": os.path.abspath(path)}
    if os.path.isfile(path):
        stat = os.stat(path)
        record["size_bytes"] = stat.st_size
        record["sha256"] = _sha256(path)
    else:
        record["missing"] = True
    return record


def build_run_manifest(
    args: Dict[str, Any],
    samples: Dict[str, list],
    *,
    timestamp: Optional[str] = None,
) -> Dict[str, Any]:
    """Build (but do not write) the run-manifest dict.

    Args:
        args: The pipeline argument dict (needs ``output``, ``run_name``,
            ``config_file``, ``data_type``).
        samples: Mapping of sample id -> list of input FASTQ paths.
        timestamp: ISO-8601 timestamp; generated (UTC) if omitted.
    """
    if timestamp is None:
        timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()

    # Keyed ``sample-<id>`` like every other output (config, samples/ dirs,
    # FASTA headers, summary tables).
    sample_inputs = {
        sample_key(sample): [_describe_input(p) for p in paths]
        for sample, paths in (samples or {}).items()
    }

    config_file = args.get("config_file")
    # Hash the resolved config so the manifest pins the exact content used, not
    # just a path that could later be edited.
    config_sha256 = _sha256(config_file) if config_file and os.path.isfile(config_file) else None

    return {
        "viralconseq_version": __version__,
        "created_utc": timestamp,
        "run_name": args.get("run_name"),
        "data_type": args.get("data_type"),
        "config_file": (os.path.abspath(args["config_file"]) if args.get("config_file") else None),
        "config_sha256": config_sha256,
        "output": (
            os.path.abspath(os.path.join(args["output"], args.get("run_name", "")))
            if args.get("output")
            else None
        ),
        "sample_count": len(samples or {}),
        "viralqc": {
            "enabled": bool(args.get("run_viralqc", True)),
            "db": (
                None
                if _is_path_sentinel(args.get("viralqc_db"))
                else os.path.abspath(str(args["viralqc_db"]))
            ),
        },
        "samples": sample_inputs,
    }


def write_run_manifest(
    args: Dict[str, Any],
    samples: Dict[str, list],
    *,
    timestamp: Optional[str] = None,
) -> str:
    """Write the run manifest into ``<output>/<run_name>/run_manifest.json``.

    Returns:
        The path to the written manifest.
    """
    manifest = build_run_manifest(args, samples, timestamp=timestamp)
    run_dir = os.path.join(args["output"], args.get("run_name", ""))
    os.makedirs(run_dir, exist_ok=True)
    manifest_path = os.path.join(run_dir, MANIFEST_FILENAME)
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    return manifest_path


def record_run_completion(
    manifest_path: str,
    *,
    status: str,
    timestamp: Optional[str] = None,
    extra: Optional[Dict[str, Any]] = None,
) -> None:
    """Patch an existing manifest with the run outcome.

    The manifest is written *before* the workflow runs (it records intent), so
    this records the actual result: ``status`` ("success"/"failed") and a finish
    timestamp. Best-effort by contract — callers should not let a failure here
    abort the analysis.

    Args:
        manifest_path: Path returned by :func:`write_run_manifest`.
        status: Outcome string, e.g. "success" or "failed".
        timestamp: ISO-8601 finish time; generated (UTC) if omitted.
        extra: Additional top-level keys to record (e.g. the path of the
            copied Snakemake log).
    """
    if timestamp is None:
        timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with open(manifest_path) as fh:
        manifest = json.load(fh)
    manifest["status"] = status
    manifest["finished_utc"] = timestamp
    for key, value in (extra or {}).items():
        manifest[key] = value
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)


def copy_snakemake_log(run_dir: str) -> Optional[str]:
    """Copy the newest Snakemake transcript into ``<run_dir>/logs/snakemake.log``.

    Snakemake writes its console transcript under ``<run_dir>/.snakemake/log/``
    with a timestamped name. Keeping a stable copy next to the results makes the
    run self-describing and survives a later ``.snakemake/`` clean-up.

    Returns:
        The path of the copy, or ``None`` when no transcript exists. Never
        raises: this is provenance, not analysis.
    """
    try:
        candidates = glob.glob(os.path.join(run_dir, ".snakemake", "log", "*.snakemake.log"))
        if not candidates:
            return None
        newest = max(candidates, key=os.path.getmtime)
        destination = os.path.join(run_dir, SNAKEMAKE_LOG_COPY)
        os.makedirs(os.path.dirname(destination), exist_ok=True)
        shutil.copyfile(newest, destination)
        return destination
    except OSError:
        return None
