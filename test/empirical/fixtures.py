"""Idempotent download / extraction / workdir helpers for empirical tests.

Every function checks sentinel files before doing work, so the suite can run
repeatedly against a warm cache without re-downloading anything. The cache
lives at ``~/.cache/viralconseq-test-data/<scenario>/`` by default; override
with the ``VIRALCONSEQ_TEST_CACHE`` environment variable.
"""

from __future__ import annotations

import csv
import os
import tarfile
import urllib.request
from pathlib import Path

from .runner import Pipeline, Scenario


def _default_cache_dir() -> Path:
    override = os.environ.get("VIRALCONSEQ_TEST_CACHE")
    if override:
        return Path(override).expanduser().resolve()
    return Path.home() / ".cache" / "viralconseq-test-data"


CACHE_DIR = _default_cache_dir()


def scenario_cache(scenario: Scenario, cache_dir: Path | None = None) -> Path:
    """Per-scenario subdir so multiple scenarios cannot collide on disk."""
    root = cache_dir or CACHE_DIR
    path = root / scenario.name
    path.mkdir(parents=True, exist_ok=True)
    return path


def ensure_tarball(scenario: Scenario, cache_dir: Path | None = None) -> Path:
    """Download + extract ``scenario.data.tarball``; skip if all sentinels exist."""
    base = scenario_cache(scenario, cache_dir)
    if all((base / rel).exists() for rel in scenario.data.expected_files):
        return base
    tarball_path = base / "data.tar.gz"
    if not tarball_path.exists() or tarball_path.stat().st_size == 0:
        urllib.request.urlretrieve(scenario.data.url, tarball_path)
    with tarfile.open(tarball_path) as tf:
        tf.extractall(base)
    tarball_path.unlink(missing_ok=True)
    missing = [rel for rel in scenario.data.expected_files if not (base / rel).exists()]
    if missing:
        raise RuntimeError(f"tarball {scenario.data.url} did not produce expected files: {missing}")
    return base


def ensure_references(scenario: Scenario, cache_dir: Path | None = None) -> Path:
    """Download each plain-URL reference; skip ones whose dest already exists."""
    base = scenario_cache(scenario, cache_dir)
    for ref in scenario.references:
        dest = base / ref.dest
        if dest.exists() and dest.stat().st_size > 0:
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(ref.url, dest)
    return base


def _write_sample_sheet(rows: list[list[str]], dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerows(rows)


def _symlink(src: Path, dst: Path) -> None:
    if dst.exists() or dst.is_symlink():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    os.symlink(src, dst)


def materialize_workdir(
    scenario: Scenario,
    pipeline: Pipeline,
    tmp_path: Path,
    cache_dir: Path | None = None,
) -> Path:
    """Build a per-test workdir: symlinks to cached data/refs plus a sample sheet."""
    base = scenario_cache(scenario, cache_dir)
    workdir = tmp_path
    _symlink(base / scenario.data.extract_to, workdir / scenario.data.extract_to)
    if scenario.references:
        ref_root = base / "references"
        if ref_root.exists():
            _symlink(ref_root, workdir / "references")
    sheet = scenario.sample_sheets[pipeline.samples_ref]
    _write_sample_sheet(sheet.rows, workdir / sheet.path)
    return workdir
