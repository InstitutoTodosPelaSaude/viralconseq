"""Empirical-data integration tests.

These tests download a real FASTQ dataset, fetch reference files, then run
each pipeline declared by every scenario manifest under
``test/empirical/scenarios/``. They assert (a) exit code 0 and (b) the presence
of a small set of headline output files. No byte-level comparisons.

Adding a new dataset / organism / pipeline flavor is a manifest change only
(see ``test/empirical/scenarios/sars_cov_2.yaml`` for the schema). The pytest
module itself does not encode anything pipeline-specific.

Prerequisite: ``viralconseq`` is on ``PATH`` and the per-rule conda envs have
been materialized once via::

    viralconseq setup --pipelines all

(or repeat ``--pipelines <name>`` for consensus-illumina and
consensus-nanopore.)

Otherwise the first invocation pays several minutes building envs.

These tests are gated behind ``@pytest.mark.empirical`` and deselected by
``pyproject.toml`` default addopts. Run them explicitly with:

    make test-empirical
"""

from __future__ import annotations

import pytest
from empirical.fixtures import ensure_references, ensure_tarball, materialize_workdir
from empirical.runner import Pipeline, Scenario, load_scenarios, run_pipeline

pytestmark = pytest.mark.empirical


def _cases() -> list[tuple[Scenario, Pipeline]]:
    return [(s, p) for s in load_scenarios() for p in s.pipelines]


def _case_id(case: tuple[Scenario, Pipeline]) -> str:
    scenario, pipeline = case
    return f"{scenario.name}::{pipeline.id}"


@pytest.mark.parametrize("case", _cases(), ids=_case_id)
def test_pipeline(case: tuple[Scenario, Pipeline], tmp_path) -> None:
    scenario, pipeline = case
    ensure_tarball(scenario)
    ensure_references(scenario)
    workdir = materialize_workdir(scenario, pipeline, tmp_path)
    proc = run_pipeline(scenario, pipeline, workdir)
    assert proc.returncode == 0, (
        f"{pipeline.id} exited {proc.returncode}\n"
        f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    )
    out_root = workdir / "results" / pipeline.id / pipeline.run_name
    missing = [rel for rel in pipeline.expected_outputs if not (out_root / rel).exists()]
    assert not missing, f"missing expected outputs under {out_root}: {missing}"
