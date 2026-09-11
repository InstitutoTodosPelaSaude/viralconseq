"""Scenario manifest loading and pipeline invocation for empirical tests.

A "scenario" is a YAML manifest under ``test/empirical/scenarios/`` that fully
describes one dataset and every pipeline that should be exercised against it.
Adding a new organism, segmented genome, or new pipeline flavor is a manifest
change — no new test code is needed.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

SCENARIOS_DIR = Path(__file__).resolve().parent / "scenarios"


@dataclass
class TarballSpec:
    url: str
    extract_to: str
    expected_files: list[str]


@dataclass
class ReferenceSpec:
    url: str
    dest: str


@dataclass
class SampleSheet:
    path: str
    rows: list[list[str]]


@dataclass
class Pipeline:
    id: str
    cli: list[str]
    samples_ref: str
    run_name: str
    extra_args: dict[str, str] = field(default_factory=dict)
    extra_flags: list[str] = field(default_factory=list)
    segmented_reference: dict[str, str] = field(default_factory=dict)
    expected_outputs: list[str] = field(default_factory=list)


@dataclass
class Scenario:
    name: str
    description: str
    data: TarballSpec
    references: list[ReferenceSpec]
    sample_sheets: dict[str, SampleSheet]
    pipelines: list[Pipeline]


def _parse_scenario(raw: dict[str, Any]) -> Scenario:
    tarball = raw["data"]["tarball"]
    data = TarballSpec(
        url=tarball["url"],
        extract_to=tarball["extract_to"],
        expected_files=list(tarball["expected_files"]),
    )
    references = [ReferenceSpec(url=r["url"], dest=r["dest"]) for r in raw.get("references", [])]
    sample_sheets = {
        kind: SampleSheet(path=spec["path"], rows=[list(row) for row in spec["rows"]])
        for kind, spec in raw.get("sample_sheets", {}).items()
    }
    pipelines = [
        Pipeline(
            id=p["id"],
            cli=list(p["cli"]),
            samples_ref=p["samples_ref"],
            run_name=p["run_name"],
            extra_args={k: str(v) for k, v in p.get("extra_args", {}).items()},
            extra_flags=list(p.get("extra_flags", [])),
            segmented_reference=dict(p.get("segmented_reference", {})),
            expected_outputs=list(p.get("expected_outputs", [])),
        )
        for p in raw["pipelines"]
    ]
    return Scenario(
        name=raw["name"],
        description=raw.get("description", ""),
        data=data,
        references=references,
        sample_sheets=sample_sheets,
        pipelines=pipelines,
    )


def load_scenarios(scenarios_dir: Path = SCENARIOS_DIR) -> list[Scenario]:
    """Read every ``*.yaml`` file under ``scenarios_dir`` and validate it."""
    if not scenarios_dir.exists():
        return []
    scenarios = []
    for path in sorted(scenarios_dir.glob("*.yaml")):
        with path.open() as fh:
            raw = yaml.safe_load(fh)
        scenarios.append(_parse_scenario(raw))
    return scenarios


def build_argv(
    pipeline: Pipeline, sample_sheet_path: str, output_dir: str, config_file: str
) -> list[str]:
    """Build the ``viralconseq`` argv for a pipeline from its manifest entry."""
    argv = list(pipeline.cli)
    argv += ["--sample-sheet", sample_sheet_path]
    argv += ["--config-file", config_file]
    argv += ["--output", output_dir]
    argv += ["--run-name", pipeline.run_name]
    for flag, value in pipeline.extra_args.items():
        argv += [flag, value]
    for segment, path in pipeline.segmented_reference.items():
        argv += ["--segmented-reference", f"{segment}={path}"]
    argv += list(pipeline.extra_flags)
    return argv


def run_pipeline(
    scenario: Scenario, pipeline: Pipeline, workdir: Path
) -> subprocess.CompletedProcess:
    """Invoke a pipeline inside ``workdir``. Returns the CompletedProcess."""
    output_dir = f"results/{pipeline.id}/"
    config_file = f"results/{pipeline.id}/config.yml"
    sample_sheet = scenario.sample_sheets[pipeline.samples_ref].path
    argv = build_argv(pipeline, sample_sheet, output_dir, config_file)
    return subprocess.run(argv, cwd=workdir, capture_output=True, text=True)
