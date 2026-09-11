"""Click CLI for ``viralconseq setup``.

Pre-builds the per-rule conda environments declared in the Snakemake
workflows under ``viralconseq/scripts/`` into a shared cache directory so
subsequent ``viralconseq consensus`` runs reuse them
instead of materializing fresh envs per working directory.

Why this exists: dynamic env creation at run-time is brittle (an upstream
conda 26.x / bioconda interaction caused ``CreateCondaEnvironmentException``
on first-run installs). Pre-warming the cache once isolates the env-build
failure surface from real pipeline runs.
"""

import logging
import os
import re
import tempfile
from pathlib import Path
from typing import List, Tuple

import click
from snakemake import snakemake

from viralconseq.config_generator import ConfigGenerator

logger = logging.getLogger(__name__)


# (pipeline, data_type) -> ".smk" filename under viralconseq/scripts/.
# Segmented variants are intentionally omitted: they share the same env
# YAMLs as their non-segmented counterparts, so adding them would only
# duplicate work.
_PIPELINE_TO_WORKFLOW = {
    "consensus-illumina": ("consensus", "illumina", "consensus_illumina.smk"),
    "consensus-nanopore": ("consensus", "nanopore", "consensus_nanopore.smk"),
}
_ALL_PIPELINES = list(_PIPELINE_TO_WORKFLOW.keys())


def _default_conda_prefix() -> str:
    """Resolve the default conda-prefix at CLI invocation time.

    Mirrors the helper in the consensus CLI so both commands share the same
    cache by default.
    """
    return os.environ.get(
        "VIRALCONSEQ_CONDA_PREFIX",
        str(Path.home() / ".cache" / "viralconseq" / "conda-envs"),
    )


def _scripts_dir() -> Path:
    """Locate ``viralconseq/scripts/`` (installed package or editable checkout)."""
    return Path(__file__).resolve().parent / "scripts"


def _collect_env_yamls(workflow_path: Path) -> List[str]:
    """Return the set of ``envs/<name>.yaml`` strings referenced by the
    top-level .smk *and* every rule module it includes.

    Used by ``--dry-run`` to print the list of envs without invoking
    Snakemake. Snakemake itself resolves these the same way.
    """
    yamls: set = set()
    # Match the env YAML name regardless of relative prefix: top-level .smk
    # uses ``"envs/foo.yaml"`` while rule modules under ``rules/`` use
    # ``"../envs/foo.yaml"``.
    pat = re.compile(r'envs/([^"\s/]+\.yaml)')
    queue = [workflow_path]
    seen: set = set()
    while queue:
        smk = queue.pop()
        if smk in seen or not smk.exists():
            continue
        seen.add(smk)
        text = smk.read_text()
        yamls.update(pat.findall(text))
        for inc in re.findall(r'include:\s*"([^"]+\.smk)"', text):
            queue.append((smk.parent / inc).resolve())
    return sorted(yamls)


def _populate_placeholders(root: Path, pipeline: str, data_type: str) -> None:
    """Touch empty input files so Snakemake's DAG build can resolve them.

    Even with ``conda_create_envs_only=True`` Snakemake walks the DAG and
    checks that rule inputs exist before short-circuiting to env creation.
    The set of files is declared in
    ``ConfigGenerator.SKELETON_PLACEHOLDERS`` so the file list and the
    skeleton config stay in lock-step.
    """
    for rel in ConfigGenerator.SKELETON_PLACEHOLDERS[pipeline][data_type]:
        target = root / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.touch()


def _expand_pipelines(selected: Tuple[str, ...]) -> List[str]:
    """Resolve ``--pipelines`` selections, expanding ``all``."""
    if not selected or "all" in selected:
        return _ALL_PIPELINES
    out = []
    for p in selected:
        if p not in _PIPELINE_TO_WORKFLOW:
            raise click.BadParameter(f"Unknown pipeline: {p}")
        if p not in out:
            out.append(p)
    return out


@click.command(name="setup")
@click.option(
    "--conda-prefix",
    default=_default_conda_prefix,
    show_default="$VIRALCONSEQ_CONDA_PREFIX or ~/.cache/viralconseq/conda-envs",
    help="Directory where per-rule conda envs are cached. Reused by every "
    "subsequent 'viralconseq consensus' run that points "
    "at the same prefix (or that picks it up from $VIRALCONSEQ_CONDA_PREFIX).",
)
@click.option(
    "--pipelines",
    type=click.Choice(_ALL_PIPELINES + ["all"]),
    multiple=True,
    default=("all",),
    show_default=True,
    help="Which pipelines to materialize envs for. Repeatable; defaults to all.",
)
@click.option(
    "--threads",
    default=4,
    show_default=True,
    type=int,
    help="Cores given to Snakemake while materializing envs.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Print the envs that would be created and exit without running conda.",
)
def setup(
    conda_prefix: str,
    pipelines: Tuple[str, ...],
    threads: int,
    dry_run: bool,
) -> None:
    """Pre-build per-rule conda envs into a shared cache.

    Run this once after ``pip install -e .`` (or after upgrading
    viralconseq) to materialize every per-rule conda env without needing
    real input data. Subsequent pipeline runs that point at the same
    ``--conda-prefix`` will reuse the cache and skip env creation
    entirely.
    """
    prefix = Path(conda_prefix).expanduser()
    selected = _expand_pipelines(pipelines)
    scripts_dir = _scripts_dir()

    click.echo(f"Conda prefix: {prefix}")
    click.echo(f"Pipelines:    {', '.join(selected)}")

    if dry_run:
        for name in selected:
            _, _, smk_name = _PIPELINE_TO_WORKFLOW[name]
            smk_path = scripts_dir / smk_name
            yamls = _collect_env_yamls(smk_path)
            click.echo(f"\n[{name}] would create {len(yamls)} envs:")
            for y in yamls:
                click.echo(f"  - {y}")
        return

    prefix.mkdir(parents=True, exist_ok=True)
    if not os.access(prefix, os.W_OK):
        raise click.ClickException(f"Conda prefix not writable: {prefix}")

    failures: List[str] = []
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        placeholder_dir = tmpdir / "placeholders"
        for name in selected:
            pipeline, data_type, smk_name = _PIPELINE_TO_WORKFLOW[name]
            smk_path = scripts_dir / smk_name
            if not smk_path.is_file():
                failures.append(name)
                click.echo(f"\n[{name}] workflow not found: {smk_path}", err=True)
                continue
            _populate_placeholders(placeholder_dir, pipeline, data_type)
            skel_path = tmpdir / f"{name}.yaml"
            ConfigGenerator.write_skeleton(
                pipeline, data_type, str(skel_path), str(placeholder_dir)
            )
            click.echo(f"\n[{name}] materializing envs into {prefix} ...")
            ok = snakemake(
                str(smk_path),
                configfiles=[str(skel_path)],
                cores=threads,
                use_conda=True,
                conda_prefix=str(prefix),
                conda_create_envs_only=True,
                targets=["all"],
            )
            if ok:
                click.echo(f"[{name}] OK")
            else:
                failures.append(name)
                click.echo(f"[{name}] FAILED", err=True)

    if failures:
        raise click.ClickException("Env creation failed for: " + ", ".join(failures))
    click.echo("\nAll envs ready. Pipeline runs against this prefix will skip env creation.")
