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
from viralconseq.constants import Clair3Models, ViralQCDatabase
from viralconseq.scripts.python import fetch_clair3_model
from viralconseq.validators import missing_clair3_model_files, missing_viralqc_database_files

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
    # ``or`` (not the get() default) so an exported-but-empty
    # VIRALCONSEQ_CONDA_PREFIX falls back instead of resolving to "".
    return os.environ.get("VIRALCONSEQ_CONDA_PREFIX", "") or str(
        Path.home() / ".cache" / "viralconseq" / "conda-envs"
    )


def _expand_models(selected: Tuple[str, ...]) -> List[str]:
    """Resolve ``--clair3-models`` values: ``all`` -> the manifest; names must be
    manifest entries (their download URLs are known only for those)."""
    names: List[str] = []
    for entry in selected:
        for name in entry.split(","):
            name = name.strip()
            if not name:
                continue
            if name == "all":
                names.extend(n for n in Clair3Models.MANIFEST if n not in names)
            elif name in Clair3Models.MANIFEST:
                if name not in names:
                    names.append(name)
            else:
                raise click.BadParameter(
                    f"unknown Clair3 model {name!r}. Known models: "
                    + ", ".join(Clair3Models.MANIFEST),
                    param_hint="--clair3-models",
                )
    return names


def _fetch_clair3_models(model_dir: Path, names: List[str]) -> List[str]:
    """Fetch every model in ``names`` that is absent or incomplete. Returns the
    names that could not be fetched."""
    failures: List[str] = []
    for name in names:
        if not missing_clair3_model_files(str(model_dir), name):
            click.echo(f"[clair3-models] {name}: already present")
            continue
        click.echo(f"[clair3-models] {name}: downloading ...")
        try:
            source = fetch_clair3_model.fetch_model(
                str(model_dir),
                name,
                list(Clair3Models.urls_for(name)),
                log=lambda message: logger.info(message),
            )
        except RuntimeError as exc:
            failures.append(name)
            click.echo(f"[clair3-models] {name}: FAILED - {exc}", err=True)
            continue
        click.echo(f"[clair3-models] {name}: OK (from {source})")
    return failures


_VIRALQC_SETUP_SMK = "viralqc_setup.smk"
_VIRALQC_DB_BANNER = (
    "  nextclade datasets + NCBI RefSeq viral BLAST reference set: about 1 GB on "
    "disk, several GB transient during the download; typically 15-60 min.\n"
    "  taxonkit's copy of the NCBI taxonomy dump is placed under $HOME/.taxonkit."
)


def _run_viralqc_db_download(db_dir: Path, prefix: Path, threads: int, workdir: Path) -> bool:
    """Drive ``scripts/viralqc_setup.smk`` to populate *db_dir*.

    Reuses the same conda prefix as the pipeline envs, so the ``viralqc.yaml``
    env built for the rules is the one used for the download. ``workdir``
    confines Snakemake's ``.snakemake/`` state to a scratch directory (the
    outputs are absolute paths); ``rerun_triggers=["mtime"]`` stops a later
    edit of the env file from silently re-downloading the databases.
    """
    db_dir.mkdir(parents=True, exist_ok=True)
    cwd = os.getcwd()
    try:
        return bool(
            snakemake(
                str(_scripts_dir() / _VIRALQC_SETUP_SMK),
                config={"viralqc_db": str(db_dir)},
                cores=threads,
                use_conda=True,
                conda_prefix=str(prefix),
                workdir=str(workdir),
                rerun_triggers=["mtime"],
                targets=["all"],
            )
        )
    finally:
        os.chdir(cwd)


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
    type=click.IntRange(min=1),
    help="Cores given to Snakemake while materializing envs and downloading databases.",
)
@click.option(
    "--viralqc-db",
    default=ViralQCDatabase.default_dir,
    show_default="$VIRALCONSEQ_VIRALQC_DB or ~/.cache/viralconseq/viralqc-db",
    help="Directory where the viralQC databases (nextclade datasets + BLAST "
    "reference set) are downloaded. 'viralconseq consensus' reads the same "
    "location by default.",
)
@click.option(
    "--skip-viralqc-db",
    is_flag=True,
    default=False,
    help="Only build conda envs; do not download the viralQC databases.",
)
@click.option(
    "--clair3-models",
    multiple=True,
    default=Clair3Models.DEFAULT_SETUP_MODELS,
    show_default=True,
    metavar="NAME",
    help="Clair3 models to download into --clair3-model-dir (repeatable, or comma-"
    "separated; 'all' fetches every model of the manifest). Names are the Clair3 "
    "model zoo names, e.g. r1041_e82_400bps_hac_v500. About 20 MB each.",
)
@click.option(
    "--clair3-model-dir",
    default=Clair3Models.default_dir,
    show_default="$VIRALCONSEQ_CLAIR3_MODELS or ~/.cache/viralconseq/clair3-models",
    help="Directory the models are downloaded into; 'viralconseq consensus nanopore' "
    "reads the same location by default.",
)
@click.option(
    "--skip-clair3-models",
    is_flag=True,
    default=False,
    help="Do not download Clair3 models.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Print the envs, databases and models that would be created and exit without "
    "running conda or downloading anything.",
)
def setup(
    conda_prefix: str,
    pipelines: Tuple[str, ...],
    threads: int,
    viralqc_db: str,
    skip_viralqc_db: bool,
    clair3_models: Tuple[str, ...],
    clair3_model_dir: str,
    skip_clair3_models: bool,
    dry_run: bool,
) -> None:
    """Pre-build per-rule conda envs into a shared cache.

    Run this once after ``pip install -e .`` (or after upgrading
    viralconseq) to materialize every per-rule conda env without needing
    real input data. Subsequent pipeline runs that point at the same
    ``--conda-prefix`` will reuse the cache and skip env creation
    entirely.

    Unless ``--skip-viralqc-db`` is given, it also downloads the viralQC
    databases (about 1 GB on disk, several GB transient, typically 15-60
    min) into ``--viralqc-db``; rerunning is a no-op once they exist.
    Unless ``--skip-clair3-models`` is given, it downloads the Clair3 models
    named by ``--clair3-models`` into ``--clair3-model-dir`` (about 20 MB
    each); a model already present and complete is not fetched again.
    """
    # Absolute paths: the database download runs Snakemake with a scratch
    # ``workdir`` (it chdirs), so relative prefixes/paths would resolve there.
    prefix = Path(conda_prefix).expanduser().absolute()
    db_dir = Path(viralqc_db).expanduser().absolute()
    model_dir = Path(clair3_model_dir).expanduser().absolute()
    selected = _expand_pipelines(pipelines)
    models = [] if skip_clair3_models else _expand_models(clair3_models)
    scripts_dir = _scripts_dir()

    click.echo(f"Conda prefix: {prefix}")
    click.echo(f"Pipelines:    {', '.join(selected)}")
    click.echo(
        f"viralQC DB:   {db_dir}{'  (skipped: --skip-viralqc-db)' if skip_viralqc_db else ''}"
    )
    click.echo(
        f"Clair3 models: {model_dir}"
        + (
            "  (skipped: --skip-clair3-models)"
            if skip_clair3_models
            else f"  ({len(models)} model(s))"
        )
    )

    if dry_run:
        for name in selected:
            _, _, smk_name = _PIPELINE_TO_WORKFLOW[name]
            smk_path = scripts_dir / smk_name
            yamls = _collect_env_yamls(smk_path)
            click.echo(f"\n[{name}] would create {len(yamls)} envs:")
            for y in yamls:
                click.echo(f"  - {y}")
        if not skip_viralqc_db:
            missing = missing_viralqc_database_files(str(db_dir))
            if not missing:
                click.echo(f"\n[viralqc-db] already present at {db_dir}; nothing to download.")
            else:
                click.echo(
                    f"\n[viralqc-db] would download the viralQC databases into {db_dir} "
                    f"(missing: {', '.join(missing)}):"
                )
                click.echo(_VIRALQC_DB_BANNER)
        if not skip_clair3_models:
            click.echo("")
            for name in models:
                missing = missing_clair3_model_files(str(model_dir), name)
                if missing:
                    click.echo(
                        f"[clair3-models] would download {name} into {model_dir} "
                        f"(missing: {', '.join(missing)})"
                    )
                else:
                    click.echo(f"[clair3-models] {name} already present at {model_dir}")
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

        if not skip_viralqc_db:
            if failures:
                click.echo("\n[viralqc-db] skipped: env creation failed above.", err=True)
            elif not missing_viralqc_database_files(str(db_dir)):
                click.echo(f"\n[viralqc-db] already present at {db_dir}; skipping download.")
            else:
                click.echo(f"\n[viralqc-db] downloading the viralQC databases into {db_dir} ...")
                click.echo(_VIRALQC_DB_BANNER)
                db_workdir = tmpdir / "viralqc_setup_workdir"
                db_workdir.mkdir()
                ok = _run_viralqc_db_download(db_dir, prefix, threads, db_workdir)
                still_missing = missing_viralqc_database_files(str(db_dir))
                if ok and not still_missing:
                    click.echo("[viralqc-db] OK")
                else:
                    failures.append("viralqc-db")
                    missing_note = (
                        f" (still missing: {', '.join(still_missing)})" if still_missing else ""
                    )
                    click.echo(
                        f"[viralqc-db] FAILED{missing_note} - see {db_dir}/logs/ and rerun "
                        f"'viralconseq setup --viralqc-db {db_dir}' (a partial "
                        f"{db_dir}/tmp_ncbi/ may remain and can be deleted)",
                        err=True,
                    )

    if not skip_clair3_models:
        click.echo("")
        model_failures = _fetch_clair3_models(model_dir, models)
        if model_failures:
            failures.append("clair3-models")
            click.echo(
                "[clair3-models] FAILED for "
                + ", ".join(model_failures)
                + " - rerun 'viralconseq setup --skip-viralqc-db --clair3-models "
                # Comma-separated: --clair3-models is multiple=True, so a
                # space-separated list is parsed as extra positional arguments
                # and click rejects the command we just told the user to run.
                + ",".join(model_failures)
                + f" --clair3-model-dir {model_dir}' (see the network errors above), or copy "
                "the model directories in from another machine",
                err=True,
            )

    if failures:
        raise click.ClickException("Setup failed for: " + ", ".join(failures))
    click.echo(
        "\nAll envs ready"
        + ("" if skip_viralqc_db else ", viralQC databases in place")
        + ("" if skip_clair3_models else ", Clair3 models in place")
        + ". Pipeline runs against this prefix will skip env creation."
    )
