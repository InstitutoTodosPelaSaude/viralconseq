"""``viralconseq create-report``: rebuild ``report.html`` for a finished run."""

import logging
from typing import Optional

import click

from viralconseq.exceptions import ViralConseqError
from viralconseq.report import create_report

logger = logging.getLogger(__name__)


@click.command(name="create-report")
@click.argument("run_dir", type=click.Path(exists=True, file_okay=False))
@click.option(
    "-o",
    "--output",
    default=None,
    type=click.Path(dir_okay=False),
    help="Where to write the page. [default: RUN_DIR/report.html]",
)
@click.option(
    "--label",
    default=None,
    help="Run name shown in the page header. [default: the run directory's name]",
)
def create_report_command(run_dir: str, output: Optional[str], label: Optional[str]) -> None:
    """Rebuild the interactive report of the run in RUN_DIR (<output>/<run_name>).

    Reads summary.tsv and the artefacts around it, takes the parameters from
    RUN_DIR/config.yml when present, and writes a single self-contained HTML
    page. Useful after --no-report, after editing a template, or to regenerate
    a page for a run produced by an older version.
    """
    try:
        path = create_report(run_dir, output=output, label=label)
    except ViralConseqError as e:
        raise click.ClickException(f"[{e.code}] {e}") from e
    click.echo(f"report written to {path}")
