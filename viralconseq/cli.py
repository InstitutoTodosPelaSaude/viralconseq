"""Top-level click CLI for viralconseq."""

import click

from viralconseq import __program__, __version__
from viralconseq.consensus_cli import consensus
from viralconseq.create_samplesheet import create_samplesheet
from viralconseq.logging_config import configure_logging
from viralconseq.report_cli import create_report_command
from viralconseq.rerun_cli import rerun
from viralconseq.setup_cli import setup


@click.group()
@click.version_option(version=__version__, prog_name=__program__)
@click.option(
    "--log-level",
    default="INFO",
    show_default=True,
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    help="Logging verbosity.",
)
@click.option(
    "--json-logs",
    is_flag=True,
    default=False,
    help="Emit structured JSON logs (one object per line) instead of text.",
)
def cli(log_level: str, json_logs: bool) -> None:
    """viralconseq infers consensus sequences from viral high-throughput sequencing data.

    \b
    Subcommands:
    * consensus           reference-guided consensus assembly (illumina/nanopore)
    * rerun               run again from a saved config (resume, --dry-run, --unlock, --set)
    * create-report       rebuild report.html for a finished run directory
    * setup               pre-build per-rule conda envs into a shared cache
    * create-samplesheet  generate a sample sheet from a sequencing directory

    Run ``viralconseq <subcommand> --help`` for subcommand-specific options.
    """
    configure_logging(level=log_level, json_logs=json_logs)


cli.add_command(consensus)
cli.add_command(rerun)
cli.add_command(create_report_command)
cli.add_command(setup)
cli.add_command(create_samplesheet)


def main() -> None:
    cli()


if __name__ == "__main__":
    main()
