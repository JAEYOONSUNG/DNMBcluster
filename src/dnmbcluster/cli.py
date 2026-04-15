"""DNMBcluster command-line interface."""
from __future__ import annotations

import sys
from pathlib import Path

import click

from . import __version__

ENGINES = ["mmseqs2", "diamond", "cd-hit", "usearch12"]


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="dnmbcluster")
def main() -> None:
    """DNMBcluster — pan-genome clustering and visualization pipeline.

    Run on a folder of GenBank files to get clustering, pan/core analysis,
    and publication-ready plots in one command.
    """


@main.command()
@click.argument(
    "input_dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--output", "-o",
    type=click.Path(path_type=Path),
    default=Path("results"),
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--tool",
    type=click.Choice(ENGINES),
    default="mmseqs2",
    show_default=True,
    help="Clustering engine. Default TBD after Geobacillus-10 benchmark (M6).",
)
@click.option(
    "--identity", type=float, default=0.5, show_default=True,
    help="Sequence identity threshold for clustering.",
)
@click.option(
    "--coverage", type=float, default=0.8, show_default=True,
    help="Alignment coverage threshold.",
)
@click.option(
    "--threads", type=int, default=0, show_default=True,
    help="Threads to use. 0 = auto-detect.",
)
def run(
    input_dir: Path,
    output: Path,
    tool: str,
    identity: float,
    coverage: float,
    threads: int,
) -> None:
    """Run the full pipeline on a folder of GenBank files."""
    click.secho("[dnmbcluster] M0 scaffold — pipeline not yet implemented.", fg="yellow")
    click.echo(f"  input  : {input_dir}")
    click.echo(f"  output : {output}")
    click.echo(f"  tool   : {tool}")
    click.echo(f"  id={identity}  cov={coverage}  threads={threads or 'auto'}")
    sys.exit(0)


@main.command("list-engines")
def list_engines() -> None:
    """List available clustering engines."""
    for engine in ENGINES:
        click.echo(engine)


if __name__ == "__main__":
    main()
