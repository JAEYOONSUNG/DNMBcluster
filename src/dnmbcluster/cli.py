"""DNMBcluster command-line interface."""
from __future__ import annotations

from pathlib import Path

import click

from . import __version__

ENGINES = ["mmseqs2", "diamond", "cd-hit", "usearch12"]
LEVELS = ["protein", "nucleotide"]


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="dnmbcluster")
def main() -> None:
    """DNMBcluster — pan-genome clustering and visualization pipeline.

    Run on a folder of GenBank files to get clustering, pan/core
    analysis, and publication-ready plots in one command.
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
    "--level",
    type=click.Choice(LEVELS),
    default="protein",
    show_default=True,
    help=(
        "Sequence level for comparative genomics. "
        "'protein' (default) clusters CDS translations; "
        "'nucleotide' clusters CDS DNA sequences (slower, more memory)."
    ),
)
@click.option(
    "--identity",
    type=float,
    default=None,
    help="Sequence identity threshold. Defaults: 0.5 for protein, 0.7 for nucleotide.",
)
@click.option(
    "--coverage", type=float, default=0.8, show_default=True,
    help="Alignment coverage threshold.",
)
@click.option(
    "--threads", type=int, default=0, show_default=True,
    help="Threads. 0 = auto-detect.",
)
def run(
    input_dir: Path,
    output: Path,
    tool: str,
    level: str,
    identity: float | None,
    coverage: float,
    threads: int,
) -> None:
    """Run the pipeline on a folder of GenBank files.

    Current implementation: M1 (GenBank parse) only. Clustering,
    presence/absence, and plots land in M2–M5.
    """
    # Defer heavy imports so `--help` stays fast and doesn't require
    # Biopython / Polars / PyArrow to be installed.
    import pyarrow.parquet as pq

    from .genbank import build_tables, parse_folder

    if identity is None:
        identity = 0.5 if level == "protein" else 0.7

    output.mkdir(parents=True, exist_ok=True)
    dnmb_dir = output / "dnmb"
    dnmb_dir.mkdir(parents=True, exist_ok=True)

    click.secho(f"[dnmbcluster] parsing GenBank files from {input_dir}", fg="cyan")
    parsed = parse_folder(input_dir, level=level, threads=threads)
    click.echo(f"[dnmbcluster]   {len(parsed)} genomes parsed")

    id_map, gene_table, genome_meta = build_tables(parsed, level)

    zstd_kwargs = {"compression": "zstd", "compression_level": 3}
    pq.write_table(id_map, dnmb_dir / "id_map.parquet", **zstd_kwargs)
    pq.write_table(gene_table, dnmb_dir / "gene_table.parquet", **zstd_kwargs)
    pq.write_table(genome_meta, dnmb_dir / "genome_meta.parquet", **zstd_kwargs)

    click.secho(
        f"[dnmbcluster]   id_map.parquet       {id_map.num_rows:>8} rows",
        fg="green",
    )
    click.secho(
        f"[dnmbcluster]   gene_table.parquet   {gene_table.num_rows:>8} rows",
        fg="green",
    )
    click.secho(
        f"[dnmbcluster]   genome_meta.parquet  {genome_meta.num_rows:>8} rows",
        fg="green",
    )
    click.secho(
        f"[dnmbcluster] M1 complete. Clustering ({tool}, id={identity}, "
        f"cov={coverage}, level={level}) not yet implemented — M2+.",
        fg="yellow",
    )


@main.command("list-engines")
def list_engines() -> None:
    """List available clustering engines."""
    for engine in ENGINES:
        click.echo(engine)


if __name__ == "__main__":
    main()
