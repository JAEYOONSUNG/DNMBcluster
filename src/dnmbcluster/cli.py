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
@click.option(
    "--max-ram", type=str, default=None,
    help="Memory cap, e.g. '8G'. Propagated to clustering tools.",
)
@click.option(
    "--parse-only", is_flag=True, default=False,
    help="Stop after GenBank parsing. Useful for debugging.",
)
def run(
    input_dir: Path,
    output: Path,
    tool: str,
    level: str,
    identity: float | None,
    coverage: float,
    threads: int,
    max_ram: str | None,
    parse_only: bool,
) -> None:
    """Run the pipeline on a folder of GenBank files."""
    # Defer heavy imports so `--help` stays fast.
    import pyarrow.parquet as pq

    from .engines import ClusterParams, EngineError, get_engine
    from .fasta import write_fasta
    from .genbank import build_tables, parse_folder

    if identity is None:
        identity = 0.5 if level == "protein" else 0.7

    max_ram_gb = _parse_ram(max_ram) if max_ram else None

    output.mkdir(parents=True, exist_ok=True)
    dnmb_dir = output / "dnmb"
    dnmb_dir.mkdir(parents=True, exist_ok=True)

    # Validate engine support for chosen level *before* parsing so we
    # fail fast on impossible combinations (SPEED.md Section 12).
    if not parse_only:
        try:
            engine = get_engine(tool)
        except EngineError as exc:
            raise click.ClickException(str(exc)) from exc
        if not engine.supports(level):  # type: ignore[arg-type]
            raise click.ClickException(
                f"Engine {tool!r} does not support --level {level!r}. "
                f"Supported: {engine.supported_levels}"
            )
        try:
            engine.check_available()
        except EngineError as exc:
            raise click.ClickException(str(exc)) from exc

    # ---------- Stage 1: GenBank parse ----------
    click.secho(f"[dnmbcluster] parsing GenBank files from {input_dir}", fg="cyan")
    parsed = parse_folder(input_dir, level=level, threads=threads)
    click.echo(f"[dnmbcluster]   {len(parsed)} genomes parsed")

    id_map, gene_table, genome_meta = build_tables(parsed, level)

    zstd_kwargs = {"compression": "zstd", "compression_level": 3}
    id_map_path = dnmb_dir / "id_map.parquet"
    gene_table_path = dnmb_dir / "gene_table.parquet"
    genome_meta_path = dnmb_dir / "genome_meta.parquet"
    pq.write_table(id_map, id_map_path, **zstd_kwargs)
    pq.write_table(gene_table, gene_table_path, **zstd_kwargs)
    pq.write_table(genome_meta, genome_meta_path, **zstd_kwargs)

    click.secho(
        f"[dnmbcluster]   id_map.parquet       {id_map.num_rows:>8} rows", fg="green",
    )
    click.secho(
        f"[dnmbcluster]   gene_table.parquet   {gene_table.num_rows:>8} rows", fg="green",
    )
    click.secho(
        f"[dnmbcluster]   genome_meta.parquet  {genome_meta.num_rows:>8} rows", fg="green",
    )

    if parse_only:
        click.secho("[dnmbcluster] --parse-only: stopping.", fg="yellow")
        return

    # ---------- Stage 2: FASTA export ----------
    fasta_path = dnmb_dir / ("proteins.faa" if level == "protein" else "cds.fna")
    n_seq, n_bytes = write_fasta(gene_table_path, fasta_path, level=level)
    click.secho(
        f"[dnmbcluster]   {fasta_path.name}   {n_seq} sequences, "
        f"{n_bytes / (1024 * 1024):.1f} MiB",
        fg="green",
    )

    # ---------- Stage 3: Clustering ----------
    click.secho(
        f"[dnmbcluster] clustering with {tool} "
        f"(id={identity}, cov={coverage}, level={level})",
        fg="cyan",
    )
    params = ClusterParams(
        identity=identity,
        coverage=coverage,
        threads=threads,
        max_ram_gb=max_ram_gb,
        level=level,  # type: ignore[arg-type]
    )
    try:
        result = engine.cluster(fasta_path, dnmb_dir, params)
    except EngineError as exc:
        raise click.ClickException(str(exc)) from exc

    click.secho(
        f"[dnmbcluster]   clusters.parquet     "
        f"{result.n_clusters:>8} clusters / {result.n_input_sequences} sequences",
        fg="green",
    )
    click.secho(
        "[dnmbcluster] M2 complete. Pan/core (M4) + plots (M5) not yet wired.",
        fg="yellow",
    )


@main.command("list-engines")
def list_engines() -> None:
    """List available clustering engines."""
    for engine in ENGINES:
        click.echo(engine)


def _parse_ram(value: str) -> float:
    """Parse a memory string like '8G', '512M' into GiB as a float."""
    s = value.strip().upper()
    if s.endswith("G"):
        return float(s[:-1])
    if s.endswith("M"):
        return float(s[:-1]) / 1024.0
    if s.endswith("K"):
        return float(s[:-1]) / (1024.0 * 1024.0)
    return float(s) / (1024.0 ** 3)


if __name__ == "__main__":
    main()
