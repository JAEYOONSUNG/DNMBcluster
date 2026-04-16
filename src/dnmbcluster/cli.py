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
    help=(
        "Sequence identity threshold. Default: 0.6 for both protein and "
        "nucleotide — a balanced mid-point suitable for species-level "
        "pan-genome work across all engines."
    ),
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
@click.option(
    "--alignment/--fast",
    "with_alignment",
    default=True,
    help=(
        "Pipeline-wide alignment enrichment switch — applies to every "
        "engine (mmseqs2, diamond, cdhit, usearch12), not just MMseqs2. "
        "--alignment (default) runs the shared bidirectional realignment "
        "stage after the engine produces cluster membership, populating "
        "pct_identity_fwd/_rev plus member/rep coverage and alignment "
        "length uniformly across engines. --fast skips that stage; "
        "clusters.parquet keeps those columns null. The realign stage "
        "uses MMseqs2 internally regardless of which engine clustered."
    ),
)
@click.option(
    "--columns",
    type=str,
    default="product",
    show_default=True,
    help=(
        "Comma-separated list of per-genome attribute columns inserted "
        "between the fixed locus_tag column and the fixed identity (%) "
        "column in comparative_genomics_N.xlsx. Pick any subset of: "
        "product, gene, protein_id, ec_number, contig, aa_length, "
        "pct_identity_rev, member_coverage, rep_coverage, alignment_length. "
        "Pass an empty string to produce minimal 2-column blocks "
        "(locus_tag + identity). Example: --columns product,gene,protein_id"
    ),
)
@click.option(
    "--phylo/--no-phylo",
    default=False,
    help=(
        "Run core-gene phylogenomics stage (single-copy core → MAFFT "
        "→ IQ-TREE fast mode → ggtree visualization). Off by default; "
        "adds 2–5 minutes on a 10-genome run. Requires mafft + iqtree "
        "on PATH."
    ),
)
@click.option(
    "--ani/--no-ani",
    default=False,
    help=(
        "Compute ANI (skani) and POCP genome-genome similarity "
        "matrices and draw pheatmap heatmaps with hierarchical "
        "clustering. Off by default; adds ~10–30 s on a 10-genome "
        "run. Requires skani on PATH."
    ),
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
    with_alignment: bool,
    columns: str,
    phylo: bool,
    ani: bool,
) -> None:
    """Run the pipeline on a folder of GenBank files."""
    # Defer heavy imports so `--help` stays fast.
    import pyarrow.parquet as pq

    from .engines import ClusterParams, EngineError, get_engine
    from .fasta import write_fasta
    from .genbank import build_tables, parse_folder
    from .r_bridge import RBridgeError  # noqa: F401

    if identity is None:
        identity = 0.6

    max_ram_gb = _parse_ram(max_ram) if max_ram else None

    output.mkdir(parents=True, exist_ok=True)
    dnmb_dir = output / "dnmb"
    inputs_dir = dnmb_dir / "inputs"
    raw_dir_root = dnmb_dir / "raw"
    processed_dir = dnmb_dir / "processed"
    for d in (inputs_dir, raw_dir_root, processed_dir):
        d.mkdir(parents=True, exist_ok=True)

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
    id_map_path = inputs_dir / "id_map.parquet"
    gene_table_path = inputs_dir / "gene_table.parquet"
    genome_meta_path = inputs_dir / "genome_meta.parquet"
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
    fasta_path = inputs_dir / ("proteins.faa" if level == "protein" else "cds.fna")
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
        with_alignment=with_alignment,
    )
    engine_raw_dir = raw_dir_root / engine.name
    try:
        result = engine.cluster(fasta_path, engine_raw_dir, processed_dir, params)
    except EngineError as exc:
        raise click.ClickException(str(exc)) from exc

    click.secho(
        f"[dnmbcluster]   clusters.parquet     "
        f"{result.n_clusters:>8} clusters / {result.n_input_sequences} sequences",
        fg="green",
    )

    # ---------- Stage 3b: membership validation + shared realignment ----------
    from .engines.realign import populate_alignment_metrics
    from .engines.validation import ValidationError, validate_clusters_table

    try:
        pre_stats = validate_clusters_table(
            result.clusters_parquet,
            check_alignment=False,
            n_input_sequences=result.n_input_sequences,
        )
    except ValidationError as exc:
        raise click.ClickException(f"post-cluster validation failed: {exc}") from exc
    click.echo(
        f"[dnmbcluster]   validation ok  ({pre_stats['n_rows']} rows / "
        f"{pre_stats['n_clusters']} clusters / {pre_stats['n_singletons']} singletons)"
    )

    if with_alignment:
        click.secho(
            "[dnmbcluster] realigning cluster members (bidirectional MMseqs2 easy-search)",
            fg="cyan",
        )
        realign_stats = populate_alignment_metrics(
            clusters_parquet=result.clusters_parquet,
            input_fasta=fasta_path,
            out_dir=raw_dir_root / "realign",
            threads=threads,
            level=level,  # type: ignore[arg-type]
            representatives_fasta=result.representatives_fasta,
        )
        click.secho(
            f"[dnmbcluster]   realign  centroids={realign_stats.n_centroids} "
            f"aligned={realign_stats.n_aligned} missing={realign_stats.n_missing}",
            fg="green",
        )
        try:
            validate_clusters_table(
                result.clusters_parquet,
                check_alignment=True,
                n_input_sequences=result.n_input_sequences,
            )
        except ValidationError as exc:
            raise click.ClickException(f"post-realign validation failed: {exc}") from exc
    else:
        click.secho(
            "[dnmbcluster]   --fast: skipping alignment enrichment; "
            "clusters.parquet alignment columns remain null",
            fg="yellow",
        )

    # ---------- Stage 4: presence/absence + pan/core + summary + Roary ----------
    from .matrix import write_presence_absence
    from .pancore import write_pan_core_curve
    from .roary_export import write_roary_csv
    from .summary import write_cluster_summary

    n_total_genomes = genome_meta.num_rows
    presence_path = processed_dir / "presence_absence.parquet"
    pan_core_path = processed_dir / "pan_core_curve.parquet"
    summary_path = processed_dir / "cluster_summary.parquet"

    pres_table = write_presence_absence(
        result.clusters_parquet, n_total_genomes, presence_path,
    )
    click.secho(
        f"[dnmbcluster]   presence_absence.parquet  {pres_table.num_rows:>6} clusters",
        fg="green",
    )

    pan_core_table = write_pan_core_curve(
        presence_path, n_total_genomes, pan_core_path, n_permutations=10, seed=0,
    )
    click.secho(
        f"[dnmbcluster]   pan_core_curve.parquet    "
        f"{pan_core_table.num_rows:>6} rows (10 permutations x {n_total_genomes} k)",
        fg="green",
    )

    summary_table = write_cluster_summary(
        result.clusters_parquet, id_map_path, presence_path, summary_path,
    )
    click.secho(
        f"[dnmbcluster]   cluster_summary.parquet   {summary_table.num_rows:>6} clusters",
        fg="green",
    )

    from .cluster_long import write_cluster_long
    cluster_long_path = processed_dir / "cluster_long.parquet"
    long_table = write_cluster_long(
        result.clusters_parquet, id_map_path, genome_meta_path,
        summary_path, cluster_long_path,
    )
    click.secho(
        f"[dnmbcluster]   cluster_long.parquet      {long_table.num_rows:>6} CDS rows",
        fg="green",
    )

    roary_csv = output / "gene_presence_absence.csv"
    roary_rtab = output / "gene_presence_absence.Rtab"
    n_cl, n_gen = write_roary_csv(
        result.clusters_parquet, id_map_path, genome_meta_path, summary_path,
        csv_out=roary_csv, rtab_out=roary_rtab,
    )
    click.secho(
        f"[dnmbcluster]   gene_presence_absence.csv  {n_cl} clusters x {n_gen} genomes",
        fg="green",
    )

    # ---------- Stage 4c: Functional categorization ----------
    from .functional import compute_functional_categories
    func_path = processed_dir / "functional_categories.parquet"
    func_table = compute_functional_categories(
        clusters_path=result.clusters_parquet,
        id_map_path=id_map_path,
        cluster_summary_path=summary_path,
        out_path=func_path,
    )
    n_classes = len(set(func_table.column("functional_class").to_pylist()))
    click.secho(
        f"[dnmbcluster]   functional_categories.parquet  "
        f"{func_table.num_rows} CDS / {n_classes} classes",
        fg="green",
    )

    # ---------- Optional: ANI + POCP genome-genome similarity ----------
    if ani:
        from .genome_similarity import run_genome_similarity
        click.secho(
            "[dnmbcluster] computing ANI (skani) + POCP matrices",
            fg="cyan",
        )
        try:
            sim = run_genome_similarity(
                clusters_path=result.clusters_parquet,
                genome_meta_path=genome_meta_path,
                raw_dir=raw_dir_root,
                processed_dir=processed_dir,
                threads=threads,
            )
            click.secho(
                f"[dnmbcluster]   ani_matrix.parquet + pocp_matrix.parquet "
                f"({sim.n_genomes}x{sim.n_genomes})",
                fg="green",
            )
        except RuntimeError as exc:
            click.secho(
                f"[dnmbcluster] ANI/POCP stage failed: {exc}", fg="red",
            )

    # ---------- Optional: Core-gene phylogenomics ----------
    if phylo:
        from .phylogenomics import run_phylogenomics
        click.secho(
            "[dnmbcluster] running core-gene phylogenomics "
            "(single-copy core → MAFFT → IQ-TREE)",
            fg="cyan",
        )
        try:
            phylo_result = run_phylogenomics(
                clusters_path=result.clusters_parquet,
                cluster_summary_path=summary_path,
                gene_table_path=gene_table_path,
                genome_meta_path=genome_meta_path,
                raw_dir=raw_dir_root,
                processed_dir=processed_dir,
                level=level,  # type: ignore[arg-type]
                threads=threads,
            )
            click.secho(
                f"[dnmbcluster]   phylo_tree.nwk  "
                f"{phylo_result.alignment_count} core clusters "
                f"/ {phylo_result.supermatrix_length} aa supermatrix "
                f"/ model {phylo_result.model}",
                fg="green",
            )
        except RuntimeError as exc:
            click.secho(
                f"[dnmbcluster] phylogenomics stage failed: {exc}", fg="red",
            )

    # BPGA-compatible single-sheet Excel (merged casting view).
    # The long-format CDS view lives in dnmb/processed/cluster_long.parquet.
    from .dnmb_excel import parse_columns_option, write_comparative_genomics_xlsx
    try:
        block_columns = parse_columns_option(columns)
    except ValueError as exc:
        raise click.ClickException(f"--columns: {exc}") from exc
    xlsx_path = output / f"comparative_genomics_{n_gen}.xlsx"
    try:
        n_cl_xlsx, n_gen_xlsx = write_comparative_genomics_xlsx(
            result.clusters_parquet, id_map_path, genome_meta_path, summary_path,
            out_path=xlsx_path,
            block_columns=block_columns,
        )
        click.secho(
            f"[dnmbcluster]   {xlsx_path.name}  "
            f"{n_cl_xlsx} clusters x {n_gen_xlsx} genomes "
            f"(merged, {2 + len(block_columns)} cols/genome)",
            fg="green",
        )
    except ImportError as exc:
        click.secho(
            f"[dnmbcluster] xlsx export skipped ({exc}). "
            f"Install xlsxwriter to enable.", fg="yellow",
        )

    # ---------- Stage 5: R visualization ----------
    from .r_bridge import run_r_plots

    click.secho(
        "[dnmbcluster] generating plots via R (DNMBcluster package)",
        fg="cyan",
    )
    try:
        plots_dir = run_r_plots(output)
        click.secho(
            f"[dnmbcluster]   plots written to {plots_dir}", fg="green",
        )
    except (RuntimeError, RBridgeError) as exc:
        click.secho(f"[dnmbcluster] R plotting failed: {exc}", fg="red")
        click.secho(
            "[dnmbcluster] pipeline continues; dataframes are still valid.",
            fg="yellow",
        )

    click.secho(
        "[dnmbcluster] pipeline complete.",
        fg="green", bold=True,
    )


@main.command("list-engines")
def list_engines() -> None:
    """List available clustering engines."""
    for engine in ENGINES:
        click.echo(engine)


@main.command("context")
@click.argument(
    "results_dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--genome", "anchor_genome",
    required=True,
    help="genome_key of the anchor CDS (as found in genome_meta.genome_key).",
)
@click.option(
    "--locus", "anchor_locus",
    required=True,
    help="locus_tag (or cds_key fallback) of the anchor CDS.",
)
@click.option(
    "--window",
    "window_bp",
    type=int, default=25000, show_default=True,
    help=(
        "Half-width of the visible window in base pairs. Default "
        "25000 → a 50 kb view centered on the anchor. Anchor is "
        "rendered at x=0 in every strain and the whole window is "
        "reflected when the anchor lies on the minus strand so "
        "every row reads left-to-right."
    ),
)
@click.option(
    "-o", "--output",
    type=click.Path(path_type=Path),
    default=None,
    help="Output PDF path. Defaults to <results_dir>/plots/context_<genome>_<locus>.pdf.",
)
def context(
    results_dir: Path,
    anchor_genome: str,
    anchor_locus: str,
    window_bp: int,
    output: Path | None,
) -> None:
    """Draw the ortholog neighborhood around a chosen (genome, locus) anchor.

    Walks the cluster_long + id_map parquet under ``results_dir`` to
    find every strain that shares the anchor's cluster, extracts the
    ±flank CDSs on the matching contig in each strain, and renders a
    stacked gggenes plot with shared-ortholog colors.
    """
    import os
    import subprocess

    from .r_bridge import find_rscript

    if output is None:
        safe = f"{anchor_genome}_{anchor_locus}".replace("/", "_")
        output = results_dir / "plots" / f"context_{safe}.pdf"
    output = output.resolve()
    results_dir = results_dir.resolve()

    rscript = find_rscript()
    r_expr = (
        "suppressPackageStartupMessages(library(DNMBcluster)); "
        f"dnmb <- DNMBcluster::load_dnmb('{results_dir}'); "
        f"DNMBcluster::context_ribbon(dnmb, anchor_genome='{anchor_genome}', "
        f"anchor_locus='{anchor_locus}', window_bp={int(window_bp)}L, "
        f"output_file='{output}')"
    )
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    try:
        subprocess.run(
            [rscript, "--vanilla", "-e", r_expr],
            check=True, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
        raise click.ClickException(
            f"context_ribbon failed:\n{stderr}"
        ) from exc
    click.secho(f"[dnmbcluster] context ribbon written to {output}", fg="green")


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
