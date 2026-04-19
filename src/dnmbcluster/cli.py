"""DNMBcluster command-line interface."""
from __future__ import annotations

import shutil
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import click

from . import __version__


# Per-stage wall-clock timings collected during a single `run` invocation.
# A list-of-tuples (not a dict) preserves stage order for the final report.
_stage_times: list[tuple[str, float]] = []


@contextmanager
def _timed(stage: str) -> Iterator[None]:
    """Time a pipeline stage and append the result to ``_stage_times``.

    The timing record is appended whether the stage succeeds or raises,
    so partial-failure runs still show where time went before the error.
    """
    start = time.monotonic()
    try:
        yield
    finally:
        _stage_times.append((stage, time.monotonic() - start))


def _print_stage_timings(pipeline_start: float) -> None:
    """Print the accumulated per-stage timings table.

    Called both on normal completion and from a finally-block on early
    exit (exception, --parse-only, KeyboardInterrupt) so a failed long
    run still reveals where time was spent before the failure.
    """
    if not _stage_times:
        return
    total = time.monotonic() - pipeline_start
    click.secho("\n[dnmbcluster] stage timings:", fg="cyan", bold=True)
    for name, dur in _stage_times:
        click.echo(f"  {name:<32s} {dur:>8.1f} s")
    click.secho(
        f"  {'TOTAL (wall-clock)':<32s} {total:>8.1f} s",
        fg="cyan", bold=True,
    )

def _preflight_disk_space(output: Path, *, min_gb: float) -> None:
    """Warn loudly if ``output``'s filesystem has less than ``min_gb`` free.

    Does not raise — some environments (CI, short runs) legitimately
    have tight but sufficient headroom. The pipeline will still fail
    mid-way if the floor is actually breached; the warning just tells
    the user where the error is going to come from.
    """
    try:
        usage = shutil.disk_usage(output)
    except OSError as exc:
        click.secho(
            f"[dnmbcluster] could not check disk space for {output}: {exc}",
            fg="yellow",
        )
        return
    free_gb = usage.free / (1024 ** 3)
    if free_gb < min_gb:
        click.secho(
            f"[dnmbcluster] WARNING: only {free_gb:.1f} GB free on "
            f"{output}'s filesystem (recommended: >= {min_gb:.0f} GB). "
            f"Clustering and realign stages may fail mid-run.",
            fg="yellow", bold=True,
        )


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
        "Run core-gene phylogenomics stage (single-copy core → FAMSA/"
        "MUSCLE5/MAFFT → ClipKIT/trimAl → IQ-TREE 3 with ModelFinder "
        "+ UFBoot → ggtree visualization). Off by default; adds 3–8 "
        "minutes on a 10-genome run. Requires at least one aligner and "
        "one iqtree binary on PATH; best tool wins automatically."
    ),
)
@click.option(
    "--generax/--no-generax",
    "generax",
    default=True,
    help=(
        "When --orthofinder-like is set, run GeneRax UndatedDTL if the "
        "'generax' binary is on PATH and replace heuristic DTL with "
        "rigorous ML reconciliation. Default: on (auto-skipped if "
        "binary absent)."
    ),
)
@click.option(
    "--foldseek-dir",
    "foldseek_dir",
    default=None,
    type=click.Path(exists=True, file_okay=True, dir_okay=True),
    help=(
        "When --orthofinder-like is set, run Foldseek structural "
        "orthology on this directory of per-protein PDBs (or an AA "
        "FASTA if --foldseek-mode=prosT5) and annotate HOG TSVs with "
        "structural support scores. Requires the 'foldseek' binary."
    ),
)
@click.option(
    "--foldseek-mode",
    "foldseek_mode",
    type=click.Choice(["pdb", "prosT5"]),
    default="pdb",
    help="Foldseek input mode: 'pdb' (default) or 'prosT5' (AA-only).",
)
@click.option(
    "--foldseek-prosT5",
    "foldseek_prosT5",
    default=None,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path to ProstT5 weights; required for --foldseek-mode=prosT5.",
)
@click.option(
    "--codeml/--no-codeml",
    "codeml",
    default=False,
    help=(
        "Run PAML codeml M0 dN/dS on single-copy-core OGs after the "
        "OrthoFinder-parity stage. Requires a codon-aware nucleotide "
        "alignment per OG at orthogroups/OG_*/aln.nuc.phy and the "
        "'codeml' binary on PATH."
    ),
)
@click.option(
    "--orthofinder-like/--no-orthofinder-like",
    "orthofinder_like",
    default=False,
    help=(
        "Run the SOTA pure-R OrthoFinder-parity stage: graph-refined "
        "OGs (Louvain or MCL), ML gene trees with bootstrap, auto "
        "partitioned-ML supermatrix, STRIDE/MAD rooting, HOGs, DTL "
        "reconciliation, and species-tree event overlay. Off by "
        "default; adds a few minutes on a 10-genome run. Writes "
        "under <output>/orthofinder_like/."
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
@click.option(
    "--annotate/--no-annotate",
    default=False,
    help=(
        "Run fast EggNOG annotation (DIAMOND blastp --fast against "
        "cached eggnog_proteins.dmnd + SQLite COG/KEGG lookup). "
        "~10x faster than eggnog-mapper. Requires the EggNOG DB "
        "at ~/.dnmb-cache/db_modules/eggnog/data/."
    ),
)
@click.option(
    "--hyphy/--no-hyphy",
    "hyphy",
    default=False,
    help=(
        "When --orthofinder-like is set, run HyPhy FEL/BUSTED/aBSREL "
        "per single-copy OG for publication-grade selection analysis "
        "(2023+ microbial-evolution standard). Requires the 'hyphy' "
        "binary on PATH and codon-aware nucleotide CDS (auto-extracted "
        "from GenBank into dnmb/inputs/cds.fna)."
    ),
)
@click.option(
    "--max",
    "max_preset",
    is_flag=True,
    default=False,
    help=(
        "ONE-FLAG MAXIMUM-INFO PRESET. Flips on every optional stage: "
        "--phylo --ani --annotate --orthofinder-like --codeml --hyphy. "
        "Individual --no-<stage> flags still override. Intended as the "
        "default invocation for reviewer-ready comparative-genomics runs."
    ),
)
@click.option(
    "--resume/--no-resume",
    default=False,
    help=(
        "Skip stages whose primary output already exists on disk. "
        "Useful for re-running with different --phylo / --ani / "
        "--columns settings without repeating the expensive "
        "clustering + realignment stages."
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
    orthofinder_like: bool,
    generax: bool,
    foldseek_dir: str | None,
    foldseek_mode: str,
    foldseek_prosT5: str | None,
    codeml: bool,
    ani: bool,
    annotate: bool,
    hyphy: bool,
    max_preset: bool,
    resume: bool,
) -> None:
    """Run the pipeline on a folder of GenBank files."""
    # --max preset: flip every optional stage ON unless the user
    # explicitly passed --no-<stage>. We detect "explicitly set" via
    # Click's parameter-source API so `--max --no-phylo` still skips
    # phylo — the only sane override semantics.
    if max_preset:
        ctx = click.get_current_context()
        def _was_explicit(name: str) -> bool:
            src = ctx.get_parameter_source(name)
            return src == click.core.ParameterSource.COMMANDLINE
        if not _was_explicit("phylo"):
            phylo = True
        if not _was_explicit("ani"):
            ani = True
        if not _was_explicit("annotate"):
            annotate = True
        if not _was_explicit("orthofinder_like"):
            orthofinder_like = True
        if not _was_explicit("codeml"):
            codeml = True
        if not _was_explicit("hyphy"):
            hyphy = True
        click.secho(
            "[dnmbcluster] --max preset: "
            f"phylo={phylo} ani={ani} annotate={annotate} "
            f"orthofinder_like={orthofinder_like} codeml={codeml} hyphy={hyphy}",
            fg="cyan",
        )

    # Reset stage-timing accumulator — module-level list persists across
    # invocations within the same process (matters for tests / scripted
    # multi-run setups).
    _stage_times.clear()
    pipeline_start = time.monotonic()

    # Wrap the entire pipeline body so the timing summary prints on any
    # exit path — normal completion, --parse-only short-circuit, engine
    # failure, or Ctrl+C. Helps diagnose which stage consumed time on a
    # run that failed 15 minutes in.
    try:
        _run_body(
            input_dir=input_dir, output=output, tool=tool, level=level,
            identity=identity, coverage=coverage, threads=threads,
            max_ram=max_ram, parse_only=parse_only,
            with_alignment=with_alignment, columns=columns,
            phylo=phylo, orthofinder_like=orthofinder_like, generax=generax,
            foldseek_dir=foldseek_dir, foldseek_mode=foldseek_mode,
            foldseek_prosT5=foldseek_prosT5,
            codeml=codeml, ani=ani, annotate=annotate, hyphy=hyphy,
            resume=resume,
        )
    finally:
        _print_stage_timings(pipeline_start)


def _run_body(
    *,
    input_dir: Path, output: Path, tool: str, level: str,
    identity: float | None, coverage: float, threads: int,
    max_ram: str | None, parse_only: bool, with_alignment: bool,
    columns: str, phylo: bool, orthofinder_like: bool, generax: bool,
    foldseek_dir: str | None, foldseek_mode: str,
    foldseek_prosT5: str | None, codeml: bool, ani: bool,
    annotate: bool, hyphy: bool, resume: bool,
) -> None:
    """Actual pipeline body — see ``run()`` for the user-facing docstring."""
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

    # Pre-flight disk space check. Engine tmp dirs, realign, eggnog and
    # phylogenomics all write into ``output`` — running out mid-pipeline
    # produces confusing partial-file errors. 2 GB is a conservative
    # floor for a ~10-genome run; large datasets should have much more.
    _preflight_disk_space(output, min_gb=2.0)

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

    # --resume shortcut: if clusters.parquet already exists in
    # processed/ from a previous run, skip the expensive parse →
    # cluster → realign stages entirely and jump straight to the
    # derived-analytics stages that are fast and cheap. Useful when
    # re-running with different --phylo / --ani / --columns settings.
    clusters_parquet = processed_dir / "clusters.parquet"
    id_map_path = inputs_dir / "id_map.parquet"
    gene_table_path = inputs_dir / "gene_table.parquet"
    genome_meta_path = inputs_dir / "genome_meta.parquet"
    fasta_path = inputs_dir / ("proteins.faa" if level == "protein" else "cds.fna")

    if resume and clusters_parquet.exists() and id_map_path.exists():
        click.secho(
            "[dnmbcluster] --resume: clusters.parquet found, "
            "skipping parse → cluster → realign",
            fg="yellow",
        )
        genome_meta = pq.read_table(genome_meta_path)
        # Downstream stages reference result.clusters_parquet and
        # genome_meta.num_rows — satisfy both without re-clustering.
        from types import SimpleNamespace
        result = SimpleNamespace(clusters_parquet=clusters_parquet)
    else:
        # ---------- Stage 1: GenBank parse ----------
        click.secho(f"[dnmbcluster] parsing GenBank files from {input_dir}", fg="cyan")
        with _timed("genbank-parse"):
            parsed = parse_folder(input_dir, level=level, threads=threads)
            click.echo(f"[dnmbcluster]   {len(parsed)} genomes parsed")

            id_map, gene_table, genome_meta = build_tables(parsed, level)

            from .io_utils import atomic_write_table
            atomic_write_table(id_map, id_map_path)
            atomic_write_table(gene_table, gene_table_path)
            atomic_write_table(genome_meta, genome_meta_path)

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
        with _timed("fasta-export"):
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
        with _timed(f"clustering-{tool}"):
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
            with _timed("realign"):
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

    # All six derived tables run back-to-back in a single _timed() block
    # so a partial failure still records the time-to-failure — matches
    # every other stage's try/finally semantics.
    from .cluster_long import write_cluster_long
    from .functional import compute_functional_categories
    cluster_long_path = processed_dir / "cluster_long.parquet"
    func_path = processed_dir / "functional_categories.parquet"
    roary_csv = output / "gene_presence_absence.csv"
    roary_rtab = output / "gene_presence_absence.Rtab"
    with _timed("derived-tables"):
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

        long_table = write_cluster_long(
            result.clusters_parquet, id_map_path, genome_meta_path,
            summary_path, cluster_long_path,
        )
        click.secho(
            f"[dnmbcluster]   cluster_long.parquet      {long_table.num_rows:>6} CDS rows",
            fg="green",
        )

        n_cl, n_gen = write_roary_csv(
            result.clusters_parquet, id_map_path, genome_meta_path, summary_path,
            csv_out=roary_csv, rtab_out=roary_rtab,
        )
        click.secho(
            f"[dnmbcluster]   gene_presence_absence.csv  {n_cl} clusters x {n_gen} genomes",
            fg="green",
        )

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

    # ---------- Optional: EggNOG fast annotation ----------
    if annotate:
        from .eggnog_fast import run_eggnog_fast
        click.secho(
            "[dnmbcluster] running fast EggNOG annotation (DIAMOND + SQLite)",
            fg="cyan",
        )
        try:
            with _timed("eggnog"):
                egg_table = run_eggnog_fast(
                    input_fasta=fasta_path,
                    clusters_path=result.clusters_parquet,
                    raw_dir=raw_dir_root,
                    processed_dir=processed_dir,
                    threads=threads,
                )
            n_hits = sum(
                1 for v in egg_table.column("eggnog_hit").to_pylist() if v
            )
            click.secho(
                f"[dnmbcluster]   eggnog_annotations.parquet  "
                f"{egg_table.num_rows} CDS / {n_hits} hits",
                fg="green",
            )
        except RuntimeError as exc:
            click.secho(
                f"[dnmbcluster] EggNOG annotation failed: {exc}", fg="red",
            )

    # ---------- Optional: ANI + POCP genome-genome similarity ----------
    if ani:
        from .genome_similarity import run_genome_similarity
        click.secho(
            "[dnmbcluster] computing ANI (skani) + POCP matrices",
            fg="cyan",
        )
        try:
            with _timed("ani+pocp"):
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
    phylo_tree_path = processed_dir / "phylo_tree.nwk"
    phylo_super_path = processed_dir / "phylo_supermatrix.fa"
    if phylo and resume and phylo_tree_path.exists() and phylo_super_path.exists():
        # --resume skip: phylo is expensive (~45 min on 10-genome runs)
        # and downstream SOTA/plots consume phylo_tree.nwk which is already
        # on disk. Re-running wastes time when recovering from a post-phylo
        # crash (e.g. SOTA/plot stage failures).
        click.secho(
            f"[dnmbcluster] --resume: phylo_tree.nwk found, skipping core-gene phylogenomics",
            fg="yellow",
        )
    elif phylo:
        from .phylogenomics import (
            _find_iqtree_binary,
            _pick_aligner,
            _pick_trimmer,
            run_phylogenomics,
        )
        try:
            _aln_tool, _ = _pick_aligner()
        except RuntimeError:
            _aln_tool = "NONE"
        _trim = _pick_trimmer()
        _trim_tool = _trim[0] if _trim else "none"
        _iq = _find_iqtree_binary() or "NONE"
        click.secho(
            f"[dnmbcluster] running core-gene phylogenomics "
            f"(aligner={_aln_tool}, trim={_trim_tool}, tree={_iq})",
            fg="cyan",
        )
        try:
            with _timed("phylo (FAMSA+ClipKIT+IQ-TREE)"):
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

    # ---------- Gene gain/loss on phylo tree (if tree exists) ----------
    phylo_tree_path = processed_dir / "phylo_tree.nwk"
    if phylo_tree_path.exists():
        from .gain_loss import compute_gain_loss
        click.secho(
            "[dnmbcluster] inferring gene gain/loss (Dollo parsimony)",
            fg="cyan",
        )
        try:
            with _timed("gain-loss (Dollo)"):
                gl_table = compute_gain_loss(
                    tree_path=phylo_tree_path,
                    presence_absence_path=presence_path,
                    genome_meta_path=genome_meta_path,
                    out_path=processed_dir / "gain_loss.parquet",
                )
            total_gain = sum(gl_table.column("n_gained").to_pylist())
            total_loss = sum(gl_table.column("n_lost").to_pylist())
            click.secho(
                f"[dnmbcluster]   gain_loss.parquet  "
                f"{gl_table.num_rows} branches / "
                f"+{total_gain} gains / -{total_loss} losses",
                fg="green",
            )
        except Exception as exc:
            click.secho(
                f"[dnmbcluster] gain/loss inference failed: {exc}", fg="red",
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
        with _timed("xlsx"):
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

    # ---------- Optional: SOTA OrthoFinder-parity stage (R) ----------
    if orthofinder_like:
        from .r_bridge import run_r_sota_phylogenomics

        # HyPhy + codeml both need a codon-source nucleotide CDS FASTA.
        # The protein-level run only writes proteins.faa, so emit a
        # sidecar cds.fna keyed by the same protein_uid via a single
        # re-walk of the GenBank folder. Skipped in nucleotide mode
        # (cds.fna already exists) or if neither hyphy nor codeml is on.
        cds_source_path: Path | None = None
        if hyphy or codeml:
            if level == "nucleotide":
                cds_source_path = fasta_path
            else:
                from .cds_fasta import write_cds_nt_fasta
                cds_source_path = inputs_dir / "cds.fna"
                # mtime cache: skip regeneration if cds.fna is newer
                # than id_map.parquet (the only upstream input that
                # governs CDS key assignment). Saves ~seconds per rerun
                # on --resume flows where only the downstream flags
                # change. If the user re-parses GenBank, id_map.parquet
                # gets rewritten → its mtime jumps → cache invalidates.
                # Also sanity-check the first byte — a run killed mid-
                # write can leave a newer-mtime but malformed file.
                cache_hit = False
                if (
                    cds_source_path.exists()
                    and cds_source_path.stat().st_size > 0
                    and cds_source_path.stat().st_mtime
                    > id_map_path.stat().st_mtime
                ):
                    try:
                        with open(cds_source_path, "rb") as _fh:
                            cache_hit = _fh.read(1) == b">"
                    except OSError:
                        cache_hit = False
                if cache_hit:
                    click.secho(
                        "[dnmbcluster]   cds.fna up-to-date "
                        "(newer than id_map.parquet) — skipping regeneration",
                        fg="yellow",
                    )
                else:
                    click.secho(
                        "[dnmbcluster] emitting nucleotide CDS sidecar "
                        "(cds.fna) for codon-aware selection analysis",
                        fg="cyan",
                    )
                    try:
                        with _timed("cds-fna-sidecar"):
                            n_cds_nt, n_missing = write_cds_nt_fasta(
                                input_dir, id_map_path, cds_source_path,
                                threads=threads,
                            )
                        # Empty-output guard: id_map might be empty in
                        # pathological runs, or every CDS could be a
                        # pseudogene. In either case HyPhy/codeml would
                        # choke on a zero-byte FASTA; short-circuit instead.
                        if n_cds_nt == 0:
                            click.secho(
                                "[dnmbcluster]   cds.fna produced 0 sequences "
                                "— disabling HyPhy / codeml for this run",
                                fg="yellow",
                            )
                            cds_source_path = None
                        else:
                            n_target = n_cds_nt + n_missing
                            miss_ratio = (
                                n_missing / n_target if n_target else 0.0
                            )
                            click.secho(
                                f"[dnmbcluster]   cds.fna  {n_cds_nt} CDS "
                                f"({n_missing} missing)",
                                fg="green",
                            )
                            # High-miss warning: >5% of id_map protein_uids
                            # couldn't be recovered — indicates upstream
                            # parse drift or malformed GenBank. HyPhy/
                            # codeml will silently drop those OGs.
                            if miss_ratio > 0.05:
                                click.secho(
                                    f"[dnmbcluster]   WARNING: "
                                    f"{miss_ratio * 100:.1f}% of protein_uids "
                                    "have no NT counterpart — HyPhy / codeml "
                                    "will skip those OGs. Check id_map/GenBank "
                                    "consistency.",
                                    fg="yellow", bold=True,
                                )
                    except Exception as exc:
                        click.secho(
                            f"[dnmbcluster] cds.fna emission failed: {exc} "
                            "— HyPhy / codeml will skip.",
                            fg="yellow",
                        )
                        cds_source_path = None

        click.secho(
            "[dnmbcluster] running SOTA OrthoFinder-parity stage "
            "(graph-refined OGs, ML gene trees, bootstrap, DTL)",
            fg="cyan",
        )
        try:
            with _timed("orthofinder-like-sota"):
                of_dir = run_r_sota_phylogenomics(
                    output, threads=max(1, threads),
                    generax=generax,
                    foldseek_dir=foldseek_dir,
                    foldseek_mode=foldseek_mode,
                    foldseek_prosT5=foldseek_prosT5,
                    codeml=codeml,
                    hyphy=hyphy,
                    cds_source=cds_source_path,
                )
            # of_dir comes back absolute from the R side; output may be
            # relative (user ran `-o result`). Resolve both before computing
            # the display path so relative_to doesn't error on mismatched
            # absolute/relative Paths.
            try:
                display_of = of_dir.resolve().relative_to(output.resolve())
            except ValueError:
                display_of = of_dir
            click.secho(
                f"[dnmbcluster]   {display_of} ready "
                "(species tree + HOGs + DTL overlay)", fg="green",
            )
        except (RuntimeError, RBridgeError) as exc:
            click.secho(
                f"[dnmbcluster] SOTA phylogenomics stage failed: {exc}",
                fg="red",
            )

    # ---------- Stage 5: R visualization ----------
    from .r_bridge import run_r_plots

    click.secho(
        "[dnmbcluster] generating plots via R (DNMBcluster package)",
        fg="cyan",
    )
    try:
        with _timed("r-plots"):
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

    # ---------- Stage 6: HTML report ----------
    from .html_report import write_html_report
    try:
        with _timed("html-report"):
            report_path = write_html_report(output)
        click.secho(
            f"[dnmbcluster]   {report_path.name}",
            fg="green",
        )
    except Exception as exc:
        click.secho(f"[dnmbcluster] HTML report failed: {exc}", fg="red")

    click.secho(
        "[dnmbcluster] pipeline complete.",
        fg="green", bold=True,
    )


@main.command("list-engines")
def list_engines() -> None:
    """List available clustering engines."""
    for engine in ENGINES:
        click.echo(engine)


@main.command("info")
@click.argument(
    "results_dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
)
def info(results_dir: Path) -> None:
    """Print summary statistics from an existing results directory."""
    import pyarrow.parquet as pq

    processed = results_dir / "dnmb" / "processed"
    inputs = results_dir / "dnmb" / "inputs"

    def _read(subdir, name):
        p = subdir / name
        if not p.exists():
            return None
        return pq.read_table(p)

    meta = _read(inputs, "genome_meta.parquet")
    clusters = _read(processed, "clusters.parquet")
    summary = _read(processed, "cluster_summary.parquet")
    func = _read(processed, "functional_categories.parquet")
    ani = _read(processed, "ani_matrix.parquet")

    if meta is None or clusters is None:
        raise click.ClickException(
            f"Not a valid DNMBcluster results dir: {results_dir}"
        )

    n_genomes = meta.num_rows
    n_cds = clusters.num_rows
    summary_pyd = summary.to_pydict() if summary else {}
    cats = summary_pyd.get("category", [])
    n_core = cats.count("core")
    n_acc = cats.count("accessory")
    n_uniq = cats.count("unique")
    n_clusters = len(cats)

    click.secho(f"\n{'='*50}", bold=True)
    click.secho(f"  DNMBcluster results: {results_dir}", bold=True)
    click.secho(f"{'='*50}\n", bold=True)

    click.echo(f"  Genomes:         {n_genomes:>10,}")
    click.echo(f"  Total CDS:       {n_cds:>10,}")
    click.echo(f"  Clusters:        {n_clusters:>10,}")
    click.echo(f"    Core:          {n_core:>10,}  ({n_core/max(n_clusters,1)*100:.1f}%)")
    click.echo(f"    Accessory:     {n_acc:>10,}  ({n_acc/max(n_clusters,1)*100:.1f}%)")
    click.echo(f"    Unique:        {n_uniq:>10,}  ({n_uniq/max(n_clusters,1)*100:.1f}%)")

    if func:
        func_pyd = func.to_pydict()
        from collections import Counter
        fc = Counter(func_pyd.get("functional_class", []))
        click.echo(f"\n  Functional classes (top 5):")
        for cls, cnt in fc.most_common(5):
            click.echo(f"    {cls:<28} {cnt:>8,}  ({cnt/max(n_cds,1)*100:.1f}%)")

    if ani:
        import pyarrow.compute as pc
        ani_vals = ani.column("ani_percent").to_pylist()
        off_diag = [v for v in ani_vals if v < 99.99]
        if off_diag:
            click.echo(f"\n  ANI (off-diagonal): {min(off_diag):.1f}% – {max(off_diag):.1f}%")

    plots = list((results_dir / "plots").glob("*.pdf")) if (results_dir / "plots").exists() else []
    click.echo(f"\n  Plots:           {len(plots):>10}")

    phylo = processed / "phylo_tree.nwk"
    click.echo(f"  Phylo tree:      {'yes' if phylo.exists() else 'no':>10}")
    click.echo()


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
