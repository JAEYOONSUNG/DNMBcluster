"""Core-gene phylogenomics: concat → MAFFT → IQ-TREE best tree.

Speed-tuned for mid-size bacterial pan-genomes. Key decisions:

- **Only single-copy core clusters** are used as markers. Multi-copy
  clusters would require paralog resolution; we punt on that for v1.
- **MAFFT ``--auto``** picks the right algorithm per cluster size;
  runs are parallelized across clusters via a ProcessPoolExecutor.
- **IQ-TREE ``--fast`` mode** with a fixed substitution model (LG+G4
  for protein, GTR+G4 for nucleotide) — skips ModelFinder entirely,
  which is the usual wall-time killer on bacterial core sets.
- **SH-aLRT + UFBoot** (1000 replicates each) for quick branch support
  without dragging the runtime to traditional ML bootstrap levels.

Optional ``trimal -automated1`` alignment trimming runs after MAFFT
when the binary is available; skipped silently otherwise so the stage
degrades gracefully on minimal installs.
"""
from __future__ import annotations

import concurrent.futures
import logging
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pyarrow.parquet as pq

log = logging.getLogger(__name__)


@dataclass
class PhyloResult:
    """Summary of a phylogenomics run."""

    core_cluster_count: int
    alignment_count: int        # clusters that survived MAFFT + trimming
    supermatrix_path: Path
    supermatrix_length: int     # total aligned columns after concat
    treefile: Path
    iqtree_log: Path
    model: str


# ---------------------------------------------------------------------------
# Core gene selection
# ---------------------------------------------------------------------------


def select_single_copy_core(
    clusters_path: Path,
    cluster_summary_path: Path,
    n_total_genomes: int,
) -> list[int]:
    """Return cluster_ids that are single-copy core across every genome.

    Definition: category == "core" AND every genome contributes
    exactly one member (no in-genome paralogs). Multi-copy core
    clusters are excluded because alignment + concatenation across
    paralogs is not well-defined for phylogenomics.
    """
    summary = pq.read_table(cluster_summary_path).to_pydict()
    core_ids = {
        cid for cid, cat in zip(summary["cluster_id"], summary["category"])
        if cat == "core"
    }
    if not core_ids:
        return []

    clusters = pq.read_table(clusters_path).to_pydict()
    cluster_ids = clusters["cluster_id"]
    genome_uids = clusters["genome_uid"]

    per_cluster_counts: dict[int, dict[int, int]] = {}
    for cid, gid in zip(cluster_ids, genome_uids):
        if cid not in core_ids:
            continue
        per_cluster_counts.setdefault(cid, {})
        per_cluster_counts[cid][gid] = per_cluster_counts[cid].get(gid, 0) + 1

    single_copy = [
        cid for cid, genome_map in per_cluster_counts.items()
        if len(genome_map) == n_total_genomes
        and all(cnt == 1 for cnt in genome_map.values())
    ]
    return sorted(single_copy)


# ---------------------------------------------------------------------------
# Per-cluster FASTA extraction
# ---------------------------------------------------------------------------


def _load_lookups(
    clusters_path: Path,
    gene_table_path: Path,
    genome_meta_path: Path,
    level: str,
) -> tuple[dict[int, str], dict[int, str], dict[int, list[int]]]:
    """Return ``(protein_uid_to_seq, genome_uid_to_key, cluster_to_uids)``
    for fast cluster-by-cluster sequence extraction.
    """
    seq_col = "translation" if level == "protein" else "nt_sequence"

    gene = pq.read_table(gene_table_path, columns=["protein_uid", seq_col]).to_pydict()
    pid_to_seq: dict[int, str] = {
        uid: seq for uid, seq in zip(gene["protein_uid"], gene[seq_col])
    }

    meta = pq.read_table(genome_meta_path, columns=["genome_uid", "genome_key"]).to_pydict()
    gid_to_key: dict[int, str] = {
        gid: gkey for gid, gkey in zip(meta["genome_uid"], meta["genome_key"])
    }

    clusters = pq.read_table(
        clusters_path,
        columns=["protein_uid", "genome_uid", "cluster_id"],
    ).to_pydict()
    cluster_to_members: dict[int, list[tuple[int, int]]] = {}
    for cid, pid, gid in zip(
        clusters["cluster_id"], clusters["protein_uid"], clusters["genome_uid"],
    ):
        cluster_to_members.setdefault(cid, []).append((pid, gid))

    return pid_to_seq, gid_to_key, cluster_to_members


def _write_cluster_fasta(
    cluster_id: int,
    members: list[tuple[int, int]],
    pid_to_seq: dict[int, str],
    gid_to_key: dict[int, str],
    out_path: Path,
) -> int:
    """Write ``>{genome_key}\\n{sequence}`` for every member of one
    single-copy core cluster. Returns the number of records written.
    """
    with open(out_path, "w") as fh:
        for pid, gid in members:
            seq = pid_to_seq.get(pid)
            gkey = gid_to_key.get(gid)
            if seq is None or gkey is None:
                continue
            fh.write(f">{gkey}\n{seq}\n")
    return len(members)


# ---------------------------------------------------------------------------
# MAFFT + optional trimAl
# ---------------------------------------------------------------------------


def _run_mafft(input_fasta: Path, output_fasta: Path, binary: str = "mafft") -> None:
    cmd = [binary, "--auto", "--thread", "1", str(input_fasta)]
    try:
        res = subprocess.run(
            cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
        raise RuntimeError(f"mafft failed on {input_fasta.name}: {stderr}") from exc
    output_fasta.write_bytes(res.stdout)


def _run_trimal(
    input_fasta: Path, output_fasta: Path, binary: str = "trimal",
) -> None:
    cmd = [binary, "-in", str(input_fasta), "-out", str(output_fasta), "-automated1"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def _align_one_cluster(args: tuple[Path, Path, Path, str, bool, str, str]) -> Path | None:
    """Worker for the ProcessPoolExecutor. Takes a tuple so it's
    pickleable. Returns the path of the final (trimmed) alignment or
    None on failure.
    """
    raw_fasta, aln_fasta, trimmed_fasta, mafft_bin, use_trimal, trimal_bin, _level = args
    try:
        _run_mafft(raw_fasta, aln_fasta, binary=mafft_bin)
        final = aln_fasta
        if use_trimal:
            _run_trimal(aln_fasta, trimmed_fasta, binary=trimal_bin)
            final = trimmed_fasta
        return final
    except Exception as exc:  # pragma: no cover - defensive
        log.warning("align cluster failed (%s): %s", raw_fasta.name, exc)
        return None


# ---------------------------------------------------------------------------
# Supermatrix concat
# ---------------------------------------------------------------------------


def _read_aligned_fasta(path: Path) -> dict[str, str]:
    """Parse a single aligned multi-FASTA into ``{header: sequence}``."""
    seqs: dict[str, list[str]] = {}
    current: str | None = None
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current is not None:
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def _concat_alignments(
    alignment_paths: list[Path],
    genome_keys: list[str],
    out_path: Path,
) -> int:
    """Concatenate aligned per-cluster FASTAs into a single supermatrix.

    Every genome key contributes exactly one supermatrix row; gaps
    fill in for any cluster that somehow ended up with fewer than
    n_genomes rows (should be impossible for single-copy core — we
    enforce it via ``select_single_copy_core`` — but we guard anyway).

    Returns the total alignment length.
    """
    per_genome: dict[str, list[str]] = {g: [] for g in genome_keys}
    total_len = 0

    for ap in alignment_paths:
        if ap is None:
            continue
        aligned = _read_aligned_fasta(ap)
        # Determine this cluster's alignment length from any record.
        try:
            this_len = len(next(iter(aligned.values())))
        except StopIteration:
            continue
        total_len += this_len
        gap_fill = "-" * this_len
        for g in genome_keys:
            per_genome[g].append(aligned.get(g, gap_fill))

    with open(out_path, "w") as out:
        for g in genome_keys:
            out.write(f">{g}\n{''.join(per_genome[g])}\n")
    return total_len


# ---------------------------------------------------------------------------
# IQ-TREE
# ---------------------------------------------------------------------------


def _find_iqtree_binary() -> str | None:
    for name in ("iqtree2", "iqtree"):
        if shutil.which(name):
            return name
    return None


def _run_iqtree(
    supermatrix: Path,
    out_dir: Path,
    level: str,
    threads: int,
) -> tuple[Path, Path, str]:
    binary = _find_iqtree_binary()
    if binary is None:
        raise RuntimeError(
            "IQ-TREE binary not found on PATH (tried iqtree2, iqtree)."
        )

    model = "LG+G4" if level == "protein" else "GTR+G4"
    prefix = out_dir / "supermatrix"
    threads_arg = str(threads) if threads > 0 else "AUTO"

    # IQ-TREE's ``--fast`` is incompatible with ultrafast bootstrap
    # (``-bb``). Since we want support values on every node for the
    # ggtree visualization, drop ``--fast`` and rely on a single-model
    # (no ModelFinder) + UFBoot + SH-aLRT recipe. On 10 bacterial
    # genomes with ~1k single-copy core clusters this is ~1-3 min.
    cmd = [
        binary,
        "-s", str(supermatrix),
        "-m", model,
        "-bb", "1000",
        "-alrt", "1000",
        "-T", threads_arg,
        "--prefix", str(prefix),
        "-redo",
    ]
    log.info("running iqtree: %s", " ".join(cmd))
    try:
        subprocess.run(
            cmd, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            cwd=str(out_dir),
        )
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
        raise RuntimeError(f"iqtree failed: {stderr}") from exc

    treefile = prefix.with_suffix(".treefile")
    logfile = prefix.with_suffix(".log")
    if not treefile.exists():
        raise RuntimeError(f"iqtree did not produce {treefile}")
    return treefile, logfile, model


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def run_phylogenomics(
    *,
    clusters_path: Path,
    cluster_summary_path: Path,
    gene_table_path: Path,
    genome_meta_path: Path,
    raw_dir: Path,
    processed_dir: Path,
    level: str = "protein",
    threads: int = 0,
    mafft_binary: str = "mafft",
    trimal_binary: str = "trimal",
    max_clusters: int | None = None,
) -> PhyloResult:
    """Run the full core-gene phylogenomics pipeline.

    Outputs land under ``raw_dir / "phylo"``:

    - ``aln/cluster_{cid}.faa``       — MAFFT raw alignments
    - ``trimmed/cluster_{cid}.faa``   — trimAl output (if available)
    - ``supermatrix.fa``              — concatenated alignment
    - ``supermatrix.treefile``        — IQ-TREE best tree
    - ``supermatrix.log``             — IQ-TREE log

    The treefile + supermatrix are also copied into ``processed_dir``
    as canonical artifacts for the R visualization layer to pick up.
    """
    phylo_dir = raw_dir / "phylo"
    aln_dir = phylo_dir / "aln"
    trimmed_dir = phylo_dir / "trimmed"
    seqs_dir = phylo_dir / "seqs"
    for d in (phylo_dir, aln_dir, trimmed_dir, seqs_dir):
        d.mkdir(parents=True, exist_ok=True)

    genome_meta = pq.read_table(
        genome_meta_path, columns=["genome_uid", "genome_key"]
    ).to_pydict()
    ordered_pairs = sorted(
        zip(genome_meta["genome_uid"], genome_meta["genome_key"]),
        key=lambda x: x[0],
    )
    genome_keys = [k for _, k in ordered_pairs]
    n_total = len(genome_keys)

    core_cluster_ids = select_single_copy_core(
        clusters_path, cluster_summary_path, n_total,
    )
    log.info("phylo: %d single-copy core clusters", len(core_cluster_ids))
    if not core_cluster_ids:
        raise RuntimeError(
            "No single-copy core clusters found — phylogenomics needs at "
            "least one cluster present exactly once in every input genome."
        )
    if max_clusters is not None and len(core_cluster_ids) > max_clusters:
        # Deterministic subsampling for huge runs — take the first N by
        # cluster_id so reruns on the same input reproduce the same tree.
        core_cluster_ids = core_cluster_ids[:max_clusters]
        log.info("phylo: capped at %d clusters via max_clusters", max_clusters)

    pid_to_seq, gid_to_key, cluster_to_members = _load_lookups(
        clusters_path, gene_table_path, genome_meta_path, level,
    )

    # Extract per-cluster FASTAs.
    raw_fastas: list[Path] = []
    for cid in core_cluster_ids:
        members = cluster_to_members.get(cid, [])
        out_fa = seqs_dir / f"cluster_{cid}.faa"
        n_written = _write_cluster_fasta(
            cid, members, pid_to_seq, gid_to_key, out_fa,
        )
        if n_written == n_total:
            raw_fastas.append(out_fa)

    log.info("phylo: extracted %d cluster FASTAs", len(raw_fastas))

    # Parallel MAFFT (+ trimAl if available).
    use_trimal = shutil.which(trimal_binary) is not None
    if not use_trimal:
        log.info("phylo: trimal not on PATH — skipping alignment trimming")

    n_workers = threads if threads > 0 else (os.cpu_count() or 1)
    n_workers = max(1, min(n_workers, 16))

    # Thread pool (not process pool) — the per-cluster work is a
    # subprocess call to mafft, so thread-level fanout avoids the
    # spawn-start-method import dance on macOS and costs nothing in
    # Python-side CPU. Each thread just blocks on waitpid().
    align_args = [
        (
            raw_fa,
            aln_dir / raw_fa.name,
            trimmed_dir / raw_fa.name,
            mafft_binary,
            use_trimal,
            trimal_binary,
            level,
        )
        for raw_fa in raw_fastas
    ]
    alignment_paths: list[Path] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as ex:
        for final in ex.map(_align_one_cluster, align_args):
            if final is not None:
                alignment_paths.append(final)
    log.info("phylo: %d alignments produced", len(alignment_paths))
    if not alignment_paths:
        raise RuntimeError("phylo: every MAFFT call failed")

    # Concat + IQ-TREE.
    supermatrix_path = phylo_dir / "supermatrix.fa"
    total_len = _concat_alignments(alignment_paths, genome_keys, supermatrix_path)
    log.info("phylo: supermatrix %d columns, %d taxa", total_len, n_total)

    treefile, logfile, model = _run_iqtree(
        supermatrix_path, phylo_dir, level, threads,
    )

    # Copy canonical artifacts into processed_dir so the R layer can
    # pick them up via a stable path.
    canonical_super = processed_dir / "phylo_supermatrix.fa"
    canonical_tree = processed_dir / "phylo_tree.nwk"
    processed_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(supermatrix_path, canonical_super)
    shutil.copyfile(treefile, canonical_tree)

    return PhyloResult(
        core_cluster_count=len(core_cluster_ids),
        alignment_count=len(alignment_paths),
        supermatrix_path=canonical_super,
        supermatrix_length=total_len,
        treefile=canonical_tree,
        iqtree_log=logfile,
        model=model,
    )
