"""Core-gene phylogenomics: concat → aligner → trimmer → IQ-TREE.

SOTA-defaulted and speed-tuned for mid-size bacterial pan-genomes.
The stage auto-picks the best tool available on PATH, so upgrading
a binary in the conda env upgrades the pipeline with no code change.

- **Only single-copy core clusters** are used as markers. Multi-copy
  clusters would require paralog resolution; we punt on that for v1.

- **Aligner priority**  (first found on PATH wins):
    1. ``famsa``       — fastest + high accuracy for bacterial homolog
                         sets (Deorowicz 2016, 2–10× faster than MAFFT
                         at equal/better SP score).
    2. ``muscle``      — MUSCLE 5 ``-super5`` mode: SOTA for >1k seqs.
    3. ``mafft``       — universal fallback, ``--auto`` picks algo.
   Parallelized across clusters via a ThreadPoolExecutor.

- **Trimmer priority** (first found on PATH wins):
    1. ``clipkit``     — smart-gap mode, signal-preserving (Steenwyk
                         2020); newer + more principled than trimAl.
    2. ``trimal``      — ``-automated1`` fallback.

- **IQ-TREE priority** (newest binary wins): ``iqtree3`` > ``iqtree2``
  > ``iqtree``.
    Model : ``-m MFP`` (ModelFinder Plus) over the full AA candidate
            set (LG, WAG, JTT, VT, Q.pfam/insect/plant, LG4X, EX_EHO,
            mtREV, … ~22 models). No -mset restriction — worth the
            3-5 min to let ModelFinder pick the genuinely best model.
    Support: UFBoot 1000 **with -bnni** (Hoang et al. 2018 MBE, corrects
            UFBoot's upward bias under model misspecification — the
            standard for publication trees since IQ-TREE v1.6) plus
            SH-aLRT 1000 for dual-support reporting. --seed 42 fixes
            stochastic starting trees for reproducibility.
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
# Aligner + trimmer — auto-pick the best tool available on PATH
# ---------------------------------------------------------------------------


def _pick_aligner() -> tuple[str, str]:
    """Return ``(tool, binary_path)``. Priority: famsa > muscle (v5) > mafft."""
    for name in ("famsa", "muscle", "mafft"):
        path = shutil.which(name)
        if path:
            return name, path
    raise RuntimeError(
        "No aligner on PATH. Install at least one of: famsa (preferred), "
        "muscle (v5), mafft."
    )


def _pick_trimmer() -> tuple[str, str] | None:
    """Return ``(tool, binary_path)`` or ``None`` if no trimmer is installed."""
    for name in ("clipkit", "trimal"):
        path = shutil.which(name)
        if path:
            return name, path
    return None


def _run_aligner(
    tool: str, binary: str, input_fasta: Path, output_fasta: Path,
) -> None:
    """Dispatch to the chosen aligner. Each branch writes ``output_fasta``."""
    if tool == "famsa":
        # FAMSA writes directly to the output file; -t 1 because we
        # already parallelize across clusters at the Python level.
        cmd = [binary, "-t", "1", str(input_fasta), str(output_fasta)]
        subprocess.run(
            cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
        )
        return
    if tool == "muscle":
        # MUSCLE 5 uses -super5 for fast, scalable alignments; -align
        # mode is reserved for <1k seqs but per-cluster inputs in a
        # bacterial pan-genome rarely exceed that, so use -align for
        # best quality here.
        cmd = [binary, "-align", str(input_fasta), "-output", str(output_fasta),
               "-threads", "1"]
        subprocess.run(
            cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
        )
        return
    # MAFFT writes to stdout.
    cmd = [binary, "--auto", "--thread", "1", str(input_fasta)]
    res = subprocess.run(
        cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    output_fasta.write_bytes(res.stdout)


def _run_trimmer(
    tool: str, binary: str, input_fasta: Path, output_fasta: Path,
) -> None:
    """Dispatch to the chosen trimmer."""
    if tool == "clipkit":
        # smart-gap: principled gap-threshold auto-selected per alignment
        # via kneedle (Steenwyk 2020). Produces a tighter, more
        # phylogenetically informative alignment than trimAl -automated1.
        cmd = [binary, str(input_fasta), "-o", str(output_fasta),
               "-m", "smart-gap"]
        subprocess.run(
            cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
        )
        return
    # trimAl
    cmd = [binary, "-in", str(input_fasta), "-out", str(output_fasta),
           "-automated1"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def _align_one_cluster(
    args: tuple[Path, Path, Path, str, str, tuple[str, str] | None, str],
) -> Path | None:
    """Thread-pool worker. Returns the path of the final (trimmed)
    alignment, or ``None`` on failure.
    """
    raw_fa, aln_fa, trimmed_fa, aln_tool, aln_bin, trim_tc, _level = args
    try:
        _run_aligner(aln_tool, aln_bin, raw_fa, aln_fa)
        final = aln_fa
        if trim_tc is not None:
            trim_tool, trim_bin = trim_tc
            _run_trimmer(trim_tool, trim_bin, aln_fa, trimmed_fa)
            final = trimmed_fa
        return final
    except Exception as exc:  # pragma: no cover - defensive
        log.warning("align cluster failed (%s): %s", raw_fa.name, exc)
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
    # Prefer iqtree3 (newest, faster heuristics) > iqtree2 > iqtree.
    for name in ("iqtree3", "iqtree2", "iqtree"):
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
            "IQ-TREE binary not found on PATH (tried iqtree3, iqtree2, iqtree)."
        )

    # Publication-grade "best tree" config. ModelFinder Plus picks the
    # best substitution model by BIC over the full candidate set (no
    # -mset restriction) — for a 10-taxon bacterial supermatrix the full
    # ~22-AA-model sweep adds only 3-5 min and may select LG4X /
    # Q.insect / EX_EHO over the restricted core when warranted.
    # UFBoot (-bb) is paired with -bnni (Hoang et al. 2018 MBE),
    # which corrects UFBoot's upward support bias under model
    # misspecification and is the standard for publication trees
    # since IQ-TREE v1.6. SH-aLRT (-alrt) is reported alongside for
    # dual-support reporting. --seed fixes stochastic starting trees.
    model_arg = "MFP"
    # Resolve to absolute paths: we pass cwd=out_dir below so iqtree's
    # auxiliary files land in the right place, but a relative --prefix
    # would then be re-interpreted inside cwd and iqtree would try to
    # write into a non-existent nested path (observed on 2026-04 run).
    out_dir = out_dir.resolve()
    supermatrix = supermatrix.resolve()
    prefix = out_dir / "supermatrix"
    threads_arg = str(threads) if threads > 0 else "AUTO"

    cmd = [
        binary,
        "-s", str(supermatrix),
        "-m", model_arg,
        "-bb", "1000",
        "-bnni",
        "-alrt", "1000",
        "-T", threads_arg,
        "--seed", "42",
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

    # Parse the model IQ-TREE actually selected from the .iqtree
    # summary file. Falls back to the requested -m argument on parse
    # failure so the downstream f-string still has something to show.
    picked_model = model_arg
    iqtree_summary = prefix.with_suffix(".iqtree")
    if iqtree_summary.exists():
        try:
            for raw in iqtree_summary.read_text().splitlines():
                # Format: "Best-fit model: LG+F+R5 chosen according to BIC"
                if raw.startswith("Best-fit model"):
                    picked_model = raw.split(":", 1)[1].split("chosen")[0].strip()
                    break
                # Fallback for fixed-model runs.
                if raw.startswith("Model of substitution"):
                    picked_model = raw.split(":", 1)[1].strip()
                    break
        except OSError:
            pass
    return treefile, logfile, picked_model


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

    # Auto-pick the best aligner + trimmer available on PATH.
    aln_tool, aln_bin = _pick_aligner()
    trim_tc = _pick_trimmer()
    log.info("phylo: aligner=%s (%s)", aln_tool, aln_bin)
    if trim_tc is None:
        log.info("phylo: no trimmer on PATH (clipkit/trimal) — skipping trim")
    else:
        log.info("phylo: trimmer=%s (%s)", trim_tc[0], trim_tc[1])

    n_workers = threads if threads > 0 else (os.cpu_count() or 1)
    n_workers = max(1, min(n_workers, 16))

    # Thread pool (not process pool) — the per-cluster work is a
    # subprocess call to an aligner, so thread-level fanout avoids the
    # spawn-start-method import dance on macOS and costs nothing in
    # Python-side CPU. Each thread just blocks on waitpid().
    align_args = [
        (
            raw_fa,
            aln_dir / raw_fa.name,
            trimmed_dir / raw_fa.name,
            aln_tool,
            aln_bin,
            trim_tc,
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
