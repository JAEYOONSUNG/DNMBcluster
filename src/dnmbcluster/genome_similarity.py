"""Genome-genome similarity: Average Nucleotide Identity (ANI) via skani
and Percentage Of Conserved Proteins (POCP) from the clustering result.

ANI uses ``skani triangle`` which is the state-of-the-art fast ANI
implementation — k-mer based, ~100x faster than pyANI for typical
bacterial collections. skani reads assembly FASTAs, not GenBank, so
this module lazily extracts a single-contig-FASTA per genome from
the original GenBank inputs referenced in ``genome_meta.file_path``.

POCP is computed directly from ``clusters.parquet`` + genome CDS
counts: for each genome pair ``(A, B)``, the conserved-protein
count for A w.r.t. B is the number of A's proteins that sit in a
cluster that also contains at least one B protein. The formula is
the standard Qin et al. 2014 ratio:

    POCP(A, B) = (C1 + C2) / (T1 + T2) * 100

where Ci is A/B's conserved count toward the other and Ti is the
total CDS count of that genome. Since our clusters already enforce
``--identity`` + ``--coverage`` thresholds at run time, membership
in a shared cluster is our homology signal — the approximation
converges to the literature definition when ``--identity`` is set
near the BLASTP POCP default of 40%.

Both outputs land as:

- ``dnmb/raw/ani/skani_triangle.tsv``  (raw skani output)
- ``dnmb/processed/ani_matrix.parquet`` (canonical long-form)
- ``dnmb/processed/pocp_matrix.parquet`` (canonical long-form)
"""
from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

log = logging.getLogger(__name__)

ANI_SPECIES_THRESHOLD: float = 95.0
POCP_GENUS_THRESHOLD: float = 60.0


@dataclass
class SimilarityResult:
    n_genomes: int
    ani_tsv: Path | None
    ani_parquet: Path
    pocp_parquet: Path


# ---------------------------------------------------------------------------
# Genome FASTA cache
# ---------------------------------------------------------------------------


def _extract_genome_fasta(
    gb_path: Path, out_path: Path, genome_key: str,
) -> None:
    """Write all contig sequences of a GenBank file to one FASTA file.

    Header = ``{genome_key}__{contig_id}`` so skani's output can be
    traced back per-contig if needed, though we only care about
    whole-genome ANI at the top level.
    """
    from Bio import SeqIO

    with open(out_path, "w") as out:
        for rec in SeqIO.parse(str(gb_path), "genbank"):
            seq_str = str(rec.seq).upper()
            if not seq_str:
                continue
            contig_id = rec.id or rec.name or "unknown"
            out.write(f">{genome_key}__{contig_id}\n")
            for i in range(0, len(seq_str), 80):
                out.write(seq_str[i : i + 80] + "\n")


def _ensure_genome_fastas(
    genome_meta_path: Path, cache_dir: Path,
) -> list[tuple[str, Path]]:
    """Materialize one FASTA per genome under ``cache_dir``.

    Returns ``[(genome_key, fasta_path), ...]`` in genome_uid order.
    Skips files that already exist (content-address via filename).
    """
    meta = pq.read_table(
        genome_meta_path, columns=["genome_uid", "genome_key", "file_path"],
    ).to_pydict()

    cache_dir.mkdir(parents=True, exist_ok=True)
    out: list[tuple[str, Path]] = []
    for uid, gkey, src in sorted(
        zip(meta["genome_uid"], meta["genome_key"], meta["file_path"])
    ):
        fasta_path = cache_dir / f"{gkey}.fna"
        if not fasta_path.exists() or fasta_path.stat().st_size == 0:
            _extract_genome_fasta(Path(src), fasta_path, gkey)
        out.append((gkey, fasta_path))
    return out


# ---------------------------------------------------------------------------
# ANI via skani
# ---------------------------------------------------------------------------


def run_skani_triangle(
    genome_fastas: list[Path], out_tsv: Path, threads: int = 0,
) -> None:
    """Run ``skani triangle --full-matrix`` on all input genomes.

    skani emits a symmetric N×N matrix in a tab-separated file —
    one header row listing genome names, then one labeled data row
    per genome. We normalize it to a long-form parquet downstream.
    """
    if shutil.which("skani") is None:
        raise RuntimeError("skani binary not found on PATH.")

    threads_arg = str(threads) if threads > 0 else "0"  # 0 = all cores
    cmd = [
        "skani", "triangle", "--full-matrix",
        "-o", str(out_tsv),
        "-t", threads_arg,
    ] + [str(p) for p in genome_fastas]

    log.info("running skani: %s", " ".join(cmd[:10]) + f" ... ({len(genome_fastas)} genomes)")
    try:
        subprocess.run(
            cmd, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
        raise RuntimeError(f"skani triangle failed: {stderr}") from exc


def _parse_skani_matrix(
    tsv_path: Path, expected_keys: list[str],
) -> dict[tuple[str, str], float]:
    """Parse a skani ``--full-matrix`` TSV into a ``(key_a, key_b)->ANI`` dict.

    skani writes the matrix with filename-stem row labels. We map
    those back to our genome_key values via the provided list.
    """
    with open(tsv_path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines:
        raise RuntimeError(f"empty skani output: {tsv_path}")

    # First line is `<N>\n` (skani triangle prints the count), then
    # the matrix. Newer skani versions may skip the count line — so
    # detect numerically.
    header_idx = 0
    try:
        int(lines[0].strip())
        header_idx = 1
    except ValueError:
        header_idx = 0

    matrix_lines = lines[header_idx:]
    # Column order == row order; skani writes each row as:
    #   genome_name  v1  v2 ... vN
    # We build a stem → genome_key index once.
    stems = [Path(ml.split("\t", 1)[0]).stem for ml in matrix_lines]
    stem_to_key = {stem: key for stem, key in zip(stems, expected_keys)}
    result: dict[tuple[str, str], float] = {}
    for i, row in enumerate(matrix_lines):
        fields = row.split("\t")
        row_stem = Path(fields[0]).stem
        row_key = stem_to_key.get(row_stem, row_stem)
        for j, val in enumerate(fields[1:]):
            col_stem = stems[j]
            col_key = stem_to_key.get(col_stem, col_stem)
            try:
                ani = float(val)
            except ValueError:
                ani = 0.0
            result[(row_key, col_key)] = ani
    return result


# ---------------------------------------------------------------------------
# POCP from clusters
# ---------------------------------------------------------------------------


def compute_pocp(
    clusters_path: Path, genome_meta_path: Path,
) -> dict[tuple[str, str], float]:
    """Compute POCP(A, B) for every ordered genome pair.

    Returns ``{(genome_key_A, genome_key_B): pocp_percent}`` for all
    pairs including diagonal (= 100.0). Diagonal is forced to 100 so
    downstream heatmap reordering by ``as.dist(100 - pocp)`` works.
    """
    meta = pq.read_table(
        genome_meta_path, columns=["genome_uid", "genome_key", "n_cds"],
    ).to_pydict()
    gid_to_key: dict[int, str] = {
        gid: gkey for gid, gkey in zip(meta["genome_uid"], meta["genome_key"])
    }
    gid_to_tcount: dict[int, int] = {
        gid: n for gid, n in zip(meta["genome_uid"], meta["n_cds"])
    }
    genome_uids = sorted(gid_to_key.keys())

    clusters = pq.read_table(
        clusters_path, columns=["cluster_id", "genome_uid"]
    ).to_pydict()
    cluster_to_genome_counts: dict[int, dict[int, int]] = {}
    for cid, gid in zip(clusters["cluster_id"], clusters["genome_uid"]):
        cluster_to_genome_counts.setdefault(cid, {})
        cluster_to_genome_counts[cid][gid] = (
            cluster_to_genome_counts[cid].get(gid, 0) + 1
        )

    # For each pair (A, B): count how many A CDSs sit in clusters
    # that also contain at least one B CDS. That's C_A_wrt_B; by
    # symmetry we also accumulate C_B_wrt_A at the same time.
    conserved: dict[tuple[int, int], int] = {}
    for cid, gmap in cluster_to_genome_counts.items():
        present = list(gmap.keys())
        for a in present:
            for b in present:
                if a == b:
                    continue
                conserved[(a, b)] = conserved.get((a, b), 0) + gmap[a]

    pocp: dict[tuple[str, str], float] = {}
    for a in genome_uids:
        key_a = gid_to_key[a]
        t_a = gid_to_tcount.get(a, 0) or 1
        for b in genome_uids:
            key_b = gid_to_key[b]
            if a == b:
                pocp[(key_a, key_b)] = 100.0
                continue
            t_b = gid_to_tcount.get(b, 0) or 1
            c_a = conserved.get((a, b), 0)
            c_b = conserved.get((b, a), 0)
            pocp[(key_a, key_b)] = (c_a + c_b) / (t_a + t_b) * 100.0
    return pocp


# ---------------------------------------------------------------------------
# Parquet writers
# ---------------------------------------------------------------------------


def _write_pairwise_parquet(
    values: dict[tuple[str, str], float],
    out_path: Path,
    value_col: str,
) -> None:
    """Persist a pairwise metric as a 3-column long-form parquet:
    ``(genome_a, genome_b, <value_col>)``.
    """
    rows = sorted(values.keys())
    genome_a = [k[0] for k in rows]
    genome_b = [k[1] for k in rows]
    vals = [values[k] for k in rows]
    table = pa.table(
        {
            "genome_a": pa.array(genome_a, type=pa.string()),
            "genome_b": pa.array(genome_b, type=pa.string()),
            value_col: pa.array(vals, type=pa.float64()),
        }
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, out_path, compression="zstd", compression_level=3)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def run_genome_similarity(
    *,
    clusters_path: Path,
    genome_meta_path: Path,
    raw_dir: Path,
    processed_dir: Path,
    threads: int = 0,
    run_ani: bool = True,
    run_pocp: bool = True,
) -> SimilarityResult:
    """Compute ANI (via skani) and POCP matrices for every genome pair.

    Paths:
      - ``raw_dir / "genomes" / {genome_key}.fna``  extracted FASTA cache
      - ``raw_dir / "ani" / skani_triangle.tsv``    raw skani matrix
      - ``processed_dir / ani_matrix.parquet``      long-form ANI
      - ``processed_dir / pocp_matrix.parquet``     long-form POCP
    """
    processed_dir.mkdir(parents=True, exist_ok=True)
    ani_parquet  = processed_dir / "ani_matrix.parquet"
    pocp_parquet = processed_dir / "pocp_matrix.parquet"
    ani_tsv: Path | None = None

    # Always read meta once for ordered genome list.
    meta = pq.read_table(
        genome_meta_path, columns=["genome_uid", "genome_key"]
    ).to_pydict()
    ordered_keys = [
        k for _, k in sorted(zip(meta["genome_uid"], meta["genome_key"]))
    ]

    if run_ani:
        cache_dir = raw_dir / "genomes"
        ani_dir = raw_dir / "ani"
        ani_dir.mkdir(parents=True, exist_ok=True)
        fastas_pairs = _ensure_genome_fastas(genome_meta_path, cache_dir)
        fasta_paths = [p for _, p in fastas_pairs]
        ani_tsv = ani_dir / "skani_triangle.tsv"
        run_skani_triangle(fasta_paths, ani_tsv, threads=threads)
        ani = _parse_skani_matrix(ani_tsv, ordered_keys)
        _write_pairwise_parquet(ani, ani_parquet, "ani_percent")
        log.info("ANI matrix written to %s", ani_parquet)

    if run_pocp:
        pocp = compute_pocp(clusters_path, genome_meta_path)
        _write_pairwise_parquet(pocp, pocp_parquet, "pocp_percent")
        log.info("POCP matrix written to %s", pocp_parquet)

    return SimilarityResult(
        n_genomes=len(ordered_keys),
        ani_tsv=ani_tsv,
        ani_parquet=ani_parquet,
        pocp_parquet=pocp_parquet,
    )
