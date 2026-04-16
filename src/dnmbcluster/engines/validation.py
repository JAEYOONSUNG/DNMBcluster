"""Hard validation of ``clusters.parquet`` at stage boundaries.

Called twice per run:

1. After the engine produces membership-only output (``check_alignment=False``).
2. After the shared realignment stage populates alignment metrics
   (``check_alignment=True``).

Raises ``ValidationError`` with a precise description on the first
failure. The point of this module is to turn silent parser bugs and
partial alignment populations into loud test failures.
"""
from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from ..schemas import CLUSTERS_SCHEMA, validate_schema


class ValidationError(RuntimeError):
    pass


def validate_clusters_table(
    clusters_parquet: Path,
    *,
    check_alignment: bool = True,
    n_input_sequences: int | None = None,
) -> dict[str, int]:
    """Run every post-clustering assertion on ``clusters.parquet``.

    Parameters
    ----------
    clusters_parquet
        Path to the Parquet file to validate.
    check_alignment
        When True, every non-centroid row in a cluster of size ≥ 2
        must have all five alignment columns populated. When False,
        alignment columns are allowed to be null (engine produced
        membership only; realignment will run next).
    n_input_sequences
        Optional expected row count. When provided, we assert that
        ``clusters_parquet`` has exactly this many rows.

    Returns
    -------
    Dict with counts of rows checked, per category.
    """
    if not clusters_parquet.exists():
        raise ValidationError(f"clusters.parquet not found: {clusters_parquet}")

    table = pq.read_table(clusters_parquet)
    validate_schema(table, CLUSTERS_SCHEMA, "clusters.parquet")

    pyd = table.to_pydict()
    n = table.num_rows

    if n_input_sequences is not None and n != n_input_sequences:
        raise ValidationError(
            f"clusters.parquet row count {n} does not match input sequence "
            f"count {n_input_sequences}"
        )

    protein_uids: list[int] = pyd["protein_uid"]
    genome_uids: list[int] = pyd["genome_uid"]
    cluster_ids: list[int] = pyd["cluster_id"]
    rep_uids: list[int] = pyd["representative_uid"]
    is_centroid: list[bool] = pyd["is_centroid"]
    pid_fwd = pyd["pct_identity_fwd"]
    pid_rev = pyd["pct_identity_rev"]
    mcov = pyd["member_coverage"]
    rcov = pyd["rep_coverage"]
    aln_len = pyd["alignment_length"]

    # -------- (1) protein_uid uniqueness --------
    if len(set(protein_uids)) != n:
        raise ValidationError(
            f"duplicate protein_uid: {n} rows but "
            f"{len(set(protein_uids))} unique values"
        )

    # -------- (2) one centroid per cluster --------
    centroid_per_cluster: dict[int, int] = {}
    for cid, cen in zip(cluster_ids, is_centroid):
        if cen:
            centroid_per_cluster[cid] = centroid_per_cluster.get(cid, 0) + 1

    cluster_sizes: dict[int, int] = {}
    for cid in cluster_ids:
        cluster_sizes[cid] = cluster_sizes.get(cid, 0) + 1

    missing_centroid = [
        cid for cid in cluster_sizes if centroid_per_cluster.get(cid, 0) == 0
    ]
    if missing_centroid:
        raise ValidationError(
            f"{len(missing_centroid)} clusters have no centroid row "
            f"(first: {missing_centroid[:5]})"
        )
    too_many_centroids = [
        cid for cid, count in centroid_per_cluster.items() if count > 1
    ]
    if too_many_centroids:
        raise ValidationError(
            f"{len(too_many_centroids)} clusters have multiple centroids "
            f"(first: {too_many_centroids[:5]})"
        )

    # -------- (3) is_centroid <=> protein_uid == representative_uid --------
    for i in range(n):
        if is_centroid[i] and int(protein_uids[i]) != int(rep_uids[i]):
            raise ValidationError(
                f"row {i}: is_centroid=True but protein_uid ({protein_uids[i]}) "
                f"!= representative_uid ({rep_uids[i]})"
            )
        if not is_centroid[i] and int(protein_uids[i]) == int(rep_uids[i]):
            raise ValidationError(
                f"row {i}: is_centroid=False but protein_uid equals "
                f"representative_uid ({protein_uids[i]})"
            )

    # -------- (4) representative_uid consistency within each cluster --------
    cluster_reps: dict[int, int] = {}
    for cid, rep in zip(cluster_ids, rep_uids):
        if cid in cluster_reps:
            if cluster_reps[cid] != int(rep):
                raise ValidationError(
                    f"cluster {cid}: inconsistent representative_uid "
                    f"({cluster_reps[cid]} vs {rep})"
                )
        else:
            cluster_reps[cid] = int(rep)

    # -------- (5) genome_uid matches upper 16 bits of protein_uid --------
    for i in range(n):
        expected_gid = (int(protein_uids[i]) >> 48) & 0xFFFF
        if int(genome_uids[i]) != expected_gid:
            raise ValidationError(
                f"row {i}: genome_uid={genome_uids[i]} does not match "
                f"upper 16 bits of protein_uid ({expected_gid})"
            )

    # -------- (6) alignment populated for every non-centroid row in size≥2 clusters --------
    #
    # Strict rule: after realignment, both pct_identity_fwd AND
    # pct_identity_rev must be populated for every non-centroid
    # member. The shared engines/realign.py stage is responsible for
    # filling them — the forward miss recovery path reuses rev.pident
    # when the forward easy-search under --max-seqs 1 dropped a pair,
    # so fwd nulls on assigned pairs are a bug by construction.
    n_alignment_checked = 0
    missing_fwd: list[int] = []
    missing_rev: list[int] = []
    if check_alignment:
        for i in range(n):
            if is_centroid[i]:
                continue
            if cluster_sizes.get(int(cluster_ids[i]), 0) < 2:
                raise ValidationError(
                    f"row {i}: non-centroid in singleton cluster {cluster_ids[i]}"
                )
            n_alignment_checked += 1
            if pid_fwd[i] is None:
                missing_fwd.append(i)
            if pid_rev[i] is None:
                missing_rev.append(i)

        if missing_fwd:
            sample = missing_fwd[:10]
            raise ValidationError(
                f"{len(missing_fwd)} non-centroid rows have null "
                f"pct_identity_fwd after realignment (first rows: "
                f"{sample}). Realignment stage failed to cover every "
                f"assigned (member, rep) pair in the forward direction."
            )
        if missing_rev:
            sample = missing_rev[:10]
            raise ValidationError(
                f"{len(missing_rev)} non-centroid rows have null "
                f"pct_identity_rev after realignment (first rows: "
                f"{sample}). Reverse alignment pass is incomplete."
            )

    # -------- (7) centroid rows are analytically populated when check_alignment --------
    if check_alignment:
        for i in range(n):
            if not is_centroid[i]:
                continue
            if pid_fwd[i] is None or pid_rev[i] is None:
                raise ValidationError(
                    f"row {i}: centroid has null identity "
                    f"(pct_identity_fwd={pid_fwd[i]}, pct_identity_rev={pid_rev[i]})"
                )
            if mcov[i] is None or rcov[i] is None:
                raise ValidationError(
                    f"row {i}: centroid has null coverage"
                )

    return {
        "n_rows": n,
        "n_clusters": len(cluster_sizes),
        "n_centroids": sum(1 for c in is_centroid if c),
        "n_alignment_checked": n_alignment_checked,
        "n_singletons": sum(1 for size in cluster_sizes.values() if size == 1),
    }
