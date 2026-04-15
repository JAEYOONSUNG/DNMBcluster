"""Pan/core curve computation from presence/absence bitmaps.

For a permutation ``π`` over the ``N`` input genomes and each
``k = 1..N``:

- **Pan(k)**  = number of clusters present in *at least one* of the
  first ``k`` genomes — i.e., clusters whose bitmap shares a bit with
  the "first k" bitmap.
- **Core(k)** = number of clusters present in *all* of the first ``k``
  genomes — i.e., clusters whose bitmap is a superset of the first-k
  bitmap.

The bitmap representation from ``matrix.py`` (list of uint64 words)
makes both tests efficient. M6 optimization will lift the outer loop
to NumPy uint64 array ops; this MVP implementation is pure Python for
correctness.

SPEED.md Section 5.
"""
from __future__ import annotations

import random
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .matrix import bits_for_first_k, is_subset
from .schemas import (
    PAN_CORE_CURVE_SCHEMA,
    PRESENCE_ABSENCE_SCHEMA,
    validate_schema,
)


def _cluster_has_any(cluster_bitmap: list[int], mask: list[int]) -> bool:
    """True if cluster_bitmap shares at least one bit with ``mask``."""
    for c, m in zip(cluster_bitmap, mask):
        if c & m:
            return True
    return False


def compute_pan_core_curve(
    presence_absence_path: Path,
    n_genomes: int,
    *,
    n_permutations: int = 10,
    seed: int = 0,
) -> pa.Table:
    """Compute pan/core curves across ``n_permutations`` random orders.

    Returns an Arrow table matching ``PAN_CORE_CURVE_SCHEMA`` with
    ``n_permutations * n_genomes`` rows.
    """
    table = pq.read_table(presence_absence_path)
    validate_schema(table, PRESENCE_ABSENCE_SCHEMA, "presence_absence.parquet")

    bitmaps_col = table.column("genome_bitmap").to_pylist()
    # Materialize as list of lists for fast Python iteration
    cluster_bitmaps: list[list[int]] = [list(b) for b in bitmaps_col]
    n_words = len(cluster_bitmaps[0]) if cluster_bitmaps else 1

    rng = random.Random(seed)

    permutation_col: list[int] = []
    k_col: list[int] = []
    pan_col: list[int] = []
    core_col: list[int] = []

    for perm_idx in range(n_permutations):
        order = list(range(n_genomes))
        rng.shuffle(order)

        for k in range(1, n_genomes + 1):
            mask = bits_for_first_k(order, k, n_words)
            pan_count = 0
            core_count = 0
            for cluster_bitmap in cluster_bitmaps:
                if _cluster_has_any(cluster_bitmap, mask):
                    pan_count += 1
                    if is_subset(mask, cluster_bitmap):
                        core_count += 1
            permutation_col.append(perm_idx)
            k_col.append(k)
            pan_col.append(pan_count)
            core_col.append(core_count)

    result = pa.table(
        {
            "permutation": pa.array(permutation_col, type=pa.uint16()),
            "k": pa.array(k_col, type=pa.uint16()),
            "pan": pa.array(pan_col, type=pa.uint32()),
            "core": pa.array(core_col, type=pa.uint32()),
        },
        schema=PAN_CORE_CURVE_SCHEMA,
    )
    validate_schema(result, PAN_CORE_CURVE_SCHEMA, "pan_core_curve.parquet")
    return result


def write_pan_core_curve(
    presence_absence_path: Path,
    n_genomes: int,
    out_path: Path,
    *,
    n_permutations: int = 10,
    seed: int = 0,
) -> pa.Table:
    table = compute_pan_core_curve(
        presence_absence_path,
        n_genomes,
        n_permutations=n_permutations,
        seed=seed,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, out_path, compression="zstd", compression_level=3)
    return table
