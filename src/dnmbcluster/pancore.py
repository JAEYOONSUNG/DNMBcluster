"""Pan/core curve computation from presence/absence bitmaps.

For a permutation ``π`` over the ``N`` input genomes and each
``k = 1..N``:

- **Pan(k)**  = number of clusters present in *at least one* of the
  first ``k`` genomes — clusters whose bitmap shares a bit with
  the "first k" bitmap.
- **Core(k)** = number of clusters present in *all* of the first ``k``
  genomes — clusters whose bitmap is a superset of the "first k" bitmap.

Implementation is vectorized via NumPy over the whole (n_clusters,
n_words) array: for each k the mask is applied once and pan/core
counts reduce to ``np.any`` / ``np.all`` over the word axis. Complexity
is ``O(n_permutations * n_genomes * n_clusters * n_words)`` wall work
but the inner 2-D ops stay in C, so even 100 000 clusters × 1000
genomes finishes in well under a second.

SPEED.md Section 5.
"""
from __future__ import annotations

import random
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from .schemas import (
    PAN_CORE_CURVE_SCHEMA,
    PRESENCE_ABSENCE_SCHEMA,
    validate_schema,
)


def _load_bitmaps(table: pa.Table) -> np.ndarray:
    """Return a ``(n_clusters, n_words)`` uint64 NumPy array.

    The ``genome_bitmap`` column is a Parquet list<uint64>. We flatten
    it to a contiguous 2-D ndarray in one pass so every downstream op
    is a pure NumPy kernel.
    """
    py_lists = table.column("genome_bitmap").to_pylist()
    n_clusters = len(py_lists)
    if n_clusters == 0:
        return np.zeros((0, 1), dtype=np.uint64)
    n_words = len(py_lists[0])
    arr = np.empty((n_clusters, n_words), dtype=np.uint64)
    for i, row in enumerate(py_lists):
        if len(row) != n_words:
            raise ValueError(
                f"bitmap row {i} has width {len(row)}, expected {n_words}"
            )
        arr[i] = row
    return arr


def _build_mask(order: np.ndarray, k: int, n_words: int) -> np.ndarray:
    """Return the first-k bitmap as a ``(n_words,)`` uint64 array."""
    mask = np.zeros(n_words, dtype=np.uint64)
    for gid in order[:k]:
        word_idx = int(gid) >> 6
        bit = int(gid) & 63
        mask[word_idx] |= np.uint64(1) << np.uint64(bit)
    return mask


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

    cluster_bitmaps = _load_bitmaps(table)  # (n_clusters, n_words)
    n_clusters, n_words = cluster_bitmaps.shape

    rng = random.Random(seed)

    # Preallocate result columns
    total_rows = n_permutations * n_genomes
    permutation_col = np.empty(total_rows, dtype=np.uint16)
    k_col = np.empty(total_rows, dtype=np.uint16)
    pan_col = np.empty(total_rows, dtype=np.uint32)
    core_col = np.empty(total_rows, dtype=np.uint32)

    row_idx = 0
    for perm_idx in range(n_permutations):
        order_list = list(range(n_genomes))
        rng.shuffle(order_list)
        order = np.asarray(order_list, dtype=np.int64)

        for k in range(1, n_genomes + 1):
            mask = _build_mask(order, k, n_words)  # (n_words,)

            if n_clusters == 0:
                pan_count = 0
                core_count = 0
            else:
                # Broadcast mask over clusters:
                #   cluster_bitmaps: (n_clusters, n_words)
                #   mask:            (n_words,)
                and_result = cluster_bitmaps & mask  # (n_clusters, n_words)

                # Pan: cluster has ANY bit in common with mask
                pan_mask = np.any(and_result != 0, axis=1)
                pan_count = int(pan_mask.sum())

                # Core: mask is a subset of cluster bitmap
                #   -> every bit of mask also set in cluster
                #   -> (cluster & mask) == mask  for every word
                core_mask = np.all(and_result == mask, axis=1)
                core_count = int(core_mask.sum())

            permutation_col[row_idx] = perm_idx
            k_col[row_idx] = k
            pan_col[row_idx] = pan_count
            core_col[row_idx] = core_count
            row_idx += 1

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
    from .io_utils import atomic_write_table
    atomic_write_table(table, out_path)
    return table
