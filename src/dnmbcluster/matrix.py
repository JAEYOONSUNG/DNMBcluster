"""Presence/absence bitmap construction.

Reads ``clusters.parquet`` and ``genome_meta.parquet``, builds one
bitmap per cluster, and writes ``presence_absence.parquet``. The
bitmap representation is a list of ``uint64`` words — single-word for
small genome counts, multi-word for larger pangenomes up to the
``uint16 genome_uid`` limit (65 535 genomes).

This module is pure-Python and correctness-first (M4). M6's
optimization pass will swap the inner loop to ``numpy.bitwise_or``
vectorized over cluster rows.

SPEED.md Sections 1 and 5.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import pyarrow as pa
import pyarrow.parquet as pq  # noqa: F401 - kept for read_table users

from .io_utils import atomic_write_table
from .schemas import (
    CLUSTERS_SCHEMA,
    PRESENCE_ABSENCE_SCHEMA,
    categorize,
    validate_schema,
)


# ---------------------------------------------------------------------------
# Bitmap helpers (pure Python — swap to numpy in M6)
# ---------------------------------------------------------------------------


def bitmap_n_words(n_genomes: int) -> int:
    """Number of uint64 words needed for ``n_genomes`` bits."""
    return max(1, math.ceil(n_genomes / 64))


def genomes_to_bitmap(genome_uids: Iterable[int], n_words: int) -> list[int]:
    """Set bits at the given ``genome_uid`` positions."""
    words = [0] * n_words
    for gid in genome_uids:
        if gid < 0:
            raise ValueError(f"genome_uid must be non-negative, got {gid}")
        word_idx = gid >> 6
        bit = gid & 63
        if word_idx >= n_words:
            raise ValueError(
                f"genome_uid {gid} exceeds bitmap width {n_words * 64}"
            )
        words[word_idx] |= 1 << bit
    return words


def popcount(words: Iterable[int]) -> int:
    """Count of set bits across all words in a bitmap."""
    return sum(w.bit_count() for w in words)


def union(a: list[int], b: list[int]) -> list[int]:
    return [x | y for x, y in zip(a, b)]


def intersect(a: list[int], b: list[int]) -> list[int]:
    return [x & y for x, y in zip(a, b)]


def is_subset(subset: list[int], full: list[int]) -> bool:
    """Return True if every bit in ``subset`` is also in ``full``."""
    return all((s & ~f) == 0 for s, f in zip(subset, full))


def bits_for_first_k(permutation: list[int], k: int, n_words: int) -> list[int]:
    """Build a bitmap with the first ``k`` bits from a permutation set."""
    words = [0] * n_words
    for gid in permutation[:k]:
        words[gid >> 6] |= 1 << (gid & 63)
    return words


# ---------------------------------------------------------------------------
# Presence/absence table builder
# ---------------------------------------------------------------------------


def build_presence_absence(
    clusters_path: Path,
    n_genomes: int,
) -> pa.Table:
    """Build the ``presence_absence.parquet`` table for one clustering run.

    Returns an Arrow table matching ``PRESENCE_ABSENCE_SCHEMA``.
    """
    clusters = pq.read_table(clusters_path)
    validate_schema(clusters, CLUSTERS_SCHEMA, "clusters.parquet")

    cluster_id_col = clusters.column("cluster_id").to_pylist()
    genome_uid_col = clusters.column("genome_uid").to_pylist()

    # Group cluster_id -> set of genome_uids + total sequence count
    cluster_to_genomes: dict[int, set[int]] = {}
    cluster_to_count: dict[int, int] = {}
    for cid, gid in zip(cluster_id_col, genome_uid_col):
        cluster_to_genomes.setdefault(cid, set()).add(gid)
        cluster_to_count[cid] = cluster_to_count.get(cid, 0) + 1

    n_words = bitmap_n_words(n_genomes)

    cluster_ids: list[int] = []
    n_genomes_list: list[int] = []
    n_sequences_list: list[int] = []
    categories: list[str] = []
    bitmaps: list[list[int]] = []

    for cid in sorted(cluster_to_genomes.keys()):
        genome_set = cluster_to_genomes[cid]
        bitmap = genomes_to_bitmap(sorted(genome_set), n_words)

        cluster_ids.append(cid)
        n_genomes_list.append(len(genome_set))
        n_sequences_list.append(cluster_to_count[cid])
        categories.append(categorize(len(genome_set), n_genomes))
        bitmaps.append(bitmap)

    result = pa.table(
        {
            "cluster_id": pa.array(cluster_ids, type=pa.uint32()),
            "n_genomes": pa.array(n_genomes_list, type=pa.uint16()),
            "n_sequences": pa.array(n_sequences_list, type=pa.uint32()),
            "category": pa.array(categories, type=pa.string()),
            "genome_bitmap": pa.array(bitmaps, type=pa.list_(pa.uint64())),
        },
        schema=PRESENCE_ABSENCE_SCHEMA,
    )
    validate_schema(result, PRESENCE_ABSENCE_SCHEMA, "presence_absence.parquet")
    return result


def write_presence_absence(
    clusters_path: Path,
    n_genomes: int,
    out_path: Path,
) -> pa.Table:
    """Build and persist ``presence_absence.parquet``."""
    table = build_presence_absence(clusters_path, n_genomes)
    atomic_write_table(table, out_path)
    return table
