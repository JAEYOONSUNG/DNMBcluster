"""Unit tests for presence/absence bitmap + pan/core curve."""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pyarrow")

import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from dnmbcluster.matrix import (  # noqa: E402
    bitmap_n_words,
    build_presence_absence,
    genomes_to_bitmap,
    is_subset,
    popcount,
    union,
    write_presence_absence,
)
from dnmbcluster.pancore import compute_pan_core_curve  # noqa: E402
from dnmbcluster.schemas import (  # noqa: E402
    CATEGORY_ACCESSORY,
    CATEGORY_CORE,
    CATEGORY_UNIQUE,
    CLUSTERS_SCHEMA,
    categorize,
)


# ---------------------------------------------------------------------------
# Bitmap primitives
# ---------------------------------------------------------------------------


def test_bitmap_n_words() -> None:
    assert bitmap_n_words(1) == 1
    assert bitmap_n_words(64) == 1
    assert bitmap_n_words(65) == 2
    assert bitmap_n_words(128) == 2
    assert bitmap_n_words(129) == 3


def test_genomes_to_bitmap_single_word() -> None:
    bm = genomes_to_bitmap([0, 3, 5], n_words=1)
    assert bm == [0b101001]
    assert popcount(bm) == 3


def test_genomes_to_bitmap_multi_word() -> None:
    bm = genomes_to_bitmap([0, 64, 65], n_words=2)
    assert bm[0] == 1
    assert bm[1] == (1 | 2)


def test_union_and_subset() -> None:
    a = [0b0011, 0]
    b = [0b0101, 0]
    assert union(a, b) == [0b0111, 0]
    assert is_subset(a, [0b0111, 0])
    assert not is_subset(a, [0b0100, 0])


def test_categorize() -> None:
    # All 10 → core
    assert categorize(10, 10) == CATEGORY_CORE
    # Strict majority short of full → accessory (no soft_core tier any more)
    assert categorize(10, 11) == CATEGORY_ACCESSORY
    assert categorize(19, 20) == CATEGORY_ACCESSORY
    assert categorize(2, 20) == CATEGORY_ACCESSORY
    # Exactly one genome → unique
    assert categorize(1, 20) == CATEGORY_UNIQUE
    # Single-genome corpus: trivially core
    assert categorize(1, 1) == CATEGORY_CORE


# ---------------------------------------------------------------------------
# build_presence_absence
# ---------------------------------------------------------------------------


def _write_clusters_parquet(path: Path, rows: list[dict]) -> None:
    cols = {
        "protein_uid": [r["protein_uid"] for r in rows],
        "genome_uid": [r["genome_uid"] for r in rows],
        "cluster_id": [r["cluster_id"] for r in rows],
        "representative_uid": [r["representative_uid"] for r in rows],
        "is_centroid": [r["is_centroid"] for r in rows],
        "pct_identity_fwd": [r.get("pct_identity_fwd") for r in rows],
        "pct_identity_rev": [r.get("pct_identity_rev") for r in rows],
        "member_coverage": [r.get("member_coverage") for r in rows],
        "rep_coverage": [r.get("rep_coverage") for r in rows],
        "alignment_length": [r.get("alignment_length") for r in rows],
    }
    table = pa.table(cols, schema=CLUSTERS_SCHEMA)
    pq.write_table(table, path)


def test_build_presence_absence_small(tmp_path: Path) -> None:
    # 3 genomes, 2 clusters: cluster 0 is core (in 0,1,2), cluster 1 only in 2.
    cluster_rows = [
        # cluster 0
        {"protein_uid": 0x0000_0000_0000_0000, "genome_uid": 0, "cluster_id": 0, "representative_uid": 0, "is_centroid": True},
        {"protein_uid": 0x0001_0000_0000_0000, "genome_uid": 1, "cluster_id": 0, "representative_uid": 0, "is_centroid": False},
        {"protein_uid": 0x0002_0000_0000_0000, "genome_uid": 2, "cluster_id": 0, "representative_uid": 0, "is_centroid": False},
        # cluster 1 (unique — only genome 2)
        {"protein_uid": 0x0002_0000_0000_0001, "genome_uid": 2, "cluster_id": 1, "representative_uid": 0x0002_0000_0000_0001, "is_centroid": True},
    ]
    clusters_path = tmp_path / "clusters.parquet"
    _write_clusters_parquet(clusters_path, cluster_rows)

    table = build_presence_absence(clusters_path, n_genomes=3)
    assert table.num_rows == 2

    rows = table.to_pydict()
    assert rows["cluster_id"] == [0, 1]
    assert rows["n_genomes"] == [3, 1]
    assert rows["n_sequences"] == [3, 1]
    assert rows["category"] == [CATEGORY_CORE, CATEGORY_UNIQUE]
    assert rows["genome_bitmap"][0] == [0b111]
    assert rows["genome_bitmap"][1] == [0b100]


# ---------------------------------------------------------------------------
# Pan/core curve
# ---------------------------------------------------------------------------


def test_pan_core_curve_monotonic(tmp_path: Path) -> None:
    # Same 3-genome / 2-cluster fixture as above
    cluster_rows = [
        {"protein_uid": 0x0000_0000_0000_0000, "genome_uid": 0, "cluster_id": 0, "representative_uid": 0, "is_centroid": True},
        {"protein_uid": 0x0001_0000_0000_0000, "genome_uid": 1, "cluster_id": 0, "representative_uid": 0, "is_centroid": False},
        {"protein_uid": 0x0002_0000_0000_0000, "genome_uid": 2, "cluster_id": 0, "representative_uid": 0, "is_centroid": False},
        {"protein_uid": 0x0002_0000_0000_0001, "genome_uid": 2, "cluster_id": 1, "representative_uid": 0x0002_0000_0000_0001, "is_centroid": True},
    ]
    clusters_path = tmp_path / "clusters.parquet"
    _write_clusters_parquet(clusters_path, cluster_rows)

    presence_path = tmp_path / "pa.parquet"
    write_presence_absence(clusters_path, 3, presence_path)

    curve = compute_pan_core_curve(presence_path, n_genomes=3, n_permutations=5, seed=42)
    rows = curve.to_pydict()
    # Invariants that must always hold:
    # - pan(k) is non-decreasing with k
    # - core(k) is non-increasing with k
    # - pan(k) >= core(k)
    # - pan(N) == total cluster count (2 in this fixture)
    by_perm: dict[int, list[tuple[int, int, int]]] = {}
    for p, k, pan, core in zip(
        rows["permutation"], rows["k"], rows["pan"], rows["core"],
    ):
        by_perm.setdefault(p, []).append((k, pan, core))

    for perm_rows in by_perm.values():
        perm_rows.sort()
        prev_pan = 0
        prev_core = 10**9
        for k, pan, core in perm_rows:
            assert pan >= prev_pan, f"pan regressed at k={k}"
            assert core <= prev_core, f"core increased at k={k}"
            assert pan >= core
            prev_pan, prev_core = pan, core
        # Final k=N
        last_k, last_pan, last_core = perm_rows[-1]
        assert last_k == 3
        assert last_pan == 2  # total clusters
        # core at k=N is 1 (only cluster 0 is present in all 3)
        assert last_core == 1
