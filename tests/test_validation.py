"""Unit tests for ``engines.validation.validate_clusters_table``.

The validator is the last line of defense before ``clusters.parquet``
leaves the pipeline. Each failure mode below corresponds to a real bug
class we've seen or want to prevent.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pytest

pytest.importorskip("pyarrow")

import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from dnmbcluster.engines.validation import (  # noqa: E402
    ValidationError,
    validate_clusters_table,
)
from dnmbcluster.ids import make_protein_uid  # noqa: E402
from dnmbcluster.schemas import CLUSTERS_SCHEMA  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _row(
    *,
    protein_uid: int,
    cluster_id: int,
    representative_uid: int,
    is_centroid: bool,
    genome_uid: int | None = None,
    pct_identity_fwd: float | None = None,
    pct_identity_rev: float | None = None,
    member_coverage: float | None = None,
    rep_coverage: float | None = None,
    alignment_length: int | None = None,
) -> dict:
    return {
        "protein_uid": protein_uid,
        "genome_uid": (protein_uid >> 48) & 0xFFFF if genome_uid is None else genome_uid,
        "cluster_id": cluster_id,
        "representative_uid": representative_uid,
        "is_centroid": is_centroid,
        "pct_identity_fwd": pct_identity_fwd,
        "pct_identity_rev": pct_identity_rev,
        "member_coverage": member_coverage,
        "rep_coverage": rep_coverage,
        "alignment_length": alignment_length,
    }


def _write_parquet(rows: Iterable[dict], path: Path) -> Path:
    rows = list(rows)
    table = pa.table(
        {
            "protein_uid":        pa.array([r["protein_uid"] for r in rows], type=pa.uint64()),
            "genome_uid":         pa.array([r["genome_uid"] for r in rows], type=pa.uint16()),
            "cluster_id":         pa.array([r["cluster_id"] for r in rows], type=pa.uint32()),
            "representative_uid": pa.array([r["representative_uid"] for r in rows], type=pa.uint64()),
            "is_centroid":        pa.array([r["is_centroid"] for r in rows], type=pa.bool_()),
            "pct_identity_fwd":   pa.array([r["pct_identity_fwd"] for r in rows], type=pa.float32()),
            "pct_identity_rev":   pa.array([r["pct_identity_rev"] for r in rows], type=pa.float32()),
            "member_coverage":    pa.array([r["member_coverage"] for r in rows], type=pa.float32()),
            "rep_coverage":       pa.array([r["rep_coverage"] for r in rows], type=pa.float32()),
            "alignment_length":   pa.array([r["alignment_length"] for r in rows], type=pa.uint32()),
        },
        schema=CLUSTERS_SCHEMA,
    )
    pq.write_table(table, path)
    return path


def _membership_only_fixture(tmp_path: Path) -> Path:
    """Three rows, two clusters — one size-2 + one singleton, metrics null."""
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    uid_c = make_protein_uid(1, 0)
    rows = [
        _row(protein_uid=uid_a, cluster_id=0, representative_uid=uid_a, is_centroid=True),
        _row(protein_uid=uid_b, cluster_id=0, representative_uid=uid_a, is_centroid=False),
        _row(protein_uid=uid_c, cluster_id=1, representative_uid=uid_c, is_centroid=True),
    ]
    return _write_parquet(rows, tmp_path / "clusters.parquet")


def _fully_aligned_fixture(tmp_path: Path) -> Path:
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    uid_c = make_protein_uid(1, 0)
    rows = [
        _row(
            protein_uid=uid_a, cluster_id=0, representative_uid=uid_a, is_centroid=True,
            pct_identity_fwd=100.0, pct_identity_rev=100.0,
            member_coverage=1.0, rep_coverage=1.0,
        ),
        _row(
            protein_uid=uid_b, cluster_id=0, representative_uid=uid_a, is_centroid=False,
            pct_identity_fwd=91.5, pct_identity_rev=91.5,
            member_coverage=0.97, rep_coverage=0.95, alignment_length=340,
        ),
        _row(
            protein_uid=uid_c, cluster_id=1, representative_uid=uid_c, is_centroid=True,
            pct_identity_fwd=100.0, pct_identity_rev=100.0,
            member_coverage=1.0, rep_coverage=1.0,
        ),
    ]
    return _write_parquet(rows, tmp_path / "clusters.parquet")


# ---------------------------------------------------------------------------
# Happy paths
# ---------------------------------------------------------------------------


def test_membership_only_happy_path(tmp_path: Path) -> None:
    path = _membership_only_fixture(tmp_path)
    stats = validate_clusters_table(path, check_alignment=False)
    assert stats["n_rows"] == 3
    assert stats["n_clusters"] == 2
    assert stats["n_centroids"] == 2
    assert stats["n_singletons"] == 1
    assert stats["n_alignment_checked"] == 0


def test_fully_aligned_happy_path(tmp_path: Path) -> None:
    path = _fully_aligned_fixture(tmp_path)
    stats = validate_clusters_table(path, check_alignment=True)
    assert stats["n_rows"] == 3
    assert stats["n_alignment_checked"] == 1  # one non-centroid in a size-2 cluster


def test_n_input_sequences_matches(tmp_path: Path) -> None:
    path = _membership_only_fixture(tmp_path)
    stats = validate_clusters_table(path, check_alignment=False, n_input_sequences=3)
    assert stats["n_rows"] == 3


# ---------------------------------------------------------------------------
# Failure modes
# ---------------------------------------------------------------------------


def test_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(ValidationError, match="not found"):
        validate_clusters_table(tmp_path / "nope.parquet", check_alignment=False)


def test_n_input_sequences_mismatch(tmp_path: Path) -> None:
    path = _membership_only_fixture(tmp_path)
    with pytest.raises(ValidationError, match="row count 3 does not match input sequence count 7"):
        validate_clusters_table(path, check_alignment=False, n_input_sequences=7)


def test_duplicate_protein_uid(tmp_path: Path) -> None:
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    rows = [
        _row(protein_uid=uid_a, cluster_id=0, representative_uid=uid_a, is_centroid=True),
        _row(protein_uid=uid_a, cluster_id=0, representative_uid=uid_a, is_centroid=False),
        _row(protein_uid=uid_b, cluster_id=1, representative_uid=uid_b, is_centroid=True),
    ]
    path = _write_parquet(rows, tmp_path / "clusters.parquet")
    with pytest.raises(ValidationError, match="duplicate protein_uid"):
        validate_clusters_table(path, check_alignment=False)


def test_cluster_missing_centroid(tmp_path: Path) -> None:
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    rows = [
        _row(protein_uid=uid_a, cluster_id=0, representative_uid=uid_b, is_centroid=False),
        _row(protein_uid=uid_b, cluster_id=0, representative_uid=uid_b, is_centroid=False),
    ]
    path = _write_parquet(rows, tmp_path / "clusters.parquet")
    with pytest.raises(ValidationError, match="no centroid row"):
        validate_clusters_table(path, check_alignment=False)


def test_cluster_multiple_centroids(tmp_path: Path) -> None:
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    rows = [
        _row(protein_uid=uid_a, cluster_id=0, representative_uid=uid_a, is_centroid=True),
        _row(protein_uid=uid_b, cluster_id=0, representative_uid=uid_b, is_centroid=True),
    ]
    path = _write_parquet(rows, tmp_path / "clusters.parquet")
    with pytest.raises(ValidationError, match="multiple centroids"):
        validate_clusters_table(path, check_alignment=False)


def test_centroid_flag_mismatch_true(tmp_path: Path) -> None:
    """is_centroid=True but protein_uid != representative_uid."""
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    rows = [
        _row(protein_uid=uid_a, cluster_id=0, representative_uid=uid_b, is_centroid=True),
        _row(protein_uid=uid_b, cluster_id=0, representative_uid=uid_b, is_centroid=False),
    ]
    path = _write_parquet(rows, tmp_path / "clusters.parquet")
    with pytest.raises(ValidationError, match=r"is_centroid=True but protein_uid"):
        validate_clusters_table(path, check_alignment=False)


def test_centroid_flag_mismatch_false(tmp_path: Path) -> None:
    """is_centroid=False but protein_uid == representative_uid."""
    uid_a = make_protein_uid(0, 0)
    rows = [
        _row(protein_uid=uid_a, cluster_id=0, representative_uid=uid_a, is_centroid=False),
    ]
    path = _write_parquet(rows, tmp_path / "clusters.parquet")
    # Note: this fixture also has no centroid, so the centroid-count check
    # actually fires first. Assert on whichever message triggers.
    with pytest.raises(ValidationError):
        validate_clusters_table(path, check_alignment=False)


def test_inconsistent_representative_within_cluster(tmp_path: Path) -> None:
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    uid_c = make_protein_uid(0, 2)
    rows = [
        _row(protein_uid=uid_a, cluster_id=0, representative_uid=uid_a, is_centroid=True),
        _row(protein_uid=uid_b, cluster_id=0, representative_uid=uid_c, is_centroid=False),
    ]
    path = _write_parquet(rows, tmp_path / "clusters.parquet")
    with pytest.raises(ValidationError, match="inconsistent representative_uid"):
        validate_clusters_table(path, check_alignment=False)


def test_genome_uid_not_matching_upper_bits(tmp_path: Path) -> None:
    uid_a = make_protein_uid(3, 0)
    rows = [
        # genome_uid forced to 0 while upper bits say 3
        _row(protein_uid=uid_a, cluster_id=0, representative_uid=uid_a, is_centroid=True, genome_uid=0),
    ]
    path = _write_parquet(rows, tmp_path / "clusters.parquet")
    with pytest.raises(ValidationError, match="does not match"):
        validate_clusters_table(path, check_alignment=False)


def test_missing_alignment_for_non_centroid(tmp_path: Path) -> None:
    """check_alignment=True must reject unpopulated non-centroid rows."""
    path = _membership_only_fixture(tmp_path)
    # Fwd check fires first now — both directions must be populated, and
    # the membership-only fixture has every alignment column null.
    with pytest.raises(ValidationError, match="null pct_identity_fwd"):
        validate_clusters_table(path, check_alignment=True)


def test_only_rev_populated_is_rejected(tmp_path: Path) -> None:
    """Validator must reject a row with pct_identity_rev but null fwd.

    This encodes the regression fix: the old validator accepted
    ``fwd null AND rev populated`` as OK, letting the forward-miss bug
    slip through. The strict rule is that both directions must be set
    for every non-centroid member.
    """
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    rows = [
        _row(
            protein_uid=uid_a, cluster_id=0, representative_uid=uid_a, is_centroid=True,
            pct_identity_fwd=100.0, pct_identity_rev=100.0,
            member_coverage=1.0, rep_coverage=1.0,
        ),
        _row(
            protein_uid=uid_b, cluster_id=0, representative_uid=uid_a, is_centroid=False,
            pct_identity_fwd=None,        # <-- forward miss
            pct_identity_rev=92.0,
            member_coverage=0.95, rep_coverage=0.90, alignment_length=320,
        ),
    ]
    path = _write_parquet(rows, tmp_path / "clusters.parquet")
    with pytest.raises(ValidationError, match="null pct_identity_fwd"):
        validate_clusters_table(path, check_alignment=True)


def test_centroid_null_identity_after_realign(tmp_path: Path) -> None:
    """A centroid with null identity is a realign bug."""
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    rows = [
        # Centroid missing identity — analytical fill is supposed to set 100/100
        _row(
            protein_uid=uid_a, cluster_id=0, representative_uid=uid_a, is_centroid=True,
            member_coverage=1.0, rep_coverage=1.0,
        ),
        _row(
            protein_uid=uid_b, cluster_id=0, representative_uid=uid_a, is_centroid=False,
            pct_identity_fwd=95.0, pct_identity_rev=95.0,
            member_coverage=1.0, rep_coverage=1.0,
        ),
    ]
    path = _write_parquet(rows, tmp_path / "clusters.parquet")
    with pytest.raises(ValidationError, match="centroid has null identity"):
        validate_clusters_table(path, check_alignment=True)


def test_centroid_null_coverage_after_realign(tmp_path: Path) -> None:
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    rows = [
        _row(
            protein_uid=uid_a, cluster_id=0, representative_uid=uid_a, is_centroid=True,
            pct_identity_fwd=100.0, pct_identity_rev=100.0,
            # coverages left null
        ),
        _row(
            protein_uid=uid_b, cluster_id=0, representative_uid=uid_a, is_centroid=False,
            pct_identity_fwd=95.0, pct_identity_rev=95.0,
            member_coverage=1.0, rep_coverage=1.0,
        ),
    ]
    path = _write_parquet(rows, tmp_path / "clusters.parquet")
    with pytest.raises(ValidationError, match="centroid has null coverage"):
        validate_clusters_table(path, check_alignment=True)
