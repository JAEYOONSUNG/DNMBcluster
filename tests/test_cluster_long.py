"""Tests for ``cluster_long.py`` — the long-format CDS view.

This view replaces the old ``long_detailed`` Excel sheet. The key
guarantees the test enforces:

1. Output schema is exactly ``CLUSTER_LONG_SCHEMA``.
2. One row per ``clusters.parquet`` row (no row drop, no duplication).
3. Alignment columns are forwarded verbatim from ``clusters.parquet``.
4. Category lookup from ``cluster_summary.parquet`` flows through.
5. genome_key resolves from genome_meta via genome_uid.
6. locus_tag falls back to cds_key when blank.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pyarrow")

import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from dnmbcluster.cluster_long import write_cluster_long  # noqa: E402
from dnmbcluster.matrix import write_presence_absence  # noqa: E402
from dnmbcluster.schemas import (  # noqa: E402
    CATEGORY_CORE,
    CATEGORY_UNIQUE,
    CLUSTER_LONG_SCHEMA,
    CLUSTERS_SCHEMA,
    GENOME_META_SCHEMA,
    ID_MAP_SCHEMA,
    validate_schema,
)
from dnmbcluster.summary import write_cluster_summary  # noqa: E402


def _write_clusters(path: Path, rows: list[dict]) -> None:
    cols = {
        "protein_uid":        [r["protein_uid"] for r in rows],
        "genome_uid":         [r["genome_uid"] for r in rows],
        "cluster_id":         [r["cluster_id"] for r in rows],
        "representative_uid": [r["representative_uid"] for r in rows],
        "is_centroid":        [r["is_centroid"] for r in rows],
        "pct_identity_fwd":   [r.get("pct_identity_fwd") for r in rows],
        "pct_identity_rev":   [r.get("pct_identity_rev") for r in rows],
        "member_coverage":    [r.get("member_coverage") for r in rows],
        "rep_coverage":       [r.get("rep_coverage") for r in rows],
        "alignment_length":   [r.get("alignment_length") for r in rows],
    }
    pq.write_table(pa.table(cols, schema=CLUSTERS_SCHEMA), path)


def _write_id_map(path: Path, rows: list[dict]) -> None:
    keys = [
        "protein_uid", "genome_uid", "gene_uid",
        "genome_key", "genome_key_source",
        "cds_key", "cds_key_source",
        "assembly_prefix", "assembly_version",
        "organism", "strain",
        "locus_tag", "protein_id", "gene",
        "product", "ec_number",
        "contig", "start", "end", "strand",
        "aa_length",
    ]
    cols = {k: [r.get(k) for r in rows] for k in keys}
    pq.write_table(pa.table(cols, schema=ID_MAP_SCHEMA), path)


def _write_genome_meta(path: Path, rows: list[dict]) -> None:
    keys = [
        "genome_uid", "genome_key", "file_path",
        "file_sha256", "file_bytes",
        "n_records", "n_cds", "n_skipped_pseudogenes",
        "organism", "strain", "definition",
        "assembly_prefix", "assembly_version",
        "total_length", "gc_percent",
    ]
    cols = {k: [r.get(k) for r in rows] for k in keys}
    pq.write_table(pa.table(cols, schema=GENOME_META_SCHEMA), path)


def test_cluster_long_end_to_end(tmp_path: Path) -> None:
    """Two genomes, two clusters: cluster 0 core (in both), cluster 1 unique."""
    clusters_path = tmp_path / "clusters.parquet"
    _write_clusters(
        clusters_path,
        [
            # cluster 0 — core
            {
                "protein_uid": 1 << 48, "genome_uid": 0, "cluster_id": 0,
                "representative_uid": 1 << 48, "is_centroid": True,
                "pct_identity_fwd": 100.0, "pct_identity_rev": 100.0,
                "member_coverage": 1.0, "rep_coverage": 1.0,
                "alignment_length": None,
            },
            {
                "protein_uid": 2 << 48, "genome_uid": 1, "cluster_id": 0,
                "representative_uid": 1 << 48, "is_centroid": False,
                "pct_identity_fwd": 92.5, "pct_identity_rev": 92.5,
                "member_coverage": 0.97, "rep_coverage": 0.95,
                "alignment_length": 340,
            },
            # cluster 1 — unique to genome 1; locus_tag is intentionally
            # blank to exercise the cds_key fallback.
            {
                "protein_uid": (2 << 48) + 1, "genome_uid": 1, "cluster_id": 1,
                "representative_uid": (2 << 48) + 1, "is_centroid": True,
                "pct_identity_fwd": 100.0, "pct_identity_rev": 100.0,
                "member_coverage": 1.0, "rep_coverage": 1.0,
                "alignment_length": None,
            },
        ],
    )

    id_map_path = tmp_path / "id_map.parquet"
    _write_id_map(
        id_map_path,
        [
            {
                "protein_uid": 1 << 48, "genome_uid": 0, "gene_uid": 0,
                "genome_key": "GCF_000000001.1", "genome_key_source": "filename_gcf",
                "cds_key": "LOC0001", "cds_key_source": "locus_tag",
                "locus_tag": "LOC0001", "gene": "dnaA",
                "product": "chromosomal replication initiator",
                "contig": "chr1", "start": 1, "end": 300, "strand": 1,
                "aa_length": 100,
            },
            {
                "protein_uid": 2 << 48, "genome_uid": 1, "gene_uid": 0,
                "genome_key": "GCF_000000002.1", "genome_key_source": "filename_gcf",
                "cds_key": "LOC1001", "cds_key_source": "locus_tag",
                "locus_tag": "LOC1001", "gene": "dnaA",
                "product": "chromosomal replication initiator",
                "contig": "chr1", "start": 1, "end": 300, "strand": 1,
                "aa_length": 100,
            },
            {
                "protein_uid": (2 << 48) + 1, "genome_uid": 1, "gene_uid": 1,
                "genome_key": "GCF_000000002.1", "genome_key_source": "filename_gcf",
                "cds_key": "cds_chr1_500_800_-", "cds_key_source": "coord",
                # locus_tag deliberately blank → must fall back to cds_key
                "locus_tag": "",
                "gene": None, "product": "hypothetical protein",
                "contig": "chr1", "start": 500, "end": 800, "strand": -1,
                "aa_length": 100,
            },
        ],
    )

    meta_path = tmp_path / "genome_meta.parquet"
    _write_genome_meta(
        meta_path,
        [
            {
                "genome_uid": 0, "genome_key": "GCF_000000001.1",
                "file_path": "/x/a.gbff",
                "file_sha256": "a" * 64, "file_bytes": 1000,
                "n_records": 1, "n_cds": 1,
            },
            {
                "genome_uid": 1, "genome_key": "GCF_000000002.1",
                "file_path": "/x/b.gbff",
                "file_sha256": "b" * 64, "file_bytes": 1000,
                "n_records": 1, "n_cds": 2,
            },
        ],
    )

    presence_path = tmp_path / "presence.parquet"
    write_presence_absence(clusters_path, n_genomes=2, out_path=presence_path)
    summary_path = tmp_path / "summary.parquet"
    write_cluster_summary(clusters_path, id_map_path, presence_path, summary_path)

    out_path = tmp_path / "cluster_long.parquet"
    table = write_cluster_long(
        clusters_path, id_map_path, meta_path, summary_path, out_path,
    )

    validate_schema(table, CLUSTER_LONG_SCHEMA, "cluster_long.parquet")
    assert table.num_rows == 3

    rows = table.to_pylist()

    # Row 1 — cluster 0 centroid in genome 0
    assert rows[0]["cluster_id"] == 0
    assert rows[0]["category"] == CATEGORY_CORE
    assert rows[0]["genome_key"] == "GCF_000000001.1"
    assert rows[0]["locus_tag"] == "LOC0001"
    assert rows[0]["gene"] == "dnaA"
    assert rows[0]["is_centroid"] is True
    assert rows[0]["pct_identity_fwd"] == pytest.approx(100.0)
    assert rows[0]["alignment_length"] is None

    # Row 2 — cluster 0 member in genome 1
    assert rows[1]["cluster_id"] == 0
    assert rows[1]["category"] == CATEGORY_CORE
    assert rows[1]["genome_key"] == "GCF_000000002.1"
    assert rows[1]["is_centroid"] is False
    assert rows[1]["pct_identity_fwd"] == pytest.approx(92.5)
    assert rows[1]["pct_identity_rev"] == pytest.approx(92.5)
    assert rows[1]["member_coverage"] == pytest.approx(0.97)
    assert rows[1]["rep_coverage"] == pytest.approx(0.95)
    assert rows[1]["alignment_length"] == 340

    # Row 3 — cluster 1 unique, locus_tag falls back to cds_key
    assert rows[2]["cluster_id"] == 1
    assert rows[2]["category"] == CATEGORY_UNIQUE
    assert rows[2]["locus_tag"] == "cds_chr1_500_800_-"
    assert rows[2]["gene"] is None
    assert rows[2]["product"] == "hypothetical protein"


def test_cluster_long_genome_uid_must_resolve(tmp_path: Path) -> None:
    """A genome_uid in clusters.parquet that's missing from genome_meta
    is a programming error — the join must raise, not silently emit blank."""
    clusters_path = tmp_path / "clusters.parquet"
    _write_clusters(
        clusters_path,
        [
            {
                "protein_uid": 1 << 48, "genome_uid": 0, "cluster_id": 0,
                "representative_uid": 1 << 48, "is_centroid": True,
                "pct_identity_fwd": 100.0, "pct_identity_rev": 100.0,
                "member_coverage": 1.0, "rep_coverage": 1.0,
            },
        ],
    )
    id_map_path = tmp_path / "id_map.parquet"
    _write_id_map(
        id_map_path,
        [
            {
                "protein_uid": 1 << 48, "genome_uid": 0, "gene_uid": 0,
                "genome_key": "G1", "genome_key_source": "filename_gcf",
                "cds_key": "L1", "cds_key_source": "locus_tag",
                "locus_tag": "L1", "gene": None, "product": None,
                "contig": "c", "start": 1, "end": 10, "strand": 1,
            },
        ],
    )
    meta_path = tmp_path / "genome_meta.parquet"
    # Empty genome_meta — genome_uid 0 is unknown here.
    _write_genome_meta(meta_path, [])

    presence_path = tmp_path / "presence.parquet"
    write_presence_absence(clusters_path, n_genomes=1, out_path=presence_path)
    summary_path = tmp_path / "summary.parquet"
    write_cluster_summary(clusters_path, id_map_path, presence_path, summary_path)

    with pytest.raises(ValueError, match="has no entry in genome_meta"):
        write_cluster_long(
            clusters_path, id_map_path, meta_path, summary_path,
            tmp_path / "cluster_long.parquet",
        )
