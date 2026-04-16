"""Roary CSV / Rtab exporter round-trip test."""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

pytest.importorskip("pyarrow")

import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from dnmbcluster.matrix import write_presence_absence  # noqa: E402
from dnmbcluster.roary_export import write_roary_csv  # noqa: E402
from dnmbcluster.schemas import (  # noqa: E402
    CLUSTERS_SCHEMA,
    GENOME_META_SCHEMA,
    ID_MAP_SCHEMA,
)
from dnmbcluster.summary import write_cluster_summary  # noqa: E402


def _write_clusters(path: Path, rows: list[dict]) -> None:
    cols = {
        "protein_uid": [r["protein_uid"] for r in rows],
        "genome_uid": [r["genome_uid"] for r in rows],
        "cluster_id": [r["cluster_id"] for r in rows],
        "representative_uid": [r["representative_uid"] for r in rows],
        "is_centroid": [r["is_centroid"] for r in rows],
        "pct_identity_fwd": [None] * len(rows),
        "pct_identity_rev": [None] * len(rows),
        "member_coverage": [None] * len(rows),
        "rep_coverage": [None] * len(rows),
        "alignment_length": [None] * len(rows),
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


def test_roary_csv_roundtrip(tmp_path: Path) -> None:
    # Two genomes, two clusters: cluster 0 present in both (core),
    # cluster 1 only in genome 1 (unique).
    clusters_path = tmp_path / "clusters.parquet"
    _write_clusters(
        clusters_path,
        [
            {"protein_uid": 1 << 48, "genome_uid": 0, "cluster_id": 0, "representative_uid": 1 << 48, "is_centroid": True},
            {"protein_uid": 2 << 48, "genome_uid": 1, "cluster_id": 0, "representative_uid": 1 << 48, "is_centroid": False},
            {"protein_uid": (2 << 48) + 1, "genome_uid": 1, "cluster_id": 1, "representative_uid": (2 << 48) + 1, "is_centroid": True},
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
                "assembly_prefix": "GCF_000000001", "assembly_version": ".1",
                "organism": "E. foo", "strain": "A",
                "locus_tag": "LOC0001", "protein_id": "WP_0001.1",
                "gene": "dnaA", "product": "chromosomal replication initiator",
                "ec_number": None,
                "contig": "chr1", "start": 1, "end": 300, "strand": 1,
                "aa_length": 100,
            },
            {
                "protein_uid": 2 << 48, "genome_uid": 1, "gene_uid": 0,
                "genome_key": "GCF_000000002.1", "genome_key_source": "filename_gcf",
                "cds_key": "LOC1001", "cds_key_source": "locus_tag",
                "assembly_prefix": "GCF_000000002", "assembly_version": ".1",
                "organism": "E. foo", "strain": "B",
                "locus_tag": "LOC1001", "protein_id": "WP_0001.1",
                "gene": "dnaA", "product": "chromosomal replication initiator",
                "ec_number": None,
                "contig": "chr1", "start": 1, "end": 300, "strand": 1,
                "aa_length": 100,
            },
            {
                "protein_uid": (2 << 48) + 1, "genome_uid": 1, "gene_uid": 1,
                "genome_key": "GCF_000000002.1", "genome_key_source": "filename_gcf",
                "cds_key": "LOC1002", "cds_key_source": "locus_tag",
                "assembly_prefix": "GCF_000000002", "assembly_version": ".1",
                "organism": "E. foo", "strain": "B",
                "locus_tag": "LOC1002", "protein_id": "WP_0002.1",
                "gene": "xyzR", "product": "hypothetical protein",
                "ec_number": None,
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
                "n_records": 1, "n_cds": 1, "n_skipped_pseudogenes": 0,
                "organism": "E. foo", "strain": "A",
                "assembly_prefix": "GCF_000000001", "assembly_version": ".1",
                "total_length": 1000, "gc_percent": 50.0,
            },
            {
                "genome_uid": 1, "genome_key": "GCF_000000002.1",
                "file_path": "/x/b.gbff",
                "file_sha256": "b" * 64, "file_bytes": 1000,
                "n_records": 1, "n_cds": 2, "n_skipped_pseudogenes": 0,
                "organism": "E. foo", "strain": "B",
                "assembly_prefix": "GCF_000000002", "assembly_version": ".1",
                "total_length": 1000, "gc_percent": 50.0,
            },
        ],
    )

    presence_path = tmp_path / "presence.parquet"
    write_presence_absence(clusters_path, n_genomes=2, out_path=presence_path)

    summary_path = tmp_path / "summary.parquet"
    write_cluster_summary(clusters_path, id_map_path, presence_path, summary_path)

    csv_out = tmp_path / "gene_presence_absence.csv"
    rtab_out = tmp_path / "gene_presence_absence.Rtab"
    n_cl, n_gen = write_roary_csv(
        clusters_path, id_map_path, meta_path, summary_path,
        csv_out=csv_out, rtab_out=rtab_out,
    )
    assert n_cl == 2
    assert n_gen == 2

    # Verify CSV
    with open(csv_out) as fh:
        reader = list(csv.reader(fh))
    header = reader[0]
    assert header[0] == "Gene"
    assert header[-2:] == ["GCF_000000001.1", "GCF_000000002.1"]
    # Row 1 is cluster 0, row 2 is cluster 1
    row_core = reader[1]
    assert row_core[0] == "dnaA"
    assert row_core[3] == "2"  # No. isolates
    assert row_core[-2] == "LOC0001"
    assert row_core[-1] == "LOC1001"

    row_unique = reader[2]
    assert row_unique[3] == "1"
    assert row_unique[-2] == ""       # genome 0 has no member
    assert row_unique[-1] == "LOC1002"

    # Verify Rtab
    with open(rtab_out) as fh:
        rtab = fh.read().splitlines()
    assert rtab[0] == "Gene\tGCF_000000001.1\tGCF_000000002.1"
    assert rtab[1].endswith("\t1\t1")  # core cluster
    assert rtab[2].endswith("\t0\t1")  # unique cluster (only in genome 1)
