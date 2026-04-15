"""Unit tests for the MMseqs2 TSV → clusters.parquet converter.

These don't need the mmseqs2 binary — they feed a synthetic TSV to
``parse_cluster_tsv_to_parquet`` and assert the schema and semantics
of the output.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pyarrow")

import pyarrow.parquet as pq  # noqa: E402

from dnmbcluster.engines.mmseqs2 import parse_cluster_tsv_to_parquet  # noqa: E402
from dnmbcluster.ids import make_protein_uid  # noqa: E402
from dnmbcluster.schemas import CLUSTERS_SCHEMA, validate_schema  # noqa: E402


def _write_tsv(path: Path, rows: list[tuple[int, int]]) -> None:
    with open(path, "w") as fh:
        for rep, member in rows:
            fh.write(f"{rep}\t{member}\n")


def test_single_singleton_cluster(tmp_path: Path) -> None:
    tsv = tmp_path / "mmseqs_out_cluster.tsv"
    uid = make_protein_uid(0, 0)
    _write_tsv(tsv, [(uid, uid)])

    out = tmp_path / "clusters.parquet"
    n_clusters = parse_cluster_tsv_to_parquet(cluster_tsv=tsv, out_parquet=out)

    assert n_clusters == 1
    table = pq.read_table(out)
    validate_schema(table, CLUSTERS_SCHEMA, "clusters.parquet")
    assert table.num_rows == 1
    assert table.column("cluster_id").to_pylist() == [0]
    assert table.column("is_centroid").to_pylist() == [True]
    # Centroid self-row has trivially perfect identity and full coverage.
    assert table.column("pct_identity_fwd").to_pylist() == [100.0]
    assert table.column("pct_identity_rev").to_pylist() == [100.0]
    assert table.column("member_coverage").to_pylist() == [1.0]
    assert table.column("rep_coverage").to_pylist() == [1.0]


def test_cluster_with_members(tmp_path: Path) -> None:
    # Two genomes, two clusters: cluster 0 has 3 members (in genome 1),
    # cluster 1 has 2 members (one in each genome).
    rep_a = make_protein_uid(1, 0)
    members_a = [make_protein_uid(1, i) for i in (0, 1, 2)]

    rep_b = make_protein_uid(0, 5)
    members_b = [make_protein_uid(0, 5), make_protein_uid(1, 10)]

    rows: list[tuple[int, int]] = []
    for m in members_a:
        rows.append((rep_a, m))
    for m in members_b:
        rows.append((rep_b, m))

    tsv = tmp_path / "mmseqs_out_cluster.tsv"
    _write_tsv(tsv, rows)

    out = tmp_path / "clusters.parquet"
    n_clusters = parse_cluster_tsv_to_parquet(cluster_tsv=tsv, out_parquet=out)
    assert n_clusters == 2

    table = pq.read_table(out)
    validate_schema(table, CLUSTERS_SCHEMA, "clusters.parquet")
    assert table.num_rows == 5

    cluster_ids = table.column("cluster_id").to_pylist()
    # First rep seen → cluster_id 0, second rep seen → cluster_id 1.
    assert cluster_ids[:3] == [0, 0, 0]
    assert cluster_ids[3:] == [1, 1]

    centroid_flags = table.column("is_centroid").to_pylist()
    assert centroid_flags == [True, False, False, True, False]

    # genome_uid is derived from the upper 16 bits of protein_uid
    genome_uids = table.column("genome_uid").to_pylist()
    assert genome_uids == [1, 1, 1, 0, 1]


def test_bidirectional_alignment_join(tmp_path: Path) -> None:
    """Forward + reverse easy-search TSVs populate both identity columns."""
    rep = make_protein_uid(0, 0)
    member = make_protein_uid(0, 1)

    tsv = tmp_path / "mmseqs_out_cluster.tsv"
    _write_tsv(tsv, [(rep, rep), (rep, member)])

    fwd = tmp_path / "fwd.tsv"
    with open(fwd, "w") as fh:
        # query member, target rep
        fh.write(f"{member}\t{rep}\t92.5\t0.95\t0.80\t180\n")

    rev = tmp_path / "rev.tsv"
    with open(rev, "w") as fh:
        # query rep, target member
        fh.write(f"{rep}\t{member}\t91.0\t0.80\t0.95\t180\n")

    out = tmp_path / "clusters.parquet"
    parse_cluster_tsv_to_parquet(
        cluster_tsv=tsv,
        alignments_fwd=fwd,
        alignments_rev=rev,
        out_parquet=out,
    )

    table = pq.read_table(out)
    # Find the member row (not the centroid)
    member_row = None
    for i, uid in enumerate(table.column("protein_uid").to_pylist()):
        if uid == member:
            member_row = i
            break
    assert member_row is not None

    pid_fwd = table.column("pct_identity_fwd").to_pylist()[member_row]
    pid_rev = table.column("pct_identity_rev").to_pylist()[member_row]
    mcov = table.column("member_coverage").to_pylist()[member_row]
    rcov = table.column("rep_coverage").to_pylist()[member_row]
    alen = table.column("alignment_length").to_pylist()[member_row]

    assert pid_fwd is not None and abs(pid_fwd - 92.5) < 1e-3
    assert pid_rev is not None and abs(pid_rev - 91.0) < 1e-3
    assert mcov is not None and abs(mcov - 0.95) < 1e-3
    assert rcov is not None and abs(rcov - 0.80) < 1e-3
    assert alen == 180


def test_dense_cluster_id_is_contiguous(tmp_path: Path) -> None:
    """Cluster IDs are dense 0..N-1 regardless of representative_uid size."""
    big_rep = make_protein_uid(1000, 12345)
    small_rep = make_protein_uid(0, 0)
    tsv = tmp_path / "t.tsv"
    _write_tsv(
        tsv,
        [
            (big_rep, big_rep),
            (small_rep, small_rep),
            (big_rep, make_protein_uid(1000, 12346)),
        ],
    )
    out = tmp_path / "out.parquet"
    parse_cluster_tsv_to_parquet(cluster_tsv=tsv, out_parquet=out)
    table = pq.read_table(out)
    assert set(table.column("cluster_id").to_pylist()) == {0, 1}
