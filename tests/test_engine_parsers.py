"""Unit tests for the DIAMOND / CD-HIT / usearch12 output parsers.

All synthetic — none of these need the actual binaries.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pyarrow")

import pyarrow.parquet as pq  # noqa: E402

from dnmbcluster.engines.cdhit import parse_cdhit_clstr  # noqa: E402
from dnmbcluster.engines.diamond import parse_diamond_cluster_tsv  # noqa: E402
from dnmbcluster.engines.usearch12 import parse_uc  # noqa: E402
from dnmbcluster.ids import make_protein_uid  # noqa: E402
from dnmbcluster.schemas import CLUSTERS_SCHEMA, validate_schema  # noqa: E402


# ---------------------------------------------------------------------------
# DIAMOND deepclust parser
# ---------------------------------------------------------------------------


def test_diamond_parser_with_header(tmp_path: Path) -> None:
    tsv = tmp_path / "diamond_cluster.tsv"
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    uid_c = make_protein_uid(1, 0)
    with open(tsv, "w") as fh:
        fh.write("#centroid\tmember\n")
        fh.write(f"{uid_a}\t{uid_a}\n")
        fh.write(f"{uid_a}\t{uid_b}\n")
        fh.write(f"{uid_c}\t{uid_c}\n")

    out = tmp_path / "clusters.parquet"
    n_clusters = parse_diamond_cluster_tsv(cluster_tsv=tsv, out_parquet=out)
    assert n_clusters == 2

    table = pq.read_table(out)
    validate_schema(table, CLUSTERS_SCHEMA, "clusters.parquet")
    assert table.num_rows == 3
    assert table.column("is_centroid").to_pylist() == [True, False, True]
    assert table.column("cluster_id").to_pylist() == [0, 0, 1]
    assert table.column("genome_uid").to_pylist() == [0, 0, 1]


# ---------------------------------------------------------------------------
# CD-HIT .clstr parser
# ---------------------------------------------------------------------------


def test_cdhit_clstr_parser(tmp_path: Path) -> None:
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    uid_c = make_protein_uid(1, 0)

    clstr = tmp_path / "cdhit_out.clstr"
    with open(clstr, "w") as fh:
        fh.write(">Cluster 0\n")
        fh.write(f"0\t350aa, >{uid_a}... *\n")
        fh.write(f"1\t340aa, >{uid_b}... at 92.15%\n")
        fh.write(">Cluster 1\n")
        fh.write(f"0\t500aa, >{uid_c}... *\n")

    out = tmp_path / "clusters.parquet"
    n_clusters = parse_cdhit_clstr(clstr_path=clstr, out_parquet=out)
    assert n_clusters == 2

    table = pq.read_table(out)
    validate_schema(table, CLUSTERS_SCHEMA, "clusters.parquet")
    assert table.num_rows == 3

    assert table.column("cluster_id").to_pylist() == [0, 0, 1]
    assert table.column("is_centroid").to_pylist() == [True, False, True]
    # All alignment columns are null from the parser; engines.realign
    # populates them uniformly across engines after clustering.
    assert table.column("pct_identity_fwd").to_pylist() == [None, None, None]
    assert table.column("pct_identity_rev").to_pylist() == [None, None, None]
    assert table.column("member_coverage").to_pylist() == [None, None, None]
    assert table.column("rep_coverage").to_pylist() == [None, None, None]
    assert table.column("representative_uid").to_pylist() == [uid_a, uid_a, uid_c]


def test_cdhit_parser_handles_nucleotide_sizes(tmp_path: Path) -> None:
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    clstr = tmp_path / "cdhit_est.clstr"
    with open(clstr, "w") as fh:
        fh.write(">Cluster 0\n")
        fh.write(f"0\t1050nt, >{uid_a}... *\n")
        fh.write(f"1\t1020nt, >{uid_b}... at 88.0%\n")

    out = tmp_path / "clusters.parquet"
    n_clusters = parse_cdhit_clstr(clstr_path=clstr, out_parquet=out)
    assert n_clusters == 1
    table = pq.read_table(out)
    assert table.num_rows == 2


# ---------------------------------------------------------------------------
# usearch12 .uc parser
# ---------------------------------------------------------------------------


def test_uc_parser_s_and_h_rows(tmp_path: Path) -> None:
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    uid_c = make_protein_uid(1, 5)

    uc = tmp_path / "clusters.uc"
    with open(uc, "w") as fh:
        # Format: 10 tab-separated columns per SPEED.md §1
        fh.write(f"S\t0\t350\t*\t.\t*\t*\t*\t{uid_a}\t*\n")
        fh.write(f"H\t0\t340\t92.1\t.\t*\t*\t=\t{uid_b}\t{uid_a}\n")
        fh.write(f"S\t1\t500\t*\t.\t*\t*\t*\t{uid_c}\t*\n")
        fh.write(f"C\t0\t2\t*\t*\t*\t*\t*\t{uid_a}\t*\n")  # summary row - ignored
        fh.write(f"C\t1\t1\t*\t*\t*\t*\t*\t{uid_c}\t*\n")

    out = tmp_path / "clusters.parquet"
    n_clusters = parse_uc(uc_path=uc, out_parquet=out)
    assert n_clusters == 2

    table = pq.read_table(out)
    validate_schema(table, CLUSTERS_SCHEMA, "clusters.parquet")
    assert table.num_rows == 3  # S, H, S — C rows dropped

    assert table.column("is_centroid").to_pylist() == [True, False, True]
    assert table.column("cluster_id").to_pylist() == [0, 0, 1]
    assert table.column("representative_uid").to_pylist() == [uid_a, uid_a, uid_c]
    # Parser emits membership only; alignment metrics populated by
    # engines.realign after clustering.
    assert table.column("pct_identity_fwd").to_pylist() == [None, None, None]
    assert table.column("pct_identity_rev").to_pylist() == [None, None, None]


def test_uc_parser_skips_comment_and_malformed(tmp_path: Path) -> None:
    uid = make_protein_uid(0, 0)
    uc = tmp_path / "clusters.uc"
    with open(uc, "w") as fh:
        fh.write("# header comment\n")
        fh.write("malformed line without enough tabs\n")
        fh.write(f"S\t0\t100\t*\t.\t*\t*\t*\t{uid}\t*\n")

    out = tmp_path / "clusters.parquet"
    n_clusters = parse_uc(uc_path=uc, out_parquet=out)
    assert n_clusters == 1
    table = pq.read_table(out)
    assert table.num_rows == 1
    assert table.column("protein_uid").to_pylist() == [uid]
