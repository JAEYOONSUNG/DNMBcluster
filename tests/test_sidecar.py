"""Tests for engine-native sidecar parquets.

Two sidecars are written today:

- ``usearch12_native.parquet`` — CIGAR + native pct identity + strand
  for every H row in the .uc file.
- ``cdhit_native.parquet`` — alignment positions (when ``-p 1`` was set)
  + native pct identity for every non-centroid member of every cluster.

DIAMOND has no sidecar by design — see ``engines/sidecar.py`` docstring.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pyarrow")

import pyarrow.parquet as pq  # noqa: E402

from dnmbcluster.engines.cdhit import parse_cdhit_clstr  # noqa: E402
from dnmbcluster.engines.sidecar import (  # noqa: E402
    CDHIT_NATIVE_SCHEMA,
    USEARCH12_NATIVE_SCHEMA,
)
from dnmbcluster.engines.usearch12 import parse_uc  # noqa: E402
from dnmbcluster.ids import make_protein_uid  # noqa: E402
from dnmbcluster.schemas import validate_schema  # noqa: E402


# ---------------------------------------------------------------------------
# usearch12 sidecar
# ---------------------------------------------------------------------------


def test_usearch12_sidecar_captures_cigar_and_strand(tmp_path: Path) -> None:
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    uid_c = make_protein_uid(1, 5)
    uid_d = make_protein_uid(1, 6)

    uc = tmp_path / "clusters.uc"
    with open(uc, "w") as fh:
        # cluster 0: centroid + one H row with CIGAR "340M"
        fh.write(f"S\t0\t350\t*\t.\t*\t*\t*\t{uid_a}\t*\n")
        fh.write(f"H\t0\t340\t92.1\t+\t*\t*\t340M\t{uid_b}\t{uid_a}\n")
        # cluster 1: centroid + one H row with strand "-" and a richer CIGAR
        fh.write(f"S\t1\t500\t*\t.\t*\t*\t*\t{uid_c}\t*\n")
        fh.write(f"H\t1\t490\t88.0\t-\t*\t*\t100M2D388M\t{uid_d}\t{uid_c}\n")

    clusters_parquet = tmp_path / "clusters.parquet"
    sidecar_parquet = tmp_path / "usearch12_native.parquet"
    n_clusters = parse_uc(
        uc_path=uc,
        out_parquet=clusters_parquet,
        sidecar_parquet=sidecar_parquet,
    )
    assert n_clusters == 2

    sidecar = pq.read_table(sidecar_parquet)
    validate_schema(sidecar, USEARCH12_NATIVE_SCHEMA, "usearch12_native.parquet")
    # Two H rows → two sidecar rows. S rows are intentionally excluded.
    assert sidecar.num_rows == 2

    rows = sidecar.to_pydict()
    assert rows["protein_uid"] == [uid_b, uid_d]
    assert rows["representative_uid"] == [uid_a, uid_c]
    assert rows["native_pct_identity"] == pytest.approx([92.1, 88.0])
    assert rows["strand"] == ["+", "-"]
    assert rows["cigar"] == ["340M", "100M2D388M"]


def test_usearch12_sidecar_handles_missing_cigar(tmp_path: Path) -> None:
    """``*`` in column 7 means CIGAR was not stored. Becomes null."""
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    uc = tmp_path / "clusters.uc"
    with open(uc, "w") as fh:
        fh.write(f"S\t0\t350\t*\t.\t*\t*\t*\t{uid_a}\t*\n")
        fh.write(f"H\t0\t340\t92.1\t.\t*\t*\t*\t{uid_b}\t{uid_a}\n")

    clusters_parquet = tmp_path / "clusters.parquet"
    sidecar_parquet = tmp_path / "usearch12_native.parquet"
    parse_uc(
        uc_path=uc,
        out_parquet=clusters_parquet,
        sidecar_parquet=sidecar_parquet,
    )

    sidecar = pq.read_table(sidecar_parquet)
    rows = sidecar.to_pydict()
    assert rows["cigar"] == [None]
    assert rows["strand"] == [None]


def test_usearch12_sidecar_optional(tmp_path: Path) -> None:
    """Caller can pass sidecar_parquet=None to skip the sidecar entirely."""
    uid_a = make_protein_uid(0, 0)
    uc = tmp_path / "clusters.uc"
    with open(uc, "w") as fh:
        fh.write(f"S\t0\t350\t*\t.\t*\t*\t*\t{uid_a}\t*\n")

    clusters_parquet = tmp_path / "clusters.parquet"
    parse_uc(uc_path=uc, out_parquet=clusters_parquet, sidecar_parquet=None)
    assert clusters_parquet.exists()
    assert not (tmp_path / "usearch12_native.parquet").exists()


# ---------------------------------------------------------------------------
# CD-HIT sidecar — default flavor (no -p 1)
# ---------------------------------------------------------------------------


def test_cdhit_sidecar_default_flavor_no_positions(tmp_path: Path) -> None:
    """Without ``-p 1``, the sidecar still captures native pct identity
    but leaves position columns null."""
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)

    clstr = tmp_path / "cdhit_out.clstr"
    with open(clstr, "w") as fh:
        fh.write(">Cluster 0\n")
        fh.write(f"0\t350aa, >{uid_a}... *\n")
        fh.write(f"1\t340aa, >{uid_b}... at 92.15%\n")

    clusters_parquet = tmp_path / "clusters.parquet"
    sidecar_parquet = tmp_path / "cdhit_native.parquet"
    parse_cdhit_clstr(
        clstr_path=clstr,
        out_parquet=clusters_parquet,
        sidecar_parquet=sidecar_parquet,
    )

    sidecar = pq.read_table(sidecar_parquet)
    validate_schema(sidecar, CDHIT_NATIVE_SCHEMA, "cdhit_native.parquet")
    assert sidecar.num_rows == 1
    rows = sidecar.to_pydict()
    assert rows["protein_uid"] == [uid_b]
    assert rows["representative_uid"] == [uid_a]
    assert rows["native_pct_identity"] == pytest.approx([92.15])
    assert rows["query_start"] == [None]
    assert rows["query_end"] == [None]
    assert rows["target_start"] == [None]
    assert rows["target_end"] == [None]
    assert rows["strand"] == [None]


# ---------------------------------------------------------------------------
# CD-HIT sidecar — -p 1 protein flavor
# ---------------------------------------------------------------------------


def test_cdhit_sidecar_p1_protein_positions(tmp_path: Path) -> None:
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)

    clstr = tmp_path / "cdhit_out.clstr"
    with open(clstr, "w") as fh:
        fh.write(">Cluster 0\n")
        fh.write(f"0\t350aa, >{uid_a}... *\n")
        fh.write(f"1\t340aa, >{uid_b}... at 1:336:5:340/92.15%\n")

    clusters_parquet = tmp_path / "clusters.parquet"
    sidecar_parquet = tmp_path / "cdhit_native.parquet"
    parse_cdhit_clstr(
        clstr_path=clstr,
        out_parquet=clusters_parquet,
        sidecar_parquet=sidecar_parquet,
    )

    sidecar = pq.read_table(sidecar_parquet)
    rows = sidecar.to_pydict()
    assert rows["protein_uid"] == [uid_b]
    assert rows["query_start"] == [1]
    assert rows["query_end"] == [336]
    assert rows["target_start"] == [5]
    assert rows["target_end"] == [340]
    assert rows["strand"] == [None]  # protein has no strand
    assert rows["native_pct_identity"] == pytest.approx([92.15])

    # Verify the canonical clusters.parquet still parses correctly with
    # the new positional format.
    table = pq.read_table(clusters_parquet)
    assert table.num_rows == 2
    assert table.column("is_centroid").to_pylist() == [True, False]


# ---------------------------------------------------------------------------
# CD-HIT sidecar — -p 1 nucleotide flavor (with strand)
# ---------------------------------------------------------------------------


def test_cdhit_sidecar_p1_nucleotide_with_strand(tmp_path: Path) -> None:
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)
    uid_c = make_protein_uid(1, 0)
    uid_d = make_protein_uid(1, 1)

    clstr = tmp_path / "cdhit_est.clstr"
    with open(clstr, "w") as fh:
        fh.write(">Cluster 0\n")
        fh.write(f"0\t1050nt, >{uid_a}... *\n")
        fh.write(f"1\t1020nt, >{uid_b}... at 1:1020:5:1024/+/92.15%\n")
        fh.write(">Cluster 1\n")
        fh.write(f"0\t800nt, >{uid_c}... *\n")
        fh.write(f"1\t790nt, >{uid_d}... at 1:790:11:800/-/88.0%\n")

    clusters_parquet = tmp_path / "clusters.parquet"
    sidecar_parquet = tmp_path / "cdhit_native.parquet"
    n_clusters = parse_cdhit_clstr(
        clstr_path=clstr,
        out_parquet=clusters_parquet,
        sidecar_parquet=sidecar_parquet,
    )
    assert n_clusters == 2

    sidecar = pq.read_table(sidecar_parquet)
    rows = sidecar.to_pydict()
    assert rows["protein_uid"] == [uid_b, uid_d]
    assert rows["strand"] == ["+", "-"]
    assert rows["query_start"] == [1, 1]
    assert rows["query_end"] == [1020, 790]
    assert rows["target_start"] == [5, 11]
    assert rows["target_end"] == [1024, 800]
    assert rows["native_pct_identity"] == pytest.approx([92.15, 88.0])


def test_cdhit_sidecar_default_flavor_nucleotide_strand(tmp_path: Path) -> None:
    """Default cd-hit-est without -p 1 prints strand but no positions:
    ``at +/92.15%``. Sidecar should pick up the strand."""
    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(0, 1)

    clstr = tmp_path / "cdhit_est.clstr"
    with open(clstr, "w") as fh:
        fh.write(">Cluster 0\n")
        fh.write(f"0\t1050nt, >{uid_a}... *\n")
        fh.write(f"1\t1020nt, >{uid_b}... at +/92.15%\n")

    clusters_parquet = tmp_path / "clusters.parquet"
    sidecar_parquet = tmp_path / "cdhit_native.parquet"
    parse_cdhit_clstr(
        clstr_path=clstr,
        out_parquet=clusters_parquet,
        sidecar_parquet=sidecar_parquet,
    )

    sidecar = pq.read_table(sidecar_parquet)
    rows = sidecar.to_pydict()
    assert rows["strand"] == ["+"]
    assert rows["query_start"] == [None]  # no positional info
    assert rows["native_pct_identity"] == pytest.approx([92.15])


def test_cdhit_sidecar_optional(tmp_path: Path) -> None:
    """Caller can pass sidecar_parquet=None to skip the sidecar entirely."""
    uid_a = make_protein_uid(0, 0)
    clstr = tmp_path / "cdhit_out.clstr"
    with open(clstr, "w") as fh:
        fh.write(">Cluster 0\n")
        fh.write(f"0\t350aa, >{uid_a}... *\n")

    clusters_parquet = tmp_path / "clusters.parquet"
    parse_cdhit_clstr(
        clstr_path=clstr,
        out_parquet=clusters_parquet,
        sidecar_parquet=None,
    )
    assert clusters_parquet.exists()
    assert not (tmp_path / "cdhit_native.parquet").exists()
