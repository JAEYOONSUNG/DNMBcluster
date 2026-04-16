"""Unit + light end-to-end tests for the shared realignment stage.

Replaces the obsolete ``test_bidirectional_alignment_join`` that lived
inside the MMseqs2 engine before the refactor in commit 41996c8.

Tests split into three tiers:

1. **Pure helpers** — ``_write_subset_fasta``, ``_load_alignments``,
   ``_build_search_cmd``. No subprocess, no Biopython needed.
2. **All-singleton short-circuit** — ``_fill_centroids_only`` exercised
   via the public entry point with a clusters table that has no
   non-centroid rows. No mmseqs invoked.
3. **End-to-end** — full ``populate_alignment_metrics`` against a tiny
   synthetic FASTA. Skipped if ``mmseqs`` is not on PATH (CI doesn't
   install bioconda binaries).
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

pytest.importorskip("pyarrow")

import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from dnmbcluster.engines.realign import (  # noqa: E402
    _build_search_cmd,
    _load_alignments,
    _write_subset_fasta,
    populate_alignment_metrics,
)
from dnmbcluster.ids import make_protein_uid  # noqa: E402
from dnmbcluster.schemas import CLUSTERS_SCHEMA  # noqa: E402


# ---------------------------------------------------------------------------
# Fixtures: tiny FASTA + clusters.parquet builders
# ---------------------------------------------------------------------------


def _write_fasta(records: dict[int, str], path: Path) -> Path:
    with open(path, "w") as fh:
        for uid, seq in records.items():
            fh.write(f">{uid}\n{seq}\n")
    return path


def _write_membership_only(rows: list[tuple[int, int, int, bool]], path: Path) -> Path:
    """Write a membership-only clusters.parquet from
    ``(protein_uid, cluster_id, representative_uid, is_centroid)`` tuples."""
    n = len(rows)
    table = pa.table(
        {
            "protein_uid":        pa.array([r[0] for r in rows], type=pa.uint64()),
            "genome_uid":         pa.array([(r[0] >> 48) & 0xFFFF for r in rows], type=pa.uint16()),
            "cluster_id":         pa.array([r[1] for r in rows], type=pa.uint32()),
            "representative_uid": pa.array([r[2] for r in rows], type=pa.uint64()),
            "is_centroid":        pa.array([r[3] for r in rows], type=pa.bool_()),
            "pct_identity_fwd":   pa.array([None] * n, type=pa.float32()),
            "pct_identity_rev":   pa.array([None] * n, type=pa.float32()),
            "member_coverage":    pa.array([None] * n, type=pa.float32()),
            "rep_coverage":       pa.array([None] * n, type=pa.float32()),
            "alignment_length":   pa.array([None] * n, type=pa.uint32()),
        },
        schema=CLUSTERS_SCHEMA,
    )
    pq.write_table(table, path)
    return path


# ---------------------------------------------------------------------------
# Tier 1 — pure helpers
# ---------------------------------------------------------------------------


def test_write_subset_fasta_filters_by_uid(tmp_path: Path) -> None:
    src = tmp_path / "all.faa"
    _write_fasta(
        {1001: "MKVLW", 1002: "MAAAY", 1003: "MGGGS"},
        src,
    )
    dst = tmp_path / "subset.faa"
    _write_subset_fasta(src, {1001, 1003}, dst)

    text = dst.read_text()
    assert ">1001" in text
    assert ">1003" in text
    assert ">1002" not in text
    assert "MKVLW" in text
    assert "MAAAY" not in text


def test_write_subset_fasta_handles_multiline_sequence(tmp_path: Path) -> None:
    src = tmp_path / "wrapped.faa"
    src.write_text(">2001\nMKVL\nWPQR\n>2002\nAAAA\n")
    dst = tmp_path / "subset.faa"
    _write_subset_fasta(src, {2001}, dst)
    text = dst.read_text()
    assert ">2001" in text
    # Both wrapped lines of 2001 must be preserved
    assert "MKVL" in text and "WPQR" in text
    assert ">2002" not in text


def test_load_alignments_parses_tsv(tmp_path: Path) -> None:
    tsv = tmp_path / "fwd.tsv"
    # query, target, pident, qcov, tcov, alnlen
    tsv.write_text(
        "100\t200\t91.5\t0.97\t0.95\t340\n"
        "101\t201\t100.0\t1.0\t1.0\t500\n"
    )
    lookup = _load_alignments(tsv)
    assert (100, 200) in lookup
    assert (101, 201) in lookup
    pid, qcov, tcov, aln = lookup[(100, 200)]
    assert pid == pytest.approx(91.5)
    assert qcov == pytest.approx(0.97)
    assert tcov == pytest.approx(0.95)
    assert aln == 340


def test_load_alignments_missing_file_returns_empty(tmp_path: Path) -> None:
    assert _load_alignments(tmp_path / "nope.tsv") == {}


def test_build_search_cmd_threshold_free(tmp_path: Path) -> None:
    cmd = _build_search_cmd(
        mmseqs_binary="mmseqs",
        query_fasta=tmp_path / "q.faa",
        target_fasta=tmp_path / "t.faa",
        out_tsv=tmp_path / "out.tsv",
        tmpdir=tmp_path / "tmp",
        threads=2,
        level="protein",
        max_seqs=1,
    )
    # Threshold-free is the whole point: identity / coverage / e-value all wide open.
    assert "--min-seq-id" in cmd and cmd[cmd.index("--min-seq-id") + 1] == "0"
    assert "-c" in cmd and cmd[cmd.index("-c") + 1] == "0"
    assert "-e" in cmd and cmd[cmd.index("-e") + 1] == "1e10"
    assert "--max-seqs" in cmd and cmd[cmd.index("--max-seqs") + 1] == "1"
    assert "--search-type" not in cmd  # protein mode


def test_build_search_cmd_nucleotide_adds_search_type(tmp_path: Path) -> None:
    cmd = _build_search_cmd(
        mmseqs_binary="mmseqs",
        query_fasta=tmp_path / "q.fna",
        target_fasta=tmp_path / "t.fna",
        out_tsv=tmp_path / "out.tsv",
        tmpdir=tmp_path / "tmp",
        threads=1,
        level="nucleotide",
        max_seqs=10,
    )
    assert "--search-type" in cmd
    assert cmd[cmd.index("--search-type") + 1] == "3"


# ---------------------------------------------------------------------------
# Tier 2 — all-singleton short-circuit (no mmseqs invoked)
# ---------------------------------------------------------------------------


def test_all_singletons_skip_mmseqs(tmp_path: Path) -> None:
    """When every row is its own centroid, no alignment work is needed."""
    if shutil.which("mmseqs") is None:
        # The which() check inside populate_alignment_metrics still runs
        # before the short-circuit, so we need mmseqs visible in PATH for
        # this test even though the binary is never invoked.
        pytest.skip("mmseqs not on PATH; cannot reach short-circuit branch")

    uid_a = make_protein_uid(0, 0)
    uid_b = make_protein_uid(1, 0)
    uid_c = make_protein_uid(2, 0)

    fasta = _write_fasta(
        {uid_a: "MKVLW", uid_b: "MAAAY", uid_c: "MGGGS"},
        tmp_path / "input.faa",
    )
    clusters = _write_membership_only(
        [
            (uid_a, 0, uid_a, True),
            (uid_b, 1, uid_b, True),
            (uid_c, 2, uid_c, True),
        ],
        tmp_path / "clusters.parquet",
    )

    result = populate_alignment_metrics(
        clusters_parquet=clusters,
        input_fasta=fasta,
        out_dir=tmp_path / "realign",
        threads=1,
    )

    assert result.n_rows == 3
    assert result.n_centroids == 3
    assert result.n_aligned == 0
    assert result.n_analytical == 3
    assert result.n_missing == 0
    assert result.fwd_tsv is None  # short-circuit didn't write any TSV

    table = pq.read_table(clusters)
    assert table.column("pct_identity_fwd").to_pylist() == [100.0, 100.0, 100.0]
    assert table.column("pct_identity_rev").to_pylist() == [100.0, 100.0, 100.0]
    assert table.column("member_coverage").to_pylist() == [1.0, 1.0, 1.0]
    assert table.column("rep_coverage").to_pylist() == [1.0, 1.0, 1.0]


# ---------------------------------------------------------------------------
# Tier 3 — end-to-end with real MMseqs2 (skipped without binary)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs not on PATH")
def test_end_to_end_populates_all_columns(tmp_path: Path) -> None:
    """A two-cluster fixture with one near-identical pair should produce
    bidirectional identity, coverages, and an alignment length for the
    non-centroid row, plus analytical 100/100/1.0/1.0 for centroids."""
    rep1 = make_protein_uid(0, 0)
    mem1 = make_protein_uid(0, 1)
    rep2 = make_protein_uid(1, 0)

    # Two ~150-aa sequences to clear the MMseqs2 k-mer prefilter.
    seq_rep1 = (
        "MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQTLSEKVQ"
        "GAQDGTRRDPRWAHIVKKILVPHHEWGNRLERMVGYHTNWQADRTLDPMHKRLNYIA"
    )
    # Member differs by 5 residues out of ~140 — enough to exercise pident < 100.
    seq_mem1 = (
        "MKVLWAALMVTFLAGCQAKVEQAVETEPEPELRQETEWQSGQRWELALGRFWDYLRWVQTLSEKVQ"
        "GAQDGSRRDPRWAHIVKKILVPHHEWGNRLERMVGYHTNWQADRTLDPMHKRLDYIA"
    )
    # Unrelated rep for cluster 2.
    seq_rep2 = (
        "MAAAYTLTFAGCSEKVELRPELDGEPYRWAQHTLGYRWDDYWLRLTLGQRWELALGRFLDDFLRWV"
        "QTLSEKVQRRDPHHEWGNRLERMVGYHTNWQADRTLDPMHKRLNYIAQDGSRRDPR"
    )

    fasta = _write_fasta(
        {rep1: seq_rep1, mem1: seq_mem1, rep2: seq_rep2},
        tmp_path / "input.faa",
    )
    clusters = _write_membership_only(
        [
            (rep1, 0, rep1, True),
            (mem1, 0, rep1, False),
            (rep2, 1, rep2, True),
        ],
        tmp_path / "clusters.parquet",
    )

    result = populate_alignment_metrics(
        clusters_parquet=clusters,
        input_fasta=fasta,
        out_dir=tmp_path / "realign",
        threads=1,
    )

    assert result.n_rows == 3
    assert result.n_centroids == 2
    assert result.n_analytical == 2
    assert result.n_missing == 0
    assert result.n_aligned >= 1
    assert result.fwd_tsv is not None and result.fwd_tsv.exists()
    assert result.rev_tsv is not None and result.rev_tsv.exists()

    table = pq.read_table(clusters)
    pid_fwd = table.column("pct_identity_fwd").to_pylist()
    pid_rev = table.column("pct_identity_rev").to_pylist()
    member_cov = table.column("member_coverage").to_pylist()
    rep_cov = table.column("rep_coverage").to_pylist()

    # Centroid rows: analytical fill
    assert pid_fwd[0] == 100.0 and pid_rev[0] == 100.0
    assert pid_fwd[2] == 100.0 and pid_rev[2] == 100.0
    assert member_cov[0] == 1.0 and rep_cov[0] == 1.0

    # Non-centroid row: real alignment metrics, both directions populated,
    # identity strictly below 100 (we mutated 4 residues).
    assert pid_fwd[1] is not None
    assert pid_rev[1] is not None
    assert 80.0 < pid_fwd[1] < 100.0
    assert 80.0 < pid_rev[1] < 100.0
    assert member_cov[1] is not None and member_cov[1] > 0.9
    assert rep_cov[1] is not None and rep_cov[1] > 0.9
    assert table.column("alignment_length").to_pylist()[1] is not None


@pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs not on PATH")
def test_derives_rep_fasta_when_none_passed(tmp_path: Path) -> None:
    """Engines like DIAMOND don't write a rep_seq.fasta. The shared
    enrichment must derive it from clusters.parquet + input FASTA."""
    rep = make_protein_uid(0, 0)
    mem = make_protein_uid(0, 1)
    seq_rep = (
        "MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQTLSEKVQ"
        "GAQDGTRRDPRWAHIVKKILVPHHEWGNRLERMVGYHTNWQADRTLDPMHKRLNYIA"
    )
    seq_mem = seq_rep.replace("MKVLW", "MKVLY")

    fasta = _write_fasta({rep: seq_rep, mem: seq_mem}, tmp_path / "input.faa")
    clusters = _write_membership_only(
        [(rep, 0, rep, True), (mem, 0, rep, False)],
        tmp_path / "clusters.parquet",
    )

    result = populate_alignment_metrics(
        clusters_parquet=clusters,
        input_fasta=fasta,
        out_dir=tmp_path / "realign",
        threads=1,
        representatives_fasta=None,  # explicit: derive it
    )

    derived = tmp_path / "realign" / "derived_rep_seq.fasta"
    assert derived.exists()
    assert ">{}".format(rep) in derived.read_text()
    assert result.n_missing == 0


def test_missing_mmseqs_binary_raises(tmp_path: Path) -> None:
    fasta = _write_fasta({1: "MKV"}, tmp_path / "in.faa")
    clusters = _write_membership_only(
        [(1, 0, 1, True)],
        tmp_path / "clusters.parquet",
    )
    with pytest.raises(RuntimeError, match="not found on PATH"):
        populate_alignment_metrics(
            clusters_parquet=clusters,
            input_fasta=fasta,
            out_dir=tmp_path / "realign",
            mmseqs_binary="mmseqs_definitely_not_a_real_binary_xyz123",
        )
