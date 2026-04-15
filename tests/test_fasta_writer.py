"""Unit tests for the FASTA writer."""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pyarrow")

import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from dnmbcluster.fasta import parse_fasta_headers, write_fasta  # noqa: E402
from dnmbcluster.ids import make_protein_uid  # noqa: E402
from dnmbcluster.schemas import gene_table_schema  # noqa: E402


def _write_gene_table(path: Path, rows: list[dict]) -> None:
    columns = {
        "protein_uid": [r["protein_uid"] for r in rows],
        "genome_uid": [r["genome_uid"] for r in rows],
        "length": [r["length"] for r in rows],
        "translation": [r["translation"] for r in rows],
    }
    table = pa.table(columns, schema=gene_table_schema("protein"))
    pq.write_table(table, path)


def test_write_fasta_basic(tmp_path: Path) -> None:
    gene_table = tmp_path / "gene_table.parquet"
    uids = [make_protein_uid(0, i) for i in range(3)]
    _write_gene_table(
        gene_table,
        [
            {"protein_uid": uids[0], "genome_uid": 0, "length": 4, "translation": "MKLA"},
            {"protein_uid": uids[1], "genome_uid": 0, "length": 6, "translation": "MHELLO"},
            {"protein_uid": uids[2], "genome_uid": 0, "length": 3, "translation": "MVA"},
        ],
    )

    fasta = tmp_path / "out.faa"
    n_seq, n_bytes = write_fasta(gene_table, fasta, level="protein")

    assert n_seq == 3
    assert n_bytes > 0
    assert fasta.exists()

    headers = parse_fasta_headers(fasta)
    # Sorted by length desc: MHELLO (6), MKLA (4), MVA (3)
    assert headers == [uids[1], uids[0], uids[2]]


def test_write_fasta_preserves_sequences(tmp_path: Path) -> None:
    gene_table = tmp_path / "gene_table.parquet"
    uid = make_protein_uid(0, 0)
    long_seq = "M" + ("A" * 250)
    _write_gene_table(
        gene_table,
        [{"protein_uid": uid, "genome_uid": 0, "length": len(long_seq), "translation": long_seq}],
    )

    fasta = tmp_path / "long.faa"
    write_fasta(gene_table, fasta, level="protein")

    content = fasta.read_text()
    # Header + wrapped body; concatenating non-header lines reconstructs the seq.
    lines = content.strip().split("\n")
    assert lines[0] == f">{uid}"
    reconstructed = "".join(lines[1:])
    assert reconstructed == long_seq
