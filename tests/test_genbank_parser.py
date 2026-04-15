"""Integration test for the GenBank parser on the Geobacillus-10 fixture.

Skipped if Biopython or PyArrow aren't installed (CI runs them via
conda env; local dev may not have them yet).
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("Bio")
pytest.importorskip("pyarrow")

from dnmbcluster.genbank import build_tables, parse_folder  # noqa: E402
from dnmbcluster.ids import unpack_protein_uid  # noqa: E402

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "geobacillus10"


@pytest.mark.skipif(not FIXTURE_DIR.exists(), reason="geobacillus10 fixture missing")
def test_parse_geobacillus10_protein() -> None:
    parsed = parse_folder(FIXTURE_DIR, level="protein", threads=1)
    assert len(parsed) == 10, "expected exactly 10 Geobacillus fixture files"

    for p in parsed:
        assert p.n_cds > 0, f"{p.file_path.name}: no CDS extracted"
        assert p.genome_key.key, f"{p.file_path.name}: empty genome_key"
        assert p.file_sha256, f"{p.file_path.name}: missing sha256"

    id_map, gene_table, genome_meta = build_tables(parsed, level="protein")

    assert id_map.num_rows == sum(p.n_cds for p in parsed)
    assert gene_table.num_rows == id_map.num_rows
    assert genome_meta.num_rows == 10

    # protein_uid must be globally unique
    uids = id_map.column("protein_uid").to_pylist()
    assert len(set(uids)) == len(uids), "duplicate protein_uid detected"

    # Composite (genome_key, cds_key) must be globally unique
    composites = list(
        zip(
            id_map.column("genome_key").to_pylist(),
            id_map.column("cds_key").to_pylist(),
        )
    )
    assert len(set(composites)) == len(composites), "duplicate (genome_key, cds_key)"

    # genome_uid from protein_uid must match the id_map column
    genome_uids = id_map.column("genome_uid").to_pylist()
    for uid, expected in zip(uids, genome_uids):
        unpacked_genome, _ = unpack_protein_uid(uid)
        assert unpacked_genome == expected

    # Every gene_table row must have a non-empty translation
    translations = gene_table.column("translation").to_pylist()
    assert all(t for t in translations), "empty translation in gene_table"


@pytest.mark.skipif(not FIXTURE_DIR.exists(), reason="geobacillus10 fixture missing")
def test_genbank_source_attribution() -> None:
    """Every CDS records which fallback rung produced its key."""
    parsed = parse_folder(FIXTURE_DIR, level="protein", threads=1)
    id_map, _, _ = build_tables(parsed, level="protein")
    from dnmbcluster.schemas import CDS_KEY_SOURCES, GENOME_KEY_SOURCES

    cds_sources = set(id_map.column("cds_key_source").to_pylist())
    assert cds_sources.issubset(set(CDS_KEY_SOURCES)), cds_sources

    genome_sources = set(id_map.column("genome_key_source").to_pylist())
    assert genome_sources.issubset(set(GENOME_KEY_SOURCES)), genome_sources
