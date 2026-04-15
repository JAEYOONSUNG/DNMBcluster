"""Tests for identifier extraction and fallback chains."""
from __future__ import annotations

from pathlib import Path

import pytest

from dnmbcluster.ids import (
    _GENOME_MAX,
    GenomeKey,
    make_protein_uid,
    resolve_cds_key,
    resolve_genome_key,
    split_assembly,
    unpack_protein_uid,
)


# ---------------------------------------------------------------------------
# uint64 packing
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "genome_uid, gene_uid",
    [(0, 0), (1, 1), (255, 0), (0, 12345), (65_535, 4_294_967_295), (42, 1_000_000)],
)
def test_protein_uid_roundtrip(genome_uid: int, gene_uid: int) -> None:
    uid = make_protein_uid(genome_uid, gene_uid)
    assert unpack_protein_uid(uid) == (genome_uid, gene_uid)


def test_protein_uid_out_of_range_genome() -> None:
    with pytest.raises(ValueError, match="genome_uid"):
        make_protein_uid(_GENOME_MAX, 0)


def test_protein_uid_out_of_range_gene() -> None:
    with pytest.raises(ValueError, match="gene_uid"):
        make_protein_uid(0, 1 << 32)


# ---------------------------------------------------------------------------
# Assembly accession splitting
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "accession, expected",
    [
        ("GCF_000009785.1", ("GCF_000009785", ".1")),
        ("GCA_030376785.2", ("GCA_030376785", ".2")),
        ("GCF_000009785", ("GCF_000009785", None)),
        ("CP187452.1", (None, None)),
        ("random_text", (None, None)),
        ("", (None, None)),
    ],
)
def test_split_assembly(accession: str, expected: tuple) -> None:
    assert split_assembly(accession) == expected


# ---------------------------------------------------------------------------
# genome_key fallback chain
# ---------------------------------------------------------------------------


def test_genome_key_manifest_wins_over_everything() -> None:
    gk = resolve_genome_key(
        file_path=Path("GCF_000000001.1.gbff"),
        dblink=["Assembly: GCF_999999999.1"],
        locus_version="CP_SOMETHING.1",
        organism="Some organism",
        strain="X",
        manifest_name="MyCustomName",
    )
    assert gk.key == "MyCustomName"
    assert gk.source == "manifest"


def test_genome_key_dblink_over_filename() -> None:
    gk = resolve_genome_key(
        file_path=Path("unrelated_name.gbff"),
        dblink=["BioProject: PRJNA123", "Assembly: GCF_000009785.1"],
    )
    assert gk.key == "GCF_000009785.1"
    assert gk.source == "dblink"
    assert gk.assembly_prefix == "GCF_000009785"
    assert gk.assembly_version == ".1"


def test_genome_key_dblink_ignores_bioproject_only() -> None:
    gk = resolve_genome_key(
        file_path=Path("GCF_030376785.1.gbff"),
        dblink=["BioProject: PRJNA123", "BioSample: SAMN456"],
    )
    # Should fall through to filename
    assert gk.key == "GCF_030376785.1"
    assert gk.source == "filename_gcf"


def test_genome_key_filename_gcf_match() -> None:
    gk = resolve_genome_key(file_path=Path("GCF_030376785.1.gbff"), dblink=[])
    assert gk.key == "GCF_030376785.1"
    assert gk.source == "filename_gcf"
    assert gk.assembly_prefix == "GCF_030376785"


def test_genome_key_filename_gca() -> None:
    gk = resolve_genome_key(file_path=Path("GCA_000009785.2.gbff"), dblink=[])
    assert gk.key == "GCA_000009785.2"
    assert gk.source == "filename_gcf"


def test_genome_key_filename_no_gcf_falls_through_to_locus_version() -> None:
    gk = resolve_genome_key(
        file_path=Path("prokka_out.gbk"),
        dblink=[],
        locus_version="CP187452.1",
    )
    assert gk.key == "CP187452.1"
    assert gk.source == "locus_version"


def test_genome_key_organism_slug() -> None:
    gk = resolve_genome_key(
        file_path=Path("unknown.gbk"),
        dblink=[],
        organism="Geobacillus sp.",
        strain="GHH01",
    )
    assert gk.key == "Geobacillus_sp_GHH01"
    assert gk.source == "organism"


def test_genome_key_organism_dedupes_strain_in_organism() -> None:
    gk = resolve_genome_key(
        file_path=Path("unknown.gbk"),
        dblink=[],
        organism="Geobacillus sp. GHH01",
        strain="GHH01",
    )
    assert gk.key == "Geobacillus_sp_GHH01"
    assert gk.source == "organism"


def test_genome_key_filename_raw_last_resort() -> None:
    gk = resolve_genome_key(file_path=Path("whatever.gbk"), dblink=[])
    assert gk.key == "whatever"
    assert gk.source == "filename_raw"
    assert gk.assembly_prefix is None


# ---------------------------------------------------------------------------
# cds_key fallback chain
# ---------------------------------------------------------------------------


def test_cds_key_prefers_locus_tag() -> None:
    key, source = resolve_cds_key(
        locus_tag="GK0001", protein_id="WP_123.1", gene="dnaA",
        gene_ordinal=1, contig="chr1", start=100, end=1000, strand=1,
    )
    assert key == "GK0001"
    assert source == "locus_tag"


def test_cds_key_falls_back_to_protein_id() -> None:
    key, source = resolve_cds_key(
        locus_tag=None, protein_id="WP_012345678.1", gene="dnaA",
        gene_ordinal=1, contig="chr1", start=100, end=1000, strand=1,
    )
    assert key == "WP_012345678.1"
    assert source == "protein_id"


def test_cds_key_falls_back_to_gene_ordinal() -> None:
    key, source = resolve_cds_key(
        locus_tag=None, protein_id=None, gene="dnaA",
        gene_ordinal=3, contig="chr1", start=100, end=1000, strand=1,
    )
    assert key == "dnaA__3"
    assert source == "gene_ordinal"


def test_cds_key_falls_back_to_coordinates() -> None:
    key, source = resolve_cds_key(
        locus_tag=None, protein_id=None, gene=None, gene_ordinal=0,
        contig="CP187452.1", start=3421, end=4560, strand=-1,
    )
    assert key == "cds_CP187452.1_3421_4560_-"
    assert source == "coord"
