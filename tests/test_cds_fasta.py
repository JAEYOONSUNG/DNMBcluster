"""Tests for the nucleotide CDS FASTA sidecar writer.

The sidecar must stay in lock-step with parse_folder / build_tables —
any divergence in the (genome_key, cds_key) resolution silently breaks
downstream codon-aware stages (HyPhy / codeml) because protein_uid
lookups fail.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("Bio")
pytest.importorskip("pyarrow")

import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from dnmbcluster.cds_fasta import (  # noqa: E402
    _resolve_genome_key_for_file,
    write_cds_nt_fasta,
)
from dnmbcluster.genbank import build_tables, parse_folder  # noqa: E402
from dnmbcluster.schemas import ID_MAP_SCHEMA  # noqa: E402

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "geobacillus10"


def _write_id_map(table: pa.Table, out: Path) -> None:
    pq.write_table(table, out)


@pytest.fixture(scope="module")
def parsed_fixture():
    """parse_folder + build_tables once for the fixture (slow)."""
    parsed = parse_folder(FIXTURE_DIR, level="protein", threads=1)
    id_map, gene_table, genome_meta = build_tables(parsed, level="protein")
    return parsed, id_map, gene_table, genome_meta


@pytest.mark.skipif(not FIXTURE_DIR.exists(), reason="geobacillus10 fixture missing")
def test_genome_key_resolver_matches_parse_folder(parsed_fixture) -> None:
    """_resolve_genome_key_for_file must produce the same key parse_folder assigns.

    This is the critical parity contract — a divergence would break every
    id_map lookup in the worker.
    """
    parsed, *_ = parsed_fixture
    for p in parsed:
        recovered = _resolve_genome_key_for_file(p.file_path)
        assert recovered == p.genome_key.key, (
            f"{p.file_path.name}: sidecar resolver={recovered!r} "
            f"parse_folder={p.genome_key.key!r}"
        )


@pytest.mark.skipif(not FIXTURE_DIR.exists(), reason="geobacillus10 fixture missing")
def test_write_cds_nt_fasta_serial(tmp_path, parsed_fixture) -> None:
    _, id_map, _, _ = parsed_fixture
    id_map_path = tmp_path / "id_map.parquet"
    _write_id_map(id_map, id_map_path)

    out = tmp_path / "cds.fna"
    n_written, n_missing = write_cds_nt_fasta(
        FIXTURE_DIR, id_map_path, out, threads=1,
    )

    assert out.exists() and out.stat().st_size > 0
    # Parity: every protein_uid in id_map must have a CDS emitted.
    assert n_missing == 0, f"{n_missing} protein_uids had no CDS match"
    assert n_written == id_map.num_rows

    # Headers must be exactly the protein_uids from id_map.
    uids_in_fasta: list[int] = []
    for line in out.read_text().splitlines():
        if line.startswith(">"):
            uids_in_fasta.append(int(line[1:]))
    assert set(uids_in_fasta) == set(id_map.column("protein_uid").to_pylist())


@pytest.mark.skipif(not FIXTURE_DIR.exists(), reason="geobacillus10 fixture missing")
def test_write_cds_nt_fasta_parallel_matches_serial(tmp_path, parsed_fixture) -> None:
    """Parallel and serial output must be byte-identical (determinism)."""
    _, id_map, _, _ = parsed_fixture
    id_map_path = tmp_path / "id_map.parquet"
    _write_id_map(id_map, id_map_path)

    serial_out = tmp_path / "cds_serial.fna"
    write_cds_nt_fasta(FIXTURE_DIR, id_map_path, serial_out, threads=1)

    parallel_out = tmp_path / "cds_parallel.fna"
    write_cds_nt_fasta(FIXTURE_DIR, id_map_path, parallel_out, threads=4)

    assert serial_out.read_bytes() == parallel_out.read_bytes()


@pytest.mark.skipif(not FIXTURE_DIR.exists(), reason="geobacillus10 fixture missing")
def test_cds_sequences_are_multiples_of_three(tmp_path, parsed_fixture) -> None:
    """CDS nt sequences must be codon-aligned (length % 3 == 0)."""
    _, id_map, _, _ = parsed_fixture
    id_map_path = tmp_path / "id_map.parquet"
    _write_id_map(id_map, id_map_path)

    out = tmp_path / "cds.fna"
    write_cds_nt_fasta(FIXTURE_DIR, id_map_path, out, threads=1)

    current_uid: int | None = None
    current_len = 0
    violations: list[tuple[int, int]] = []
    for line in out.read_text().splitlines():
        if line.startswith(">"):
            if current_uid is not None and current_len % 3 != 0:
                violations.append((current_uid, current_len))
            current_uid = int(line[1:])
            current_len = 0
        else:
            current_len += len(line)
    if current_uid is not None and current_len % 3 != 0:
        violations.append((current_uid, current_len))

    assert not violations, f"non-codon-aligned CDS: {violations[:5]}"


def test_empty_lookup_yields_empty_output(tmp_path) -> None:
    """id_map with no matching genome_key → empty FASTA, warnings logged."""
    if not FIXTURE_DIR.exists():
        pytest.skip("geobacillus10 fixture missing")

    # Build a fake id_map whose genome_key matches nothing on disk.
    empty_rows = {name.name: [] for name in ID_MAP_SCHEMA}
    # Need at least one row to exercise the grouping loop; make it a
    # genome_key that cannot possibly match any fixture file.
    empty_rows["protein_uid"] = [1]
    empty_rows["genome_uid"] = [0]
    empty_rows["gene_uid"] = [0]
    empty_rows["genome_key"] = ["NO_SUCH_GENOME"]
    empty_rows["genome_key_source"] = ["synthetic"]
    empty_rows["cds_key"] = ["no_such_cds"]
    empty_rows["cds_key_source"] = ["synthetic"]
    empty_rows["assembly_prefix"] = [None]
    empty_rows["assembly_version"] = [None]
    empty_rows["organism"] = [None]
    empty_rows["strain"] = [None]
    empty_rows["locus_tag"] = [None]
    empty_rows["protein_id"] = [None]
    empty_rows["gene"] = [None]
    empty_rows["product"] = [None]
    empty_rows["ec_number"] = [None]
    empty_rows["contig"] = ["x"]
    empty_rows["start"] = [0]
    empty_rows["end"] = [0]
    empty_rows["strand"] = [0]
    empty_rows["aa_length"] = [None]

    table = pa.table(empty_rows, schema=ID_MAP_SCHEMA)
    id_map_path = tmp_path / "id_map.parquet"
    _write_id_map(table, id_map_path)

    out = tmp_path / "cds.fna"
    n_written, n_missing = write_cds_nt_fasta(
        FIXTURE_DIR, id_map_path, out, threads=1,
    )

    assert out.exists()
    assert out.stat().st_size == 0
    assert n_written == 0
    assert n_missing == 1
