"""GenBank folder parser — walks input, extracts CDS features per file.

M1 implementation uses Biopython's ``SeqIO.parse`` for correctness.
Speed optimization (custom byte-level parser) is deferred to M6 if
benchmarks warrant it — SPEED.md Section 3.

Parallelizes across input files with ``ProcessPoolExecutor``; each file
is independent (SPEED.md Section 7, Level 1 parallelism).
"""
from __future__ import annotations

import hashlib
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, Optional

import pyarrow as pa

from .ids import GenomeKey, make_protein_uid, resolve_cds_key, resolve_genome_key
from .schemas import (
    GENOME_META_SCHEMA,
    ID_MAP_SCHEMA,
    gene_table_schema,
    validate_schema,
)

log = logging.getLogger(__name__)

GENBANK_SUFFIXES = (".gb", ".gbk", ".gbff")
SUPPORTED_LEVELS = ("protein", "nucleotide")


# ---------------------------------------------------------------------------
# Filesystem helpers
# ---------------------------------------------------------------------------


def find_genbank_files(input_dir: Path) -> list[Path]:
    """Return sorted GenBank files directly under ``input_dir``.

    Only looks at the top level (not recursive) — nested directories
    usually mean "multiple datasets" and should be explicit.
    """
    files: list[Path] = []
    for path in input_dir.iterdir():
        if path.is_file() and path.suffix.lower() in GENBANK_SUFFIXES:
            files.append(path)
    return sorted(files, key=lambda p: p.name)


def _sha256_and_size(path: Path, chunk: int = 1 << 20) -> tuple[str, int]:
    """Single-pass SHA256 + byte count of a file."""
    h = hashlib.sha256()
    total = 0
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
            total += len(block)
    return h.hexdigest(), total


# ---------------------------------------------------------------------------
# Per-file parse result
# ---------------------------------------------------------------------------


@dataclass
class ParsedGenome:
    """All information extracted from one GenBank file.

    `rows` are CDS-level dicts; protein_uid is not yet assigned — it
    requires the run-global genome_uid which is only known at aggregation.
    """

    file_path: Path
    file_sha256: str
    file_bytes: int
    genome_key: GenomeKey
    organism: Optional[str]
    strain: Optional[str]
    definition: Optional[str]
    n_records: int
    n_cds: int
    n_skipped_pseudogenes: int
    total_length: int
    gc_count: int
    rows: list[dict[str, Any]] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


def _iter_records(file_path: Path) -> Iterator[Any]:
    """Yield Biopython SeqRecord objects from a GenBank file.

    Biopython import is deferred so the CLI can import this module
    without failing when Biopython is absent (e.g., for `--help`).
    """
    from Bio import SeqIO  # type: ignore[import-untyped]

    with open(file_path) as handle:
        yield from SeqIO.parse(handle, "genbank")


def parse_genbank_file(
    file_path: Path,
    *,
    level: str = "protein",
    manifest_name: Optional[str] = None,
) -> ParsedGenome:
    """Parse one GenBank file end-to-end.

    Extracts CDS features only — never `gene`, `tRNA`, `rRNA`, or
    `misc_feature`. For ``level='protein'`` uses ``/translation`` and
    skips pseudogenes. For ``level='nucleotide'`` extracts the genomic
    span and reverse-complements on the minus strand.
    """
    if level not in SUPPORTED_LEVELS:
        raise ValueError(f"unknown --level: {level!r}. Expected one of {SUPPORTED_LEVELS}")

    file_sha256, file_bytes = _sha256_and_size(file_path)

    rows: list[dict[str, Any]] = []
    n_records = 0
    n_cds = 0
    n_skipped_pseudogenes = 0
    total_length = 0
    gc_count = 0

    organism: Optional[str] = None
    strain: Optional[str] = None
    definition: Optional[str] = None
    dblink_entries: list[str] = []
    first_locus_version: Optional[str] = None

    # Per-file counters for deterministic fallback identifiers
    gene_ordinals: dict[str, int] = {}
    locus_tag_occurrences: dict[str, int] = {}

    for record in _iter_records(file_path):
        n_records += 1

        seq_len = len(record.seq) if record.seq is not None else 0
        total_length += seq_len
        if seq_len > 0:
            seq_str = str(record.seq).upper()
            gc_count += seq_str.count("G") + seq_str.count("C")

        if first_locus_version is None:
            first_locus_version = record.id or record.name or None

        if definition is None:
            desc = getattr(record, "description", None)
            # Biopython fills description with a placeholder "." when the
            # DEFINITION line is empty; treat that as missing so reports
            # fall back to organism/strain cleanly.
            if desc and desc.strip() not in ("", "."):
                definition = desc.strip().rstrip(".")

        annotations = record.annotations or {}
        dblink_entries.extend(annotations.get("dblink") or [])
        if organism is None:
            organism = annotations.get("organism") or None

        for feat in record.features:
            if feat.type == "source" and strain is None:
                strain_vals = feat.qualifiers.get("strain")
                if strain_vals:
                    strain = strain_vals[0]
                break  # only need first source

        contig = record.id or record.name or file_path.stem

        for feat in record.features:
            if feat.type != "CDS":
                continue

            quals = feat.qualifiers
            is_pseudogene = "pseudo" in quals or "pseudogene" in quals

            if level == "protein":
                translation_vals = quals.get("translation")
                if not translation_vals or is_pseudogene:
                    n_skipped_pseudogenes += 1
                    continue
                sequence = "".join(translation_vals[0].split())
                if not sequence:
                    n_skipped_pseudogenes += 1
                    continue
            else:  # nucleotide
                if record.seq is None:
                    n_skipped_pseudogenes += 1
                    continue
                try:
                    nt = feat.extract(record.seq)
                except Exception as exc:  # malformed location
                    log.warning("%s: skipping malformed CDS: %s", file_path.name, exc)
                    n_skipped_pseudogenes += 1
                    continue
                sequence = str(nt).upper()
                if not sequence:
                    n_skipped_pseudogenes += 1
                    continue

            n_cds += 1

            locus_tag = (quals.get("locus_tag") or [None])[0]
            protein_id = (quals.get("protein_id") or [None])[0]
            gene_name = (quals.get("gene") or [None])[0]
            product = (quals.get("product") or [None])[0]
            ec_number = (quals.get("EC_number") or [None])[0]

            # Defensive handling of duplicate locus_tags in a single file
            # (broken Prokka reruns) — disambiguate with an ordinal suffix.
            if locus_tag:
                seen = locus_tag_occurrences.get(locus_tag, 0)
                locus_tag_occurrences[locus_tag] = seen + 1
                if seen > 0:
                    log.warning(
                        "%s: duplicate locus_tag %r — suffixing __dup%d",
                        file_path.name, locus_tag, seen,
                    )
                    locus_tag = f"{locus_tag}__dup{seen}"

            # gene_ordinal supports the gene-name fallback rung
            gene_key_for_ordinal = gene_name or ""
            gene_ordinals[gene_key_for_ordinal] = (
                gene_ordinals.get(gene_key_for_ordinal, 0) + 1
            )
            gene_ordinal = gene_ordinals[gene_key_for_ordinal]

            try:
                start = int(feat.location.start)
                end = int(feat.location.end)
            except (TypeError, AttributeError):
                start = 0
                end = 0
            strand = int(feat.location.strand or 0) if feat.location is not None else 0

            cds_key, cds_key_source = resolve_cds_key(
                locus_tag=locus_tag,
                protein_id=protein_id,
                gene=gene_name,
                gene_ordinal=gene_ordinal,
                contig=contig,
                start=start,
                end=end,
                strand=strand,
            )

            aa_length = len(sequence) if level == "protein" else max(len(sequence) // 3, 0)

            rows.append(
                {
                    "cds_key": cds_key,
                    "cds_key_source": cds_key_source,
                    "locus_tag": locus_tag,
                    "protein_id": protein_id,
                    "gene": gene_name,
                    "product": product,
                    "ec_number": ec_number,
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "sequence": sequence,
                    "aa_length": aa_length,
                }
            )

    genome_key = resolve_genome_key(
        file_path=file_path,
        dblink=dblink_entries,
        locus_version=first_locus_version,
        organism=organism,
        strain=strain,
        manifest_name=manifest_name,
    )

    return ParsedGenome(
        file_path=file_path,
        file_sha256=file_sha256,
        file_bytes=file_bytes,
        genome_key=genome_key,
        organism=organism,
        strain=strain,
        definition=definition,
        n_records=n_records,
        n_cds=n_cds,
        n_skipped_pseudogenes=n_skipped_pseudogenes,
        total_length=total_length,
        gc_count=gc_count,
        rows=rows,
    )


# ---------------------------------------------------------------------------
# Directory-level parallel driver
# ---------------------------------------------------------------------------


def parse_folder(
    input_dir: Path,
    *,
    level: str = "protein",
    threads: int = 0,
    manifest: Optional[dict[str, str]] = None,
) -> list[ParsedGenome]:
    """Parse all GenBank files in ``input_dir`` in parallel.

    Returns results sorted by filename for deterministic `genome_uid`
    assignment in `build_tables`.
    """
    files = find_genbank_files(input_dir)
    if not files:
        raise FileNotFoundError(
            f"No GenBank files (*.gb, *.gbk, *.gbff) in {input_dir}"
        )

    n_workers = threads if threads > 0 else min(len(files), 8)

    results: list[ParsedGenome] = []
    if n_workers <= 1 or len(files) == 1:
        for path in files:
            mname = manifest.get(path.name) if manifest else None
            results.append(parse_genbank_file(path, level=level, manifest_name=mname))
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            future_to_path = {
                pool.submit(
                    parse_genbank_file,
                    path,
                    level=level,
                    manifest_name=(manifest.get(path.name) if manifest else None),
                ): path
                for path in files
            }
            for future in as_completed(future_to_path):
                results.append(future.result())

    results.sort(key=lambda p: p.file_path.name)
    return results


# ---------------------------------------------------------------------------
# Aggregation — assign uids and emit Arrow tables
# ---------------------------------------------------------------------------


def build_tables(
    parsed: list[ParsedGenome],
    level: str,
) -> tuple[pa.Table, pa.Table, pa.Table]:
    """Assign ``genome_uid`` / ``protein_uid`` and produce Arrow tables.

    Returns ``(id_map, gene_table, genome_meta)``. All three are
    schema-validated before return.
    """
    if level not in SUPPORTED_LEVELS:
        raise ValueError(f"unknown --level: {level!r}")

    sequence_column = "translation" if level == "protein" else "nt_sequence"

    # Build columnar dicts directly so we can hand them to pa.table(...).
    id_map_cols: dict[str, list[Any]] = {
        "protein_uid": [], "genome_uid": [], "gene_uid": [],
        "genome_key": [], "genome_key_source": [],
        "cds_key": [], "cds_key_source": [],
        "assembly_prefix": [], "assembly_version": [],
        "organism": [], "strain": [],
        "locus_tag": [], "protein_id": [], "gene": [],
        "product": [], "ec_number": [],
        "contig": [], "start": [], "end": [], "strand": [],
        "aa_length": [],
    }
    gene_cols: dict[str, list[Any]] = {
        "protein_uid": [], "genome_uid": [], "length": [], sequence_column: [],
    }
    meta_cols: dict[str, list[Any]] = {
        "genome_uid": [], "genome_key": [], "file_path": [],
        "file_sha256": [], "file_bytes": [],
        "n_records": [], "n_cds": [], "n_skipped_pseudogenes": [],
        "organism": [], "strain": [], "definition": [],
        "assembly_prefix": [], "assembly_version": [],
        "total_length": [], "gc_percent": [],
    }

    seen_genome_keys: dict[str, Path] = {}
    next_genome_uid = 0

    for parsed_genome in parsed:
        gkey = parsed_genome.genome_key.key
        if gkey in seen_genome_keys:
            raise ValueError(
                f"Duplicate genome_key {gkey!r}: "
                f"{seen_genome_keys[gkey].name} and {parsed_genome.file_path.name} "
                f"resolved to the same key. Use --manifest to disambiguate."
            )
        seen_genome_keys[gkey] = parsed_genome.file_path

        if next_genome_uid >= (1 << 16):
            raise ValueError(
                "Exceeded 65535 input genomes in a single run "
                "(uint16 genome_uid limit)."
            )
        genome_uid = next_genome_uid
        next_genome_uid += 1

        meta_cols["genome_uid"].append(genome_uid)
        meta_cols["genome_key"].append(gkey)
        meta_cols["file_path"].append(str(parsed_genome.file_path))
        meta_cols["file_sha256"].append(parsed_genome.file_sha256)
        meta_cols["file_bytes"].append(parsed_genome.file_bytes)
        meta_cols["n_records"].append(parsed_genome.n_records)
        meta_cols["n_cds"].append(parsed_genome.n_cds)
        meta_cols["n_skipped_pseudogenes"].append(parsed_genome.n_skipped_pseudogenes)
        meta_cols["organism"].append(parsed_genome.organism)
        meta_cols["strain"].append(parsed_genome.strain)
        meta_cols["definition"].append(parsed_genome.definition)
        meta_cols["assembly_prefix"].append(parsed_genome.genome_key.assembly_prefix)
        meta_cols["assembly_version"].append(parsed_genome.genome_key.assembly_version)
        meta_cols["total_length"].append(parsed_genome.total_length)
        if parsed_genome.total_length > 0:
            meta_cols["gc_percent"].append(
                round(parsed_genome.gc_count * 100.0 / parsed_genome.total_length, 2)
            )
        else:
            meta_cols["gc_percent"].append(None)

        for gene_uid, row in enumerate(parsed_genome.rows):
            protein_uid = make_protein_uid(genome_uid, gene_uid)

            id_map_cols["protein_uid"].append(protein_uid)
            id_map_cols["genome_uid"].append(genome_uid)
            id_map_cols["gene_uid"].append(gene_uid)
            id_map_cols["genome_key"].append(gkey)
            id_map_cols["genome_key_source"].append(parsed_genome.genome_key.source)
            id_map_cols["cds_key"].append(row["cds_key"])
            id_map_cols["cds_key_source"].append(row["cds_key_source"])
            id_map_cols["assembly_prefix"].append(parsed_genome.genome_key.assembly_prefix)
            id_map_cols["assembly_version"].append(parsed_genome.genome_key.assembly_version)
            id_map_cols["organism"].append(parsed_genome.organism)
            id_map_cols["strain"].append(parsed_genome.strain)
            id_map_cols["locus_tag"].append(row["locus_tag"])
            id_map_cols["protein_id"].append(row["protein_id"])
            id_map_cols["gene"].append(row["gene"])
            id_map_cols["product"].append(row["product"])
            id_map_cols["ec_number"].append(row["ec_number"])
            id_map_cols["contig"].append(row["contig"])
            id_map_cols["start"].append(row["start"])
            id_map_cols["end"].append(row["end"])
            id_map_cols["strand"].append(row["strand"])
            id_map_cols["aa_length"].append(row["aa_length"])

            gene_cols["protein_uid"].append(protein_uid)
            gene_cols["genome_uid"].append(genome_uid)
            gene_cols["length"].append(len(row["sequence"]))
            gene_cols[sequence_column].append(row["sequence"])

    id_map_table = pa.table(id_map_cols, schema=ID_MAP_SCHEMA)
    gene_table_ = pa.table(gene_cols, schema=gene_table_schema(level))
    genome_meta_table = pa.table(meta_cols, schema=GENOME_META_SCHEMA)

    validate_schema(id_map_table, ID_MAP_SCHEMA, "id_map.parquet")
    validate_schema(gene_table_, gene_table_schema(level), "gene_table.parquet")
    validate_schema(genome_meta_table, GENOME_META_SCHEMA, "genome_meta.parquet")

    return id_map_table, gene_table_, genome_meta_table
