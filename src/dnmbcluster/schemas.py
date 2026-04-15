"""PyArrow schemas for DNMBcluster intermediate Parquet files.

Every reader and writer validates against these at stage boundaries.
See SPEED.md Section 1 and 6 for the design rationale.
"""
from __future__ import annotations

from typing import Final

import pyarrow as pa

# ---------------------------------------------------------------------------
# Enum value sets (stored as plain pa.string() — application-level validation)
# ---------------------------------------------------------------------------

GENOME_KEY_SOURCES: Final[tuple[str, ...]] = (
    "manifest",
    "dblink",
    "filename_gcf",
    "locus_version",
    "organism",
    "filename_raw",
)

CDS_KEY_SOURCES: Final[tuple[str, ...]] = (
    "locus_tag",
    "protein_id",
    "gene_ordinal",
    "coord",
)

# ---------------------------------------------------------------------------
# id_map.parquet — the canonical identifier table
# ---------------------------------------------------------------------------

ID_MAP_SCHEMA: Final[pa.Schema] = pa.schema([
    # hot-path integer keys
    pa.field("protein_uid", pa.uint64(), nullable=False),
    pa.field("genome_uid", pa.uint16(), nullable=False),
    pa.field("gene_uid", pa.uint32(), nullable=False),

    # correctness-layer composite key
    pa.field("genome_key", pa.string(), nullable=False),
    pa.field("genome_key_source", pa.string(), nullable=False),
    pa.field("cds_key", pa.string(), nullable=False),
    pa.field("cds_key_source", pa.string(), nullable=False),

    # reannotation detection
    pa.field("assembly_prefix", pa.string(), nullable=True),
    pa.field("assembly_version", pa.string(), nullable=True),

    # original GenBank attributes (non-unique; display + dereplication)
    pa.field("organism", pa.string(), nullable=True),
    pa.field("strain", pa.string(), nullable=True),
    pa.field("locus_tag", pa.string(), nullable=True),
    pa.field("protein_id", pa.string(), nullable=True),
    pa.field("gene", pa.string(), nullable=True),
    pa.field("product", pa.string(), nullable=True),
    pa.field("ec_number", pa.string(), nullable=True),
    pa.field("contig", pa.string(), nullable=False),
    pa.field("start", pa.uint32(), nullable=False),
    pa.field("end", pa.uint32(), nullable=False),
    pa.field("strand", pa.int8(), nullable=False),
    pa.field("aa_length", pa.uint32(), nullable=True),
])

# ---------------------------------------------------------------------------
# gene_table.parquet — sequence payload keyed by protein_uid
# Column set differs by --level (never both, to avoid bloat).
# ---------------------------------------------------------------------------


def gene_table_schema(level: str) -> pa.Schema:
    """Return the gene_table schema for the given `--level`."""
    base = [
        pa.field("protein_uid", pa.uint64(), nullable=False),
        pa.field("genome_uid", pa.uint16(), nullable=False),
        pa.field("length", pa.uint32(), nullable=False),
    ]
    if level == "protein":
        base.append(pa.field("translation", pa.string(), nullable=False))
    elif level == "nucleotide":
        base.append(pa.field("nt_sequence", pa.string(), nullable=False))
    else:
        raise ValueError(f"unknown --level: {level!r}")
    return pa.schema(base)


# ---------------------------------------------------------------------------
# genome_meta.parquet — one row per input GenBank file
# ---------------------------------------------------------------------------

CLUSTERS_SCHEMA: Final[pa.Schema] = pa.schema([
    pa.field("protein_uid", pa.uint64(), nullable=False),
    pa.field("genome_uid", pa.uint16(), nullable=False),
    pa.field("cluster_id", pa.uint32(), nullable=False),
    pa.field("representative_uid", pa.uint64(), nullable=False),
    pa.field("is_centroid", pa.bool_(), nullable=False),
    pa.field("pct_identity", pa.float32(), nullable=True),
])

GENOME_META_SCHEMA: Final[pa.Schema] = pa.schema([
    pa.field("genome_uid", pa.uint16(), nullable=False),
    pa.field("genome_key", pa.string(), nullable=False),
    pa.field("file_path", pa.string(), nullable=False),
    pa.field("file_sha256", pa.string(), nullable=False),
    pa.field("file_bytes", pa.uint64(), nullable=False),
    pa.field("n_records", pa.uint32(), nullable=False),
    pa.field("n_cds", pa.uint32(), nullable=False),
    pa.field("n_skipped_pseudogenes", pa.uint32(), nullable=True),
    pa.field("organism", pa.string(), nullable=True),
    pa.field("strain", pa.string(), nullable=True),
    pa.field("assembly_prefix", pa.string(), nullable=True),
    pa.field("assembly_version", pa.string(), nullable=True),
    pa.field("total_length", pa.uint64(), nullable=True),
    pa.field("gc_percent", pa.float32(), nullable=True),
])


def validate_schema(table: pa.Table, expected: pa.Schema, label: str) -> None:
    """Fail fast if `table`'s schema does not match `expected`.

    Compares field names and types, ignoring metadata and field order.
    """
    if len(table.schema) != len(expected):
        raise ValueError(
            f"{label}: column count mismatch "
            f"(expected {len(expected)}, got {len(table.schema)})"
        )
    for i, field in enumerate(expected):
        actual = table.schema.field(i)
        if actual.name != field.name:
            raise ValueError(
                f"{label}: column {i} name mismatch "
                f"(expected {field.name!r}, got {actual.name!r})"
            )
        if not actual.type.equals(field.type):
            raise ValueError(
                f"{label}: column {field.name!r} type mismatch "
                f"(expected {field.type}, got {actual.type})"
            )
