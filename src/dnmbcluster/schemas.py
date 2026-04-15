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

# ---------------------------------------------------------------------------
# presence_absence.parquet — one row per cluster, bitmap of genome membership
# ---------------------------------------------------------------------------
#
# genome_bitmap is a list of uint64 words of length ceil(n_genomes / 64);
# bit i of word (i // 64) represents genome_uid i. Single-word lists are
# used for n_genomes <= 64 (the common case).
#
# category follows the standard pan-genome definition:
#   core      — present in all n_total genomes
#   soft_core — present in >= 95% * n_total
#   shell     — present in [2, 0.95 * n_total)
#   cloud     — present in exactly 1 genome

PRESENCE_ABSENCE_SCHEMA: Final[pa.Schema] = pa.schema([
    pa.field("cluster_id", pa.uint32(), nullable=False),
    pa.field("n_genomes", pa.uint16(), nullable=False),
    pa.field("n_sequences", pa.uint32(), nullable=False),
    pa.field("category", pa.string(), nullable=False),
    pa.field("genome_bitmap", pa.list_(pa.uint64()), nullable=False),
])

# ---------------------------------------------------------------------------
# pan_core_curve.parquet — cumulative curves across random permutations
# ---------------------------------------------------------------------------

PAN_CORE_CURVE_SCHEMA: Final[pa.Schema] = pa.schema([
    pa.field("permutation", pa.uint16(), nullable=False),
    pa.field("k", pa.uint16(), nullable=False),
    pa.field("pan", pa.uint32(), nullable=False),
    pa.field("core", pa.uint32(), nullable=False),
])

# ---------------------------------------------------------------------------
# cluster_summary.parquet — per-cluster metadata + category + representative
# ---------------------------------------------------------------------------

CLUSTER_SUMMARY_SCHEMA: Final[pa.Schema] = pa.schema([
    pa.field("cluster_id", pa.uint32(), nullable=False),
    pa.field("n_genomes", pa.uint16(), nullable=False),
    pa.field("n_sequences", pa.uint32(), nullable=False),
    pa.field("category", pa.string(), nullable=False),
    pa.field("representative_uid", pa.uint64(), nullable=False),
    pa.field("representative_gene", pa.string(), nullable=True),
    pa.field("representative_product", pa.string(), nullable=True),
    pa.field("representative_locus_tag", pa.string(), nullable=True),
    pa.field("representative_protein_id", pa.string(), nullable=True),
])

CATEGORY_CORE: Final[str] = "core"
CATEGORY_SOFT_CORE: Final[str] = "soft_core"
CATEGORY_SHELL: Final[str] = "shell"
CATEGORY_CLOUD: Final[str] = "cloud"
CATEGORIES: Final[tuple[str, ...]] = (
    CATEGORY_CORE, CATEGORY_SOFT_CORE, CATEGORY_SHELL, CATEGORY_CLOUD,
)


def categorize(n_genomes: int, n_total: int, soft_core_frac: float = 0.95) -> str:
    if n_genomes >= n_total:
        return CATEGORY_CORE
    if n_genomes >= int(round(soft_core_frac * n_total)):
        return CATEGORY_SOFT_CORE
    if n_genomes == 1:
        return CATEGORY_CLOUD
    return CATEGORY_SHELL

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
