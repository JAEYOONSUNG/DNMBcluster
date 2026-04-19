"""Engine-native sidecar schemas + writers.

The unified ``clusters.parquet`` deliberately throws away engine-specific
fields so every engine produces identical columns. That's the right call
for downstream code, but it means information unique to each clustering
algorithm — usearch12's CIGAR strings, CD-HIT's per-member alignment
positions — would otherwise be lost between the raw text output on disk
and the canonical Parquet view.

This module preserves those fields in a per-engine **sidecar parquet**
written next to ``clusters.parquet``:

    out_dir/
        clusters.parquet                # canonical, identical schema
        usearch12_native.parquet        # only when --tool usearch12
        cdhit_native.parquet            # only when --tool cd-hit

Both sidecars are keyed by ``protein_uid`` so they join back to the
canonical table on a single column. Centroid rows are intentionally
excluded — there is no member-vs-rep alignment for a sequence aligned
to itself, and the ``is_centroid`` column on the canonical table makes
that distinction available without duplication.

DIAMOND has no sidecar. ``deepclust`` emits only the 2-column
``centroid<TAB>member`` TSV with no per-pair metrics, and the
``--aln-out`` flag was empirically unreliable in DIAMOND 2.1.x — see
the audit notes in commit 41996c8.
"""
from __future__ import annotations

from pathlib import Path
from typing import Final

import pyarrow as pa
import pyarrow.parquet as pq


# ---------------------------------------------------------------------------
# usearch12 native sidecar
# ---------------------------------------------------------------------------
#
# Source columns from the ``.uc`` format (10 tab-separated fields):
#   col 3: percent identity (float, H rows only)
#   col 4: strand (+/-/.)
#   col 7: compressed alignment / CIGAR string (e.g. "342M", "100M2D40M")
#
# Only H rows are emitted; S/C/N rows have no alignment data.

USEARCH12_NATIVE_SCHEMA: Final[pa.Schema] = pa.schema([
    pa.field("protein_uid", pa.uint64(), nullable=False),
    pa.field("representative_uid", pa.uint64(), nullable=False),
    pa.field("native_pct_identity", pa.float32(), nullable=True),
    pa.field("strand", pa.string(), nullable=True),
    pa.field("cigar", pa.string(), nullable=True),
])


# ---------------------------------------------------------------------------
# CD-HIT native sidecar
# ---------------------------------------------------------------------------
#
# Available only when CD-HIT was run with ``-p 1``, which expands the
# ``.clstr`` member line from
#
#     1\t340aa, >uid... at 92.15%
#
# to
#
#     1\t340aa, >uid... at 1:336:5:340/92.15%             (protein)
#     1\t1020nt, >uid... at 1:1020:5:1024/+/92.15%        (nucleotide)
#
# query_start / query_end are the member's aligned region; target_start /
# target_end are the centroid's. All four are 1-based inclusive coordinates
# as written by CD-HIT. ``strand`` is None for protein clustering.

CDHIT_NATIVE_SCHEMA: Final[pa.Schema] = pa.schema([
    pa.field("protein_uid", pa.uint64(), nullable=False),
    pa.field("representative_uid", pa.uint64(), nullable=False),
    pa.field("native_pct_identity", pa.float32(), nullable=True),
    pa.field("query_start", pa.uint32(), nullable=True),
    pa.field("query_end", pa.uint32(), nullable=True),
    pa.field("target_start", pa.uint32(), nullable=True),
    pa.field("target_end", pa.uint32(), nullable=True),
    pa.field("strand", pa.string(), nullable=True),
])


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------


def write_usearch12_native(
    *,
    rows: list[dict],
    out_path: Path,
) -> int:
    """Write the usearch12 native sidecar.

    ``rows`` is a list of dicts with keys matching ``USEARCH12_NATIVE_SCHEMA``.
    Returns the number of rows written.
    """
    table = pa.table(
        {
            "protein_uid":         pa.array([r["protein_uid"] for r in rows], type=pa.uint64()),
            "representative_uid":  pa.array([r["representative_uid"] for r in rows], type=pa.uint64()),
            "native_pct_identity": pa.array([r["native_pct_identity"] for r in rows], type=pa.float32()),
            "strand":              pa.array([r["strand"] for r in rows], type=pa.string()),
            "cigar":               pa.array([r["cigar"] for r in rows], type=pa.string()),
        },
        schema=USEARCH12_NATIVE_SCHEMA,
    )
    from ..io_utils import atomic_write_table
    atomic_write_table(table, out_path)
    return table.num_rows


def write_cdhit_native(
    *,
    rows: list[dict],
    out_path: Path,
) -> int:
    """Write the CD-HIT native sidecar.

    ``rows`` is a list of dicts with keys matching ``CDHIT_NATIVE_SCHEMA``.
    Returns the number of rows written.
    """
    table = pa.table(
        {
            "protein_uid":         pa.array([r["protein_uid"] for r in rows], type=pa.uint64()),
            "representative_uid":  pa.array([r["representative_uid"] for r in rows], type=pa.uint64()),
            "native_pct_identity": pa.array([r["native_pct_identity"] for r in rows], type=pa.float32()),
            "query_start":         pa.array([r["query_start"] for r in rows], type=pa.uint32()),
            "query_end":           pa.array([r["query_end"] for r in rows], type=pa.uint32()),
            "target_start":        pa.array([r["target_start"] for r in rows], type=pa.uint32()),
            "target_end":          pa.array([r["target_end"] for r in rows], type=pa.uint32()),
            "strand":              pa.array([r["strand"] for r in rows], type=pa.string()),
        },
        schema=CDHIT_NATIVE_SCHEMA,
    )
    from ..io_utils import atomic_write_table
    atomic_write_table(table, out_path)
    return table.num_rows
