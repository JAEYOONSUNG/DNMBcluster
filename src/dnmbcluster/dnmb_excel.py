"""BPGA-compatible comparative genomics Excel exporter.

The Excel artifact contains a single ``merged`` sheet — the BPGA
``casting_merge`` view where each genome contributes a side-by-side
block of columns and each row is one cluster. Block layout is
user-configurable via ``block_columns``:

    [locus_tag]  +  [caller-chosen middle columns]  +  [identity (%)]

``locus_tag`` is always the first column of each genome block (the
identifier anchor) and the forward identity is always the last (the
primary metric). Everything between the two is picked by the CLI's
``--columns`` flag — see ``parse_columns_option`` below for the
whitelist. Default is a single ``product`` column, matching the
original BPGA triple.

Cells where a cluster has no member from a given genome are written as
**blank**, not the BPGA legacy ``"Not assigned"`` sentinel — blank is
unambiguous in Excel formulas and avoids polluting numeric columns.
Numeric columns (identities, coverages, aa_length, alignment_length)
are written as real floats/ints so Excel handles them as numbers.
"""
from __future__ import annotations

from pathlib import Path
from typing import Final

import pandas as pd
import pyarrow.parquet as pq

from .schemas import (
    CLUSTER_SUMMARY_SCHEMA,
    CLUSTERS_SCHEMA,
    GENOME_META_SCHEMA,
    ID_MAP_SCHEMA,
    validate_schema,
)


# ---------------------------------------------------------------------------
# Column whitelist + origin routing
# ---------------------------------------------------------------------------

#: Middle columns users can ask for via --columns. Mapped to the
#: parquet file they originate in — everything here is loaded from
#: disk on demand, avoiding unnecessary wide reads when the caller
#: only wants the default ``product`` column.
_MIDDLE_COLUMNS_FROM_ID_MAP: Final[frozenset[str]] = frozenset({
    "product", "gene", "protein_id", "ec_number", "contig", "aa_length",
})
_MIDDLE_COLUMNS_FROM_CLUSTERS: Final[frozenset[str]] = frozenset({
    "pct_identity_rev", "member_coverage", "rep_coverage", "alignment_length",
})
_ALLOWED_MIDDLE_COLUMNS: Final[frozenset[str]] = (
    _MIDDLE_COLUMNS_FROM_ID_MAP | _MIDDLE_COLUMNS_FROM_CLUSTERS
)

#: Fixed columns — cannot be duplicated by --columns.
_FIXED_FIRST: Final[str] = "locus_tag"   # identifier anchor (with cds_key fallback)
_FIXED_LAST: Final[str]  = "pct_identity_fwd"

#: Display-time label overrides. Parquet column names on disk are
#: kept literal; these rewrites apply only to the bottom header row.
_LABEL_REWRITE: Final[dict[str, str]] = {
    "pct_identity_fwd": "identity (%)",
    "pct_identity_rev": "identity_rev (%)",
    "member_coverage":  "member_cov",
    "rep_coverage":     "rep_cov",
}

#: Columns that should be rendered with a 2-decimal numeric format.
_TWO_DECIMAL_COLUMNS: Final[frozenset[str]] = frozenset({
    "pct_identity_fwd", "pct_identity_rev", "member_coverage", "rep_coverage",
})
#: Columns that are numeric integers (no decimal format).
_INTEGER_COLUMNS: Final[frozenset[str]] = frozenset({
    "aa_length", "alignment_length",
})


def parse_columns_option(raw: str) -> list[str]:
    """Parse the comma-separated ``--columns`` CLI string.

    Returns the cleaned, de-duplicated list of middle column names in
    user-supplied order. Raises ``ValueError`` with an actionable
    message on any of:

    - unknown column name (not in ``_ALLOWED_MIDDLE_COLUMNS``)
    - duplicate entry
    - inclusion of the fixed ``locus_tag`` or ``pct_identity_fwd`` columns
    """
    items = [c.strip() for c in (raw or "").split(",") if c.strip()]
    seen: set[str] = set()
    cleaned: list[str] = []
    for col in items:
        if col in (_FIXED_FIRST, _FIXED_LAST):
            raise ValueError(
                f"{col!r} is a fixed block column and cannot be listed "
                f"in --columns (locus_tag is always first, identity (%) "
                f"is always last)."
            )
        if col == "locus_or_fallback":
            # Internal alias; reject so users can't stumble into it.
            raise ValueError(
                f"{col!r} is an internal alias; use {_FIXED_FIRST!r} instead."
            )
        if col not in _ALLOWED_MIDDLE_COLUMNS:
            allowed = ", ".join(sorted(_ALLOWED_MIDDLE_COLUMNS))
            raise ValueError(
                f"unknown column {col!r}. Allowed: {allowed}"
            )
        if col in seen:
            raise ValueError(f"duplicate column {col!r} in --columns.")
        seen.add(col)
        cleaned.append(col)
    return cleaned


# ---------------------------------------------------------------------------
# Data loading + join
# ---------------------------------------------------------------------------


def _load_joined(
    clusters_path: Path,
    id_map_path: Path,
    genome_meta_path: Path,
    cluster_summary_path: Path,
    *,
    extra_id_map_cols: list[str],
    extra_cluster_cols: list[str],
) -> tuple[pd.DataFrame, list[str], dict[str, str], pd.DataFrame]:
    """Read Parquet inputs, validate, and join clusters + id_map + meta.

    Returns ``(joined_df, ordered_genome_keys, genome_label_map, summary_df)``.
    Only the columns requested via ``extra_*`` are loaded beyond the
    mandatory base set so wide runs (lots of selected middle columns)
    stay roughly pay-as-you-go in memory.
    """
    clusters_full = pq.read_table(clusters_path)
    id_map_full = pq.read_table(id_map_path)
    genome_meta = pq.read_table(genome_meta_path)
    summary = pq.read_table(cluster_summary_path)

    validate_schema(clusters_full, CLUSTERS_SCHEMA, "clusters.parquet")
    validate_schema(id_map_full, ID_MAP_SCHEMA, "id_map.parquet")
    validate_schema(genome_meta, GENOME_META_SCHEMA, "genome_meta.parquet")
    validate_schema(summary, CLUSTER_SUMMARY_SCHEMA, "cluster_summary.parquet")

    base_cluster_cols = [
        "protein_uid", "cluster_id", "is_centroid", "pct_identity_fwd",
    ]
    cluster_cols = base_cluster_cols + [
        c for c in extra_cluster_cols if c not in base_cluster_cols
    ]
    clusters_df = clusters_full.select(cluster_cols).to_pandas()

    base_id_map_cols = [
        "protein_uid", "genome_key", "locus_tag", "cds_key",
    ]
    id_map_cols = base_id_map_cols + [
        c for c in extra_id_map_cols if c not in base_id_map_cols
    ]
    id_map_df = id_map_full.select(id_map_cols).to_pandas()

    joined = clusters_df.merge(id_map_df, on="protein_uid", how="left")
    joined["locus_or_fallback"] = joined["locus_tag"].where(
        joined["locus_tag"].notna() & (joined["locus_tag"] != ""),
        joined["cds_key"],
    )

    meta_df = genome_meta.to_pandas().sort_values("genome_uid")
    ordered_genome_keys = meta_df["genome_key"].tolist()

    label_map: dict[str, str] = {}
    for _, row in meta_df.iterrows():
        gkey = row["genome_key"]
        definition = row.get("definition") if "definition" in row.index else None
        organism = row.get("organism")
        strain = row.get("strain")
        if isinstance(definition, str) and definition.strip():
            label_map[gkey] = definition.strip()
        elif isinstance(organism, str) and organism.strip():
            if (
                isinstance(strain, str) and strain.strip()
                and strain.strip() not in organism
            ):
                label_map[gkey] = f"{organism.strip()} {strain.strip()}"
            else:
                label_map[gkey] = organism.strip()
        else:
            label_map[gkey] = gkey

    return joined, ordered_genome_keys, label_map, summary.to_pandas()


# ---------------------------------------------------------------------------
# Pivot + merge helpers
# ---------------------------------------------------------------------------


def _pivot(
    joined: pd.DataFrame,
    genome_keys: list[str],
    value_col: str,
    *,
    aggfunc: str = "first",
    round_decimals: int | None = None,
) -> pd.DataFrame:
    """Pivot ``joined`` into ``cluster_id x genome_key`` on ``value_col``."""
    df = joined[["cluster_id", "genome_key", value_col]]

    if aggfunc == "concat":
        def agg_fn(series: pd.Series):
            vals = [str(v) for v in series.dropna().tolist() if str(v)]
            return "; ".join(vals) if vals else pd.NA
    elif aggfunc == "first":
        def agg_fn(series: pd.Series):
            non_null = series.dropna()
            return non_null.iloc[0] if len(non_null) else pd.NA
    else:
        raise ValueError(f"unknown aggfunc {aggfunc!r}")

    pivot = df.pivot_table(
        index="cluster_id",
        columns="genome_key",
        values=value_col,
        aggfunc=agg_fn,
        observed=True,
    )
    for g in genome_keys:
        if g not in pivot.columns:
            pivot[g] = pd.NA
    pivot = pivot[genome_keys]

    if round_decimals is not None:
        pivot = pivot.round(round_decimals)

    pivot.index.name = "cluster_id"
    return pivot.reset_index()


def _aggfunc_for(parquet_col: str) -> str:
    """locus_tag gets a ``"; "``-join across multiple CDSs per genome;
    every other column takes the first non-null."""
    if parquet_col == "locus_or_fallback":
        return "concat"
    return "first"


def _round_for(parquet_col: str) -> int | None:
    return 2 if parquet_col in _TWO_DECIMAL_COLUMNS else None


def _build_merged(
    pivots: dict[str, pd.DataFrame],
    genome_keys: list[str],
    block_parquet_order: list[str],
) -> pd.DataFrame:
    """Interleave N pivots (one per block column) per genome.

    ``block_parquet_order`` holds the parquet column names in the
    order they should appear in each genome block — e.g.
    ``["locus_or_fallback", "product", "gene", "pct_identity_fwd"]``.

    ``pivot_table`` with ``observed=True`` drops cluster_ids that have
    no non-null value in the pivoted column, so different pivots can
    end up with different row counts. Reindex every pivot to the
    union of cluster_ids before interleaving so the resulting frame
    is well-formed regardless of which middle columns were picked.
    """
    all_cids: pd.Index = pd.Index([], dtype="uint32")
    for col in block_parquet_order:
        all_cids = all_cids.union(pivots[col]["cluster_id"], sort=True)

    canonical = pd.DataFrame({"cluster_id": all_cids.astype("uint32")})

    aligned: dict[str, pd.DataFrame] = {}
    for col in block_parquet_order:
        reindexed = (
            pivots[col]
            .set_index("cluster_id")
            .reindex(canonical["cluster_id"])
            .reset_index()
        )
        aligned[col] = reindexed

    frames: list[pd.DataFrame] = [canonical]
    for g in genome_keys:
        block = pd.DataFrame(
            {
                f"{g}__{col}": aligned[col][g].values
                for col in block_parquet_order
            }
        )
        frames.append(block)
    return pd.concat(frames, axis=1)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def write_comparative_genomics_xlsx(
    clusters_path: Path,
    id_map_path: Path,
    genome_meta_path: Path,
    cluster_summary_path: Path,
    *,
    out_path: Path,
    block_columns: list[str] | None = None,
) -> tuple[int, int]:
    """Write the single-sheet BPGA ``merged`` Excel artifact.

    ``block_columns`` is the user-supplied middle column list (from
    the ``--columns`` CLI flag, already validated by
    ``parse_columns_option``). The on-disk per-genome block is
    ``[locus_tag]  +  block_columns  +  [pct_identity_fwd]``.

    Returns ``(n_clusters, n_genomes)``.
    """
    if block_columns is None:
        block_columns = ["product"]

    # Decide which raw parquet files each requested column comes from.
    extra_id_map_cols = [
        c for c in block_columns if c in _MIDDLE_COLUMNS_FROM_ID_MAP
    ]
    extra_cluster_cols = [
        c for c in block_columns if c in _MIDDLE_COLUMNS_FROM_CLUSTERS
    ]

    joined, genome_keys, genome_labels, summary_df = _load_joined(
        clusters_path, id_map_path, genome_meta_path, cluster_summary_path,
        extra_id_map_cols=extra_id_map_cols,
        extra_cluster_cols=extra_cluster_cols,
    )

    # The internal parquet-column order for each genome block.
    # locus_tag comes from the joined ``locus_or_fallback`` helper
    # column so missing locus_tags fall back to cds_key automatically.
    block_parquet_order: list[str] = ["locus_or_fallback"]
    block_parquet_order.extend(block_columns)
    block_parquet_order.append(_FIXED_LAST)

    # Build one pivot per block column, keyed by the parquet name.
    pivots: dict[str, pd.DataFrame] = {}
    for col in block_parquet_order:
        pivots[col] = _pivot(
            joined, genome_keys, col,
            aggfunc=_aggfunc_for(col),
            round_decimals=_round_for(col),
        )

    merged = _build_merged(pivots, genome_keys, block_parquet_order)

    # Attach summary columns (gene/annotation/category/n_genomes).
    summary_slim = summary_df[[
        "cluster_id",
        "representative_gene",
        "representative_product",
        "category",
        "n_genomes",
    ]].rename(
        columns={
            "representative_gene": "gene",
            "representative_product": "annotation",
        }
    )
    summary_slim["gene"] = summary_slim["gene"].replace("", pd.NA)

    merged_sheet = summary_slim.merge(merged, on="cluster_id", how="right")

    # Column layout metadata used for borders / labels / widths.
    n_meta = 5  # cluster_id, gene, annotation, category, n_genomes
    block_width = len(block_parquet_order)
    n_genomes = len(genome_keys)
    border_col_idxs = {n_meta + block_width * i for i in range(n_genomes)}

    # Columns that need numeric formatting in the xlsx.
    two_dec_sheet_cols: set[str] = set()
    int_sheet_cols: set[str] = set()
    for g in genome_keys:
        for col in block_parquet_order:
            sheet_name = f"{g}__{col}"
            if col in _TWO_DECIMAL_COLUMNS:
                two_dec_sheet_cols.add(sheet_name)
            elif col in _INTEGER_COLUMNS:
                int_sheet_cols.add(sheet_name)

    header_row_top = 0
    header_row_bot = 1
    data_start = 2
    n_rows = len(merged_sheet)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        merged_sheet.to_excel(
            writer, sheet_name="merged", index=False, na_rep="",
            startrow=data_start, header=False,
        )

        ws = writer.sheets["merged"]
        wb = writer.book

        # --- formats -------------------------------------------------
        cell_center = {"align": "center"}
        cell_fmt           = wb.add_format(cell_center)
        cell_num_fmt       = wb.add_format({**cell_center, "num_format": "0.00"})
        cell_int_fmt       = wb.add_format({**cell_center, "num_format": "0"})
        cell_border_fmt    = wb.add_format({**cell_center, "left": 1})
        cell_border_num_fmt = wb.add_format(
            {**cell_center, "num_format": "0.00", "left": 1}
        )
        cell_border_int_fmt = wb.add_format(
            {**cell_center, "num_format": "0", "left": 1}
        )

        header_base = {
            "bold": True,
            "align": "center",
            "valign": "vcenter",
            "bg_color": "#F2F2F2",
            "top": 1,
            "bottom": 1,
        }
        header_fmt          = wb.add_format(header_base)
        header_border_fmt   = wb.add_format({**header_base, "left": 1})
        header_group_fmt    = wb.add_format(
            {**header_base, "left": 1, "text_wrap": True}
        )

        # --- row 0: definition + accession per genome block ----------
        for i, gkey in enumerate(genome_keys):
            first_col = n_meta + block_width * i
            last_col  = first_col + block_width - 1
            definition = genome_labels.get(gkey, gkey)
            label = f"{definition}\n({gkey})"
            ws.merge_range(
                header_row_top, first_col, header_row_top, last_col,
                label, header_group_fmt,
            )
        for col_idx in range(n_meta):
            ws.write_blank(header_row_top, col_idx, None, header_fmt)

        # --- row 1: per-column labels --------------------------------
        for col_idx, col in enumerate(merged_sheet.columns):
            if "__" in col and col_idx >= n_meta:
                raw = col.split("__", 1)[1]
                # Display locus_or_fallback as plain "locus_tag"
                if raw == "locus_or_fallback":
                    raw = "locus_tag"
                label = _LABEL_REWRITE.get(raw, raw)
            else:
                label = col
            fmt = header_border_fmt if col_idx in border_col_idxs else header_fmt
            ws.write(header_row_bot, col_idx, label, fmt)

        ws.freeze_panes(data_start, 1)
        ws.set_row(header_row_top, 36)
        ws.set_row(header_row_bot, 18)

        # --- column widths + per-column formats ----------------------
        def _content_width(column_name: str, label_text: str) -> int:
            non_null = merged_sheet[column_name].dropna()
            values_len = (
                non_null.astype(str).str.len().max() if len(non_null) else 0
            )
            max_len = max(int(values_len), len(label_text))
            return min(max_len + 2, 40)

        # Per-block width is keyed off the locus_tag column (the first
        # column of each genome block) so every column in a block has
        # the same width and the block reads as a tidy grid.
        block_widths: dict[int, int] = {}
        for i in range(n_genomes):
            locus_col_name = merged_sheet.columns[n_meta + block_width * i]
            block_widths[i] = _content_width(locus_col_name, "locus_tag")

        for col_idx, col in enumerate(merged_sheet.columns):
            is_genome_col = col_idx >= n_meta
            has_border = col_idx in border_col_idxs
            is_two_dec = col in two_dec_sheet_cols
            is_int     = col in int_sheet_cols

            if is_genome_col:
                block_idx = (col_idx - n_meta) // block_width
                width = block_widths[block_idx]
                if has_border and is_two_dec:
                    fmt = cell_border_num_fmt
                elif has_border and is_int:
                    fmt = cell_border_int_fmt
                elif has_border:
                    fmt = cell_border_fmt
                elif is_two_dec:
                    fmt = cell_num_fmt
                elif is_int:
                    fmt = cell_int_fmt
                else:
                    fmt = cell_fmt
                ws.set_column(col_idx, col_idx, width, fmt)
            else:
                width = _content_width(col, str(col))
                ws.set_column(col_idx, col_idx, width)

    return n_rows, len(genome_keys)
