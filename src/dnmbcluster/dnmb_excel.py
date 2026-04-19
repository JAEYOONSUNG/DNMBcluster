"""BPGA-compatible comparative genomics Excel exporter.

The Excel artifact contains two sheets that share the BPGA
``casting_merge`` layout (each genome gets a side-by-side block,
each row is one cluster). The two sheets correspond to the two
directions of the MMseqs2 reciprocal realign:

* ``forward`` — *member (query) aligned against seed (target)*:
  per-member ``coverage`` (``member_cov``) + ``identity (%)``
  (``identity_member`` = matches / member_length). Answers "for
  each member, how much of it aligns to the seed, and how much of
  its length is identical?"
* ``reverse`` — *seed (query) aligned against member (target)*:
  ``coverage`` (``rep_cov``) + ``identity (%)`` (``identity_seed``
  = matches / seed_length). Answers "when realigned the other way,
  how much of the seed covers each member and how much of the seed
  length is identical?"

Raw MMseqs2 ``pct_identity_fwd`` / ``pct_identity_rev`` are essentially
symmetric (``pident = matches / alnlen``) so they hide length
asymmetry. Multiplying by the direction-appropriate coverage converts
them into sequence-length-normalized identities that differ whenever
member and seed have different lengths:

    identity_member = pident × member_coverage   (sheet 1)
    identity_seed   = pident × rep_coverage      (sheet 2)

Cells where identity == 100% (seed self-alignment / perfect matches)
are rendered **bold** so the anchor rows pop visually.

Block layout per sheet:

    [locus_tag]  +  [user middle cols]  +  [direction coverage]  +  [identity (%)]

``locus_tag`` is always first; the trailing identity is always last.
The direction-appropriate coverage (``member_cov`` for merged,
``rep_cov`` for reverse) is inserted automatically just before the
identity. Other middle columns come from the CLI's ``--columns`` flag
— see ``parse_columns_option`` for the whitelist.

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
    "member_coverage", "rep_coverage", "alignment_length",
})
_ALLOWED_MIDDLE_COLUMNS: Final[frozenset[str]] = (
    _MIDDLE_COLUMNS_FROM_ID_MAP | _MIDDLE_COLUMNS_FROM_CLUSTERS
)

#: Fixed columns — cannot be duplicated by --columns.
_FIXED_FIRST: Final[str] = "locus_tag"   # identifier anchor (with cds_key fallback)
#: Trailing columns (derived at load time — see ``_load_joined``):
#:   forward → identity_member = matches / member_length = pident × member_cov
#:             (member aligned against seed)
#:   reverse → identity_seed   = matches / seed_length   = pident × rep_cov
#:             (seed aligned against member)
_FIXED_LAST: Final[str]     = "identity_member"
_FIXED_LAST_REV: Final[str] = "identity_seed"
#: Direction-appropriate coverage auto-inserted just before the identity.
_COV_COL: Final[str]     = "member_coverage"  # merged sheet
_COV_COL_REV: Final[str] = "rep_coverage"     # reverse sheet

#: Display-time label overrides. Parquet column names on disk are
#: kept literal; these rewrites apply only to the bottom header row.
#: Each sheet is already direction-scoped (merged = forward only,
#: reverse = reverse only), so the trailing columns drop their
#: redundant ``_member`` / ``_seed`` suffix and read as plain
#: ``coverage`` / ``identity (%)``.
_LABEL_REWRITE: Final[dict[str, str]] = {
    "pct_identity_fwd": "pident (%)",
    "pct_identity_rev": "pident_rev (%)",
    "identity_seed":    "identity (%)",
    "identity_member":  "identity (%)",
    "member_coverage":  "coverage",
    "rep_coverage":     "coverage",
}

#: Per-column width overrides (Excel "characters" units). Locus tags,
#: coverages, and identities are short fields — fixing them here
#: prevents the block from being stretched to the locus_tag concat
#: width (which was the default before and left a lot of dead space).
_NARROW_WIDTHS: Final[dict[str, int]] = {
    "locus_or_fallback": 14,   # holds "B9K58_00005"-style tags + a little margin
    "member_coverage":   10,   # header "coverage" is 8 chars
    "rep_coverage":      10,
    "identity_member":   14,   # header "identity (%)" is 12 chars
    "identity_seed":     14,
}

#: Columns that should be rendered with a 2-decimal numeric format.
_TWO_DECIMAL_COLUMNS: Final[frozenset[str]] = frozenset({
    "pct_identity_fwd", "pct_identity_rev",
    "identity_seed", "identity_member",
    "member_coverage", "rep_coverage",
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
        if col in (_FIXED_FIRST, _FIXED_LAST, _FIXED_LAST_REV):
            raise ValueError(
                f"{col!r} is a fixed block column and cannot be listed "
                f"in --columns (locus_tag is always first; identity (%) "
                f"is the trailing column — identity_member on the "
                f"'forward' sheet and identity_seed on the 'reverse' sheet)."
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

    # Both identity directions + coverages are loaded unconditionally
    # — the writer derives length-adjusted identities from them and
    # exposes coverages as default middle columns.
    base_cluster_cols = [
        "protein_uid", "cluster_id", "is_centroid",
        "pct_identity_fwd", "pct_identity_rev",
        "member_coverage", "rep_coverage",
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

    # Derive length-adjusted identities:
    #   identity_seed   = pident × rep_cov    = matches / seed_length
    #   identity_member = pident × member_cov = matches / member_length
    # Raw MMseqs2 pident is alignment-intrinsic (matches/alnlen) and
    # therefore symmetric across query/target swap — multiplying by the
    # respective coverage converts it to a sequence-length-normalized
    # identity that reflects member/seed length differences. Cap at 100
    # since floating-point rounding can push the product to ~100.0001.
    pid = joined["pct_identity_fwd"].astype("float32")
    rcov = joined["rep_coverage"].astype("float32")
    mcov = joined["member_coverage"].astype("float32")
    joined["identity_seed"]   = (pid * rcov).clip(upper=100.0).astype("float32")
    joined["identity_member"] = (pid * mcov).clip(upper=100.0).astype("float32")

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


def _write_one_sheet(
    writer: pd.ExcelWriter,
    sheet_name: str,
    *,
    fixed_last_col: str,
    coverage_col: str,
    joined: pd.DataFrame,
    genome_keys: list[str],
    genome_labels: dict[str, str],
    summary_slim: pd.DataFrame,
    block_columns: list[str],
) -> int:
    """Render one ``merged``-style sheet with direction-specific trailing cols.

    ``forward`` and ``reverse`` differ only in the two direction-
    specific trailing columns (``coverage_col`` just before
    ``fixed_last_col``). Everything else — header layout, block widths,
    formats — is driven off the shared ``joined`` / ``summary_slim``
    frames so each sheet costs only one extra pivot per direction col.
    Returns the data row count.
    """
    # Filter out the "other" direction's coverage if the caller passed
    # it in block_columns — we always auto-insert this sheet's own
    # coverage, never the opposite one.
    other_cov = _COV_COL_REV if coverage_col == _COV_COL else _COV_COL
    user_mids = [c for c in block_columns if c not in (coverage_col, other_cov)]

    block_parquet_order: list[str] = ["locus_or_fallback"]
    block_parquet_order.extend(user_mids)
    block_parquet_order.append(coverage_col)
    block_parquet_order.append(fixed_last_col)

    pivots: dict[str, pd.DataFrame] = {}
    for col in block_parquet_order:
        pivots[col] = _pivot(
            joined, genome_keys, col,
            aggfunc=_aggfunc_for(col),
            round_decimals=_round_for(col),
        )

    merged = _build_merged(pivots, genome_keys, block_parquet_order)
    merged_sheet = summary_slim.merge(merged, on="cluster_id", how="right")

    n_meta = 5  # cluster_id, gene, annotation, category, n_genomes
    block_width = len(block_parquet_order)
    n_genomes = len(genome_keys)
    border_col_idxs = {n_meta + block_width * i for i in range(n_genomes)}

    two_dec_sheet_cols: set[str] = set()
    int_sheet_cols: set[str] = set()
    for g in genome_keys:
        for col in block_parquet_order:
            col_sheet_name = f"{g}__{col}"
            if col in _TWO_DECIMAL_COLUMNS:
                two_dec_sheet_cols.add(col_sheet_name)
            elif col in _INTEGER_COLUMNS:
                int_sheet_cols.add(col_sheet_name)

    header_row_top = 0
    header_row_bot = 1
    data_start = 2

    merged_sheet.to_excel(
        writer, sheet_name=sheet_name, index=False, na_rep="",
        startrow=data_start, header=False,
    )

    ws = writer.sheets[sheet_name]
    wb = writer.book

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

    for col_idx, col in enumerate(merged_sheet.columns):
        if "__" in col and col_idx >= n_meta:
            raw = col.split("__", 1)[1]
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

    def _content_width(column_name: str, label_text: str) -> int:
        non_null = merged_sheet[column_name].dropna()
        values_len = (
            non_null.astype(str).str.len().max() if len(non_null) else 0
        )
        max_len = max(int(values_len), len(label_text))
        return min(max_len + 2, 40)

    def _col_width(parquet_col: str, sheet_col_name: str) -> int:
        """Width for one per-genome cell column.

        Narrow-typed columns (locus_tag, coverage, identity) get a
        compact hard-coded width so the block doesn't inherit the
        widest content — usually a semicolon-joined locus_tag. Other
        middle columns (product, gene, …) still auto-size to content.
        """
        if parquet_col in _NARROW_WIDTHS:
            return _NARROW_WIDTHS[parquet_col]
        label = _LABEL_REWRITE.get(parquet_col, parquet_col)
        return _content_width(sheet_col_name, label)

    for col_idx, col in enumerate(merged_sheet.columns):
        is_genome_col = col_idx >= n_meta
        has_border = col_idx in border_col_idxs
        is_two_dec = col in two_dec_sheet_cols
        is_int     = col in int_sheet_cols

        if is_genome_col:
            parquet_col = col.split("__", 1)[1]
            width = _col_width(parquet_col, col)
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

    # Bold the seed self-alignment cells. The centroid row always
    # scores 100% identity in its own genome column — making those
    # cells bold gives readers an immediate visual anchor for
    # "this is the seed's own value" amid all the realigned members.
    bold_100_fmt = wb.add_format({"bold": True, "num_format": "0.00"})
    n_data_rows = len(merged_sheet)
    if n_data_rows:
        last_row = data_start + n_data_rows - 1
        for col_idx, col in enumerate(merged_sheet.columns):
            if not col.endswith(f"__{fixed_last_col}"):
                continue
            ws.conditional_format(
                data_start, col_idx, last_row, col_idx,
                {"type": "cell", "criteria": "==",
                 "value": 100, "format": bold_100_fmt},
            )

    return n_data_rows


def write_comparative_genomics_xlsx(
    clusters_path: Path,
    id_map_path: Path,
    genome_meta_path: Path,
    cluster_summary_path: Path,
    *,
    out_path: Path,
    block_columns: list[str] | None = None,
) -> tuple[int, int]:
    """Write the two-sheet BPGA Excel artifact.

    ``block_columns`` is the user-supplied middle column list (from
    the ``--columns`` CLI flag, already validated by
    ``parse_columns_option``). The on-disk per-genome block is
    ``[locus_tag]  +  block_columns  +  [coverage]  +  [identity]``
    where ``coverage`` is ``member_cov`` on ``forward`` and ``rep_cov``
    on ``reverse`` (auto-inserted per sheet direction), and
    ``identity`` is the matching length-adjusted value.

    Returns ``(n_clusters, n_genomes)``.
    """
    if block_columns is None:
        # User-facing default: just the product annotation. The
        # direction-appropriate coverage + identity are auto-inserted
        # per sheet so the user gets member_cov + identity_member in
        # forward and rep_cov + identity_seed in reverse without
        # needing to opt in.
        block_columns = ["product"]

    # Decide which raw parquet files each requested column comes from.
    # Both coverage columns are always loaded (needed for the derived
    # length-adjusted identities), so no need to add them here.
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

    # Summary columns attached to both sheets (gene/annotation/category/n_genomes).
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

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        n_rows = _write_one_sheet(
            writer, "forward",
            fixed_last_col=_FIXED_LAST,
            coverage_col=_COV_COL,
            joined=joined,
            genome_keys=genome_keys,
            genome_labels=genome_labels,
            summary_slim=summary_slim,
            block_columns=block_columns,
        )
        _write_one_sheet(
            writer, "reverse",
            fixed_last_col=_FIXED_LAST_REV,
            coverage_col=_COV_COL_REV,
            joined=joined,
            genome_keys=genome_keys,
            genome_labels=genome_labels,
            summary_slim=summary_slim,
            block_columns=block_columns,
        )

    return n_rows, len(genome_keys)
