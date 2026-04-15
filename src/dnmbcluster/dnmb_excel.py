"""BPGA-compatible comparative genomics Excel exporter.

BPGA's canonical output is a wide multi-factor table where each input
genome gets one column per "factor" (locus_tag, product, identity) and
each row is one cluster (Gene_number). DNMBcluster mirrors that format
directly so BPGA users can drop their existing downstream analyses
onto the new pipeline.

We also produce a simple 3-column long format table because it is the
most interoperable shape for custom scripts: one row per (cluster,
genome, locus_tag) triple, trivially filterable in Excel / pandas /
dplyr.

Both formats land in one multi-sheet xlsx file:

  comparative_genomics_N.xlsx
    ├── long              (cluster_id, genome_key, locus_tag, product, pct_identity)
    ├── locus_tag         (wide pivot: cluster_id x genome_key -> locus_tag)
    ├── product           (wide pivot: cluster_id x genome_key -> product)
    ├── identity          (wide pivot: cluster_id x genome_key -> pct_identity_fwd)
    └── merged            (BPGA casting_merge: per-genome 3 columns side-by-side)

Missing cells are filled with "Not assigned" to match BPGA convention.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from .schemas import (
    CLUSTER_SUMMARY_SCHEMA,
    CLUSTERS_SCHEMA,
    GENOME_META_SCHEMA,
    ID_MAP_SCHEMA,
    validate_schema,
)

NOT_ASSIGNED = "Not assigned"


def _load_joined(
    clusters_path: Path,
    id_map_path: Path,
    genome_meta_path: Path,
    cluster_summary_path: Path,
) -> tuple[pd.DataFrame, list[str], pd.DataFrame]:
    """Read Parquet inputs, validate, and return the per-CDS long table.

    Returns ``(joined_df, ordered_genome_keys, cluster_summary_df)``.
    """
    clusters = pq.read_table(clusters_path)
    id_map = pq.read_table(id_map_path)
    genome_meta = pq.read_table(genome_meta_path)
    summary = pq.read_table(cluster_summary_path)

    validate_schema(clusters, CLUSTERS_SCHEMA, "clusters.parquet")
    validate_schema(id_map, ID_MAP_SCHEMA, "id_map.parquet")
    validate_schema(genome_meta, GENOME_META_SCHEMA, "genome_meta.parquet")
    validate_schema(summary, CLUSTER_SUMMARY_SCHEMA, "cluster_summary.parquet")

    clusters_df = clusters.select(
        ["protein_uid", "cluster_id", "is_centroid", "pct_identity_fwd"]
    ).to_pandas()
    id_map_df = id_map.select(
        ["protein_uid", "genome_key", "locus_tag", "cds_key", "product", "gene", "protein_id"]
    ).to_pandas()

    joined = clusters_df.merge(id_map_df, on="protein_uid", how="left")
    joined["locus_or_fallback"] = joined["locus_tag"].where(
        joined["locus_tag"].notna() & (joined["locus_tag"] != ""),
        joined["cds_key"],
    )

    meta_df = genome_meta.to_pandas().sort_values("genome_uid")
    ordered_genome_keys = meta_df["genome_key"].tolist()

    return joined, ordered_genome_keys, summary.to_pandas()


def _pivot(
    joined: pd.DataFrame,
    genome_keys: list[str],
    value_col: str,
    *,
    aggfunc: str = "first",
    fmt: str | None = None,
) -> pd.DataFrame:
    """Pivot ``joined`` into cluster_id x genome_key on ``value_col``.

    Multiple CDSs from the same genome in the same cluster are joined
    with ``"; "`` separation (equivalent to BPGA's ``paste(collapse=", ")``).
    """
    df = joined.copy()
    if fmt == "round2":
        df[value_col] = df[value_col].round(2)
    df[value_col] = df[value_col].astype(str).replace({"nan": pd.NA, "<NA>": pd.NA, "": pd.NA})

    if aggfunc == "concat":
        def joiner(series: pd.Series) -> str:
            vals = [v for v in series.dropna().astype(str).tolist() if v]
            return "; ".join(vals) if vals else NOT_ASSIGNED
        agg = joiner
    else:
        def first_non_null(series: pd.Series) -> str:
            for v in series:
                if pd.notna(v) and v != "":
                    return str(v)
            return NOT_ASSIGNED
        agg = first_non_null

    pivot = df.pivot_table(
        index="cluster_id",
        columns="genome_key",
        values=value_col,
        aggfunc=agg,
        observed=True,
    )
    # Ensure every genome_key column exists in the expected order
    for g in genome_keys:
        if g not in pivot.columns:
            pivot[g] = NOT_ASSIGNED
    pivot = pivot[genome_keys]
    pivot = pivot.fillna(NOT_ASSIGNED)
    pivot.index.name = "cluster_id"
    return pivot.reset_index()


def _build_merged(
    locus_pivot: pd.DataFrame,
    product_pivot: pd.DataFrame,
    identity_pivot: pd.DataFrame,
    genome_keys: list[str],
) -> pd.DataFrame:
    """Interleave the three factor pivots per genome (BPGA casting_merge).

    Column layout:
        cluster_id
        genome1_locus_tag  genome1_product  genome1_identity
        genome2_locus_tag  genome2_product  genome2_identity
        ...
    """
    frames = []
    frames.append(locus_pivot[["cluster_id"]])
    for g in genome_keys:
        block = pd.DataFrame(
            {
                f"{g}__locus_tag":   locus_pivot[g].values,
                f"{g}__product":     product_pivot[g].values,
                f"{g}__pct_identity": identity_pivot[g].values,
            }
        )
        frames.append(block)
    return pd.concat(frames, axis=1)


def _build_long(joined: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Long-format tables: a strict 3-column view and a detailed view.

    - ``long`` — exactly ``cluster_id, genome_key, locus_tag``. One row
      per (cluster, CDS) pair. The minimal shape for filtering in Excel
      or downstream scripts.
    - ``long_detailed`` — same rows, plus ``product``, ``pct_identity``,
      ``is_centroid``, ``gene`` for richer inspection.
    """
    detailed = joined[
        [
            "cluster_id", "genome_key", "locus_or_fallback",
            "gene", "product", "pct_identity_fwd", "is_centroid",
        ]
    ].copy()
    detailed.columns = [
        "cluster_id", "genome_key", "locus_tag",
        "gene", "product", "pct_identity", "is_centroid",
    ]
    detailed["gene"] = detailed["gene"].fillna("")
    detailed["product"] = detailed["product"].fillna(NOT_ASSIGNED)
    detailed["pct_identity"] = detailed["pct_identity"].round(2)
    detailed = detailed.sort_values(
        ["cluster_id", "genome_key", "locus_tag"]
    ).reset_index(drop=True)

    minimal = detailed[["cluster_id", "genome_key", "locus_tag"]].copy()

    return minimal, detailed


def write_comparative_genomics_xlsx(
    clusters_path: Path,
    id_map_path: Path,
    genome_meta_path: Path,
    cluster_summary_path: Path,
    *,
    out_path: Path,
) -> tuple[int, int]:
    """Write the multi-sheet BPGA-compatible Excel file.

    Returns ``(n_clusters, n_genomes)``.
    """
    joined, genome_keys, summary_df = _load_joined(
        clusters_path, id_map_path, genome_meta_path, cluster_summary_path,
    )

    long_df, long_detailed_df = _build_long(joined)

    locus_pivot = _pivot(joined, genome_keys, "locus_or_fallback", aggfunc="concat")
    product_pivot = _pivot(joined, genome_keys, "product", aggfunc="first")
    identity_pivot = _pivot(joined, genome_keys, "pct_identity_fwd", aggfunc="first", fmt="round2")

    # Attach summary columns (gene name, annotation, category) to each pivot
    summary_slim = summary_df[
        ["cluster_id", "representative_gene", "representative_product", "category", "n_genomes"]
    ].rename(
        columns={
            "representative_gene": "gene",
            "representative_product": "annotation",
        }
    )
    summary_slim["gene"] = summary_slim["gene"].fillna("").replace("", pd.NA)
    summary_slim["gene"] = summary_slim["gene"].where(
        summary_slim["gene"].notna(),
        summary_slim["cluster_id"].apply(lambda c: f"group_{c}"),
    )
    summary_slim["annotation"] = summary_slim["annotation"].fillna(NOT_ASSIGNED)

    def with_summary(pivot: pd.DataFrame) -> pd.DataFrame:
        return summary_slim.merge(pivot, on="cluster_id", how="right")

    locus_sheet = with_summary(locus_pivot)
    product_sheet = with_summary(product_pivot)
    identity_sheet = with_summary(identity_pivot)
    merged_sheet = with_summary(
        _build_merged(locus_pivot, product_pivot, identity_pivot, genome_keys)
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        long_df.to_excel(writer, sheet_name="long", index=False)
        long_detailed_df.to_excel(writer, sheet_name="long_detailed", index=False)
        locus_sheet.to_excel(writer, sheet_name="locus_tag", index=False)
        product_sheet.to_excel(writer, sheet_name="product", index=False)
        identity_sheet.to_excel(writer, sheet_name="identity", index=False)
        merged_sheet.to_excel(writer, sheet_name="merged", index=False)

        # Freeze panes + auto column widths
        for sheet_name, df in [
            ("long", long_df),
            ("long_detailed", long_detailed_df),
            ("locus_tag", locus_sheet),
            ("product", product_sheet),
            ("identity", identity_sheet),
            ("merged", merged_sheet),
        ]:
            ws = writer.sheets[sheet_name]
            ws.freeze_panes(1, 1)
            for col_idx, col in enumerate(df.columns):
                if len(df) == 0:
                    max_len = len(str(col))
                else:
                    series = df[col].astype(str).fillna("")
                    values_len = series.str.len().max()
                    if pd.isna(values_len):
                        values_len = 0
                    max_len = max(int(values_len), len(str(col)))
                ws.set_column(col_idx, col_idx, min(max_len + 2, 40))

    return len(locus_pivot), len(genome_keys)
