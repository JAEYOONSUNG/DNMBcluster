"""Long-format per-CDS view of the clustering result.

This module produces ``dnmb/cluster_long.parquet`` — one row per CDS,
with the cluster it joined, the cluster's pan-genome category, the
CDS metadata (gene / product / locus_tag) inlined, the genome label,
and **every alignment metric** from ``clusters.parquet`` carried
through unchanged.

It replaces the ``long_detailed`` sheet that used to live inside
``comparative_genomics_N.xlsx``. Parquet is a much better home for it:

- 33 116 rows (Geobacillus-10) is small for parquet, large for an
  Excel sheet (Excel column-format/freeze-pane gets sluggish around
  20 k rows in many viewers).
- Real ``float32`` identity columns can be filtered numerically in
  pandas / polars / DuckDB / R without the string-coercion gymnastics
  Excel forces.
- Centroid-relative metrics are exactly what ortholog QC scripts want;
  putting them next to the CDS metadata removes a join.

The Excel exporter still produces the wide ``merged`` sheet for users
who specifically want to look at one cluster across all genomes
side-by-side; the long-format slice lives here instead.
"""
from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .schemas import (
    CLUSTER_LONG_SCHEMA,
    CLUSTER_SUMMARY_SCHEMA,
    CLUSTERS_SCHEMA,
    GENOME_META_SCHEMA,
    ID_MAP_SCHEMA,
    validate_schema,
)


def build_cluster_long(
    clusters_path: Path,
    id_map_path: Path,
    genome_meta_path: Path,
    cluster_summary_path: Path,
) -> pa.Table:
    """Join clusters + id_map + genome_meta + summary into a long table.

    The output has exactly one row per row in ``clusters.parquet`` (i.e.
    one row per input CDS), preserving the alignment columns verbatim
    from clusters.parquet so callers see the same numbers downstream
    code already trusts.
    """
    clusters = pq.read_table(clusters_path)
    id_map = pq.read_table(id_map_path)
    genome_meta = pq.read_table(genome_meta_path)
    summary = pq.read_table(cluster_summary_path)

    validate_schema(clusters, CLUSTERS_SCHEMA, "clusters.parquet")
    validate_schema(id_map, ID_MAP_SCHEMA, "id_map.parquet")
    validate_schema(genome_meta, GENOME_META_SCHEMA, "genome_meta.parquet")
    validate_schema(summary, CLUSTER_SUMMARY_SCHEMA, "cluster_summary.parquet")

    clusters_pyd = clusters.to_pydict()
    n_rows = clusters.num_rows

    # protein_uid -> (locus_tag, gene, product) lookup
    id_map_pyd = id_map.to_pydict()
    pid_to_idx = {uid: i for i, uid in enumerate(id_map_pyd["protein_uid"])}
    locus_col = id_map_pyd.get("locus_tag", [None] * id_map.num_rows)
    gene_col = id_map_pyd.get("gene", [None] * id_map.num_rows)
    product_col = id_map_pyd.get("product", [None] * id_map.num_rows)
    cds_key_col = id_map_pyd.get("cds_key", [None] * id_map.num_rows)

    # genome_uid -> genome_key lookup
    gm_pyd = genome_meta.to_pydict()
    gid_to_key = {
        gid: key for gid, key in zip(gm_pyd["genome_uid"], gm_pyd["genome_key"])
    }

    # cluster_id -> category lookup (from cluster_summary)
    summary_pyd = summary.to_pydict()
    cid_to_category = {
        cid: cat
        for cid, cat in zip(summary_pyd["cluster_id"], summary_pyd["category"])
    }

    out_locus: list[str | None] = []
    out_gene: list[str | None] = []
    out_product: list[str | None] = []
    out_genome_key: list[str] = []
    out_category: list[str] = []

    for i in range(n_rows):
        pid = clusters_pyd["protein_uid"][i]
        gid = clusters_pyd["genome_uid"][i]
        cid = clusters_pyd["cluster_id"][i]

        idx = pid_to_idx.get(pid)
        if idx is None:
            out_locus.append(None)
            out_gene.append(None)
            out_product.append(None)
        else:
            # Fall back to cds_key when locus_tag is missing/empty so the
            # column is never blank for non-NCBI assemblies.
            lt = locus_col[idx]
            if lt is None or lt == "":
                lt = cds_key_col[idx]
            out_locus.append(lt)
            out_gene.append(gene_col[idx])
            out_product.append(product_col[idx])

        gkey = gid_to_key.get(gid)
        if gkey is None:
            raise ValueError(
                f"row {i}: genome_uid={gid} has no entry in genome_meta"
            )
        out_genome_key.append(gkey)
        out_category.append(cid_to_category.get(cid, ""))

    table = pa.table(
        {
            "cluster_id":       pa.array(clusters_pyd["cluster_id"], type=pa.uint32()),
            "category":         pa.array(out_category, type=pa.string()),
            "genome_uid":       pa.array(clusters_pyd["genome_uid"], type=pa.uint16()),
            "genome_key":       pa.array(out_genome_key, type=pa.string()),
            "protein_uid":      pa.array(clusters_pyd["protein_uid"], type=pa.uint64()),
            "locus_tag":        pa.array(out_locus, type=pa.string()),
            "gene":             pa.array(out_gene, type=pa.string()),
            "product":          pa.array(out_product, type=pa.string()),
            "is_centroid":      pa.array(clusters_pyd["is_centroid"], type=pa.bool_()),
            "pct_identity_fwd": pa.array(clusters_pyd["pct_identity_fwd"], type=pa.float32()),
            "pct_identity_rev": pa.array(clusters_pyd["pct_identity_rev"], type=pa.float32()),
            "member_coverage":  pa.array(clusters_pyd["member_coverage"], type=pa.float32()),
            "rep_coverage":     pa.array(clusters_pyd["rep_coverage"], type=pa.float32()),
            "alignment_length": pa.array(clusters_pyd["alignment_length"], type=pa.uint32()),
        },
        schema=CLUSTER_LONG_SCHEMA,
    )
    validate_schema(table, CLUSTER_LONG_SCHEMA, "cluster_long.parquet")
    return table


def write_cluster_long(
    clusters_path: Path,
    id_map_path: Path,
    genome_meta_path: Path,
    cluster_summary_path: Path,
    out_path: Path,
) -> pa.Table:
    """Build and persist ``dnmb/cluster_long.parquet``."""
    table = build_cluster_long(
        clusters_path, id_map_path, genome_meta_path, cluster_summary_path,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, out_path, compression="zstd", compression_level=3)
    return table
