"""Per-cluster summary: category + representative metadata.

Joins ``clusters.parquet`` with ``id_map.parquet`` and
``presence_absence.parquet`` to produce a one-row-per-cluster view
that downstream reporting / Roary export can consume in a single read.
"""
from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .schemas import (
    CLUSTER_SUMMARY_SCHEMA,
    CLUSTERS_SCHEMA,
    ID_MAP_SCHEMA,
    PRESENCE_ABSENCE_SCHEMA,
    validate_schema,
)


def build_cluster_summary(
    clusters_path: Path,
    id_map_path: Path,
    presence_absence_path: Path,
) -> pa.Table:
    """Return the per-cluster summary table (one row per cluster)."""
    clusters = pq.read_table(clusters_path)
    id_map = pq.read_table(id_map_path)
    presence = pq.read_table(presence_absence_path)

    validate_schema(clusters, CLUSTERS_SCHEMA, "clusters.parquet")
    validate_schema(id_map, ID_MAP_SCHEMA, "id_map.parquet")
    validate_schema(presence, PRESENCE_ABSENCE_SCHEMA, "presence_absence.parquet")

    # Build lookups keyed by protein_uid. These are small (≤ few million
    # rows even at large scale) and PyArrow-native dict conversion is
    # fast enough for M4 MVP. Polars streaming join would be the hot-
    # path choice after M6 profiling.
    id_map_pydict = id_map.to_pydict()
    protein_uid_to_idx = {
        uid: i for i, uid in enumerate(id_map_pydict["protein_uid"])
    }
    gene_col = id_map_pydict.get("gene", [None] * id_map.num_rows)
    product_col = id_map_pydict.get("product", [None] * id_map.num_rows)
    locus_tag_col = id_map_pydict.get("locus_tag", [None] * id_map.num_rows)
    protein_id_col = id_map_pydict.get("protein_id", [None] * id_map.num_rows)

    # Pull representative_uid per cluster_id from the clusters table.
    # Every row for a cluster carries the same representative_uid (set at
    # engine-parse time), so `setdefault` on first encounter is sufficient.
    clusters_pydict = clusters.to_pydict()
    cluster_id_col = clusters_pydict["cluster_id"]
    rep_col = clusters_pydict["representative_uid"]

    rep_uid_per_cluster: dict[int, int] = {}
    for cid, rep_uid in zip(cluster_id_col, rep_col):
        rep_uid_per_cluster.setdefault(cid, rep_uid)

    presence_pydict = presence.to_pydict()
    pres_cluster_ids = presence_pydict["cluster_id"]
    pres_n_genomes = presence_pydict["n_genomes"]
    pres_n_sequences = presence_pydict["n_sequences"]
    pres_category = presence_pydict["category"]

    out_cluster_id: list[int] = []
    out_n_genomes: list[int] = []
    out_n_sequences: list[int] = []
    out_category: list[str] = []
    out_rep_uid: list[int] = []
    out_rep_gene: list[str | None] = []
    out_rep_product: list[str | None] = []
    out_rep_locus_tag: list[str | None] = []
    out_rep_protein_id: list[str | None] = []

    for cid, n_g, n_s, cat in zip(
        pres_cluster_ids, pres_n_genomes, pres_n_sequences, pres_category,
    ):
        rep_uid = rep_uid_per_cluster.get(cid)
        if rep_uid is None:
            # Shouldn't happen: presence_absence is derived from clusters.
            continue
        row_idx = protein_uid_to_idx.get(rep_uid)
        out_cluster_id.append(cid)
        out_n_genomes.append(n_g)
        out_n_sequences.append(n_s)
        out_category.append(cat)
        out_rep_uid.append(rep_uid)
        if row_idx is None:
            out_rep_gene.append(None)
            out_rep_product.append(None)
            out_rep_locus_tag.append(None)
            out_rep_protein_id.append(None)
        else:
            out_rep_gene.append(gene_col[row_idx])
            out_rep_product.append(product_col[row_idx])
            out_rep_locus_tag.append(locus_tag_col[row_idx])
            out_rep_protein_id.append(protein_id_col[row_idx])

    result = pa.table(
        {
            "cluster_id": pa.array(out_cluster_id, type=pa.uint32()),
            "n_genomes": pa.array(out_n_genomes, type=pa.uint16()),
            "n_sequences": pa.array(out_n_sequences, type=pa.uint32()),
            "category": pa.array(out_category, type=pa.string()),
            "representative_uid": pa.array(out_rep_uid, type=pa.uint64()),
            "representative_gene": pa.array(out_rep_gene, type=pa.string()),
            "representative_product": pa.array(out_rep_product, type=pa.string()),
            "representative_locus_tag": pa.array(out_rep_locus_tag, type=pa.string()),
            "representative_protein_id": pa.array(out_rep_protein_id, type=pa.string()),
        },
        schema=CLUSTER_SUMMARY_SCHEMA,
    )
    validate_schema(result, CLUSTER_SUMMARY_SCHEMA, "cluster_summary.parquet")
    return result


def write_cluster_summary(
    clusters_path: Path,
    id_map_path: Path,
    presence_absence_path: Path,
    out_path: Path,
) -> pa.Table:
    table = build_cluster_summary(
        clusters_path, id_map_path, presence_absence_path,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, out_path, compression="zstd", compression_level=3)
    return table
