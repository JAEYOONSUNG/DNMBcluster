"""Roary-compatible ``gene_presence_absence.csv`` exporter.

Roary's CSV format is the de-facto pan-genome output standard —
Scoary, Phandango, Piggy, and Coinfinder all consume it. Producing it
means DNMBcluster outputs are a drop-in for that ecosystem regardless
of which clustering engine was used.

Format: 14 fixed metadata columns + one column per genome whose cell
value is the locus_tag(s) of that genome's CDS(es) in the cluster
(semicolon-separated when multiple).

We also emit the ``.Rtab`` binary companion (tab-separated 0/1 matrix)
that Scoary and downstream tools expect.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pyarrow.parquet as pq

from .schemas import (
    CLUSTER_SUMMARY_SCHEMA,
    CLUSTERS_SCHEMA,
    GENOME_META_SCHEMA,
    ID_MAP_SCHEMA,
    validate_schema,
)

ROARY_METADATA_COLUMNS: tuple[str, ...] = (
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Order with Fragment",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
)


def write_roary_csv(
    clusters_path: Path,
    id_map_path: Path,
    genome_meta_path: Path,
    cluster_summary_path: Path,
    *,
    csv_out: Path,
    rtab_out: Path,
) -> tuple[int, int]:
    """Write Roary ``gene_presence_absence.csv`` and ``.Rtab``.

    Returns ``(n_clusters, n_genomes)``.
    """
    clusters = pq.read_table(clusters_path)
    id_map = pq.read_table(id_map_path)
    genome_meta = pq.read_table(genome_meta_path)
    summary = pq.read_table(cluster_summary_path)

    validate_schema(clusters, CLUSTERS_SCHEMA, "clusters.parquet")
    validate_schema(id_map, ID_MAP_SCHEMA, "id_map.parquet")
    validate_schema(genome_meta, GENOME_META_SCHEMA, "genome_meta.parquet")
    validate_schema(summary, CLUSTER_SUMMARY_SCHEMA, "cluster_summary.parquet")

    # Ordered list of genome column names — use genome_key as the header.
    meta_pydict = genome_meta.to_pydict()
    genome_uids = meta_pydict["genome_uid"]
    genome_keys = meta_pydict["genome_key"]
    # Order by genome_uid for deterministic column order
    genome_order: list[tuple[int, str]] = sorted(
        zip(genome_uids, genome_keys), key=lambda x: x[0],
    )
    genome_uid_to_col: dict[int, int] = {
        gid: i for i, (gid, _) in enumerate(genome_order)
    }
    column_headers = [gkey for _, gkey in genome_order]
    n_genomes = len(column_headers)

    # Build a lookup: protein_uid -> locus_tag (or fallback)
    id_pydict = id_map.to_pydict()
    protein_uid_to_locus: dict[int, str] = {}
    for uid, locus_tag, cds_key in zip(
        id_pydict["protein_uid"],
        id_pydict["locus_tag"],
        id_pydict["cds_key"],
    ):
        protein_uid_to_locus[uid] = locus_tag if locus_tag else cds_key

    # Group clusters: for each cluster_id, for each genome_uid, collect
    # locus_tags of member CDSs.
    clusters_pydict = clusters.to_pydict()
    cluster_to_genome_loci: dict[int, dict[int, list[str]]] = {}
    cluster_to_lengths: dict[int, list[int]] = {}

    # Need per-member lengths; we read them from the id_map aa_length column.
    uid_to_length: dict[int, int] = {}
    for uid, aa_len in zip(id_pydict["protein_uid"], id_pydict["aa_length"]):
        if aa_len is not None:
            # Roary reports "group size nuc" — multiply by 3 for AA length.
            uid_to_length[uid] = int(aa_len) * 3

    for cid, puid, guid in zip(
        clusters_pydict["cluster_id"],
        clusters_pydict["protein_uid"],
        clusters_pydict["genome_uid"],
    ):
        cluster_to_genome_loci.setdefault(cid, {}).setdefault(guid, []).append(
            protein_uid_to_locus.get(puid, str(puid))
        )
        length = uid_to_length.get(puid)
        if length is not None:
            cluster_to_lengths.setdefault(cid, []).append(length)

    # Cluster summary keyed by cluster_id
    summary_pydict = summary.to_pydict()
    summary_by_cid: dict[int, dict] = {}
    for i, cid in enumerate(summary_pydict["cluster_id"]):
        summary_by_cid[cid] = {
            "n_genomes": summary_pydict["n_genomes"][i],
            "n_sequences": summary_pydict["n_sequences"][i],
            "category": summary_pydict["category"][i],
            "gene": summary_pydict["representative_gene"][i],
            "product": summary_pydict["representative_product"][i],
        }

    csv_out.parent.mkdir(parents=True, exist_ok=True)
    rtab_out.parent.mkdir(parents=True, exist_ok=True)

    header = list(ROARY_METADATA_COLUMNS) + column_headers
    sorted_cluster_ids = sorted(summary_by_cid.keys())

    n_clusters = 0
    with open(csv_out, "w", newline="") as csv_fh, open(rtab_out, "w") as rtab_fh:
        writer = csv.writer(csv_fh, quoting=csv.QUOTE_MINIMAL)
        writer.writerow(header)
        rtab_fh.write("Gene\t" + "\t".join(column_headers) + "\n")

        for cid in sorted_cluster_ids:
            summary_row = summary_by_cid[cid]
            gene = summary_row.get("gene") or f"group_{cid}"
            annotation = summary_row.get("product") or ""
            n_iso = summary_row["n_genomes"]
            n_seq = summary_row["n_sequences"]
            avg_seq = f"{n_seq / n_iso:.2f}" if n_iso else "0"

            lengths = cluster_to_lengths.get(cid, [])
            if lengths:
                min_len = min(lengths)
                max_len = max(lengths)
                avg_len = sum(lengths) / len(lengths)
            else:
                min_len = max_len = 0
                avg_len = 0.0

            metadata = [
                gene,
                "",                 # Non-unique Gene name
                annotation,
                str(n_iso),
                str(n_seq),
                avg_seq,
                "",                 # Genome Fragment
                "",                 # Order within Fragment
                "",                 # Accessory Fragment
                "",                 # Accessory Order with Fragment
                "",                 # QC
                str(min_len),
                str(max_len),
                f"{avg_len:.2f}",
            ]

            genome_cells: list[str] = []
            rtab_cells: list[str] = []
            loci_by_genome = cluster_to_genome_loci.get(cid, {})
            for gid, _ in genome_order:
                if gid in loci_by_genome:
                    genome_cells.append(";".join(loci_by_genome[gid]))
                    rtab_cells.append("1")
                else:
                    genome_cells.append("")
                    rtab_cells.append("0")

            writer.writerow(metadata + genome_cells)
            rtab_fh.write(gene + "\t" + "\t".join(rtab_cells) + "\n")
            n_clusters += 1

    return n_clusters, n_genomes
