"""Lightweight functional categorization from GenBank product annotations.

Maps each CDS's ``product`` field to one of ~15 broad functional
categories via keyword matching. This is NOT a rigorous COG/KEGG
classification — it's a fast, zero-dependency first-pass that works
on any input immediately because the data is already in
``id_map.parquet``. For a proper functional annotation, run
eggnog-mapper externally and feed the output to the dedicated
eggnog parser (TBD).

The category scheme follows a simplified COG-inspired grouping that
captures the major functional themes visible in typical bacterial
genomes: metabolism, transport, regulation, replication, cell wall,
mobile elements, hypotheticals, and a catch-all "Other".

Outputs:

- ``dnmb/processed/functional_categories.parquet`` — one row per CDS
  with ``protein_uid``, ``cluster_id``, ``category``, ``functional_class``
- R visualization: stacked bar of functional classes by pan-genome
  category (core / accessory / unique)
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Keyword → functional class mapping
# ---------------------------------------------------------------------------

# Order matters: first match wins. Patterns are compiled once at import
# time for speed over the full CDS set (~30k rows on Geobacillus-10).

_RULES: list[tuple[str, re.Pattern]] = [
    ("Hypothetical",     re.compile(r"hypothetical|uncharacterized|DUF\d|unknown function", re.I)),
    ("Mobile element",   re.compile(r"transpos|integrase|phage|insertion.?seq|IS\d|recombinase.*(site|phage)", re.I)),
    ("Transcription",    re.compile(r"transcription|sigma.?factor|anti.?sigma|RNA polymerase|repressor|activator.*transcri", re.I)),
    ("Translation",      re.compile(r"ribosom|tRNA|translation|elongation factor|aminoacyl.*tRNA|peptide chain release", re.I)),
    ("DNA replication",  re.compile(r"DNA polymer|DNA gyrase|helicase|topoisomerase|DNA ligase|primase|replicat|recombination|DNA repair|exonuclease|endonuclease", re.I)),
    ("Cell wall",        re.compile(r"peptidoglycan|murein|penicillin|transpeptidase|cell wall|lipopolysaccharide|LPS|teichoic", re.I)),
    ("Membrane/Transport", re.compile(r"ABC transporter|permease|channel|porin|efflux|MFS|symporter|antiporter|import|export.*protein|secretion|type.?[IVX]+.*secret", re.I)),
    ("Signal transduction", re.compile(r"sensor.*kinase|response regulator|two.?component|signal.*transduc|histidine kinase|GGDEF|diguanylate", re.I)),
    ("Amino acid metabolism", re.compile(r"amino.?acid|aminotransferase|deaminase|decarboxylase|synthase.*amino|tryptophan|glutamate|aspartate|lysine|methionine|cysteine|serine|threonine|alanine|valine|leucine|isoleucine|proline|arginine|histidine.*bio", re.I)),
    ("Carbohydrate metabolism", re.compile(r"glycos|sugar|glucose|galactose|xylose|manno|fructose|PTS|phosphotransferase.*sugar|cellulose|starch|amylase|glucan|glycolytic|pentose", re.I)),
    ("Energy metabolism", re.compile(r"cytochrome|oxidoreductase|dehydrogenase|NADH|ATP synthase|electron.*transfer|ferredoxin|oxidase.*terminal|respir|ubiquinone|menaquinone", re.I)),
    ("Lipid metabolism",  re.compile(r"lipid|fatty.?acid|acyl|phospholipid|lipase|esterase.*lip|beta.?oxidation|CoA.*fatty", re.I)),
    ("Cofactor/Vitamin",  re.compile(r"cobalamin|biotin|thiamin|riboflavin|folate|NAD|FAD|coenzyme|iron.?sulfur|molybdo|heme.*bio", re.I)),
    ("Stress response",   re.compile(r"chaperone|heat.?shock|cold.?shock|stress|catalase|superoxide|peroxidase|thioredoxin|universal stress|Clp|GroE|DnaK|DnaJ", re.I)),
    ("Other",             re.compile(r".*")),  # catch-all
]


def classify_product(product: str | None) -> str:
    """Return the first matching functional class for a product string."""
    if not product or not product.strip():
        return "Hypothetical"
    for cls, pat in _RULES:
        if pat.search(product):
            return cls
    return "Other"


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def compute_functional_categories(
    *,
    clusters_path: Path,
    id_map_path: Path,
    cluster_summary_path: Path,
    out_path: Path,
) -> pa.Table:
    """Classify every CDS by product-keyword functional category.

    Writes ``functional_categories.parquet`` with columns:
    ``protein_uid, cluster_id, pan_category, functional_class``.

    Returns the Arrow table.
    """
    clusters = pq.read_table(
        clusters_path, columns=["protein_uid", "cluster_id"]
    ).to_pydict()
    id_map = pq.read_table(
        id_map_path, columns=["protein_uid", "product"]
    ).to_pydict()
    summary = pq.read_table(
        cluster_summary_path, columns=["cluster_id", "category"]
    ).to_pydict()

    cid_to_cat: dict[int, str] = {
        cid: cat for cid, cat in zip(summary["cluster_id"], summary["category"])
    }
    uid_to_product: dict[int, str | None] = {
        uid: prod for uid, prod in zip(id_map["protein_uid"], id_map["product"])
    }

    protein_uids = clusters["protein_uid"]
    cluster_ids = clusters["cluster_id"]
    n = len(protein_uids)

    pan_categories: list[str] = []
    functional_classes: list[str] = []
    for i in range(n):
        pan_categories.append(cid_to_cat.get(cluster_ids[i], ""))
        functional_classes.append(
            classify_product(uid_to_product.get(protein_uids[i]))
        )

    table = pa.table({
        "protein_uid":     pa.array(protein_uids, type=pa.uint64()),
        "cluster_id":      pa.array(cluster_ids, type=pa.uint32()),
        "pan_category":    pa.array(pan_categories, type=pa.string()),
        "functional_class": pa.array(functional_classes, type=pa.string()),
    })

    from .io_utils import atomic_write_table
    atomic_write_table(table, out_path)
    log.info(
        "functional categories: %d CDS, %d classes",
        n, len(set(functional_classes)),
    )
    return table
