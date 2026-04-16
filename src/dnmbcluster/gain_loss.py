"""Gene gain/loss inference on a phylogenetic tree via Dollo parsimony.

Maps the binary presence/absence of each cluster onto the core-gene
phylogeny and counts per-branch gain and loss events under the Dollo
model (a gene can be gained once on the tree but lost independently
on any number of branches). This is the standard approach for
accessory-genome dynamics in bacterial pan-genome papers.

Requirements:
- ``dnmb/processed/phylo_tree.nwk`` (from ``--phylo`` stage)
- ``dnmb/processed/presence_absence.parquet``
- ``dnmb/inputs/genome_meta.parquet``

Outputs:
- ``dnmb/processed/gain_loss.parquet`` — one row per internal branch
  with ``parent_node, child_node, n_gained, n_lost`` counts.

The R visualization layer picks this up and annotates the ggtree
phylo_tree plot with +N / -N labels on each branch.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

log = logging.getLogger(__name__)


def _parse_newick_simple(nwk_path: Path):
    """Parse a Newick tree via ete3 or Biopython's Phylo."""
    try:
        from Bio import Phylo
        import io
        tree = Phylo.read(str(nwk_path), "newick")
        return tree
    except ImportError:
        raise RuntimeError(
            "Biopython required for gain/loss inference (Bio.Phylo)"
        )


def _get_tip_names(tree) -> list[str]:
    """Return leaf names from a Biopython Phylo tree."""
    return [tip.name for tip in tree.get_terminals() if tip.name]


def _dollo_parsimony(
    tree, cluster_presence: dict[int, set[str]],
) -> list[dict]:
    """Run Dollo parsimony for each cluster over the tree.

    Returns a list of dicts ``{branch: (parent_name, child_name),
    gains: int, losses: int}`` aggregated over all clusters.
    """
    from Bio import Phylo

    # Label internal nodes with unique IDs. IQ-TREE writes support
    # values (e.g. "100/100") as node labels, which aren't unique
    # across branches — replace them with deterministic N0, N1, ...
    # while the original support string is only used by ggtree.
    internal_id = 0
    for clade in tree.find_clades(order="level"):
        is_tip = clade.is_terminal()
        if not is_tip:
            clade.name = f"N{internal_id}"
            internal_id += 1

    # Build branch list: (parent, child) for every edge
    branches: list[tuple[str, str]] = []
    def _walk(parent, clade):
        if parent is not None:
            branches.append((parent.name, clade.name))
        for child in clade.clades:
            _walk(clade, child)
    _walk(None, tree.root)

    # Per-branch gain/loss counters
    gain_count: dict[tuple[str, str], int] = {b: 0 for b in branches}
    loss_count: dict[tuple[str, str], int] = {b: 0 for b in branches}

    # For each cluster, do a bottom-up + top-down Dollo inference.
    # Dollo: gene gained at the LCA of all tips that have it;
    # lost independently on every branch leading to a tip that
    # lacks it below the LCA.
    all_tips = set(_get_tip_names(tree))

    for cid, present_tips in cluster_presence.items():
        present_in_tree = present_tips & all_tips
        if not present_in_tree or present_in_tree == all_tips:
            continue  # skip core / absent clusters

        # Find LCA of present tips
        if len(present_in_tree) == 1:
            lca = tree.find_any(name=list(present_in_tree)[0])
        else:
            tips_list = [tree.find_any(name=t) for t in present_in_tree]
            tips_list = [t for t in tips_list if t is not None]
            if len(tips_list) < 2:
                continue
            lca = tree.common_ancestor(tips_list)

        if lca is None:
            continue

        # Mark gain at the branch leading TO the LCA
        lca_name = lca.name

        # Find the parent of LCA
        path = tree.get_path(lca)
        if path and len(path) >= 2:
            parent_of_lca = path[-2].name if hasattr(path[-2], 'name') else None
        else:
            parent_of_lca = None

        if parent_of_lca:
            branch_key = (parent_of_lca, lca_name)
            if branch_key in gain_count:
                gain_count[branch_key] += 1

        # Mark losses: for every tip below the LCA that LACKS the
        # gene, trace up to the LCA and mark the first branch as a loss.
        for tip in lca.get_terminals():
            if tip.name not in present_in_tree:
                # Find the branch from lca downward that leads to this tip
                tip_path = tree.get_path(tip)
                # The branch right below the LCA on the way to this tip
                for i, node in enumerate(tip_path):
                    if node == lca and i + 1 < len(tip_path):
                        child_of_lca = tip_path[i + 1]
                        branch_key = (lca_name, child_of_lca.name)
                        if branch_key in loss_count:
                            loss_count[branch_key] += 1
                        break

    rows = []
    for (parent, child) in branches:
        rows.append({
            "parent_node": parent,
            "child_node": child,
            "n_gained": gain_count[(parent, child)],
            "n_lost": loss_count[(parent, child)],
        })
    return rows


def compute_gain_loss(
    *,
    tree_path: Path,
    presence_absence_path: Path,
    genome_meta_path: Path,
    out_path: Path,
) -> pa.Table:
    """Infer gene gain/loss events and write ``gain_loss.parquet``."""
    tree = _parse_newick_simple(tree_path)

    pa_table = pq.read_table(presence_absence_path)
    meta = pq.read_table(
        genome_meta_path, columns=["genome_uid", "genome_key"]
    ).to_pydict()
    gid_to_key = {
        gid: gkey for gid, gkey in zip(meta["genome_uid"], meta["genome_key"])
    }

    # Build cluster_id → set of genome_keys
    pa_pyd = pa_table.to_pydict()
    cluster_presence: dict[int, set[str]] = {}
    for cid, bitmap in zip(pa_pyd["cluster_id"], pa_pyd["genome_bitmap"]):
        present_gids: set[int] = set()
        for word_idx, word in enumerate(bitmap):
            for bit in range(64):
                if word & (1 << bit):
                    present_gids.add(word_idx * 64 + bit)
        cluster_presence[cid] = {
            gid_to_key[g] for g in present_gids if g in gid_to_key
        }

    rows = _dollo_parsimony(tree, cluster_presence)

    table = pa.table({
        "parent_node": pa.array([r["parent_node"] for r in rows], type=pa.string()),
        "child_node":  pa.array([r["child_node"] for r in rows], type=pa.string()),
        "n_gained":    pa.array([r["n_gained"] for r in rows], type=pa.uint32()),
        "n_lost":      pa.array([r["n_lost"] for r in rows], type=pa.uint32()),
    })

    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, out_path, compression="zstd", compression_level=3)
    log.info("gain/loss: %d branches, %d total gains, %d total losses",
             len(rows),
             sum(r["n_gained"] for r in rows),
             sum(r["n_lost"] for r in rows))
    return table
