"""Gene gain/loss inference via Fitch parsimony on a phylogenetic tree.

For each accessory cluster (not core, not absent in all), reconstructs
ancestral presence/absence states at every internal node using Fitch's
maximum parsimony algorithm, then counts per-branch gain (0→1) and
loss (1→0) events.

Fitch parsimony is preferred over Dollo for bacterial pan-genomes
because HGT allows the same gene family to be gained independently
on multiple branches — Dollo's single-gain constraint underestimates
this. Wagner parsimony (gain == loss cost) is the special case we use
here.

Requirements:
- ``dnmb/processed/phylo_tree.nwk``
- ``dnmb/processed/presence_absence.parquet``
- ``dnmb/inputs/genome_meta.parquet``

Outputs:
- ``dnmb/processed/gain_loss.parquet`` — one row per branch
  with ``parent_node, child_node, n_gained, n_lost``.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

log = logging.getLogger(__name__)


def _parse_tree(nwk_path: Path):
    """Parse Newick via Biopython, assign unique internal node names."""
    from Bio import Phylo
    tree = Phylo.read(str(nwk_path), "newick")
    internal_id = 0
    for clade in tree.find_clades(order="level"):
        if not clade.is_terminal():
            clade.name = f"N{internal_id}"
            internal_id += 1
    return tree


def _fitch_parsimony_one_cluster(tree, tip_states: dict[str, int]):
    """Run Fitch parsimony for one binary character (presence/absence).

    Returns ``{(parent_name, child_name): (gain, loss)}`` for every
    branch. gain=1 means 0→1 transition, loss=1 means 1→0.

    Bottom-up pass: each node gets a SET of possible states.
    Top-down pass: each node gets a single resolved state.
    """
    # Bottom-up: compute state sets
    state_sets: dict[str, set[int]] = {}

    def _bottom_up(clade):
        if clade.is_terminal():
            state_sets[clade.name] = {tip_states.get(clade.name, 0)}
            return
        child_sets = []
        for child in clade.clades:
            _bottom_up(child)
            child_sets.append(state_sets[child.name])
        # Fitch rule: intersection if non-empty, else union
        inter = child_sets[0]
        for s in child_sets[1:]:
            inter = inter & s
        if inter:
            state_sets[clade.name] = inter
        else:
            union = set()
            for s in child_sets:
                union = union | s
            state_sets[clade.name] = union

    _bottom_up(tree.root)

    # Top-down: resolve ambiguous states
    resolved: dict[str, int] = {}

    def _top_down(clade, parent_state: int | None):
        ss = state_sets[clade.name]
        if parent_state is not None and parent_state in ss:
            resolved[clade.name] = parent_state
        else:
            # Pick 0 (absent) as default when ambiguous — conservative
            # for gain counting (avoids inflating gains).
            resolved[clade.name] = min(ss)
        for child in clade.clades:
            _top_down(child, resolved[clade.name])

    _top_down(tree.root, None)

    # Count transitions per branch
    transitions: dict[tuple[str, str], tuple[int, int]] = {}

    def _count(parent, clade):
        if parent is not None:
            p_state = resolved[parent.name]
            c_state = resolved[clade.name]
            gain = 1 if p_state == 0 and c_state == 1 else 0
            loss = 1 if p_state == 1 and c_state == 0 else 0
            transitions[(parent.name, clade.name)] = (gain, loss)
        for child in clade.clades:
            _count(clade, child)

    _count(None, tree.root)
    return transitions


def compute_gain_loss(
    *,
    tree_path: Path,
    presence_absence_path: Path,
    genome_meta_path: Path,
    out_path: Path,
) -> pa.Table:
    """Infer gene gain/loss events via Fitch parsimony."""
    tree = _parse_tree(tree_path)

    pa_table = pq.read_table(presence_absence_path)
    meta = pq.read_table(
        genome_meta_path, columns=["genome_uid", "genome_key"]
    ).to_pydict()
    gid_to_key = {
        gid: gkey for gid, gkey in zip(meta["genome_uid"], meta["genome_key"])
    }
    all_tips = {t.name for t in tree.get_terminals()}

    pa_pyd = pa_table.to_pydict()

    # Build cluster_id → {genome_key: 1} presence map
    cluster_presence: dict[int, dict[str, int]] = {}
    for cid, bitmap in zip(pa_pyd["cluster_id"], pa_pyd["genome_bitmap"]):
        if not bitmap:
            continue
        # Vectorized set-bit extraction: view uint64 words as little-endian
        # bytes, then unpackbits with bitorder='little' so bit position b
        # within word w lands at index w*64 + b — matching the original
        # word_idx*64+bit mapping.
        bm = np.asarray(bitmap, dtype="<u8")
        bits = np.unpackbits(bm.view(np.uint8), bitorder="little")
        present_gids = np.flatnonzero(bits).tolist()
        present_keys = {
            gid_to_key[g] for g in present_gids if g in gid_to_key
        }
        # Skip clusters present in ALL or NONE of the tree tips
        in_tree = present_keys & all_tips
        if not in_tree or in_tree == all_tips:
            continue
        cluster_presence[cid] = {k: (1 if k in in_tree else 0) for k in all_tips}

    log.info("gain_loss: %d accessory clusters for Fitch parsimony", len(cluster_presence))

    # Aggregate over all clusters
    total_gain: dict[tuple[str, str], int] = {}
    total_loss: dict[tuple[str, str], int] = {}

    for cid, tip_states in cluster_presence.items():
        transitions = _fitch_parsimony_one_cluster(tree, tip_states)
        for branch, (g, l) in transitions.items():
            total_gain[branch] = total_gain.get(branch, 0) + g
            total_loss[branch] = total_loss.get(branch, 0) + l

    # Build output table — one row per branch
    branches = sorted(total_gain.keys())
    rows_parent = [b[0] for b in branches]
    rows_child = [b[1] for b in branches]
    rows_gain = [total_gain[b] for b in branches]
    rows_loss = [total_loss[b] for b in branches]

    table = pa.table({
        "parent_node": pa.array(rows_parent, type=pa.string()),
        "child_node":  pa.array(rows_child, type=pa.string()),
        "n_gained":    pa.array(rows_gain, type=pa.uint32()),
        "n_lost":      pa.array(rows_loss, type=pa.uint32()),
    })

    from .io_utils import atomic_write_table
    atomic_write_table(table, out_path)
    total_g = sum(rows_gain)
    total_l = sum(rows_loss)
    log.info("gain/loss (Fitch): %d branches, +%d gains, -%d losses",
             len(branches), total_g, total_l)
    return table
