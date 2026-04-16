"""Fast EggNOG annotation via DIAMOND direct + SQLite lookup.

Bypasses eggnog-mapper's Python overhead entirely. Steps:

1. DIAMOND blastp ``--fast`` against the cached ``eggnog_proteins.dmnd``
   (~9 GB) — takes ~2 min on 33k bacterial proteins vs eggnog-mapper's
   ~30 min.
2. Parse top hits (one per query).
3. Batch-lookup each hit's ``(COG_categories, kegg_ko, kegg_cog)`` from
   ``eggnog.db`` (SQLite, indexed).
4. Write ``dnmb/processed/eggnog_annotations.parquet`` with per-CDS
   COG category + KEGG ko.

Requires the EggNOG database cached at a known path (same as
DNMBsuite's ``~/.dnmb-cache/db_modules/eggnog/data/``).
"""
from __future__ import annotations

import logging
import os
import shutil
import sqlite3
import subprocess
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

log = logging.getLogger(__name__)

_DEFAULT_DB_DIRS = [
    Path.home() / ".dnmb-cache" / "db_modules" / "eggnog" / "data",
    Path.home() / ".local" / "share" / "eggnog-mapper",
]


def _find_eggnog_db(override: Path | None = None) -> tuple[Path, Path]:
    """Return ``(dmnd_path, sqlite_path)`` or raise."""
    candidates = [override] if override else _DEFAULT_DB_DIRS
    for d in candidates:
        if d is None:
            continue
        dmnd = d / "eggnog_proteins.dmnd"
        db = d / "eggnog.db"
        if dmnd.exists() and db.exists():
            return dmnd, db
    raise RuntimeError(
        "EggNOG database not found. Expected eggnog_proteins.dmnd + "
        "eggnog.db in one of: "
        + ", ".join(str(d) for d in _DEFAULT_DB_DIRS)
        + ". Download via: download_eggnog_data.py -y --data_dir <path>"
    )


def _run_diamond_search(
    query_fasta: Path,
    dmnd_db: Path,
    out_tsv: Path,
    threads: int,
) -> None:
    """DIAMOND blastp --fast against the EggNOG protein DB."""
    if shutil.which("diamond") is None:
        raise RuntimeError("diamond not found on PATH")

    threads_arg = str(threads) if threads > 0 else str(os.cpu_count() or 1)
    # Default outfmt 6 = standard BLAST tabular (12 cols). Don't
    # pass field names — older DIAMOND + eggnog-format dmnd DBs
    # silently produce 0 hits when custom field specs conflict with
    # the internal DB format.
    cmd = [
        "diamond", "blastp",
        "--db", str(dmnd_db),
        "--query", str(query_fasta),
        "--out", str(out_tsv),
        "--outfmt", "6",
        "--max-target-seqs", "1",
        "--evalue", "1e-5",
        "--threads", threads_arg,
    ]
    log.info("eggnog-fast DIAMOND: %s", " ".join(cmd[:10]))
    try:
        subprocess.run(
            cmd, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
        raise RuntimeError(f"diamond blastp failed: {stderr}") from exc


def _parse_diamond_hits(tsv_path: Path) -> dict[int, str]:
    """Parse DIAMOND outfmt 6 → {query_uid: hit_name}."""
    hits: dict[int, str] = {}
    with open(tsv_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                qid = int(parts[0])
            except ValueError:
                continue
            hits[qid] = parts[1]
    return hits


def _batch_lookup_cog_kegg(
    hit_names: list[str], db_path: Path,
) -> dict[str, tuple[str, str]]:
    """Batch-lookup ``(COG_category_letter, kegg_ko)`` from eggnog.db.

    The ``prots.ogs`` field contains a comma-separated list of ortholog
    group assignments like ``"COG3209@1,COG3209@2,1TR8F@1239"``. We
    extract the first ``COG\\d+`` ID, then look up its single-letter
    functional category from the ``og`` table. ``kegg_ko`` comes from
    ``prots.kegg_ko`` directly (often empty for non-KEGG organisms).

    Returns ``{hit_name: (cog_letter, kegg_ko)}``.
    """
    import re
    _COG_RE = re.compile(r"(COG\d+)@")

    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()

    # Step 1: get ogs + kegg_ko from prots
    name_to_ogs: dict[str, tuple[str, str]] = {}
    batch_size = 500
    names = list(set(hit_names))

    for i in range(0, len(names), batch_size):
        batch = names[i : i + batch_size]
        placeholders = ",".join("?" * len(batch))
        cur.execute(
            f"SELECT name, ogs, kegg_ko FROM prots "
            f"WHERE name IN ({placeholders})",
            batch,
        )
        for row in cur.fetchall():
            name_to_ogs[row[0]] = (row[1] or "", row[2] or "")

    # Step 2: extract unique COG IDs and batch-lookup their categories
    all_cog_ids: set[str] = set()
    name_to_first_cog: dict[str, str] = {}
    for name, (ogs, _) in name_to_ogs.items():
        match = _COG_RE.search(ogs)
        if match:
            cog_id = match.group(1)
            name_to_first_cog[name] = cog_id
            all_cog_ids.add(cog_id)

    cog_to_cat: dict[str, str] = {}
    cog_list = list(all_cog_ids)
    for i in range(0, len(cog_list), batch_size):
        batch = cog_list[i : i + batch_size]
        placeholders = ",".join("?" * len(batch))
        cur.execute(
            f"SELECT og, COG_categories FROM og "
            f"WHERE og IN ({placeholders}) AND level = '1'",
            batch,
        )
        for row in cur.fetchall():
            cog_to_cat[row[0]] = row[1] or ""
        # Fallback: level '2' if level '1' didn't match
        missing = [c for c in batch if c not in cog_to_cat]
        if missing:
            placeholders2 = ",".join("?" * len(missing))
            cur.execute(
                f"SELECT og, COG_categories FROM og "
                f"WHERE og IN ({placeholders2})",
                missing,
            )
            for row in cur.fetchall():
                if row[0] not in cog_to_cat:
                    cog_to_cat[row[0]] = row[1] or ""

    conn.close()

    # Step 3: assemble final result
    result: dict[str, tuple[str, str]] = {}
    for name in names:
        ogs_kegg = name_to_ogs.get(name)
        if ogs_kegg is None:
            result[name] = ("", "")
            continue
        kegg_ko = ogs_kegg[1]
        cog_id = name_to_first_cog.get(name, "")
        cog_cat = cog_to_cat.get(cog_id, "")
        result[name] = (cog_cat, kegg_ko)

    return result


def run_eggnog_fast(
    *,
    input_fasta: Path,
    clusters_path: Path,
    raw_dir: Path,
    processed_dir: Path,
    threads: int = 0,
    eggnog_data_dir: Path | None = None,
) -> pa.Table:
    """Run the fast EggNOG annotation pipeline.

    Outputs ``dnmb/processed/eggnog_annotations.parquet`` with columns:
    ``protein_uid, cluster_id, eggnog_hit, cog_category, kegg_ko``.
    """
    dmnd_db, sqlite_db = _find_eggnog_db(eggnog_data_dir)
    log.info("eggnog-fast: using DB at %s", dmnd_db.parent)

    eggnog_raw = raw_dir / "eggnog"
    eggnog_raw.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)

    hits_tsv = eggnog_raw / "diamond_hits.tsv"
    _run_diamond_search(input_fasta, dmnd_db, hits_tsv, threads)

    hits = _parse_diamond_hits(hits_tsv)
    log.info("eggnog-fast: %d queries → %d hits", 0, len(hits))

    unique_targets = list(set(hits.values()))
    cog_kegg = _batch_lookup_cog_kegg(unique_targets, sqlite_db)

    clusters = pq.read_table(
        clusters_path, columns=["protein_uid", "cluster_id"]
    ).to_pydict()

    protein_uids = clusters["protein_uid"]
    cluster_ids = clusters["cluster_id"]
    n = len(protein_uids)

    eggnog_hits: list[str | None] = []
    cog_cats: list[str] = []
    kegg_kos: list[str] = []

    for i in range(n):
        uid = protein_uids[i]
        hit_name = hits.get(uid)
        if hit_name is None:
            eggnog_hits.append(None)
            cog_cats.append("")
            kegg_kos.append("")
        else:
            eggnog_hits.append(hit_name)
            ck = cog_kegg.get(hit_name, ("", ""))
            cog_cats.append(ck[0])
            kegg_kos.append(ck[1])

    table = pa.table({
        "protein_uid":  pa.array(protein_uids, type=pa.uint64()),
        "cluster_id":   pa.array(cluster_ids, type=pa.uint32()),
        "eggnog_hit":   pa.array(eggnog_hits, type=pa.string()),
        "cog_category": pa.array(cog_cats, type=pa.string()),
        "kegg_ko":      pa.array(kegg_kos, type=pa.string()),
    })

    out_path = processed_dir / "eggnog_annotations.parquet"
    pq.write_table(table, out_path, compression="zstd", compression_level=3)
    log.info("eggnog-fast: wrote %s (%d rows)", out_path, n)
    return table
