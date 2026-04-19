"""Nucleotide CDS FASTA export for codon-aware downstream stages.

HyPhy FEL/BUSTED/aBSREL and codeml M0 need codon alignments, which the
R `back_translate_og_codons()` helper builds by matching AA alignment
columns to nucleotide CDS sequences keyed by ``protein_uid``.

The core pipeline writes ``proteins.faa`` (translations) when
``--level protein``; this module emits a *parallel* ``cds.fna`` keyed
by the same ``protein_uid`` by reparsing the GenBank folder and joining
against ``id_map.parquet`` via ``(genome_uid, cds_key)``.

Why reparse rather than carrying ``nt_sequence`` in ``gene_table``?
The gene table stores exactly one sequence per CDS (selected by
``--level``) to keep the hot-path parquet slim — a dual schema would
bloat every protein-level run. This module is the opt-in sidecar.
"""
from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pyarrow.parquet as pq

from .genbank import _iter_records, find_genbank_files
from .ids import resolve_cds_key, resolve_genome_key  # noqa: F401 (resolve_genome_key used by _resolve_genome_key_for_file)

log = logging.getLogger(__name__)


# Worker-global lookup table populated by ``_init_worker`` at pool
# startup. Avoids pickling the full lookup dict per task — critical on
# large runs where the dict can be hundreds of MB.
_WORKER_LOOKUP: dict[str, dict[str, int]] | None = None


def _init_worker(lookup_by_genome: dict[str, dict[str, int]]) -> None:
    global _WORKER_LOOKUP
    _WORKER_LOOKUP = lookup_by_genome


def _process_one_file(file_path: Path) -> tuple[str, bytes, int]:
    """Worker: emit FASTA bytes for CDS in this file.

    Single pass over the GenBank file — collects genome_key metadata
    (DBLINK / organism / strain / first locus_version) alongside CDS
    candidates, then resolves genome_key and filters against
    ``_WORKER_LOOKUP`` once at the end. Walking the file only once
    halves BioPython parse cost for large WGS assemblies.

    Returns ``(file_name, fasta_bytes, n_written)`` so the driver can
    order outputs by filename (determinism) and sum written counts.
    """
    assert _WORKER_LOOKUP is not None, "worker initializer not run"

    # Genome-key inputs (accumulated across all records, just like
    # parse_folder — see _resolve_genome_key_for_file's docstring).
    dblink_entries: list[str] = []
    organism: str | None = None
    strain: str | None = None
    first_locus_version: str | None = None

    # CDS candidates collected during the walk. We defer emission until
    # after genome_key is known so that files whose key isn't in the
    # id_map don't pay the nt-extraction cost.
    candidates: list[tuple[str, str]] = []  # (cds_key, nt_sequence)

    gene_ordinals: dict[str, int] = {}
    locus_tag_occurrences: dict[str, int] = {}

    try:
        for record in _iter_records(file_path):
            ann = record.annotations or {}
            dblink_entries.extend(ann.get("dblink") or [])
            if organism is None:
                organism = ann.get("organism") or None
            if first_locus_version is None:
                first_locus_version = record.id or record.name or None
            if strain is None:
                for feat in record.features:
                    if feat.type == "source":
                        strain_vals = feat.qualifiers.get("strain")
                        if strain_vals:
                            strain = strain_vals[0]
                        break

            if record.seq is None:
                continue
            contig = record.id or record.name or file_path.stem
            for feat in record.features:
                if feat.type != "CDS":
                    continue
                quals = feat.qualifiers
                is_pseudo = "pseudo" in quals or "pseudogene" in quals
                if is_pseudo or not quals.get("translation"):
                    # Must mirror parse_folder's protein-level skip
                    # rule so cds_key ordinals align.
                    continue

                locus_tag = (quals.get("locus_tag") or [None])[0]
                protein_id = (quals.get("protein_id") or [None])[0]
                gene_name = (quals.get("gene") or [None])[0]

                if locus_tag:
                    seen = locus_tag_occurrences.get(locus_tag, 0)
                    locus_tag_occurrences[locus_tag] = seen + 1
                    if seen > 0:
                        locus_tag = f"{locus_tag}__dup{seen}"

                gene_key_for_ordinal = gene_name or ""
                gene_ordinals[gene_key_for_ordinal] = (
                    gene_ordinals.get(gene_key_for_ordinal, 0) + 1
                )
                gene_ordinal = gene_ordinals[gene_key_for_ordinal]

                try:
                    start = int(feat.location.start)
                    end = int(feat.location.end)
                except (TypeError, AttributeError):
                    start, end = 0, 0
                strand = (
                    int(feat.location.strand or 0)
                    if feat.location is not None else 0
                )

                cds_key, _ = resolve_cds_key(
                    locus_tag=locus_tag,
                    protein_id=protein_id,
                    gene=gene_name,
                    gene_ordinal=gene_ordinal,
                    contig=contig,
                    start=start,
                    end=end,
                    strand=strand,
                )
                try:
                    nt = str(feat.extract(record.seq)).upper()
                except Exception as exc:
                    log.warning(
                        "cds_fasta: skip malformed CDS in %s: %s",
                        file_path.name, exc,
                    )
                    continue
                if not nt:
                    continue
                candidates.append((cds_key, nt))
    except Exception as exc:
        log.warning("cds_fasta: failed to parse %s: %s", file_path.name, exc)
        return (file_path.name, b"", 0)

    gk = resolve_genome_key(
        file_path=file_path,
        dblink=dblink_entries,
        locus_version=first_locus_version,
        organism=organism,
        strain=strain,
    )
    genome_key = gk.key
    per_file_uids = _WORKER_LOOKUP.get(genome_key) or {}
    if not per_file_uids:
        log.warning(
            "cds_fasta: no id_map rows for genome_key=%s (%s)",
            genome_key, file_path.name,
        )
        return (file_path.name, b"", 0)

    buf = bytearray()
    written = 0
    for cds_key, nt in candidates:
        puid = per_file_uids.get(cds_key)
        if puid is None:
            continue
        buf.extend(f">{puid}\n".encode("ascii"))
        for j in range(0, len(nt), 80):
            buf.extend(nt[j : j + 80].encode("ascii"))
            buf.append(0x0A)  # '\n'
        written += 1

    return (file_path.name, bytes(buf), written)


def write_cds_nt_fasta(
    input_dir: Path,
    id_map_path: Path,
    out_path: Path,
    threads: int = 0,
) -> tuple[int, int]:
    """Write a nucleotide CDS FASTA keyed by ``protein_uid``.

    Headers are ``>{protein_uid}`` (decimal, matching proteins.faa).
    Only CDS that passed protein-level filtering (i.e. appear in
    ``id_map.parquet``) are emitted — pseudogenes and translation-less
    CDS are silently dropped to keep the AA/NT key spaces aligned.

    Parallelized via ProcessPoolExecutor so per-file GenBank reparse +
    ``feat.extract`` runs on multiple cores — matches parse_folder's
    parallelism. ``threads=0`` (default) auto-picks up to 8 workers.

    Returns ``(n_sequences_written, n_missing)``. ``n_missing`` counts
    protein_uids present in id_map whose GenBank CDS could not be
    resolved (usually a malformed location); typically 0.
    """
    id_map = pq.read_table(
        id_map_path,
        columns=["protein_uid", "genome_key", "cds_key"],
    ).to_pydict()

    # Pre-group by genome_key so each worker only sees its own slice.
    # The full lookup is shared once per worker process via the
    # initializer — avoids N-way pickling which would blow up on large
    # runs (5M rows × N files easily >1 GB serialized).
    lookup_by_genome: dict[str, dict[str, int]] = {}
    for puid, gkey, cdskey in zip(
        id_map["protein_uid"], id_map["genome_key"], id_map["cds_key"],
    ):
        lookup_by_genome.setdefault(gkey, {})[cdskey] = int(puid)
    n_total = sum(len(v) for v in lookup_by_genome.values())

    files = find_genbank_files(input_dir)
    if not files:
        raise RuntimeError(f"no GenBank files under {input_dir}")

    n_workers = threads if threads > 0 else min(len(files), 8)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Collect results indexed by filename so the written FASTA is
    # byte-identical across reruns regardless of worker completion order.
    chunks: dict[str, bytes] = {}
    total_written = 0

    n_total_files = len(files)
    # Log every ~10% of files (at least 1) so progress is visible on
    # large assemblies without spamming small runs.
    log_every = max(1, n_total_files // 10)
    completed = 0

    def _log_progress() -> None:
        if completed == n_total_files or completed % log_every == 0:
            log.info(
                "cds_fasta: %d/%d files processed", completed, n_total_files,
            )

    if n_workers <= 1 or len(files) == 1:
        _init_worker(lookup_by_genome)
        for fp in files:
            name, buf, w = _process_one_file(fp)
            chunks[name] = buf
            total_written += w
            completed += 1
            _log_progress()
    else:
        with ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=_init_worker,
            initargs=(lookup_by_genome,),
        ) as pool:
            futures = {pool.submit(_process_one_file, fp): fp for fp in files}
            for fut in as_completed(futures):
                name, buf, w = fut.result()
                chunks[name] = buf
                total_written += w
                completed += 1
                _log_progress()

    with open(out_path, "wb") as out:
        for fp in files:
            out.write(chunks.get(fp.name, b""))

    n_missing = n_total - total_written
    log.info(
        "cds_fasta: wrote %d CDS to %s (%d in id_map, %d missing)",
        total_written, out_path, n_total, n_missing,
    )
    return total_written, max(0, n_missing)


def _resolve_genome_key_for_file(file_path: Path) -> str:
    """Recover the genome_key.key string parse_folder assigns to a file.

    Mirrors parse_folder's scan exactly: accumulates DBLINK across ALL
    records (multi-contig WGS assemblies often scatter BioProject /
    BioSample across contigs), takes the first-seen locus_version /
    organism / strain, then calls ``resolve_genome_key`` with the same
    argument shape. Must not short-circuit — an early break on the
    first record would under-collect DBLINK and diverge the key,
    silently breaking every id_map lookup downstream.
    """
    dblink_entries: list[str] = []
    organism = None
    strain = None
    first_locus_version = None
    try:
        for rec in _iter_records(file_path):
            ann = rec.annotations or {}
            dblink_entries.extend(ann.get("dblink") or [])
            if organism is None:
                organism = ann.get("organism") or None
            if first_locus_version is None:
                first_locus_version = rec.id or rec.name or None
            if strain is None:
                for feat in rec.features:
                    if feat.type == "source":
                        strain_vals = feat.qualifiers.get("strain")
                        if strain_vals:
                            strain = strain_vals[0]
                        break
    except Exception:
        return file_path.stem
    # NOTE: if CLI ever grows manifest support (parse_folder already
    # accepts `manifest=`), thread manifest_name through here too —
    # resolve_genome_key keys on it and the sidecar will silently drift.
    gk = resolve_genome_key(
        file_path=file_path,
        dblink=dblink_entries,
        locus_version=first_locus_version,
        organism=organism,
        strain=strain,
    )
    return gk.key
