"""Shared bidirectional alignment enrichment for clustering engines.

Every clustering backend (MMseqs2, DIAMOND, CD-HIT, usearch12) produces
membership-only ``clusters.parquet`` at first — the five alignment
columns ``pct_identity_fwd``, ``pct_identity_rev``, ``member_coverage``,
``rep_coverage``, ``alignment_length`` start life null. This module
then runs one shared enrichment pass to populate them uniformly across
engines.

### Why a shared enrichment stage

1. **Data uniformity.** Without it, MMseqs2 output has identity and
   coverage while DIAMOND / CD-HIT / usearch12 output has nulls. The
   unified schema is meaningless if only one engine populates it.

2. **Threshold-free alignment.** The previous MMseqs2 pipeline used
   ``easy-search --min-seq-id 0.5 -c 0.8`` for the enrichment pass,
   which filtered out cluster pairs that linclust's k-mer heuristic
   assigned but that full Smith-Waterman couldn't score above the
   threshold. Roughly 1.4% of non-centroid rows silently received null
   alignment metrics. Here we drop all thresholds (``--min-seq-id 0
   -c 0 -e 1000``) so every assigned (member, rep) pair is reported.

3. **Centroid rows computed analytically.** A centroid aligned against
   itself is trivially 100% identity and full coverage. We fill those
   rows directly and only align non-centroid members, halving the
   alignment workload vs the naive all-vs-rep pass.

4. **Representative FASTA derivation.** Some engines (notably DIAMOND
   deepclust) never write a rep_seq.fasta. We derive it from
   ``clusters.parquet`` + the input FASTA so enrichment works
   uniformly regardless of engine.

``--fast`` is the only opt-out: pipeline-wide skip of this stage.
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pyarrow as pa
import pyarrow.csv as pa_csv
import pyarrow.parquet as pq

from ..schemas import CLUSTERS_SCHEMA, validate_schema

log = logging.getLogger(__name__)

PerfectAlignment = (100.0, 1.0, 1.0)  # (pident, member_cov, rep_cov)


@dataclass
class RealignResult:
    """Summary of a realignment run."""

    n_rows: int
    n_centroids: int
    n_aligned: int                 # rows populated from the alignment TSVs
    n_analytical: int              # rows filled by centroid self-rule
    n_missing: int                 # rows still null after the pass
    fwd_tsv: Path | None
    rev_tsv: Path | None


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def populate_alignment_metrics(
    clusters_parquet: Path,
    input_fasta: Path,
    out_dir: Path,
    *,
    threads: int = 0,
    level: str = "protein",
    representatives_fasta: Path | None = None,
    mmseqs_binary: str = "mmseqs",
) -> RealignResult:
    """Enrich ``clusters.parquet`` in place with bidirectional alignment metrics.

    Parameters
    ----------
    clusters_parquet
        Path to the engine's existing ``clusters.parquet``. The
        cluster membership columns (``protein_uid``, ``cluster_id``,
        ``representative_uid``, ``is_centroid``) must already be
        populated; the five alignment columns may be null.
    input_fasta
        Integer-header FASTA fed to the clustering engine. Every
        ``protein_uid`` in ``clusters.parquet`` must appear as a header
        in this file.
    out_dir
        Directory for intermediate alignment TSVs and tmpdirs. Cleaned
        on success.
    threads
        Thread count; 0 = auto-detect.
    level
        ``protein`` or ``nucleotide`` — only affects MMseqs2's
        ``--search-type 3`` switch for nucleotide mode.
    representatives_fasta
        Optional pre-built centroid FASTA. When None, we derive it from
        ``clusters.parquet`` + ``input_fasta``.
    mmseqs_binary
        MMseqs2 binary name on PATH; override for testing.
    """
    if shutil.which(mmseqs_binary) is None:
        raise RuntimeError(
            f"{mmseqs_binary!r} not found on PATH; required for the shared "
            f"realignment stage. Pass --fast to skip alignment enrichment."
        )

    out_dir.mkdir(parents=True, exist_ok=True)

    # ---------- Load existing clusters.parquet ----------
    table = pq.read_table(clusters_parquet)
    validate_schema(table, CLUSTERS_SCHEMA, "clusters.parquet")
    pydict = table.to_pydict()

    protein_uids: list[int] = pydict["protein_uid"]
    rep_uids: list[int] = pydict["representative_uid"]
    is_centroid: list[bool] = pydict["is_centroid"]
    n_rows = len(protein_uids)

    # Member rows are everything that is not its own centroid.
    member_indices: list[int] = [
        i for i, cen in enumerate(is_centroid) if not cen
    ]
    centroid_indices: list[int] = [
        i for i, cen in enumerate(is_centroid) if cen
    ]
    log.info(
        "realign: %d rows (%d centroids, %d non-centroid members)",
        n_rows, len(centroid_indices), len(member_indices),
    )

    # ---------- Build representative FASTA (if needed) ----------
    if representatives_fasta is None or not representatives_fasta.exists():
        rep_uid_set = {int(rep_uids[i]) for i in centroid_indices}
        representatives_fasta = out_dir / "derived_rep_seq.fasta"
        _write_subset_fasta(input_fasta, rep_uid_set, representatives_fasta)
        log.info("realign: derived rep FASTA with %d sequences", len(rep_uid_set))

    # ---------- Build non-centroid member query FASTA ----------
    members_fasta = out_dir / "realign_members.faa"
    member_uid_set = {int(protein_uids[i]) for i in member_indices}
    _write_subset_fasta(input_fasta, member_uid_set, members_fasta)
    log.info("realign: wrote %d member sequences to %s", len(member_uid_set), members_fasta.name)

    # If there are no non-centroid members (everything is a singleton)
    # we can skip alignment entirely — every row is filled analytically.
    if not member_indices:
        return _fill_centroids_only(
            table, pydict, centroid_indices, n_rows, clusters_parquet,
        )

    # ---------- Run forward + reverse MMseqs2 easy-search ----------
    align_tmp = out_dir / "realign_tmp"
    align_tmp.mkdir(parents=True, exist_ok=True)

    fwd_tsv = out_dir / "realign_fwd.tsv"
    rev_tsv = out_dir / "realign_rev.tsv"

    threads = threads if threads > 0 else (os.cpu_count() or 1)
    env = _clean_env()

    fwd_cmd = _build_search_cmd(
        mmseqs_binary=mmseqs_binary,
        query_fasta=members_fasta,
        target_fasta=representatives_fasta,
        out_tsv=fwd_tsv,
        tmpdir=align_tmp / "fwd",
        threads=threads,
        level=level,
        max_seqs=1,
    )
    log.info("realign fwd: %s", " ".join(fwd_cmd))
    _run(fwd_cmd, env, "realign easy-search (forward)")

    rev_cmd = _build_search_cmd(
        mmseqs_binary=mmseqs_binary,
        query_fasta=representatives_fasta,
        target_fasta=members_fasta,
        out_tsv=rev_tsv,
        tmpdir=align_tmp / "rev",
        threads=threads,
        level=level,
        max_seqs=1_000_000,  # generous cap — every member might hit its rep
    )
    log.info("realign rev: %s", " ".join(rev_cmd))
    _run(rev_cmd, env, "realign easy-search (reverse)")

    shutil.rmtree(align_tmp, ignore_errors=True)

    fwd_lookup = _load_alignments(fwd_tsv)
    rev_lookup = _load_alignments(rev_tsv)

    # ---------- Merge metrics into columnar lists ----------
    pid_fwd: list[float | None] = [None] * n_rows
    pid_rev: list[float | None] = [None] * n_rows
    member_cov: list[float | None] = [None] * n_rows
    rep_cov: list[float | None] = [None] * n_rows
    aln_len: list[int | None] = [None] * n_rows

    n_aligned = 0
    n_analytical = 0

    for i in centroid_indices:
        pid_fwd[i] = 100.0
        pid_rev[i] = 100.0
        member_cov[i] = 1.0
        rep_cov[i] = 1.0
        aln_len[i] = None  # centroid self-alignment length is trivially the seq length; leave null
        n_analytical += 1

    for i in member_indices:
        m = int(protein_uids[i])
        r = int(rep_uids[i])

        fwd = fwd_lookup.get((m, r))
        if fwd is not None:
            pid_fwd[i] = fwd[0]
            member_cov[i] = fwd[1]
            rep_cov[i] = fwd[2]
            aln_len[i] = fwd[3]

        rev = rev_lookup.get((r, m))
        if rev is not None:
            pid_rev[i] = rev[0]
            if member_cov[i] is None:
                # reverse qcov is rep's coverage; reverse tcov is member's
                rep_cov[i] = rev[1]
                member_cov[i] = rev[2]
            if aln_len[i] is None:
                aln_len[i] = rev[3]

        if pid_fwd[i] is not None or pid_rev[i] is not None:
            n_aligned += 1

    # ---------- Fallback: direct pairwise alignment for short peptides ----------
    # MMseqs2's k-mer prefilter can't seed on sequences shorter than ~30 aa.
    # These pairs get a full Biopython Gotoh alignment instead — slow per
    # pair but the total count is tiny (usually <10 for a 10-genome run).
    still_missing = [
        i for i in member_indices
        if pid_fwd[i] is None and pid_rev[i] is None
    ]
    n_rescued = 0
    if still_missing:
        log.info(
            "realign: %d pair(s) unaligned by MMseqs2; running Biopython fallback",
            len(still_missing),
        )
        rescued = _biopython_fallback(
            pairs=[(int(protein_uids[i]), int(rep_uids[i])) for i in still_missing],
            input_fasta=input_fasta,
        )
        for i, (m, r) in zip(
            still_missing,
            [(int(protein_uids[i]), int(rep_uids[i])) for i in still_missing],
        ):
            metrics = rescued.get((m, r))
            if metrics is not None:
                pid_fwd[i] = metrics[0]
                pid_rev[i] = metrics[0]  # symmetric global alignment
                member_cov[i] = metrics[1]
                rep_cov[i] = metrics[2]
                aln_len[i] = metrics[3]
                n_aligned += 1
                n_rescued += 1
        log.info("realign: Biopython fallback rescued %d/%d pairs",
                 n_rescued, len(still_missing))

    n_missing = sum(
        1 for i in member_indices
        if pid_fwd[i] is None and pid_rev[i] is None
    )

    # ---------- Write updated clusters.parquet atomically ----------
    result_table = pa.table(
        {
            "protein_uid":        pydict["protein_uid"],
            "genome_uid":         pydict["genome_uid"],
            "cluster_id":         pydict["cluster_id"],
            "representative_uid": pydict["representative_uid"],
            "is_centroid":        pydict["is_centroid"],
            "pct_identity_fwd":   pa.array(pid_fwd, type=pa.float32()),
            "pct_identity_rev":   pa.array(pid_rev, type=pa.float32()),
            "member_coverage":    pa.array(member_cov, type=pa.float32()),
            "rep_coverage":       pa.array(rep_cov, type=pa.float32()),
            "alignment_length":   pa.array(aln_len, type=pa.uint32()),
        },
        schema=CLUSTERS_SCHEMA,
    )
    validate_schema(result_table, CLUSTERS_SCHEMA, "clusters.parquet")
    _atomic_write_parquet(result_table, clusters_parquet)

    return RealignResult(
        n_rows=n_rows,
        n_centroids=len(centroid_indices),
        n_aligned=n_aligned,
        n_analytical=n_analytical,
        n_missing=n_missing,
        fwd_tsv=fwd_tsv,
        rev_tsv=rev_tsv,
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_subset_fasta(
    input_fasta: Path,
    wanted_uids: set[int],
    out_path: Path,
) -> None:
    """Stream ``input_fasta`` and write records whose header int is wanted.

    Headers are the decimal ``protein_uid`` produced by ``fasta.py``.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    keep = False
    written = 0
    with open(input_fasta) as fh, open(out_path, "w") as out:
        for line in fh:
            if line.startswith(">"):
                try:
                    uid = int(line[1:].split()[0])
                except ValueError:
                    keep = False
                    continue
                keep = uid in wanted_uids
                if keep:
                    out.write(line)
                    written += 1
            elif keep:
                out.write(line)
    if written != len(wanted_uids):
        log.warning(
            "subset fasta wrote %d records; expected %d wanted uids",
            written, len(wanted_uids),
        )


def _build_search_cmd(
    *,
    mmseqs_binary: str,
    query_fasta: Path,
    target_fasta: Path,
    out_tsv: Path,
    tmpdir: Path,
    threads: int,
    level: str,
    max_seqs: int,
) -> list[str]:
    tmpdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        mmseqs_binary, "easy-search",
        str(query_fasta), str(target_fasta),
        str(out_tsv), str(tmpdir),
        # Drop every threshold so linclust-assigned pairs that score
        # below 0.5 identity still produce an alignment record.
        "--min-seq-id", "0",
        "-c", "0",
        "--cov-mode", "0",
        "-e", "1e10",                # huge e-value ceiling
        "--format-output", "query,target,pident,qcov,tcov,alnlen",
        "--format-mode", "0",
        "--max-seqs", str(max_seqs),
        "--threads", str(threads),
    ]
    if level == "nucleotide":
        cmd += ["--search-type", "3"]
    return cmd


def _clean_env() -> dict[str, str]:
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    return env


def _run(cmd: Iterable[str], env: dict[str, str], label: str) -> None:
    try:
        subprocess.run(
            list(cmd), check=True, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
        raise RuntimeError(
            f"{label} failed (exit {exc.returncode}):\n{stderr}"
        ) from exc


def _load_alignments(
    path: Path,
) -> dict[tuple[int, int], tuple[float, float, float, int]]:
    """Load an easy-search TSV into a ``(query, target) -> metrics`` dict."""
    if not path.exists():
        return {}
    table = pa_csv.read_csv(
        path,
        read_options=pa_csv.ReadOptions(
            column_names=["query", "target", "pident", "qcov", "tcov", "alnlen"],
        ),
        parse_options=pa_csv.ParseOptions(delimiter="\t"),
        convert_options=pa_csv.ConvertOptions(
            column_types={
                "query":  pa.uint64(),
                "target": pa.uint64(),
                "pident": pa.float32(),
                "qcov":   pa.float32(),
                "tcov":   pa.float32(),
                "alnlen": pa.uint32(),
            },
        ),
    )
    q = table.column("query").to_pylist()
    t = table.column("target").to_pylist()
    pid = table.column("pident").to_pylist()
    qc = table.column("qcov").to_pylist()
    tc = table.column("tcov").to_pylist()
    al = table.column("alnlen").to_pylist()
    lookup: dict[tuple[int, int], tuple[float, float, float, int]] = {}
    for qi, ti, pi, qci, tci, ali in zip(q, t, pid, qc, tc, al):
        lookup[(int(qi), int(ti))] = (pi, qci, tci, ali)
    return lookup


def _atomic_write_parquet(table: pa.Table, path: Path) -> None:
    """Write Parquet to a tempfile and atomically rename into place."""
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    pq.write_table(table, tmp_path, compression="zstd", compression_level=3)
    tmp_path.replace(path)


def _biopython_fallback(
    *,
    pairs: list[tuple[int, int]],
    input_fasta: Path,
) -> dict[tuple[int, int], tuple[float, float, float, int]]:
    """Compute pairwise global alignment for short peptides MMseqs2 can't seed.

    Returns ``{(member_uid, rep_uid): (pident, member_cov, rep_cov, alnlen)}``.
    Uses Biopython's ``PairwiseAligner`` with the BLOSUM62 matrix, affine
    gap penalties, and global mode — equivalent to Needleman-Wunsch.
    """
    try:
        from Bio import Align
        from Bio.Align import substitution_matrices
    except ImportError as exc:
        log.warning("Biopython unavailable for fallback alignment: %s", exc)
        return {}

    # Collect all uids we need from the FASTA in one pass
    wanted: set[int] = set()
    for m, r in pairs:
        wanted.add(m)
        wanted.add(r)

    sequences: dict[int, str] = {}
    header_uid: int | None = None
    with open(input_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                try:
                    uid = int(line[1:].split()[0])
                except ValueError:
                    header_uid = None
                    continue
                header_uid = uid if uid in wanted else None
                if header_uid is not None:
                    sequences[header_uid] = ""
            elif header_uid is not None:
                sequences[header_uid] += line.strip()

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    try:
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    except Exception:  # pragma: no cover - defensive
        aligner.match_score = 2
        aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1

    results: dict[tuple[int, int], tuple[float, float, float, int]] = {}
    for m, r in pairs:
        ms, rs = sequences.get(m), sequences.get(r)
        if not ms or not rs:
            continue
        try:
            alignment = aligner.align(ms, rs)[0]
        except Exception as exc:  # pragma: no cover
            log.warning("pairwise fallback failed for (%d,%d): %s", m, r, exc)
            continue

        fmt = alignment.format().splitlines()
        # PairwiseAligner.format returns three lines: query, match, target
        q_line = fmt[0]
        t_line = fmt[2]
        aligned_cols = min(len(q_line), len(t_line))
        matches = 0
        q_nongap = 0
        t_nongap = 0
        for qc, tc in zip(q_line[:aligned_cols], t_line[:aligned_cols]):
            if qc != "-":
                q_nongap += 1
            if tc != "-":
                t_nongap += 1
            if qc != "-" and tc != "-" and qc == tc:
                matches += 1

        aln_len = aligned_cols
        pident = (matches / aln_len * 100.0) if aln_len else 0.0
        member_cov = (q_nongap / len(ms)) if ms else 0.0
        rep_cov = (t_nongap / len(rs)) if rs else 0.0
        results[(m, r)] = (pident, member_cov, rep_cov, aln_len)

    return results


def _fill_centroids_only(
    table: pa.Table,
    pydict: dict,
    centroid_indices: list[int],
    n_rows: int,
    out_path: Path,
) -> RealignResult:
    """Short-circuit: every row is a centroid (all singleton clusters)."""
    pid_fwd = [100.0 if i in set(centroid_indices) else None for i in range(n_rows)]
    pid_rev = pid_fwd[:]
    member_cov = [1.0 if i in set(centroid_indices) else None for i in range(n_rows)]
    rep_cov = member_cov[:]
    aln_len = [None] * n_rows

    result_table = pa.table(
        {
            "protein_uid":        pydict["protein_uid"],
            "genome_uid":         pydict["genome_uid"],
            "cluster_id":         pydict["cluster_id"],
            "representative_uid": pydict["representative_uid"],
            "is_centroid":        pydict["is_centroid"],
            "pct_identity_fwd":   pa.array(pid_fwd, type=pa.float32()),
            "pct_identity_rev":   pa.array(pid_rev, type=pa.float32()),
            "member_coverage":    pa.array(member_cov, type=pa.float32()),
            "rep_coverage":       pa.array(rep_cov, type=pa.float32()),
            "alignment_length":   pa.array(aln_len, type=pa.uint32()),
        },
        schema=CLUSTERS_SCHEMA,
    )
    _atomic_write_parquet(result_table, out_path)

    return RealignResult(
        n_rows=n_rows,
        n_centroids=len(centroid_indices),
        n_aligned=0,
        n_analytical=len(centroid_indices),
        n_missing=0,
        fwd_tsv=None,
        rev_tsv=None,
    )
