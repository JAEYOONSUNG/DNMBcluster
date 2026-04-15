"""MMseqs2 clustering engine.

Three-stage workflow by default:

1. ``mmseqs easy-linclust`` — fast linear-time clustering.
2. Forward ``mmseqs easy-search`` (member as query, seed as target) —
   populates ``pct_identity_fwd``, ``member_coverage``, ``rep_coverage``,
   ``alignment_length``.
3. Reverse ``mmseqs easy-search`` (seed as query, member as target) —
   populates ``pct_identity_rev``.

Two independent alignments give bidirectional identity. Smith-Waterman
is symmetric for untrimmed global alignment but the heuristic pipelines
(k-mer seeding, gapped extension, pre-filter thresholds) can introduce
direction-dependent identity, and reporting both is the only way to
catch that for length-disparate sequences.

Users who want pure speed pass ``--fast`` (``with_alignment=False``) to
skip both alignment passes; those fields are then null.

SPEED.md Sections 4, 12, and the length / bidirectional identity /
coverage requirement added 2026-04-15.
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
from pathlib import Path

import pyarrow as pa
import pyarrow.csv as pa_csv
import pyarrow.parquet as pq

from ..fasta import parse_fasta_headers
from ..schemas import CLUSTERS_SCHEMA, validate_schema
from .base import ClusterEngine, ClusterParams, ClusterResult, EngineError

log = logging.getLogger(__name__)


class MMseqs2Engine(ClusterEngine):
    name = "mmseqs2"
    supported_levels = ("protein", "nucleotide")

    def __init__(self, binary: str = "mmseqs") -> None:
        self.binary = binary

    def check_available(self) -> None:
        if shutil.which(self.binary) is None:
            raise EngineError(
                f"mmseqs2 binary {self.binary!r} not found on PATH. "
                f"Install via `conda install -c bioconda mmseqs2` or run inside "
                f"the DNMBcluster Docker image."
            )

    def cluster(
        self,
        input_fasta: Path,
        out_dir: Path,
        params: ClusterParams,
    ) -> ClusterResult:
        self.check_available()
        if not self.supports(params.level):
            raise EngineError(f"mmseqs2 does not support level={params.level!r}")

        out_dir.mkdir(parents=True, exist_ok=True)
        work_prefix = out_dir / "mmseqs_out"
        tmpdir = out_dir / "mmseqs_tmp"
        tmpdir.mkdir(parents=True, exist_ok=True)

        env = _mmseqs_env()
        threads = params.threads if params.threads > 0 else (os.cpu_count() or 1)

        # ---------- Stage 1: fast clustering ----------
        cluster_cmd = self._build_cluster_command(
            input_fasta, work_prefix, tmpdir, params, threads,
        )
        log.info("running mmseqs2 linclust: %s", " ".join(cluster_cmd))
        _run(cluster_cmd, env, "mmseqs2 easy-linclust")

        cluster_tsv = work_prefix.with_name(work_prefix.name + "_cluster.tsv")
        rep_fasta = work_prefix.with_name(work_prefix.name + "_rep_seq.fasta")
        if not cluster_tsv.exists():
            raise EngineError(f"mmseqs2 did not produce {cluster_tsv}")

        # ---------- Stage 2 + 3: bidirectional alignment pass (optional) ----------
        alignments_fwd: Path | None = None
        alignments_rev: Path | None = None
        if params.with_alignment and rep_fasta.exists():
            align_tmpdir = out_dir / "mmseqs_align_tmp"
            align_tmpdir.mkdir(parents=True, exist_ok=True)

            alignments_fwd = out_dir / "mmseqs_alignments_fwd.tsv"
            fwd_cmd = self._build_search_command(
                query_fasta=input_fasta,
                target_fasta=rep_fasta,
                out_tsv=alignments_fwd,
                tmpdir=align_tmpdir / "fwd",
                params=params,
                threads=threads,
                max_seqs=1,
            )
            log.info("mmseqs2 fwd easy-search: %s", " ".join(fwd_cmd))
            try:
                _run(fwd_cmd, env, "mmseqs2 easy-search (forward)")
            except EngineError as exc:
                log.warning("fwd alignment pass failed: %s", exc)
                alignments_fwd = None

            # Reverse: every rep vs every input protein. Each rep gets
            # many hits (all members of its cluster), so max_seqs must
            # be large enough to cover the largest cluster. We pass a
            # generous default — MMseqs2 will simply return fewer when
            # clusters are smaller.
            alignments_rev = out_dir / "mmseqs_alignments_rev.tsv"
            rev_cmd = self._build_search_command(
                query_fasta=rep_fasta,
                target_fasta=input_fasta,
                out_tsv=alignments_rev,
                tmpdir=align_tmpdir / "rev",
                params=params,
                threads=threads,
                max_seqs=100_000,
            )
            log.info("mmseqs2 rev easy-search: %s", " ".join(rev_cmd))
            try:
                _run(rev_cmd, env, "mmseqs2 easy-search (reverse)")
            except EngineError as exc:
                log.warning("rev alignment pass failed: %s", exc)
                alignments_rev = None

            shutil.rmtree(align_tmpdir, ignore_errors=True)

        # ---------- Convert to unified clusters.parquet ----------
        clusters_parquet = out_dir / "clusters.parquet"
        n_clusters = parse_cluster_tsv_to_parquet(
            cluster_tsv=cluster_tsv,
            alignments_fwd=alignments_fwd,
            alignments_rev=alignments_rev,
            out_parquet=clusters_parquet,
        )

        shutil.rmtree(tmpdir, ignore_errors=True)

        n_input = len(parse_fasta_headers(input_fasta))
        return ClusterResult(
            clusters_parquet=clusters_parquet,
            representatives_fasta=rep_fasta if rep_fasta.exists() else None,
            engine=self.name,
            params=params,
            n_input_sequences=n_input,
            n_clusters=n_clusters,
        )

    # ------------------------------------------------------------------
    # Command construction
    # ------------------------------------------------------------------

    def _build_cluster_command(
        self,
        input_fasta: Path,
        work_prefix: Path,
        tmpdir: Path,
        params: ClusterParams,
        threads: int,
    ) -> list[str]:
        cmd: list[str] = [
            self.binary,
            "easy-linclust",
            str(input_fasta),
            str(work_prefix),
            str(tmpdir),
            "--min-seq-id", f"{params.identity}",
            "-c", f"{params.coverage}",
            "--cov-mode", "0",
            "--threads", str(threads),
        ]
        if params.max_ram_gb and params.max_ram_gb > 0:
            limit = max(int(params.max_ram_gb * 0.8), 1)
            cmd += ["--split-memory-limit", f"{limit}G"]
        if params.level == "nucleotide":
            cmd += ["--search-type", "3"]
        return cmd

    def _build_search_command(
        self,
        *,
        query_fasta: Path,
        target_fasta: Path,
        out_tsv: Path,
        tmpdir: Path,
        params: ClusterParams,
        threads: int,
        max_seqs: int,
    ) -> list[str]:
        tmpdir.mkdir(parents=True, exist_ok=True)
        cmd: list[str] = [
            self.binary,
            "easy-search",
            str(query_fasta),
            str(target_fasta),
            str(out_tsv),
            str(tmpdir),
            "--min-seq-id", f"{params.identity}",
            "-c", f"{params.coverage}",
            "--cov-mode", "0",
            "--format-output", "query,target,pident,qcov,tcov,alnlen",
            "--format-mode", "0",
            "--max-seqs", str(max_seqs),
            "--threads", str(threads),
        ]
        if params.level == "nucleotide":
            cmd += ["--search-type", "3"]
        return cmd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _mmseqs_env() -> dict[str, str]:
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    return env


def _run(cmd: list[str], env: dict[str, str], label: str) -> None:
    try:
        subprocess.run(
            cmd, check=True, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
        raise EngineError(f"{label} failed (exit {exc.returncode}):\n{stderr}") from exc


# ---------------------------------------------------------------------------
# TSV → unified clusters.parquet
# ---------------------------------------------------------------------------


def _load_alignments(path: Path) -> dict[tuple[int, int], tuple[float, float, float, int]]:
    """Load an easy-search TSV and return a ``(query, target) -> metrics`` dict.

    Metrics tuple is ``(pident, qcov, tcov, alnlen)``.
    """
    aln_table = pa_csv.read_csv(
        path,
        read_options=pa_csv.ReadOptions(
            column_names=["query", "target", "pident", "qcov", "tcov", "alnlen"],
        ),
        parse_options=pa_csv.ParseOptions(delimiter="\t"),
        convert_options=pa_csv.ConvertOptions(
            column_types={
                "query": pa.uint64(),
                "target": pa.uint64(),
                "pident": pa.float32(),
                "qcov": pa.float32(),
                "tcov": pa.float32(),
                "alnlen": pa.uint32(),
            },
        ),
    )
    q_col = aln_table.column("query").to_pylist()
    t_col = aln_table.column("target").to_pylist()
    pid_col = aln_table.column("pident").to_pylist()
    qcov_col = aln_table.column("qcov").to_pylist()
    tcov_col = aln_table.column("tcov").to_pylist()
    alnlen_col = aln_table.column("alnlen").to_pylist()
    lookup: dict[tuple[int, int], tuple[float, float, float, int]] = {}
    for q, t, pid, qc, tc, al in zip(
        q_col, t_col, pid_col, qcov_col, tcov_col, alnlen_col,
    ):
        lookup[(int(q), int(t))] = (pid, qc, tc, al)
    return lookup


def parse_cluster_tsv_to_parquet(
    *,
    cluster_tsv: Path,
    out_parquet: Path,
    alignments_fwd: Path | None = None,
    alignments_rev: Path | None = None,
    input_fasta: Path | None = None,  # unused; kept for API compat
) -> int:
    """Convert MMseqs2 output into the unified DNMB schema.

    - ``cluster_tsv``: ``_cluster.tsv`` from easy-linclust (two-column).
    - ``alignments_fwd``: optional easy-search output with member as
      query and seed as target. Used to populate ``pct_identity_fwd``
      and the two coverage columns.
    - ``alignments_rev``: optional easy-search output with seed as
      query and member as target. Used to populate ``pct_identity_rev``.

    Returns the number of distinct clusters.
    """
    cluster_table = pa_csv.read_csv(
        cluster_tsv,
        read_options=pa_csv.ReadOptions(
            column_names=["representative_uid", "protein_uid"],
        ),
        parse_options=pa_csv.ParseOptions(delimiter="\t"),
        convert_options=pa_csv.ConvertOptions(
            column_types={
                "representative_uid": pa.uint64(),
                "protein_uid": pa.uint64(),
            },
        ),
    )
    representatives = cluster_table.column("representative_uid").to_pylist()
    members = cluster_table.column("protein_uid").to_pylist()
    n_rows = len(representatives)

    dense_map: dict[int, int] = {}
    cluster_id_list: list[int] = []
    for rep in representatives:
        if rep not in dense_map:
            dense_map[rep] = len(dense_map)
        cluster_id_list.append(dense_map[rep])

    genome_uid_list = [(uid >> 48) & 0xFFFF for uid in members]
    is_centroid_list = [int(m) == int(r) for m, r in zip(members, representatives)]

    pct_id_fwd: list[float | None] = [None] * n_rows
    pct_id_rev: list[float | None] = [None] * n_rows
    member_cov: list[float | None] = [None] * n_rows
    rep_cov: list[float | None] = [None] * n_rows
    aln_len: list[int | None] = [None] * n_rows

    fwd_lookup = _load_alignments(alignments_fwd) if alignments_fwd and alignments_fwd.exists() else {}
    rev_lookup = _load_alignments(alignments_rev) if alignments_rev and alignments_rev.exists() else {}

    for i, (member, rep) in enumerate(zip(members, representatives)):
        m_int, r_int = int(member), int(rep)

        # Forward direction: query = member, target = rep
        fwd = fwd_lookup.get((m_int, r_int))
        if fwd is not None:
            pct_id_fwd[i] = fwd[0]
            member_cov[i] = fwd[1]
            rep_cov[i] = fwd[2]
            aln_len[i] = fwd[3]

        # Reverse direction: query = rep, target = member
        rev = rev_lookup.get((r_int, m_int))
        if rev is not None:
            pct_id_rev[i] = rev[0]
            if member_cov[i] is None:
                # Reverse qcov is rep's coverage; reverse tcov is member's.
                rep_cov[i] = rev[1]
                member_cov[i] = rev[2]
                aln_len[i] = rev[3] if aln_len[i] is None else aln_len[i]

        # Centroid self-alignment: perfect by definition.
        if m_int == r_int:
            if pct_id_fwd[i] is None:
                pct_id_fwd[i] = 100.0
            if pct_id_rev[i] is None:
                pct_id_rev[i] = 100.0
            if member_cov[i] is None:
                member_cov[i] = 1.0
            if rep_cov[i] is None:
                rep_cov[i] = 1.0

    result = pa.table(
        {
            "protein_uid": pa.array(members, type=pa.uint64()),
            "genome_uid": pa.array(genome_uid_list, type=pa.uint16()),
            "cluster_id": pa.array(cluster_id_list, type=pa.uint32()),
            "representative_uid": pa.array(representatives, type=pa.uint64()),
            "is_centroid": pa.array(is_centroid_list, type=pa.bool_()),
            "pct_identity_fwd": pa.array(pct_id_fwd, type=pa.float32()),
            "pct_identity_rev": pa.array(pct_id_rev, type=pa.float32()),
            "member_coverage": pa.array(member_cov, type=pa.float32()),
            "rep_coverage": pa.array(rep_cov, type=pa.float32()),
            "alignment_length": pa.array(aln_len, type=pa.uint32()),
        },
        schema=CLUSTERS_SCHEMA,
    )
    validate_schema(result, CLUSTERS_SCHEMA, "clusters.parquet")

    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(result, out_parquet, compression="zstd", compression_level=3)
    return len(dense_map)
