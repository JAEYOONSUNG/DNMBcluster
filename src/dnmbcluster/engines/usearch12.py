"""usearch12 clustering engine.

Uses ``cluster_fast`` (or ``cluster_mt`` when ≥8 threads are available
— SPEED.md Section 4). Alphabet is auto-detected by usearch12 from the
FASTA content, so the same command works for both protein and
nucleotide modes.

The UC file format (``-uc``) is a 10-column TSV; see SPEED.md Section 1
for the column description. We only keep rows of type ``S`` (seed /
centroid) and ``H`` (hit / member) — ``C`` and ``N`` rows are
redundant or noise.
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from ..fasta import parse_fasta_headers
from ..schemas import CLUSTERS_SCHEMA, validate_schema
from .base import ClusterEngine, ClusterParams, ClusterResult, EngineError

log = logging.getLogger(__name__)


class Usearch12Engine(ClusterEngine):
    name = "usearch12"
    supported_levels = ("protein", "nucleotide")

    def __init__(self, binary: str = "usearch12") -> None:
        self.binary = binary

    def check_available(self) -> None:
        if shutil.which(self.binary) is None:
            raise EngineError(
                f"{self.binary!r} not found on PATH. In the DNMBcluster "
                f"Docker image this is compiled from rcedgar/usearch12 at "
                f"build time; on bare metal, build from "
                f"https://github.com/rcedgar/usearch12 and put the binary "
                f"on PATH."
            )

    def cluster(
        self,
        input_fasta: Path,
        out_dir: Path,
        params: ClusterParams,
    ) -> ClusterResult:
        self.check_available()
        out_dir.mkdir(parents=True, exist_ok=True)

        uc_path = out_dir / "usearch12_clusters.uc"
        rep_path = out_dir / "usearch12_rep_seq.fasta"

        threads = params.threads if params.threads > 0 else (os.cpu_count() or 1)
        subcommand = "-cluster_mt" if threads >= 8 else "-cluster_fast"

        cmd: list[str] = [
            self.binary,
            subcommand, str(input_fasta),
            "-id", f"{params.identity}",
            "-uc", str(uc_path),
            "-centroids", str(rep_path),
        ]
        if subcommand == "-cluster_fast":
            # -sort length is beneficial for cluster_fast; cluster_mt
            # picks its own ordering.
            cmd += ["-sort", "length"]
        if threads > 1 and subcommand == "-cluster_mt":
            cmd += ["-threads", str(threads)]

        log.info("running usearch12: %s", " ".join(cmd))
        try:
            subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except subprocess.CalledProcessError as exc:
            stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
            raise EngineError(f"usearch12 failed (exit {exc.returncode}):\n{stderr}") from exc

        if not uc_path.exists():
            raise EngineError(f"usearch12 did not produce {uc_path}")

        clusters_parquet = out_dir / "clusters.parquet"
        n_clusters = parse_uc(
            uc_path=uc_path,
            out_parquet=clusters_parquet,
        )

        n_input = len(parse_fasta_headers(input_fasta))
        return ClusterResult(
            clusters_parquet=clusters_parquet,
            representatives_fasta=rep_path if rep_path.exists() else None,
            engine=self.name,
            params=params,
            n_input_sequences=n_input,
            n_clusters=n_clusters,
        )


def parse_uc(*, uc_path: Path, out_parquet: Path) -> int:
    """Parse usearch12 ``.uc`` output into the unified schema.

    Columns of interest (0-indexed):

    - 0: record type (S=centroid, H=hit, C=cluster summary, N=no hit)
    - 1: 0-based cluster id
    - 3: percent identity (H rows only)
    - 8: query label (our integer protein_uid)
    - 9: target label (centroid protein_uid, '*' for S rows)
    """
    protein_uids: list[int] = []
    cluster_ids: list[int] = []
    representative_uids: list[int] = []
    pct_identities: list[float | None] = []
    is_centroid_flags: list[bool] = []

    # First pass: collect S rows so we can later resolve representatives
    # for H rows in a single allocation.
    centroid_for_cluster: dict[int, int] = {}

    with open(uc_path) as fh:
        # Two-pass over the file: first find centroids, then emit rows.
        # A .uc file is usually small enough that this is fine; for very
        # large runs we could stream once and defer centroid resolution,
        # but it's not the hot path.
        lines = fh.readlines()

    for line in lines:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 10:
            continue
        rec_type = parts[0]
        if rec_type != "S":
            continue
        try:
            cluster_id = int(parts[1])
            uid = int(parts[8])
        except ValueError:
            continue
        centroid_for_cluster[cluster_id] = uid

    for line in lines:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 10:
            continue
        rec_type = parts[0]
        if rec_type not in {"S", "H"}:
            continue
        try:
            cluster_id = int(parts[1])
            uid = int(parts[8])
        except ValueError:
            continue

        rep_uid = centroid_for_cluster.get(cluster_id)
        if rep_uid is None:
            # Orphan H row — usearch12 sometimes emits H before S in
            # unsorted output. Use the query as its own rep.
            rep_uid = uid

        if rec_type == "S":
            pct: float | None = None
            is_centroid = True
        else:
            try:
                pct = float(parts[3])
            except ValueError:
                pct = None
            is_centroid = False

        protein_uids.append(uid)
        cluster_ids.append(cluster_id)
        representative_uids.append(rep_uid)
        pct_identities.append(pct)
        is_centroid_flags.append(is_centroid)

    # usearch12 cluster_id is already 0-based dense; keep as-is but
    # re-rank defensively in case the input had gaps.
    unique_clusters = sorted(set(cluster_ids))
    rank = {c: i for i, c in enumerate(unique_clusters)}
    cluster_ids = [rank[c] for c in cluster_ids]

    genome_uids = [(u >> 48) & 0xFFFF for u in protein_uids]
    n_rows = len(protein_uids)
    null_f32: list[float | None] = [None] * n_rows
    null_u32: list[int | None] = [None] * n_rows

    result = pa.table(
        {
            "protein_uid": pa.array(protein_uids, type=pa.uint64()),
            "genome_uid": pa.array(genome_uids, type=pa.uint16()),
            "cluster_id": pa.array(cluster_ids, type=pa.uint32()),
            "representative_uid": pa.array(representative_uids, type=pa.uint64()),
            "is_centroid": pa.array(is_centroid_flags, type=pa.bool_()),
            # .uc column 3 is a single pct_identity; stored as forward.
            "pct_identity_fwd": pa.array(pct_identities, type=pa.float32()),
            "pct_identity_rev": pa.array(null_f32, type=pa.float32()),
            "member_coverage": pa.array(null_f32, type=pa.float32()),
            "rep_coverage": pa.array(null_f32, type=pa.float32()),
            "alignment_length": pa.array(null_u32, type=pa.uint32()),
        },
        schema=CLUSTERS_SCHEMA,
    )
    validate_schema(result, CLUSTERS_SCHEMA, "clusters.parquet")

    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(result, out_parquet, compression="zstd", compression_level=3)
    return len(unique_clusters)
