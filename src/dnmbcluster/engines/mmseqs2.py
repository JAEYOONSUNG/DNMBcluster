"""MMseqs2 clustering engine.

Produces cluster membership via ``mmseqs easy-linclust``. The shared
``realign.py`` stage populates the alignment metric columns uniformly
across every engine, so this module is now membership-only — it does
not run the bidirectional ``easy-search`` pass itself.

SPEED.md Sections 4 and 12.
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
        raw_dir: Path,
        processed_dir: Path,
        params: ClusterParams,
    ) -> ClusterResult:
        self.check_available()
        if not self.supports(params.level):
            raise EngineError(f"mmseqs2 does not support level={params.level!r}")

        raw_dir.mkdir(parents=True, exist_ok=True)
        processed_dir.mkdir(parents=True, exist_ok=True)
        work_prefix = raw_dir / "mmseqs_out"
        tmpdir = raw_dir / "mmseqs_tmp"
        tmpdir.mkdir(parents=True, exist_ok=True)

        cmd = self._build_cluster_command(
            input_fasta, work_prefix, tmpdir, params,
        )
        log.info("running mmseqs2 linclust: %s", " ".join(cmd))

        env = os.environ.copy()
        env.setdefault("OMP_NUM_THREADS", "1")
        env.setdefault("OPENBLAS_NUM_THREADS", "1")
        env.setdefault("MKL_NUM_THREADS", "1")

        try:
            subprocess.run(
                cmd, check=True, env=env,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
        except subprocess.CalledProcessError as exc:
            stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
            raise EngineError(
                f"mmseqs2 easy-linclust failed (exit {exc.returncode}):\n{stderr}"
            ) from exc

        cluster_tsv = work_prefix.with_name(work_prefix.name + "_cluster.tsv")
        rep_fasta = work_prefix.with_name(work_prefix.name + "_rep_seq.fasta")
        if not cluster_tsv.exists():
            raise EngineError(f"mmseqs2 did not produce {cluster_tsv}")

        clusters_parquet = processed_dir / "clusters.parquet"
        n_clusters = parse_cluster_tsv_to_parquet(
            cluster_tsv=cluster_tsv,
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

    def _build_cluster_command(
        self,
        input_fasta: Path,
        work_prefix: Path,
        tmpdir: Path,
        params: ClusterParams,
    ) -> list[str]:
        threads = params.threads if params.threads > 0 else (os.cpu_count() or 1)
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


# ---------------------------------------------------------------------------
# TSV → unified clusters.parquet (membership only)
# ---------------------------------------------------------------------------


def parse_cluster_tsv_to_parquet(
    *,
    cluster_tsv: Path,
    out_parquet: Path,
) -> int:
    """Convert MMseqs2 ``_cluster.tsv`` into membership-only ``clusters.parquet``.

    Alignment metric columns (``pct_identity_fwd``, ``pct_identity_rev``,
    ``member_coverage``, ``rep_coverage``, ``alignment_length``) are
    written as null. The shared ``realign.py`` stage populates them
    after this module returns.

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

    null_f32: list[float | None] = [None] * n_rows
    null_u32: list[int | None] = [None] * n_rows

    result = pa.table(
        {
            "protein_uid":        pa.array(members, type=pa.uint64()),
            "genome_uid":         pa.array(genome_uid_list, type=pa.uint16()),
            "cluster_id":         pa.array(cluster_id_list, type=pa.uint32()),
            "representative_uid": pa.array(representatives, type=pa.uint64()),
            "is_centroid":        pa.array(is_centroid_list, type=pa.bool_()),
            "pct_identity_fwd":   pa.array(null_f32, type=pa.float32()),
            "pct_identity_rev":   pa.array(null_f32, type=pa.float32()),
            "member_coverage":    pa.array(null_f32, type=pa.float32()),
            "rep_coverage":       pa.array(null_f32, type=pa.float32()),
            "alignment_length":   pa.array(null_u32, type=pa.uint32()),
        },
        schema=CLUSTERS_SCHEMA,
    )
    validate_schema(result, CLUSTERS_SCHEMA, "clusters.parquet")

    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(result, out_parquet, compression="zstd", compression_level=3)
    return len(dense_map)
