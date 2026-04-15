"""MMseqs2 clustering engine (``easy-linclust`` by default).

MMseqs2 is the speed-first default engine for DNMBcluster. Linclust is
near-linear in input size and scales to billions of sequences;
``easy-cluster`` is reserved for a future ``--sensitive`` mode.

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
from .base import ClusterEngine, ClusterParams, ClusterResult, EngineError, Level

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

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

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

        cmd = self._build_command(input_fasta, work_prefix, tmpdir, params)
        log.info("running mmseqs2: %s", " ".join(cmd))

        env = os.environ.copy()
        # Prevent OMP/BLAS oversubscription — SPEED.md Section 7.
        env.setdefault("OMP_NUM_THREADS", "1")
        env.setdefault("OPENBLAS_NUM_THREADS", "1")
        env.setdefault("MKL_NUM_THREADS", "1")

        try:
            subprocess.run(
                cmd,
                check=True,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except subprocess.CalledProcessError as exc:
            stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
            raise EngineError(f"mmseqs2 failed (exit {exc.returncode}):\n{stderr}") from exc

        cluster_tsv = work_prefix.with_name(work_prefix.name + "_cluster.tsv")
        if not cluster_tsv.exists():
            raise EngineError(
                f"mmseqs2 did not produce expected output {cluster_tsv}"
            )

        clusters_parquet = out_dir / "clusters.parquet"
        n_clusters = parse_cluster_tsv_to_parquet(
            cluster_tsv=cluster_tsv,
            out_parquet=clusters_parquet,
            input_fasta=input_fasta,
        )

        rep_fasta = work_prefix.with_name(work_prefix.name + "_rep_seq.fasta")

        # Wipe MMseqs scratch to avoid cache bloat — SPEED.md Section 8.
        shutil.rmtree(tmpdir, ignore_errors=True)

        n_input = len(parse_fasta_headers(input_fasta))
        return ClusterResult(
            clusters_parquet=clusters_parquet,
            representatives_fasta=rep_fasta if rep_fasta.exists() else None,
            engine=self.name,
            params=params,
            n_input_sequences=n_input,
            n_clusters=n_clusters,
            tmpdir=None,
        )

    # ------------------------------------------------------------------
    # Command construction
    # ------------------------------------------------------------------

    def _build_command(
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
            # MMseqs2 understands suffixes like "8G"; pick 80% to leave
            # headroom for Polars/NumPy in the same process.
            limit = max(int(params.max_ram_gb * 0.8), 1)
            cmd += ["--split-memory-limit", f"{limit}G"]

        if params.level == "nucleotide":
            # --search-type 3 forces nucleotide DNA mode end-to-end.
            cmd += ["--search-type", "3"]

        return cmd


# ---------------------------------------------------------------------------
# TSV → unified clusters.parquet
# ---------------------------------------------------------------------------


def parse_cluster_tsv_to_parquet(
    *,
    cluster_tsv: Path,
    out_parquet: Path,
    input_fasta: Path | None = None,
) -> int:
    """Convert MMseqs2 ``_cluster.tsv`` into the unified DNMB schema.

    The MMseqs2 TSV has two columns: representative ID and member ID,
    both being the integer ``protein_uid`` we wrote in the FASTA. We
    need to:

    1. Parse as uint64.
    2. Assign a dense 0-based ``cluster_id`` per unique representative.
    3. Derive ``genome_uid`` from the upper 16 bits of ``protein_uid``.
    4. Flag the centroid row (``protein_uid == representative_uid``).

    MMseqs2 does not expose per-member identity in ``easy-linclust`` TSV,
    so ``pct_identity`` is left null.

    Returns the number of distinct clusters.
    """
    # Load as two uint64 columns. polars would be faster for >10M rows
    # but Arrow CSV reader is already vectorized and avoids adding a
    # polars dependency at this layer.
    read_opts = pa_csv.ReadOptions(
        column_names=["representative_uid", "protein_uid"],
    )
    parse_opts = pa_csv.ParseOptions(delimiter="\t")
    convert_opts = pa_csv.ConvertOptions(
        column_types={
            "representative_uid": pa.uint64(),
            "protein_uid": pa.uint64(),
        },
    )
    table = pa_csv.read_csv(
        cluster_tsv,
        read_options=read_opts,
        parse_options=parse_opts,
        convert_options=convert_opts,
    )

    representatives = table.column("representative_uid")
    protein_uids = table.column("protein_uid")

    # Assign dense cluster_id: one per unique representative, in the
    # order they first appear in the TSV.
    dense_map: dict[int, int] = {}
    cluster_id_list: list[int] = []
    for rep in representatives.to_pylist():
        if rep not in dense_map:
            dense_map[rep] = len(dense_map)
        cluster_id_list.append(dense_map[rep])

    genome_uid_list = [
        (uid >> 48) & 0xFFFF for uid in protein_uids.to_pylist()
    ]
    is_centroid_list = [
        int(p) == int(r)
        for p, r in zip(protein_uids.to_pylist(), representatives.to_pylist())
    ]

    result = pa.table(
        {
            "protein_uid": protein_uids,
            "genome_uid": pa.array(genome_uid_list, type=pa.uint16()),
            "cluster_id": pa.array(cluster_id_list, type=pa.uint32()),
            "representative_uid": representatives,
            "is_centroid": pa.array(is_centroid_list, type=pa.bool_()),
            "pct_identity": pa.array(
                [None] * len(cluster_id_list),
                type=pa.float32(),
            ),
        },
        schema=CLUSTERS_SCHEMA,
    )
    validate_schema(result, CLUSTERS_SCHEMA, "clusters.parquet")

    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(
        result,
        out_parquet,
        compression="zstd",
        compression_level=3,
    )
    return len(dense_map)
