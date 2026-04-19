"""DIAMOND deepclust clustering engine (protein-only).

DIAMOND is the fastest protein clustering tool at DNMBcluster's target
scales — up to 82× CD-HIT in published benchmarks. It is **protein
only**; nucleotide level raises EngineError and the CLI pre-validates
so we never reach engine code with the wrong level.

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


class DiamondEngine(ClusterEngine):
    name = "diamond"
    supported_levels = ("protein",)

    def __init__(self, binary: str = "diamond") -> None:
        self.binary = binary

    def check_available(self) -> None:
        if shutil.which(self.binary) is None:
            raise EngineError(
                f"diamond binary {self.binary!r} not found on PATH. "
                f"Install via `conda install -c bioconda diamond>=2.1` or "
                f"run inside the DNMBcluster Docker image."
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
            raise EngineError(
                f"diamond is protein-only. Use --tool mmseqs2 or usearch12 "
                f"for nucleotide clustering."
            )

        raw_dir.mkdir(parents=True, exist_ok=True)
        processed_dir.mkdir(parents=True, exist_ok=True)
        cluster_tsv = raw_dir / "diamond_cluster.tsv"
        rep_fasta = raw_dir / "diamond_rep_seq.fasta"
        tmpdir = raw_dir / "diamond_tmp"
        tmpdir.mkdir(parents=True, exist_ok=True)

        threads = params.threads if params.threads > 0 else (os.cpu_count() or 1)
        approx_id = int(params.identity * 100)
        member_cover = int(params.coverage * 100)

        cmd: list[str] = [
            self.binary,
            "deepclust",
            "-d", str(input_fasta),
            "-o", str(cluster_tsv),
            "--approx-id", str(approx_id),
            "--member-cover", str(member_cover),
            "-p", str(threads),
            "--header",  # emit a header line so our parser can skip it
            "--tmpdir", str(tmpdir),
        ]
        if params.max_ram_gb and params.max_ram_gb > 0:
            limit = max(int(params.max_ram_gb * 0.75), 1)
            cmd += ["-M", f"{limit}G"]

        log.info("running diamond: %s", " ".join(cmd))
        env = os.environ.copy()
        env.setdefault("OMP_NUM_THREADS", "1")
        env.setdefault("OPENBLAS_NUM_THREADS", "1")

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
            raise EngineError(f"diamond failed (exit {exc.returncode}):\n{stderr}") from exc

        if not cluster_tsv.exists():
            raise EngineError(f"diamond did not produce expected output {cluster_tsv}")

        clusters_parquet = processed_dir / "clusters.parquet"
        n_clusters = parse_diamond_cluster_tsv(
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


def parse_diamond_cluster_tsv(
    *,
    cluster_tsv: Path,
    out_parquet: Path,
) -> int:
    """Parse DIAMOND ``deepclust`` output into the unified schema.

    DIAMOND writes a two-column TSV: ``centroid<TAB>member``. With
    ``--header``, the first line is ``#centroid<TAB>member``, which the
    Arrow CSV reader treats as a normal header.
    """
    read_opts = pa_csv.ReadOptions(skip_rows_after_names=0)
    parse_opts = pa_csv.ParseOptions(delimiter="\t")
    convert_opts = pa_csv.ConvertOptions(
        column_types={
            "#centroid": pa.uint64(),
            "centroid": pa.uint64(),
            "member": pa.uint64(),
        },
        include_columns=None,
    )
    table = pa_csv.read_csv(
        cluster_tsv,
        read_options=read_opts,
        parse_options=parse_opts,
        convert_options=convert_opts,
    )

    # DIAMOND may or may not prefix the header with '#'. Normalize.
    col_names = table.column_names
    if "#centroid" in col_names:
        rep_col = table.column("#centroid")
    elif "centroid" in col_names:
        rep_col = table.column("centroid")
    else:
        rep_col = table.column(0)
    member_col = table.column("member") if "member" in col_names else table.column(1)

    dense_map: dict[int, int] = {}
    cluster_id_list: list[int] = []
    for rep in rep_col.to_pylist():
        if rep not in dense_map:
            dense_map[rep] = len(dense_map)
        cluster_id_list.append(dense_map[rep])

    member_uids = member_col.to_pylist()
    rep_uids = rep_col.to_pylist()
    genome_uid_list = [(u >> 48) & 0xFFFF for u in member_uids]
    is_centroid_list = [int(m) == int(r) for m, r in zip(member_uids, rep_uids)]

    n_rows = len(cluster_id_list)
    null_f32: list[float | None] = [None] * n_rows
    null_u32: list[int | None] = [None] * n_rows
    result = pa.table(
        {
            "protein_uid": pa.array(member_uids, type=pa.uint64()),
            "genome_uid": pa.array(genome_uid_list, type=pa.uint16()),
            "cluster_id": pa.array(cluster_id_list, type=pa.uint32()),
            "representative_uid": pa.array(rep_uids, type=pa.uint64()),
            "is_centroid": pa.array(is_centroid_list, type=pa.bool_()),
            "pct_identity_fwd": pa.array(null_f32, type=pa.float32()),
            "pct_identity_rev": pa.array(null_f32, type=pa.float32()),
            "member_coverage": pa.array(null_f32, type=pa.float32()),
            "rep_coverage": pa.array(null_f32, type=pa.float32()),
            "alignment_length": pa.array(null_u32, type=pa.uint32()),
        },
        schema=CLUSTERS_SCHEMA,
    )
    validate_schema(result, CLUSTERS_SCHEMA, "clusters.parquet")

    from ..io_utils import atomic_write_table
    atomic_write_table(result, out_parquet)
    return len(dense_map)
