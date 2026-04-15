"""CD-HIT clustering engine.

Dispatches to ``cd-hit`` for protein and ``cd-hit-est`` for nucleotide
(they are separate binaries with slightly different flags and
word-size rules). CD-HIT's ``.clstr`` output needs regex parsing — the
format is text-oriented, not tabular.

SPEED.md Sections 4 and 12.
"""
from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from ..fasta import parse_fasta_headers
from ..schemas import CLUSTERS_SCHEMA, validate_schema
from .base import ClusterEngine, ClusterParams, ClusterResult, EngineError

log = logging.getLogger(__name__)


# Word size rule of thumb per CD-HIT manual:
#   protein:    0.7..1.0 → -n 5,   0.6..0.7 → -n 4,   0.5..0.6 → -n 3,   0.4..0.5 → -n 2
#   nucleotide: 0.95..1.0 → -n 10, 0.9..0.95 → -n 8,  0.88..0.9 → -n 7,  0.85..0.88 → -n 6,
#               0.8..0.85 → -n 5, 0.75..0.8 → -n 4
def _word_size(identity: float, level: str) -> int:
    if level == "protein":
        if identity >= 0.7:
            return 5
        if identity >= 0.6:
            return 4
        if identity >= 0.5:
            return 3
        return 2
    # nucleotide (cd-hit-est)
    if identity >= 0.95:
        return 10
    if identity >= 0.90:
        return 8
    if identity >= 0.88:
        return 7
    if identity >= 0.85:
        return 6
    if identity >= 0.80:
        return 5
    return 4


class CDHitEngine(ClusterEngine):
    name = "cd-hit"
    supported_levels = ("protein", "nucleotide")

    def __init__(
        self,
        protein_binary: str = "cd-hit",
        nucleotide_binary: str = "cd-hit-est",
    ) -> None:
        self.protein_binary = protein_binary
        self.nucleotide_binary = nucleotide_binary

    def _binary_for_level(self, level: str) -> str:
        return self.protein_binary if level == "protein" else self.nucleotide_binary

    def check_available(self) -> None:
        # Check only the binary that would be used by default (protein).
        # The CLI re-checks for the chosen level at cluster() time.
        if shutil.which(self.protein_binary) is None:
            raise EngineError(
                f"{self.protein_binary!r} not found on PATH. Install via "
                f"`conda install -c bioconda cd-hit` or run inside the "
                f"DNMBcluster Docker image."
            )

    def cluster(
        self,
        input_fasta: Path,
        out_dir: Path,
        params: ClusterParams,
    ) -> ClusterResult:
        binary = self._binary_for_level(params.level)
        if shutil.which(binary) is None:
            raise EngineError(
                f"{binary!r} not found on PATH (needed for --level {params.level})."
            )
        if not self.supports(params.level):
            raise EngineError(f"cd-hit does not support level={params.level!r}")

        out_dir.mkdir(parents=True, exist_ok=True)
        rep_output = out_dir / "cdhit_out"
        clstr_output = out_dir / "cdhit_out.clstr"

        threads = params.threads if params.threads > 0 else (os.cpu_count() or 1)
        word_size = _word_size(params.identity, params.level)
        # CD-HIT memory flag is MB. 0 = unlimited. We pass a hard cap
        # when max_ram_gb is set to avoid OOMs on bounded containers.
        if params.max_ram_gb and params.max_ram_gb > 0:
            mem_mb = str(int(params.max_ram_gb * 1024 * 0.8))
        else:
            mem_mb = "0"

        cmd: list[str] = [
            binary,
            "-i", str(input_fasta),
            "-o", str(rep_output),
            "-c", f"{params.identity}",
            "-n", str(word_size),
            "-aS", f"{params.coverage}",
            "-M", mem_mb,
            "-T", str(threads),
            "-d", "0",  # preserve full headers (our integer protein_uid)
        ]

        log.info("running cd-hit: %s", " ".join(cmd))
        try:
            subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except subprocess.CalledProcessError as exc:
            stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
            raise EngineError(f"cd-hit failed (exit {exc.returncode}):\n{stderr}") from exc

        if not clstr_output.exists():
            raise EngineError(f"cd-hit did not produce {clstr_output}")

        clusters_parquet = out_dir / "clusters.parquet"
        n_clusters = parse_cdhit_clstr(
            clstr_path=clstr_output,
            out_parquet=clusters_parquet,
        )

        n_input = len(parse_fasta_headers(input_fasta))
        return ClusterResult(
            clusters_parquet=clusters_parquet,
            representatives_fasta=rep_output if rep_output.exists() else None,
            engine=self.name,
            params=params,
            n_input_sequences=n_input,
            n_clusters=n_clusters,
        )


# Matches a .clstr member line:
#   0\t350aa, >12345... *                          (centroid)
#   1\t348aa, >12346... at 92.15%                  (member w/ identity)
_MEMBER_RE = re.compile(
    r"^\s*\d+\s+\d+(?:aa|nt),\s+>(\S+)\.\.\.\s+(?:\*|at\s+([\d.]+)%)\s*$"
)
_CLUSTER_HEADER_RE = re.compile(r"^>Cluster\s+(\d+)\s*$")


def parse_cdhit_clstr(
    *,
    clstr_path: Path,
    out_parquet: Path,
) -> int:
    """Parse a CD-HIT ``.clstr`` file into the unified cluster schema.

    The ``.clstr`` format looks like::

        >Cluster 0
        0\t350aa, >12345... *
        1\t348aa, >12346... at 92.15%
        >Cluster 1
        0\t400aa, >23456... *
    """
    protein_uids: list[int] = []
    cluster_ids: list[int] = []
    pct_identities: list[float | None] = []
    representative_uids: list[int] = []
    is_centroid_flags: list[bool] = []

    current_cluster: int | None = None
    current_members: list[tuple[int, float | None, bool]] = []

    def _flush() -> None:
        if current_cluster is None or not current_members:
            return
        # Find centroid
        rep_uid: int | None = None
        for uid, _pct, is_cen in current_members:
            if is_cen:
                rep_uid = uid
                break
        if rep_uid is None:
            rep_uid = current_members[0][0]
        for uid, pct, is_cen in current_members:
            protein_uids.append(uid)
            cluster_ids.append(current_cluster)
            pct_identities.append(pct)
            representative_uids.append(rep_uid)
            is_centroid_flags.append(is_cen)

    with open(clstr_path) as fh:
        for line in fh:
            header_match = _CLUSTER_HEADER_RE.match(line)
            if header_match:
                _flush()
                current_cluster = int(header_match.group(1))
                current_members = []
                continue

            member_match = _MEMBER_RE.match(line)
            if member_match:
                raw_id = member_match.group(1)
                try:
                    uid = int(raw_id)
                except ValueError:
                    # Header wasn't our integer protein_uid — skip defensively.
                    log.warning("cd-hit non-integer header %r; skipping", raw_id)
                    continue
                pct_str = member_match.group(2)
                pct = float(pct_str) if pct_str else None
                is_cen = pct_str is None
                current_members.append((uid, pct, is_cen))

    _flush()

    genome_uids = [(u >> 48) & 0xFFFF for u in protein_uids]

    result = pa.table(
        {
            "protein_uid": pa.array(protein_uids, type=pa.uint64()),
            "genome_uid": pa.array(genome_uids, type=pa.uint16()),
            "cluster_id": pa.array(cluster_ids, type=pa.uint32()),
            "representative_uid": pa.array(representative_uids, type=pa.uint64()),
            "is_centroid": pa.array(is_centroid_flags, type=pa.bool_()),
            "pct_identity": pa.array(pct_identities, type=pa.float32()),
        },
        schema=CLUSTERS_SCHEMA,
    )
    validate_schema(result, CLUSTERS_SCHEMA, "clusters.parquet")

    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(result, out_parquet, compression="zstd", compression_level=3)
    return len(set(cluster_ids))
