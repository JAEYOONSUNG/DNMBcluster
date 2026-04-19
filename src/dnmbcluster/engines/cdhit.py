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
from dataclasses import dataclass
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from ..fasta import parse_fasta_headers
from ..schemas import CLUSTERS_SCHEMA, validate_schema
from .base import ClusterEngine, ClusterParams, ClusterResult, EngineError
from .sidecar import write_cdhit_native

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
        raw_dir: Path,
        processed_dir: Path,
        params: ClusterParams,
    ) -> ClusterResult:
        binary = self._binary_for_level(params.level)
        if shutil.which(binary) is None:
            raise EngineError(
                f"{binary!r} not found on PATH (needed for --level {params.level})."
            )
        if not self.supports(params.level):
            raise EngineError(f"cd-hit does not support level={params.level!r}")

        raw_dir.mkdir(parents=True, exist_ok=True)
        processed_dir.mkdir(parents=True, exist_ok=True)
        rep_output = raw_dir / "cdhit_out"
        clstr_output = raw_dir / "cdhit_out.clstr"

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
            "-p", "1",  # emit per-member alignment positions in .clstr
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

        clusters_parquet = processed_dir / "clusters.parquet"
        sidecar_parquet = processed_dir / "cdhit_native.parquet"
        n_clusters = parse_cdhit_clstr(
            clstr_path=clstr_output,
            out_parquet=clusters_parquet,
            sidecar_parquet=sidecar_parquet,
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


# Matches a .clstr member line in any of CD-HIT's four output flavors:
#
#   0\t350aa, >12345... *                                 (centroid)
#   1\t348aa, >12346... at 92.15%                         (default protein)
#   1\t348aa, >12346... at 1:336:5:340/92.15%             (-p 1 protein)
#   1\t1020nt, >12346... at +/92.15%                      (default nucleotide)
#   1\t1020nt, >12346... at 1:1020:5:1024/+/92.15%        (-p 1 nucleotide)
#
# Centroid groups: id only.
# Default groups: id + pct.
# -p 1 groups: id + qstart + qend + tstart + tend + (strand?) + pct.
_CENTROID_RE = re.compile(
    r"^\s*\d+\s+\d+(?:aa|nt),\s+>(\S+)\.\.\.\s+\*\s*$"
)
_MEMBER_DEFAULT_RE = re.compile(
    r"^\s*\d+\s+\d+(?:aa|nt),\s+>(\S+)\.\.\.\s+at\s+(?:([+\-])/)?([\d.]+)%\s*$"
)
_MEMBER_POSITIONAL_RE = re.compile(
    r"^\s*\d+\s+\d+(?:aa|nt),\s+>(\S+)\.\.\.\s+at\s+"
    r"(\d+):(\d+):(\d+):(\d+)"      # qstart:qend:tstart:tend
    r"(?:/([+\-]))?"                # optional strand
    r"/([\d.]+)%\s*$"
)
_CLUSTER_HEADER_RE = re.compile(r"^>Cluster\s+(\d+)\s*$")


@dataclass
class _Member:
    uid: int
    pct: float | None
    is_centroid: bool
    query_start: int | None = None
    query_end: int | None = None
    target_start: int | None = None
    target_end: int | None = None
    strand: str | None = None


def parse_cdhit_clstr(
    *,
    clstr_path: Path,
    out_parquet: Path,
    sidecar_parquet: Path | None = None,
) -> int:
    """Parse a CD-HIT ``.clstr`` file into the unified cluster schema.

    The ``.clstr`` format looks like::

        >Cluster 0
        0\t350aa, >12345... *
        1\t348aa, >12346... at 92.15%
        >Cluster 1
        0\t400aa, >23456... *

    With ``-p 1`` member lines also include alignment positions, which
    are preserved in the sidecar parquet when ``sidecar_parquet`` is
    given. See ``engines/sidecar.py`` for the rationale.
    """


    protein_uids: list[int] = []
    cluster_ids: list[int] = []
    pct_identities: list[float | None] = []
    representative_uids: list[int] = []
    is_centroid_flags: list[bool] = []
    sidecar_rows: list[dict] = []

    current_cluster: int | None = None
    current_members: list[_Member] = []

    def _flush() -> None:
        if current_cluster is None or not current_members:
            return
        rep_uid: int | None = None
        for m in current_members:
            if m.is_centroid:
                rep_uid = m.uid
                break
        if rep_uid is None:
            rep_uid = current_members[0].uid
        for m in current_members:
            protein_uids.append(m.uid)
            cluster_ids.append(current_cluster)
            pct_identities.append(m.pct)
            representative_uids.append(rep_uid)
            is_centroid_flags.append(m.is_centroid)
            if not m.is_centroid:
                sidecar_rows.append(
                    {
                        "protein_uid": m.uid,
                        "representative_uid": rep_uid,
                        "native_pct_identity": m.pct,
                        "query_start": m.query_start,
                        "query_end": m.query_end,
                        "target_start": m.target_start,
                        "target_end": m.target_end,
                        "strand": m.strand,
                    }
                )

    def _parse_member(line: str) -> _Member | None:
        # Centroid (no identity, no positions)
        cm = _CENTROID_RE.match(line)
        if cm:
            try:
                return _Member(uid=int(cm.group(1)), pct=None, is_centroid=True)
            except ValueError:
                log.warning("cd-hit non-integer header %r; skipping", cm.group(1))
                return None

        # -p 1 positional flavor: try first because the default regex
        # would also match an unanchored prefix of these lines.
        pm = _MEMBER_POSITIONAL_RE.match(line)
        if pm:
            try:
                uid = int(pm.group(1))
            except ValueError:
                log.warning("cd-hit non-integer header %r; skipping", pm.group(1))
                return None
            return _Member(
                uid=uid,
                pct=float(pm.group(7)),
                is_centroid=False,
                query_start=int(pm.group(2)),
                query_end=int(pm.group(3)),
                target_start=int(pm.group(4)),
                target_end=int(pm.group(5)),
                strand=pm.group(6),
            )

        dm = _MEMBER_DEFAULT_RE.match(line)
        if dm:
            try:
                uid = int(dm.group(1))
            except ValueError:
                log.warning("cd-hit non-integer header %r; skipping", dm.group(1))
                return None
            return _Member(
                uid=uid,
                pct=float(dm.group(3)),
                is_centroid=False,
                strand=dm.group(2),
            )
        return None

    with open(clstr_path) as fh:
        for line in fh:
            header_match = _CLUSTER_HEADER_RE.match(line)
            if header_match:
                _flush()
                current_cluster = int(header_match.group(1))
                current_members = []
                continue

            member = _parse_member(line)
            if member is not None:
                current_members.append(member)

    _flush()

    genome_uids = [(u >> 48) & 0xFFFF for u in protein_uids]
    n_rows = len(protein_uids)
    null_f32: list[float | None] = [None] * n_rows
    null_u32: list[int | None] = [None] * n_rows

    # CD-HIT ``at X%`` identity is preserved in the native .clstr file
    # on disk. The unified clusters.parquet schema is populated by the
    # shared engines.realign stage so every engine produces identical
    # metrics — per-engine identity variants are intentionally dropped
    # here to avoid a mixed-source pct_identity_fwd column.
    _ = pct_identities  # preserved only for debugging / tests

    result = pa.table(
        {
            "protein_uid": pa.array(protein_uids, type=pa.uint64()),
            "genome_uid": pa.array(genome_uids, type=pa.uint16()),
            "cluster_id": pa.array(cluster_ids, type=pa.uint32()),
            "representative_uid": pa.array(representative_uids, type=pa.uint64()),
            "is_centroid": pa.array(is_centroid_flags, type=pa.bool_()),
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

    if sidecar_parquet is not None:
        write_cdhit_native(rows=sidecar_rows, out_path=sidecar_parquet)

    return len(set(cluster_ids))
