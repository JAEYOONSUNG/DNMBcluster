"""Abstract clustering engine interface.

Every engine (MMseqs2, DIAMOND deepclust, CD-HIT, usearch12) accepts a
DNMB FASTA file (integer-header, see fasta.py) and emits
``clusters.parquet`` matching ``schemas.CLUSTERS_SCHEMA``. The unified
schema decouples downstream code from the specific clustering backend.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

Level = Literal["protein", "nucleotide"]


@dataclass
class ClusterParams:
    """Normalized clustering parameters passed to every engine."""

    identity: float = 0.5
    coverage: float = 0.8
    threads: int = 0
    max_ram_gb: float | None = None
    level: Level = "protein"
    #: If True, engines that support a separate alignment pass will run it
    #: and populate ``member_coverage`` / ``rep_coverage`` / ``alignment_length``
    #: / bidirectional ``pct_identity`` in clusters.parquet. If False, only
    #: the fast clustering step runs and those fields are null.
    with_alignment: bool = True


@dataclass
class ClusterResult:
    """Handle to a completed clustering run."""

    clusters_parquet: Path
    representatives_fasta: Path | None
    engine: str
    params: ClusterParams
    n_input_sequences: int = 0
    n_clusters: int = 0
    tmpdir: Path | None = None
    extra: dict[str, object] = field(default_factory=dict)


class EngineError(RuntimeError):
    """Raised when a clustering engine fails or is unavailable."""


class ClusterEngine(ABC):
    """Base class for all clustering backends."""

    #: Short unique identifier used by the CLI `--tool` option.
    name: str = ""

    #: Sequence levels this engine supports. DIAMOND overrides to
    #: ``("protein",)`` only; others default to both.
    supported_levels: tuple[Level, ...] = ("protein", "nucleotide")

    def supports(self, level: Level) -> bool:
        return level in self.supported_levels

    def check_available(self) -> None:
        """Raise `EngineError` if the backend binary is not installed.

        Subclasses implement this by shelling out to a version probe.
        The CLI calls it before any parsing so we fail fast on missing
        dependencies rather than mid-pipeline.
        """

    @abstractmethod
    def cluster(
        self,
        input_fasta: Path,
        raw_dir: Path,
        processed_dir: Path,
        params: ClusterParams,
    ) -> ClusterResult:
        """Cluster the input FASTA.

        ``raw_dir`` is the engine's scratch + native output directory:
        anything the underlying binary writes (``.uc``, ``.clstr``,
        ``deepclust`` TSV, MMseqs2 work files, scratch tmpdirs) lives
        here untouched, so users can audit the engine's own bookkeeping.

        ``processed_dir`` is where the engine writes the **canonical**
        Parquet artifacts every downstream stage consumes: at minimum
        ``clusters.parquet`` matching ``schemas.CLUSTERS_SCHEMA``, and
        any engine-native sidecar (e.g. ``cdhit_native.parquet``).
        """
