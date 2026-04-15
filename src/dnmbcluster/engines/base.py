"""Abstract clustering engine interface.

All concrete engines (MMseqs2, DIAMOND deepclust, CD-HIT, usearch12) implement
this interface and produce output in the unified DNMB schema:

    cluster_id:int64, member_id:str, genome_id:str,
    representative_id:str, is_centroid:bool,
    pct_identity:float32|NA, member_length:int32
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ClusterResult:
    """Handle to a completed clustering run."""

    clusters_tsv: Path
    """Unified cluster membership table (one row per input protein)."""

    representatives_fasta: Path
    """FASTA of cluster representatives."""

    engine: str
    params: dict = field(default_factory=dict)


class ClusterEngine(ABC):
    """Base class for clustering backends."""

    name: str

    @abstractmethod
    def cluster(
        self,
        input_fasta: Path,
        out_dir: Path,
        identity: float = 0.5,
        coverage: float = 0.8,
        threads: int = 0,
    ) -> ClusterResult:
        """Cluster the given protein FASTA and emit unified output."""
