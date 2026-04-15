"""Pluggable clustering engine backends.

All engines share the `ClusterEngine` interface defined in ``base.py``
and emit ``clusters.parquet`` in the unified schema from
``..schemas.CLUSTERS_SCHEMA``.
"""
from __future__ import annotations

from .base import ClusterEngine, ClusterParams, ClusterResult, EngineError


def get_engine(name: str) -> ClusterEngine:
    """Return a clustering engine instance by its ``--tool`` name."""
    if name == "mmseqs2":
        from .mmseqs2 import MMseqs2Engine
        return MMseqs2Engine()
    if name == "diamond":
        from .diamond import DiamondEngine
        return DiamondEngine()
    if name == "cd-hit":
        from .cdhit import CDHitEngine
        return CDHitEngine()
    if name == "usearch12":
        from .usearch12 import Usearch12Engine
        return Usearch12Engine()
    raise EngineError(f"unknown engine {name!r}")


__all__ = [
    "ClusterEngine",
    "ClusterParams",
    "ClusterResult",
    "EngineError",
    "get_engine",
]
