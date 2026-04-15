"""Pluggable clustering engine backends.

All engines share the `ClusterEngine` interface defined in ``base.py``
and emit ``clusters.parquet`` in the unified schema from
``..schemas.CLUSTERS_SCHEMA``.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from .base import ClusterEngine, ClusterParams, ClusterResult, EngineError

if TYPE_CHECKING:
    pass


def get_engine(name: str) -> ClusterEngine:
    """Return a clustering engine instance by its ``--tool`` name."""
    if name == "mmseqs2":
        from .mmseqs2 import MMseqs2Engine
        return MMseqs2Engine()
    # Stubs for engines not yet implemented — M3.
    if name in {"diamond", "cd-hit", "usearch12"}:
        raise EngineError(f"engine {name!r} not yet implemented (coming in M3)")
    raise EngineError(f"unknown engine {name!r}")


__all__ = [
    "ClusterEngine",
    "ClusterParams",
    "ClusterResult",
    "EngineError",
    "get_engine",
]
