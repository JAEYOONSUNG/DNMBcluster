"""Shared I/O helpers for crash-safe persistence.

Every parquet in DNMBcluster is downstream input for another stage —
corrupt files from a killed mid-write kill ``--resume`` flows. The
helpers here write to a sibling ``.tmp`` file then ``os.replace()``
into place. ``os.replace`` is atomic on POSIX (and on Windows for
same-filesystem renames), so a reader either sees the old file or
the fully-flushed new one, never a half-written intermediate.
"""
from __future__ import annotations

import os
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


def atomic_write_table(
    table: pa.Table,
    out_path: Path,
    *,
    compression: str = "zstd",
    compression_level: int = 3,
    **kwargs,
) -> None:
    """Atomically write ``table`` to ``out_path`` as parquet.

    Writes to ``{out_path}.tmp``, flushes, then ``os.replace``s into
    place. On failure mid-write, the temp file is removed and the
    original (if any) is untouched.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    try:
        pq.write_table(
            table,
            tmp_path,
            compression=compression,
            compression_level=compression_level,
            **kwargs,
        )
        os.replace(tmp_path, out_path)
    except BaseException:
        # Best-effort cleanup of the temp file on any failure path
        # (exception, KeyboardInterrupt). The rename is the "commit
        # point" — if we never reach it, the destination stays as it
        # was before this call.
        try:
            tmp_path.unlink(missing_ok=True)
        except OSError:
            pass
        raise
