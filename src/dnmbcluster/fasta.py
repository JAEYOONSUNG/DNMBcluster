"""FASTA writer for the clustering stage.

Reads `gene_table.parquet` and emits a single protein (or nucleotide)
FASTA file where each header is the decimal `protein_uid`. Short
integer headers mean:

- Clustering tools parse headers ~2× faster.
- FASTA on disk is ~60 MB smaller at 3M-protein scale.
- The unified clusters.parquet parser gets pure integer columns —
  no string hashing in the hot path.

SPEED.md Sections 1 and 4.
"""
from __future__ import annotations

from pathlib import Path
from typing import Literal

import pyarrow.parquet as pq

from .schemas import gene_table_schema, validate_schema

Level = Literal["protein", "nucleotide"]


def _sequence_column(level: Level) -> str:
    return "translation" if level == "protein" else "nt_sequence"


def write_fasta(
    gene_table_path: Path,
    output_path: Path,
    *,
    level: Level = "protein",
    sort_by_length_desc: bool = True,
) -> tuple[int, int]:
    """Write a FASTA file from ``gene_table.parquet``.

    Headers are ``>{protein_uid}`` (decimal). Sequences are wrapped to
    80 columns. Returns ``(n_sequences, n_bytes_written)``.

    If ``sort_by_length_desc`` is True (default), sequences are written
    in descending length order. This helps usearch12 ``cluster_fast``
    and is neutral for MMseqs2 / DIAMOND / CD-HIT.
    """
    table = pq.read_table(gene_table_path)
    validate_schema(table, gene_table_schema(level), "gene_table.parquet")

    seq_col_name = _sequence_column(level)
    protein_uids = table.column("protein_uid").to_pylist()
    sequences = table.column(seq_col_name).to_pylist()

    if sort_by_length_desc:
        order = sorted(
            range(len(sequences)),
            key=lambda i: len(sequences[i]),
            reverse=True,
        )
    else:
        order = range(len(sequences))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_seq = 0
    n_bytes = 0
    with open(output_path, "w") as out:
        for i in order:
            header = f">{protein_uids[i]}\n"
            seq = sequences[i]
            # 80-column wrapping — cheap and preserves compatibility with
            # older parsers. Most modern tools tolerate any line length,
            # but the overhead is negligible.
            wrapped_lines = [seq[j : j + 80] for j in range(0, len(seq), 80)]
            body = "\n".join(wrapped_lines) + "\n"
            out.write(header)
            out.write(body)
            n_seq += 1
            n_bytes += len(header) + len(body)

    return n_seq, n_bytes


def parse_fasta_headers(path: Path) -> list[int]:
    """Return the list of ``protein_uid`` values from a DNMB FASTA.

    Used by tests and by engine parsers that need to cross-check that
    their output covers the input set.
    """
    uids: list[int] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                uids.append(int(line[1:].split()[0]))
    return uids
