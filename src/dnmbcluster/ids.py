"""Identifier extraction and fallback-chain resolution.

Implements the two-level `(genome_key, cds_key)` composite primary key
described in SPEED.md Section 1, plus the compact integer `protein_uid`
hot-path alias.

No single GenBank field is globally unique in a pan-genome context:

- `locus_tag` is reassigned on reannotation and absent in old files.
- `protein_id` (RefSeq WP_) is intentionally shared across closely
  related strains — it collides in the pan-genome setting.
- Filename alone is brittle across filesystems and not meaningful.

The fix is a composite key with fallback chains on each half and a
`_source` enum recording which rung was used. This module owns the
fallback logic; `genbank.py` calls it once per GenBank file / CDS.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional


# ---------------------------------------------------------------------------
# Assembly accession regex (RefSeq GCF_ / GenBank GCA_, 9 digits, optional .N)
# ---------------------------------------------------------------------------

ASSEMBLY_RE_STRICT = re.compile(r"^(GC[FA]_\d{9})(\.\d+)?$")
ASSEMBLY_RE_ANYWHERE = re.compile(r"(GC[FA]_\d{9})(\.\d+)?")


def split_assembly(accession: str) -> tuple[Optional[str], Optional[str]]:
    """Split ``'GCF_000009785.1'`` into ``('GCF_000009785', '.1')``.

    Returns ``(None, None)`` if the string is not an assembly accession.
    """
    m = ASSEMBLY_RE_STRICT.match(accession)
    if not m:
        return None, None
    return m.group(1), (m.group(2) or None)


# ---------------------------------------------------------------------------
# genome_key — two-level primary key, upper half
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GenomeKey:
    """Result of `resolve_genome_key` — key, source, and assembly parts."""

    key: str
    source: str  # one of GENOME_KEY_SOURCES in schemas.py
    assembly_prefix: Optional[str] = None
    assembly_version: Optional[str] = None


def _from_dblink(dblink_entries: Iterable[str]) -> Optional[str]:
    """Extract an assembly accession from GenBank DBLINK entries.

    DBLINK entries look like ``'BioProject: PRJNA123'`` or
    ``'Assembly: GCF_000009785.1'``; we only care about the Assembly line.
    """
    for entry in dblink_entries:
        if "Assembly:" not in entry:
            continue
        tail = entry.split("Assembly:", 1)[1].strip()
        m = ASSEMBLY_RE_ANYWHERE.search(tail)
        if m:
            prefix = m.group(1)
            version = m.group(2) or ""
            return f"{prefix}{version}"
    return None


def _from_filename_gcf(file_path: Path) -> Optional[str]:
    """Extract GCF/GCA assembly accession from a GenBank filename."""
    stem = file_path.name
    for ext in (".gbff", ".gbk", ".gb", ".gz"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
    m = ASSEMBLY_RE_ANYWHERE.search(stem)
    if not m:
        return None
    return f"{m.group(1)}{m.group(2) or ''}"


def _slugify_organism(organism: str, strain: Optional[str]) -> Optional[str]:
    """Build a filesystem-safe slug from organism + strain."""
    if not organism:
        return None
    parts = [organism]
    if strain and strain not in organism:
        parts.append(strain)
    raw = "_".join(parts)
    slug = re.sub(r"[^A-Za-z0-9]+", "_", raw).strip("_")
    return slug or None


def resolve_genome_key(
    *,
    file_path: Path,
    dblink: Iterable[str] = (),
    locus_version: Optional[str] = None,
    organism: Optional[str] = None,
    strain: Optional[str] = None,
    manifest_name: Optional[str] = None,
) -> GenomeKey:
    """Apply the 6-step genome_key fallback chain (SPEED.md Section 1).

    Priority, highest first:

    1. ``manifest_name``       → ``source='manifest'``
    2. DBLINK ``Assembly:``    → ``source='dblink'``
    3. filename GCF/GCA match  → ``source='filename_gcf'``
    4. LOCUS/VERSION of first record → ``source='locus_version'``
    5. organism + strain slug  → ``source='organism'``
    6. filename stem (raw)     → ``source='filename_raw'``

    Always returns a `GenomeKey`; never raises.
    """
    if manifest_name:
        prefix, version = split_assembly(manifest_name)
        return GenomeKey(manifest_name, "manifest", prefix, version)

    dblink_acc = _from_dblink(dblink)
    if dblink_acc:
        prefix, version = split_assembly(dblink_acc)
        return GenomeKey(dblink_acc, "dblink", prefix, version)

    filename_acc = _from_filename_gcf(file_path)
    if filename_acc:
        prefix, version = split_assembly(filename_acc)
        return GenomeKey(filename_acc, "filename_gcf", prefix, version)

    if locus_version:
        return GenomeKey(locus_version, "locus_version")

    slug = _slugify_organism(organism or "", strain)
    if slug:
        return GenomeKey(slug, "organism")

    return GenomeKey(file_path.stem, "filename_raw")


# ---------------------------------------------------------------------------
# cds_key — two-level primary key, lower half (scoped within one genome)
# ---------------------------------------------------------------------------


def resolve_cds_key(
    *,
    locus_tag: Optional[str],
    protein_id: Optional[str],
    gene: Optional[str],
    gene_ordinal: int,
    contig: str,
    start: int,
    end: int,
    strand: int,
) -> tuple[str, str]:
    """Apply the 4-step cds_key fallback chain.

    Priority, highest first:

    1. ``locus_tag``                        → ``source='locus_tag'``
    2. ``protein_id``                       → ``source='protein_id'``
    3. ``f'{gene}__{ordinal}'``             → ``source='gene_ordinal'``
    4. ``f'cds_{contig}_{start}_{end}_{strand}'`` → ``source='coord'``

    Returns ``(key, source)``.
    """
    if locus_tag:
        return locus_tag, "locus_tag"
    if protein_id:
        return protein_id, "protein_id"
    if gene:
        return f"{gene}__{gene_ordinal}", "gene_ordinal"
    strand_char = "+" if strand >= 0 else "-"
    return f"cds_{contig}_{start}_{end}_{strand_char}", "coord"


# ---------------------------------------------------------------------------
# protein_uid — compact integer hot-path key
# ---------------------------------------------------------------------------

# Layout: [genome_uid:16 | unused:16 | gene_uid:32]
_GENOME_SHIFT = 48
_GENE_MASK = 0xFFFFFFFF
_GENOME_MASK = 0xFFFF
_GENOME_MAX = 1 << 16
_GENE_MAX = 1 << 32


def make_protein_uid(genome_uid: int, gene_uid: int) -> int:
    """Pack ``(genome_uid:uint16, gene_uid:uint32)`` into a uint64.

    See SPEED.md Section 1: this is a run-scoped hot-path alias for the
    correctness-layer ``(genome_key, cds_key)`` composite. The middle 16
    bits are reserved for future flags and must be zero.
    """
    if not (0 <= genome_uid < _GENOME_MAX):
        raise ValueError(
            f"genome_uid {genome_uid} out of range [0, {_GENOME_MAX})"
        )
    if not (0 <= gene_uid < _GENE_MAX):
        raise ValueError(
            f"gene_uid {gene_uid} out of range [0, {_GENE_MAX})"
        )
    return (genome_uid << _GENOME_SHIFT) | gene_uid


def unpack_protein_uid(protein_uid: int) -> tuple[int, int]:
    """Inverse of `make_protein_uid`."""
    return (protein_uid >> _GENOME_SHIFT) & _GENOME_MASK, protein_uid & _GENE_MASK
