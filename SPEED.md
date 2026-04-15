# DNMBcluster — Speed & Data-Integrity Design

This document is load-bearing. Every M1–M7 implementation decision must
respect the rules here. If a rule conflicts with reality (benchmarked),
update this doc first — never silently diverge.

## North-star targets

| Workload | Wall time target | Peak RSS target |
|---|---|---|
| 10 Geobacillus genomes (~30K proteins) | **< 30 s** end-to-end | < 2 GB |
| 100 genomes (~300K proteins) | **< 5 min** end-to-end | < 8 GB |
| 1000 genomes (~3M proteins) | **< 45 min** end-to-end | < 32 GB |

End-to-end = GenBank parse + cluster + dataframe + plots.

If a build regresses any target by >20%, CI fails.

---

## 1. Identifier scheme — two-level composite + compact integer hot path

The identifier layer has two jobs:

1. **Correctness layer** — uniquely name every CDS across every genome in
   a run, without silent collisions and without trusting any single
   GenBank field. This is the `(genome_key, cds_key)` composite, the
   field consensus across Roary / Panaroo / PPanGGOLiN / get_homologues /
   anvi'o.
2. **Speed layer** — make all hot-path joins work on `uint64` integers.
   Strings exist only at parse time and at final-report time.

Both layers coexist in `id_map.parquet`. The integer `protein_uid` is
the hot-path key; the `(genome_key, cds_key)` pair is the human-readable
primary key stored alongside it.

### Why neither field alone works (motivation)

- **`locus_tag`**: non-redundant within a file, but reassigned on
  reannotation. Absent in pre-2005 GenBanks and some Prokka edge cases.
  Can duplicate within a broken-rerun file.
- **`protein_id`** (WP_xxxxxxx.x): RefSeq deliberately **shares** WP_
  accessions across closely related strains with identical protein
  sequences. Useless as a global key in pan-genome context; this is
  why BPGA had to prefix it with a per-run `strain_number`.
- **`gene`**: duplicated by paralogs.
- **Filename**: brittle across filesystems, not meaningful.

→ Any single field is wrong. The fix is `(genome_key, cds_key)` with a
fallback chain on each half and a `_source` enum recording which rung
of the chain was used.

### `genome_key` extraction priority (per GenBank file)

```
1. --manifest entry                  user-supplied name in optional TSV
2. DBLINK Assembly accession         "GCF_000009785.1" (RefSeq) / GCA_ (GenBank)
3. Filename stem matching GCF/GCA regex    "GCF_030376785"   (Prokka-produced GBK has no DBLINK)
4. LOCUS / VERSION of first record   "CP187452.1"   (single-replicon only)
5. Organism + strain slug            "Geobacillus_sp_strain_X"
6. Filename stem (raw)               last-resort
```

Rules:
- **Assembly accession includes the version suffix** (`.1`, `.2`). Two
  versions of the same assembly → two different `genome_key`s. This is
  how we detect reannotation and force a rebuild instead of silently
  mixing old and new.
- `assembly_prefix` (the `GCF_000009785` without `.1`) is stored as a
  **separate column** so we can warn the user: *"genome_key changed from
  `.1` to `.2`, this is a new run."*
- Every genome records `genome_key_source:enum` — one of `manifest /
  dblink / filename_gcf / locus_version / organism / filename_raw` — so
  cluster inspection never leaves the user wondering how a label was
  derived.
- One GenBank file = one genome, regardless of how many LOCUS records
  it contains (chromosome + plasmids + WGS contigs). We never key the
  genome by per-record VERSION.

### `cds_key` fallback chain (per CDS, scoped inside a genome)

```
1. locus_tag                         "GK0001"
2. protein_id                        "WP_012345678.1"   (only as within-genome fallback — safe because scoped)
3. gene + ordinal                    "dnaA__3"          (when the same gene symbol recurs)
4. synthesized coordinate key        "cds_CP187452.1_3421_4560_-1"   (pseudogenes, legacy GB)
```

Every CDS records `cds_key_source:enum` — one of `locus_tag /
protein_id / gene_ordinal / coord`. Downstream QC can filter out
`coord`-derived entries if desired.

Duplicate `locus_tag` within a single file (broken Prokka rerun) is
detected at parse time and resolved by appending an ordinal suffix
(`GK0001__dup2`) with a warning logged to `manifest.json`.

### Compact integer hot-path IDs (unchanged from previous version)

```
genome_uid  : uint16   (≤ 65 535 genomes per run)
gene_uid    : uint32   (per-genome CDS index)
protein_uid : uint64 = (genome_uid << 48) | gene_uid
```

- Polars/Arrow joins on `int64` are ~3–5× faster than on `Utf8`.
- FASTA headers shrink from ~35 bytes to ~15 bytes; on 3M-protein input
  that's ~60 MB less for the clustering tool to parse.
- Clustering `.uc` / `_cluster.tsv` parsing becomes pure integer work.
- Presence/absence matrix keys are native `int64` — no string hashing.

The integer IDs are **run-scoped, not run-stable**. Two runs on the
same input will produce the same `genome_key` but potentially different
`genome_uid` (depends on file discovery order). The composite
`(genome_key, cds_key)` is the stable primary key; `protein_uid` is a
fast alias that lives only for the duration of a run.

### `id_map.parquet` — full schema

Written once after GenBank parsing; all downstream joins consult it.

```
── hot-path integer keys ─────────────────────────────────────────────
protein_uid            : uint64    primary hot-path key
genome_uid             : uint16
gene_uid               : uint32

── correctness-layer composite key ───────────────────────────────────
genome_key             : Utf8      NOT NULL  (primary human-readable key)
genome_key_source      : enum      NOT NULL  {manifest, dblink, filename_gcf,
                                               locus_version, organism, filename_raw}
cds_key                : Utf8      NOT NULL
cds_key_source         : enum      NOT NULL  {locus_tag, protein_id,
                                               gene_ordinal, coord}

── reannotation detection ────────────────────────────────────────────
assembly_prefix        : Utf8      nullable  (e.g. "GCF_000009785" sans version)
assembly_version       : Utf8      nullable  (e.g. ".1")

── original GenBank attributes (non-unique, for display / dereplication) ──
organism               : Utf8
strain                 : Utf8
locus_tag              : Utf8      nullable
protein_id             : Utf8      nullable  (WP_*; useful cross-genome as a
                                               "same identical protein" hint)
gene                   : Utf8      nullable
product                : Utf8      nullable
ec_number              : Utf8      nullable
contig                 : Utf8
start                  : uint32
end                    : uint32
strand                 : int8
aa_length              : uint32
```

`protein_id` is deliberately kept as a **non-key attribute**. It's the
right place for RefSeq WP_ accessions: it lets downstream queries find
"all CDSs across genomes that are the exact same RefSeq protein" in one
Polars filter — effectively a free dereplication hint — without
corrupting the primary key.

### Label format at report time

Cluster members in user-facing output (plots, HTML, Excel) are labelled
`{genome_key}::{cds_key}`, e.g.:

```
GCF_000009785.1::GK0001
GCF_030376785.1::dnaA__1
CP187452.1::cds_CP187452.1_3421_4560_-1
```

Stable, meaningful, collision-free, and traceable back to the original
annotation through `cds_key_source`.

---

## 2. Vectorized tabular layer (Polars + Arrow only)

### Non-negotiable rules

1. **No pandas in hot paths.** pandas is allowed only for `to_excel`
   export at the very end, and only via `polars → pandas → xlsx`.
2. **No Python loops on DataFrames >1000 rows.** If you need per-row
   logic, it's a Polars expression, a `.map_batches`, or a NumPy kernel.
3. **Lazy where possible.** Use `pl.scan_parquet` / `scan_csv` so query
   optimizer can push down projections and predicates.
4. **Parquet with ZSTD-3** for all intermediate tables. Row-group size
   = 128 MB default, tune if profiling says otherwise.
5. **Arrow, not CSV, is the wire format** between Python and R. R reads
   via `arrow::read_parquet` (also ~5× faster than `read.table`).

### Core dataframes (all Polars-schema-enforced)

```
gene_table.parquet         protein_uid, genome_uid, gene_uid, start, end,
                           strand, length, aa_length, Mw, pI, product,
                           translation
clusters.parquet           protein_uid, cluster_id, representative_uid,
                           is_centroid, pct_identity
genome_meta.parquet        genome_uid, genome_accession, organism, n_cds,
                           gc_percent, assembly_length
presence_absence.parquet   cluster_id, genome_bitmap (uint64[]), n_genomes
pan_core_curve.parquet     permutation, k, pan, core
cluster_summary.parquet    cluster_id, n_genomes, category (core/shell/cloud),
                           rep_product, rep_protein_id
```

Schemas live in `src/dnmbcluster/schemas.py` as `pa.schema(...)`
constants. Every writer validates before writing; every reader validates
on load. `schema.equals(expected, check_metadata=False)` is the gate.

---

## 3. GenBank parser — custom fast path

Biopython `SeqIO.parse` on GenBank is known-slow (pure Python, regex-heavy).
Measured baseline: ~2 s per typical bacterial genome (~3–6 MB .gbff).
At 1000 genomes that's 30+ minutes wasted before clustering even starts.

### Implementation strategy

1. **First pass:** try `pyrodigal-gv` / `pyhmmer` style C-backed parse. If
   there's no good C parser for GenBank CDS extraction, write our own.
2. **Custom parser** (expected winning path):
   - Byte-level line iteration (`open(..., 'rb')`, `.split(b'\n')`).
   - State machine: scan for `CDS ` blocks, pull `/locus_tag=`,
     `/protein_id=`, `/product=`, `/translation="..."` (multi-line).
   - Skip the entire ORIGIN section — we never need genomic sequence
     in the hot path. (Separate pass for `nt_seq` only if requested.)
   - Target: ≥20× Biopython on our Geobacillus fixture.
3. **Parallelism:** `ProcessPoolExecutor(max_workers=nproc)` one worker
   per file. GenBank files are independent.
4. **Output:** directly emit Parquet row batches via `pyarrow.RecordBatch`
   — do not stage through Python lists/dicts of per-gene dicts.

### Validation

Every parsed file emits a `sha256` of the source bytes + a `schema_hash`
of the output RecordBatch schema into `manifest.json`. Re-runs that find
a matching hash skip parsing entirely (cache hit).

---

## 4. Clustering engine invocation rules

### Shared flags

- **Threads:** always pass `--threads $(nproc)` (or the CLI `--threads`
  override). Never let the tool default.
- **Input FASTA pre-sort:** sort by length descending before invocation.
  UCLUST/cluster_fast benefit directly; others benefit from better
  memory locality.
- **Tmpdir on tmpfs:** Docker `--tmpfs /tmp:size=8g`. Clustering tools
  write GB of scratch; disk-backed /tmp is a 2–3× slowdown.
- **Header format:** `>protein_uid` where protein_uid is the decimal
  int64. Clustering tools never see locus_tags.
- **No stdout capture** for tool output; always `-o FILE`. Python reads
  the file afterward via `pl.scan_csv`.

### Per-engine specifics

| Engine | Default command | Key flag |
|---|---|---|
| MMseqs2 | `easy-linclust` | `--min-seq-id 0.5 -c 0.8 --cov-mode 0 --split-memory-limit 80%` |
| DIAMOND | `diamond deepclust` | `--approx-id 50 --member-cover 80 -M 75%RAM --header` |
| CD-HIT | `cd-hit` | `-c 0.5 -n 3 -aS 0.8 -d 0 -M 0 -T $(nproc)` |
| usearch12 | `cluster_mt` if ≥8 cores else `cluster_fast` | `-id 0.5 -sort length -uc $OUT` |

`--split-memory-limit 80%` / `-M 75%RAM` lets the tool self-manage RAM
and avoid OS OOM — critical for reproducibility.

### Output parsing

All parsers emit the unified schema in one Polars query (no row loops):

```python
# MMseqs2 example (representative example — real code in engines/mmseqs2.py)
lf = (
    pl.scan_csv(tsv, separator="\t", has_header=False,
                new_columns=["rep_str", "member_str"])
      .with_columns(
          pl.col("rep_str").cast(pl.UInt64).alias("representative_uid"),
          pl.col("member_str").cast(pl.UInt64).alias("protein_uid"),
      )
      .drop(["rep_str", "member_str"])
      .with_columns(
          (pl.col("representative_uid") == pl.col("protein_uid")).alias("is_centroid"),
      )
      # cluster_id = dense rank of representative_uid
      .with_columns(
          pl.col("representative_uid").rank("dense").cast(pl.UInt32).alias("cluster_id"),
      )
)
lf.sink_parquet(out / "clusters.parquet", compression="zstd")
```

One streaming query, no materialization until sink. Same pattern for
every engine.

---

## 5. Presence/absence + pan/core — bitmap arithmetic

### Representation

For N genomes ≤ 64: **one `uint64` per cluster = bitmap of genome
membership.** Core = `bitmap == all_ones_mask`. Singleton cloud =
`popcount(bitmap) == 1`.

For N genomes > 64: **`np.ndarray[uint64]` of ceil(N/64) words per
cluster.** Stored in Parquet as `list[uint64]`.

### Pan/core curve (the traditional O(N²) trap)

Naive: for each permutation π, for k in 1..N, intersect/union
cumulative bitmap. N=1000, 100 perms → 100K bitmap ops, each is a few
uint64 word ops with popcount. **NumPy vectorized ⇒ well under 1 s.**

```python
# cluster_bitmaps: np.ndarray[uint64, shape=(n_clusters, words)]
# perm:            np.ndarray[uint16, shape=(N,)]
# Build cumulative union/intersection in-place with np.bitwise_or.accumulate
# over bits selected by perm — vectorized over clusters.
```

Never use Python loops over clusters or genomes in this stage.

### Bootstrap parallelism

100 permutations → `joblib.Parallel(n_jobs=nproc, prefer="threads")`
because bitmap ops release the GIL (NumPy C code). Avoid process pools
here — the bitmap array is large and copying hurts.

---

## 6. Data integrity guarantees

### Manifest

`results/manifest.json` records, for every stage:

```json
{
  "stage": "genbank_parse",
  "inputs": [{"path": "genbank/foo.gbff", "sha256": "...", "bytes": 4821033}],
  "outputs": [{"path": "gene_table.parquet", "sha256": "...", "rows": 3982}],
  "schema_hash": "arrow-fingerprint-hex",
  "started": "...", "finished": "...", "engine_version": "mmseqs2=15.6f452"
}
```

On re-run, if all input hashes match a prior entry, skip the stage.
`--force` disables cache.

### Schema validation

`pyarrow.Schema.equals(expected, check_metadata=False)` at every
reader/writer boundary. Missing column → fail fast with the exact column
name and expected type. No silent coercion.

### Property tests

`tests/test_integrity.py` uses `hypothesis` to generate fake GenBank
records and assert:

- Parser round-trips: parse → emit → re-parse → identical RecordBatch.
- ID-map invertibility: for every protein_uid in clusters.parquet, the
  id_map has exactly one matching row.
- Presence/absence consistency: `bitmap[genome_uid]` set ⇒ there exists
  at least one row in clusters.parquet with that `(cluster_id, genome_uid)`.

---

## 7. Parallelism layout — no oversubscription

Three levels of parallelism exist in the pipeline. They must not
multiply each other uncontrollably.

```
Level 1 (stage-outer): ProcessPool over genomes — GenBank parse only
Level 2 (tool-inner):  clustering tool's --threads
Level 3 (post-cluster): NumPy/Polars thread pool
```

Rules:

- Only **one level is active at a time**. The pipeline runs stages
  sequentially, and each stage picks the right parallelism level for
  its workload. No nested pools.
- Global `DNMBCLUSTER_THREADS` env var is the single knob; CLI
  `--threads N` sets it. Stages consult this, not `os.cpu_count()`
  directly.
- `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `POLARS_MAX_THREADS` is
  the one we let scale — prevents NumPy from spawning a thread pool
  inside a Polars worker, which causes the classic 200% CPU thrash.

---

## 8. Caching and resumability — bounded, never unbounded

Cache must never crash a run by filling the disk. Every cache operation
is bounded, observable, and evictable.

### Layout

- Content-addressed cache at `results/.cache/<sha256>/`
- Key = sha256 of (input genome set hash + engine params + tool version
  + DNMBcluster version + SPEED.md hash)
- Stage outputs are **hardlinked** (not symlinked, not copied) into
  `results/` after a successful build. Hardlink = single inode, zero
  extra bytes, survives cache eviction of the source.
- Each cache entry is a directory with `meta.json`:
  ```json
  {"key": "...", "stage": "genbank_parse", "size_bytes": 483922013,
   "created": "...", "last_access": "...", "hit_count": 3}
  ```

### Size cap and eviction

- **Default cap: `min(10 GB, 25% of free disk on results volume)`.**
  Tuneable via `--cache-max 5G` or `DNMBCLUSTER_CACHE_MAX=5G`.
- **LRU eviction** by `last_access`, triggered when: (a) insertion would
  exceed cap, or (b) `dnmbcluster cache clean` invoked manually.
- Entries currently referenced by a running pipeline are locked
  (flock on `entry/.lock`) and never evicted mid-run.
- `dnmbcluster cache stats` — shows total size, per-stage breakdown,
  LRU order.
- `dnmbcluster cache clean [--older-than 7d] [--stage mmseqs2]` —
  manual eviction, scoped.
- `--no-cache` — disable reads and writes entirely for a run.
- `--force` — disable reads (but still write), for stress testing.
- `--resume` (default) — read + write.

### Pre-flight disk check

Before stage start:
```
free = statvfs(output_dir).f_bavail * f_frsize
need = estimate(stage, n_genomes)           # empirical per-stage curves
if free < need * 1.5:
    abort with: "Need ~{need} GB free, have {free} GB on {path}. "
                "Run `dnmbcluster cache clean` or free disk."
```
Per-stage `estimate()` uses linear models fitted from M6 benchmark
runs, stored in `src/dnmbcluster/capacity.py`.

### Intermediate file discipline

- Every stage writes to `$TMPDIR/<stage>/` first, then atomically moves
  into the cache on success. Partial outputs on failure are deleted by
  the signal handler, never orphaned.
- Clustering-tool scratch (`mmseqs tmp`, `diamond tmpdir`, CD-HIT temp)
  goes to `$TMPDIR`, never into the cache directory. Wiped on stage exit
  regardless of outcome.
- Polars `sink_parquet(..., compression="zstd", compression_level=3)`
  for everything in `results/`. Level-3 is the Pareto sweet spot; going
  to level 9 saves ~15% size but costs ~3× CPU.
- For ephemeral intermediates (kept only during one stage), use
  `compression="lz4"` — ~2× faster write, ~1.5× larger. Wins when the
  file is deleted within seconds.
- NEVER write uncompressed Parquet.

---

## 9. R side (visualization)

Mirror the Python rules on the R side:

- **`arrow::read_parquet`** everywhere. Never `read.table`, never
  `read.csv` for files >1 MB.
- **`data.table`** for any tabular manipulation the R code needs.
  `dplyr` acceptable only for plot-layer sugar.
- Plotting is the hot path on the R side: `ggplot2` is fine for most
  plots, but pan/core curves with thousands of points should use
  `geom_line` with pre-aggregated bootstrap summaries (never plot raw
  per-permutation points).
- `ggtree` is slow on large trees — pre-compute tree layouts once,
  cache layout, re-plot from cache.
- Rscript invocation: `--vanilla` to skip Rprofile/Renviron; `env
  R_COMPILE_PKGS=0` to skip JIT.

---

## 10. Docker runtime tuning

- `--tmpfs /tmp:size=${DNMB_TMPFS:-8g},mode=1777` recommended. Document
  how to raise it in README for large runs.
- `--cpus` — if the user caps CPUs, `DNMBCLUSTER_THREADS` must honor
  it. Detect via `/sys/fs/cgroup/cpu.max` not `nproc` (which lies inside
  containers).
- `--memory` cap → propagate to clustering-tool memory flags
  (`--split-memory-limit`, `-M`) as (cap × 0.8).
- `OMP_PROC_BIND=false` — avoids MMseqs2/DIAMOND trying to pin threads
  that conflict with the container scheduler.
- Build-time: `ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1`.

---

## 11. Disk & memory budget enforcement

The pipeline must **fail fast with a clear message** rather than crash
or swap. Every stage declares an expected high-water mark; the runner
checks it against available resources before starting.

### Memory

- `--max-ram 8G` CLI flag. `DNMBCLUSTER_MAX_RAM` env var equivalent.
- Default: 80% of container memory limit (from `/sys/fs/cgroup/memory.max`),
  falling back to 80% of `/proc/meminfo` MemAvailable on bare metal.
- Propagated to every clustering tool:
  - MMseqs2 `--split-memory-limit {max_ram * 0.8}`
  - DIAMOND `-M {max_ram * 0.8}`
  - CD-HIT `-M {max_ram_mb}`
  - usearch12 — no RAM flag; enforced by OS cgroup
- Polars respects `POLARS_STREAMING` for memory-bounded execution when
  a stage operates on something larger than RAM.

### Disk

- Track cumulative bytes written per run in `manifest.json`.
- Every writer calls `budget.check(estimated_bytes)` before writing.
  If the check fails: raise a clean `BudgetExceeded` with a recovery
  hint (`try --cache-max 5G` or `free disk on {path}`).
- Abort cleanly via signal handler (`SIGINT`, `SIGTERM`, `SIGXFSZ`):
  delete partial outputs, unlock cache entries, print trace.

### Python-side memory hygiene

- Streaming I/O by default (`pl.scan_parquet`, `sink_parquet`).
- Never build giant Python lists before `pl.DataFrame(rows)`. Always
  build `pyarrow.RecordBatch`es and concatenate with zero-copy.
- Explicitly `del` large intermediates inside stages; don't rely on GC.
- Run the pipeline with `MALLOC_TRIM_THRESHOLD_=131072` to let glibc
  return freed heap to the OS sooner.
- On Linux, set `MALLOC_ARENA_MAX=2` — NumPy/BLAS default arenas waste
  RSS at scale.
- **jemalloc or mimalloc** as `LD_PRELOAD` in the Docker image for the
  Python process. Benchmark in M6; choose the one that wins on the
  1000-genome profile.

### R-side memory hygiene

- `Rscript --vanilla --max-mem-size={max_ram}` where supported.
- Use `arrow::open_dataset()` for lazy Parquet scans; call `collect()`
  only inside the plot function that needs the data.
- Between plots, `rm(list = ls()); gc(full = TRUE)`.

---

## 12. Sequence level — protein default, nucleotide optional

Comparative genomics in DNMBcluster operates on **CDS features from
GenBank**. The `--level` flag controls what representation is extracted
and clustered:

| `--level` | Extracts | Clustered representation | Default ID threshold |
|---|---|---|---|
| `protein` (default) | CDS `/translation` | amino acid | 0.5 (50%) |
| `nucleotide` | CDS nucleotide (genomic span, reverse-complemented on `-` strand) | DNA | 0.7 (70%) |

Both modes still operate strictly on CDS features — never on `gene`,
`tRNA`, `rRNA`, `ncRNA`, or `misc_feature`. Non-CDS features are
ignored at parse time. (A future flag `--include-rna` could lift this,
but is out of scope for v1.)

### Why a lower protein ID threshold than nucleotide

Protein sequences are more conserved; 50% AA ≈ 70% NT at ortholog
distance in bacteria. These defaults match BPGA / Panaroo / Roary
conventions.

### Engine dispatch

The clustering engine adapter must route nucleotide mode correctly:

| Engine | `protein` | `nucleotide` |
|---|---|---|
| **MMseqs2** | `easy-linclust --min-seq-id 0.5 -c 0.8` | `easy-linclust --min-seq-id 0.7 -c 0.8 --search-type 3` |
| **DIAMOND** | `diamond deepclust` | **NOT SUPPORTED** — DIAMOND is protein-only. CLI errors out with a clean message and suggests MMseqs2 or usearch12 |
| **CD-HIT** | `cd-hit -c 0.5 -n 3` | **`cd-hit-est` binary**, `-c 0.7 -n 5` (different binary, different word-size rules) |
| **usearch12** | `cluster_fast -id 0.5` on protein | `cluster_fast -id 0.7` on nucleotide (same command, alphabet auto-detected) |

Engine adapters inherit a `supports(level) -> bool` method; the CLI
validates `--level` against `--tool` up front and errors out before
any parsing begins. No partial runs.

### Parser impact

- `protein` mode: extract CDS `/translation` only. Skip ORIGIN entirely
  (no genomic sequence in RAM). Pseudogenes with no `/translation` are
  skipped with a `skipped_pseudogene` counter.
- `nucleotide` mode: extract CDS nucleotide span from the LOCUS sequence
  (requires reading ORIGIN), reverse-complement on `-` strand. This is
  ~3× slower than protein mode and uses more RAM; document the
  tradeoff in the CLI help.
- Both modes share the same identifier scheme (Section 1); the only
  difference is the `translation` column vs a `nt_sequence` column in
  `gene_table.parquet`.

### Schema impact

`gene_table.parquet` gets one of these two columns depending on mode
(never both, to avoid bloat):

```
-- --level protein (default)
translation : Utf8

-- --level nucleotide
nt_sequence : Utf8
```

The mode is recorded in `manifest.json` and in the Parquet file
metadata so downstream code can branch on it without inspection.

---

## Open questions (Codex rescue pass will weigh in)

1. Is there a faster GenBank parser than a hand-rolled one? (pyhmmer
   ecosystem? genomepy? biotite?)
2. Should `clusters.parquet` use row-group-level statistics for
   per-cluster predicate pushdown, or flat? (depends on downstream
   access pattern)
3. Roaring bitmaps vs `np.ndarray[uint64]` — crossover point?
4. Is there any gain in offloading the Polars → Arrow → R handoff
   through shared memory (`pyarrow.plasma` / Arrow IPC fd passing)
   instead of on-disk Parquet?
5. For the 1000-genome target: does `diamond deepclust` scale linearly,
   or do we need sharding + MCL post-merge?

These get resolved in M6 benchmarks and this document updated.
