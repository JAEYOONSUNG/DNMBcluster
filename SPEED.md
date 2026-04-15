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

## 1. Compact integer identifier scheme (single biggest lever)

**Rule:** everything past the GenBank parser works on integer IDs, never
strings. Strings exist only in the parser (input) and the final report
layer (output).

### Scheme

```
genome_uid  : uint16   (supports 65 535 genomes per run)
gene_uid    : uint32   (per-genome CDS index, supports 4.3B)
protein_uid : uint64 = (genome_uid << 48) | gene_uid
```

### Why

- Polars/Arrow joins on `int64` are ~3–5× faster than on `Utf8`.
- FASTA headers shrink from ~35 bytes (`lcl|WP_012345678.1|Geobacillus...`)
  to ~15 bytes (`>g3p00001234`). On 3M-protein input that's ~60 MB less for
  the clustering tool to parse.
- `.uc` / `_cluster.tsv` parsing becomes pure integer work.
- Presence/absence matrix keys are native int64 — no string hashing
  per lookup.

### Round-trip mapping

Sidecar `id_map.parquet` columns (written once after GenBank parse):

```
protein_uid:uint64  genome_uid:uint16  gene_uid:uint32
genome_accession:Utf8   (e.g., GCF_000009785.1)
locus_tag:Utf8          (e.g., GK0001)
protein_id:Utf8         (e.g., WP_012345678.1)
product:Utf8
```

All downstream joins that need human-readable labels go through
`id_map.parquet` **once** at report time, never in the hot path.

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

## 8. Caching and resumability

- Content-addressed cache at `results/.cache/<sha256>/`
- Key = sha256 of (input genome set hash + engine params + tool version)
- Stage outputs are symlinked into `results/` after successful build
- `--force` bypasses cache
- `--resume` (default) reuses cache hits

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

- `--tmpfs /tmp:size=8g,mode=1777` always set in recommended `docker
  run`. Document in README.
- `--cpus` — if the user caps CPUs, `DNMBCLUSTER_THREADS` must honor
  it. Detect via `/sys/fs/cgroup/cpu.max` not `nproc` (which lies inside
  containers).
- `OMP_PROC_BIND=false` — avoids MMseqs2/DIAMOND trying to pin threads
  that conflict with the container scheduler.
- Build-time: `ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1`.

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
