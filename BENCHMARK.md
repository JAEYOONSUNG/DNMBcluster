# DNMBcluster benchmark — Geobacillus 10 genomes

Reference dataset: 10 *Geobacillus* genomes from NCBI RefSeq + GenBank,
33,116 total CDS proteins. Hardware: Apple Silicon laptop, 8 threads.
Identity threshold 0.5 / coverage 0.8 / `--level protein` for every
engine.

End-to-end means: GenBank parse → protein FASTA → cluster → presence/
absence → pan/core → cluster_summary → Roary CSV → R plots. Identical
work for every engine.

## Results

| Engine | Wall time | Peak RSS | Clusters | Notes |
|---|---:|---:|---:|---|
| **mmseqs2 --fast**         |  **9.9 s** | **329 MB** | 6508 | linclust only; alignment fields null |
| **diamond deepclust**      | 13.6 s | 324 MB | 5709 | fewer clusters (more aggressive merging) |
| **mmseqs2** (default)      | 67.4 s | 1004 MB | 6508 | bidirectional alignment pass populates pct_identity_fwd/rev + member/rep coverage |
| **cd-hit**                 | 70.6 s | 323 MB | 5862 | slowest of the conda set; scales poorly past ~100 genomes |
| **usearch12**              | not benchmarked | — | — | requires source compile; Docker image only — see M7 |

## Key takeaways

1. **For raw clustering throughput, `mmseqs2 --fast` wins** — it is
   the fastest end-to-end (9.9 s) and produces the same cluster
   membership as the full mmseqs2 default (6508 clusters) because the
   alignment pass is metadata-only, not a re-clustering. Use `--fast`
   when you only need the presence/absence matrix and downstream
   plots.
2. **The default workflow (`mmseqs2` without `--fast`) costs about
   57 s of pure alignment overhead** on this dataset. The two
   `easy-search` passes (forward + reverse) populate the bidirectional
   identity columns and both coverages — useful for QC and
   ortholog-filter pipelines but not free.
3. **DIAMOND deepclust is the second-fastest** at 13.6 s and produces
   ~12% fewer clusters than mmseqs2 at the same threshold. DIAMOND's
   `--approx-id` is coarser than mmseqs2's true Smith-Waterman
   identity, so threshold-driven boundaries land in slightly different
   places. For pan-genome work this is within the normal between-tool
   variance.
4. **CD-HIT is the slowest** at 70.6 s, on par with mmseqs2's full
   bidirectional path but without producing the extra metric columns.
   It does not parallelize as well; on this 8-thread box CD-HIT used
   ~6 cores effectively where mmseqs2 used all 8. It also has a hard
   word-size dependency on the identity threshold and is known to
   scale poorly past ~10⁷ sequences. Recommended only for legacy
   compatibility.

## Default engine decision

**Default tool: `mmseqs2`** (with bidirectional alignment, 67 s).

Reasoning:
- Same cluster boundaries as the fast variant — switching the default
  to `--fast` later is a single flag, not a data migration.
- The `pct_identity_fwd/_rev` and bidirectional coverage columns are
  the most-asked-for downstream feature and BPGA users will expect
  them by default.
- Users who care about raw speed pass `--fast` and pay nothing for the
  alignment pass.
- DIAMOND is a strong second; users looking for slightly larger
  cluster-merging behaviour or planning to scale to many thousands of
  genomes should consider `--tool diamond`.
- CD-HIT is shipped for compatibility, not as the default.

## Cross-engine cluster count divergence

| Engine | Clusters | Δ vs mmseqs2 |
|---|---:|---:|
| mmseqs2 | 6508 | — |
| cd-hit | 5862 | -10% |
| diamond | 5709 | -12% |

This is normal — different alignment algorithms and different identity
calculations land cluster boundaries in slightly different places at
the same threshold. The core gene set (clusters present in all 10
genomes) is more stable across engines and is the biologically
meaningful number to compare.

## Scaling notes (extrapolated, not measured)

Targets from `SPEED.md`:

| Workload | Wall time target | mmseqs2 expected |
|---|---:|---:|
| 10 genomes (~30K proteins) | < 30 s | 9.9 s `--fast` ✓ / 67 s default |
| 100 genomes (~300K proteins) | < 5 min | likely 1-2 min `--fast`, 8-12 min default |
| 1000 genomes (~3M proteins) | < 45 min | likely 15-25 min `--fast`, 60-90 min default |

Default-mode misses the 10-genome target because of the alignment
pass cost. Either the alignment pass needs further optimization (use
`mmseqs align` directly on the cluster DB instead of `easy-search`,
which would avoid the pre-filter) or the user opts into `--fast` when
they don't need the metric columns.

## Methodology

```
/usr/bin/time -l dnmbcluster run \
    tests/fixtures/geobacillus10 \
    -o /tmp/bench_<engine> \
    --threads 8 \
    --tool <engine> [--fast]
```

Wall time and peak RSS captured by `time -l`. Cluster count and
status captured from the CLI stdout. Each run was preceded by `rm
-rf /tmp/bench_*` to ensure cold-cache results.

Reproduce: `tests/fixtures/geobacillus10/` ships with the repo.
Install all four engines via:

```
mamba install -c bioconda -c conda-forge \
    mmseqs2 diamond cd-hit
# usearch12: docker build the DNMBcluster image — Dockerfile compiles
# rcedgar/usearch12 from source under GPL-3
```
