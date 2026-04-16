# DNMBcluster

Dockerized pan-genome clustering pipeline with integrated visualization.
A modern, Mac/Linux-native, single-command replacement for BPGA
(Bacterial Pan Genome Analysis). Takes a folder of GenBank files,
returns Parquet dataframes, a Roary-compatible CSV, and publication
plots.

**Status**: M0–M7 complete. Verified end-to-end on 10 *Geobacillus*
genomes (33,116 proteins → 6,508 clusters) in 65 s including plots.

## What it does

1. Walks a folder of `.gb` / `.gbk` / `.gbff` files.
2. Extracts CDS features (translations by default, nucleotide with
   `--level nucleotide`).
3. Clusters with your choice of engine: **MMseqs2** (default), **DIAMOND
   deepclust**, **CD-HIT**, or **usearch12**.
4. Produces a unified `clusters.parquet` with **bidirectional** identity
   (member→seed + seed→member) and **both-direction** coverage, plus
   an engine-native sidecar (CIGAR for usearch12, alignment positions
   for CD-HIT `-p 1`) when the engine emits extras.
5. Computes presence/absence bitmaps, pan/core curves, per-cluster
   summaries, and a long-format `cluster_long.parquet` that joins every
   CDS to its cluster + category + alignment metrics in one file.
6. Emits a Roary-compatible `gene_presence_absence.csv` for the whole
   Roary/Scoary/Phandango ecosystem, plus a BPGA-style single-sheet
   `comparative_genomics_N.xlsx` (merged-view with DEFINITION/accession
   headers, vertical separators every 3 columns, centered cells).
7. Generates **17 publication-quality plots** via an internal R package
   (absorbed from [BPGAconverter](https://github.com/JAEYOONSUNG/BPGAconverter)).
   See the [Visualizations](#visualizations) section below.
8. **Optional phylogenomics** (`--phylo`): single-copy core extraction
   → MAFFT alignment → optional trimAl → concatenation → IQ-TREE best
   tree with UFBoot + SH-aLRT support → ggtree visualization.
9. **Optional ANI + POCP similarity** (`--ani`): skani triangle ANI
   matrix + cluster-based POCP, both rendered as pheatmap heatmaps
   with hierarchical clustering.
10. **Context-ribbon subcommand** (`dnmbcluster context`): pick a
    genome + locus_tag and get a stacked gggenes neighborhood plot
    with shared-cluster colors showing ortholog context across every
    strain that carries the anchor cluster.

## Install

### Docker (recommended)

```bash
docker pull ghcr.io/jaeyoonsung/dnmbcluster:latest
```

### From source (local dev)

```bash
git clone https://github.com/JAEYOONSUNG/DNMBcluster.git
cd DNMBcluster

# Create a conda env with all four clustering engines + Python + R
mamba env create -n dnmb -f env.yml

# Install the Python package in editable mode
~/miniforge3/envs/dnmb/bin/pip install -e .

# Install the R visualization package
~/miniforge3/envs/dnmb/bin/R CMD INSTALL R
```

`usearch12` is not on bioconda — the Docker image compiles it from
[rcedgar/usearch12](https://github.com/rcedgar/usearch12) at build
time. For local dev without Docker, build it yourself and put the
binary on `$PATH` as `usearch12`.

## Quick start

Put your GenBank files in a folder, then run one command:

```bash
docker run --rm \
  -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) \
  --tmpfs /tmp:size=8g,mode=1777 \
  -v $PWD:/data \
  ghcr.io/jaeyoonsung/dnmbcluster:latest \
  run /data/my_genbanks
```

Same command, local install:

```bash
dnmbcluster run ./my_genbanks
```

Outputs land in `./results/` by default.

## The identity threshold — the most important parameter

`--identity` is the single knob that changes everything downstream —
cluster count, core gene set, pan/core curves, flower plot values, and
the Roary CSV. Pick it carefully for your biological question.

### Defaults

| `--level` | Default `--identity` | Default `--coverage` | Reasoning |
|---|---:|---:|---|
| `protein` (default) | **0.6** (60%) | 0.8 (80%) | Mid-point suitable for species-level pan-genome work across all supported engines; balanced between ortholog capture and paralog separation |
| `nucleotide` | **0.6** (60%) | 0.8 (80%) | Same single default so engine/level swaps don't silently change the clustering threshold; tighten to 0.7–0.8 for within-species |

### When to change it

| Scenario | Suggested `--identity` | Why |
|---|---:|---|
| Within-species core genome (e.g. same species, different strains) | 0.90 – 0.95 | Tighter grouping; splits out allelic variants |
| Species-level pan-genome (default BPGA use case) | **0.50 – 0.70** | Broad enough to capture orthologs, tight enough to separate paralog families |
| Genus-level comparison (mixed species) | 0.40 – 0.50 | Needed for cross-species orthology |
| Family-level or deeper | 0.30 – 0.40 | Below this, clustering tools lose sensitivity — use MMseqs2 or DIAMOND, not CD-HIT |

### CD-HIT word-size rule

CD-HIT imposes a hard mapping from `--identity` to the internal word
size (`-n`):

| Protein identity | CD-HIT word size |
|---|---:|
| ≥ 0.70 | 5 |
| 0.60 – 0.70 | 4 |
| 0.50 – 0.60 | 3 |
| 0.40 – 0.50 | 2 (unstable; prefer MMseqs2 below 0.5) |

DNMBcluster picks the word size automatically from your `--identity`,
but if you go below 0.4 on CD-HIT the tool becomes unreliable — switch
engines.

### How `--coverage` interacts with `--identity`

`--coverage` is the fraction of the shorter sequence that must be
aligned for the pair to be clustered. It prevents long divergent
proteins from being merged just because a short domain happens to
match. Default **0.8** (80%) is the right choice for almost every
pan-genome use case. Lower it to 0.5 only if you explicitly want
domain-level families.

### The bidirectional alignment pass

Running `mmseqs2` without `--fast` triggers two separate alignment
passes:

1. **Forward** — each member aligned as query against its cluster
   representative as target. Gives `pct_identity_fwd`,
   `member_coverage`, `rep_coverage`, `alignment_length`.
2. **Reverse** — each representative aligned as query against all its
   members as targets. Gives `pct_identity_rev`.

For length-disparate sequences the heuristic aligner can produce
slightly different alignments depending on which sequence is the
query, so reporting both directions catches asymmetries that a
single-direction score would hide. If you only need cluster
membership, pass `--fast` to skip both alignment passes — the extra
columns land as `NULL` but wall time drops ~6.8× on the Geobacillus-10
benchmark.

## Full CLI reference

```
dnmbcluster run INPUT_DIR [OPTIONS]

Arguments:
  INPUT_DIR                   Folder containing *.gb / *.gbk / *.gbff files

Options:
  -o, --output PATH           Output directory (default: ./results)
  --tool [mmseqs2|diamond|cd-hit|usearch12]
                              Clustering engine (default: mmseqs2)
  --level [protein|nucleotide]
                              Sequence level (default: protein)
  --identity FLOAT            Identity threshold 0.0-1.0 (default: 0.6 for both levels)
  --coverage FLOAT            Alignment coverage threshold (default: 0.8)
  --threads INT               Thread count, 0=auto-detect (default: 0)
  --max-ram SIZE              Memory cap, e.g. "8G" (default: container limit × 0.8)
  --alignment / --fast        Pipeline-wide bidirectional alignment pass (default: --alignment)
  --columns TEXT              Comma-separated middle block columns in comparative_genomics_N.xlsx
                              (default: "product"; see Comparative Excel section)
  --parse-only                Stop after GenBank parsing
  -h, --help                  Show this message and exit
```

## Examples

**Default species-level pan-genome (10 Geobacillus genomes):**
```bash
dnmbcluster run ./geobacillus10
```
Runs with `mmseqs2`, protein level, identity 0.5, coverage 0.8,
bidirectional alignment on. Outputs in `./results/`.

**Maximum speed, same clustering:**
```bash
dnmbcluster run ./geobacillus10 --fast
```
Drops the alignment pass for ~6.8× speedup. Cluster membership is
identical; only the identity/coverage metric columns are null.

**Strain-level (tight threshold):**
```bash
dnmbcluster run ./strains --identity 0.95
```

**Nucleotide clustering of CDS DNA:**
```bash
dnmbcluster run ./genomes --level nucleotide --identity 0.7
```

**Explicit engine choice + thread cap + memory cap:**
```bash
dnmbcluster run ./genomes \
    --tool diamond \
    --identity 0.5 --coverage 0.8 \
    --threads 16 --max-ram 32G
```

**Parse only (debug GenBank ingestion):**
```bash
dnmbcluster run ./genomes --parse-only -o /tmp/parse_check
```

**Legacy CD-HIT path (for BPGA-style reproduction):**
```bash
dnmbcluster run ./genomes --tool cd-hit --identity 0.5
```

## Comparative genomics Excel — flexible column layout

The `comparative_genomics_N.xlsx` merged sheet gives every genome its
own side-by-side block of columns. Two columns per block are **fixed**:

- **First column of every block**: `locus_tag` (falls back to `cds_key`
  if the CDS has no locus_tag in the GenBank source)
- **Last column of every block**: `identity (%)` (= `pct_identity_fwd`
  from the realignment stage)

Everything between them is controlled by the `--columns` CLI flag.
The default produces a 3-column block `locus_tag | product | identity (%)`,
matching the BPGA triple. Add more attributes by listing them comma-
separated:

```bash
dnmbcluster run ./genomes \
    --columns product,gene,protein_id,pct_identity_rev
```

yields a 6-column block per genome:

```
locus_tag | product | gene | protein_id | identity_rev (%) | identity (%)
```

Empty string produces the minimal 2-column block `locus_tag | identity (%)`.

### Selectable middle columns

| Column | Source | Type | Agg rule |
|---|---|---|---|
| `product` | id_map (default) | string | first non-null |
| `gene` | id_map | string | first non-null |
| `protein_id` | id_map | string | first non-null |
| `ec_number` | id_map | string | first non-null |
| `contig` | id_map | string | first non-null |
| `aa_length` | id_map | uint32 | first non-null |
| `pct_identity_rev` | clusters | float32 (0-100) | first non-null, 2dp |
| `member_coverage` | clusters | float32 (0-1) | first non-null, 2dp |
| `rep_coverage` | clusters | float32 (0-1) | first non-null, 2dp |
| `alignment_length` | clusters | uint32 | first non-null |

Unknown column names, duplicates, and the two fixed columns are
rejected at the CLI layer before any data is loaded, so typos fail
fast with a clear error message.

Column widths, borders, and numeric formats all adjust automatically
to the block width. The DEFINITION + accession two-line merged header
still spans the entire block regardless of how many middle columns
you added.

## Output structure

The output directory is split into three zones: `inputs/` (engine-
agnostic genome parsing artifacts, shared across reruns), `raw/`
(per-engine native output untouched), and `processed/` (canonical
Parquet every downstream consumer joins on). Top-level files are the
user-facing reports.

```
results/
├── dnmb/
│   ├── inputs/                        # Engine-agnostic, rerun-stable
│   │   ├── id_map.parquet             # canonical identifiers, 1 row per CDS
│   │   ├── gene_table.parquet         # protein/DNA sequences keyed by protein_uid
│   │   ├── genome_meta.parquet        # 1 row per input file (sha256, GC%, n_cds,
│   │   │                              #   organism, strain, DEFINITION, ...)
│   │   └── proteins.faa (or cds.fna)  # integer-header FASTA fed to the engine
│   ├── raw/                           # Untouched engine-native output
│   │   ├── cd-hit/ | diamond/ | mmseqs2/ | usearch12/
│   │   │   └── (engine's own text / TSV / clstr / uc files)
│   │   └── realign/                   # Shared MMseqs2 realignment stage
│   │       ├── realign_fwd.tsv
│   │       └── realign_rev.tsv
│   └── processed/                     # Canonical Parquet dataframes
│       ├── clusters.parquet           # 1 row per CDS: cluster_id, is_centroid,
│       │                              #   pct_identity_fwd/_rev, member_coverage,
│       │                              #   rep_coverage, alignment_length
│       ├── cdhit_native.parquet       # (CD-HIT only) per-member alignment
│       │                              #   positions from -p 1
│       ├── usearch12_native.parquet   # (usearch12 only) CIGAR + strand
│       ├── presence_absence.parquet   # 1 row per cluster, uint64[] genome bitmap
│       ├── pan_core_curve.parquet     # N_permutations × N_genomes rows
│       ├── cluster_summary.parquet    # 1 row per cluster: category, rep metadata
│       └── cluster_long.parquet       # 1 row per CDS: everything joined
│                                      #   (cluster_id, category, genome_key,
│                                      #   locus_tag, gene, product, is_centroid,
│                                      #   pct_identity_fwd/_rev, coverages, aln_len)
├── comparative_genomics_N.xlsx        # BPGA-style merged sheet, 1 row per cluster
├── gene_presence_absence.csv          # Roary-compatible; Scoary/Phandango direct input
├── gene_presence_absence.Rtab         # Roary binary 0/1 companion
└── plots/                             # 13 R-rendered publication PDFs
    ├── flower_plot.pdf
    ├── pan_core_plot.pdf              # + Heaps' law gamma in subtitle
    ├── category_bar.pdf
    ├── presence_absence_heatmap.pdf
    ├── cluster_size_dist.pdf
    ├── identity_distribution.pdf
    ├── genome_jaccard_heatmap.pdf     # dual-metric, two color scales
    ├── coverage_scatter.pdf
    ├── gene_content_mds.pdf
    ├── centroid_length_distribution.pdf
    ├── core_gene_conservation.pdf
    ├── cumulative_pan_contribution.pdf
    └── singleton_top_n.pdf
```

## Visualizations

Default pipeline runs render **13 plots** from `dnmb/processed/` +
`dnmb/inputs/`. Opt-in stages (`--phylo`, `--ani`) add up to **4 more**.
The standalone `dnmbcluster context` subcommand produces per-anchor
ortholog ribbon plots on demand.

| Category | Plots |
|---|---|
| **Overview** | `flower_plot`, `category_bar`, `cluster_size_dist`, `ortho_euler` |
| **Rarefaction** | `pan_core_plot` (with Heaps' law γ in subtitle) |
| **Matrix / relatedness** | `presence_absence_heatmap`, `genome_jaccard_heatmap`, `gene_content_mds` |
| **Per-genome contribution** | `cumulative_pan_contribution`, `singleton_top_n` |
| **QC** | `identity_distribution`, `coverage_scatter` |
| **Biology** | `centroid_length_distribution`, `core_gene_conservation` |
| **Opt-in — ANI/POCP** (`--ani`) | `ani_heatmap`, `pocp_heatmap` |
| **Opt-in — phylogenomics** (`--phylo`) | `phylo_tree` (IQ-TREE best tree via ggtree) |
| **Standalone subcommand** | `context_ribbon` (gggenes, per-anchor) |

The R layer is also importable directly — `library(DNMBcluster);
load_dnmb("results")` returns the Parquet dataframes as tibbles if you
want to script your own analyses on top.

### Phylogenomics (`--phylo`)

Runs only when the user passes `--phylo`. Under the hood:

1. Select single-copy core clusters — category `core` **and** exactly
   one member per genome (multi-copy core is skipped).
2. Extract per-cluster FASTAs from `gene_table.parquet`.
3. Align each cluster with MAFFT `--auto` in a thread pool; optional
   trimAl `-automated1` trimming when the binary is available.
4. Concatenate aligned clusters into a single supermatrix.
5. IQ-TREE with a fixed model (`LG+G4` protein / `GTR+G4` nucleotide)
   + UFBoot 1000 + SH-aLRT 1000. ModelFinder is intentionally skipped
   to keep wall time under 5 min on 10-genome runs.
6. Newick tree lands at `dnmb/processed/phylo_tree.nwk`; ggtree renders
   `plots/phylo_tree.pdf` with support-colored internal nodes.

Scales: ~3 min on Geobacillus-10 with 300 core clusters capped via the
internal `max_clusters=300` default; remove the cap by editing the
call in `cli.py` if you want the full core set.

### ANI + POCP (`--ani`)

Runs only when `--ani` is set. Under the hood:

1. Lazily extract one contig FASTA per genome into
   `dnmb/raw/genomes/{genome_key}.fna` (cached; only touched if missing).
2. Run `skani triangle --full-matrix` once across all genomes → raw
   TSV at `dnmb/raw/ani/skani_triangle.tsv`.
3. Parse the matrix into long-form `dnmb/processed/ani_matrix.parquet`.
4. Compute POCP directly from `clusters.parquet` + `genome_meta.n_cds`
   using the Qin et al. 2014 ratio. Our approximation replaces the
   classical BLASTP 40%/50% threshold with "members share a cluster",
   so the POCP value reflects whatever `--identity` the clustering
   run used — document this in your methods.
5. `plots/ani_heatmap.pdf` and `plots/pocp_heatmap.pdf` render via
   pheatmap with average-linkage hierarchical clustering on
   `1 - metric/100`. Species boundary at ANI ≥ 95%, genus boundary at
   POCP ≥ 60% are called out in the title.

### Context ribbon (`dnmbcluster context`)

Standalone subcommand that works against an existing results
directory. Example:

```bash
dnmbcluster context ./results \
    --genome GCF_000009785.1 \
    --locus GK_RS07590 \
    --flank 5
```

Finds the anchor CDS, walks every strain that shares the anchor's
cluster, extracts the ±`flank` flanking CDSs on each strain's matching
contig, and renders a stacked gggenes plot. Shared-cluster genes
across strains get the same fill color (your visual ribbon), the
anchor homolog in every strain is marked with a ★, and strain-unique
genes fade to grey. Writes to
`<results>/plots/context_<genome>_<locus>.pdf` by default.

### Identifier scheme

DNMBcluster uses a **two-level composite primary key** for every CDS:

```
{genome_key}::{cds_key}     e.g.  "GCF_000009785.1::GK_RS07590"
```

The `genome_key` comes from a 6-step fallback chain (DBLINK Assembly
accession → filename `GCF_*` → LOCUS/VERSION → organism slug →
filename stem) and the `cds_key` from a 4-step chain (`locus_tag` →
`protein_id` → `{gene}__{ordinal}` → coordinate hash). Both fallback
rungs are recorded in `genome_key_source` / `cds_key_source` columns so
you can see exactly how a label was derived. RefSeq WP_ protein
accessions are deliberately NOT used as the primary key because they
are shared across closely related strains by design.

A compact `protein_uid:uint64` (run-scoped) is used as the hot-path
join key. The stable primary key is the `(genome_key, cds_key)` pair.
See `SPEED.md` Section 1 for the full design.

## Analyzable scope

What the pipeline is designed to do, and what it's not.

### Input scale

| Parameter | Supported range | Hard cap | Notes |
|---|---|---|---|
| Genomes per run | 2 – 65,535 | `uint16 genome_uid` | Single-word genome bitmap up to 64; multi-word beyond, no code-level limit until 65,535 |
| CDS per genome | 1 – 4,294,967,295 | `uint32 gene_uid` | Packed into the lower 32 bits of `protein_uid` |
| Total CDS across run | ~10M benchmarked | unbounded in schema | Wall time grows roughly linearly; memory grows with clusters × genomes |
| Sequence level | `protein` (default), `nucleotide` | — | DIAMOND is protein-only; `--level nucleotide --tool diamond` is rejected up-front |
| Identity threshold | 0.0 – 1.0 | engine-specific (see below) | MMseqs2 / DIAMOND effective range ≥ 0.3 ; CD-HIT ≥ 0.4 (hard word-size floor); usearch12 ≥ 0.3 |
| Coverage threshold | 0.0 – 1.0 | — | Applied as `--cov-mode 0` (coverage of shorter sequence) |

### What each engine covers

| Engine | protein | nucleotide | Notes |
|---|:---:|:---:|---|
| **MMseqs2** (default) | yes | yes | Fastest via linclust; default choice for bacterial pan-genomes |
| **DIAMOND deepclust** | yes | no | Up to 82× CD-HIT on large protein sets; rejects `--level nucleotide` |
| **CD-HIT** | yes (`cd-hit`) | yes (`cd-hit-est`) | Dispatches to the right binary by level; identity minimum 0.4 (protein) |
| **usearch12** | yes | yes | Source-only build; bundled in Docker, manual install for local dev |

### What you get out (post-run)

- **11 core Parquet dataframes** under `dnmb/inputs/` and `dnmb/processed/`, all schema-validated at stage boundaries.
- **1 BPGA-style Excel** with merged-cell DEFINITION headers + centered cells + 3-column block borders.
- **2 Roary-compatible files** for direct Scoary / Phandango / panGOLIN consumption.
- **13 publication PDFs** covering overview, rarefaction, relatedness, per-genome contribution, QC, and biology (list above).
- **Uniform schema across all engines** — `clusters.parquet` columns and types are identical regardless of `--tool`, so downstream code never branches on the clustering backend.

### Tested engine × identity matrix

The pipeline is end-to-end tested on the Geobacillus-10 fixture
across every combination of:

```
engines    : mmseqs2, diamond, cd-hit
identities : 0.3, 0.5, 0.7, 0.9   (cd-hit skips 0.3 per its 0.4 floor)
```

All 11 runs pass structural, schema, file-presence, validation, and
plotting checks. Cluster counts scale monotonically with `--identity`
within each engine; cross-engine cluster count spreads fall within
the expected algorithmic-variation range.

## Caveats & known limits

Things that will bite if you don't know about them.

### Pipeline-wide

- **MMseqs2 is a hard runtime dependency even for non-MMseqs2 engines**
  when `--alignment` (default) is on. The shared
  `engines/realign.py` stage uses MMseqs2 `easy-search` to populate
  the bidirectional alignment metric columns uniformly across every
  engine. Pass `--fast` to skip realignment entirely (alignment
  columns land null; cluster membership is unaffected).
- **Non-centroid members in clusters ≥ 2 must have both
  `pct_identity_fwd` and `pct_identity_rev` populated after
  realignment.** The shared validator enforces this as a hard
  assertion and fails the run loudly if any row slips through. A
  short-peptide Biopython fallback covers the rare pairs MMseqs2's
  k-mer prefilter can't seed.
- **`genome_uid` is 16-bit.** Single run ceiling is 65,535 genomes.
  Larger collections need to be split.
- **`protein_uid` is run-scoped**, not stable across runs. Use the
  `(genome_key, cds_key)` composite from `id_map.parquet` as the
  stable cross-run key.
- **`cluster_id` is dense and run-scoped**, re-emitted fresh per
  engine. Don't compare `cluster_id` values across different runs or
  different `--tool` choices.

### Identifiers

- **Duplicate `genome_key`s crash early.** If two input files resolve
  to the same composite key (commonly: two copies of the same
  assembly under different filenames), `build_tables` raises before
  clustering starts — disambiguate via a manifest or rename the
  files.
- **`locus_tag` is NOT unique across an assembly for all Prokka
  reruns.** The parser catches duplicates in a single file and
  appends `__dupN` to keep identifiers unique, emitting a warning.
- **RefSeq `WP_` protein accessions are intentionally not used as a
  primary key** — RefSeq shares them across related strains by
  design, so using them would collapse orthologs into the same row.

### Engine-specific

- **CD-HIT identity minimum is 0.4** (protein). Below that, the
  word-size heuristic becomes unstable — the pipeline keeps passing
  the threshold through but CD-HIT itself will refuse or produce
  unreliable clusters. Use MMseqs2 or DIAMOND for `--identity < 0.4`.
- **DIAMOND is protein-only.** `--tool diamond --level nucleotide`
  is rejected at the CLI layer before parsing starts.
- **DIAMOND has no engine-native sidecar.** `deepclust` emits a
  bare 2-column TSV and its `--aln-out` mode is unreliable in the
  2.1.x series, so `diamond_native.parquet` is intentionally not
  produced. All alignment metrics come from the shared realign
  stage.
- **usearch12 is under GPL-3** and has no bioconda package. The
  Docker image compiles it from source at build time; for local dev
  you need to build it yourself and put the binary on `$PATH` as
  `usearch12`.

### Downstream / reporting

- **Category definitions are strict 3-tier**: `core` =
  present in all genomes, `unique` = present in exactly one,
  `accessory` = everything in between. The earlier Roary-style
  core/soft_core/shell/cloud split has been collapsed — if your
  downstream analysis depends on a 4-tier breakdown, reclassify off
  `n_genomes` in `presence_absence.parquet`.
- **R plotting requires Rscript on `$PATH`** (or `$DNMBCLUSTER_RSCRIPT`).
  In a conda setup, the R package must be installed into the same
  env as the python CLI via `R CMD INSTALL R`; otherwise the Python
  pipeline will still succeed but the plot stage will warn and skip.
- **Heaps' law γ needs at least 2 data points on the pan curve.**
  A single-genome run produces no meaningful openness estimate and
  the subtitle falls back to empty.
- **`gene_content_mds.pdf` needs ≥ 3 genomes.** Smaller runs emit a
  warning and skip that specific plot.

## Benchmark snapshot

10 *Geobacillus* genomes (33,116 proteins), Apple Silicon, 8 threads,
identity 0.5, coverage 0.8:

| Engine | Wall | RSS | Clusters |
|---|---:|---:|---:|
| `mmseqs2 --fast` | **9.9 s** | 329 MB | 6508 |
| `diamond` | 13.6 s | 324 MB | 5709 |
| `mmseqs2` default (bidirectional aln) | 67.4 s | 1004 MB | 6508 |
| `cd-hit` | 70.6 s | 323 MB | 5862 |

Full methodology and per-engine analysis in
[BENCHMARK.md](BENCHMARK.md). `usearch12` is bundled in the Docker
image but not benchmarked in the local conda path (it requires a
source compile from `rcedgar/usearch12` under GPL-3).

## Roadmap

- [x] **M0** repo scaffold
- [x] **M1** GenBank parser + parallel protein extraction
- [x] **M2** MMseqs2 engine + unified cluster TSV
- [x] **M3** DIAMOND, CD-HIT, usearch12 engines
- [x] **M4** DNMB-native Parquet dataframes + Roary-compatible CSV
- [x] **M5** R visualization layer (absorbed from BPGAconverter)
- [x] **M6** Geobacillus-10 benchmark — see [BENCHMARK.md](BENCHMARK.md)
- [x] **M7** Multi-arch CI (amd64/arm64) + GHCR release workflow

## Repository layout

```
DNMBcluster/
├── src/dnmbcluster/        # Python pipeline (CLI, parser, engines, dataframes)
│   ├── cli.py              # Click entry point
│   ├── genbank.py          # Biopython GenBank parser
│   ├── ids.py              # genome_key + cds_key fallback chains
│   ├── schemas.py          # PyArrow schemas for every stage
│   ├── fasta.py            # integer-header FASTA writer
│   ├── engines/            # mmseqs2, diamond, cdhit, usearch12
│   ├── matrix.py           # presence/absence bitmap builder
│   ├── pancore.py          # numpy-vectorized pan/core curves
│   ├── summary.py          # per-cluster summary
│   ├── roary_export.py     # gene_presence_absence.csv / .Rtab
│   └── r_bridge.py         # subprocess bridge to the R package
├── R/                      # DNMBcluster R package (visualization)
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   └── R/
│       ├── load_dnmb.R     # arrow::read_parquet loader
│       ├── flower_plot.R
│       ├── pan_core_plot.R             # + Heaps' law gamma
│       ├── category_bar.R
│       ├── presence_absence_heatmap.R
│       ├── cluster_size_dist.R
│       ├── identity_distribution.R
│       ├── genome_jaccard_heatmap.R    # dual-metric (ggnewscale)
│       ├── coverage_scatter.R
│       ├── gene_content_mds.R
│       ├── centroid_length_distribution.R
│       ├── core_gene_conservation.R
│       ├── cumulative_pan_contribution.R
│       ├── singleton_top_n.R
│       └── run_dnmb_plot.R # entry point called from Python
├── tests/
│   ├── fixtures/
│   │   └── geobacillus10/  # 10-genome benchmark GenBank set
│   └── test_*.py           # 85 unit + integration tests
├── Dockerfile              # multi-stage: usearch12 builder + conda runtime
├── docker-entrypoint.sh    # gosu + tini, HOST_UID auto-detect
├── env.yml                 # conda env spec
├── pyproject.toml          # Python package spec
├── SPEED.md                # speed + integrity design doc (load-bearing)
├── BENCHMARK.md            # Geobacillus-10 timing results
└── .github/workflows/      # ci.yml (pytest + buildx), release.yml (GHCR)
```

## License

MIT for DNMBcluster itself. Bundled third-party tools retain their
original licenses:

- **usearch12** — GPL-3.0 (compiled in Docker image from source)
- **DIAMOND** — GPL-3.0
- **CD-HIT** — GPL-2.0
- **MMseqs2** — MIT
- **Biopython** — Biopython License (BSD-like)
- **R packages (arrow, ggplot2, dplyr, tidyr, tibble, magrittr, scales)** — MIT / GPL-2+

Invoking GPL tools as subprocesses from DNMBcluster (MIT) is mere
aggregation and does not affect the pipeline code's license.
