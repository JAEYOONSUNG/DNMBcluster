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
   (member→seed + seed→member) and **both-direction** coverage.
5. Computes presence/absence bitmaps, pan/core curves, and per-cluster
   summaries.
6. Emits a Roary-compatible `gene_presence_absence.csv` for the whole
   Roary/Scoary/Phandango ecosystem.
7. Generates flower, pan-core, and category plots via an internal R
   package (absorbed from [BPGAconverter](https://github.com/JAEYOONSUNG/BPGAconverter)).

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
| `protein` (default) | **0.5** (50%) | 0.8 (80%) | BPGA / Roary / Panaroo convention for bacterial pan-genome work; captures orthologs across species boundaries |
| `nucleotide` | **0.7** (70%) | 0.8 (80%) | Nucleotide sequences are less conserved than protein translations at ortholog distance; 50% AA ≈ 70% NT in bacteria |

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
  --identity FLOAT            Identity threshold 0.0-1.0 (default: 0.5 protein / 0.7 nucleotide)
  --coverage FLOAT            Alignment coverage threshold (default: 0.8)
  --threads INT               Thread count, 0=auto-detect (default: 0)
  --max-ram SIZE              Memory cap, e.g. "8G" (default: container limit × 0.8)
  --alignment / --fast        MMseqs2 bidirectional alignment pass (default: --alignment)
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

## Output structure

```
results/
├── dnmb/                              # DNMB-native Parquet (fastest I/O)
│   ├── id_map.parquet                 # canonical identifiers, 1 row per CDS
│   ├── gene_table.parquet             # protein/DNA sequences keyed by protein_uid
│   ├── genome_meta.parquet            # 1 row per input file (sha256, GC%, n_cds, ...)
│   ├── clusters.parquet               # 1 row per CDS: cluster_id, is_centroid,
│   │                                  #   pct_identity_fwd/_rev, member_coverage,
│   │                                  #   rep_coverage, alignment_length
│   ├── presence_absence.parquet       # 1 row per cluster, uint64[] genome bitmap
│   ├── pan_core_curve.parquet         # N_permutations × N_genomes rows
│   ├── cluster_summary.parquet        # 1 row per cluster: category, rep metadata
│   ├── proteins.faa (or cds.fna)      # integer-header FASTA fed to the engine
│   ├── mmseqs_out_*.*                 # raw engine intermediates (can be deleted)
│   └── mmseqs_alignments_{fwd,rev}.tsv
├── gene_presence_absence.csv          # Roary-compatible; Scoary/Phandango direct input
├── gene_presence_absence.Rtab         # Roary binary 0/1 companion
└── plots/
    ├── flower_plot.pdf                # core/accessory/unique petal plot
    ├── pan_core_plot.pdf              # pan vs core genome curves
    └── category_bar.pdf               # stacked cluster categories per genome
```

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
│       ├── pan_core_plot.R
│       ├── category_bar.R
│       └── run_dnmb_plot.R # entry point called from Python
├── tests/
│   ├── fixtures/
│   │   └── geobacillus10/  # 10-genome benchmark GenBank set
│   └── test_*.py           # 48 unit + integration tests
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
