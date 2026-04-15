# DNMBcluster

Dockerized pan-genome clustering pipeline with integrated visualization.

**Status**: pre-alpha (M0 scaffold). Not yet runnable.

## Goal

A single-command, Docker-shipped replacement for BPGA (Bacterial Pan Genome
Analysis) that runs on macOS and Linux:

1. Takes a folder of GenBank files as input.
2. Runs one of four protein clustering engines: **MMseqs2** (default for speed),
   **DIAMOND deepclust**, **CD-HIT**, or **usearch12** (BPGA-compat).
3. Unifies cluster outputs into a DNMB-native Parquet/TSV dataframe schema.
4. Generates flower plot, pan/core curve, phylogeny, and a single HTML report.

The R visualization layer from
[BPGAconverter](https://github.com/JAEYOONSUNG/BPGAconverter) is absorbed
directly into DNMBcluster as an internal R package (`R/`).

## Quick start

Put your GenBank files (`.gb` / `.gbk` / `.gbff`) in a folder and point
DNMBcluster at it:

```bash
docker run --rm \
  -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) \
  --tmpfs /tmp:size=8g,mode=1777 \
  -v $PWD:/data \
  ghcr.io/jaeyoonsung/dnmbcluster:latest \
  run /data/my_genbanks
```

Outputs land in `./results/` by default: `dnmb/*.parquet` (7 tables
including the unified `clusters.parquet` with bidirectional identity
and coverage), `gene_presence_absence.csv` (Roary-compatible), and
`plots/*.pdf` (flower, pan-core, category bar).

### Options

- `--tool {mmseqs2,diamond,cd-hit,usearch12}` — clustering engine.
  Default is `mmseqs2`. See [BENCHMARK.md](BENCHMARK.md) for per-engine
  tradeoffs.
- `--fast` — skip the bidirectional MMseqs2 alignment pass.
  `clusters.parquet` alignment columns land as null but end-to-end
  wall time drops ~6.8× on 10 genomes.
- `--level {protein,nucleotide}` — default protein. Nucleotide mode
  clusters CDS DNA sequences instead of translations.
- `--identity FLOAT --coverage FLOAT` — sequence identity and
  coverage thresholds. Default 0.5 / 0.8 for protein, 0.7 / 0.8 for
  nucleotide.
- `--threads N` — default auto-detect.
- `--max-ram 8G` — propagated to clustering tools as a memory cap.
- `--parse-only` — stop after the GenBank parsing stage.

## Roadmap

- [x] **M0** repo scaffold
- [x] **M1** GenBank parser + parallel protein extraction
- [x] **M2** MMseqs2 engine + unified cluster TSV
- [x] **M3** Remaining engines (DIAMOND, CD-HIT, usearch12)
- [x] **M4** DNMB-native Parquet dataframes + Roary-compatible CSV
- [x] **M5** R visualization layer (absorbed from BPGAconverter)
- [x] **M6** Geobacillus-10 benchmark (see [BENCHMARK.md](BENCHMARK.md))
- [ ] **M7** Multi-arch CI (amd64/arm64) + GHCR release

## Benchmark snapshot (10 Geobacillus genomes, 33K proteins)

| Engine | Wall | RSS | Clusters |
|---|---:|---:|---:|
| `mmseqs2 --fast` | **9.9 s** | 329 MB | 6508 |
| `diamond` | 13.6 s | 324 MB | 5709 |
| `mmseqs2` default (bidirectional aln) | 67.4 s | 1004 MB | 6508 |
| `cd-hit` | 70.6 s | 323 MB | 5862 |

Full methodology and per-engine analysis in [BENCHMARK.md](BENCHMARK.md).

## Repository layout

```
DNMBcluster/
├── src/dnmbcluster/        # Python pipeline (orchestration + clustering)
│   ├── cli.py              # Click CLI entry point
│   ├── engines/            # pluggable clustering backends
│   └── ...
├── R/                      # DNMBcluster R package (visualization)
├── tests/
│   └── fixtures/
│       └── geobacillus10/  # benchmark GenBank set
├── Dockerfile              # multi-stage (usearch12 builder + conda runtime)
├── env.yml                 # conda env specification
├── docker-entrypoint.sh    # gosu + tini user-id handling
└── pyproject.toml
```

## License

MIT for DNMBcluster itself. Bundled third-party tools retain their original
licenses (usearch12: GPL-3.0, DIAMOND: GPL-3.0, CD-HIT: GPL-2.0,
MMseqs2: MIT).
