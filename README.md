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

## Quick start (planned — M6+)

```bash
docker run --rm \
  -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) \
  -v $PWD:/data \
  ghcr.io/jaeyoonsung/dnmbcluster:latest \
  run /data/my_genbanks
```

## Roadmap

- [x] **M0** repo scaffold
- [ ] **M1** GenBank parser + parallel protein extraction
- [ ] **M2** MMseqs2 engine + unified cluster TSV
- [ ] **M3** Remaining engines (DIAMOND, CD-HIT, usearch12)
- [ ] **M4** DNMB-native Parquet dataframes + Roary-compatible CSV
- [ ] **M5** R visualization layer (absorbed from BPGAconverter)
- [ ] **M6** Speed optimization + benchmark on 10 Geobacillus genomes
- [ ] **M7** Multi-arch CI (amd64/arm64) + GHCR release

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
