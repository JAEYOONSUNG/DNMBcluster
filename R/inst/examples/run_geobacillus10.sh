#!/usr/bin/env bash
# =====================================================================
# Full DNMB pipeline on the real Geobacillus-10 fixture — --max preset
# ---------------------------------------------------------------------
#   inputs  : tests/fixtures/geobacillus10/*.gbff  (10 genomes, 87 MB)
#   engine  : mmseqs2  id=0.6 cov=0.8 protein-level, bidirectional realign
#
#   `dnmbcluster run --max` flips on every optional stage in one flag:
#     --phylo              single-copy core → FAMSA → ClipKIT →
#                          IQ-TREE 3 ModelFinder + UFBoot + -bnni
#                          (Dollo gain/loss auto-runs on the tree)
#     --ani                skani + POCP genome-genome matrices
#     --annotate           fast EggNOG (DIAMOND --fast + SQLite)
#     --orthofinder-like   SOTA R stage: graph-refined OGs + ML gene
#                          trees + DTL reconciliation + optional
#                          GeneRax + Foldseek
#     --codeml             PAML codeml M0 dN/dS per SC-OG
#     --hyphy              HyPhy FEL + BUSTED + aBSREL per SC-OG
#                          (auto-extracts cds.fna from GenBank)
#
#   `dnmbcluster run` auto-invokes the R plotting side at the end, so
#   this is still single-command.
# ---------------------------------------------------------------------
# Requires: conda env `dnmb-test`
#   mmseqs2, famsa/mafft, clipkit/trimal, iqtree3, skani, diamond,
#   codeml, hyphy, generax, foldseek, and the EggNOG cache at
#   ~/.dnmb-cache/db_modules/eggnog/.
#
# Typical runtime on an M-series laptop: ~25–40 min with --max
# (7.5 min clustering + ~8 min IQ-TREE + ~1 min ANI + ~2 min EggNOG
#  + ~10–20 min OrthoFinder-parity + HyPhy + codeml).
# =====================================================================

set -eo pipefail   # -u trips conda's activate.d on unset GFORTRAN

REPO="/Users/JaeYoon/Dropbox/0.Personal folder/5. Bioinformatics/DNMBcluster"
FIXTURE="${REPO}/tests/fixtures/geobacillus10"
OUT_ROOT="${DNMB_GEO_OUT:-$(mktemp -d -t dnmb_geobacillus_XXXX)}"

source /Users/JaeYoon/miniforge3/etc/profile.d/conda.sh
conda activate dnmb-test
# mafft lives at /usr/local/bin/mafft on this host; conda env doesn't
# ship it. Keep /usr/local/bin after the conda env so both resolve.
export PATH="$PATH:/usr/local/bin"

echo "[geo] inputs      : $FIXTURE  ($(ls "$FIXTURE" | wc -l | tr -d ' ') .gbff files)"
echo "[geo] output root : $OUT_ROOT"
echo "[geo] tools       : mmseqs=$(which mmseqs) famsa=$(which famsa 2>/dev/null || echo -) clipkit=$(which clipkit 2>/dev/null || echo -) iqtree3=$(which iqtree3) skani=$(which skani) hyphy=$(which hyphy) codeml=$(which codeml)"
echo

/usr/bin/time -p dnmbcluster run "$FIXTURE" \
    --output "$OUT_ROOT" \
    --tool mmseqs2 \
    --level protein \
    --identity 0.6 \
    --coverage 0.8 \
    --threads 0 \
    --alignment \
    --max

echo
echo "--- Parquet artefacts ---"
find "$OUT_ROOT/dnmb" -name '*.parquet' | sort
echo
echo "--- OrthoFinder-parity outputs ---"
find "$OUT_ROOT/orthofinder_like" -maxdepth 2 -type f 2>/dev/null | sort
echo
echo "--- Plots (PDF) ---"
ls "$OUT_ROOT/plots" 2>/dev/null | sort
echo
echo "[geo] report : $OUT_ROOT/report.html"
echo "[geo] xlsx   : $OUT_ROOT/comparative_genomics_10.xlsx"
