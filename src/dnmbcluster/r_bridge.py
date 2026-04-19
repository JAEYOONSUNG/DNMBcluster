"""Thin bridge to the DNMBcluster R package.

The R package lives in ``R/`` of this repo. It is installed once at
Docker build time via ``R CMD INSTALL /app/R`` and thereafter is just
``library(DNMBcluster)`` from any Rscript process. This module invokes
Rscript in a subprocess with ``--vanilla`` to skip Rprofile/Renviron
for determinism (SPEED.md Section 9).
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


class RBridgeError(RuntimeError):
    pass


def find_rscript() -> str:
    """Return the path to Rscript, preferring ``$DNMBCLUSTER_RSCRIPT``."""
    override = os.environ.get("DNMBCLUSTER_RSCRIPT")
    if override:
        return override
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise RBridgeError(
            "Rscript not found on PATH. Install R (conda-forge r-base) or "
            "run inside the DNMBcluster Docker image."
        )
    return rscript


def run_r_plots(results_dir: Path) -> Path:
    """Invoke ``DNMBcluster::run_dnmb_plot`` on a results directory.

    Returns the path to the plots output directory on success.
    """
    rscript = find_rscript()
    results_dir = results_dir.resolve()
    plots_dir = results_dir / "plots"

    # Build a small R one-liner; escape Windows-style backslashes defensively
    # even though DNMBcluster targets macOS/Linux.
    path_literal = str(results_dir).replace("\\", "/")
    r_expr = (
        'suppressPackageStartupMessages(library(DNMBcluster)); '
        f'DNMBcluster::run_dnmb_plot("{path_literal}")'
    )

    env = os.environ.copy()
    env.setdefault("R_COMPILE_PKGS", "0")
    env.setdefault("OMP_NUM_THREADS", "1")

    # Inherit parent stdout/stderr so R plotting progress appears live
    # rather than buffering until completion. run_dnmb_plot is fast (a
    # few seconds) on 10-genome inputs so live output is not noisy; on
    # larger inputs the live stream is critical UX feedback.
    try:
        subprocess.run(
            [rscript, "--vanilla", "-e", r_expr],
            check=True,
            env=env,
        )
    except subprocess.CalledProcessError as exc:
        raise RBridgeError(
            f"Rscript failed (exit {exc.returncode}) — see R output above."
        ) from exc

    if not plots_dir.exists():
        raise RBridgeError(
            f"R completed but plots directory does not exist: {plots_dir}"
        )
    return plots_dir


def run_r_sota_phylogenomics(
    results_dir: Path,
    out_subdir: str = "orthofinder_like",
    threads: int = 4,
    workers: int = 1,
    bootstrap: int = 100,
    og_bootstrap_reps: int = 100,
    min_dup_support: float = 50,
    plot: bool = True,
    report: bool = True,
    generax: bool = True,
    foldseek_dir: str | None = None,
    foldseek_mode: str = "pdb",
    foldseek_prosT5: str | None = None,
    codeml: bool = False,
    hyphy: bool = False,
    cds_source: Path | None = None,
) -> Path:
    """Invoke ``DNMBcluster::run_phylogenomics_sota`` on a results dir.

    Turns on every accuracy-oriented option built into the R package:
    graph-refined OGs, ML gene trees with bootstrap, auto partitioned
    ML supermatrix with site-bootstrap, support-gated STRIDE / HOG /
    DTL, and a species-tree DTL-event overlay PDF.

    Selection analysis: pass ``hyphy=True`` with a ``cds_source`` path
    (single nucleotide FASTA keyed by the same headers as the AA
    alignments) to run HyPhy FEL+BUSTED+aBSREL per single-copy OG.
    """
    rscript = find_rscript()
    results_dir = results_dir.resolve()
    out_dir = results_dir / out_subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    results_lit = str(results_dir).replace("\\", "/")
    out_lit = str(out_dir).replace("\\", "/")
    generax_lit = "TRUE" if generax else "FALSE"
    fsk_dir_lit = (
        f'"{str(Path(foldseek_dir).resolve()).replace(chr(92), "/")}"'
        if foldseek_dir else "NULL"
    )
    fsk_p5_lit = (
        f'"{str(Path(foldseek_prosT5).resolve()).replace(chr(92), "/")}"'
        if foldseek_prosT5 else "NULL"
    )
    # HyPhy + codeml delegation: when cds_source is provided,
    # run_phylogenomics_sota() runs back_translate → codon MSA → HyPhy
    # FEL+BUSTED+aBSREL (+ optional codeml M0) in one pass. The
    # selection_workers default equals the main `workers` arg inside
    # the R function; we bump threads allocation via the R-side split.
    if cds_source is not None:
        cds_lit = f'"{str(Path(cds_source).resolve()).replace(chr(92), "/")}"'
    else:
        cds_lit = "NULL"
    hyphy_lit = "\"auto\"" if hyphy else "FALSE"
    codeml_lit = "\"auto\"" if codeml else "FALSE"

    r_expr = (
        "suppressPackageStartupMessages(library(DNMBcluster)); "
        f'dnmb <- DNMBcluster::load_dnmb("{results_lit}"); '
        f'DNMBcluster::run_phylogenomics_sota('
        f'dnmb, "{out_lit}", '
        f"threads = {int(threads)}L, workers = {int(workers)}L, "
        f"bootstrap = {int(bootstrap)}L, "
        f"og_bootstrap_reps = {int(og_bootstrap_reps)}L, "
        f"min_dup_support = {float(min_dup_support)}, "
        f"plot = {'TRUE' if plot else 'FALSE'}, "
        f"report = {'TRUE' if report else 'FALSE'}, "
        f"generax = {generax_lit}, "
        f"foldseek_dir = {fsk_dir_lit}, "
        f'foldseek_mode = "{foldseek_mode}", '
        f"foldseek_prosT5 = {fsk_p5_lit}, "
        f"hyphy = {hyphy_lit}, "
        f"codeml = {codeml_lit}, "
        f"cds_source = {cds_lit}"
        ")"
    )

    env = os.environ.copy()
    env.setdefault("R_COMPILE_PKGS", "0")

    # Inherit parent stdout/stderr. SOTA phylogenomics runs 10-20 min
    # with many per-stage progress prints (OG refinement, gene trees,
    # bootstrap, DTL); buffering until exit meant the user saw a frozen
    # terminal for the whole stage. Live streaming matters here much
    # more than post-mortem error capture.
    try:
        subprocess.run(
            [rscript, "--vanilla", "-e", r_expr],
            check=True, env=env,
        )
    except subprocess.CalledProcessError as exc:
        raise RBridgeError(
            f"SOTA phylogenomics Rscript failed (exit {exc.returncode}) "
            "— see R output above."
        ) from exc

    tree_path = out_dir / "species_tree_rooted.nwk"
    if not tree_path.exists():
        raise RBridgeError(
            f"SOTA phylogenomics completed but tree missing: {tree_path}"
        )
    return out_dir
