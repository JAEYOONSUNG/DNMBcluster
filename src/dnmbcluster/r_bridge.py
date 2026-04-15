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

    try:
        result = subprocess.run(
            [rscript, "--vanilla", "-e", r_expr],
            check=True,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
        stdout = exc.stdout.decode(errors="replace") if exc.stdout else ""
        raise RBridgeError(
            f"Rscript failed (exit {exc.returncode}):\n"
            f"STDOUT:\n{stdout}\nSTDERR:\n{stderr}"
        ) from exc

    if not plots_dir.exists():
        raise RBridgeError(
            f"R completed but plots directory does not exist: {plots_dir}"
        )
    return plots_dir
