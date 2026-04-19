#' State-of-the-art phylogenomics pipeline preset
#'
#' One-call wrapper around `run_orthofinder_like()` that turns on the
#' full set of accuracy-oriented options built into DNMBcluster:
#'
#' \itemize{
#'   \item Length-normalized graph refinement of upstream OGs
#'         (`refine_graph = TRUE`, Louvain by default; MCL if the
#'         `MCL` package is installed).
#'   \item ML gene trees with bootstrap support
#'         (`method = "ml"`, `og_bootstrap_reps = 100`).
#'   \item Auto species-tree strategy: partitioned ML supermatrix when
#'         enough single-copy OGs exist, STAG otherwise
#'         (`species_tree_method = "auto"`, `supermatrix_partitioned = TRUE`).
#'   \item Site-bootstrap species-tree support (`bootstrap = 100`).
#'   \item Support-gated duplication calling for STRIDE rooting and
#'         HOG cutting (`min_dup_support = 50`).
#'   \item Full DTL (S/D/T/loss) reconciliation and per-edge loss
#'         aggregation (`dtl = TRUE`).
#'   \item Writes a species-tree overlay PDF with S/D/T pie charts and
#'         loss-weighted edge widths (`dtl_events_overlay.pdf`).
#' }
#'
#' Everything is still a pure-R pipeline — no external GeneRax/ALE is
#' invoked. For rigorous cost-based DTL, take the reconciled per-OG
#' trees and run them through `export_og_bundles()` + GeneRax/ALE
#' separately.
#'
#' @inheritParams run_orthofinder_like
#' @param bootstrap Species-tree site-bootstrap replicates.
#'   Default 100.
#' @param og_bootstrap_reps Per-OG gene-tree bootstrap replicates.
#'   Default 100.
#' @param min_dup_support Bootstrap threshold (%) below which a
#'   duplication is ignored by STRIDE, HOG cutting, and DTL D/T
#'   calling. Default 50.
#' @param long_branch_quantile DTL long-branch gate for transfer
#'   candidates. Default 0.9 (top-decile sibling edge required).
#' @param plot If TRUE, writes `dtl_events_overlay.pdf`. Default TRUE.
#' @param report If TRUE and the `rmarkdown` package is available,
#'   renders the OrthoFinder-parity HTML report at the end.
#' @param generax If `"auto"` (default), runs `run_generax()` when the
#'   `generax` binary is on PATH and overwrites the heuristic
#'   `res$dtl` events with rigorous UndatedDTL results. `TRUE` forces
#'   and errors if binary is missing; `FALSE` never runs it.
#' @param foldseek_dir Optional path to a directory containing
#'   per-protein PDB files (or an AA FASTA when `foldseek_mode =
#'   "prosT5"`). When set and the `foldseek` binary is on PATH,
#'   per-HOG structural support is computed and joined onto the
#'   HOG TSVs. Default NULL (disabled).
#' @param foldseek_mode `"pdb"` (default) or `"prosT5"`; see
#'   `run_foldseek_hog()`.
#' @param foldseek_prosT5 Path to ProstT5 weights when
#'   `foldseek_mode = "prosT5"`.
#' @param hyphy If `"auto"` (default), runs modern HyPhy selection tests
#'   (FEL + BUSTED[S] + aBSREL) when the `hyphy` binary is on PATH AND
#'   `cds_source` is supplied. `TRUE` errors if the binary is missing;
#'   `FALSE` never runs it. 2025 SOTA for per-OG / per-site / per-branch
#'   dN/dS.
#' @param cds_source Path to a nucleotide CDS FASTA (or per-genome dir
#'   of FASTAs) used by `back_translate_og_codons()` to build codon
#'   alignments. Required when `hyphy` or `codeml` is on. NULL disables.
#' @param codeml If `"auto"` (default), also runs `run_codeml_m0()` as
#'   a fast baseline per-OG ω when the `codeml` binary is on PATH. Kept
#'   as a sensitivity check alongside HyPhy per 2024–25 reviewer norms
#'   (Jeffares et al. MBE 2023).
#' @param selection_methods Subset of `c("fel","busted","absrel")` to
#'   run under HyPhy. Default all three.
#' @param selection_workers Parallel HyPhy workers. Default
#'   `max(1, workers)`.
#' @return Named list from `run_orthofinder_like()`, extended with
#'   `overlay_pdf` (if `plot = TRUE`) and `report_html` (if `report = TRUE`).
#' @export
run_phylogenomics_sota <- function(dnmb,
                                    out_dir,
                                    min_seqs = 3L,
                                    max_seqs = 500L,
                                    threads = 4L,
                                    workers = 1L,
                                    bootstrap = 100L,
                                    og_bootstrap_reps = 100L,
                                    min_dup_support = 50,
                                    long_branch_quantile = 0.9,
                                    outgroup = NULL,
                                    xeno_z = 3.0,
                                    plot = TRUE,
                                    report = TRUE,
                                    generax = "auto",
                                    foldseek_dir = NULL,
                                    foldseek_mode = c("pdb", "prosT5"),
                                    foldseek_prosT5 = NULL,
                                    hyphy = "auto",
                                    cds_source = NULL,
                                    codeml = "auto",
                                    selection_methods = c("fel", "busted", "absrel"),
                                    selection_workers = NULL,
                                    verbose = TRUE) {
  foldseek_mode <- match.arg(foldseek_mode)
  selection_methods <- match.arg(selection_methods, several.ok = TRUE)
  if (is.null(selection_workers)) selection_workers <- max(1L, workers)
  refine_method <- if (requireNamespace("MCL", quietly = TRUE)) "mcl" else "louvain"
  tree_method <- if (nzchar(Sys.which("iqtree2")) || nzchar(Sys.which("iqtree"))) "iqtree" else "ml"
  if (verbose) message("[sota] refine_method=", refine_method,
                        ", tree_method=", tree_method,
                        ", bootstrap=", bootstrap,
                        ", og_bootstrap_reps=", og_bootstrap_reps,
                        ", min_dup_support=", min_dup_support)

  res <- run_orthofinder_like(
    dnmb,
    out_dir,
    method                  = tree_method,
    min_seqs                = min_seqs,
    max_seqs                = max_seqs,
    use_multi_copy          = FALSE,
    outgroup                = outgroup,
    xeno_z                  = xeno_z,
    threads                 = threads,
    verbose                 = verbose,
    refine_graph            = TRUE,
    refine_method           = refine_method,
    bootstrap               = bootstrap,
    force                   = FALSE,
    workers                 = workers,
    species_tree_method     = "auto",
    trim_gaps               = TRUE,
    max_gap_frac            = 0.5,
    aa_model                = "LG",
    gamma_k                 = 4L,
    og_bootstrap_reps       = og_bootstrap_reps,
    min_dup_support         = min_dup_support,
    supermatrix_partitioned = TRUE,
    dtl                     = TRUE
  )

  # Preflight: `generax = TRUE` is documented as a hard requirement — if
  # the binary is missing, stop() rather than silently degrading. Auto
  # mode only runs when the binary is present, and FALSE never runs it.
  if (identical(generax, TRUE) && !nzchar(Sys.which("generax"))) {
    stop("run_phylogenomics_sota(generax = TRUE) but 'generax' is not on PATH. ",
         "Install it (bioconda: generax) or set generax = \"auto\"/FALSE.")
  }
  want_generax <- identical(generax, TRUE) ||
    (identical(generax, "auto") && nzchar(Sys.which("generax")))
  if (want_generax) {
    if (verbose) message("[sota] running GeneRax UndatedDTL reconciliation...")
    bundle <- tryCatch(
      export_og_bundles(dnmb, res$og_result, out_dir,
                        tools = "generax",
                        species_tree = file.path(out_dir, "species_tree_rooted.nwk"),
                        aa_model = "LG"),
      error = function(e) { message("[sota] export_og_bundles failed: ",
                                    conditionMessage(e)); NULL })
    if (!is.null(bundle)) {
      grx_out <- file.path(out_dir, "generax")
      grx <- tryCatch(run_generax(bundle$generax, grx_out, force = TRUE,
                                   threads = threads),
                      error = function(e) { message("[sota] generax failed: ",
                                                     conditionMessage(e)); NULL })
      if (!is.null(grx)) {
        # Keep GeneRax results in their own slot. The heuristic
        # `res$dtl` uses a per-node schema (cluster_id, g_node,
        # g_n_leaves, sp_lca, event) that the overlay PDF consumes;
        # GeneRax output is per-family totals (cluster_id, event,
        # count) and per-family summary, so overwriting res$dtl makes
        # the object, the on-disk dtl_*.tsv, and the overlay disagree.
        res$generax <- grx
        if (verbose) message("[sota] GeneRax events stored in res$generax; ",
                              "overlay still uses heuristic res$dtl schema.")
      }
    }
  }

  if (!is.null(foldseek_dir) && nzchar(Sys.which("foldseek")) &&
      length(res$hogs) && nrow(res$hogs)) {
    if (verbose) message("[sota] running Foldseek structural orthology...")
    hog_long <- .sota_collect_hog_table(out_dir)
    fsk_out <- file.path(out_dir, "foldseek")
    fsk <- tryCatch(run_foldseek_hog(dnmb, hog_long, foldseek_dir, fsk_out,
                                      mode = foldseek_mode,
                                      prosT5_weights = foldseek_prosT5,
                                      threads = threads),
                    error = function(e) { message("[sota] foldseek failed: ",
                                                   conditionMessage(e)); NULL })
    if (!is.null(fsk)) {
      res$foldseek <- fsk
      ann_ok <- tryCatch({
        annotate_hogs_with_foldseek(out_dir, fsk$per_hog); TRUE
      }, error = function(e) {
        message("[sota] annotate_hogs_with_foldseek failed: ",
                conditionMessage(e)); FALSE
      })
      if (isTRUE(ann_ok) && verbose)
        message("[sota] annotated HOG TSVs with Foldseek support.")
    }
  }

  # --- Selection analysis: HyPhy (SOTA) + optional codeml M0 baseline. ---
  if (identical(hyphy, TRUE) && !nzchar(Sys.which("hyphy"))) {
    stop("run_phylogenomics_sota(hyphy = TRUE) but 'hyphy' is not on PATH. ",
         "Install it (bioconda: hyphy) or set hyphy = \"auto\"/FALSE.")
  }
  want_hyphy <- identical(hyphy, TRUE) ||
                (identical(hyphy, "auto") && nzchar(Sys.which("hyphy")))
  want_codeml <- identical(codeml, TRUE) ||
                 (identical(codeml, "auto") && nzchar(Sys.which("codeml")))

  if ((want_hyphy || want_codeml) && !is.null(cds_source) &&
      !is.null(res$og_result) && nrow(res$og_result)) {
    if (verbose) message("[sota] back-translating SC-OG AA alignments to codon MSA...")
    has_tree <- !is.na(res$og_result$tree_path)
    sc_flag  <- if ("single_copy" %in% names(res$og_result))
                  isTRUE(any(res$og_result$single_copy, na.rm = TRUE)) else FALSE
    sc_ok <- if (sc_flag) {
      res$og_result[has_tree & !is.na(res$og_result$single_copy) &
                      res$og_result$single_copy, , drop = FALSE]
    } else {
      # No single-copy column or none flagged — fall back to every OG
      # that has an alignment + tree.
      res$og_result[has_tree, , drop = FALSE]
    }
    bt <- tryCatch(
      back_translate_og_codons(sc_ok, cds_source = cds_source,
                               write_phylip = TRUE, verbose = FALSE),
      error = function(e) {
        message("[sota] back_translate_og_codons failed: ", conditionMessage(e))
        NULL
      })
    if (!is.null(bt)) {
      # `which()` drops NAs; bt$ok is logical so this filters to successful rows.
      good <- bt[which(bt$ok), , drop = FALSE]
      og_dirs <- unique(dirname(good$nuc_path[!is.na(good$nuc_path)]))
      res$selection <- list(back_translate = bt, og_dirs = og_dirs)

      if (length(og_dirs)) {
        if (want_hyphy) {
          if (verbose) message("[sota] running HyPhy batch (",
                                paste(selection_methods, collapse = "+"),
                                ") on ", length(og_dirs), " OGs...")
          hy <- tryCatch(
            run_hyphy_batch(og_dirs,
                            work_dir = file.path(out_dir, "hyphy"),
                            methods  = selection_methods,
                            workers  = selection_workers,
                            threads  = max(1L, threads %/% max(1L, selection_workers))),
            error = function(e) {
              message("[sota] run_hyphy_batch failed: ", conditionMessage(e))
              NULL
            })
          res$selection$hyphy <- hy
        }
        if (want_codeml) {
          if (verbose) message("[sota] running codeml M0 baseline on ",
                                length(og_dirs), " OGs...")
          cml <- tryCatch(
            run_codeml_m0(og_dirs,
                          aln_name  = "aln.nuc.phy",
                          tree_name = "tree.nwk",
                          work_dir  = file.path(out_dir, "codeml")),
            error = function(e) {
              message("[sota] run_codeml_m0 failed: ", conditionMessage(e))
              NULL
            })
          res$selection$codeml <- cml
        }
      } else if (verbose) {
        message("[sota] no OGs passed back-translation; skipping selection.")
      }
    }
  } else if (want_hyphy && is.null(cds_source) && verbose) {
    message("[sota] hyphy requested but cds_source = NULL; ",
            "selection analysis skipped (DNMBcluster does not keep CDS).")
  }

  if (isTRUE(plot) && !is.null(res$dtl) && nrow(res$dtl$events)) {
    overlay <- file.path(out_dir, "dtl_events_overlay.pdf")
    tryCatch({
      plot_dtl_branch_events(res$dtl$events, res$species_tree,
                              losses_tbl = res$dtl$losses,
                              out_pdf = overlay)
      res$overlay_pdf <- overlay
      if (verbose) message("[sota] wrote DTL overlay -> ", overlay)
    }, error = function(e) {
      message("[sota] overlay plot failed: ", conditionMessage(e))
    })
  }

  if (isTRUE(report) && requireNamespace("rmarkdown", quietly = TRUE)) {
    res$report_html <- tryCatch(
      render_orthofinder_report(out_dir,
                                 run_title = "DNMBcluster SOTA pipeline"),
      error = function(e) {
        message("[sota] report render failed: ", conditionMessage(e))
        NULL
      })
  }

  if (verbose) {
    message("[sota] done. outputs in ", out_dir)
  }
  invisible(res)
}


# Flatten per-node HOG TSVs into one (hog_id, protein_uid) table for
# downstream structural-support scoring.
.sota_collect_hog_table <- function(out_dir) {
  files <- list.files(file.path(out_dir, "HOGs"),
                      pattern = "^N\\d+_.*\\.tsv$", full.names = TRUE)
  rows <- list()
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    body  <- lines[!startsWith(lines, "#") & !startsWith(lines, "hog_id")]
    if (!length(body)) next
    parts <- strsplit(body, "\t", fixed = TRUE)
    for (p in parts) {
      if (length(p) < 4L) next
      leaves <- strsplit(p[4], ",", fixed = TRUE)[[1]]
      uids <- suppressWarnings(as.integer(sub("_g.*$", "",
                                               sub("^p", "", leaves))))
      rows[[length(rows) + 1L]] <- tibble::tibble(
        hog_id = p[1], protein_uid = uids
      )
    }
  }
  dplyr::bind_rows(rows)
}
