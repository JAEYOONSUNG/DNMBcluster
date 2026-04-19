#' Run the full OrthoFinder-parity workflow on a DNMB results directory
#'
#' One-shot wrapper that chains the four phase functions:
#' \enumerate{
#'   \item `per_og_trees()`             — per-OG alignment + gene trees
#'   \item `stag_species_tree()`        — species tree from gene trees
#'   \item `stride_root()`              — duplication-consistent rooting
#'   \item `cut_hogs()`                 — HOG tables per internal node
#'   \item `reconcile_relationships()`  — ortholog / paralog / xenolog labels
#' }
#'
#' Phase E-a (tool exports) is a separate call because not every user
#' wants GeneRax/PAML/MCScanX bundles on every run.
#'
#' All output lands under `out_dir`:
#' \preformatted{
#'   orthogroups/OG_<id>/{aln.fasta,tree.nwk}
#'   orthogroup_trees.tsv
#'   single_copy_OGs.tsv
#'   species_tree_rooted.nwk
#'   HOGs/N<node>_<clade>.tsv
#'   relationships.parquet
#' }
#'
#' @param dnmb Result of `load_dnmb()`.
#' @param out_dir Destination directory.
#' @param method Tree method: `"nj"` (fast distance), `"ml"` (phangorn
#'   NNI, accurate), or `"iqtree"` (IQ-TREE `--fast` + UFBoot, fastest
#'   ML for large OGs when the `iqtree2`/`iqtree` binary is on PATH;
#'   auto-falls back to phangorn ML otherwise).
#' @param min_seqs Skip OGs smaller than this.
#' @param max_seqs Subsample OGs larger than this.
#' @param use_multi_copy Include multi-copy OGs in STAG aggregation.
#' @param outgroup Force-rooting outgroup (overrides STRIDE).
#' @param xeno_z Xenolog-candidate Z-score cutoff.
#' @param threads Thread count for engines.
#' @param verbose Echo per-phase progress.
#' @param refine_graph If TRUE, pre-refine DNMB's identity-based OGs with
#'   `refine_ogs_graph()` (length-normalized edges + community split)
#'   before building gene trees. Closes OrthoFinder's RBH-reweighting gap.
#' @param refine_method Passed to `refine_ogs_graph()` when refine_graph=TRUE.
#' @param bootstrap Species-tree bootstrap replicates for STAG node support.
#' @param force If FALSE (default), reuse existing per-OG alignments and
#'   trees on disk; if TRUE, rebuild everything. Passed to
#'   `per_og_trees()` for resume-friendly reruns.
#' @param workers Parallel workers for per-OG tree building. Requires
#'   the `future.apply` Suggests package when > 1.
#' @param species_tree_method Which species-tree builder to use:
#'   `"stag"` (default — distance aggregation, handles multi-copy),
#'   `"supermatrix"` (concatenated SC-OG alignment + ML/NJ), or
#'   `"auto"` (supermatrix if >= 20 SC-OGs, else STAG).
#' @param trim_gaps Drop high-gap alignment columns before tree building
#'   (passed to `per_og_trees`). Default TRUE.
#' @param max_gap_frac Column-gap threshold for `trim_gaps`. Default 0.5.
#' @param aa_model Amino-acid model for ML branch. Default `"LG"`.
#' @param gamma_k Γ rate categories for ML. Default 4 (0 disables).
#' @param og_bootstrap_reps Per-OG bootstrap replicates (0 disables).
#'   Support percentages land on `tree$node.label` in each OG tree.
#' @param min_dup_support When `og_bootstrap_reps > 0`, gene-tree
#'   duplication nodes with bootstrap support below this percentage are
#'   ignored by STRIDE rooting and HOG cutting. Default 0.
#' @param supermatrix_partitioned If TRUE and ML supermatrix is used,
#'   fit each SC-OG as its own partition via `phangorn::pmlPart`. Falls
#'   back to concatenated ML on error. Default FALSE.
#' @param dtl If TRUE, run per-OG DTL reconciliation via `reconcile_dtl()`
#'   and write `dtl_per_og.tsv` + `dtl_events.tsv` to `out_dir`. Default
#'   FALSE — the reconciliation is a heuristic LCA pass, not a cost-space
#'   search; for rigorous DTL use `export_og_bundles()` + GeneRax/ALE.
#' @return Named list with keys `og_result`, `species_tree`, `hogs`,
#'   `relationships`, `dtl` (if `dtl = TRUE`), and `refined` (if `refine_graph`).
#' @export
run_orthofinder_like <- function(dnmb,
                                 out_dir,
                                 method = c("nj", "ml", "iqtree"),
                                 min_seqs = 3L,
                                 max_seqs = 500L,
                                 use_multi_copy = FALSE,
                                 outgroup = NULL,
                                 xeno_z = 3.0,
                                 threads = 2L,
                                 verbose = TRUE,
                                 refine_graph = FALSE,
                                 refine_method = c("louvain", "mcl", "walktrap"),
                                 bootstrap = 0L,
                                 force = FALSE,
                                 workers = 1L,
                                 species_tree_method = c("stag", "supermatrix", "auto"),
                                 trim_gaps = TRUE,
                                 max_gap_frac = 0.5,
                                 aa_model = "LG",
                                 gamma_k = 4L,
                                 og_bootstrap_reps = 0L,
                                 min_dup_support = 0,
                                 supermatrix_partitioned = FALSE,
                                 dtl = FALSE) {
  species_tree_method <- match.arg(species_tree_method)
  refine_method <- match.arg(refine_method)
  method <- match.arg(method)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  t_start <- Sys.time()
  refined <- NULL
  if (refine_graph) {
    if (verbose) message("[0/5] Refining OGs via length-normalized graph...")
    refined <- refine_ogs_graph(dnmb, out_dir = out_dir,
                                 method = refine_method, threads = threads)
    # Rewrite dnmb$clusters so downstream phases see the split OGs.
    # refined_id is character → hash to new integer cluster ids.
    lvl <- unique(refined$refined_id)
    new_cid <- stats::setNames(seq_along(lvl), lvl)
    remap <- tibble::tibble(
      protein_uid = refined$protein_uid,
      new_cluster = as.integer(new_cid[refined$refined_id])
    )
    dnmb$clusters <- dnmb$clusters %>%
      dplyr::left_join(remap, by = "protein_uid") %>%
      dplyr::mutate(cluster_id = dplyr::coalesce(new_cluster, cluster_id)) %>%
      dplyr::select(-new_cluster)
  }

  if (verbose) message(sprintf("[%d/%d] Building per-OG trees...",
                                if (refine_graph) 1L else 1L,
                                if (refine_graph) 5L else 4L))
  og_result <- per_og_trees(dnmb, out_dir,
                            min_seqs = min_seqs, max_seqs = max_seqs,
                            method = method, threads = threads,
                            verbose = verbose, force = force,
                            workers = workers,
                            trim_gaps = trim_gaps,
                            max_gap_frac = max_gap_frac,
                            model = aa_model,
                            gamma_k = gamma_k,
                            bootstrap_reps = og_bootstrap_reps)

  n_sc <- sum(og_result$single_copy, na.rm = TRUE)
  st_method <- species_tree_method
  if (st_method == "auto") {
    st_method <- if (n_sc >= 20L) "supermatrix" else "stag"
    if (verbose) message(sprintf("  [auto] %d SC-OGs -> species_tree_method=%s",
                                  n_sc, st_method))
  }

  if (st_method == "supermatrix" && n_sc == 0L) {
    if (verbose) message("  supermatrix requested but no SC-OGs -- falling back to STAG (multi-copy).")
    st_method <- "stag"
  }

  if (st_method == "supermatrix") {
    if (verbose) message("[2] Supermatrix species tree (concatenated SC-OGs)...")
    # supermatrix_species_tree uses phangorn internally (no iqtree path);
    # map iqtree -> ml so SOTA (tree_method="iqtree") doesn't error here.
    sm_method <- if (identical(method, "iqtree")) "ml" else method
    sp_unrooted <- supermatrix_species_tree(
      og_result, dnmb, method = sm_method,
      model = aa_model, gamma_k = gamma_k,
      partitioned = isTRUE(supermatrix_partitioned),
      bootstrap = bootstrap
    )
  } else {
    if (verbose) message("[2] STAG + STRIDE species tree...")
    effective_multi <- use_multi_copy || n_sc == 0
    if (effective_multi && !use_multi_copy && verbose) {
      message("  no single-copy OGs present -- falling back to use_multi_copy=TRUE")
    }
    sp_unrooted <- stag_species_tree(og_result, dnmb,
                                      use_multi_copy = effective_multi,
                                      bootstrap = bootstrap)
  }
  species_tree <- stride_root(sp_unrooted, og_result, dnmb, outgroup = outgroup,
                              min_dup_support = min_dup_support)
  ape::write.tree(species_tree, file.path(out_dir, "species_tree_rooted.nwk"))

  if (verbose) message("[3/4] Cutting HOGs...")
  hogs <- cut_hogs(species_tree, og_result, dnmb, out_dir,
                   min_dup_support = min_dup_support)

  if (verbose) message("[4/4] Reconciling ortholog/paralog/xenolog labels...")
  relationships <- reconcile_relationships(og_result, species_tree, dnmb,
                                            out_dir = out_dir,
                                            xeno_z = xeno_z)

  dtl_result <- NULL
  if (isTRUE(dtl)) {
    if (verbose) message("[5] DTL reconciliation (per-OG LCA + species-overlap)...")
    dtl_result <- .run_dtl_batch(og_result, species_tree, dnmb, out_dir, verbose)
  }

  # Render the curated figure set. Each is a composite that subsumes
  # earlier sidecar PDFs; individual `plot_*` functions remain exported
  # for users who want them directly. Failures are non-fatal — data
  # artifacts are already on disk.
  tryCatch(
    plot_orthofinder_overview(out_dir, verbose = verbose),
    error = function(e) if (verbose)
      message("  [orthofinder_overview] skipped: ", conditionMessage(e))
  )
  tryCatch(
    plot_pangenome_landscape(out_dir, verbose = verbose),
    error = function(e) if (verbose)
      message("  [pangenome_landscape] skipped: ", conditionMessage(e))
  )
  tryCatch(
    plot_duplication_burden(out_dir, verbose = verbose),
    error = function(e) if (verbose)
      message("  [duplication_burden] skipped: ", conditionMessage(e))
  )
  if (!is.null(dtl_result) && nrow(dtl_result$events)) {
    tryCatch(
      plot_dtl_summary(out_dir, verbose = verbose),
      error = function(e) if (verbose)
        message("  [dtl_summary] skipped: ", conditionMessage(e))
    )
    tryCatch(
      plot_dtl_top_ogs(out_dir, verbose = verbose),
      error = function(e) if (verbose)
        message("  [dtl_top_ogs] skipped: ", conditionMessage(e))
    )
  }
  tryCatch({
    grid_pdf <- file.path(out_dir, "per_og_tree_grid.pdf")
    grid_plot <- per_og_tree_grid(og_result)
    ggplot2::ggsave(grid_pdf, grid_plot, width = 13, height = 9,
                    dpi = 200, limitsize = FALSE, bg = "#F4F1EA")
    if (verbose) message("  [per_og_tree_grid] wrote ", grid_pdf)
  }, error = function(e) if (verbose)
    message("  [per_og_tree_grid] skipped: ", conditionMessage(e)))

  if (verbose) {
    elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)
    extra <- if (!is.null(dtl_result)) sprintf(", DTL S/D/T=%d/%d/%d",
                                                sum(dtl_result$per_og$n_S),
                                                sum(dtl_result$per_og$n_D),
                                                sum(dtl_result$per_og$n_T)) else ""
    message(sprintf("[run_orthofinder_like] done in %ss. OGs=%d, SC-OGs=%d, HOG-nodes=%d, relationships=%d%s",
                    elapsed,
                    sum(!is.na(og_result$tree_path)),
                    sum(og_result$single_copy),
                    nrow(hogs),
                    nrow(relationships),
                    extra))
  }

  invisible(list(
    og_result     = og_result,
    species_tree  = species_tree,
    hogs          = hogs,
    relationships = relationships,
    dtl           = dtl_result,
    refined       = refined
  ))
}

.run_dtl_batch <- function(og_result, species_tree, dnmb, out_dir, verbose) {
  ok <- og_result[!is.na(og_result$tree_path), , drop = FALSE]
  if (!nrow(ok)) return(NULL)
  uid_to_key <- stats::setNames(as.character(dnmb$genome_meta$genome_key),
                                 as.character(dnmb$genome_meta$genome_uid))
  per_og <- vector("list", nrow(ok))
  events <- vector("list", nrow(ok))
  losses <- vector("list", nrow(ok))
  for (i in seq_len(nrow(ok))) {
    cid <- ok$cluster_id[i]
    tr <- tryCatch(ape::read.tree(ok$tree_path[i]), error = function(e) NULL)
    if (is.null(tr) || length(tr$tip.label) < 2L) next
    if (!ape::is.rooted(tr)) tr <- mad_root(tr)
    guid_str <- sub(".*_g", "", tr$tip.label)
    sp_keys <- unname(uid_to_key[guid_str])
    sp_keys[is.na(sp_keys)] <- guid_str[is.na(sp_keys)]
    tip_to_species <- stats::setNames(sp_keys, tr$tip.label)
    res <- tryCatch(reconcile_dtl(tr, species_tree,
                                   tip_to_species = tip_to_species),
                    error = function(e) NULL)
    if (is.null(res)) next
    per_og[[i]] <- tibble::tibble(
      cluster_id = cid,
      n_S = res$summary$n_S,
      n_D = res$summary$n_D,
      n_T = res$summary$n_T,
      n_loss = res$summary$n_loss
    )
    events[[i]] <- tibble::tibble(cluster_id = cid, res$events)
    if (!is.null(res$losses) && nrow(res$losses))
      losses[[i]] <- tibble::tibble(cluster_id = cid, res$losses)
  }
  per_og_tbl <- dplyr::bind_rows(per_og)
  events_tbl <- dplyr::bind_rows(events)
  losses_tbl <- dplyr::bind_rows(losses)
  utils::write.table(per_og_tbl, file.path(out_dir, "dtl_per_og.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(events_tbl, file.path(out_dir, "dtl_events.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  if (nrow(losses_tbl))
    utils::write.table(losses_tbl, file.path(out_dir, "dtl_losses.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  if (verbose) message(sprintf("  [dtl] %d OGs processed; S=%d D=%d T=%d loss=%d",
                               nrow(per_og_tbl),
                               sum(per_og_tbl$n_S),
                               sum(per_og_tbl$n_D),
                               sum(per_og_tbl$n_T),
                               sum(per_og_tbl$n_loss)))
  list(per_og = per_og_tbl, events = events_tbl, losses = losses_tbl)
}
