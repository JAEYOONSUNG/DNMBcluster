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
#' @param method Tree method: `"nj"` (fast) or `"ml"` (phangorn NNI).
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
#' @return Named list with keys `og_result`, `species_tree`, `hogs`,
#'   `relationships`, and (if refine_graph) `refined`.
#' @export
run_orthofinder_like <- function(dnmb,
                                 out_dir,
                                 method = c("nj", "ml"),
                                 min_seqs = 3L,
                                 max_seqs = 500L,
                                 use_multi_copy = FALSE,
                                 outgroup = NULL,
                                 xeno_z = 3.0,
                                 threads = 2L,
                                 verbose = TRUE,
                                 refine_graph = FALSE,
                                 refine_method = c("louvain", "mcl", "walktrap"),
                                 bootstrap = 0L) {
  refine_method <- match.arg(refine_method)
  method <- match.arg(method)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  t_start <- Sys.time()
  refined <- NULL
  if (refine_graph) {
    if (verbose) message("[0/5] Refining OGs via length-normalized graph…")
    refined <- refine_ogs_graph(dnmb, out_dir = out_dir,
                                 method = refine_method, threads = threads)
    # Rewrite dnmb$clusters so downstream phases see the split OGs.
    # refined_id is character → hash to new integer cluster ids.
    lvl <- unique(refined$refined_id)
    new_cid <- setNames(seq_along(lvl), lvl)
    remap <- tibble::tibble(
      protein_uid = refined$protein_uid,
      new_cluster = as.integer(new_cid[refined$refined_id])
    )
    dnmb$clusters <- dnmb$clusters %>%
      dplyr::left_join(remap, by = "protein_uid") %>%
      dplyr::mutate(cluster_id = dplyr::coalesce(new_cluster, cluster_id)) %>%
      dplyr::select(-new_cluster)
  }

  if (verbose) message(sprintf("[%d/%d] Building per-OG trees…",
                                if (refine_graph) 1L else 1L,
                                if (refine_graph) 5L else 4L))
  og_result <- per_og_trees(dnmb, out_dir,
                            min_seqs = min_seqs, max_seqs = max_seqs,
                            method = method, threads = threads,
                            verbose = verbose)

  if (verbose) message("[2] STAG + STRIDE species tree…")
  n_sc <- sum(og_result$single_copy, na.rm = TRUE)
  effective_multi <- use_multi_copy || n_sc == 0
  if (effective_multi && !use_multi_copy && verbose) {
    message("  no single-copy OGs present → falling back to use_multi_copy=TRUE")
  }
  sp_unrooted <- stag_species_tree(og_result, dnmb,
                                    use_multi_copy = effective_multi,
                                    bootstrap = bootstrap)
  species_tree <- stride_root(sp_unrooted, og_result, dnmb, outgroup = outgroup)
  ape::write.tree(species_tree, file.path(out_dir, "species_tree_rooted.nwk"))

  if (verbose) message("[3/4] Cutting HOGs…")
  hogs <- cut_hogs(species_tree, dnmb, out_dir)

  if (verbose) message("[4/4] Reconciling ortholog/paralog/xenolog labels…")
  relationships <- reconcile_relationships(og_result, species_tree, dnmb,
                                            out_dir = out_dir,
                                            xeno_z = xeno_z)

  if (verbose) {
    elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)
    message(sprintf("[run_orthofinder_like] done in %ss. OGs=%d, SC-OGs=%d, HOG-nodes=%d, relationships=%d",
                    elapsed,
                    sum(!is.na(og_result$tree_path)),
                    sum(og_result$single_copy),
                    nrow(hogs),
                    nrow(relationships)))
  }

  invisible(list(
    og_result     = og_result,
    species_tree  = species_tree,
    hogs          = hogs,
    relationships = relationships,
    refined       = refined
  ))
}
