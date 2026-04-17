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
#' @return Named list with keys `og_result`, `species_tree`, `hogs`,
#'   `relationships`. All phase outputs are also written to `out_dir`.
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
                                 verbose = TRUE) {
  method <- match.arg(method)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  t_start <- Sys.time()
  if (verbose) message("[1/4] Building per-OG trees…")
  og_result <- per_og_trees(dnmb, out_dir,
                            min_seqs = min_seqs, max_seqs = max_seqs,
                            method = method, threads = threads,
                            verbose = verbose)

  if (verbose) message("[2/4] STAG + STRIDE species tree…")
  n_sc <- sum(og_result$single_copy, na.rm = TRUE)
  effective_multi <- use_multi_copy || n_sc == 0
  if (effective_multi && !use_multi_copy && verbose) {
    message("  no single-copy OGs present → falling back to use_multi_copy=TRUE")
  }
  sp_unrooted <- stag_species_tree(og_result, dnmb,
                                    use_multi_copy = effective_multi)
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
    relationships = relationships
  ))
}
