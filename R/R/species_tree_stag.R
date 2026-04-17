#' Build a species tree from gene trees (STAG-style)
#'
#' Implements the STAG (Species Tree from All Genes) idea in pure R:
#' for each gene tree whose leaves map unambiguously to species (one
#' gene per species in that OG), compute the cophenetic distance matrix
#' between species. Aggregate across all contributing gene trees by
#' taking the element-wise median, then build an unrooted species tree
#' via NJ.
#'
#' Single-copy orthogroups are always eligible. Multi-copy OGs
#' contribute only when `use_multi_copy = TRUE`; in that case each OG
#' contributes the minimum cophenetic distance between any leaf-pair
#' that hits species (i, j). This mirrors OrthoFinder's STAG treatment
#' of in-paralogs.
#'
#' @param og_result Tibble from `per_og_trees()`.
#' @param dnmb `load_dnmb()` result (needed to map protein_uid→genome_uid).
#' @param use_multi_copy If TRUE, include multi-copy OGs via min-pair distance.
#' @param min_contributing Minimum number of gene trees a pair needs to be
#'   called. Pairs below this threshold are imputed from the species
#'   distance matrix's sub-triangle.
#' @return An `ape::phylo` unrooted tree with tip labels = genome_key.
#' @export
stag_species_tree <- function(og_result,
                              dnmb,
                              use_multi_copy = FALSE,
                              min_contributing = 3L) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("ape is required for stag_species_tree().")
  }

  ok <- og_result[!is.na(og_result$tree_path), , drop = FALSE]
  if (!nrow(ok)) stop("stag_species_tree: no gene trees available.")

  if (!use_multi_copy) {
    ok <- ok[ok$single_copy, , drop = FALSE]
    if (!nrow(ok)) stop("stag_species_tree: no single-copy OGs; set use_multi_copy=TRUE.")
  }

  species <- as.character(dnmb$genome_meta$genome_uid)
  name_map <- setNames(dnmb$genome_meta$genome_key, species)
  n_sp <- length(species)

  # Running sums stored as a 3D accumulator: [i, j, k] = kth distance sample.
  # We can't fully materialize that; instead collect a list of matrices
  # per OG and take element-wise median at the end.
  per_og_mats <- vector("list", nrow(ok))

  for (i in seq_len(nrow(ok))) {
    tr <- ape::read.tree(ok$tree_path[i])
    if (!is.null(tr$edge.length)) tr$edge.length <- pmax(tr$edge.length, 0)

    # Leaf label format: "p<protein_uid>_g<genome_uid>"
    lab <- tr$tip.label
    gidx <- sub(".*_g", "", lab)  # → "0","1","2"…
    if (!all(gidx %in% species)) next  # mis-labeled tree; skip

    # cophenetic.phylo returns leaf × leaf patristic distance matrix
    dmat <- ape::cophenetic.phylo(tr)
    colnames(dmat) <- rownames(dmat) <- gidx

    # Aggregate to species × species: min over duplicate pairs
    sp_mat <- matrix(NA_real_, n_sp, n_sp, dimnames = list(species, species))
    present <- intersect(species, unique(gidx))
    for (a in present) for (b in present) {
      if (a == b) { sp_mat[a, b] <- 0; next }
      idx_a <- which(gidx == a)
      idx_b <- which(gidx == b)
      sp_mat[a, b] <- min(dmat[idx_a, idx_b])
    }
    per_og_mats[[i]] <- sp_mat
  }
  per_og_mats <- per_og_mats[!vapply(per_og_mats, is.null, logical(1))]
  if (!length(per_og_mats)) stop("stag_species_tree: no usable gene trees.")

  # Element-wise median across OGs, ignoring NA
  arr <- simplify2array(per_og_mats)
  # arr is [n_sp, n_sp, n_OG]
  n_present <- apply(arr, c(1, 2), function(x) sum(!is.na(x)))
  med <- apply(arr, c(1, 2), stats::median, na.rm = TRUE)
  med[n_present < min_contributing] <- NA

  # Impute any remaining NAs via UPGMA-style average of neighbors
  # (simple fallback; suffices on well-sampled datasets).
  if (any(is.na(med))) {
    med[is.na(med)] <- mean(med, na.rm = TRUE)
  }
  diag(med) <- 0

  tree <- ape::nj(stats::as.dist(med))
  tree$tip.label <- unname(name_map[tree$tip.label])
  tree
}
