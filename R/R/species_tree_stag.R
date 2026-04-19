#' Build a species tree from gene trees (STAG-style, with bootstrap)
#'
#' Implements STAG (Species Tree from All Genes) in pure R. For each
#' gene tree, computes the cophenetic distance between species present
#' in that tree. Aggregates across all contributing gene trees and
#' builds a species tree via NJ.
#'
#' Features added beyond the 0.0.2 MVP, to close the gap with
#' OrthoFinder-3's STAG:
#' \itemize{
#'   \item \strong{Weighted aggregation}. Each OG can contribute with a
#'         weight equal to `sqrt(n_members)` (coverage signal without
#'         letting paralog-heavy families dominate — pure n_members
#'         caused huge over-representation). Disable with
#'         `weighted = FALSE` for unweighted median.
#'   \item \strong{Outlier trimming}. For each species pair, drop
#'         distance samples outside `[Q1 - k·IQR, Q3 + k·IQR]` before
#'         aggregation. `trim_iqr = 1.5` (default Tukey rule) matches
#'         STAG's outlier filter.
#'   \item \strong{Bootstrap support}. `bootstrap = N` resamples the OG
#'         set with replacement N times and rebuilds the tree per
#'         replicate. The final tree's internal-node labels are the
#'         percentage of bootstrap trees in which the corresponding
#'         bipartition appears.
#' }
#'
#' @param og_result Tibble from `per_og_trees()`.
#' @param dnmb `load_dnmb()` result.
#' @param use_multi_copy If TRUE, include multi-copy OGs via min-pair distance.
#' @param min_contributing Minimum gene trees a pair needs to be called.
#' @param weighted Logical; weight each OG's contribution by n_members.
#' @param trim_iqr Numeric k for IQR outlier rejection (1.5 = Tukey);
#'   set to NULL or 0 to disable.
#' @param bootstrap Integer; number of bootstrap replicates for node support.
#'   0 disables.
#' @param bootstrap_seed RNG seed for reproducibility.
#' @return An `ape::phylo` unrooted tree with tip labels = genome_key.
#'   If `bootstrap > 0`, `$node.label` carries bipartition support (0-100).
#' @export
stag_species_tree <- function(og_result,
                              dnmb,
                              use_multi_copy  = FALSE,
                              min_contributing = 3L,
                              weighted        = TRUE,
                              trim_iqr        = 1.5,
                              bootstrap       = 0L,
                              bootstrap_seed  = 42L) {
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
  name_map <- stats::setNames(dnmb$genome_meta$genome_key, species)

  # Collect per-OG distance matrices + weights (member count)
  per_og <- .collect_og_distances(ok, species)

  # Point estimate tree from full set
  main_tree <- .aggregate_and_nj(per_og$mats, per_og$weights, species,
                                  name_map,
                                  min_contributing = min_contributing,
                                  weighted = weighted,
                                  trim_iqr = trim_iqr)

  if (bootstrap > 0L) {
    set.seed(bootstrap_seed)
    n <- length(per_og$mats)
    reps <- vector("list", bootstrap)
    for (b in seq_len(bootstrap)) {
      idx <- sample.int(n, n, replace = TRUE)
      reps[[b]] <- .aggregate_and_nj(per_og$mats[idx], per_og$weights[idx],
                                      species, name_map,
                                      min_contributing = min_contributing,
                                      weighted = weighted,
                                      trim_iqr = trim_iqr)
    }
    # Clade-by-clade support: prop.clades returns #reps each node appears in.
    support <- ape::prop.clades(main_tree, reps, rooted = FALSE)
    support <- round(100 * support / bootstrap)
    main_tree$node.label <- support
  }
  main_tree
}


# ---- internal helpers ------------------------------------------------------

.collect_og_distances <- function(ok, species) {
  n_sp <- length(species)
  mats <- vector("list", nrow(ok))
  weights <- numeric(nrow(ok))
  k <- 0L
  for (i in seq_len(nrow(ok))) {
    tr <- ape::read.tree(ok$tree_path[i])
    if (!is.null(tr$edge.length)) tr$edge.length <- pmax(tr$edge.length, 0)
    lab <- tr$tip.label
    gidx <- sub(".*_g", "", lab)
    if (!all(gidx %in% species)) next

    dmat <- ape::cophenetic.phylo(tr)
    colnames(dmat) <- rownames(dmat) <- gidx

    sp_mat <- matrix(NA_real_, n_sp, n_sp, dimnames = list(species, species))
    present <- intersect(species, unique(gidx))
    for (a in present) for (b in present) {
      if (a == b) { sp_mat[a, b] <- 0; next }
      ia <- which(gidx == a); ib <- which(gidx == b)
      sp_mat[a, b] <- min(dmat[ia, ib])
    }
    k <- k + 1L
    mats[[k]] <- sp_mat
    # sqrt dampens paralog-rich OGs; a 100-member fusion family would
    # otherwise outweigh 25 single-copy OGs combined.
    weights[k] <- sqrt(ok$n_members[i])
  }
  list(mats = mats[seq_len(k)], weights = weights[seq_len(k)])
}

.aggregate_and_nj <- function(mats, weights, species, name_map,
                              min_contributing, weighted, trim_iqr) {
  if (!length(mats)) stop(".aggregate_and_nj: empty mats.")
  arr <- simplify2array(mats)
  if (length(dim(arr)) == 2) {
    dim(arr) <- c(dim(arr), 1L)
    dimnames(arr) <- c(dimnames(mats[[1]]), list(NULL))
  }

  n_present <- apply(arr, c(1, 2), function(x) sum(!is.na(x)))

  # Auto-scale the contributing-tree threshold: when the caller asks for
  # e.g. 3 but only 2 gene trees exist, every pair would fail the gate
  # and every species would be dropped. Cap the threshold at the number
  # of trees so small OG counts still produce a tree.
  n_trees <- length(mats)
  effective_min <- min(as.integer(min_contributing), as.integer(n_trees))
  if (effective_min < 1L) effective_min <- 1L

  # Aggregate per cell: trim_iqr → weighted median
  med <- apply(arr, c(1, 2), function(samples) {
    keep <- !is.na(samples)
    if (!any(keep)) return(NA_real_)
    x <- samples[keep]
    w <- weights[keep]
    if (!is.null(trim_iqr) && trim_iqr > 0 && length(x) >= 4L) {
      q <- stats::quantile(x, c(0.25, 0.75), names = FALSE)
      iqr <- q[2] - q[1]
      lo <- q[1] - trim_iqr * iqr
      hi <- q[2] + trim_iqr * iqr
      keep2 <- x >= lo & x <= hi
      if (any(keep2)) { x <- x[keep2]; w <- w[keep2] }
    }
    if (weighted && length(x) >= 1L && sum(w) > 0) {
      .weighted_median(x, w)
    } else {
      stats::median(x)
    }
  })
  med[n_present < effective_min] <- NA
  diag(med) <- 0

  # Drop species that have no callable distance to any other species.
  # This avoids silently imputing fabricated distances with the global
  # mean, which destroys the topology when a genome is only present
  # in a handful of OGs.
  off_diag_na <- is.na(med)
  diag(off_diag_na) <- FALSE
  fully_na <- rowSums(!off_diag_na) == 1L  # only self-comparison left
  if (any(fully_na)) {
    dropped <- rownames(med)[fully_na]
    warning("stag_species_tree: dropping species with no callable pair: ",
            paste(unname(name_map[dropped]), collapse = ", "))
    med <- med[!fully_na, !fully_na, drop = FALSE]
    if (nrow(med) < 3) stop("stag_species_tree: <3 species remain after NA drop.")
  }

  # Any remaining NAs are species pairs that never co-occurred in any
  # OG. NJ cannot consume NAs, so we interpolate each gap from the
  # shortest path through other species in the observed-distance graph
  # (Floyd–Warshall with NA-tolerant additions). This preserves
  # topology far better than replacing with the global mean.
  if (anyNA(med)) {
    med <- .floyd_fill(med)
  }

  tree <- ape::nj(stats::as.dist(med))
  tree$tip.label <- unname(name_map[tree$tip.label])
  tree
}

.floyd_fill <- function(d) {
  n <- nrow(d)
  m <- d
  m[is.na(m)] <- Inf
  for (k in seq_len(n)) {
    m <- pmin(m, outer(m[, k], m[k, ], "+"))
  }
  if (any(is.infinite(m))) {
    stop(".floyd_fill: tree graph disconnected; cannot reconstruct pairwise distances.")
  }
  dimnames(m) <- dimnames(d)
  m
}

.weighted_median <- function(x, w) {
  ord <- order(x)
  x <- x[ord]; w <- w[ord]
  cw <- cumsum(w) / sum(w)
  # Smallest x whose cumulative weight >= 0.5
  x[which(cw >= 0.5)[1]]
}
