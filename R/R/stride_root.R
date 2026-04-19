#' Root a species tree by duplication-consistency (STRIDE-style)
#'
#' For each gene tree containing at least one duplication, the
#' duplication node partitions species into two subsets. A valid root
#' for the species tree cannot lie on any edge that would make the
#' duplication's two children inconsistent with each other — i.e. the
#' two subsets must both be "clades" in the species tree with respect
#' to the candidate root.
#'
#' This implementation scores each edge in the unrooted species tree by
#' the fraction of duplication events it is consistent with. The edge
#' maximizing this score becomes the root. When no duplications are
#' present (or none are resolvable), falls back to midpoint rooting.
#'
#' A duplication is detected in a gene tree when a node's two children
#' each contain at least one overlapping species (the "species overlap"
#' heuristic — same criterion OrthoFinder uses).
#'
#' @param species_tree Unrooted `ape::phylo` (e.g. from `stag_species_tree()`).
#' @param og_result Tibble from `per_og_trees()`.
#' @param dnmb `load_dnmb()` result.
#' @param outgroup Optional character vector of `genome_key` values to
#'   force as the outgroup. Overrides STRIDE scoring.
#' @param min_dup_support Skip duplication evidence from gene-tree nodes
#'   whose bootstrap support (`tree$node.label`, 0–100) is below this
#'   threshold. Requires per-OG bootstrap replicates to have been run.
#'   Default 0 (use every duplication).
#' @return A rooted `ape::phylo` with the selected root.
#' @export
stride_root <- function(species_tree, og_result, dnmb, outgroup = NULL,
                        min_dup_support = 0) {
  if (!requireNamespace("ape", quietly = TRUE)) stop("ape required.")

  if (!is.null(outgroup) && length(outgroup)) {
    hits <- intersect(outgroup, species_tree$tip.label)
    if (length(hits)) return(ape::root(species_tree, outgroup = hits, resolve.root = TRUE))
    warning("stride_root: outgroup not found on tree; falling back to STRIDE.")
  }

  ok <- og_result[!is.na(og_result$tree_path), , drop = FALSE]
  dup_partitions <- list()
  for (i in seq_len(nrow(ok))) {
    tr <- ape::read.tree(ok$tree_path[i])
    if (!is.null(tr$edge.length)) tr$edge.length <- pmax(tr$edge.length, 0)
    parts <- .detect_duplications(tr, dnmb, min_support = min_dup_support)
    if (length(parts)) dup_partitions <- c(dup_partitions, parts)
  }

  if (!length(dup_partitions)) {
    message("[stride_root] No duplications detected; using MAD rooting (with midpoint fallback).")
    return(.mad_root(species_tree))
  }

  # Score each edge by fraction of duplication partitions that are
  # compatible with rooting on that edge.
  scores <- .score_edges(species_tree, dup_partitions)
  best_edge <- which.max(scores)
  rooted <- .root_on_edge(species_tree, best_edge)
  rooted
}


# ---- helpers ---------------------------------------------------------------

.detect_duplications <- function(gene_tree, dnmb, min_support = 0) {
  tips <- gene_tree$tip.label
  gidx <- sub(".*_g", "", tips)
  if (length(tips) < 4) return(list())

  name_map <- stats::setNames(dnmb$genome_meta$genome_key, as.character(dnmb$genome_meta$genome_uid))
  # for each internal node, get the leaf sets of each child subtree
  n_int <- gene_tree$Nnode
  n_tip <- length(tips)
  partitions <- list()

  nlab <- gene_tree$node.label
  have_support <- !is.null(nlab) && length(nlab) == n_int && min_support > 0

  for (nd in (n_tip + 1):(n_tip + n_int)) {
    children <- gene_tree$edge[gene_tree$edge[, 1] == nd, 2]
    if (length(children) != 2) next
    if (have_support) {
      s <- suppressWarnings(as.numeric(nlab[nd - n_tip]))
      if (!is.na(s) && s < min_support) next
    }
    sp_sets <- lapply(children, function(ch) {
      if (ch <= n_tip) {
        unique(gidx[ch])
      } else {
        desc <- .descendant_tips(gene_tree, ch)
        unique(gidx[desc])
      }
    })
    shared <- intersect(sp_sets[[1]], sp_sets[[2]])
    # Species-overlap heuristic (OrthoFinder STRIDE): a node is a
    # duplication iff its two child subtrees share >= 1 species.
    if (length(shared) >= 1) {
      partitions[[length(partitions) + 1]] <- list(
        A = unname(name_map[sp_sets[[1]]]),
        B = unname(name_map[sp_sets[[2]]])
      )
    }
  }
  partitions
}

.descendant_tips <- function(tr, node) {
  n_tip <- length(tr$tip.label)
  stack <- node
  tips <- integer()
  while (length(stack)) {
    cur <- stack[1]; stack <- stack[-1]
    ch <- tr$edge[tr$edge[, 1] == cur, 2]
    leaf <- ch[ch <= n_tip]
    tips <- c(tips, leaf)
    stack <- c(stack, ch[ch > n_tip])
  }
  tips
}

.score_edges <- function(sp_tree, partitions) {
  n_edge <- nrow(sp_tree$edge)
  scores <- numeric(n_edge)
  for (e in seq_len(n_edge)) {
    bip <- .edge_bipartition(sp_tree, e)
    # Count partitions whose {A, B} matches (or subsets) the bipartition
    s <- 0
    for (p in partitions) {
      a_in_1 <- all(p$A %in% bip$side1)
      b_in_2 <- all(p$B %in% bip$side2)
      a_in_2 <- all(p$A %in% bip$side2)
      b_in_1 <- all(p$B %in% bip$side1)
      if ((a_in_1 && b_in_2) || (a_in_2 && b_in_1)) s <- s + 1
    }
    scores[e] <- s
  }
  scores
}

.edge_bipartition <- function(tr, edge_idx) {
  # Split leaves by which side of the edge they fall on
  n_tip <- length(tr$tip.label)
  parent <- tr$edge[edge_idx, 1]
  child  <- tr$edge[edge_idx, 2]
  if (child <= n_tip) {
    side1 <- tr$tip.label[child]
    side2 <- setdiff(tr$tip.label, side1)
  } else {
    desc <- .descendant_tips(tr, child)
    side1 <- tr$tip.label[desc]
    side2 <- setdiff(tr$tip.label, side1)
  }
  list(side1 = side1, side2 = side2)
}

.root_on_edge <- function(tr, edge_idx) {
  # Root by designating the subtree below `edge_idx` as outgroup
  n_tip <- length(tr$tip.label)
  child <- tr$edge[edge_idx, 2]
  if (child <= n_tip) {
    og <- tr$tip.label[child]
  } else {
    og <- tr$tip.label[.descendant_tips(tr, child)]
  }
  ape::root(tr, outgroup = og, resolve.root = TRUE)
}

.midpoint_root <- function(tr) {
  if (is.null(tr$edge.length) || !length(tr$edge.length)) {
    tr$edge.length <- rep(1, nrow(tr$edge))
  }
  if (requireNamespace("phangorn", quietly = TRUE)) {
    return(phangorn::midpoint(tr))
  }
  d <- ape::cophenetic.phylo(tr)
  mx <- which(d == max(d), arr.ind = TRUE)[1, ]
  tip_a <- rownames(d)[mx[1]]
  ape::root(tr, outgroup = tip_a, resolve.root = TRUE)
}

#' Root a tree using Minimum Ancestor Deviation (MAD)
#'
#' Implements the MAD algorithm of Tria, Landan & Dagan (2017) for
#' rootless binary species trees. For every internal edge, MAD finds
#' the root position that minimizes a weighted root-to-tip deviation
#' across tip pairs spanning the edge, then picks the edge with the
#' smallest residual. Falls back to midpoint rooting when edge lengths
#' are missing or degenerate.
#'
#' Only spanning pairs are scored — same-side pairs contribute a
#' constant that does not affect the argmin across edges.
#'
#' @param tr Unrooted `ape::phylo`.
#' @return Rooted `ape::phylo` with the MAD-optimal edge as root.
#' @export
mad_root <- function(tr) .mad_root(tr)

.mad_root <- function(tr) {
  if (!requireNamespace("ape", quietly = TRUE)) stop("ape required.")
  if (is.null(tr$edge.length) || !length(tr$edge.length) ||
      all(is.na(tr$edge.length)) || max(tr$edge.length, na.rm = TRUE) <= 0) {
    return(.midpoint_root(tr))
  }
  n_tip <- length(tr$tip.label)
  DN <- tryCatch(ape::dist.nodes(tr), error = function(e) NULL)
  if (is.null(DN)) return(.midpoint_root(tr))

  best_score <- Inf
  best_edge  <- NA_integer_
  best_rho   <- NA_real_
  tip_seq    <- seq_len(n_tip)

  for (e in seq_len(nrow(tr$edge))) {
    u <- tr$edge[e, 1]; v <- tr$edge[e, 2]
    L <- tr$edge.length[e]
    if (is.na(L) || L <= 0) next

    v_tips <- .descendant_tips(tr, v)
    if (v <= n_tip) v_tips <- v
    u_tips <- setdiff(tip_seq, v_tips)
    if (!length(v_tips) || !length(u_tips)) next

    d_ui <- DN[u, u_tips]
    d_vj <- DN[v, v_tips]
    c_mat <- outer(d_ui, d_vj, function(a, b) L + b - a)
    d_mat <- outer(d_ui, d_vj, function(a, b) a + L + b)
    d_mat[d_mat == 0] <- NA_real_   # self-pairs shouldn't occur but guard
    w <- 1 / d_mat^2
    w[is.na(w)] <- 0
    denom <- sum(w)
    if (denom == 0) next
    rho_star <- 0.5 * sum(c_mat * w, na.rm = TRUE) / denom
    rho_star <- max(0, min(L, rho_star))
    dev <- (2 * rho_star - c_mat) / d_mat
    rms <- sqrt(mean(dev^2, na.rm = TRUE))
    if (!is.finite(rms)) next
    if (rms < best_score) {
      best_score <- rms; best_edge <- e; best_rho <- rho_star
    }
  }

  if (is.na(best_edge)) return(.midpoint_root(tr))
  .root_on_edge(tr, best_edge)
}
