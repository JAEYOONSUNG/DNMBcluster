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
#' @return A rooted `ape::phylo` with the selected root.
#' @export
stride_root <- function(species_tree, og_result, dnmb, outgroup = NULL) {
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
    parts <- .detect_duplications(tr, dnmb)
    if (length(parts)) dup_partitions <- c(dup_partitions, parts)
  }

  if (!length(dup_partitions)) {
    message("[stride_root] No duplications detected; using midpoint rooting.")
    return(.midpoint_root(species_tree))
  }

  # Score each edge by fraction of duplication partitions that are
  # compatible with rooting on that edge.
  scores <- .score_edges(species_tree, dup_partitions)
  best_edge <- which.max(scores)
  rooted <- .root_on_edge(species_tree, best_edge)
  rooted
}


# ---- helpers ---------------------------------------------------------------

.detect_duplications <- function(gene_tree, dnmb) {
  tips <- gene_tree$tip.label
  gidx <- sub(".*_g", "", tips)
  if (length(tips) < 4) return(list())

  name_map <- setNames(dnmb$genome_meta$genome_key, as.character(dnmb$genome_meta$genome_uid))
  # for each internal node, get the leaf sets of each child subtree
  n_int <- gene_tree$Nnode
  n_tip <- length(tips)
  partitions <- list()

  for (nd in (n_tip + 1):(n_tip + n_int)) {
    children <- gene_tree$edge[gene_tree$edge[, 1] == nd, 2]
    if (length(children) != 2) next
    sp_sets <- lapply(children, function(ch) {
      if (ch <= n_tip) {
        unique(gidx[ch])
      } else {
        desc <- .descendant_tips(gene_tree, ch)
        unique(gidx[desc])
      }
    })
    shared <- intersect(sp_sets[[1]], sp_sets[[2]])
    # Duplication if the two subtrees share species
    if (length(shared) >= 2) {
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
  if (requireNamespace("phangorn", quietly = TRUE)) {
    return(phangorn::midpoint(tr))
  }
  # Manual midpoint fallback
  d <- ape::cophenetic.phylo(tr)
  mx <- which(d == max(d), arr.ind = TRUE)[1, ]
  tip_a <- rownames(d)[mx[1]]
  ape::root(tr, outgroup = tip_a, resolve.root = TRUE)
}
