#' Cut Hierarchical Orthogroups (HOGs) from gene trees + a rooted species tree
#'
#' Implements duplication-aware HOG inference in the style of
#' OrthoFinder. For each internal species-tree node N and each gene
#' tree G:
#' \enumerate{
#'   \item Annotate every internal G-node as a duplication (D) or a
#'         speciation (S) using the species-overlap heuristic — a node
#'         is D iff its two child subtrees share at least one species.
#'   \item Reconcile every G-node to the species-tree LCA of its leaf
#'         species.
#'   \item A HOG at N is a maximal subtree of G such that (a) all its
#'         leaves belong to genomes under N and (b) no duplication
#'         G-node inside it reconciles to a species-tree node at or
#'         above N. Equivalently, cut every such duplication edge and
#'         emit each resulting component as a HOG.
#' }
#'
#' This is the difference from the 0.0.2 implementation, which only
#' restricted the flat OG table to descendant-genome members and
#' ignored gene-tree duplication events entirely.
#'
#' Output layout under `out_dir`:
#' \preformatted{
#'   HOGs/
#'     N001_<clade>.tsv   # per-internal-node HOG table
#'     ...
#' }
#'
#' Each HOG table row:
#' `hog_id, cluster_id, n_members, member_leaves, member_genome_keys`.
#'
#' @param rooted_species_tree `ape::phylo` rooted tree.
#' @param og_result Tibble from `per_og_trees()`; must have column
#'   `tree_path` pointing at per-OG gene trees.
#' @param dnmb `load_dnmb()` result.
#' @param out_dir Directory under which `HOGs/` is written.
#' @param min_members Drop HOGs with fewer than this many leaves.
#' @param min_dup_support Minimum bootstrap support (0–100) required to
#'   treat a gene-tree node as a duplication. Requires per-OG bootstrap
#'   replicates to have populated `tree$node.label`. Default 0.
#' @return Tibble with columns `node_id, clade_label, n_leaves,
#'   n_hogs, n_members, path`.
#' @export
cut_hogs <- function(rooted_species_tree, og_result, dnmb, out_dir,
                     min_members = 2L,
                     min_dup_support = 0) {
  if (!requireNamespace("ape", quietly = TRUE)) stop("ape required.")
  dir.create(file.path(out_dir, "HOGs"), recursive = TRUE, showWarnings = FALSE)

  tr_sp  <- rooted_species_tree
  n_tip  <- length(tr_sp$tip.label)
  n_int  <- tr_sp$Nnode
  internal <- (n_tip + 1L):(n_tip + n_int)

  desc_keys_by_node <- lapply(internal, function(nd)
    tr_sp$tip.label[.descendants(tr_sp, nd)])
  names(desc_keys_by_node) <- as.character(internal)

  uid2key <- stats::setNames(dnmb$genome_meta$genome_key,
                      as.character(dnmb$genome_meta$genome_uid))

  ok <- og_result[!is.na(og_result$tree_path), , drop = FALSE]
  if (!nrow(ok)) {
    warning("cut_hogs: no gene trees in og_result; returning empty.")
    return(tibble::tibble(node_id = integer(), clade_label = character(),
                          n_leaves = integer(), n_hogs = integer(),
                          n_members = integer(), path = character()))
  }

  # Per-node accumulator: list of {cluster_id, leaves, genome_keys}
  per_node <- vector("list", length(internal))
  names(per_node) <- as.character(internal)

  for (i in seq_len(nrow(ok))) {
    tr <- tryCatch(ape::read.tree(ok$tree_path[i]),
                   error = function(e) NULL)
    if (is.null(tr) || length(tr$tip.label) < 2) next
    if (!is.null(tr$edge.length)) tr$edge.length <- pmax(tr$edge.length, 0)

    gtips <- tr$tip.label
    gidx  <- sub(".*_g", "", gtips)
    gkeys <- unname(uid2key[gidx])
    keep  <- !is.na(gkeys)
    if (!any(keep)) next
    if (!all(keep)) {
      tr <- ape::keep.tip(tr, gtips[keep])
      gtips <- tr$tip.label
      gidx  <- sub(".*_g", "", gtips)
      gkeys <- unname(uid2key[gidx])
    }
    if (length(gtips) < 2) next

    # per_og_trees() produces unrooted NJ/ML gene trees. HOG logic
    # requires a rooted traversal: midpoint-root each gene tree up
    # front. (OrthoFinder uses STRIDE-style rooting of gene trees too;
    # midpoint is a pragmatic pure-R default that closes the bug.)
    tr <- .hog_root_gene_tree(tr)
    # Rooting can reorder tip.label; recompute the species map.
    gtips <- tr$tip.label
    gidx  <- sub(".*_g", "", gtips)
    gkeys <- unname(uid2key[gidx])

    ann <- .annotate_gene_tree(tr, gkeys, tr_sp,
                                min_support = min_dup_support)

    cid <- ok$cluster_id[i]
    for (N in internal) {
      desc_keys <- desc_keys_by_node[[as.character(N)]]
      hog_roots <- .find_hog_roots(tr, ann, N, desc_keys, tr_sp)
      if (!length(hog_roots)) next
      for (hr in hog_roots) {
        members <- if (hr <= length(gtips)) hr else .descendant_tips_hog(tr, hr)
        if (length(members) < min_members) next
        leaves <- gtips[members]
        keys   <- gkeys[members]
        per_node[[as.character(N)]][[length(per_node[[as.character(N)]]) + 1]] <- list(
          cluster_id   = cid,
          leaves       = leaves,
          genome_keys  = keys
        )
      }
    }
  }

  rows <- list()
  for (N in internal) {
    hogs <- per_node[[as.character(N)]]
    if (!length(hogs)) next
    desc_keys <- desc_keys_by_node[[as.character(N)]]
    clade_label <- if (N == n_tip + 1L) "root" else paste(utils::head(desc_keys, 3), collapse = "+")
    safe_label  <- gsub("[^A-Za-z0-9._-]+", "_", clade_label)
    path <- file.path(out_dir, "HOGs",
                      sprintf("N%03d_%s.tsv", N - n_tip, safe_label))
    con <- file(path, "w")
    cat("# node_id=", N, "\n", sep = "", file = con)
    cat("# clade_leaves=", paste(desc_keys, collapse = ","), "\n",
        sep = "", file = con)
    cat("hog_id\tcluster_id\tn_members\tmember_leaves\tmember_genome_keys\n",
        file = con)
    for (j in seq_along(hogs)) {
      h <- hogs[[j]]
      hog_id <- sprintf("N%03d.HOG%05d", N - n_tip, j)
      cat(hog_id, "\t",
          h$cluster_id, "\t",
          length(h$leaves), "\t",
          paste(h$leaves, collapse = ","), "\t",
          paste(h$genome_keys, collapse = ","), "\n",
          sep = "", file = con)
    }
    close(con)
    rows[[length(rows) + 1]] <- tibble::tibble(
      node_id      = N - n_tip,
      clade_label  = clade_label,
      n_leaves     = length(desc_keys),
      n_hogs       = length(hogs),
      n_members    = sum(vapply(hogs, function(h) length(h$leaves), integer(1))),
      path         = path
    )
  }
  dplyr::bind_rows(rows)
}


# ---- helpers ---------------------------------------------------------------

.hog_root_gene_tree <- function(tr) {
  if (!is.null(attr(tr, "order")) && isTRUE(ape::is.rooted(tr))) return(tr)
  if (requireNamespace("phangorn", quietly = TRUE)) {
    return(tryCatch(phangorn::midpoint(tr), error = function(e) .manual_midpoint(tr)))
  }
  .manual_midpoint(tr)
}

.manual_midpoint <- function(tr) {
  d <- ape::cophenetic.phylo(tr)
  mx <- which(d == max(d), arr.ind = TRUE)[1, ]
  tip_a <- rownames(d)[mx[1]]
  ape::root(tr, outgroup = tip_a, resolve.root = TRUE)
}

.descendants <- function(tr, node) {
  n_tip <- length(tr$tip.label)
  if (node <= n_tip) return(node)
  stack <- node
  tips <- integer()
  while (length(stack)) {
    cur <- stack[1]; stack <- stack[-1]
    ch <- tr$edge[tr$edge[, 1] == cur, 2]
    leaf <- ch[ch <= n_tip]
    tips <- c(tips, leaf)
    stack <- c(stack, ch[ch > n_tip])
  }
  sort(unique(tips))
}

.descendant_tips_hog <- function(tr, node) {
  n_tip <- length(tr$tip.label)
  if (node <= n_tip) return(node)
  stack <- node
  tips <- integer()
  while (length(stack)) {
    cur <- stack[1]; stack <- stack[-1]
    ch <- tr$edge[tr$edge[, 1] == cur, 2]
    leaf <- ch[ch <= n_tip]
    tips <- c(tips, leaf)
    stack <- c(stack, ch[ch > n_tip])
  }
  sort(unique(tips))
}

# Annotate a rooted gene tree with: children per node, post-order,
# species set per node, is_duplication, species-tree LCA node.
.annotate_gene_tree <- function(tr, gkeys, tr_sp, min_support = 0) {
  n_tip <- length(tr$tip.label)
  n_int <- tr$Nnode
  total <- n_tip + n_int

  nlab <- tr$node.label
  support_vec <- rep(NA_real_, total)
  if (!is.null(nlab) && length(nlab) == n_int) {
    s <- suppressWarnings(as.numeric(nlab))
    support_vec[(n_tip + 1L):total] <- s
  }

  children <- vector("list", total)
  parents  <- tr$edge[, 1]
  childv   <- tr$edge[, 2]
  for (e in seq_along(parents)) {
    p <- parents[e]
    children[[p]] <- c(children[[p]], childv[e])
  }

  sp_set <- vector("list", total)
  for (v in seq_len(n_tip)) sp_set[[v]] <- gkeys[v]

  # Iterative post-order starting at root (n_tip + 1 for rooted phylo)
  root <- n_tip + 1L
  order_vec <- integer(total)
  ptr <- 0L
  stack <- list(list(node = root, state = 0L))
  while (length(stack)) {
    top <- stack[[length(stack)]]
    v <- top$node
    ch <- children[[v]]
    if (top$state == 0L) {
      stack[[length(stack)]]$state <- 1L
      for (c in ch) stack[[length(stack) + 1L]] <- list(node = c, state = 0L)
    } else {
      ptr <- ptr + 1L
      order_vec[ptr] <- v
      stack[[length(stack)]] <- NULL
    }
  }
  order_vec <- order_vec[seq_len(ptr)]

  for (v in order_vec) {
    if (v <= n_tip) next
    ch <- children[[v]]
    sp_set[[v]] <- unique(unlist(sp_set[ch]))
  }

  is_dup <- logical(total)
  for (v in order_vec) {
    if (v <= n_tip) next
    ch <- children[[v]]
    if (length(ch) < 2) next
    if (min_support > 0) {
      s <- support_vec[v]
      if (!is.na(s) && s < min_support) next
    }
    any_ov <- FALSE
    for (a in seq_along(ch)) for (b in seq_len(a - 1L)) {
      if (length(intersect(sp_set[[ch[a]]], sp_set[[ch[b]]]))) {
        any_ov <- TRUE; break
      }
      if (any_ov) break
    }
    is_dup[v] <- any_ov
  }

  # Species-tree LCA reconciliation per node
  lca_node <- integer(total)
  n_sp_tip <- length(tr_sp$tip.label)
  for (v in seq_len(total)) {
    sp <- sp_set[[v]]
    sp <- sp[sp %in% tr_sp$tip.label]
    if (!length(sp)) { lca_node[v] <- NA_integer_; next }
    if (length(sp) == 1L) {
      lca_node[v] <- which(tr_sp$tip.label == sp)[1]
    } else {
      lca_node[v] <- ape::getMRCA(tr_sp, sp)
    }
  }

  list(children = children, sp_set = sp_set, is_dup = is_dup,
       lca_node = lca_node, root = root, n_tip = n_tip)
}

# Given species-tree node N and gene-tree annotation, return HOG roots
# at level N: maximal gene-tree nodes whose descendant species all lie
# under N and whose enclosing structure is not a duplication at/above N.
.find_hog_roots <- function(tr, ann, N, desc_keys, tr_sp) {
  children <- ann$children
  sp_set   <- ann$sp_set
  is_dup   <- ann$is_dup
  lca_node <- ann$lca_node
  root     <- ann$root

  # Ancestor-or-equal test in species tree: is lca_v at or above N?
  # lca_v is at/above N iff the set of species under lca_v ⊇ desc_keys(N).
  n_sp_tip <- length(tr_sp$tip.label)
  anc_or_eq_to_N <- function(lca_v) {
    if (is.na(lca_v)) return(FALSE)
    if (lca_v <= n_sp_tip) {
      # leaf node — only ancestor-or-equal to N if N is that same leaf
      return(lca_v %in% which(tr_sp$tip.label %in% desc_keys) && length(desc_keys) == 1L)
    }
    lca_desc <- tr_sp$tip.label[.descendants(tr_sp, lca_v)]
    all(desc_keys %in% lca_desc)
  }

  hog_roots <- integer()
  stack <- root
  while (length(stack)) {
    v <- stack[length(stack)]; stack <- stack[-length(stack)]
    sp <- sp_set[[v]]
    any_under <- any(sp %in% desc_keys)
    if (!any_under) next
    all_under <- all(sp %in% desc_keys)
    if (!all_under) {
      stack <- c(stack, children[[v]])
      next
    }
    if (is_dup[v] && length(children[[v]]) > 0 && anc_or_eq_to_N(lca_node[v])) {
      # Duplication reconciles at/above N — cut here, descend.
      stack <- c(stack, children[[v]])
      next
    }
    hog_roots <- c(hog_roots, v)
  }
  hog_roots
}
