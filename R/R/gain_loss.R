#' Build a species × cluster presence/absence matrix from a DNMB result
#'
#' Rows are species (ordered to match `species_tree$tip.label` when
#' supplied), columns are cluster IDs. Values are 0/1.
#'
#' @param dnmb Result of `load_dnmb()`.
#' @param species_tree Optional rooted `ape::phylo` used to fix row order.
#' @return Integer matrix with species on rows, clusters on columns.
#' @export
presence_matrix_from_dnmb <- function(dnmb, species_tree = NULL) {
  stopifnot(is.list(dnmb), "clusters" %in% names(dnmb),
            "genome_meta" %in% names(dnmb))
  uid2key <- stats::setNames(dnmb$genome_meta$genome_key,
                      as.character(dnmb$genome_meta$genome_uid))
  tab <- dnmb$clusters %>%
    dplyr::distinct(cluster_id, genome_uid)
  tab$genome_key <- unname(uid2key[as.character(tab$genome_uid)])
  tab <- tab[!is.na(tab$genome_key), , drop = FALSE]

  sp_levels <- if (!is.null(species_tree)) species_tree$tip.label
                else sort(unique(tab$genome_key))
  cl_levels <- sort(unique(tab$cluster_id))

  m <- matrix(0L, nrow = length(sp_levels), ncol = length(cl_levels),
              dimnames = list(sp_levels, as.character(cl_levels)))
  ri <- match(tab$genome_key, sp_levels)
  ci <- match(tab$cluster_id, cl_levels)
  ok <- !is.na(ri) & !is.na(ci)
  m[cbind(ri[ok], ci[ok])] <- 1L
  m
}

#' Reconstruct gene gains and losses along a rooted species tree
#'
#' Given a rooted species tree and a species × HOG (or cluster) presence
#' matrix, infers the ancestral state at every internal node and counts
#' gains/losses on every branch. Three algorithms are supported:
#'
#' \describe{
#'   \item{\code{"dollo"}}{Dollo parsimony: one gain per HOG (at the
#'     MRCA of presence tips), unlimited losses on lineages lacking it.}
#'   \item{\code{"fitch"}}{Binary Fitch parsimony: minimum state changes
#'     without the "no reversal" constraint.}
#'   \item{\code{"ml"}}{Per-HOG ML via \code{ape::ace} with an ARD
#'     two-state model; ancestral state = 1 if P(state = 1) >=
#'     \code{min_ml_prob}.}
#' }
#'
#' Dollo and Fitch are computed vectorised across all HOG columns at
#' once. For ML, identical presence columns are de-duplicated before
#' calling \code{ape::ace}, and — if \code{workers > 1} and
#' \code{future.apply} is installed — unique patterns are processed in
#' parallel.
#'
#' @param species_tree Rooted `ape::phylo`.
#' @param presence Integer/logical matrix (species on rows, HOGs on cols).
#'   Row names must be species matching `species_tree$tip.label`.
#' @param method One of `"dollo"`, `"fitch"`, `"ml"`.
#' @param min_ml_prob ML threshold for state = 1 call. Default 0.5.
#' @param gain_weights Optional numeric vector (length = `ncol(presence)`).
#' @param return_states If TRUE, also return the full node × HOG state
#'   matrix (needed for `ancestral_pan_core_sizes`).
#' @param workers If > 1 and `future.apply` is installed, run ML ace
#'   calls in parallel on unique presence patterns.
#' @return A list with class `dnmb_gain_loss`, fields `method`, `edges`,
#'   `tree`, `n_hogs`, optionally `states`. Gains at the root are silent
#'   (no edge to attribute them to) — HOGs whose MRCA is the root
#'   contribute 0 to the gains tally. Use Fitch or ML if the pangenome
#'   contains many such families.
#' @export
reconstruct_gain_loss <- function(species_tree,
                                   presence,
                                   method = c("dollo", "fitch", "ml"),
                                   min_ml_prob = 0.5,
                                   gain_weights = NULL,
                                   return_states = FALSE,
                                   workers = 1L) {
  method <- match.arg(method)
  if (!requireNamespace("ape", quietly = TRUE)) stop("ape required.")
  tr <- species_tree
  if (!ape::is.rooted(tr)) stop("species_tree must be rooted.")
  if (!is.matrix(presence)) stop("presence must be a species x HOG matrix.")
  if (is.null(rownames(presence))) stop("presence needs species rownames.")

  n_tip <- length(tr$tip.label)
  n_hog <- ncol(presence)

  aligned <- matrix(0L, nrow = n_tip, ncol = n_hog,
                    dimnames = list(tr$tip.label, colnames(presence)))
  tip_order <- match(tr$tip.label, rownames(presence))
  if (anyNA(tip_order)) {
    warning("reconstruct_gain_loss: species missing from presence matrix: ",
            paste(tr$tip.label[is.na(tip_order)], collapse = ", "),
            " -- treating as absent on those lineages.")
  }
  ok <- !is.na(tip_order)
  aligned[ok, ] <- as.integer(presence[tip_order[ok], , drop = FALSE])

  if (is.null(gain_weights)) gain_weights <- rep(1, n_hog)
  if (length(gain_weights) != n_hog)
    stop("gain_weights must have length ncol(presence).")

  topo <- .gl_topology(tr)

  states_mat <- switch(
    method,
    dollo = .dollo_matrix(tr, aligned, topo),
    fitch = .fitch_matrix(aligned, topo),
    ml    = .ml_matrix(tr, aligned, topo, min_ml_prob, workers)
  )

  # Edge-level gains/losses (weighted)
  parent_states <- states_mat[tr$edge[, 1], , drop = FALSE]
  child_states  <- states_mat[tr$edge[, 2], , drop = FALSE]
  gain_hit <- (parent_states == 0L) & (child_states == 1L)
  loss_hit <- (parent_states == 1L) & (child_states == 0L)
  if (all(gain_weights == 1)) {
    edge_gains  <- rowSums(gain_hit)
    edge_losses <- rowSums(loss_hit)
  } else {
    wmat <- matrix(gain_weights, nrow = nrow(gain_hit), ncol = ncol(gain_hit),
                   byrow = TRUE)
    edge_gains  <- rowSums(gain_hit * wmat)
    edge_losses <- rowSums(loss_hit * wmat)
  }

  child_is_tip <- tr$edge[, 2] <= n_tip
  child_label  <- ifelse(child_is_tip,
                          tr$tip.label[pmin(tr$edge[, 2], n_tip)],
                          NA_character_)
  edges <- tibble::tibble(
    edge_index   = seq_len(nrow(tr$edge)),
    parent       = tr$edge[, 1],
    child        = tr$edge[, 2],
    child_is_tip = child_is_tip,
    child_label  = child_label,
    gains        = edge_gains,
    losses       = edge_losses,
    net_delta    = edge_gains - edge_losses
  )

  out <- list(method = method, edges = edges, tree = tr, n_hogs = n_hog)
  if (return_states) {
    rownames(states_mat) <- c(tr$tip.label, paste0("N", seq_len(tr$Nnode)))
    colnames(states_mat) <- colnames(presence)
    out$states <- states_mat
  }
  class(out) <- c("dnmb_gain_loss", "list")
  out
}

#' Ancestral pan / core size per species-tree node
#'
#' Rolls up a gain/loss reconstruction (with `return_states = TRUE`)
#' into per-node pan-genome and core-genome sizes.
#'
#' @param gl Result of `reconstruct_gain_loss(..., return_states = TRUE)`.
#' @return Tibble with columns `node, is_tip, label, pan_size, core_size`.
#' @export
ancestral_pan_core_sizes <- function(gl) {
  if (is.null(gl$states))
    stop("ancestral_pan_core_sizes: run reconstruct_gain_loss(return_states = TRUE).")
  tr <- gl$tree
  states <- gl$states
  n_tip <- length(tr$tip.label)
  total <- n_tip + tr$Nnode

  pan <- rowSums(states == 1L)

  # tip_under[v, t] = 1 iff tip t is under v
  tip_under <- .tip_under_matrix(tr)
  tip_states <- states[seq_len(n_tip), , drop = FALSE]
  core <- integer(total)
  for (v in seq_len(total)) {
    d <- which(tip_under[v, ])
    if (!length(d)) next
    if (length(d) == 1L) {
      core[v] <- pan[v]
    } else {
      core[v] <- sum(colSums(tip_states[d, , drop = FALSE] == 1L) == length(d))
    }
  }

  tibble::tibble(
    node      = seq_len(total),
    is_tip    = seq_len(total) <= n_tip,
    label     = c(tr$tip.label, paste0("N", seq_len(tr$Nnode))),
    pan_size  = as.integer(pan),
    core_size = core
  )
}

#' Flag HOGs whose presence pattern is HGT-like
#'
#' Runs a single vectorised Dollo pass plus a single vectorised Fitch
#' pass, then compares the implied gain/loss counts for each HOG. A
#' presence pattern that Fitch can explain with fewer events than Dollo
#' — especially when Fitch requires more than one gain — is flagged as
#' an HGT candidate.
#'
#' @param species_tree Rooted `ape::phylo`.
#' @param presence Species × HOG presence matrix.
#' @param min_fitch_gains Flag HOGs whose Fitch reconstruction requires
#'   at least this many gains. Default 2.
#' @param min_loss_diff Alternatively flag HOGs whose Dollo loss count
#'   exceeds the Fitch loss count by this much. Default 2.
#' @return Tibble with one row per HOG and Boolean `hgt_candidate`.
#' @export
detect_hgt_candidates <- function(species_tree, presence,
                                   min_fitch_gains = 2L,
                                   min_loss_diff   = 2L) {
  if (!requireNamespace("ape", quietly = TRUE)) stop("ape required.")
  tr <- species_tree
  n_tip <- length(tr$tip.label)
  n_hog <- ncol(presence)
  hog_ids <- if (!is.null(colnames(presence))) colnames(presence)
             else as.character(seq_len(n_hog))

  # Align presence rows to tip.label
  aligned <- matrix(0L, nrow = n_tip, ncol = n_hog,
                    dimnames = list(tr$tip.label, hog_ids))
  tip_order <- match(tr$tip.label, rownames(presence))
  ok <- !is.na(tip_order)
  aligned[ok, ] <- as.integer(presence[tip_order[ok], , drop = FALSE])

  topo <- .gl_topology(tr)
  sd_mat <- .dollo_matrix(tr, aligned, topo)
  sf_mat <- .fitch_matrix(aligned, topo)
  e <- tr$edge

  dg <- colSums(sd_mat[e[, 1], , drop = FALSE] == 0L & sd_mat[e[, 2], , drop = FALSE] == 1L)
  dl <- colSums(sd_mat[e[, 1], , drop = FALSE] == 1L & sd_mat[e[, 2], , drop = FALSE] == 0L)
  fg <- colSums(sf_mat[e[, 1], , drop = FALSE] == 0L & sf_mat[e[, 2], , drop = FALSE] == 1L)
  fl <- colSums(sf_mat[e[, 1], , drop = FALSE] == 1L & sf_mat[e[, 2], , drop = FALSE] == 0L)

  tibble::tibble(
    hog_id        = hog_ids,
    dollo_gains   = as.integer(dg),
    dollo_losses  = as.integer(dl),
    fitch_gains   = as.integer(fg),
    fitch_losses  = as.integer(fl),
    hgt_candidate = fg >= min_fitch_gains | (dl - fl) >= min_loss_diff
  )
}

# ---- internal helpers ------------------------------------------------------

.gl_topology <- function(tr) {
  n_tip <- length(tr$tip.label)
  total <- n_tip + tr$Nnode
  children <- vector("list", total)
  parent   <- integer(total)
  for (i in seq_len(nrow(tr$edge))) {
    p <- tr$edge[i, 1]; c <- tr$edge[i, 2]
    children[[p]] <- c(children[[p]], c)
    parent[c] <- p
  }
  po <- .gl_post_order(tr, children)
  list(children = children, parent = parent, post_order = po,
       root = n_tip + 1L, n_tip = n_tip, total = total)
}

.gl_post_order <- function(tr, children = NULL) {
  n_tip <- length(tr$tip.label)
  total <- n_tip + tr$Nnode
  if (is.null(children)) {
    children <- vector("list", total)
    for (i in seq_len(nrow(tr$edge))) {
      p <- tr$edge[i, 1]; c <- tr$edge[i, 2]
      children[[p]] <- c(children[[p]], c)
    }
  }
  out <- integer(total); ptr <- 0L
  stack <- list(list(node = n_tip + 1L, state = 0L))
  while (length(stack)) {
    top <- stack[[length(stack)]]
    if (top$state == 0L) {
      stack[[length(stack)]]$state <- 1L
      for (c in children[[top$node]])
        stack[[length(stack) + 1L]] <- list(node = c, state = 0L)
    } else {
      ptr <- ptr + 1L
      out[ptr] <- top$node
      stack[[length(stack)]] <- NULL
    }
  }
  out[seq_len(ptr)]
}

# Compute logical matrix [total × n_tip] where [v, t] is TRUE iff tip t
# is in the subtree rooted at v.
.tip_under_matrix <- function(tr) {
  n_tip <- length(tr$tip.label)
  total <- n_tip + tr$Nnode
  m <- matrix(FALSE, nrow = total, ncol = n_tip)
  for (t in seq_len(n_tip)) m[t, t] <- TRUE
  children <- vector("list", total)
  for (i in seq_len(nrow(tr$edge))) {
    p <- tr$edge[i, 1]; c <- tr$edge[i, 2]
    children[[p]] <- c(children[[p]], c)
  }
  for (v in .gl_post_order(tr, children)) {
    if (v <= n_tip) next
    for (c in children[[v]]) m[v, ] <- m[v, ] | m[c, ]
  }
  m
}

# Vectorised Dollo parsimony over all HOG columns.
# Returns integer matrix [total × n_hog] with ancestral states.
.dollo_matrix <- function(tr, presence, topo) {
  n_tip <- topo$n_tip; total <- topo$total
  n_hog <- ncol(presence)
  tip_under <- .tip_under_matrix(tr)   # total × n_tip

  # has_any_present[v, k] = any descendant tip of v is present in HOG k.
  has_any_present <- (tip_under %*% presence) > 0L
  mode(has_any_present) <- "integer"
  # all_tips_of_v_present[v, k] = every descendant tip of v is present
  # — used by nothing here, but handy.

  # For each HOG, find MRCA and mark in_subtree vector.
  in_subtree <- matrix(0L, nrow = total, ncol = n_hog)
  # Precompute descendant nodes under each node (including node itself)
  # as a logical matrix [total × total]. Memory: for N=200 tips+nodes
  # that's 40k bits — fine. For very large species trees we still get
  # manageable memory (10k² = 100M bits ≈ 12.5 MB).
  node_under <- .node_under_matrix(tr, topo)

  for (k in seq_len(n_hog)) {
    pt <- which(presence[, k] == 1L)
    if (!length(pt)) next
    if (length(pt) == 1L) {
      # Single-tip gain: mrca is the tip itself.
      in_subtree[pt, k] <- 1L
      next
    }
    mrca_k <- ape::getMRCA(tr, tr$tip.label[pt])
    if (is.null(mrca_k)) next
    in_subtree[, k] <- as.integer(node_under[mrca_k, ])
  }
  state <- in_subtree * as.integer(has_any_present)
  storage.mode(state) <- "integer"
  state
}

# For each node n, node_under[n, v] = 1 iff v is in subtree rooted at n
# (including v == n).
.node_under_matrix <- function(tr, topo) {
  total <- topo$total
  m <- diag(1L, total)
  for (v in topo$post_order) {
    if (v <= topo$n_tip) next
    for (c in topo$children[[v]]) m[v, ] <- m[v, ] | m[c, ]
  }
  storage.mode(m) <- "integer"
  m
}

# Vectorised Fitch parsimony over all HOG columns using bitmask encoding
# (1 = {0}, 2 = {1}, 3 = {0,1}).
.fitch_matrix <- function(presence, topo) {
  n_tip <- topo$n_tip; total <- topo$total
  n_hog <- ncol(presence)
  bits <- matrix(0L, nrow = total, ncol = n_hog)
  # Tips
  bits[seq_len(n_tip), ] <- ifelse(presence == 1L, 2L, 1L)

  # Up-pass
  for (v in topo$post_order) {
    if (v <= n_tip) next
    ch <- topo$children[[v]]
    if (!length(ch)) next
    if (length(ch) == 1L) {
      bits[v, ] <- bits[ch, ]
      next
    }
    inter <- bits[ch[1], ]
    uni   <- bits[ch[1], ]
    for (j in seq.int(2L, length(ch))) {
      inter <- bitwAnd(inter, bits[ch[j], ])
      uni   <- bitwOr(uni,   bits[ch[j], ])
    }
    bits[v, ] <- ifelse(inter != 0L, inter, uni)
  }

  # Down-pass: pick concrete 0/1 states.
  state <- matrix(0L, nrow = total, ncol = n_hog)
  root <- topo$root
  state[root, ] <- ifelse(bitwAnd(bits[root, ], 1L) != 0L, 0L, 1L)
  pre_order <- rev(topo$post_order)
  for (v in pre_order) {
    if (v == root) next
    p <- topo$parent[v]
    parent_mask <- ifelse(state[p, ] == 0L, 1L, 2L)
    keep_parent <- bitwAnd(bits[v, ], parent_mask) != 0L
    alt_state   <- ifelse(bitwAnd(bits[v, ], 1L) != 0L, 0L, 1L)
    state[v, ] <- ifelse(keep_parent, state[p, ], alt_state)
  }
  storage.mode(state) <- "integer"
  state
}

# ML via ape::ace, with column de-duplication and optional parallelism.
.ml_matrix <- function(tr, presence, topo, min_prob, workers) {
  n_tip <- topo$n_tip; total <- topo$total
  n_hog <- ncol(presence)

  # Dedupe: hash each column to find unique patterns.
  keys <- apply(presence, 2L, paste0, collapse = "")
  uniq <- unique(keys)
  key_to_idx <- stats::setNames(seq_along(uniq), uniq)
  rep_cols <- match(uniq, keys)

  do_one <- function(col_idx) {
    p <- presence[, col_idx]
    st <- integer(total)
    st[seq_len(n_tip)] <- as.integer(p)
    if (length(unique(p)) < 2L) {
      st[] <- as.integer(p[1])
      return(st)
    }
    x <- factor(p, levels = c(0L, 1L))
    res <- tryCatch(
      ape::ace(x, tr, type = "discrete", model = "ARD"),
      error = function(e) NULL
    )
    if (is.null(res)) {
      # fall back to Fitch for this column
      sub_topo <- topo
      sub <- .fitch_matrix(matrix(p, ncol = 1L), sub_topo)
      return(sub[, 1])
    }
    la <- res$lik.anc
    col1 <- if ("1" %in% colnames(la)) "1" else ncol(la)
    for (i in seq_len(nrow(la))) {
      v <- n_tip + i
      p1 <- la[i, col1]
      st[v] <- if (!is.na(p1) && p1 >= min_prob) 1L else 0L
    }
    st
  }

  use_parallel <- workers > 1L && requireNamespace("future.apply", quietly = TRUE)
  if (use_parallel) {
    if (!requireNamespace("future", quietly = TRUE)) use_parallel <- FALSE
  }
  if (use_parallel) {
    prev_plan <- future::plan()
    if (inherits(prev_plan, "sequential")) {
      future::plan(future::multisession, workers = workers)
      on.exit(future::plan(prev_plan), add = TRUE)
    }
    uniq_states <- future.apply::future_lapply(
      rep_cols, do_one,
      future.seed = TRUE,
      future.packages = c("ape")
    )
  } else {
    uniq_states <- lapply(rep_cols, do_one)
  }

  state <- matrix(0L, nrow = total, ncol = n_hog)
  idx_for_col <- key_to_idx[keys]
  for (k in seq_len(n_hog)) {
    state[, k] <- uniq_states[[idx_for_col[k]]]
  }
  state
}
