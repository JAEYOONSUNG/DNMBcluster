#' Reconcile a gene tree against a rooted species tree (LCA + heuristic DTL)
#'
#' Implements a pragmatic DTL-style reconciliation:
#'
#' \enumerate{
#'   \item Map every gene-tree node to the LCA of its descendant species
#'         in the species tree ("LCA mapping").
#'   \item Label every internal gene-tree node as \bold{S} (speciation)
#'         or \bold{D} (duplication) using the species-overlap
#'         heuristic: overlap → D, disjoint → S.
#'   \item Flag duplications as transfer-candidates (\bold{T}) when one
#'         child's species set is NOT a descendant of the other child's
#'         map in the species tree, which is Ranger-DTL's classical
#'         heuristic signature of a horizontal transfer.
#'   \item For every speciation, count losses on the species-tree edges
#'         that a child subtree fails to descend into.
#' }
#'
#' This is a deliberately lightweight, pure-R implementation — it does
#' NOT search the cost space of every possible labeling the way
#' Ranger-DTL / ALE do. For rigorous DTL inference, export gene and
#' species trees via `export_og_bundles()` and run GeneRax/ALE
#' externally. This routine is fast enough to run on every OG during
#' `run_orthofinder_like()` and produces a per-node events table that
#' is useful for downstream summaries, HGT enrichment, and plotting.
#'
#' @param gene_tree Rooted `ape::phylo`.
#' @param species_tree Rooted `ape::phylo`.
#' @param tip_to_species Named character vector mapping
#'   `gene_tree$tip.label` → `species_tree$tip.label`. If NULL, the
#'   function attempts to parse the `p<uid>_g<key>` label scheme used
#'   elsewhere in DNMBcluster.
#' @param count_losses If TRUE, also return a per-species-edge loss
#'   count tibble.
#' @param min_dup_support When the gene tree carries `node.label`
#'   bootstrap percentages, only classify a node as D/T if support is
#'   at least this value. Defaults to 0 (no gating).
#' @param long_branch_quantile Relative branch-length gate for T
#'   calls: when both child edges of a candidate T exist, classify
#'   only if at least one child edge length exceeds this quantile of
#'   all gene-tree edge lengths (default 0.9). Set to 1 to disable
#'   the gate. Has no effect on trees without edge lengths.
#' @return A list with fields:
#'   - `events`: tibble with one row per gene-tree internal node
#'     (`g_node, g_n_leaves, sp_lca, event, species_child_map`).
#'   - `losses`: tibble (if `count_losses`) with per-species-edge
#'     losses. Columns `sp_parent, sp_child, loss_count`.
#'   - `summary`: one-row tibble (`n_S, n_D, n_T, n_loss`).
#' @export
reconcile_dtl <- function(gene_tree,
                           species_tree,
                           tip_to_species = NULL,
                           count_losses = TRUE,
                           min_dup_support = 0,
                           long_branch_quantile = 0.9) {
  if (!requireNamespace("ape", quietly = TRUE)) stop("ape required.")
  if (!ape::is.rooted(species_tree)) stop("species_tree must be rooted.")
  if (!ape::is.rooted(gene_tree))    stop("gene_tree must be rooted.")

  if (is.null(tip_to_species)) {
    tip_to_species <- .dtl_infer_tip_map(gene_tree, species_tree)
  }
  missing <- setdiff(gene_tree$tip.label, names(tip_to_species))
  if (length(missing)) {
    stop("reconcile_dtl: tip_to_species missing gene-tree tips: ",
         paste(utils::head(missing, 5), collapse = ", "),
         if (length(missing) > 5) ", ..." else "")
  }
  sp_keys <- unname(tip_to_species[gene_tree$tip.label])
  keep <- sp_keys %in% species_tree$tip.label
  if (!any(keep)) stop("reconcile_dtl: no gene-tree tips map to a species-tree tip.")
  if (!all(keep)) {
    gene_tree <- ape::keep.tip(gene_tree, gene_tree$tip.label[keep])
    sp_keys   <- unname(tip_to_species[gene_tree$tip.label])
  }

  n_gtip <- length(gene_tree$tip.label)
  n_gint <- gene_tree$Nnode
  total_g <- n_gtip + n_gint

  g_children <- vector("list", total_g)
  for (i in seq_len(nrow(gene_tree$edge))) {
    p <- gene_tree$edge[i, 1]; c <- gene_tree$edge[i, 2]
    g_children[[p]] <- c(g_children[[p]], c)
  }
  g_post <- .dtl_post_order(gene_tree, g_children)

  # Species set per gene-tree node
  sp_set <- vector("list", total_g)
  for (t in seq_len(n_gtip)) sp_set[[t]] <- sp_keys[t]
  for (v in g_post) {
    if (v <= n_gtip) next
    sp_set[[v]] <- unique(unlist(sp_set[g_children[[v]]]))
  }

  # LCA map: each gene node → species-tree node
  n_stip <- length(species_tree$tip.label)
  lca_map <- integer(total_g)
  for (v in seq_len(total_g)) {
    ss <- sp_set[[v]]
    ss <- ss[ss %in% species_tree$tip.label]
    if (!length(ss)) { lca_map[v] <- NA_integer_; next }
    if (length(ss) == 1L) lca_map[v] <- which(species_tree$tip.label == ss)[1]
    else lca_map[v] <- ape::getMRCA(species_tree, ss)
  }

  # Species-tree descendants-of-each-node lookup (for loss + T check)
  s_children <- vector("list", n_stip + species_tree$Nnode)
  for (i in seq_len(nrow(species_tree$edge))) {
    p <- species_tree$edge[i, 1]; c <- species_tree$edge[i, 2]
    s_children[[p]] <- c(s_children[[p]], c)
  }
  desc_of <- .dtl_desc_sets(species_tree, s_children)

  # Event labels per internal gene node
  event <- rep(NA_character_, total_g)
  loss_acc <- integer(nrow(species_tree$edge))

  node_support <- .dtl_node_support(gene_tree, n_gtip)
  edge_len_by_child <- if (!is.null(gene_tree$edge.length) &&
                            length(gene_tree$edge.length) == nrow(gene_tree$edge)) {
    stats::setNames(gene_tree$edge.length, gene_tree$edge[, 2])
  } else {
    stats::setNames(numeric(0), character(0))
  }
  long_thr <- if (length(edge_len_by_child) && long_branch_quantile < 1) {
    stats::quantile(edge_len_by_child, long_branch_quantile, na.rm = TRUE)
  } else NA_real_

  for (v in g_post) {
    if (v <= n_gtip) next
    ch <- g_children[[v]]
    if (length(ch) < 2L) { event[v] <- "S"; next }
    s1 <- sp_set[[ch[1]]]; s2 <- sp_set[[ch[2]]]
    shared <- intersect(s1, s2)
    if (length(shared)) {
      if (min_dup_support > 0 && !is.na(node_support[v]) &&
          node_support[v] < min_dup_support) {
        event[v] <- "S"
        next
      }
      m1 <- lca_map[ch[1]]; m2 <- lca_map[ch[2]]
      if (!is.na(m1) && !is.na(m2)) {
        if (!(m1 %in% desc_of[[m2]] || m2 %in% desc_of[[m1]])) {
          # T candidate: require a long-branch sibling to reduce
          # false positives from poorly resolved gene trees.
          if (is.na(long_thr)) {
            event[v] <- "T"
          } else {
            e1 <- edge_len_by_child[as.character(ch[1])]
            e2 <- edge_len_by_child[as.character(ch[2])]
            any_long <- any(c(e1, e2) >= long_thr, na.rm = TRUE)
            event[v] <- if (any_long) "T" else "D"
          }
        } else {
          event[v] <- "D"
        }
      } else {
        event[v] <- "D"
      }
    } else {
      event[v] <- "S"
      if (count_losses && !is.na(lca_map[v])) {
        m_v <- lca_map[v]
        if (m_v > n_stip) {
          # Losses: for each species-tree child of m_v, if no gene-tree
          # child covers its descendants, count a loss.
          sc <- s_children[[m_v]]
          for (sc_i in sc) {
            sc_desc_tips <- desc_of[[sc_i]]
            sc_tip_set <- species_tree$tip.label[sc_desc_tips[sc_desc_tips <= n_stip]]
            covered <- any(vapply(ch, function(g_c) {
                              any(sp_set[[g_c]] %in% sc_tip_set)
                            }, logical(1)))
            if (!covered) {
              e_idx <- which(species_tree$edge[, 1] == m_v &
                              species_tree$edge[, 2] == sc_i)
              if (length(e_idx)) loss_acc[e_idx] <- loss_acc[e_idx] + 1L
            }
          }
        }
      }
    }
  }

  internal_nodes <- which(!is.na(event))
  events_tbl <- tibble::tibble(
    g_node     = internal_nodes,
    g_n_leaves = vapply(internal_nodes,
                         function(v) length(sp_set[[v]]), integer(1)),
    sp_lca     = lca_map[internal_nodes],
    event      = event[internal_nodes]
  )

  summary_tbl <- tibble::tibble(
    n_S = sum(event == "S", na.rm = TRUE),
    n_D = sum(event == "D", na.rm = TRUE),
    n_T = sum(event == "T", na.rm = TRUE),
    n_loss = sum(loss_acc)
  )

  out <- list(events = events_tbl, summary = summary_tbl,
              tree = gene_tree, species_tree = species_tree)
  if (count_losses) {
    out$losses <- tibble::tibble(
      sp_parent  = species_tree$edge[, 1],
      sp_child   = species_tree$edge[, 2],
      loss_count = loss_acc
    )
  }
  class(out) <- c("dnmb_dtl", "list")
  out
}

# ---- helpers ---------------------------------------------------------------

.dtl_post_order <- function(tr, children) {
  n_tip <- length(tr$tip.label)
  out <- integer(n_tip + tr$Nnode); ptr <- 0L
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

.dtl_desc_sets <- function(tr, children) {
  n_tip <- length(tr$tip.label)
  total <- n_tip + tr$Nnode
  out <- vector("list", total)
  po <- .dtl_post_order(tr, children)
  for (v in po) {
    if (v <= n_tip) { out[[v]] <- v; next }
    out[[v]] <- c(v, unlist(lapply(children[[v]], function(c) out[[c]])))
  }
  out
}

# Local null-coalesce (avoids taking hard dep on rlang).
`%||%` <- function(a, b) if (is.null(a)) b else a

# Parse integer support from gene_tree$node.label; NA when absent.
.dtl_node_support <- function(tr, n_tip) {
  total <- n_tip + tr$Nnode
  out <- rep(NA_real_, total)
  if (is.null(tr$node.label) || !length(tr$node.label)) return(out)
  sup <- suppressWarnings(as.numeric(tr$node.label))
  out[(n_tip + 1L):total] <- sup
  out
}

# Try to derive tip_to_species from DNMB's p<uid>_g<key> scheme.
.dtl_infer_tip_map <- function(gene_tree, species_tree) {
  keys <- sub(".*_g", "", gene_tree$tip.label)
  if (!all(nzchar(keys))) stop("reconcile_dtl: could not parse species key from tips.")
  stats::setNames(keys, gene_tree$tip.label)
}


#' Aggregate DTL event counts onto the species tree and plot
#'
#' Takes a per-OG events table produced by `.run_dtl_batch()` (or a
#' list of individual `reconcile_dtl()` results) and summarizes per
#' species-tree node: total speciations, duplications, transfer
#' candidates, and losses (via the paired loss tables). Draws a
#' species tree with node pies (S/D/T) and edge widths proportional
#' to loss counts.
#'
#' @param events_tbl Tibble with columns `cluster_id, g_node, sp_lca,
#'   event, ...` as produced by `.run_dtl_batch()`. The `sp_lca`
#'   column indexes nodes/tips of `species_tree`.
#' @param species_tree The rooted `ape::phylo` used for reconciliation.
#' @param losses_tbl Optional tibble of per-species-edge losses
#'   (`sp_parent, sp_child, loss_count`). If supplied, edges are
#'   drawn with widths proportional to cumulative loss count.
#' @param node_scale Scalar applied to node-pie radius (user units).
#'   Default 0.6.
#' @return Invisibly returns a tibble keyed by `sp_node` with columns
#'   `n_S, n_D, n_T, n_loss_in` (losses on the edge leading to the node).
#' @export
plot_dtl_events <- function(events_tbl,
                            species_tree,
                            losses_tbl = NULL,
                            node_scale = 0.6) {
  if (!requireNamespace("ape", quietly = TRUE)) stop("ape required.")
  stopifnot(inherits(species_tree, "phylo"))

  n_tip <- length(species_tree$tip.label)
  n_int <- species_tree$Nnode
  total <- n_tip + n_int

  counts <- matrix(0L, nrow = total, ncol = 3,
                   dimnames = list(NULL, c("S", "D", "T")))
  if (!is.null(events_tbl) && nrow(events_tbl)) {
    ev <- events_tbl[!is.na(events_tbl$sp_lca) &
                       events_tbl$event %in% c("S", "D", "T"), , drop = FALSE]
    for (i in seq_len(nrow(ev))) {
      counts[ev$sp_lca[i], ev$event[i]] <-
        counts[ev$sp_lca[i], ev$event[i]] + 1L
    }
  }

  loss_in <- integer(total)
  if (!is.null(losses_tbl) && nrow(losses_tbl)) {
    for (i in seq_len(nrow(losses_tbl))) {
      loss_in[losses_tbl$sp_child[i]] <-
        loss_in[losses_tbl$sp_child[i]] + losses_tbl$loss_count[i]
    }
  }

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(2, 1, 2, 1))

  edge_widths <- rep(1, nrow(species_tree$edge))
  if (any(loss_in > 0L)) {
    child_idx <- species_tree$edge[, 2]
    w <- loss_in[child_idx]
    if (max(w) > 0) edge_widths <- 1 + 3 * (w / max(w))
  }

  ape::plot.phylo(species_tree, edge.width = edge_widths,
                  label.offset = 0.02, show.tip.label = TRUE,
                  no.margin = FALSE)
  xy <- ape::nodelabels(text = NULL, frame = "n", adj = c(0.5, 0.5),
                        pch = NA, node = seq_len(total))
  pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)

  for (v in seq_len(total)) {
    tot <- sum(counts[v, ])
    if (tot == 0L) next
    r <- node_scale * sqrt(tot / max(1L, max(rowSums(counts))))
    x <- pp$xx[v]; y <- pp$yy[v]
    .dtl_draw_pie(x, y, r, counts[v, ], cols = c(S = "#2E7D32",
                                                  D = "#E65100",
                                                  T = "#C62828"))
    graphics::text(x, y, labels = tot, cex = 0.55, col = "white")
  }

  graphics::legend("topleft",
                   legend = c("Speciation", "Duplication", "Transfer"),
                   fill = c("#2E7D32", "#E65100", "#C62828"),
                   bty = "n", cex = 0.8)
  if (any(loss_in > 0L)) {
    graphics::mtext(sprintf("edge width ~ loss count (max = %d)", max(loss_in)),
                    side = 1, line = 0.5, cex = 0.7)
  }

  invisible(tibble::tibble(
    sp_node   = seq_len(total),
    n_S       = counts[, "S"],
    n_D       = counts[, "D"],
    n_T       = counts[, "T"],
    n_loss_in = loss_in
  ))
}

.dtl_draw_pie <- function(x, y, r, counts, cols) {
  total <- sum(counts)
  if (!total) return(invisible(NULL))
  theta0 <- 0
  for (nm in names(counts)) {
    k <- counts[[nm]]
    if (!k) next
    theta1 <- theta0 + 2 * pi * (k / total)
    ang <- seq(theta0, theta1, length.out = max(8L, ceiling((theta1 - theta0) * 20)))
    graphics::polygon(x + c(0, r * cos(ang)),
                      y + c(0, r * sin(ang)),
                      col = cols[[nm]], border = NA)
    theta0 <- theta1
  }
  graphics::symbols(x, y, circles = r, inches = FALSE,
                    add = TRUE, fg = "grey30")
}
