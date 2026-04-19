#' Branch-annotated DTL event overlay (ALE / GeneRax convention)
#'
#' Redraws the DTL reconciliation summary in the form used throughout
#' the ALE / GeneRax literature: events are placed on **branches**, not
#' at nodes, because duplications, transfers and losses occur *along a
#' lineage* between divergences. Published figures use a stacked bar at
#' the branch midpoint with a fixed S / D / T / L color key for
#' aggregated counts across many OGs — this function follows that
#' convention and pairs it with a companion panel showing per-branch
#' event totals.
#'
#' Event-to-edge mapping:
#'
#' \itemize{
#'   \item **D**, **T** — placed on the edge *leading into* the
#'         species node recorded as their LCA (`sp_lca == child`).
#'   \item **loss** — already edge-keyed in `dtl_losses.tsv`
#'         (`sp_parent -> sp_child`).
#'   \item **S** — shown for reference only; speciations track the tree
#'         topology itself and add little information at this zoom
#'         level. Rendered as a lighter tint behind D / T / L.
#' }
#'
#' Rendering uses `ggtree` for the phylogram (when available) combined
#' with a custom branch-midpoint `geom_rect` layer via
#' `ggnewscale::new_scale_fill()` so the tree and stacked-bar fill
#' scales do not collide. `ggtreeExtra::geom_fruit` is *not* used —
#' its layout is for tip-aligned rings, not branch annotations. For
#' sites where `ggtree` is unavailable, a manual ggplot tree (same
#' layout used elsewhere in the package) is used as a fallback.
#'
#' @param events_tbl `data.frame` with columns `cluster_id`, `sp_lca`
#'   and `event` (one of `S`/`D`/`T`).
#' @param species_tree A rooted `ape::phylo`.
#' @param losses_tbl Optional per-OG loss table with columns
#'   `sp_parent`, `sp_child`, `loss_count`.
#' @param out_pdf Output path. If `NULL`, returns the patchwork object.
#' @param width,height PDF dimensions (inches).
#' @param bar_scale Scalar multiplier for stacked-bar heights.
#' @return Invisibly: the written path, or the patchwork object if
#'   `out_pdf` is NULL.
#' @export
plot_dtl_branch_events <- function(events_tbl,
                                    species_tree,
                                    losses_tbl = NULL,
                                    out_pdf = NULL,
                                    width = 15,
                                    height = 8.3,
                                    bar_scale = 0.9) {
  stopifnot(inherits(species_tree, "phylo"))
  for (pkg in c("ggplot2", "ape")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("[plot_dtl_branch_events] %s required", pkg))
  }
  tree <- ape::ladderize(species_tree)

  layout  <- .dtl_tree_layout(tree)
  edges   <- .dtl_edge_counts(tree, layout, events_tbl, losses_tbl)
  bars    <- .dtl_stacked_bar_df(edges, bar_scale)

  pal_loc <- .dnmb_anvio_pal
  ev_cols <- c(D = pal_loc$red, T = pal_loc$navy, L = pal_loc$text)
  n_ogs <- length(unique(events_tbl$cluster_id[
    events_tbl$event %in% c("S", "D", "T")]))

  attr(edges, "tree") <- tree
  tree_panel <- .dtl_tree_panel(tree, layout, edges, bars, ev_cols)
  side_panel <- .dtl_side_bar_panel(edges, ev_cols)

  design <- patchwork::wrap_plots(tree_panel, side_panel,
                                    ncol = 2, widths = c(2.85, 1.2)) +
    patchwork::plot_layout(guides = "collect", axis_titles = "collect") +
    patchwork::plot_annotation(
      title = "DTL reconciliation events per species-tree branch",
      subtitle = sprintf(
        paste0("%s OGs | totals S=%s D=%s T=%s L=%s | ",
               "left: log10(count+1) stacks on each branch midpoint | ",
               "right: raw D / T / L counts per branch"),
        formatC(n_ogs, format = "d", big.mark = ","),
        formatC(sum(edges$n_S), format = "d", big.mark = ","),
        formatC(sum(edges$n_D), format = "d", big.mark = ","),
        formatC(sum(edges$n_T), format = "d", big.mark = ","),
        formatC(sum(edges$n_L), format = "d", big.mark = ",")
      ),
      caption = "events placed on the edge leading into their sp_lca (ALE / GeneRax convention)",
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", size = 15,
                                                colour = "#1A2238"),
        plot.subtitle = ggplot2::element_text(size = 10, colour = "grey35"),
        plot.caption  = ggplot2::element_text(size = 8,  colour = "grey45",
                                                hjust = 1)
      )
    ) &
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = "#F4F1EA", colour = NA),
      panel.background = ggplot2::element_rect(fill = "#F4F1EA", colour = NA),
      legend.position  = "right",
      legend.box       = "vertical",
      legend.direction = "vertical",
      legend.justification = "top",
      legend.box.just  = "left",
      legend.spacing.y = ggplot2::unit(0.55, "cm"),
      legend.margin    = ggplot2::margin(2, 2, 2, 2),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.key.size  = ggplot2::unit(0.32, "cm"),
      legend.text      = ggplot2::element_text(size = 7.5, colour = "#1A2238",
                                                 margin = ggplot2::margin(l = 2)),
      legend.title     = ggplot2::element_text(size = 8.5, face = "bold",
                                                 colour = "#1A2238",
                                                 lineheight = 1.0,
                                                 margin = ggplot2::margin(b = 4))
    )

  if (is.null(out_pdf)) return(invisible(design))
  ggplot2::ggsave(out_pdf, design, width = width, height = height,
                   dpi = 300, limitsize = FALSE, bg = "#F4F1EA")
  invisible(out_pdf)
}


# ---- tree coordinates --------------------------------------------------

.dtl_tree_layout <- function(tree) {
  n_tip  <- length(tree$tip.label)
  n_node <- tree$Nnode
  total  <- n_tip + n_node

  xpos <- numeric(total); ypos <- numeric(total)
  ypos[seq_len(n_tip)] <- seq_len(n_tip)

  root <- n_tip + 1L
  visited <- logical(total)
  queue <- c(root); visited[root] <- TRUE; xpos[root] <- 0
  while (length(queue)) {
    parent <- queue[1]; queue <- queue[-1]
    rows <- which(tree$edge[, 1] == parent)
    for (ri in rows) {
      child <- tree$edge[ri, 2]
      blen <- tree$edge.length[ri]
      if (is.null(blen) || is.na(blen)) blen <- 1
      xpos[child] <- xpos[parent] + blen
      if (!visited[child]) { queue <- c(queue, child); visited[child] <- TRUE }
    }
  }
  for (k in rev(seq_len(n_node))) {
    nid <- n_tip + k
    desc <- which(tree$edge[, 1] == nid)
    ypos[nid] <- mean(ypos[tree$edge[desc, 2]])
  }

  list(n_tip = n_tip, n_node = n_node, total = total,
       xpos = xpos, ypos = ypos, max_x = max(xpos))
}

.dtl_edge_counts <- function(tree, layout, events_tbl, losses_tbl) {
  edge_df <- data.frame(
    parent = tree$edge[, 1],
    child  = tree$edge[, 2],
    x_p    = layout$xpos[tree$edge[, 1]],
    x_c    = layout$xpos[tree$edge[, 2]],
    y      = layout$ypos[tree$edge[, 2]]
  )
  edge_df$x_mid <- (edge_df$x_p + edge_df$x_c) / 2
  edge_df$label <- ifelse(edge_df$child <= layout$n_tip,
                           tree$tip.label[edge_df$child],
                           sprintf("N%03d",
                                    edge_df$child - layout$n_tip))

  counts <- matrix(0L, nrow = layout$total, ncol = 4L,
                    dimnames = list(NULL, c("S", "D", "T", "L")))
  if (!is.null(events_tbl) && nrow(events_tbl)) {
    ev <- events_tbl[!is.na(events_tbl$sp_lca) &
                       events_tbl$event %in% c("S", "D", "T"), , drop = FALSE]
    if (nrow(ev)) {
      tab <- table(ev$sp_lca, ev$event)
      idx <- as.integer(rownames(tab))
      for (ev_type in intersect(colnames(tab), c("S", "D", "T")))
        counts[idx, ev_type] <- counts[idx, ev_type] + tab[, ev_type]
    }
  }
  if (!is.null(losses_tbl) && nrow(losses_tbl)) {
    la <- stats::aggregate(loss_count ~ sp_child, data = losses_tbl,
                            FUN = sum)
    counts[la$sp_child, "L"] <- counts[la$sp_child, "L"] + la$loss_count
  }

  edge_df$n_S <- counts[edge_df$child, "S"]
  edge_df$n_D <- counts[edge_df$child, "D"]
  edge_df$n_T <- counts[edge_df$child, "T"]
  edge_df$n_L <- counts[edge_df$child, "L"]
  edge_df$n_total <- edge_df$n_S + edge_df$n_D +
                      edge_df$n_T + edge_df$n_L
  edge_df
}

.dtl_stacked_bar_df <- function(edge_df, bar_scale) {
  # S is now encoded on the branch itself (edge colour); stacked bars
  # only carry the "reconciliation" events D / T / L so the bars show
  # actionable signal rather than being dominated by topology-S.
  rows <- do.call(rbind, lapply(seq_len(nrow(edge_df)), function(i) {
    counts_i <- c(D = edge_df$n_D[i],
                   T = edge_df$n_T[i], L = edge_df$n_L[i])
    if (sum(counts_i) == 0) return(NULL)
    order_keys <- c("L", "T", "D")
    h <- log10(counts_i[order_keys] + 1)
    cum <- cumsum(h)
    data.frame(
      edge_i = i,
      event  = factor(order_keys, levels = c("D", "T", "L")),
      x_mid  = edge_df$x_mid[i],
      y      = edge_df$y[i],
      y_low  = c(0, utils::head(cum, -1L)),
      y_high = cum,
      raw    = as.integer(counts_i[order_keys])
    )
  }))
  if (is.null(rows) || !nrow(rows)) return(NULL)
  max_stack <- max(rows$y_high, 1)
  rows$y_scale <- (0.80 * bar_scale) / max_stack
  rows$y0 <- rows$y + rows$y_low  * rows$y_scale
  rows$y1 <- rows$y + rows$y_high * rows$y_scale
  rows
}


# ---- panel A: tree with branch-midpoint stacks ------------------------

.dtl_tree_panel <- function(tree, layout, edge_df, bar_rows, ev_cols) {
  use_ggtree <- requireNamespace("ggtree", quietly = TRUE) &&
                requireNamespace("ggnewscale", quietly = TRUE)
  if (use_ggtree) {
    return(.dtl_tree_panel_ggtree(tree, layout, edge_df, bar_rows, ev_cols))
  }
  .dtl_tree_panel_manual(tree, layout, edge_df, bar_rows, ev_cols)
}

.dtl_tree_panel_ggtree <- function(tree, layout, edge_df, bar_rows, ev_cols) {
  n_tip <- layout$n_tip; n_node <- layout$n_node
  max_x <- layout$max_x

  # ggtree is used purely for its ladderized coordinates (size = 0 hides
  # its default edges so we can draw our own event-aware branches).
  base <- ggtree::ggtree(tree, ladderize = TRUE, size = 0, colour = NA)
  td <- base$data
  node_map <- data.frame(node = td$node, x = td$x, y = td$y)

  edge_p <- merge(edge_df[, c("parent", "child", "label",
                                "n_S", "n_D", "n_T", "n_L", "n_total")],
                   node_map, by.x = "parent", by.y = "node",
                   suffixes = c("", ".p"))
  names(edge_p)[names(edge_p) == "x"] <- "x_p"
  names(edge_p)[names(edge_p) == "y"] <- "y_p"
  edge_p <- merge(edge_p, node_map, by.x = "child", by.y = "node")
  names(edge_p)[names(edge_p) == "x"] <- "x_c"
  names(edge_p)[names(edge_p) == "y"] <- "y_c"
  edge_p$x_mid <- (edge_p$x_p + edge_p$x_c) / 2

  # Biologically meaningful edge colour: duplication density (D per S).
  # A branch where gene families duplicate often is distinct from a
  # large-S backbone with no duplication — raw S alone does not show
  # that.
  edge_p$dup_rate <- edge_p$n_D / pmax(edge_p$n_S, 1L)

  # Rank branches by raw D count for top-K highlighting.
  dup_order <- order(-edge_p$n_D)
  top_k <- min(3L, sum(edge_p$n_D > 0))
  edge_p$is_top <- FALSE
  if (top_k > 0) edge_p$is_top[dup_order[seq_len(top_k)]] <- TRUE

  # Realign bar-row coordinates to ggtree's rendered y ordering.
  bar_new <- NULL
  if (!is.null(bar_rows) && nrow(bar_rows)) {
    bar_new <- bar_rows
    bar_new$x_mid <- edge_p$x_mid[match(bar_rows$edge_i,
                                           seq_len(nrow(edge_df)))]
    bar_new$y     <- edge_p$y_c[match(bar_rows$edge_i,
                                         seq_len(nrow(edge_df)))]
    bar_new$y0 <- bar_new$y + bar_new$y_low  * bar_new$y_scale
    bar_new$y1 <- bar_new$y + bar_new$y_high * bar_new$y_scale
  }

  parents_u <- unique(tree$edge[, 1])
  seg_v <- do.call(rbind, lapply(parents_u, function(pnode) {
    children <- tree$edge[tree$edge[, 1] == pnode, 2]
    x  <- node_map$x[match(pnode, node_map$node)]
    ys <- node_map$y[match(children, node_map$node)]
    data.frame(x = x, xend = x, y = min(ys), yend = max(ys))
  }))

  tip_map <- node_map[node_map$node <= n_tip, , drop = FALSE]
  tip_map$label <- tree$tip.label[tip_map$node]
  tip_map$x_end <- max_x * 1.02
  node_lab <- node_map[node_map$node > n_tip, , drop = FALSE]
  node_lab$label <- sprintf("N%03d", node_lab$node - n_tip)

  # Per-tip lineage rollup: sum D/T/L along the root->tip path. Shown
  # as a miniature stacked bar anchored past the tip label for easy
  # per-species reading ("which lineages saw the most dup/loss total").
  lineage_rollup <- do.call(rbind, lapply(tip_map$node, function(tip) {
    path <- numeric(0)
    cur  <- tip
    while (length(cur) && length(which(tree$edge[, 2] == cur))) {
      path <- c(path, cur)
      cur  <- tree$edge[which(tree$edge[, 2] == cur), 1]
    }
    ei <- which(edge_df$child %in% path)
    data.frame(
      node = tip,
      y    = tip_map$y[match(tip, tip_map$node)],
      n_D  = sum(edge_df$n_D[ei]),
      n_T  = sum(edge_df$n_T[ei]),
      n_L  = sum(edge_df$n_L[ei])
    )
  }))
  max_roll <- max(rowSums(lineage_rollup[, c("n_D", "n_T", "n_L")]), 1)

  # Clade shading: colour the two subtrees below the root with alternating
  # light bands so the viewer sees the main dichotomy at a glance.
  root <- n_tip + 1L
  root_children <- tree$edge[tree$edge[, 1] == root, 2]
  clade_bands <- do.call(rbind, lapply(seq_along(root_children),
                                         function(i) {
    rc <- root_children[i]
    leaves <- if (rc <= n_tip) rc else
      ape::extract.clade(tree, rc)$tip.label
    if (is.character(leaves))
      leaves <- tip_map$node[match(leaves, tip_map$label)]
    ys <- tip_map$y[match(leaves, tip_map$node)]
    data.frame(
      ymin = min(ys) - 0.46,
      ymax = max(ys) + 0.46,
      fill = if (i %% 2 == 1) "#F3F0EA" else "#EAE5D6"
    )
  }))

  # Scale bar geometry.
  scale_len <- signif(max_x / 5, 1)
  scale_x0  <- max_x * 0.02
  scale_x1  <- scale_x0 + scale_len
  scale_y   <- -0.15

  bar_halfw <- max_x * 0.010
  max_S <- max(edge_p$n_S, 1)

  # Summary label under each branch: raw D/T/L numbers (skip zeros).
  # Only emit labels for the top-K busiest branches so dense species trees
  # don't turn into a wall of overlapping text. Everything else is still
  # readable via the D/T/L side heatmap.
  totals_all <- edge_p[edge_p$n_D + edge_p$n_T + edge_p$n_L > 0 &
                          !(edge_p$is_top &
                             edge_p$n_D == max(edge_p$n_D, na.rm = TRUE)),
                        , drop = FALSE]
  top_label_k <- min(8L, nrow(totals_all))
  if (top_label_k > 0) {
    totals_df <- totals_all[order(-(totals_all$n_D + totals_all$n_T +
                                      totals_all$n_L)), ][seq_len(top_label_k), ,
                                                            drop = FALSE]
    totals_df$label_full <- vapply(seq_len(nrow(totals_df)), function(i) {
      parts <- c(
        if (totals_df$n_D[i] > 0)
          sprintf("D %s",
                   formatC(totals_df$n_D[i], format = "d", big.mark = ",")),
        if (totals_df$n_T[i] > 0)
          sprintf("T %s",
                   formatC(totals_df$n_T[i], format = "d", big.mark = ",")),
        if (totals_df$n_L[i] > 0)
          sprintf("L %s",
                   formatC(totals_df$n_L[i], format = "d", big.mark = ","))
      )
      paste(parts, collapse = " | ")
    }, character(1))
  } else {
    totals_df <- totals_all[0, , drop = FALSE]
    totals_df$label_full <- character(0)
  }

  # ---- assemble ---------------------------------------------------------
  p <- base +
    # Clade bands — alternating light rectangles for the two root subtrees
    # so the viewer sees the main dichotomy at a glance.
    ggplot2::geom_rect(
      data = clade_bands,
      ggplot2::aes(xmin = -Inf, xmax = Inf,
                    ymin = .data$ymin, ymax = .data$ymax,
                    fill = .data$fill),
      inherit.aes = FALSE, alpha = 0.55, colour = NA) +
    ggplot2::scale_fill_identity() +
    # Dotted tip-extensions so italic names all left-align.
    ggplot2::geom_segment(
      data = tip_map,
      ggplot2::aes(x = .data$x, xend = .data$x_end,
                    y = .data$y, yend = .data$y),
      linetype = 3, linewidth = 0.3, colour = "grey70",
      inherit.aes = FALSE) +
    # Vertical connectors (parent -> child y).
    ggplot2::geom_segment(
      data = seg_v,
      ggplot2::aes(x = .data$x, xend = .data$xend,
                    y = .data$y, yend = .data$yend),
      linewidth = 0.6, colour = "grey45", lineend = "round",
      inherit.aes = FALSE) +
    # Horizontal branch segments, coloured by **duplication density**
    # (D per 100 speciations) — a meaningful rate rather than a raw count
    # that just tracks how many OGs pass through a lineage. Width scales
    # with total events so the busy branches still stand out.
    ggplot2::geom_segment(
      data = edge_p,
      ggplot2::aes(x = .data$x_p, xend = .data$x_c,
                    y = .data$y_c, yend = .data$y_c,
                    colour = .data$dup_rate * 100,
                    linewidth = .data$n_total + 1),
      lineend = "round", inherit.aes = FALSE) +
    ggplot2::scale_colour_gradientn(
      name = "duplication\nper 100 S",
      colours = c("#D9DEE3", "#9A8E76", "#9E3A38", "#6E1E20"),
      values  = scales::rescale(c(0, 0.15, 0.55, 1)),
      guide = ggplot2::guide_colorbar(order = 3, direction = "vertical",
                                        barwidth = 0.32, barheight = 3.6,
                                        title.position = "top",
                                        title.hjust = 0)
    ) +
    # Linewidth legend suppressed — the glossary strip at the bottom of
    # the tree panel already explains "thickness ~ total events".
    ggplot2::scale_linewidth(
      range = c(0.9, 4.0), trans = "log10", guide = "none"
    )

  if (!is.null(bar_new) && nrow(bar_new)) {
    bar_new$xmin <- bar_new$x_mid - bar_halfw
    bar_new$xmax <- bar_new$x_mid + bar_halfw
    # Faint drop-shadow rectangle behind the stack for separation.
    shadow <- do.call(rbind, lapply(split(bar_new, bar_new$edge_i),
                                     function(sub) {
      data.frame(
        xmin = sub$xmin[1] - bar_halfw * 0.25,
        xmax = sub$xmax[1] + bar_halfw * 0.25,
        ymin = min(sub$y0) - 0.04,
        ymax = max(sub$y1) + 0.04
      )
    }))
    p <- p +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_rect(
        data = shadow,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                      ymin = .data$ymin, ymax = .data$ymax),
        fill = "#FFFFFF", colour = "grey85", linewidth = 0.2,
        inherit.aes = FALSE) +
      ggplot2::geom_rect(
        data = bar_new,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                      ymin = .data$y0,   ymax = .data$y1,
                      fill = .data$event),
        colour = "white", linewidth = 0.2,
        inherit.aes = FALSE) +
      ggplot2::scale_fill_manual(
        values = ev_cols, name = "event",
        guide  = ggplot2::guide_legend(
          order = 1, direction = "vertical", ncol = 1,
          keywidth  = ggplot2::unit(0.32, "cm"),
          keyheight = ggplot2::unit(0.32, "cm"),
          override.aes = list(colour = "white", size = 3.2))
      )
  }

  # Root gets a gold star; other internal nodes keep the navy badge.
  inner_lab <- node_lab[node_lab$node != root, , drop = FALSE]
  root_lab  <- node_lab[node_lab$node == root, , drop = FALSE]
  # Clade size (tip count in each internal node's subtree) — rendered as
  # a small italic label below the node badge. Helps the viewer see how
  # big a clade is without counting tips manually.
  if (nrow(inner_lab)) {
    inner_lab$clade_n <- vapply(inner_lab$node, function(nd) {
      tryCatch(length(ape::extract.clade(tree, nd)$tip.label),
               error = function(e) NA_integer_)
    }, integer(1))
  }
  p <- p +
    # Internal node badges: filled circle with the node id printed in
    # white inside — avoids the label-overlap problem we had before.
    ggplot2::geom_point(
      data = inner_lab,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, fill = "#1A2238", colour = "white",
      size = 5.2, stroke = 0.4, inherit.aes = FALSE) +
    ggplot2::geom_text(
      data = inner_lab,
      ggplot2::aes(x = .data$x, y = .data$y,
                    label = sub("^N0*", "",
                                 as.character(.data$label))),
      size = 2.3, colour = "white", fontface = "bold",
      inherit.aes = FALSE) +
    ggplot2::geom_text(
      data = inner_lab[!is.na(inner_lab$clade_n), , drop = FALSE],
      ggplot2::aes(x = .data$x + max_x * 0.013, y = .data$y,
                    label = sprintf("n=%d", .data$clade_n)),
      size = 1.9, colour = "grey45", fontface = "italic",
      hjust = 0, vjust = 1.8,
      inherit.aes = FALSE) +
    # Root marker — gold five-pointed star (shape = 8 filled via stroke).
    ggplot2::geom_point(
      data = root_lab,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 8, size = 5.5, stroke = 1.2, colour = "#B38315",
      inherit.aes = FALSE) +
    ggplot2::geom_point(
      data = root_lab,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = 2.2, fill = "#E0B33A", colour = "#7A5D1F",
      stroke = 0.4, inherit.aes = FALSE) +
    ggplot2::geom_text(
      data = root_lab,
      ggplot2::aes(x = .data$x, y = .data$y - 0.55,
                    label = "root"),
      size = 2.4, fontface = "bold", colour = "#7A5D1F",
      inherit.aes = FALSE) +
    # Tip markers at the end of the dotted extension.
    ggplot2::geom_point(
      data = tip_map,
      ggplot2::aes(x = .data$x_end, y = .data$y),
      shape = 21, fill = "#6E1E20", colour = "white",
      size = 2.5, stroke = 0.35, inherit.aes = FALSE)

  # The tip whose root-to-tip path accumulates the most D/T/L events is
  # drawn in a stronger accent colour + bold so the viewer's eye can
  # find the "most affected" lineage without reading every rollup number.
  # Split into two geom_text calls with fixed aesthetics to avoid needing
  # an identity scale that would clash with other colour scales.
  tip_rollup_total <- lineage_rollup$n_D + lineage_rollup$n_T +
                       lineage_rollup$n_L
  hot_tip_node <- if (length(tip_rollup_total) && any(tip_rollup_total > 0))
                    lineage_rollup$node[which.max(tip_rollup_total)]
                  else NA_integer_
  is_hot <- !is.na(hot_tip_node) & tip_map$node == hot_tip_node
  p <- p +
    ggplot2::geom_text(
      data = tip_map[!is_hot, , drop = FALSE],
      ggplot2::aes(x = .data$x_end, y = .data$y, label = .data$label),
      hjust = 0, size = 3.3, nudge_x = max_x * 0.015,
      fontface = "italic", colour = "#1A2238",
      inherit.aes = FALSE)
  if (any(is_hot)) {
    p <- p + ggplot2::geom_text(
      data = tip_map[is_hot, , drop = FALSE],
      ggplot2::aes(x = .data$x_end, y = .data$y, label = .data$label),
      hjust = 0, size = 3.5, nudge_x = max_x * 0.015,
      fontface = "bold.italic", colour = "#9E3A38",
      inherit.aes = FALSE)
  }

  # Gold diamonds + rank number over the top-3 duplication branches.
  top_df <- edge_p[edge_p$is_top, , drop = FALSE]
  if (nrow(top_df)) {
    top_df$rank <- rank(-top_df$n_D, ties.method = "first")
    p <- p +
      ggplot2::geom_point(
        data = top_df,
        ggplot2::aes(x = .data$x_mid, y = .data$y_c + 0.38),
        shape = 23, fill = "#E0B33A", colour = "#7A5D1F",
        stroke = 0.4, size = 3.1, inherit.aes = FALSE) +
      ggplot2::geom_text(
        data = top_df,
        ggplot2::aes(x = .data$x_mid, y = .data$y_c + 0.38,
                      label = .data$rank),
        size = 2.1, fontface = "bold", colour = "#3C2B06",
        inherit.aes = FALSE)

  }

  # Per-tip lineage rollup: miniature horizontal stacks to the right of
  # each italic tip label. Anchored at a fixed x-offset past the label
  # block so all tips align regardless of label width.
  rollup_x0 <- max_x * 1.18
  rollup_len <- max_x * 0.12
  roll_long <- do.call(rbind, lapply(c("D", "T", "L"), function(k) {
    col <- paste0("n_", k)
    data.frame(
      node  = lineage_rollup$node,
      y     = lineage_rollup$y,
      event = factor(k, levels = c("D", "T", "L")),
      n     = lineage_rollup[[col]]
    )
  }))
  roll_long$frac <- roll_long$n / max_roll
  roll_long$width <- roll_long$frac * rollup_len
  roll_long <- roll_long[order(roll_long$node,
                                 as.integer(roll_long$event)), , drop = FALSE]
  roll_long$xmin <- NA_real_; roll_long$xmax <- NA_real_
  for (nd in unique(roll_long$node)) {
    idx <- which(roll_long$node == nd)
    cumw <- c(0, cumsum(roll_long$width[idx]))
    roll_long$xmin[idx] <- rollup_x0 + utils::head(cumw, -1)
    roll_long$xmax[idx] <- rollup_x0 + utils::tail(cumw, -1)
  }
  roll_long <- roll_long[roll_long$n > 0, , drop = FALSE]
  if (nrow(roll_long)) {
    rollup_totals <- data.frame(
      y     = lineage_rollup$y,
      total = lineage_rollup$n_D + lineage_rollup$n_T + lineage_rollup$n_L
    )
    rollup_totals <- rollup_totals[rollup_totals$total > 0, , drop = FALSE]
    rollup_totals$x <- rollup_x0 + rollup_len + max_x * 0.008
    rollup_totals$label <- formatC(rollup_totals$total, format = "d",
                                     big.mark = ",")
    p <- p +
      ggplot2::geom_rect(
        data = roll_long,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                      ymin = .data$y - 0.24, ymax = .data$y + 0.24,
                      fill = .data$event),
        colour = "white", linewidth = 0.2, inherit.aes = FALSE) +
      ggplot2::geom_text(
        data = rollup_totals,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        hjust = 0, size = 2.5, colour = "grey25",
        family = "mono", inherit.aes = FALSE)
    # No new_scale_fill needed — ev_cols reuses the existing event scale.

    # Axis label above the rollup column.
    p <- p +
      ggplot2::annotate(
        "text", x = rollup_x0 + rollup_len / 2, y = n_tip + 0.70,
        label = "lineage rollup\n(root -> tip)",
        size = 2.6, fontface = "bold", colour = "grey30",
        lineheight = 0.95) +
      ggplot2::annotate(
        "segment",
        x = rollup_x0, xend = rollup_x0 + rollup_len,
        y = 0.55, yend = 0.55,
        linewidth = 0.35, colour = "grey55") +
      ggplot2::annotate(
        "text", x = rollup_x0 + rollup_len / 2, y = 0.30,
        label = sprintf("0 -> %s (max across tips)",
                         formatC(max_roll, format = "d", big.mark = ",")),
        size = 2.2, colour = "grey45", fontface = "italic")
  }

  # Scale bar for branch lengths.
  p <- p +
    ggplot2::annotate(
      "segment", x = scale_x0, xend = scale_x1,
      y = scale_y, yend = scale_y,
      linewidth = 0.7, colour = "#1A2238") +
    ggplot2::annotate(
      "text", x = (scale_x0 + scale_x1) / 2, y = scale_y - 0.3,
      label = sprintf("branch length = %g", scale_len),
      size = 2.5, colour = "grey25")

  # Compact symbol glossary row placed to the right of the scale bar,
  # explaining non-obvious markers (root star / top-dup diamonds). Each
  # symbol gets its own left-aligned label to avoid the label overlap
  # that happens when three centered labels share a horizontal line.
  gloss_y  <- scale_y - 0.05
  gloss_x1 <- scale_x1 + max_x * 0.12
  gloss_x2 <- gloss_x1 + max_x * 0.18
  gloss_x3 <- gloss_x2 + max_x * 0.28
  label_pad <- max_x * 0.018
  p <- p +
    ggplot2::annotate(
      "point", x = gloss_x1, y = gloss_y,
      shape = 8, size = 4.4, stroke = 1.0, colour = "#B38315") +
    ggplot2::annotate(
      "point", x = gloss_x1, y = gloss_y,
      shape = 21, size = 1.7, fill = "#E0B33A",
      colour = "#7A5D1F", stroke = 0.35) +
    ggplot2::annotate(
      "text", x = gloss_x1 + label_pad, y = gloss_y,
      label = "root", hjust = 0, size = 2.3,
      colour = "grey25") +
    ggplot2::annotate(
      "point", x = gloss_x2, y = gloss_y,
      shape = 23, fill = "#E0B33A", colour = "#7A5D1F",
      stroke = 0.4, size = 3.1) +
    ggplot2::annotate(
      "text", x = gloss_x2 + label_pad, y = gloss_y,
      label = "top-3 duplication branch", hjust = 0, size = 2.3,
      colour = "grey25") +
    ggplot2::annotate(
      "text", x = gloss_x3, y = gloss_y,
      label = "thickness ~ total events  |  colour = dup density",
      hjust = 0, size = 2.3, colour = "grey25",
      fontface = "italic")

  # Inline per-branch total labels removed: the D/T/L side heatmap on
  # the right panel already carries the raw counts for *every* branch in
  # a cleaner tabular form, so duplicating them on the tree just crowded
  # the branch stacks and the node badges without adding information.
  # The TOP DUPLICATION callout below still marks the rank-1 branch.

  # ---- global-proportion inset (top-left, above the canopy) ------------
  # Horizontal stacked bar showing D / T / L share of all reconciliation
  # events across the tree. Lets the viewer anchor individual branch
  # differences against the whole-tree mix at a glance.
  tot_D <- sum(edge_p$n_D); tot_T <- sum(edge_p$n_T); tot_L <- sum(edge_p$n_L)
  tot_all <- tot_D + tot_T + tot_L
  if (tot_all > 0) {
    gi_x0 <- max_x * 0.03
    gi_x1 <- max_x * 0.42
    gi_y  <- n_tip + 1.55
    gi_h  <- 0.22
    shares <- c(D = tot_D, T = tot_T, L = tot_L) / tot_all
    cuts   <- c(0, cumsum(shares))
    gi_df  <- data.frame(
      xmin = gi_x0 + cuts[-4] * (gi_x1 - gi_x0),
      xmax = gi_x0 + cuts[-1] * (gi_x1 - gi_x0),
      ymin = gi_y - gi_h, ymax = gi_y + gi_h,
      event = factor(c("D", "T", "L"), levels = c("D", "T", "L"))
    )
    # Narrow slices get a compact "<letter>%" label (e.g. "L14%") instead
    # of "L 14%"; the space makes the text too wide to fit inside thin
    # slices and the "%" gets clipped against the next segment.
    gi_labels <- data.frame(
      x     = (gi_df$xmin + gi_df$xmax) / 2,
      y     = gi_y,
      label = ifelse(
        shares >= 0.18,
        sprintf("%s %s%%", c("D", "T", "L"),
                 formatC(round(shares * 100), format = "d")),
        sprintf("%s%s%%", c("D", "T", "L"),
                 formatC(round(shares * 100), format = "d"))),
      bright = shares > 0.08
    )
    gi_labels <- gi_labels[gi_labels$bright, , drop = FALSE]
    p <- p +
      ggplot2::annotate(
        "text", x = gi_x0, y = gi_y + gi_h + 0.35,
        label = "global event share",
        hjust = 0, size = 2.6, fontface = "bold", colour = "#1A2238") +
      ggplot2::annotate(
        "text", x = gi_x1, y = gi_y + gi_h + 0.35,
        label = sprintf("n = %s events",
                         formatC(tot_all, format = "d", big.mark = ",")),
        hjust = 1, size = 2.2, colour = "grey40", fontface = "italic") +
      ggplot2::geom_rect(
        data = gi_df,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                      ymin = .data$ymin, ymax = .data$ymax,
                      fill = .data$event),
        colour = "white", linewidth = 0.3, inherit.aes = FALSE) +
      ggplot2::geom_text(
        data = gi_labels,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        size = 2.5, colour = "white", fontface = "bold",
        inherit.aes = FALSE)
  }

  p +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.06))) +
    ggplot2::scale_y_continuous(
      limits = c(-0.8, n_tip + 2.4), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = "#F4F1EA", colour = NA),
      panel.background = ggplot2::element_rect(fill = "#F4F1EA", colour = NA),
      legend.position   = "right",
      legend.direction  = "vertical",
      legend.box        = "vertical",
      legend.box.just   = "left",
      legend.justification = "top",
      legend.spacing.y  = ggplot2::unit(0.02, "cm"),
      legend.spacing.x  = ggplot2::unit(0.05, "cm"),
      legend.margin     = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.key.size   = ggplot2::unit(0.28, "cm"),
      legend.text       = ggplot2::element_text(size = 6.5, colour = "#1A2238"),
      legend.title      = ggplot2::element_text(size = 7, face = "bold",
                                                   colour = "#1A2238",
                                                   margin = ggplot2::margin(b = 1)),
      plot.margin       = ggplot2::margin(10, 4, 10, 14)
    ) +
    ggplot2::coord_cartesian(clip = "off")
}

.dtl_tree_panel_manual <- function(tree, layout, edge_df, bar_rows, ev_cols) {
  xpos <- layout$xpos; ypos <- layout$ypos
  n_tip <- layout$n_tip; n_node <- layout$n_node
  total <- layout$total; max_x <- layout$max_x

  parents_u <- unique(tree$edge[, 1])
  seg_v <- do.call(rbind, lapply(parents_u, function(p) {
    cy <- ypos[tree$edge[tree$edge[, 1] == p, 2]]
    data.frame(x = xpos[p], xend = xpos[p],
                y = min(cy), yend = max(cy))
  }))

  tip_df <- data.frame(x = xpos[seq_len(n_tip)], y = ypos[seq_len(n_tip)],
                        label = tree$tip.label)
  node_df <- data.frame(
    x = xpos[(n_tip + 1L):total], y = ypos[(n_tip + 1L):total],
    label = sprintf("N%03d", seq_len(n_node))
  )
  bar_halfw <- max_x * 0.008

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg_v,
      ggplot2::aes(x = .data$x, xend = .data$xend,
                    y = .data$y, yend = .data$yend),
      linewidth = 0.45, colour = "grey40", lineend = "round") +
    ggplot2::geom_segment(data = edge_df,
      ggplot2::aes(x = .data$x_p, xend = .data$x_c,
                    y = .data$y, yend = .data$y,
                    linewidth = .data$n_total + 1),
      colour = "grey22", lineend = "round") +
    ggplot2::scale_linewidth(range = c(0.55, 3.4), trans = "log10",
      name = "total events\non branch",
      guide = ggplot2::guide_legend(order = 2,
        override.aes = list(colour = "grey22")))

  if (!is.null(bar_rows) && nrow(bar_rows)) {
    bar_rows$xmin <- bar_rows$x_mid - bar_halfw
    bar_rows$xmax <- bar_rows$x_mid + bar_halfw
    p <- p +
      ggplot2::geom_rect(data = bar_rows,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                      ymin = .data$y0,   ymax = .data$y1,
                      fill = .data$event),
        colour = "white", linewidth = 0.18) +
      ggplot2::scale_fill_manual(values = ev_cols, name = "event",
        guide = ggplot2::guide_legend(order = 1,
          override.aes = list(colour = "white")))
  }

  p +
    ggplot2::geom_point(data = tip_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, fill = "#9A8E76", colour = "white",
      size = 2.2, stroke = 0.35) +
    ggplot2::geom_text(data = tip_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = 0, size = 3.1, nudge_x = max_x * 0.015,
      fontface = "italic", colour = "#1A2238") +
    ggplot2::geom_point(data = node_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, fill = "#1F4858", colour = "white",
      size = 2.6, stroke = 0.35) +
    ggplot2::geom_text(data = node_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 2.3, nudge_y = 0.35, colour = "#1A2238") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.28))) +
    ggplot2::scale_y_continuous(
      limits = c(0.5, n_tip + 0.9), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = "#F4F1EA", colour = NA),
      panel.background = ggplot2::element_rect(fill = "#F4F1EA", colour = NA),
      legend.position  = "right",
      legend.box       = "vertical",
      legend.key.size  = ggplot2::unit(0.4, "cm"),
      legend.text      = ggplot2::element_text(size = 8),
      legend.title     = ggplot2::element_text(size = 9, face = "bold"),
      plot.margin      = ggplot2::margin(10, 6, 10, 14)
    )
}


# ---- panel B: side bar chart of D / T / L per branch ------------------

.dtl_side_bar_panel <- function(edge_df, ev_cols) {
  tree_obj <- attr(edge_df, "tree")
  n_tip <- if (!is.null(tree_obj)) length(tree_obj$tip.label) else 0L

  # Align rows 1-to-1 with the tree panel: only show tip-edges, and use
  # ggtree's tip y-values directly. Internal-edge events are still visible
  # on the tree panel's stacked bars; the TOTAL row at the bottom still
  # aggregates events across all edges (tips + internal).
  edge_df$is_tip <- if (!is.null(tree_obj)) edge_df$child <= n_tip else FALSE
  totals_all <- list(
    D = sum(edge_df$n_D), T = sum(edge_df$n_T), L = sum(edge_df$n_L)
  )
  # Reuse the tree panel's ggtree y-values so rows align 1-to-1. The tree
  # is already ladderized by plot_dtl_branch_events(); ggtree with its
  # default ladderize=TRUE is idempotent on an already-ladderized tree.
  tree_y <- if (requireNamespace("ggtree", quietly = TRUE) && !is.null(tree_obj)) {
    td <- ggtree::ggtree(tree_obj, ladderize = TRUE)$data
    ymap <- stats::setNames(td$y, td$node)
    ymap[as.character(edge_df$child)]
  } else edge_df$y
  edge_df$tree_y <- tree_y
  edge_df <- edge_df[edge_df$is_tip, , drop = FALSE]
  edge_df <- edge_df[order(edge_df$tree_y), , drop = FALSE]
  edge_df$y <- edge_df$tree_y

  # Top-3 duplication-rank branches — the same rows highlighted with gold
  # diamonds on the tree. Tying them back in the heatmap lets the viewer
  # jump between panels.
  top_rank_rows <- head(order(-edge_df$n_D), 3L)
  top_rank_rows <- top_rank_rows[edge_df$n_D[top_rank_rows] > 0]
  top_rank_df <- if (length(top_rank_rows)) {
    data.frame(
      y    = edge_df$y[top_rank_rows],
      rank = seq_along(top_rank_rows)
    )
  } else data.frame(y = numeric(0), rank = integer(0))

  # Clade bands — alternating light shading on the two root subtrees,
  # mirroring the tree panel's dichotomy cue. We build them from the
  # rank-based y so bands sit exactly one row wide.
  clade_rect_df <- data.frame(ymin = numeric(0), ymax = numeric(0),
                               fill = character(0))
  if (!is.null(tree_obj)) {
    ltr <- ape::ladderize(tree_obj)
    root <- n_tip + 1L
    root_children <- ltr$edge[ltr$edge[, 1] == root, 2]
    for (i in seq_along(root_children)) {
      rc <- root_children[i]
      if (rc <= n_tip) {
        children_nodes <- rc
      } else {
        clade_tips   <- ape::extract.clade(ltr, rc)$tip.label
        clade_tip_ids <- match(clade_tips, ltr$tip.label)
        # Also include every internal edge whose child's descendant set
        # lies entirely inside the clade (i.e. its edge belongs to the
        # clade's subtree).
        children_nodes <- c(clade_tip_ids, rc)
        # Walk down the subtree to collect all internal child-node ids.
        stack <- rc
        while (length(stack)) {
          cur <- stack[1]; stack <- stack[-1]
          kids <- ltr$edge[ltr$edge[, 1] == cur, 2]
          children_nodes <- c(children_nodes, kids)
          stack <- c(stack, kids[kids > n_tip])
        }
        children_nodes <- unique(children_nodes)
      }
      ys <- edge_df$y[edge_df$child %in% children_nodes]
      if (length(ys) > 0 && all(is.finite(ys))) {
        clade_rect_df <- rbind(clade_rect_df, data.frame(
          ymin = min(ys) - 0.5, ymax = max(ys) + 0.5,
          fill = if (i %% 2 == 1) "#F3F0EA" else "#EAE5D6",
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  long <- do.call(rbind, lapply(c("D", "T", "L"), function(k) {
    col <- paste0("n_", k)
    data.frame(branch = edge_df$label,
                y      = edge_df$y,
                event  = factor(k, levels = c("D", "T", "L")),
                n      = edge_df[[col]])
  }))

  # Append a TOTAL row below the branches. Totals span *all* edges (tips
  # and internal), captured before we filtered edge_df to tip rows above.
  # total_y sits inside the shared y range (see tree panel: -0.8..n_tip+2.4)
  # so the row aligns with the tree panel's tip y coordinates.
  total_y <- -0.3
  total_row <- data.frame(
    branch = "TOTAL",
    y      = total_y,
    event  = factor(c("D", "T", "L"), levels = c("D", "T", "L")),
    n      = c(totals_all$D, totals_all$T, totals_all$L)
  )
  long <- rbind(long, total_row)

  y_breaks <- c(edge_df$y, total_y)
  # ASCII "TOTAL" avoids MBCS conversion failures when the user's locale
  # cannot encode U+2211 (Sigma) on the default PDF device.
  y_labels <- c(edge_df$label, "TOTAL")
  y_ord    <- order(y_breaks)
  y_breaks <- y_breaks[y_ord]; y_labels <- y_labels[y_ord]

  if (!sum(long$n, na.rm = TRUE)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
               label = "no D / T / L events",
               colour = "grey55", fontface = "italic", size = 4) +
             ggplot2::theme(
               plot.background = ggplot2::element_rect(fill = "#F4F1EA",
                                                         colour = NA)))
  }

  long$log_n <- log10(long$n + 1)
  long$text  <- ifelse(long$n > 0,
                        formatC(long$n, format = "d", big.mark = ","),
                        "\u00b7")
  long$text_colour <- ifelse(long$log_n > log10(max(long$n, 10) + 1) * 0.55,
                              "white", "#1A2238")

  # Continuous x lets us freely place decorations (dominance dot, column
  # header bar, row-total sparkline) at arbitrary positions beyond the
  # main 3-column matrix. Core tiles sit at x = 1 / 2 / 3.
  long$ev_x <- c(D = 1, T = 2, L = 3)[as.character(long$event)]

  # Hotspot flag: tiles above the 95th percentile of non-zero counts per
  # event type get a gold border — visually ties back to the gold-diamond
  # top-K markers in the tree panel.
  is_total <- long$branch == "TOTAL"
  long$hotspot <- FALSE
  for (ev in c("D", "T", "L")) {
    idx <- which(long$event == ev & long$n > 0 & !is_total)
    if (length(idx) >= 4) {
      thr <- stats::quantile(long$n[idx], 0.95, names = FALSE)
      long$hotspot[idx] <- long$n[idx] >= thr
    }
  }
  hotspots <- long[long$hotspot, , drop = FALSE]

  # Alternating row striping for readability across many rows.
  uniq_y <- sort(unique(long$y[!is_total]))
  stripes <- data.frame(
    y = uniq_y,
    band = as.integer(seq_along(uniq_y) %% 2 == 0)
  )
  stripes <- stripes[stripes$band == 1L, , drop = FALSE]

  total_band <- data.frame(y = total_y)

  # ---- decorations: dominance / column headers / row sparkline ---------
  per_branch <- split(long[!is_total, , drop = FALSE],
                       factor(long$branch[!is_total],
                               unique(long$branch[!is_total])))
  dominance <- do.call(rbind, lapply(per_branch, function(df) {
    if (sum(df$n) == 0) return(NULL)
    i <- which.max(df$n)
    data.frame(y = df$y[1],
                dom_ev = factor(as.character(df$event[i]),
                                 levels = c("D", "T", "L")),
                dom_n  = df$n[i],
                total  = sum(df$n),
                stringsAsFactors = FALSE)
  }))

  # Column headers: colored pill with the event letter + total count above.
  # Single tight row instead of a bar chart above the matrix — the matrix
  # and TOTAL row already convey magnitude; this header just names the
  # columns consistently.
  col_totals_df <- data.frame(
    event  = factor(c("D", "T", "L"), levels = c("D", "T", "L")),
    ev_x   = 1:3,
    n      = c(totals_all$D, totals_all$T, totals_all$L),
    stringsAsFactors = FALSE
  )
  # Keep headers inside the shared y range (n_tip + 2.4 max) so rows align
  # with the tree panel. n_tip + 1.1 / n_tip + 1.9 leaves room at the top.
  header_pill_y  <- n_tip + 1.1
  header_count_y <- n_tip + 1.9
  col_totals_df$label  <- formatC(col_totals_df$n, format = "d",
                                    big.mark = ",")

  # Row-total sparkline (on the right of the matrix).
  dom_x    <- 0.30
  spark_x0 <- 3.55
  spark_x1 <- 4.12
  if (!is.null(dominance) && nrow(dominance)) {
    max_rt <- max(dominance$total, 1)
    dominance$spark_xmax <- spark_x0 +
      (dominance$total / max_rt) * (spark_x1 - spark_x0)
  }

  p <- ggplot2::ggplot(long,
                        ggplot2::aes(x = .data$ev_x, y = .data$y,
                                      fill = .data$log_n))

  # Clade bands (mirror the tree panel's root dichotomy shading). Each
  # rectangle is added as its own layer with a hardcoded fill so we don't
  # consume the fill scale slot used by the tile gradient.
  if (nrow(clade_rect_df)) {
    for (i in seq_len(nrow(clade_rect_df))) {
      cb <- clade_rect_df[i, , drop = FALSE]
      p <- p + ggplot2::geom_rect(
        data = cb, inherit.aes = FALSE,
        ggplot2::aes(xmin = -Inf, xmax = Inf,
                      ymin = .data$ymin, ymax = .data$ymax),
        fill = cb$fill, colour = NA, alpha = 0.55)
    }
  } else if (nrow(stripes)) {
    # Fallback: plain alternating grey stripes when we can't compute clades.
    p <- p + ggplot2::geom_rect(
      data = stripes, inherit.aes = FALSE,
      ggplot2::aes(xmin = -Inf, xmax = Inf,
                    ymin = .data$y - 0.44, ymax = .data$y + 0.44),
      fill = "#F1EEE8", colour = NA)
  }

  p <- p +
    # TOTAL row band: deeper navy tint + a solid rule on its top edge so
    # the row reads as a distinct summary strip rather than just another
    # data row.
    ggplot2::geom_rect(
      data = total_band, inherit.aes = FALSE,
      ggplot2::aes(xmin = -Inf, xmax = Inf,
                    ymin = .data$y - 0.48, ymax = .data$y + 0.48),
      fill = "#1A2238", colour = NA, alpha = 0.14) +
    ggplot2::geom_segment(
      data = total_band, inherit.aes = FALSE,
      ggplot2::aes(x = -Inf, xend = Inf,
                    y = .data$y + 0.48, yend = .data$y + 0.48),
      colour = "#1A2238", linewidth = 0.45) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6,
                        width = 0.9, height = 0.8) +
    ggplot2::geom_tile(
      data = hotspots, inherit.aes = FALSE,
      ggplot2::aes(x = .data$ev_x, y = .data$y),
      fill = NA, colour = "#E0B33A", linewidth = 0.85,
      width = 0.9, height = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$text, colour = .data$text_colour),
      size = 2.9, family = "mono") +
    ggplot2::scale_fill_gradientn(
      name   = "log10(count+1)",
      colours = c("#F3EEE6", "#E9D8AE", "#D09B60", "#9E3A38", "#6E1E20"),
      values  = scales::rescale(c(0, 0.15, 0.45, 0.75, 1)),
      guide   = ggplot2::guide_colorbar(direction = "vertical",
                                          barwidth = 0.32, barheight = 3.6,
                                          title.position = "top",
                                          title.hjust = 0)
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_hline(yintercept = total_y + 0.55,
                         colour = "grey70", linewidth = 0.3,
                         linetype = 2)

  # Switch to event-colour fill for the decoration layers.
  if (requireNamespace("ggnewscale", quietly = TRUE)) {
    p <- p + ggnewscale::new_scale_fill()

    # Column header: colored pill with event letter in white, plus the
    # column's total count above. Fixed-size so layout is stable.
    p <- p +
      ggplot2::geom_point(
        data = col_totals_df, inherit.aes = FALSE,
        ggplot2::aes(x = .data$ev_x, y = header_pill_y,
                      fill = .data$event),
        shape = 21, size = 9, stroke = 0.5, colour = "white") +
      ggplot2::geom_text(
        data = col_totals_df, inherit.aes = FALSE,
        ggplot2::aes(x = .data$ev_x, y = header_pill_y,
                      label = as.character(.data$event)),
        size = 4.2, fontface = "bold", colour = "white") +
      ggplot2::geom_text(
        data = col_totals_df, inherit.aes = FALSE,
        ggplot2::aes(x = .data$ev_x, y = header_count_y,
                      label = .data$label),
        size = 2.6, fontface = "bold", colour = "#1A2238")

    # Dominance dot per active branch row — no caption; meaning goes in
    # the side-panel subtitle.
    if (!is.null(dominance) && nrow(dominance)) {
      p <- p +
        ggplot2::geom_point(
          data = dominance, inherit.aes = FALSE,
          ggplot2::aes(x = dom_x, y = .data$y, fill = .data$dom_ev),
          shape = 21, size = 2.4, stroke = 0.35, colour = "white")
    }

    p <- p + ggplot2::scale_fill_manual(values = ev_cols, guide = "none")
  }

  # Row-total sparkline — background rail + filled bar.
  if (!is.null(dominance) && nrow(dominance)) {
    rail_df <- data.frame(y = dominance$y)
    p <- p +
      ggplot2::geom_rect(
        data = rail_df, inherit.aes = FALSE,
        ggplot2::aes(xmin = spark_x0, xmax = spark_x1,
                      ymin = .data$y - 0.09, ymax = .data$y + 0.09),
        fill = "grey88", colour = NA) +
      ggplot2::geom_rect(
        data = dominance, inherit.aes = FALSE,
        ggplot2::aes(xmin = spark_x0, xmax = .data$spark_xmax,
                      ymin = .data$y - 0.22, ymax = .data$y + 0.22),
        fill = "#6E1E20", colour = "white", linewidth = 0.2)
  }

  # Top-3 duplication rank markers — gold diamonds to the very left of
  # each corresponding heatmap row, ties back to the tree panel.
  rank_x <- -0.18
  if (nrow(top_rank_df)) {
    p <- p +
      ggplot2::geom_point(
        data = top_rank_df, inherit.aes = FALSE,
        ggplot2::aes(x = rank_x, y = .data$y),
        shape = 23, fill = "#E0B33A", colour = "#7A5D1F",
        stroke = 0.4, size = 2.6) +
      ggplot2::geom_text(
        data = top_rank_df, inherit.aes = FALSE,
        ggplot2::aes(x = rank_x, y = .data$y, label = .data$rank),
        size = 1.9, fontface = "bold", colour = "#3C2B06")
  }

  # Per-row text labels, coloured by row type: tips stay navy and bold,
  # internal-node rows (N0*) go a lighter slate + regular weight, TOTAL
  # stays bold navy. Done via geom_text so the three styles can coexist
  # — axis.text.y.right in ggplot only accepts a single style.
  label_x <- spark_x1 + 0.20
  row_lab_df <- data.frame(
    y      = c(edge_df$y, total_y),
    label  = c(edge_df$label, "TOTAL"),
    kind   = c(ifelse(edge_df$is_tip, "tip", "node"), "total"),
    stringsAsFactors = FALSE
  )
  row_lab_df$colour <- c(tip   = "#1A2238",
                          node  = "#525866",
                          total = "#1A2238")[row_lab_df$kind]
  row_lab_df$face <- c(tip   = "bold",
                        node  = "plain",
                        total = "bold")[row_lab_df$kind]
  for (kind_i in c("tip", "node", "total")) {
    sub_df <- row_lab_df[row_lab_df$kind == kind_i, , drop = FALSE]
    if (!nrow(sub_df)) next
    p <- p + ggplot2::geom_text(
      data = sub_df, inherit.aes = FALSE,
      ggplot2::aes(x = label_x, y = .data$y, label = .data$label),
      hjust = 0, size = 2.6, family = "mono",
      fontface = sub_df$face[1], colour = sub_df$colour[1])
  }

  # Final scales & theme. Rows use the tree panel's tip y-values so rows
  # line up 1-to-1; y range matches the tree panel (-0.8 .. n_tip + 2.4).
  y_bot <- -0.8
  y_top <- n_tip + 2.4
  p +
    ggplot2::scale_x_continuous(
      limits = c(-0.35, spark_x1 + 1.45),
      expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      limits = c(y_bot, y_top),
      expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.background    = ggplot2::element_rect(fill = "#F4F1EA", colour = NA),
      panel.background   = ggplot2::element_rect(fill = "#F4F1EA", colour = NA),
      legend.position    = "right",
      legend.direction   = "vertical",
      legend.box         = "vertical",
      legend.justification = "top",
      legend.spacing.y   = ggplot2::unit(0.02, "cm"),
      legend.spacing.x   = ggplot2::unit(0.05, "cm"),
      legend.margin      = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin  = ggplot2::margin(0, 0, 0, 0),
      legend.key.size    = ggplot2::unit(0.28, "cm"),
      legend.text        = ggplot2::element_text(size = 6.5),
      legend.title       = ggplot2::element_text(size = 7, face = "bold",
                                                   colour = "#1A2238",
                                                   margin = ggplot2::margin(b = 1)),
      plot.margin        = ggplot2::margin(10, 4, 10, 4)
    ) +
    ggplot2::coord_cartesian(clip = "off")
}
