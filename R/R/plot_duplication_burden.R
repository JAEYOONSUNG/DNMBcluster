#' Duplication burden overview
#'
#' Two-panel summary of where duplications concentrate:
#'
#' \itemize{
#'   \item **Left** — species tree with per-branch duplication counts.
#'         Each branch is labeled with the number of `D` events whose
#'         `sp_lca` lands on that node (tip branch = intra-species
#'         duplication; internal branch = ancestral duplication).
#'   \item **Right** — top-K OGs ranked by `n_D` (horizontal bar chart).
#'         Answers "which gene families are driving the duplication
#'         signal?" — usually a small set of paralog-rich regulatory
#'         / transporter superfamilies.
#' }
#'
#' Relies on the already-written DTL artifacts; if a run was produced
#' with `dtl = FALSE`, `dtl_events.tsv` / `dtl_per_og.tsv` won't exist
#' and this returns `NULL` with a warning.
#'
#' Intentionally complements `dtl_events_overlay.pdf` (which shows all
#' S/D/T/L as pie charts) by zooming in on just the D signal, which is
#' the most actionable for users curating gene families.
#'
#' @param out_dir `run_orthofinder_like()` directory.
#' @param out_pdf Destination PDF. Defaults to
#'   `<out_dir>/duplication_burden.pdf`.
#' @param top_n Number of top OGs to display on the bar panel.
#' @param verbose Echo progress.
#' @return Path to the written PDF (invisibly) or `NULL` on failure.
#' @export
plot_duplication_burden <- function(out_dir,
                                     out_pdf = NULL,
                                     top_n = 25L,
                                     verbose = TRUE) {
  ev_tsv <- file.path(out_dir, "dtl_events.tsv")
  po_tsv <- file.path(out_dir, "dtl_per_og.tsv")
  sp_nwk <- file.path(out_dir, "species_tree_rooted.nwk")
  if (!file.exists(ev_tsv) || !file.exists(po_tsv) || !file.exists(sp_nwk)) {
    warning("[duplication_burden] missing DTL artifacts or species tree -- skipping")
    return(NULL)
  }
  if (!requireNamespace("ape", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    warning("[duplication_burden] requires ape + ggplot2 + patchwork; skipping")
    return(NULL)
  }

  pal <- .dnmb_anvio_pal

  events <- utils::read.table(ev_tsv, sep = "\t", header = TRUE,
                               stringsAsFactors = FALSE)
  per_og <- utils::read.table(po_tsv, sep = "\t", header = TRUE,
                               stringsAsFactors = FALSE)
  tree <- ape::ladderize(ape::read.tree(sp_nwk))

  d_events <- events[events$event == "D", , drop = FALSE]
  if (!nrow(d_events)) {
    if (is.null(out_pdf)) out_pdf <- file.path(out_dir, "duplication_burden.pdf")
    empty <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0.62,
                         label = "No duplication events detected",
                         size = 5.4, colour = pal$text, fontface = "bold") +
      ggplot2::annotate("text", x = 0, y = 0.30,
                         label = sprintf("All %d OGs reconcile to S only.",
                                          nrow(per_og)),
                         size = 3.8, colour = pal$subtitle, fontface = "italic") +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::theme_void() +
      .dnmb_panel_card(pal, c(40, 40, 40, 40)) +
      patchwork::plot_annotation(
        title = "Duplication burden",
        caption = "no D events in dtl_events.tsv -- nothing to rank",
        theme = ggplot2::theme(
          plot.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
          plot.title    = ggplot2::element_text(face = "bold", size = 14,
                                                  colour = pal$text),
          plot.caption  = ggplot2::element_text(size = 8.5,
                                                  colour = pal$subtitle,
                                                  hjust = 1)
        )
      )
    ggplot2::ggsave(out_pdf, empty, width = 11, height = 4.2,
                     limitsize = FALSE, bg = pal$bg)
    if (verbose) message(sprintf(
      "[duplication_burden] wrote %s (empty: 0 D events)", out_pdf))
    return(invisible(out_pdf))
  }

  n_tip <- length(tree$tip.label)
  d_counts <- stats::aggregate(cluster_id ~ sp_lca, data = d_events, FUN = length)
  names(d_counts)[2] <- "n_D"
  # Join back to the full node range so nodes with 0 D's still render
  # as labels (important for the tree panel — missing labels look like
  # a parse error).
  all_nodes <- data.frame(sp_lca = seq_len(n_tip + tree$Nnode))
  d_counts <- merge(all_nodes, d_counts, by = "sp_lca", all.x = TRUE)
  d_counts$n_D[is.na(d_counts$n_D)] <- 0L

  tree_panel <- .build_dup_tree_panel(tree, d_counts, pal) +
                 .dnmb_panel_card(pal, c(8, 10, 8, 10))

  # Top-N OGs by n_D
  og_sorted <- per_og[order(-per_og$n_D), , drop = FALSE]
  og_top <- utils::head(og_sorted[og_sorted$n_D > 0, , drop = FALSE], top_n)
  if (nrow(og_top)) {
    og_top$label <- sprintf("OG_%07d", og_top$cluster_id)
    og_top$label <- factor(og_top$label, levels = rev(og_top$label))
    bar_panel <- ggplot2::ggplot(og_top,
                                  ggplot2::aes(x = .data$n_D, y = .data$label)) +
      ggplot2::geom_col(fill = pal$red, width = 0.78) +
      ggplot2::geom_text(ggplot2::aes(label = .data$n_D),
                          hjust = -0.2, size = 2.8, colour = pal$text) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.14))
      ) +
      ggplot2::labs(
        title = sprintf("Top %d duplicated OGs", nrow(og_top)),
        x = "duplications per OG (n_D)", y = NULL
      ) +
      .dnmb_anvio_theme(base_size = 10) +
      ggplot2::theme(
        plot.title         = ggplot2::element_text(face = "bold", size = 12,
                                                    colour = pal$text),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.text.y        = ggplot2::element_text(family = "mono", size = 8,
                                                    colour = pal$text)
      ) +
      .dnmb_panel_card(pal, c(8, 10, 8, 10))
  } else {
    bar_panel <- patchwork::plot_spacer()
  }

  combined <- patchwork::wrap_plots(tree_panel, bar_panel,
                                     widths = c(1, 1), nrow = 1) +
    patchwork::plot_annotation(
      title = "Duplication burden",
      subtitle = sprintf(
        "%d total duplications across %d OGs (mean %.2f / OG with D > 0)",
        sum(per_og$n_D),
        sum(per_og$n_D > 0),
        if (sum(per_og$n_D > 0)) mean(per_og$n_D[per_og$n_D > 0]) else 0
      ),
      caption = sprintf(
        "left: where duplications fall on the species tree | right: which OGs drive the signal\n%s",
        .dnmb_footer()),
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
        plot.title    = ggplot2::element_text(face = "bold", size = 16,
                                                colour = pal$text,
                                                family = "serif"),
        plot.subtitle = ggplot2::element_text(size = 10, colour = pal$subtitle),
        plot.caption  = ggplot2::element_text(size = 8, colour = pal$subtitle,
                                                hjust = 1, lineheight = 1.4)
      )
    ) &
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = pal$bg,
                                                             colour = NA))

  if (is.null(out_pdf)) out_pdf <- file.path(out_dir, "duplication_burden.pdf")
  ggplot2::ggsave(out_pdf, combined,
                   width = max(10, 6 + n_tip * 0.2),
                   height = max(6, top_n * 0.22 + 2),
                   limitsize = FALSE, bg = pal$bg)
  if (verbose) message(sprintf(
    "[duplication_burden] wrote %s (%d D events, %d OGs with D>0)",
    out_pdf, sum(per_og$n_D), sum(per_og$n_D > 0)
  ))
  invisible(out_pdf)
}


# Species tree with per-branch duplication counts labeled at each node.
# Manual ggplot build (vs. ggtree) to sidestep the known S7/patchwork
# aesthetics bug — shares the approach used by plot_hog_occupancy.
.build_dup_tree_panel <- function(tree, d_counts, pal = .dnmb_anvio_pal) {
  n_tip <- length(tree$tip.label)
  n_node <- tree$Nnode
  total <- n_tip + n_node
  xpos <- numeric(total); ypos <- numeric(total)
  tip_y <- stats::setNames(seq_len(n_tip), tree$tip.label)
  ypos[seq_len(n_tip)] <- tip_y[tree$tip.label]
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

  seg_h <- data.frame(
    x = xpos[tree$edge[, 1]], xend = xpos[tree$edge[, 2]],
    y = ypos[tree$edge[, 2]], yend = ypos[tree$edge[, 2]]
  )
  parents <- unique(tree$edge[, 1])
  seg_v <- do.call(rbind, lapply(parents, function(p) {
    cy <- ypos[tree$edge[tree$edge[, 1] == p, 2]]
    data.frame(x = xpos[p], xend = xpos[p], y = min(cy), yend = max(cy))
  }))

  # D count per node — size & color mapped. Zero-count nodes render
  # small + pale so they don't visually dominate.
  node_df <- data.frame(
    x = xpos, y = ypos, sp_lca = seq_len(total)
  )
  node_df <- merge(node_df, d_counts, by = "sp_lca", all.x = TRUE)
  node_df$n_D[is.na(node_df$n_D)] <- 0L

  tip_df <- data.frame(
    x = xpos[seq_len(n_tip)], y = ypos[seq_len(n_tip)],
    label = tree$tip.label
  )
  max_x <- max(xpos)
  leader_x <- max_x * 1.04
  tip_lead <- data.frame(
    x = tip_df$x, xend = leader_x,
    y = tip_df$y, yend = tip_df$y
  )
  tip_df$leader_x <- leader_x

  ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg_h,
                           ggplot2::aes(x = .data$x, xend = .data$xend,
                                         y = .data$y, yend = .data$yend),
                           linewidth = 0.5, colour = pal$text) +
    ggplot2::geom_segment(data = seg_v,
                           ggplot2::aes(x = .data$x, xend = .data$xend,
                                         y = .data$y, yend = .data$yend),
                           linewidth = 0.5, colour = pal$text) +
    ggplot2::geom_point(data = node_df,
                         ggplot2::aes(x = .data$x, y = .data$y,
                                       size = .data$n_D,
                                       fill = .data$n_D),
                         shape = 21, colour = pal$bg, stroke = 0.4) +
    ggplot2::geom_text(data = node_df[node_df$n_D > 0, , drop = FALSE],
                        ggplot2::aes(x = .data$x, y = .data$y,
                                      label = .data$n_D),
                        size = 2.6, nudge_y = 0.32,
                        colour = pal$text, fontface = "bold") +
    ggplot2::geom_segment(data = tip_lead,
                           ggplot2::aes(x = .data$x, xend = .data$xend,
                                         y = .data$y, yend = .data$yend),
                           linetype = "dotted", linewidth = 0.25,
                           colour = pal$rule) +
    ggplot2::geom_text(data = tip_df,
                        ggplot2::aes(x = .data$leader_x, y = .data$y,
                                      label = .data$label),
                        hjust = 0, size = 3.2, nudge_x = max_x * 0.015,
                        fontface = "italic", colour = pal$text) +
    ggplot2::scale_size_area(max_size = 12, guide = "none") +
    ggplot2::scale_fill_gradientn(
      name = "n_D",
      colours = c(pal$surface, pal$orange, pal$red, pal$red_strong),
      values  = scales::rescale(c(0, 0.25, 0.65, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 6, barheight = 0.32,
                                        title.position = "top",
                                        title.hjust = 0.5)
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.5))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.5, n_tip + 0.7), expand = c(0, 0)
    ) +
    ggplot2::labs(title = "Duplications per species-tree branch",
                   subtitle = "node size + colour = duplication count at LCA") +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 12,
                                              colour = pal$text),
      plot.subtitle = ggplot2::element_text(size = 9, colour = pal$subtitle),
      plot.margin   = ggplot2::margin(8, 8, 8, 12),
      legend.position = "bottom",
      legend.title    = ggplot2::element_text(size = 8.5, face = "bold",
                                                colour = pal$text),
      legend.text     = ggplot2::element_text(size = 8, colour = pal$text)
    )
}
