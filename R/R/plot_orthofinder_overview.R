#' Multi-panel overview of a `run_orthofinder_like()` directory
#'
#' One PDF that tells a first-time reader what came out of the run:
#'
#' \itemize{
#'   \item **A** rooted species tree (branch lengths + internal node IDs).
#'   \item **B** OG size distribution (log-x histogram).
#'   \item **C** HOG counts per internal species-tree node.
#'   \item **D** Per-branch loss totals from `dtl_losses.tsv` when
#'         `dtl = TRUE` was supplied; placeholder otherwise.
#' }
#'
#' Anvi'o palette / theme via `.dnmb_anvio_pal` so the panel composes
#' cleanly with the rest of the plotting family.
#'
#' @param out_dir Directory produced by `run_orthofinder_like()`.
#' @param out_pdf Destination PDF.
#' @param verbose Echo progress.
#' @return Path to the written PDF (invisibly) or `NULL` on failure.
#' @export
plot_orthofinder_overview <- function(out_dir,
                                       out_pdf = NULL,
                                       verbose = TRUE) {
  sp_nwk   <- file.path(out_dir, "species_tree_rooted.nwk")
  ogs_tsv  <- file.path(out_dir, "orthogroup_trees.tsv")
  hog_dir  <- file.path(out_dir, "HOGs")
  loss_tsv <- file.path(out_dir, "dtl_losses.tsv")

  if (!file.exists(sp_nwk)) {
    warning("[orthofinder_overview] missing species tree -- skipping")
    return(NULL)
  }
  for (pkg in c("ape", "ggplot2", "patchwork")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(sprintf("[orthofinder_overview] requires %s; skipping", pkg))
      return(NULL)
    }
  }

  pal  <- .dnmb_anvio_pal
  tree <- ape::ladderize(ape::read.tree(sp_nwk))

  og_df <- if (file.exists(ogs_tsv))
    utils::read.table(ogs_tsv, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE)
    else NULL
  n_og_trees <- if (!is.null(og_df))
    sum(!is.na(og_df$tree_path)) else 0L
  n_sc <- if (!is.null(og_df) && "single_copy" %in% names(og_df))
    sum(og_df$single_copy %in% TRUE) else 0L
  n_hog_files <- length(list.files(hog_dir, pattern = "\\.tsv$"))
  total_losses <- if (file.exists(loss_tsv)) {
    lt <- utils::read.table(loss_tsv, sep = "\t", header = TRUE,
                             stringsAsFactors = FALSE)
    if (nrow(lt)) sum(lt$loss_count, na.rm = TRUE) else 0L
  } else NA_integer_

  header <- .ovw_header(length(tree$tip.label), n_og_trees, n_sc,
                         n_hog_files, total_losses, pal)

  p_tree  <- .build_overview_tree_panel(tree, pal) +
             ggplot2::labs(tag = "A") + .dnmb_panel_card(pal) +
             .dnmb_tag_theme(pal)
  p_sizes <- .build_og_size_panel(ogs_tsv, pal)    +
             ggplot2::labs(tag = "B") + .dnmb_panel_card(pal) +
             .dnmb_tag_theme(pal)
  p_hogs  <- .build_hog_node_panel(hog_dir, tree, pal) +
             ggplot2::labs(tag = "C") + .dnmb_panel_card(pal) +
             .dnmb_tag_theme(pal)
  p_loss  <- .build_loss_panel(loss_tsv, tree, pal)    +
             ggplot2::labs(tag = "D") + .dnmb_panel_card(pal) +
             .dnmb_tag_theme(pal)

  grid <- patchwork::wrap_plots(
    p_tree, p_sizes,
    p_hogs, p_loss,
    nrow = 2, ncol = 2,
    heights = c(1, 0.95)
  )

  combined <- patchwork::wrap_elements(full = header) /
              patchwork::wrap_elements(full = grid) +
    patchwork::plot_layout(heights = c(0.18, 1.95)) +
    patchwork::plot_annotation(
      title = "OrthoFinder-parity run overview",
      caption = sprintf(
        "header: at-a-glance totals | A: rooted species tree | B: OG size distribution | C: HOGs per node | D: per-branch losses\n%s",
        .dnmb_footer()),
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
        plot.title    = ggplot2::element_text(face = "bold", size = 17,
                                                colour = pal$text,
                                                family = "serif"),
        plot.caption  = ggplot2::element_text(size = 8,
                                                colour = pal$subtitle,
                                                hjust = 1, lineheight = 1.4),
        plot.margin   = ggplot2::margin(12, 14, 10, 14)
      )
    ) &
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA)
    )

  if (is.null(out_pdf)) out_pdf <- file.path(out_dir, "orthofinder_overview.pdf")
  ggplot2::ggsave(
    out_pdf, combined,
    width  = max(12, 7 + length(tree$tip.label) * 0.15),
    height = 10.4,
    limitsize = FALSE,
    bg = pal$bg
  )
  if (verbose) message(sprintf(
    "[orthofinder_overview] wrote %s", out_pdf
  ))
  invisible(out_pdf)
}


# ---- internals ---------------------------------------------------------

.ovw_node_label_layer <- function(node_df, pal, max_x) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    return(ggrepel::geom_text_repel(
      data = node_df,
      mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 2.6, fontface = "bold", colour = pal$navy,
      segment.colour = pal$rule, segment.size = 0.25,
      box.padding = 0.35, point.padding = 0.18,
      min.segment.length = 0, max.overlaps = Inf, seed = 1L
    ))
  }
  ggplot2::geom_text(
    data = node_df,
    mapping = ggplot2::aes(x = .data$x,
                              y = .data$y + .data$nudge_y,
                              label = .data$label,
                              vjust = .data$vjust),
    size = 2.5, colour = pal$navy, fontface = "bold"
  )
}

.ovw_header <- function(n_tips, n_ogs, n_sc, n_hogs, n_loss, pal) {
  keys <- c("T", "OG", "SC", "H")
  vals <- c(n_tips, n_ogs, n_sc, n_hogs)
  labs <- c("Tree tips", "OGs w/ trees", "Single-copy", "HOG tables")
  fills <- c(pal$navy, pal$blue_cool, pal$orange, pal$red)
  if (!is.na(n_loss)) {
    keys  <- c(keys, "L")
    vals  <- c(vals, n_loss)
    labs  <- c(labs, "DTL losses")
    fills <- c(fills, pal$red_strong)
  }
  stat_df <- data.frame(key = keys, value = vals, label = labs,
                         fill = fills, stringsAsFactors = FALSE)
  .dnmb_stat_band(stat_df, pal)
}

.build_overview_tree_panel <- function(tree, pal) {
  n_tip  <- length(tree$tip.label)
  n_node <- tree$Nnode
  total  <- n_tip + n_node

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

  tip_df <- data.frame(
    x = xpos[seq_len(n_tip)], y = ypos[seq_len(n_tip)],
    label = tree$tip.label
  )
  node_df <- data.frame(
    x = xpos[(n_tip + 1L):total],
    y = ypos[(n_tip + 1L):total],
    label = sprintf("N%03d", seq_len(n_node))
  )
  max_x <- max(xpos)
  leader_x <- max_x * 1.02
  tip_lead <- data.frame(
    x = tip_df$x, xend = leader_x,
    y = tip_df$y, yend = tip_df$y
  )
  tip_df$leader_x <- leader_x
  # Stagger labels: alternate above/below the node so densely packed
  # internal nodes near the root don't overlap.
  ord <- order(node_df$x, node_df$y)
  side <- rep_len(c(1, -1), n_node)[order(ord)]
  node_df$nudge_y <- side * 0.36
  node_df$vjust   <- ifelse(side > 0, 0, 1)

  ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg_h,
                           ggplot2::aes(x = .data$x, xend = .data$xend,
                                         y = .data$y, yend = .data$yend),
                           linewidth = 0.55, colour = pal$text) +
    ggplot2::geom_segment(data = seg_v,
                           ggplot2::aes(x = .data$x, xend = .data$xend,
                                         y = .data$y, yend = .data$yend),
                           linewidth = 0.55, colour = pal$text) +
    ggplot2::geom_point(data = node_df,
                         ggplot2::aes(x = .data$x, y = .data$y),
                         shape = 21, fill = pal$navy, colour = pal$bg,
                         size = 3.4, stroke = 0.5) +
    .ovw_node_label_layer(node_df, pal, max_x) +
    ggplot2::geom_segment(data = tip_lead,
                           ggplot2::aes(x = .data$x, xend = .data$xend,
                                         y = .data$y, yend = .data$yend),
                           linetype = "dotted", linewidth = 0.25,
                           colour = pal$rule) +
    ggplot2::geom_text(data = tip_df,
                        ggplot2::aes(x = .data$leader_x, y = .data$y,
                                      label = .data$label),
                        hjust = 0, size = 3.3, nudge_x = max_x * 0.015,
                        fontface = "italic", colour = pal$text) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.55))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.5, n_tip + 0.7), expand = c(0, 0)
    ) +
    ggplot2::labs(title = "Rooted species tree",
                   subtitle = "branch lengths as computed; N### = internal node id") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
      plot.title       = ggplot2::element_text(face = "bold", size = 12,
                                                 colour = pal$text),
      plot.subtitle    = ggplot2::element_text(size = 9, colour = pal$subtitle),
      plot.margin      = ggplot2::margin(8, 6, 6, 14)
    )
}

.build_og_size_panel <- function(ogs_tsv, pal) {
  if (!file.exists(ogs_tsv)) {
    return(.placeholder_panel("OG size distribution",
                               "orthogroup_trees.tsv not found", pal))
  }
  df <- utils::read.table(ogs_tsv, sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE)
  df <- df[!is.na(df$n_members) & df$n_members > 0, , drop = FALSE]
  if (!nrow(df)) {
    return(.placeholder_panel("OG size distribution", "no OGs found", pal))
  }
  med <- stats::median(df$n_members)
  mx  <- max(df$n_members)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$n_members)) +
    ggplot2::geom_histogram(bins = 40, fill = pal$navy,
                             colour = pal$bg, linewidth = 0.18, alpha = 0.95) +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept = med, linetype = 2,
                         colour = pal$red, linewidth = 0.5) +
    ggplot2::annotate("text", x = med, y = Inf,
                       label = sprintf("median = %g", med),
                       hjust = -0.12, vjust = 1.6, size = 3.1,
                       colour = pal$red, fontface = "italic") +
    ggplot2::labs(
      x = "members per OG (log scale)", y = "OG count",
      title = "OG size distribution",
      subtitle = sprintf("%d OGs \u2022 largest = %d members",
                          nrow(df), mx)
    ) +
    .dnmb_anvio_theme(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, colour = pal$subtitle),
      panel.grid.minor = ggplot2::element_blank()
    )
}

.build_hog_node_panel <- function(hog_dir, tree, pal) {
  hog_files <- list.files(hog_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (!length(hog_files)) {
    return(.placeholder_panel("HOGs per internal node",
                               "HOGs/ directory is empty", pal))
  }
  rows <- lapply(hog_files, function(p) {
    n <- max(0L, length(readLines(p)) - 3L)
    node_id <- as.integer(sub(".*N(\\d+)_.*", "\\1", basename(p)))
    data.frame(node_id = node_id, n_hogs = n)
  })
  hd <- do.call(rbind, rows)
  hd <- stats::aggregate(n_hogs ~ node_id, data = hd, FUN = sum)
  hd$label <- sprintf("N%03d", hd$node_id)
  hd$label <- factor(hd$label, levels = hd$label[order(hd$node_id)])

  ggplot2::ggplot(hd, ggplot2::aes(x = .data$label, y = .data$n_hogs)) +
    ggplot2::geom_col(fill = pal$blue_cool, width = 0.74) +
    ggplot2::geom_text(ggplot2::aes(label = .data$n_hogs),
                        vjust = -0.45, size = 2.8, colour = pal$text) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.18))
    ) +
    ggplot2::labs(
      x = "species-tree node (N###)", y = "HOGs at node",
      title = "HOG counts per internal node",
      subtitle = sprintf("%d HOG tables \u2022 total = %d HOGs",
                          length(hog_files), sum(hd$n_hogs))
    ) +
    .dnmb_anvio_theme(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle    = ggplot2::element_text(size = 9, colour = pal$subtitle),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_text(angle = 40, hjust = 1, size = 8)
    )
}

.build_loss_panel <- function(loss_tsv, tree, pal) {
  if (!file.exists(loss_tsv)) {
    return(.placeholder_panel("Per-branch losses",
                               "dtl_losses.tsv not found (run with dtl = TRUE)",
                               pal))
  }
  df <- utils::read.table(loss_tsv, sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE)
  if (!nrow(df)) {
    return(.placeholder_panel("Per-branch losses", "no loss rows", pal))
  }
  agg <- stats::aggregate(loss_count ~ sp_child, data = df, FUN = sum)
  agg <- agg[agg$loss_count > 0, , drop = FALSE]
  if (!nrow(agg)) {
    return(.placeholder_panel("Per-branch losses", "all loss_count == 0", pal))
  }
  n_tip <- length(tree$tip.label)
  agg$label <- ifelse(agg$sp_child <= n_tip,
                      tree$tip.label[agg$sp_child],
                      sprintf("N%03d", agg$sp_child - n_tip))
  agg <- agg[order(-agg$loss_count), , drop = FALSE]
  agg$label <- factor(agg$label, levels = rev(agg$label))

  ggplot2::ggplot(agg, ggplot2::aes(x = .data$loss_count, y = .data$label)) +
    ggplot2::geom_col(fill = pal$red, width = 0.74) +
    ggplot2::geom_text(ggplot2::aes(label = .data$loss_count),
                        hjust = -0.25, size = 2.8, colour = pal$text) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.16))
    ) +
    ggplot2::labs(
      x = "losses on this branch", y = NULL,
      title = "Per-branch gene losses (DTL)",
      subtitle = sprintf("%d branches with \u22651 loss \u2022 total = %d",
                          nrow(agg), sum(agg$loss_count))
    ) +
    .dnmb_anvio_theme(base_size = 11) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle      = ggplot2::element_text(size = 9, colour = pal$subtitle),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = 8, family = "mono",
                                                   colour = pal$text)
    )
}

.placeholder_panel <- function(title, msg, pal = .dnmb_anvio_pal) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg,
                       size = 4, colour = pal$subtitle, fontface = "italic") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
      plot.title       = ggplot2::element_text(face = "bold", size = 12,
                                                 colour = pal$text),
      plot.margin      = ggplot2::margin(10, 10, 10, 10)
    )
}
