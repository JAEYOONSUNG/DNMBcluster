#' Root-HOG occupancy / copy-number heatmap
#'
#' Side-by-side panel: rooted species tree (left) + HOG copy-number
#' heatmap (right). Core HOGs render as solid horizontal bands across
#' all tips; accessory HOGs as sparse columns. HOGs with prevalence
#' below `min_prevalence` are dropped to keep the long singleton tail
#' from washing out the figure.
#'
#' Anvi'o palette / theme via `.dnmb_anvio_pal` for visual consistency
#' with the rest of the plotting family.
#'
#' @param out_dir `run_orthofinder_like()` output directory.
#' @param out_pdf Destination PDF.
#' @param min_prevalence Minimum leaf count for a HOG to be plotted.
#' @param sort_by `"prevalence"` (core-like first) or `"cluster_id"`.
#' @param width,height Override auto-scaled PDF dimensions.
#' @param verbose Echo progress.
#' @return Path to the written PDF (invisibly), or `NULL` on failure.
#' @export
plot_hog_occupancy_heatmap <- function(out_dir,
                                       out_pdf = NULL,
                                       min_prevalence = 2L,
                                       sort_by = c("prevalence", "cluster_id"),
                                       width = NULL,
                                       height = NULL,
                                       verbose = TRUE) {
  sort_by <- match.arg(sort_by)
  pal <- .dnmb_anvio_pal

  sp_nwk <- file.path(out_dir, "species_tree_rooted.nwk")
  hog_tsv <- file.path(out_dir, "HOGs", "N001_root.tsv")
  if (!file.exists(sp_nwk) || !file.exists(hog_tsv)) {
    warning("[hog_occupancy] missing ", sp_nwk, " or ", hog_tsv, " -- skipping")
    return(NULL)
  }
  for (pkg in c("ape", "ggplot2", "patchwork")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(sprintf("[hog_occupancy] requires %s; skipping", pkg))
      return(NULL)
    }
  }

  tree <- ape::ladderize(ape::read.tree(sp_nwk))
  tip_order <- tree$tip.label

  hogs_df <- utils::read.table(hog_tsv, sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE,
                                comment.char = "#", quote = "")
  if (!nrow(hogs_df)) {
    warning("[hog_occupancy] empty HOG table -- skipping")
    return(NULL)
  }

  exploded <- do.call(rbind, lapply(seq_len(nrow(hogs_df)), function(i) {
    keys <- strsplit(hogs_df$member_genome_keys[i], ",", fixed = TRUE)[[1]]
    if (!length(keys)) return(NULL)
    tab <- table(keys)
    data.frame(
      hog_id = hogs_df$hog_id[i],
      cluster_id = hogs_df$cluster_id[i],
      leaf = names(tab),
      copies = as.integer(tab),
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(exploded) || !nrow(exploded)) {
    warning("[hog_occupancy] no HOG/leaf pairs parsed -- skipping")
    return(NULL)
  }

  prevalence <- stats::aggregate(leaf ~ hog_id, data = exploded,
                                  FUN = function(x) length(unique(x)))
  names(prevalence)[2] <- "n_leaves"
  keep_hogs <- prevalence$hog_id[prevalence$n_leaves >= min_prevalence]
  exploded <- exploded[exploded$hog_id %in% keep_hogs, , drop = FALSE]
  if (!nrow(exploded)) {
    warning(sprintf("[hog_occupancy] no HOGs with prevalence >= %d -- skipping",
                    min_prevalence))
    return(NULL)
  }

  if (sort_by == "prevalence") {
    hog_order <- prevalence[prevalence$hog_id %in% keep_hogs, , drop = FALSE]
    hog_order <- hog_order[order(-hog_order$n_leaves, hog_order$hog_id), ]
    hog_levels <- hog_order$hog_id
  } else {
    ord <- unique(exploded[order(exploded$cluster_id), c("hog_id")])
    hog_levels <- ord
  }

  exploded$leaf   <- factor(exploded$leaf,   levels = tip_order)
  exploded$hog_id <- factor(exploded$hog_id, levels = hog_levels)
  exploded <- exploded[!is.na(exploded$leaf), , drop = FALSE]

  exploded$copies_cap <- pmin(exploded$copies, 5L)
  copies_breaks <- c(1, 2, 3, 4, 5)
  copies_labels <- c("1", "2", "3", "4", ">=5")
  copies_colors <- c(pal$blue_cool, pal$orange, pal$red,
                     pal$red_strong, "#7A1D1D")

  tree_plot <- .build_tree_panel(tree, tip_order, pal)

  hm <- ggplot2::ggplot(exploded,
                        ggplot2::aes(x = .data$hog_id, y = .data$leaf,
                                      fill = factor(.data$copies_cap,
                                                     levels = copies_breaks))) +
    ggplot2::geom_tile(width = 1, height = 0.9) +
    ggplot2::scale_fill_manual(
      name = "copies", values = stats::setNames(copies_colors, copies_breaks),
      labels = copies_labels, drop = FALSE, na.value = pal$bg
    ) +
    ggplot2::scale_y_discrete(limits = tip_order, position = "right") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::labs(x = sprintf("root HOGs (n=%d, prevalence >= %d)",
                               length(hog_levels), min_prevalence),
                   y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
      axis.text.x      = ggplot2::element_blank(),
      axis.ticks.x     = ggplot2::element_blank(),
      axis.title.x     = ggplot2::element_text(size = 9, colour = pal$subtitle,
                                                 margin = ggplot2::margin(t = 6)),
      axis.text.y      = ggplot2::element_text(size = 9, colour = pal$text,
                                                 face = "italic"),
      panel.grid       = ggplot2::element_blank(),
      legend.position  = "bottom",
      legend.key.size  = ggplot2::unit(0.36, "cm"),
      legend.title     = ggplot2::element_text(size = 9, face = "bold",
                                                 colour = pal$text),
      legend.text      = ggplot2::element_text(size = 8.5, colour = pal$text),
      plot.margin      = ggplot2::margin(8, 12, 6, 4)
    )

  combined <- patchwork::wrap_plots(
    tree_plot, hm,
    widths = c(1, 3.5),
    nrow = 1
  ) +
    patchwork::plot_annotation(
      title = "Root-HOG occupancy & copy number",
      subtitle = sprintf("%d leaves x %d HOGs (>= %d leaves)",
                         length(tip_order), length(hog_levels), min_prevalence),
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
        plot.title    = ggplot2::element_text(face = "bold", size = 14,
                                                colour = pal$text),
        plot.subtitle = ggplot2::element_text(size = 10, colour = pal$subtitle)
      )
    ) &
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA)
    )

  if (isFALSE(out_pdf)) return(invisible(combined))
  if (is.null(out_pdf)) out_pdf <- file.path(out_dir, "hog_occupancy_heatmap.pdf")
  if (is.null(width))   width  <- max(8, 4 + length(hog_levels) * 0.01)
  if (is.null(height))  height <- max(4, length(tip_order) * 0.28 + 1.4)
  width  <- min(width, 40)
  height <- min(height, 30)

  ggplot2::ggsave(out_pdf, combined, width = width, height = height,
                   limitsize = FALSE, bg = pal$bg)
  if (verbose) message(sprintf("[hog_occupancy] wrote %s (%d HOGs x %d leaves)",
                                out_pdf, length(hog_levels), length(tip_order)))
  invisible(out_pdf)
}


.build_tree_panel <- function(tree, tip_order, pal = .dnmb_anvio_pal) {
  n_tip <- length(tree$tip.label)
  n_node <- tree$Nnode
  total <- n_tip + n_node
  xpos <- numeric(total); ypos <- numeric(total)
  tip_y <- stats::setNames(seq_along(tip_order), tip_order)
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
  max_x <- max(xpos)

  ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg_h,
                           ggplot2::aes(x = .data$x, xend = .data$xend,
                                         y = .data$y, yend = .data$yend),
                           linewidth = 0.55, colour = pal$text) +
    ggplot2::geom_segment(data = seg_v,
                           ggplot2::aes(x = .data$x, xend = .data$xend,
                                         y = .data$y, yend = .data$yend),
                           linewidth = 0.55, colour = pal$text) +
    ggplot2::geom_text(data = tip_df,
                        ggplot2::aes(x = .data$x, y = .data$y,
                                      label = .data$label),
                        hjust = 0, size = 3.1, nudge_x = max_x * 0.02,
                        colour = pal$text, fontface = "italic") +
    ggplot2::scale_y_continuous(
      limits = c(0.5, length(tip_order) + 0.5),
      breaks = seq_along(tip_order), labels = tip_order,
      expand = c(0, 0)
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.45))) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
      plot.margin      = ggplot2::margin(8, 4, 6, 6)
    )
}
