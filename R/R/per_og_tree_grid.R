#' Faceted grid of per-orthogroup gene trees
#'
#' Plots the first `top_n` trees from `per_og_trees()` output as a small-
#' multiple ggtree grid. Useful as a QC thumbnail and as a figure in
#' comparative-genomics write-ups.
#'
#' Each panel shows the gene tree with genome-label tips coloured by
#' genome, a scale bar, and a header strip carrying the OG id, member
#' count and a single-copy flag when applicable.
#'
#' @param og_result Tibble returned by `per_og_trees()`.
#' @param top_n How many orthogroups to show (default 24, a 6x4 grid).
#' @param ncol Grid column count.
#' @param prefer Selection order: `"largest"` (most members), `"single_copy"`,
#'   or `"random"`.
#' @param max_tip_label Max tips below which tip labels are drawn.
#'   Larger trees only get tip points. Default 30.
#' @return A ggplot / patchwork object. Caller saves via `ggsave()`.
#' @export
per_og_tree_grid <- function(og_result,
                             top_n = 24L,
                             ncol  = 4L,
                             prefer = c("largest", "single_copy", "random"),
                             max_tip_label = 30L) {
  prefer <- match.arg(prefer)
  for (pkg in c("ggtree", "patchwork", "ggplot2", "ape")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("per_og_tree_grid: '%s' not installed.", pkg))
  }
  pal <- .dnmb_anvio_pal

  ok <- og_result[!is.na(og_result$tree_path), , drop = FALSE]
  if (!nrow(ok)) stop("per_og_tree_grid: no successful trees in og_result.")

  pick <- switch(prefer,
    largest     = ok[order(-ok$n_members), , drop = FALSE],
    single_copy = ok[ok$single_copy, , drop = FALSE],
    random      = ok[sample.int(nrow(ok)), , drop = FALSE]
  )
  pick <- pick[seq_len(min(top_n, nrow(pick))), , drop = FALSE]

  # Stable per-genome palette across panels. Tips are <prot>_g<genome_uid>.
  all_genomes <- unique(unlist(lapply(pick$tree_path, function(p) {
    tips <- tryCatch(ape::read.tree(p)$tip.label, error = function(e) character())
    sub(".*_g", "", tips)
  })))
  pal_seq <- c(pal$navy, pal$red, pal$orange, pal$blue_cool, pal$red_strong,
               "#5C5470", "#586F4C", "#8C6B3F", "#6E3A4A", "#2E5A53",
               "#7C403E", "#3E5970", "#B27839", "#586073", "#7E6748")
  g_colours <- stats::setNames(
    rep(pal_seq, length.out = length(all_genomes)),
    all_genomes
  )

  plots <- lapply(seq_len(nrow(pick)), function(i) {
    tr <- tryCatch(ape::read.tree(pick$tree_path[i]), error = function(e) NULL)
    if (is.null(tr) || !length(tr$tip.label))
      return(.empty_grid_panel(pick$cluster_id[i], pal))
    if (!is.null(tr$edge.length)) tr$edge.length <- pmax(tr$edge.length, 0)
    .single_og_panel(tr, pick[i, , drop = FALSE], g_colours, pal,
                     max_tip_label)
  })

  n_leg <- length(all_genomes)
  leg_ncol <- min(n_leg, 8L)

  n_sc <- sum(isTRUE(pick$single_copy) | pick$single_copy %in% TRUE)
  subtitle_txt <- sprintf(
    "tips colored by genome - SC = single-copy OG (n=%d / %d) - scale bar = substitutions/site",
    n_sc, nrow(pick))

  patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      title    = sprintf("Per-OG gene trees (top %d)", nrow(pick)),
      subtitle = subtitle_txt,
      caption  = .dnmb_footer(),
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
        plot.title      = ggplot2::element_text(face = "bold", size = 16,
                                                  colour = pal$text,
                                                  family = "serif"),
        plot.subtitle   = ggplot2::element_text(size = 9.5,
                                                  colour = pal$subtitle),
        plot.caption    = ggplot2::element_text(size = 8,
                                                  colour = pal$subtitle,
                                                  hjust = 1)
      )
    ) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      legend.position  = "bottom",
      legend.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
      legend.key       = ggplot2::element_rect(fill = pal$bg, colour = NA),
      legend.title     = ggplot2::element_text(size = 9, face = "bold",
                                                 colour = pal$text),
      legend.text      = ggplot2::element_text(size = 8, colour = pal$text,
                                                 family = "mono"),
      legend.margin    = ggplot2::margin(4, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0)
    ) &
    ggplot2::guides(colour = ggplot2::guide_legend(
      title = "genome",
      ncol  = leg_ncol,
      byrow = TRUE,
      override.aes = list(size = 3, label = "", shape = 16)
    ))
}

.single_og_panel <- function(tr, meta, g_colours, pal, max_tip_label) {
  td <- data.frame(
    label  = tr$tip.label,
    genome = sub(".*_g", "", tr$tip.label),
    stringsAsFactors = FALSE
  )
  td$short <- sub("_g.*$", "", td$label)

  show_labels <- length(tr$tip.label) <= max_tip_label
  tree_depth  <- if (!is.null(tr$edge.length)) sum(tr$edge.length) else length(tr$tip.label)
  x_pad       <- if (show_labels) 0.35 else 0.08

  g <- ggtree::ggtree(tr, ladderize = TRUE,
                      colour = pal$rule, linewidth = 0.35)
  g <- ggtree::`%<+%`(g, td)
  g <- g +
    ggtree::geom_tippoint(ggplot2::aes(colour = .data$genome),
                           size = 1.25) +
    ggplot2::scale_colour_manual(values = g_colours, na.value = pal$neutral,
                                   drop = FALSE, name = "genome")

  if (show_labels) {
    g <- g +
      ggtree::geom_tiplab(ggplot2::aes(label = .data$short,
                                         colour = .data$genome),
                           size = 2.1, offset = tree_depth * 0.015,
                           show.legend = FALSE)
  }

  is_sc   <- isTRUE(meta$single_copy)
  title   <- sprintf("OG_%07d", meta$cluster_id)
  subtitle <- sprintf("n=%d%s", meta$n_members,
                       if (is_sc) "  \u2022  SINGLE-COPY" else "")
  title_col    <- if (is_sc) pal$navy else pal$text
  subtitle_col <- if (is_sc) pal$red  else pal$subtitle

  g +
    ggtree::geom_treescale(width = NULL, fontsize = 2.4,
                             color = pal$subtitle, linesize = 0.3,
                             offset = 0.4) +
    ggplot2::coord_cartesian(clip = "off",
                              xlim = c(0, tree_depth * (1 + x_pad))) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      plot.title        = ggplot2::element_text(face = "bold", size = 9.5,
                                                  colour = title_col,
                                                  margin = ggplot2::margin(b = 0)),
      plot.subtitle     = ggplot2::element_text(size = 7.5,
                                                  colour = subtitle_col,
                                                  face = if (is_sc) "bold" else "plain",
                                                  margin = ggplot2::margin(b = 2)),
      plot.margin       = ggplot2::margin(4, 6, 4, 6)
    )
}

.empty_grid_panel <- function(cid, pal) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0,
                       label = sprintf("OG_%07d\n(tree unavailable)", cid),
                       size  = 3, colour = pal$subtitle, fontface = "italic") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA)
    )
}
