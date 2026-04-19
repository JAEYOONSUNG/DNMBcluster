#' dN/dS selection-analysis visualizations (anvi'o-style)
#'
#' Publication-ready plotting helpers for the selection results returned by
#' [run_hyphy_batch()], [run_hyphy_fel()], [run_hyphy_absrel()], and
#' [run_codeml_m0()]. Styled to match DNMBcluster's anvi'o-inspired
#' pangenome look (see `phylo_circular()` / `circos_pangenome()`):
#' navy / orange / red triad on `#FAFAFA` background, italic OG/species
#' labels, tight serif-leaning typography, 300 dpi output.
#'
#' @name plot_selection
NULL


# ----------------- shared anvi'o palette & theme -------------------------

.dnmb_anvio_pal <- list(
  bg          = "#FFFFFF",   # white page background
  surface     = "#FAFAF7",   # near-white card surface (subtle warmth)
  text        = "#1A2238",   # ink navy-black
  subtitle    = "#525866",   # muted slate
  grid        = "#E4E0D6",   # soft sand grid
  rule        = "#9A8E76",   # bronze rule
  navy        = "#1F4858",   # deep teal-navy
  orange      = "#C97B3F",   # burnt sienna
  red         = "#9E3A38",   # cordovan
  red_strong  = "#6E1E20",   # oxblood
  blue_cool   = "#56789A",   # dusty slate-blue
  neutral     = "#8A7E6C"    # warm taupe
)

#' Build a horizontal stat-card header band.
#'
#' Internal helper shared across composite figures (dtl_summary,
#' orthofinder_overview, pangenome_landscape). Accepts a named list of
#' value/colour/label specs and returns a ggplot ready to drop into a
#' patchwork composition above the main panels.
#'
#' @param stats `data.frame` with columns `key`, `value`, `label`, `fill`.
#' @noRd
.dnmb_stat_band <- function(stats, pal = .dnmb_anvio_pal) {
  stats$x <- seq_len(nrow(stats))
  ggplot2::ggplot(stats) +
    # Card body
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$x - 0.46,
                                      xmax = .data$x + 0.46,
                                      ymin = 0.06, ymax = 1,
                                      fill = .data$fill),
                        colour = NA, alpha = 0.10) +
    # Card thin border
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$x - 0.46,
                                      xmax = .data$x + 0.46,
                                      ymin = 0.06, ymax = 1),
                        colour = pal$rule, fill = NA, linewidth = 0.18) +
    # Top accent rule (the colored bar that marks the stat)
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$x - 0.46,
                                      xmax = .data$x + 0.46,
                                      ymin = 0.92, ymax = 1.0,
                                      fill = .data$fill),
                        colour = NA) +
    # Watermark key letter (large, very faded) sitting behind the value
    ggplot2::geom_text(ggplot2::aes(x = .data$x + 0.32, y = 0.50,
                                      label = .data$key,
                                      colour = .data$fill),
                        size = 16, fontface = "bold", alpha = 0.10,
                        family = "serif", hjust = 1, vjust = 0.5) +
    # The value (large serif numeral)
    ggplot2::geom_text(ggplot2::aes(x = .data$x, y = 0.58,
                                      label = format(.data$value,
                                                       big.mark = ","),
                                      colour = .data$fill),
                        size = 8.5, fontface = "bold", family = "serif") +
    # Label
    ggplot2::geom_text(ggplot2::aes(x = .data$x, y = 0.20,
                                      label = .data$label),
                        size = 3.1, colour = pal$subtitle) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.25)) +
    ggplot2::scale_y_continuous(limits = c(-0.02, 1.06), expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
      plot.margin      = ggplot2::margin(6, 10, 6, 10)
    )
}

.dnmb_presence_pal <- function() {
  c(core = "#2C5F7A", accessory = "#F2A766", unique = "#D06461")
}

.dnmb_footer <- function(extra = NULL) {
  v <- tryCatch(as.character(utils::packageVersion("DNMBcluster")),
                 error = function(e) NA_character_)
  parts <- c("DNMBcluster",
             if (!is.na(v)) sprintf("v%s", v),
             format(Sys.Date()),
             extra)
  paste(parts, collapse = "  \u00B7  ")
}

.dnmb_tag_theme <- function(pal = .dnmb_anvio_pal) {
  ggplot2::theme(
    plot.tag          = ggplot2::element_text(face = "bold", size = 12,
                                                 colour = pal$navy),
    plot.tag.position = c(0.012, 0.985)
  )
}

.dnmb_panel_card <- function(pal = .dnmb_anvio_pal, margin = c(10, 12, 10, 12)) {
  ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = pal$surface,
                                                colour = pal$rule,
                                                linewidth = 0.35),
    panel.background = ggplot2::element_rect(fill = pal$surface, colour = NA),
    legend.background = ggplot2::element_rect(fill = pal$surface, colour = NA),
    legend.key       = ggplot2::element_rect(fill = pal$surface, colour = NA),
    plot.margin      = ggplot2::margin(margin[1], margin[2], margin[3], margin[4])
  )
}

.dnmb_anvio_theme <- function(base_size = 12, legend_pos = "right") {
  pal <- .dnmb_anvio_pal
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.grid.major  = ggplot2::element_line(colour = pal$grid, linewidth = 0.25),
      panel.grid.minor  = ggplot2::element_blank(),
      axis.text         = ggplot2::element_text(colour = pal$text, size = base_size - 3),
      axis.title        = ggplot2::element_text(colour = pal$text, size = base_size - 1),
      plot.title        = ggplot2::element_text(colour = pal$text, size = 14,
                                                face  = "bold"),
      plot.subtitle     = ggplot2::element_text(colour = pal$subtitle, size = 9),
      plot.caption      = ggplot2::element_text(colour = pal$subtitle, size = 7,
                                                hjust = 1),
      legend.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
      legend.key        = ggplot2::element_rect(fill = pal$bg, colour = NA),
      legend.position   = legend_pos,
      legend.key.size   = ggplot2::unit(0.35, "cm"),
      legend.text       = ggplot2::element_text(size = 7, colour = pal$text),
      legend.title      = ggplot2::element_text(size = 8, face = "bold",
                                                colour = pal$text),
      strip.text        = ggplot2::element_text(colour = pal$text, size = 9,
                                                face  = "bold"),
      plot.margin       = ggplot2::margin(10, 14, 10, 14)
    )
}


# ----------------- public plotting API -----------------------------------

#' @rdname plot_selection
#' @description `plot_busted_volcano()` — per-OG episodic-selection
#'   volcano plot. X: log10 omega of the positive rate class. Y:
#'   -log10(q-value) after BH correction. Point size scales with the
#'   positive-class weight. Top OGs are labelled with italic `OG_xxxxxxx`.
#' @param busted Tibble returned by [run_hyphy_busted()] or the `busted`
#'   slot of [run_hyphy_batch()]. Must have `p_value`, `omega_positive`,
#'   `weight_positive`.
#' @param q_threshold BH cutoff used to flag selected OGs. Default 0.05.
#' @param label_top_n If > 0, annotate the top-N most significant OGs.
#' @param output_file Optional PDF/PNG path (300 dpi, 9x6.5 in).
#' @return A ggplot object.
#' @export
plot_busted_volcano <- function(busted, q_threshold = 0.05,
                                 label_top_n = 10L, output_file = NULL) {
  pal <- .dnmb_anvio_pal
  if (is.null(busted) || !nrow(busted))
    return(.plot_empty("No BUSTED results available"))
  d <- busted
  if (!"q_value" %in% names(d))
    d$q_value <- stats::p.adjust(d$p_value, method = "BH")
  d <- d[!is.na(d$p_value), , drop = FALSE]
  if (!nrow(d)) return(.plot_empty("No BUSTED p-values available"))
  d$log10_omega <- ifelse(!is.na(d$omega_positive) & d$omega_positive > 0,
                           log10(d$omega_positive), NA_real_)
  d$neg_log10_q <- -log10(pmax(d$q_value, .Machine$double.eps))
  d$selected <- !is.na(d$q_value) & d$q_value <= q_threshold
  d$class <- factor(ifelse(d$selected, "selected", "ns"),
                     levels = c("ns", "selected"))

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$log10_omega,
                                        y = .data$neg_log10_q,
                                        colour = .data$class,
                                        fill   = .data$class)) +
    ggplot2::geom_vline(xintercept = 0, colour = pal$rule, linetype = 2,
                         linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(q_threshold),
                         colour = pal$rule, linetype = 2, linewidth = 0.3) +
    ggplot2::annotate("text",
      x = -Inf, y = -log10(q_threshold),
      label = sprintf("  q = %g", q_threshold),
      hjust = 0, vjust = -0.4, size = 2.8,
      colour = pal$subtitle, fontface = "italic") +
    ggplot2::annotate("text",
      x = 0, y = Inf, label = "\u03c9 = 1", hjust = -0.15, vjust = 1.6,
      colour = pal$subtitle, size = 2.8, fontface = "italic") +
    ggplot2::geom_point(ggplot2::aes(size = .data$weight_positive),
                         shape = 21, stroke = 0.3, alpha = 0.85) +
    ggplot2::scale_fill_manual(
      values = c(ns = pal$neutral, selected = pal$red),
      labels = c(ns = sprintf("q > %g", q_threshold),
                  selected = sprintf("q \u2264 %g", q_threshold)),
      name = NULL, guide = ggplot2::guide_legend(override.aes = list(size = 3.5))) +
    ggplot2::scale_colour_manual(
      values = c(ns = pal$neutral, selected = pal$red_strong),
      guide = "none") +
    ggplot2::scale_size_continuous(range = c(1.3, 5.2),
                                    name = "positive\nweight") +
    ggplot2::labs(
      x = expression(log[10] ~ omega[positive]),
      y = expression(-log[10] ~ italic(q)[BH]),
      title = "BUSTED[S] episodic-selection volcano",
      subtitle = sprintf("%d OGs tested \u2022 %d significant at q \u2264 %g",
                          nrow(d), sum(d$selected), q_threshold),
      caption = "HyPhy v2.5 \u2022 BUSTED with synonymous-rate variation (BUSTED[S])"
    ) +
    .dnmb_anvio_theme(base_size = 12)

  if (label_top_n > 0L &&
      requireNamespace("ggrepel", quietly = TRUE) &&
      any(d$selected)) {
    top <- d[order(d$q_value), , drop = FALSE]
    top <- utils::head(top[top$selected, , drop = FALSE], label_top_n)
    if (nrow(top))
      p <- p + ggrepel::geom_text_repel(
        data = top,
        ggplot2::aes(label = sprintf("OG_%07d", .data$cluster_id)),
        size = 2.8, fontface = "italic", colour = pal$text,
        segment.colour = pal$rule, segment.size = 0.25,
        min.segment.length = 0.1, box.padding = 0.45,
        point.padding = 0.25, force = 2.5,
        nudge_y = 0.25, direction = "both",
        max.overlaps = 30, seed = 7L)
  }
  .maybe_write(p, output_file, width = 9.5, height = 7)
  p
}


#' @rdname plot_selection
#' @description `plot_fel_sites()` — per-site α (synonymous) vs β
#'   (non-synonymous) scatter for one OG's FEL run. Points above the
#'   diagonal with p <= `alpha` are sites under pervasive positive
#'   selection; below the diagonal with p <= `alpha` are negative
#'   (purifying) sites. The diagonal shows α=β (neutrality).
#' @param fel_json Path to a HyPhy FEL `.json` result file. Alternately,
#'   a parsed list already read with [jsonlite::fromJSON].
#' @param alpha Site-level significance cutoff. Default 0.1 (HyPhy default).
#' @param title Optional title override.
#' @param output_file Optional PDF/PNG path (300 dpi, 7.5x7 in).
#' @return A ggplot object.
#' @export
plot_fel_sites <- function(fel_json, alpha = 0.1,
                            title = NULL, output_file = NULL) {
  pal <- .dnmb_anvio_pal
  if (is.character(fel_json) && !file.exists(fel_json))
    return(.plot_empty("FEL JSON not found"))
  if (!requireNamespace("jsonlite", quietly = TRUE))
    return(.plot_empty("jsonlite package required"))
  j <- if (is.list(fel_json)) fel_json else
         jsonlite::fromJSON(fel_json, simplifyVector = TRUE)
  mle <- j[["MLE"]]
  if (is.null(mle) || is.null(mle$content))
    return(.plot_empty("FEL JSON has no MLE$content"))
  mat <- mle$content[["0"]]
  if (!is.matrix(mat)) mat <- do.call(rbind, mat)
  if (!length(mat)) return(.plot_empty("FEL JSON MLE matrix empty"))
  d <- tibble::tibble(
    site    = seq_len(nrow(mat)),
    alpha_v = as.numeric(mat[, 1L]),
    beta_v  = as.numeric(mat[, 2L]),
    pval    = as.numeric(mat[, 5L])
  )
  d$class <- dplyr::case_when(
    d$pval <= alpha & d$beta_v > d$alpha_v ~ "positive",
    d$pval <= alpha & d$alpha_v > d$beta_v ~ "negative",
    TRUE                                    ~ "neutral"
  )
  d$class <- factor(d$class, levels = c("negative", "neutral", "positive"))
  d$neg_log_p <- -log10(pmax(d$pval, .Machine$double.eps))

  lim <- max(c(d$alpha_v, d$beta_v, 1), na.rm = TRUE)
  n_pos <- sum(d$class == "positive", na.rm = TRUE)
  n_neg <- sum(d$class == "negative", na.rm = TRUE)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$alpha_v, y = .data$beta_v)) +
    # Light tinted regions for purifying (below diagonal) / adaptive
    # (above diagonal).
    ggplot2::annotate("polygon",
      x = c(0, lim, lim, 0), y = c(0, lim, 0, 0),
      fill = pal$navy, alpha = 0.03) +
    ggplot2::annotate("polygon",
      x = c(0, lim, 0, 0), y = c(0, lim, lim, 0),
      fill = pal$red, alpha = 0.04) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2,
                          colour = pal$rule, linewidth = 0.35) +
    ggplot2::annotate("text",
      x = lim * 0.08, y = lim * 0.92,
      label = "adaptive   (\u03b2 > \u03b1)", hjust = 0,
      colour = pal$red, size = 3.1, fontface = "italic") +
    ggplot2::annotate("text",
      x = lim * 0.92, y = lim * 0.08,
      label = "purifying   (\u03b1 > \u03b2)", hjust = 1,
      colour = pal$navy, size = 3.1, fontface = "italic") +
    ggplot2::geom_point(ggplot2::aes(fill = .data$class,
                                      size = .data$neg_log_p),
                         shape = 21, colour = "#1F2E4A33",
                         stroke = 0.25, alpha = 0.85) +
    ggplot2::scale_fill_manual(
      values = c(negative = pal$navy,
                  neutral  = pal$neutral,
                  positive = pal$red),
      labels = c(negative = sprintf("purifying  (n=%d)", n_neg),
                  neutral  = "neutral",
                  positive = sprintf("adaptive   (n=%d)", n_pos)),
      name = NULL,
      guide = ggplot2::guide_legend(override.aes = list(size = 3.5))) +
    ggplot2::scale_size_continuous(range = c(0.9, 4.5),
                                    name = expression(-log[10]~italic(p))) +
    ggplot2::coord_fixed(xlim = c(0, lim), ylim = c(0, lim)) +
    ggplot2::labs(
      x = expression(alpha ~ "  (synonymous rate)"),
      y = expression(beta  ~ "  (non-synonymous rate)"),
      title = title %||% "HyPhy FEL per-site rates",
      subtitle = sprintf(
        "%d codon sites \u2022 %d adaptive / %d purifying at p \u2264 %g",
        nrow(d), n_pos, n_neg, alpha),
      caption = "Diagonal = neutral expectation (\u03b1 = \u03b2) \u2022 above \u2192 adaptive"
    ) +
    .dnmb_anvio_theme(base_size = 12)
  .maybe_write(p, output_file, width = 7.5, height = 7)
  p
}


#' @rdname plot_selection
#' @description `plot_absrel_tree()` — gene tree with branches coloured by
#'   aBSREL corrected p-value. Selected branches (p <= 0.05) are drawn
#'   thick red; everything else in muted navy. Uses base graphics via
#'   [ape::plot.phylo] because `ggtree` is a Bioconductor Suggest.
#' @param tree An `ape::phylo` gene tree.
#' @param absrel_json Path to the aBSREL JSON for the same tree.
#' @param output_file Optional PDF path (300 dpi, 8x6 in).
#' @return Invisibly, the `phylo` object. This renders on the active
#'   device (or a PDF if `output_file` is set).
#' @export
plot_absrel_tree <- function(tree, absrel_json, output_file = NULL) {
  pal <- .dnmb_anvio_pal
  stopifnot(inherits(tree, "phylo"))
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning("[plot_absrel_tree] jsonlite required; skipping.")
    return(invisible(tree))
  }
  j <- jsonlite::fromJSON(absrel_json, simplifyVector = FALSE)
  ba <- j[["branch attributes"]][[1L]]
  if (is.null(ba)) {
    warning("[plot_absrel_tree] no branch attributes in JSON.")
    return(invisible(tree))
  }
  lab <- c(tree$tip.label,
           if (!is.null(tree$node.label)) tree$node.label else
             character(tree$Nnode))
  pvec <- vapply(seq_along(lab), function(i) {
    nm <- lab[i]
    v  <- ba[[nm]][["Corrected P-value"]]
    if (is.null(v)) NA_real_ else as.numeric(v)
  }, numeric(1))
  edge_p   <- pvec[tree$edge[, 2L]]
  sel      <- !is.na(edge_p) & edge_p <= 0.05
  edge_col <- ifelse(sel, pal$red_strong, pal$navy)
  edge_lwd <- ifelse(sel, 2.6, 1.0)

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(output_file, width = 8, height = 6, bg = pal$bg)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  op <- graphics::par(mar = c(1, 1, 2.2, 1), bg = pal$bg,
                       col.main = pal$text, family = "")
  on.exit(graphics::par(op), add = TRUE)
  ape::plot.phylo(tree, edge.color = edge_col, edge.width = edge_lwd,
                  cex = 0.75, font = 3, tip.color = pal$text,
                  no.margin = FALSE)
  graphics::title(
    main = sprintf("aBSREL \u2022 %d / %d branches under selection (p \u2264 0.05)",
                    sum(sel, na.rm = TRUE), length(edge_p)),
    col.main = pal$text, cex.main = 1.0, font.main = 2)
  invisible(tree)
}


#' @rdname plot_selection
#' @description `plot_codeml_vs_hyphy()` — scatter of codeml M0 ω vs
#'   HyPhy BUSTED positive-class ω across OGs. Points on the diagonal
#'   agree; OGs far above the diagonal are flagged by BUSTED as
#'   episodically selected but hidden by M0's one-ratio average.
#' @param codeml_tbl Tibble from [run_codeml_m0()] (needs `cluster_id`,
#'   `omega`).
#' @param busted_tbl Tibble from [run_hyphy_busted()] (needs
#'   `cluster_id`, `omega_positive`).
#' @param output_file Optional PDF/PNG path (300 dpi, 7.5x7 in).
#' @return A ggplot.
#' @export
plot_codeml_vs_hyphy <- function(codeml_tbl, busted_tbl, output_file = NULL) {
  pal <- .dnmb_anvio_pal
  if (is.null(codeml_tbl) || is.null(busted_tbl))
    return(.plot_empty("codeml or BUSTED input is NULL"))
  cml <- codeml_tbl[, c("cluster_id", "omega"), drop = FALSE]
  bst <- busted_tbl[, c("cluster_id", "omega_positive", "p_value"),
                    drop = FALSE]
  d <- dplyr::inner_join(cml, bst, by = "cluster_id")
  d <- d[!is.na(d$omega) & !is.na(d$omega_positive), , drop = FALSE]
  if (!nrow(d)) return(.plot_empty("no OGs with both codeml and BUSTED \u03c9"))
  d$sig   <- !is.na(d$p_value) & d$p_value <= 0.05
  d$class <- factor(ifelse(d$sig, "episodic", "concordant"),
                     levels = c("concordant", "episodic"))

  # Zoom to the active region instead of naive coord_fixed(): use a
  # robust 95th percentile + 10% headroom so a handful of large outliers
  # don't push the cloud into the corner.
  q_hi <- max(
    stats::quantile(c(d$omega, d$omega_positive), 0.95, na.rm = TRUE),
    1.2, na.rm = TRUE
  )
  lim <- min(
    max(c(d$omega, d$omega_positive, 1), na.rm = TRUE),
    q_hi * 1.6
  )
  lim <- max(lim, 1.5)
  # Clip to lim but leave a small safety margin so the point glyph
  # radius doesn't get cut at the axis edge.
  clip_to <- lim * 0.985
  d_clip <- d
  d_clip$omega          <- pmin(d_clip$omega, clip_to)
  d_clip$omega_positive <- pmin(d_clip$omega_positive, clip_to)
  p <- ggplot2::ggplot(d_clip, ggplot2::aes(x = .data$omega,
                                             y = .data$omega_positive)) +
    ggplot2::annotate("rect",
      xmin = 0, xmax = 1, ymin = 1, ymax = lim,
      fill = pal$red, alpha = 0.06) +
    ggplot2::annotate("text",
      x = 0.5, y = lim * 0.95,
      label = "BUSTED-only\nepisodic zone",
      hjust = 0.5, vjust = 1,
      colour = pal$red_strong, size = 2.9, fontface = "italic",
      lineheight = 0.9) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2,
                          colour = pal$rule, linewidth = 0.35) +
    ggplot2::annotate("text",
      x = lim * 0.82, y = lim * 0.78,
      label = "concordant", angle = 45,
      hjust = 0.5, vjust = -0.5,
      colour = pal$subtitle, size = 2.9, fontface = "italic") +
    ggplot2::geom_vline(xintercept = 1, colour = pal$rule,
                         linetype = 3, linewidth = 0.25) +
    ggplot2::geom_hline(yintercept = 1, colour = pal$rule,
                         linetype = 3, linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$class),
                         shape = 21, size = 2.6, stroke = 0.3,
                         colour = "#1F2E4A66", alpha = 0.9) +
    ggplot2::scale_fill_manual(
      values = c(concordant = pal$navy, episodic = pal$red),
      labels = c(concordant = "concordant",
                  episodic   = sprintf("BUSTED episodic (p \u2264 0.05, n=%d)",
                                        sum(d$sig))),
      name = NULL,
      guide = ggplot2::guide_legend(override.aes = list(size = 3.5))) +
    ggplot2::coord_cartesian(xlim = c(-lim * 0.02, lim),
                              ylim = c(-lim * 0.02, lim),
                              expand = FALSE) +
    ggplot2::labs(
      x = expression("codeml M0 " ~ omega),
      y = expression("BUSTED positive-class " ~ omega),
      title = "codeml M0 vs BUSTED \u03c9 per OG",
      subtitle = sprintf(
        "%d OGs \u2022 diagonal = agreement \u2022 upper-left = BUSTED-only signal",
        nrow(d)),
      caption = "Wolf 2021: pervasive M0 \u03c9 can mask branch-specific episodic selection"
    ) +
    .dnmb_anvio_theme(base_size = 12)
  .maybe_write(p, output_file, width = 7.5, height = 7)
  p
}


#' @rdname plot_selection
#' @description `plot_selection_overview()` — four-panel publication
#'   figure: (1) BUSTED volcano, (2) FEL positive/negative site counts
#'   per OG, (3) BUSTED q-value histogram with BH cutoff, (4) FEL ω
#'   median distribution. Requires `patchwork`.
#' @param selection Named list with `busted`, `fel`, `absrel`, `combined`.
#' @param output_file Optional PDF/PNG path (300 dpi, 14x10 in).
#' @return A patchwork / ggplot object.
#' @export
plot_selection_overview <- function(selection, output_file = NULL) {
  pal <- .dnmb_anvio_pal
  if (is.null(selection))
    return(.plot_empty("No selection results"))

  panel_v <- plot_busted_volcano(selection$busted, label_top_n = 4L)
  panel_b <- .panel_fel_counts(selection$fel)
  panel_q <- .panel_busted_qhist(selection$busted)
  panel_o <- .panel_fel_omega(selection$fel)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(
      panel_v, panel_q,
      panel_b, panel_o,
      nrow = 2, ncol = 2,
      heights = c(1, 0.95)
    ) +
      patchwork::plot_annotation(
        title = "HyPhy selection summary",
        subtitle = "BUSTED[S] + FEL \u2022 anvi'o-style DNMB display",
        caption = "DNMBcluster \u2022 HyPhy v2.5 \u2022 BH-adjusted q-values",
        tag_levels = "A",
        theme = ggplot2::theme(
          plot.background   = ggplot2::element_rect(fill = pal$bg, colour = NA),
          plot.title        = ggplot2::element_text(size = 16, face = "bold",
                                                     colour = pal$text),
          plot.subtitle     = ggplot2::element_text(size = 10,
                                                     colour = pal$subtitle),
          plot.caption      = ggplot2::element_text(size = 8,
                                                     colour = pal$subtitle)
        )
      ) &
      ggplot2::theme(
        plot.tag = ggplot2::element_text(face = "bold", size = 14,
                                          colour = pal$text),
        plot.tag.position = c(0.01, 0.99))
  } else {
    combined <- panel_v
  }
  .maybe_write(combined, output_file, width = 13, height = 11)
  combined
}


#' @rdname plot_selection
#' @description `plot_selection_circular()` — anvi'o-style circular
#'   summary: the gene / OG tree in the centre, with concentric data
#'   rings for each HyPhy test. Ring inside-out:
#'   (1) BUSTED significance dot (navy safe / orange marginal / red q<=0.05),
#'   (2) BUSTED positive-class ω tile (YlOrRd gradient),
#'   (3) FEL positive-site count dot (size-scaled navy),
#'   (4) aBSREL min corrected p-value tile (YlGnBu reversed).
#'   If `tree` is omitted, OGs are arranged radially by `cluster_id`.
#' @param selection Named list with `busted`, `fel`, `absrel`.
#' @param tree Optional `ape::phylo` OG tree (tip labels = cluster_id or
#'   `OG_xxxxxxx`). When NULL, OGs fan around an invisible trunk.
#' @param q_threshold BH cutoff used by the inner ring. Default 0.05.
#' @param output_file Optional PDF/PNG path (300 dpi, 12x12 in).
#' @return A ggplot object, or an empty placeholder if `ggtree` /
#'   `ggnewscale` are unavailable.
#' @export
plot_selection_circular <- function(selection, tree = NULL,
                                     q_threshold = 0.05, output_file = NULL) {
  pal <- .dnmb_anvio_pal
  if (is.null(selection)) return(.plot_empty("No selection results"))
  if (!requireNamespace("ggtree", quietly = TRUE) ||
      !requireNamespace("ggnewscale", quietly = TRUE)) {
    warning("[plot_selection_circular] ggtree + ggnewscale required; ",
            "falling back to plot_selection_overview().")
    return(plot_selection_overview(selection, output_file = output_file))
  }

  bst <- selection$busted
  fel <- selection$fel
  abs <- selection$absrel
  if (is.null(bst) || !nrow(bst))
    return(.plot_empty("No BUSTED results; cannot render circular view"))

  # Build per-OG payload ------------------------------------------------
  bst <- bst
  if (!"q_value" %in% names(bst))
    bst$q_value <- stats::p.adjust(bst$p_value, method = "BH")
  payload <- tibble::tibble(
    cluster_id     = bst$cluster_id,
    label          = sprintf("OG_%07d", bst$cluster_id),
    busted_q       = bst$q_value,
    omega_positive = bst$omega_positive,
    weight_positive = bst$weight_positive
  )
  payload$busted_class <- dplyr::case_when(
    is.na(payload$busted_q)                      ~ "missing",
    payload$busted_q <= q_threshold              ~ "selected",
    payload$busted_q <= min(0.25, q_threshold*5) ~ "marginal",
    TRUE                                          ~ "ns"
  )
  payload$busted_class <- factor(payload$busted_class,
    levels = c("selected", "marginal", "ns", "missing"))

  if (!is.null(fel) && nrow(fel)) {
    fel_sub <- fel[, c("cluster_id", "n_positive", "n_negative",
                        "median_omega"), drop = FALSE]
    names(fel_sub)[-1] <- paste0("fel_", names(fel_sub)[-1])
    payload <- dplyr::left_join(payload, fel_sub, by = "cluster_id")
  } else {
    payload$fel_n_positive  <- NA_integer_
    payload$fel_n_negative  <- NA_integer_
    payload$fel_median_omega <- NA_real_
  }
  if (!is.null(abs) && nrow(abs)) {
    abs_sub <- abs[, c("cluster_id", "n_selected", "min_corrected_p"),
                   drop = FALSE]
    names(abs_sub)[-1] <- paste0("absrel_", names(abs_sub)[-1])
    payload <- dplyr::left_join(payload, abs_sub, by = "cluster_id")
  } else {
    payload$absrel_n_selected      <- NA_integer_
    payload$absrel_min_corrected_p <- NA_real_
  }

  # Assemble the center tree (or synthetic fan) -------------------------
  phylo <- .ensure_circular_tree(tree, payload)
  meta  <- payload[match(phylo$tip.label, payload$label), , drop = FALSE]
  # ggtree's %<+% treats the first column as the join key and renames it
  # to "label", so reorder and drop any leftover duplicate.
  meta$.tip <- phylo$tip.label
  meta <- meta[, c(".tip",
                    setdiff(colnames(meta), c(".tip", "label"))),
               drop = FALSE]

  max_x <- max(ggtree::fortify(phylo)$x)
  r1 <- max_x * 1.20   # BUSTED class dot
  r2 <- max_x * 1.35   # BUSTED omega gradient
  r3 <- max_x * 1.55   # FEL positive sites
  r4 <- max_x * 1.75   # aBSREL min corrected p

  p <- ggtree::ggtree(phylo, layout = "circular", size = 0.45,
                       ladderize = TRUE, colour = pal$text) +
    ggplot2::xlim(-0.2, max_x * 2.35)
  p <- ggtree::`%<+%`(p, as.data.frame(meta))

  # Ring 1: BUSTED class dot
  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(x = r1, colour = .data$busted_class,
                    size   = .data$weight_positive),
      shape = 16, alpha = 0.9) +
    ggplot2::scale_colour_manual(
      name  = "BUSTED class",
      values = c(selected = pal$red_strong, marginal = pal$orange,
                  ns = pal$navy, missing = pal$neutral),
      drop = FALSE) +
    ggplot2::scale_size_continuous(name = "pos weight",
                                    range = c(1.0, 3.5),
                                    guide = ggplot2::guide_legend(order = 2))

  # Ring 2: BUSTED omega tile
  p <- p + ggnewscale::new_scale_color() +
    ggplot2::geom_point(
      ggplot2::aes(x = r2,
                    color = log10(pmax(.data$omega_positive, 1e-3))),
      shape = 15, size = 2.9) +
    ggplot2::scale_color_distiller(
      name = expression(log[10] ~ omega[pos]),
      palette = "YlOrRd", direction = 1, na.value = "grey85")

  # Ring 3: FEL positive-site dot
  if (any(!is.na(meta$fel_n_positive))) {
    p <- p + ggnewscale::new_scale_color() +
      ggplot2::geom_point(
        ggplot2::aes(x = r3, size = .data$fel_n_positive),
        colour = pal$navy, alpha = 0.85, shape = 16) +
      ggplot2::scale_size_continuous(
        name = "FEL +sites", range = c(0.8, 4.5),
        guide = ggplot2::guide_legend(order = 3))
  }

  # Ring 4: aBSREL min p
  if (any(!is.na(meta$absrel_min_corrected_p))) {
    p <- p + ggnewscale::new_scale_color() +
      ggplot2::geom_point(
        ggplot2::aes(x = r4,
                      color = -log10(pmax(.data$absrel_min_corrected_p,
                                           .Machine$double.eps))),
        shape = 15, size = 2.9) +
      ggplot2::scale_color_distiller(
        name = expression(-log[10] ~ italic(p)[aBSREL]),
        palette = "YlGnBu", direction = 1, na.value = "grey85")
  }

  # Ring-level annotation labels placed at 12 o'clock (y = 0 in tree
  # coords, which after circular layout corresponds to angle 0 / top).
  ring_ann <- data.frame(
    x = c(r1, r2, r3, r4),
    y = 0,
    label = c("BUSTED class",
              "log\u2081\u2080 \u03c9\u208a",
              "FEL +sites",
              "\u2212log\u2081\u2080 p aBSREL")
  )

  p <- p +
    ggtree::geom_tiplab2(
      ggplot2::aes(label = .data$label),
      fontface = "italic", size = 2.4, colour = pal$text,
      align = TRUE, linesize = 0.1, linetype = 3,
      offset = max_x * 0.85, show.legend = FALSE) +
    ggplot2::geom_text(
      data = ring_ann, inherit.aes = FALSE,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      nudge_y = -0.5,
      colour = pal$subtitle, size = 2.5,
      fontface = "italic") +
    ggplot2::labs(
      title = "HyPhy selection \u2022 circular summary",
      subtitle = paste0("Concentric rings per OG: BUSTED class \u2022 ",
                         "log\u2081\u2080 \u03c9\u208a \u2022 FEL +sites \u2022 ",
                         "\u2212log\u2081\u2080 p aBSREL"),
      caption = "DNMBcluster \u2022 HyPhy v2.5 \u2022 Benjamini\u2013Hochberg FDR") +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
      plot.title       = ggplot2::element_text(size = 14, face = "bold",
                                                 colour = pal$text),
      plot.subtitle    = ggplot2::element_text(size = 8,
                                                 colour = pal$subtitle),
      plot.caption     = ggplot2::element_text(size = 7,
                                                 colour = pal$subtitle),
      legend.position  = "right",
      legend.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
      legend.key       = ggplot2::element_rect(fill = pal$bg, colour = NA),
      legend.key.size  = ggplot2::unit(0.45, "cm"),
      legend.spacing.y = ggplot2::unit(0.35, "cm"),
      legend.text      = ggplot2::element_text(size = 8, colour = pal$text),
      legend.title     = ggplot2::element_text(size = 9, face = "bold",
                                                 colour = pal$text),
      plot.margin      = ggplot2::margin(14, 14, 14, 14))

  .maybe_write(p, output_file, width = 13, height = 12)
  p
}


# ----------------- internals ---------------------------------------------

.panel_fel_counts <- function(fel) {
  pal <- .dnmb_anvio_pal
  if (is.null(fel) || !nrow(fel)) return(.plot_empty("No FEL results"))
  d <- tidyr::pivot_longer(
    dplyr::mutate(fel, positive = n_positive, negative = n_negative),
    c("positive", "negative"),
    names_to = "class", values_to = "count")
  d <- d[!is.na(d$count) & d$count > 0, , drop = FALSE]
  if (!nrow(d)) return(.plot_empty("No FEL sites passed p \u2264 0.1"))
  ids <- sort(unique(d$cluster_id))
  d$og <- factor(sprintf("OG_%07d", d$cluster_id),
                  levels = sprintf("OG_%07d", ids))
  d$class <- factor(d$class, levels = c("negative", "positive"))

  # Thin x-axis labels when the set is too dense for readability.
  n_og <- length(ids)
  step <- if (n_og <= 12) 1L else if (n_og <= 25) 2L else
            max(1L, ceiling(n_og / 12))
  keep_levels <- levels(d$og)[seq(1L, length(levels(d$og)), by = step)]
  x_label_fun <- function(x) ifelse(x %in% keep_levels, x, "")

  ggplot2::ggplot(d, ggplot2::aes(x = .data$og, y = .data$count,
                                   fill = .data$class)) +
    ggplot2::geom_col(position = "stack", width = 0.8,
                       colour = "#1F2E4A22", linewidth = 0.15) +
    ggplot2::scale_x_discrete(labels = x_label_fun) +
    ggplot2::scale_fill_manual(
      values = c(negative = pal$navy, positive = pal$red),
      labels = c(negative = "purifying", positive = "adaptive"),
      name = NULL) +
    ggplot2::labs(
      x = NULL, y = "FEL site count",
      title = "FEL per-OG site tally",
      subtitle = sprintf(
        "%d OGs with \u22651 significant site%s",
        n_og,
        if (step > 1L) sprintf(" \u2022 every %d\u1d57\u02b0 label shown", step)
        else "")) +
    .dnmb_anvio_theme(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1,
                                            vjust = 1,
                                            face = "italic", size = 7))
}

.panel_busted_qhist <- function(bst) {
  pal <- .dnmb_anvio_pal
  if (is.null(bst) || !nrow(bst))
    return(.plot_empty("No BUSTED results"))
  d <- bst
  if (!"q_value" %in% names(d))
    d$q_value <- stats::p.adjust(d$p_value, method = "BH")
  d <- d[!is.na(d$q_value), , drop = FALSE]
  if (!nrow(d)) return(.plot_empty("No BUSTED q-values"))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$q_value)) +
    ggplot2::geom_histogram(
      fill = pal$navy, colour = pal$bg, linewidth = 0.3,
      bins = 25, alpha = 0.9) +
    ggplot2::geom_vline(xintercept = 0.05, colour = pal$red_strong,
                         linetype = 2, linewidth = 0.4) +
    ggplot2::annotate("text", x = 0.05, y = Inf,
                       label = "q = 0.05", hjust = -0.1, vjust = 1.6,
                       colour = pal$red_strong, size = 3,
                       fontface = "italic") +
    ggplot2::labs(
      x = expression(italic(q)[BH]), y = "OG count",
      title = "BUSTED q-value distribution",
      subtitle = sprintf("%d OGs \u2022 %d below q = 0.05",
                          nrow(d), sum(d$q_value <= 0.05))) +
    .dnmb_anvio_theme(base_size = 11)
}

.panel_fel_omega <- function(fel) {
  pal <- .dnmb_anvio_pal
  if (is.null(fel) || !nrow(fel)) return(.plot_empty("No FEL results"))
  d <- fel[!is.na(fel$median_omega), , drop = FALSE]
  if (!nrow(d)) return(.plot_empty("No FEL median omega"))
  d$class <- ifelse(!is.na(d$n_positive) & d$n_positive > 0,
                     "adaptive", "conserved")
  d$class <- factor(d$class, levels = c("conserved", "adaptive"))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$median_omega,
                                   fill = .data$class)) +
    ggplot2::geom_histogram(
      bins = 22, colour = pal$bg, linewidth = 0.3, alpha = 0.9,
      position = "identity") +
    ggplot2::geom_vline(xintercept = 1, colour = pal$rule,
                         linetype = 2, linewidth = 0.35) +
    ggplot2::scale_fill_manual(
      values = c(conserved = pal$navy, adaptive = pal$red),
      name = NULL) +
    ggplot2::scale_x_continuous(trans = "log10",
                                 labels = function(x)
                                   formatC(x, format = "fg", digits = 2)) +
    ggplot2::labs(
      x = expression("FEL median " ~ omega ~ "  (log scale)"),
      y = "OG count",
      title = "FEL median \u03c9 per OG",
      subtitle = sprintf(
        "%d OGs \u2022 %d with \u22651 adaptive site \u2022 dashed = neutrality",
        nrow(d), sum(d$class == "adaptive"))) +
    .dnmb_anvio_theme(base_size = 11)
}

.ensure_circular_tree <- function(tree, payload) {
  # Map / build a phylo with tip labels = payload$label.
  wanted <- payload$label
  if (inherits(tree, "phylo")) {
    renamed <- .match_tree_labels(tree, payload)
    keep <- intersect(renamed$tip.label, wanted)
    if (length(keep) >= 2L) {
      drop <- setdiff(renamed$tip.label, keep)
      if (length(drop)) renamed <- ape::drop.tip(renamed, drop)
      if (length(renamed$tip.label) >= 2L) return(renamed)
    }
  }
  # Fallback: synthetic star-tree fan ordered by cluster_id.
  ord <- order(payload$cluster_id)
  labs <- wanted[ord]
  nw <- paste0("(", paste(labs, collapse = ","), ");")
  ape::read.tree(text = nw)
}

.match_tree_labels <- function(tree, payload) {
  # If tree tips already look like OG_xxxxxxx leave them; otherwise try
  # to coerce integer cluster_ids to that form.
  if (all(grepl("^OG_", tree$tip.label))) return(tree)
  as_int <- suppressWarnings(as.integer(tree$tip.label))
  if (all(!is.na(as_int))) {
    tree$tip.label <- sprintf("OG_%07d", as_int)
  }
  tree
}

.plot_empty <- function(msg) {
  pal <- .dnmb_anvio_pal
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg,
                       size = 4.5, colour = pal$subtitle,
                       fontface = "italic") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA))
}

.maybe_write <- function(p, output_file, width = 9, height = 6.5) {
  if (is.null(output_file)) return(invisible(NULL))
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(output_file, p, width = width, height = height,
                  dpi = 300, bg = .dnmb_anvio_pal$bg)
  invisible(output_file)
}

# null-coalescing operator (DNMBcluster doesn't export its own).
`%||%` <- function(a, b) if (is.null(a)) b else a
