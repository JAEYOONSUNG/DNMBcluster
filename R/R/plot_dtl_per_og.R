#' Per-OG DTL figures — one PDF per question
#'
#' `plot_duplication_burden()` focuses on duplications (D). This family
#' of plots widens the lens to cover the four reconciliation event
#' classes that `reconcile_dtl()` writes to `dtl_per_og.tsv`:
#'
#' \itemize{
#'   \item **S** — speciation events along the gene tree
#'   \item **D** — duplications (already surfaced by `plot_duplication_burden`)
#'   \item **T** — horizontal transfer candidates
#'   \item **loss** — gene losses along the species tree
#' }
#'
#' Each entry point writes its own PDF so reviewers can pick the
#' representative figure(s) to merge later:
#'
#' \describe{
#'   \item{`plot_dtl_event_histogram()`}{log-scale faceted histograms
#'     per event type, with per-facet mean / median / max stats.}
#'   \item{`plot_dtl_dt_coupling()`}{D \eqn{\times} T hex-bin density
#'     scatter with a 1:1 reference and the top-2% OGs overlaid.}
#'   \item{`plot_dtl_concentration()`}{Lorenz-style cumulative curves
#'     (one per event type) with Gini-like concentration metric.}
#'   \item{`plot_dtl_top_ogs()`}{horizontal-bar triptych of top-N OGs
#'     by D, T and loss.}
#'   \item{`plot_dtl_per_og_summary()`}{convenience wrapper — calls all
#'     four so they render in one go.}
#' }
#'
#' All figures use the same `.dnmb_anvio_pal` look as
#' [plot_busted_volcano] so the selection / reconciliation panels read
#' as one family in a write-up.
#'
#' @param out_dir `run_orthofinder_like()` directory containing
#'   `dtl_per_og.tsv` (produced when `dtl = TRUE`).
#' @param out_pdf Destination PDF. If `NULL`, a sensible default is
#'   written next to the TSV.
#' @param top_n Number of OGs to show on the bar panels (default 15).
#' @param verbose Echo progress.
#' @return Path to the written PDF (invisibly), or `NULL` on failure.
#' @name plot_dtl_per_og
NULL

# ---- shared loader -----------------------------------------------------

.dtl_per_og_load <- function(out_dir) {
  po_tsv <- file.path(out_dir, "dtl_per_og.tsv")
  if (!file.exists(po_tsv)) {
    warning("[dtl_per_og] dtl_per_og.tsv not found -- run with dtl = TRUE")
    return(NULL)
  }
  df <- utils::read.table(po_tsv, sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE)
  if (!nrow(df)) {
    warning("[dtl_per_og] dtl_per_og.tsv is empty")
    return(NULL)
  }
  df
}

.dtl_need <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1L),
                           quietly = TRUE)]
  if (length(missing)) {
    warning(sprintf("[dtl_per_og] required package(s) missing: %s",
                    paste(missing, collapse = ", ")))
    return(FALSE)
  }
  TRUE
}

.dtl_ev_cols <- function(pal) {
  c(S = pal$navy, D = pal$red, T = pal$orange, loss = pal$blue_cool)
}

# ---- event histogram ---------------------------------------------------

#' @rdname plot_dtl_per_og
#' @export
plot_dtl_event_histogram <- function(out_dir,
                                      out_pdf = NULL,
                                      verbose = TRUE) {
  if (!.dtl_need("ggplot2")) return(NULL)
  df <- .dtl_per_og_load(out_dir); if (is.null(df)) return(NULL)
  pal <- .dnmb_anvio_pal
  ev_cols <- .dtl_ev_cols(pal)

  p <- .panel_dtl_hist(df, pal, ev_cols) +
    ggplot2::labs(
      title    = "Per-OG DTL event distribution",
      subtitle = sprintf(
        "%d OGs \u2022 totals S=%d D=%d T=%d loss=%d",
        nrow(df), sum(df$n_S), sum(df$n_D),
        sum(df$n_T), sum(df$n_loss))
    )
  if (is.null(out_pdf))
    out_pdf <- file.path(out_dir, "dtl_event_histogram.pdf")
  ggplot2::ggsave(out_pdf, p, width = 11, height = 4.6,
                   dpi = 300, bg = pal$bg, limitsize = FALSE)
  if (verbose) message(sprintf("[dtl_event_histogram] wrote %s", out_pdf))
  invisible(out_pdf)
}

# ---- D vs T scatter ----------------------------------------------------

#' @rdname plot_dtl_per_og
#' @export
plot_dtl_dt_coupling <- function(out_dir,
                                  out_pdf = NULL,
                                  verbose = TRUE) {
  if (!.dtl_need("ggplot2")) return(NULL)
  df <- .dtl_per_og_load(out_dir); if (is.null(df)) return(NULL)
  pal <- .dnmb_anvio_pal
  p <- .panel_dtl_couple(df, pal)
  if (is.null(out_pdf))
    out_pdf <- file.path(out_dir, "dtl_dt_coupling.pdf")
  ggplot2::ggsave(out_pdf, p, width = 8.5, height = 7.5,
                   dpi = 300, bg = pal$bg, limitsize = FALSE)
  if (verbose) message(sprintf("[dtl_dt_coupling] wrote %s", out_pdf))
  invisible(out_pdf)
}

# ---- Lorenz / concentration -------------------------------------------

#' @rdname plot_dtl_per_og
#' @export
plot_dtl_concentration <- function(out_dir,
                                    out_pdf = NULL,
                                    verbose = TRUE) {
  if (!.dtl_need("ggplot2")) return(NULL)
  df <- .dtl_per_og_load(out_dir); if (is.null(df)) return(NULL)
  pal <- .dnmb_anvio_pal
  ev_cols <- .dtl_ev_cols(pal)
  p <- .panel_dtl_lorenz(df, pal, ev_cols)
  if (is.null(out_pdf))
    out_pdf <- file.path(out_dir, "dtl_concentration.pdf")
  ggplot2::ggsave(out_pdf, p, width = 7.5, height = 7.5,
                   dpi = 300, bg = pal$bg, limitsize = FALSE)
  if (verbose) message(sprintf("[dtl_concentration] wrote %s", out_pdf))
  invisible(out_pdf)
}

# ---- top-OG bar triptych ----------------------------------------------

#' @rdname plot_dtl_per_og
#' @export
plot_dtl_top_ogs <- function(out_dir,
                              out_pdf = NULL,
                              top_n = 15L,
                              verbose = TRUE) {
  if (!.dtl_need(c("ggplot2", "patchwork"))) return(NULL)
  df <- .dtl_per_og_load(out_dir); if (is.null(df)) return(NULL)
  pal <- .dnmb_anvio_pal
  ev_cols <- .dtl_ev_cols(pal)

  # If there is literally nothing to rank, render a single empty-state card
  # instead of three giant "No OGs..." panels — much less visually loud.
  if (sum(df$n_D, df$n_T, df$n_loss, na.rm = TRUE) == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0.65,
                         label = "No duplication, transfer, or loss events",
                         size = 5.2, colour = pal$text, fontface = "bold") +
      ggplot2::annotate("text", x = 0, y = 0.30,
                         label = sprintf("All %d OGs in this run reconcile to S only.",
                                          nrow(df)),
                         size = 3.8, colour = pal$subtitle,
                         fontface = "italic") +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::theme_void() +
      .dnmb_panel_card(pal, c(40, 40, 40, 40))
    if (is.null(out_pdf))
      out_pdf <- file.path(out_dir, "dtl_top_ogs.pdf")
    ggplot2::ggsave(out_pdf, p, width = 11, height = 4.2,
                     dpi = 300, bg = pal$bg, limitsize = FALSE)
    if (verbose) message(sprintf("[dtl_top_ogs] empty-state -> %s", out_pdf))
    return(invisible(out_pdf))
  }

  p <- .panel_dtl_top(df, pal, ev_cols, top_n) +
    patchwork::plot_annotation(
      title = sprintf("Top %d OGs by duplication / transfer / loss", top_n),
      subtitle = sprintf(
        "%d OGs in pool \u2022 totals D=%d T=%d loss=%d",
        nrow(df), sum(df$n_D), sum(df$n_T), sum(df$n_loss)),
      caption = sprintf(
        "bars sorted within event class \u2022 colours match other DTL panels\n%s",
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
  if (is.null(out_pdf))
    out_pdf <- file.path(out_dir, "dtl_top_ogs.pdf")
  ggplot2::ggsave(out_pdf, p, width = 13, height = 5.5,
                   dpi = 300, bg = pal$bg, limitsize = FALSE)
  if (verbose) message(sprintf("[dtl_top_ogs] wrote %s", out_pdf))
  invisible(out_pdf)
}

# ---- convenience: all four --------------------------------------------

#' @rdname plot_dtl_per_og
#' @export
plot_dtl_per_og_summary <- function(out_dir,
                                     top_n = 15L,
                                     verbose = TRUE) {
  out <- list(
    histogram     = plot_dtl_event_histogram(out_dir, verbose = verbose),
    dt_coupling   = plot_dtl_dt_coupling(out_dir,    verbose = verbose),
    concentration = plot_dtl_concentration(out_dir,  verbose = verbose),
    top_ogs       = plot_dtl_top_ogs(out_dir, top_n = top_n,
                                      verbose = verbose)
  )
  invisible(out)
}

# ---- shared panel builders --------------------------------------------

.panel_dtl_hist <- function(df, pal, ev_cols) {
  long <- data.frame(
    cluster_id = rep(df$cluster_id, 4L),
    event = rep(c("S", "D", "T", "loss"), each = nrow(df)),
    count = c(df$n_S, df$n_D, df$n_T, df$n_loss)
  )
  long$event <- factor(long$event, levels = c("S", "D", "T", "loss"))
  long <- long[!is.na(long$count), , drop = FALSE]

  stats_df <- do.call(rbind, lapply(levels(long$event), function(ev) {
    v <- long$count[long$event == ev]
    data.frame(
      event  = ev,
      mean   = mean(v, na.rm = TRUE),
      median = stats::median(v, na.rm = TRUE),
      max    = max(v, na.rm = TRUE),
      n_nz   = sum(v > 0, na.rm = TRUE)
    )
  }))
  stats_df$event <- factor(stats_df$event, levels = levels(long$event))
  stats_df$label <- sprintf("mean %.2f \u2022 median %g\nmax %g \u2022 n\u2260 0 = %d",
                             stats_df$mean, stats_df$median,
                             stats_df$max,  stats_df$n_nz)

  ggplot2::ggplot(long, ggplot2::aes(x = .data$count + 1L,
                                      fill = .data$event)) +
    ggplot2::geom_histogram(bins = 28, colour = pal$bg,
                             linewidth = 0.2, alpha = 0.9) +
    ggplot2::scale_x_log10(
      breaks = c(1, 2, 5, 10, 25, 50, 100, 250, 500),
      labels = function(x) x - 1L) +
    ggplot2::scale_y_continuous(
      trans = "pseudo_log",
      breaks = c(0, 1, 10, 100, 1000, 10000),
      labels = function(x) formatC(x, format = "d", big.mark = ",")
    ) +
    ggplot2::scale_fill_manual(values = ev_cols, guide = "none") +
    ggplot2::facet_wrap(~ event, nrow = 1L, scales = "free_y") +
    ggplot2::geom_text(
      data = stats_df, inherit.aes = FALSE,
      ggplot2::aes(x = Inf, y = Inf, label = .data$label),
      hjust = 1.05, vjust = 1.4, size = 2.7,
      colour = pal$subtitle, fontface = "italic", lineheight = 0.95
    ) +
    ggplot2::labs(
      x = "events per OG (log, +1 offset)", y = "OG count"
    ) +
    .dnmb_anvio_theme(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 9.5, colour = pal$subtitle),
      strip.text = ggplot2::element_text(face = "bold", size = 10,
                                           colour = pal$text),
      panel.spacing.x = ggplot2::unit(0.8, "lines")
    )
}

.panel_dtl_couple <- function(df, pal) {
  d <- df
  d$total <- d$n_D + d$n_T + d$n_loss
  d$highlight <- d$total >= stats::quantile(d$total, 0.98,
                                             na.rm = TRUE, names = FALSE) &
                 d$total > 0
  use_hex <- max(d$n_D, 0, na.rm = TRUE) > 0 &&
             max(d$n_T, 0, na.rm = TRUE) > 0 &&
             requireNamespace("hexbin", quietly = TRUE)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$n_D + 1L,
                                         y = .data$n_T + 1L)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2,
                          colour = pal$rule, linewidth = 0.3)
  if (use_hex) {
    p <- p +
      ggplot2::geom_hex(bins = 30, alpha = 0.9) +
      ggplot2::scale_fill_gradientn(
        name = "OG density",
        colours = c(pal$bg, "#CCD8DE", pal$navy, pal$red_strong),
        values  = scales::rescale(c(0, 0.1, 0.5, 1)),
        trans   = "log",
        breaks  = c(1, 5, 20, 100, 500),
        guide   = ggplot2::guide_colorbar(barwidth = 0.4, barheight = 6)
      )
  } else {
    p <- p + ggplot2::geom_jitter(width = 0.05, height = 0.05,
                                    colour = pal$navy, alpha = 0.35,
                                    size = 0.9)
  }
  p <- p +
    ggplot2::geom_point(
      data = d[d$highlight, , drop = FALSE],
      shape = 21, fill = pal$red, colour = "#1F2E4A", stroke = 0.35,
      size = 2.3, alpha = 0.95
    )

  top <- d[d$highlight, , drop = FALSE]
  top <- top[order(-top$total), , drop = FALSE]
  top <- utils::head(top, 8L)
  if (nrow(top) && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = top,
      ggplot2::aes(label = sprintf("OG_%07d", .data$cluster_id)),
      size = 2.6, fontface = "italic", colour = pal$text,
      segment.colour = pal$rule, segment.size = 0.25,
      min.segment.length = 0.1, box.padding = 0.4, max.overlaps = 20,
      seed = 11L
    )
  }

  p +
    ggplot2::scale_x_log10(breaks = c(1, 2, 5, 10, 25, 100, 500),
                            labels = function(x) x - 1L) +
    ggplot2::scale_y_log10(breaks = c(1, 2, 5, 10, 25, 100, 500),
                            labels = function(x) x - 1L) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = "Duplication \u00d7 Transfer coupling",
      subtitle = sprintf(
        "%d OGs \u2022 diagonal = equal D/T load \u2022 red points = top 2%% by D+T+loss",
        nrow(d)),
      x = "duplications (n_D, log+1)", y = "transfers (n_T, log+1)"
    ) +
    .dnmb_anvio_theme(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 9, colour = pal$subtitle),
      legend.position = "right"
    )
}

.panel_dtl_lorenz <- function(df, pal, ev_cols) {
  build <- function(vec, ev) {
    v <- sort(vec[!is.na(vec) & vec > 0], decreasing = TRUE)
    if (!length(v)) return(NULL)
    data.frame(
      event = ev,
      rank_frac = seq_along(v) / length(v),
      cum_frac  = cumsum(v) / sum(v)
    )
  }
  curves <- rbind(
    build(df$n_S, "S"),
    build(df$n_D, "D"),
    build(df$n_T, "T"),
    build(df$n_loss, "loss")
  )
  if (is.null(curves) || !nrow(curves)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = "Concentration"))
  }
  curves$event <- factor(curves$event, levels = c("S", "D", "T", "loss"))

  gini <- function(dd) {
    if (!nrow(dd)) return(NA_real_)
    ord <- order(dd$rank_frac)
    x <- c(0, dd$rank_frac[ord])
    y <- c(0, dd$cum_frac[ord])
    2 * sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2) - 1
  }
  gini_df <- do.call(rbind, lapply(levels(curves$event), function(ev) {
    sub <- curves[curves$event == ev, , drop = FALSE]
    data.frame(event = ev, gini = gini(sub))
  }))
  gini_df$event <- factor(gini_df$event, levels = levels(curves$event))
  gini_df$label <- sprintf("%s: concentration %.2f",
                            gini_df$event, gini_df$gini)

  ggplot2::ggplot(curves,
                   ggplot2::aes(x = .data$rank_frac, y = .data$cum_frac,
                                 colour = .data$event)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                          linetype = 2, linewidth = 0.3,
                          colour = pal$rule) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = ev_cols, name = NULL) +
    ggplot2::scale_x_continuous(labels = scales::percent) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::annotate("text",
      x = 0.02, y = 0.97, hjust = 0, vjust = 1, size = 3,
      lineheight = 1.05, colour = pal$subtitle, fontface = "italic",
      label = paste(gini_df$label, collapse = "\n")
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = "Event concentration (Lorenz curves)",
      subtitle = "y = cumulative % of events vs. x = cumulative % of OGs (desc).\nCloser to the top-left = signal concentrated in fewer OGs.",
      x = "OG rank (cumulative %)", y = "event share (cumulative %)"
    ) +
    .dnmb_anvio_theme(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 9, colour = pal$subtitle,
                                              lineheight = 1.05),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.8, "cm")
    )
}

.panel_dtl_top <- function(df, pal, ev_cols, top_n) {
  make_bar <- function(ev_col, ev_label, fill) {
    sub <- df[, c("cluster_id", ev_col), drop = FALSE]
    names(sub)[2] <- "n"
    sub <- sub[order(-sub$n), , drop = FALSE]
    sub <- sub[sub$n > 0, , drop = FALSE]
    sub <- utils::head(sub, top_n)
    if (!nrow(sub)) {
      return(
        ggplot2::ggplot() + ggplot2::theme_void() +
          ggplot2::annotate("text", x = 0, y = 0,
            label = sprintf("No OGs with %s > 0", ev_label),
            size = 3.4, colour = pal$subtitle, fontface = "italic") +
          ggplot2::labs(title = sprintf("Top by %s", ev_label)) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 11,
                                                 colour = pal$text),
            plot.background = ggplot2::element_rect(fill = pal$bg,
                                                      colour = NA))
      )
    }
    sub$label <- sprintf("OG_%07d", sub$cluster_id)
    sub$label <- factor(sub$label, levels = rev(sub$label))
    ggplot2::ggplot(sub, ggplot2::aes(x = .data$n, y = .data$label)) +
      ggplot2::geom_col(fill = fill, width = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = .data$n),
                          hjust = -0.25, size = 2.9,
                          colour = pal$text) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.15))
      ) +
      ggplot2::labs(
        title = sprintf("Top %d OGs by %s", nrow(sub), ev_label),
        x = sprintf("n_%s", ev_label), y = NULL
      ) +
      .dnmb_anvio_theme(base_size = 10) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 11,
                                             colour = pal$text),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(family = "mono", size = 7.5,
                                              colour = pal$text)
      )
  }
  bar_D    <- make_bar("n_D",    "D",    ev_cols["D"])    +
              .dnmb_panel_card(pal, c(8, 10, 8, 10))
  bar_T    <- make_bar("n_T",    "T",    ev_cols["T"])    +
              .dnmb_panel_card(pal, c(8, 10, 8, 10))
  bar_loss <- make_bar("n_loss", "loss", ev_cols["loss"]) +
              .dnmb_panel_card(pal, c(8, 10, 8, 10))
  patchwork::wrap_plots(bar_D, bar_T, bar_loss, nrow = 1L)
}
