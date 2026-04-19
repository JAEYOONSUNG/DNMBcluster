#' One-page DTL summary composite
#'
#' Stitches the four DTL views into a single editorial-quality figure
#' so reviewers can read the reconciliation story without flipping
#' between sidecar PDFs:
#'
#' \itemize{
#'   \item **Top** branch-annotated DTL overlay
#'         (`plot_dtl_branch_events()`) showing where on the species
#'         tree D / T / L events accumulate.
#'   \item **Bottom-left** per-OG event histogram (`.panel_dtl_hist`).
#'   \item **Bottom-mid** D x T coupling scatter (`.panel_dtl_couple`).
#'   \item **Bottom-right** Lorenz concentration curves (`.panel_dtl_lorenz`).
#' }
#'
#' Inputs are read from `out_dir` so the function can be called as a
#' single line after `run_orthofinder_like(..., dtl = TRUE)`.
#'
#' @param out_dir Directory containing `dtl_per_og.tsv`, `dtl_events.tsv`,
#'   `dtl_losses.tsv` and `species_tree_rooted.nwk`.
#' @param out_pdf Destination PDF.
#' @param width,height PDF dimensions (inches).
#' @param verbose Echo progress.
#' @return Path to the written PDF (invisibly).
#' @export
plot_dtl_summary <- function(out_dir,
                              out_pdf = NULL,
                              width   = 16,
                              height  = 14,
                              verbose = TRUE) {
  for (pkg in c("ggplot2", "patchwork", "ape")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(sprintf("[dtl_summary] %s required; skipping", pkg))
      return(NULL)
    }
  }

  ev_tsv  <- file.path(out_dir, "dtl_events.tsv")
  po_tsv  <- file.path(out_dir, "dtl_per_og.tsv")
  ls_tsv  <- file.path(out_dir, "dtl_losses.tsv")
  sp_nwk  <- file.path(out_dir, "species_tree_rooted.nwk")
  if (!file.exists(ev_tsv) || !file.exists(po_tsv) || !file.exists(sp_nwk)) {
    warning("[dtl_summary] missing dtl_*.tsv or species tree -- skipping")
    return(NULL)
  }

  pal     <- .dnmb_anvio_pal
  ev_cols <- .dtl_ev_cols(pal)

  events <- utils::read.table(ev_tsv, sep = "\t", header = TRUE,
                               stringsAsFactors = FALSE)
  per_og <- utils::read.table(po_tsv, sep = "\t", header = TRUE,
                               stringsAsFactors = FALSE)
  losses <- if (file.exists(ls_tsv))
    utils::read.table(ls_tsv, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE) else NULL
  tree   <- ape::ladderize(ape::read.tree(sp_nwk))

  # Header stat-band: at-a-glance totals.
  header <- .panel_dtl_header(per_og, pal, ev_cols)

  # Top: branch-annotated overlay (returns its own patchwork).
  top <- plot_dtl_branch_events(events, tree, losses_tbl = losses,
                                 out_pdf = NULL)

  # Bottom row: three per-OG distribution panels with shared header.
  hist_p   <- .panel_dtl_hist(per_og, pal, ev_cols)   +
              ggplot2::labs(title = "Event distribution per OG") +
              .dnmb_panel_card(pal, c(8, 10, 8, 10))
  couple_p <- .panel_dtl_couple(per_og, pal)          +
              ggplot2::labs(title = "Duplication x Transfer") +
              .dnmb_panel_card(pal, c(8, 10, 8, 10))
  lorenz_p <- .panel_dtl_lorenz(per_og, pal, ev_cols) +
              ggplot2::labs(title = "Event concentration") +
              .dnmb_panel_card(pal, c(8, 10, 8, 10))

  bottom <- patchwork::wrap_plots(hist_p, couple_p, lorenz_p,
                                    nrow = 1, widths = c(1.4, 1, 1))

  combined <- patchwork::wrap_elements(full = header) /
              patchwork::wrap_elements(full = top) /
              patchwork::wrap_elements(full = bottom) +
    patchwork::plot_layout(heights = c(0.18, 1.55, 1)) +
    patchwork::plot_annotation(
      title    = "DTL reconciliation summary",
      caption  = sprintf(
        "header: at-a-glance totals | middle: per-branch event overlay | bottom: per-OG distributions, coupling, concentration\n%s",
        .dnmb_footer()),
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
        plot.title      = ggplot2::element_text(face = "bold", size = 18,
                                                  colour = pal$text,
                                                  family = "serif"),
        plot.caption    = ggplot2::element_text(size = 8,
                                                  colour = pal$subtitle,
                                                  hjust = 1, lineheight = 1.4)
      )
    ) &
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA)
    )

  if (is.null(out_pdf))
    out_pdf <- file.path(out_dir, "dtl_summary.pdf")
  ggplot2::ggsave(out_pdf, combined, width = width, height = height,
                   dpi = 300, limitsize = FALSE, bg = pal$bg)
  if (verbose) message(sprintf("[dtl_summary] wrote %s", out_pdf))
  invisible(out_pdf)
}


.panel_dtl_header <- function(per_og, pal, ev_cols) {
  stat_df <- data.frame(
    key   = c("#", "S", "D", "T", "L"),
    value = c(nrow(per_og), sum(per_og$n_S), sum(per_og$n_D),
              sum(per_og$n_T), sum(per_og$n_loss)),
    label = c("Orthogroups", "Speciations", "Duplications",
              "Transfers", "Losses"),
    fill  = c(pal$subtitle,
              unname(ev_cols["S"]), unname(ev_cols["D"]),
              unname(ev_cols["T"]), unname(ev_cols["loss"])),
    stringsAsFactors = FALSE
  )
  .dnmb_stat_band(stat_df, pal)
}
