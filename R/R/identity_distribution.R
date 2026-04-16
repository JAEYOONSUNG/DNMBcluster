#' Member-vs-representative identity distribution
#'
#' Histogram of `pct_identity_fwd` across every non-centroid cluster
#' member, faceted by pan-genome category. A healthy clustering run
#' shows a tight, high-identity peak for core clusters, a broader
#' shoulder for accessory, and a dispersed tail for unique (where
#' singleton rows are excluded anyway because they have no partner).
#'
#' Use this plot to QC the realignment stage — if the core histogram
#' has significant mass below the clustering identity threshold, the
#' threshold may be too loose; if it's piled against 100%, the
#' threshold is too strict.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return A ggplot object.
#' @export
identity_distribution <- function(dnmb, output_file = NULL) {
  cat_lookup <- dnmb$cluster_summary %>%
    dplyr::select(cluster_id, category)

  df <- dnmb$clusters %>%
    dplyr::filter(!is_centroid, !is.na(pct_identity_fwd)) %>%
    dplyr::select(cluster_id, pct_identity_fwd) %>%
    dplyr::left_join(cat_lookup, by = "cluster_id") %>%
    dplyr::mutate(
      category = factor(category, levels = c("core", "accessory", "unique"))
    )

  palette <- c(
    core      = "#2C5F7A",
    accessory = "#F2A766",
    unique    = "#D06461"
  )

  # Summary stats to anchor each panel's reader intuition.
  summary_df <- df %>%
    dplyr::group_by(category) %>%
    dplyr::summarise(
      n      = dplyr::n(),
      median = median(pct_identity_fwd),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = sprintf("n=%s  median=%.1f%%", format(n, big.mark = ","), median)
    )

  plot <- ggplot2::ggplot(df, ggplot2::aes(x = pct_identity_fwd, fill = category)) +
    ggplot2::geom_histogram(binwidth = 2, boundary = 0, color = "white", linewidth = 0.2) +
    ggplot2::geom_vline(
      data = summary_df,
      ggplot2::aes(xintercept = median),
      linetype = "dashed",
      color = "grey30"
    ) +
    ggplot2::geom_text(
      data = summary_df,
      ggplot2::aes(x = Inf, y = Inf, label = label),
      hjust = 1.05, vjust = 1.5,
      size = 3.2,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 10)
    ) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::facet_wrap(~ category, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = "Member vs representative identity",
      x     = "pct_identity_fwd (%)",
      y     = "Pair count",
      fill  = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position  = "none",
      strip.background = ggplot2::element_rect(fill = "#F2F2F2"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, plot, width = 10, height = 8, dpi = 300)
    message("identity_distribution written to: ", output_file)
  }

  invisible(plot)
}
