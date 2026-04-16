#' Core gene conservation histogram
#'
#' For every core cluster (present in all input genomes), computes
#' the median member-vs-representative identity across its
#' non-centroid members and histograms the result. A well-behaved
#' core block produces a single tight peak near 100% identity;
#' the lower tail picks out core genes under relaxed selection
#' (chaperones, surface proteins, phage remnants) and candidates
#' for clustering artifacts (paralog fusions, truncated CDSs).
#'
#' A dashed vertical line marks the overall median so the bulk of
#' the distribution is easy to locate.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return A ggplot object.
#' @export
core_gene_conservation <- function(dnmb, output_file = NULL) {
  core_clusters <- dnmb$cluster_summary %>%
    dplyr::filter(category == "core") %>%
    dplyr::select(cluster_id)

  df <- dnmb$clusters %>%
    dplyr::filter(!is_centroid, !is.na(pct_identity_fwd)) %>%
    dplyr::inner_join(core_clusters, by = "cluster_id") %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(
      median_identity = stats::median(pct_identity_fwd),
      n_members       = dplyr::n(),
      .groups         = "drop"
    )

  if (nrow(df) == 0L) {
    warning("core_gene_conservation: no core clusters with non-centroid members")
    return(invisible(NULL))
  }

  overall_median <- stats::median(df$median_identity)
  n_core <- nrow(df)

  plot <- ggplot2::ggplot(df, ggplot2::aes(x = median_identity)) +
    ggplot2::geom_histogram(
      binwidth = 1, boundary = 100,
      fill = "#2C5F7A", color = "white", linewidth = 0.2
    ) +
    ggplot2::geom_vline(
      xintercept = overall_median,
      linetype = "dashed", color = "#C65A11", linewidth = 0.7
    ) +
    ggplot2::annotate(
      "text",
      x = overall_median, y = Inf,
      label = sprintf("overall median = %.1f%%", overall_median),
      hjust = -0.05, vjust = 1.6,
      color = "#C65A11", size = 3.4
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, 10)) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(
      title    = sprintf(
        "Core gene conservation  (%s core clusters)",
        format(n_core, big.mark = ",")
      ),
      subtitle = "Per-cluster median member-vs-representative identity",
      x = "Median pct_identity_fwd (%)",
      y = "Core cluster count"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, plot, width = 10, height = 6, dpi = 300)
    message("core_gene_conservation written to: ", output_file)
  }

  invisible(plot)
}
