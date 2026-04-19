#' Cluster size distribution histogram
#'
#' How many clusters fall into each n_genomes bucket, colored by
#' pan-genome category. Complements the category bar by showing the
#' full width of the distribution — the characteristic U-shape of
#' most pan-genomes (tall core spike on the right, tall unique spike
#' on the left, sparse accessory valley in the middle) is immediately
#' visible.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return A ggplot object.
#' @export
cluster_size_dist <- function(dnmb, output_file = NULL) {
  df <- dnmb$cluster_summary %>%
    dplyr::count(n_genomes, category, name = "n_clusters") %>%
    dplyr::mutate(
      category = factor(category, levels = c("core", "accessory", "unique"))
    )

  palette <- .dnmb_presence_pal()

  max_n <- max(df$n_genomes)

  plot <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = n_genomes, y = n_clusters, fill = category)
  ) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::geom_text(
      ggplot2::aes(label = format(n_clusters, big.mark = ",")),
      vjust = -0.4,
      size  = 3
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_x_continuous(breaks = seq(1, max_n, 1)) +
    ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(
      title = "Cluster size distribution",
      x     = "Number of genomes per cluster",
      y     = "Cluster count",
      fill  = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position   = "bottom",
      panel.grid.minor  = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, plot, width = 10, height = 6, dpi = 300)
    message("cluster_size_dist written to: ", output_file)
  }

  invisible(plot)
}
