#' Stacked bar plot of cluster categories per genome
#'
#' For each genome, counts how many core / soft_core / shell / cloud
#' clusters contain at least one CDS from that genome, and draws a
#' stacked bar. Complements the flower plot with an ordered,
#' quantitative view.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return A ggplot object.
#' @export
category_bar_plot <- function(dnmb, output_file = NULL) {
  # cluster_id -> category lookup
  cat_lookup <- dnmb$cluster_summary %>%
    dplyr::select(cluster_id, category)

  per_genome_cat <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(cat_lookup, by = "cluster_id") %>%
    dplyr::count(genome_uid, category, name = "n_clusters") %>%
    dplyr::left_join(
      dnmb$genome_meta %>% dplyr::select(genome_uid, genome_key),
      by = "genome_uid"
    ) %>%
    dplyr::mutate(
      category = factor(category, levels = c("core", "soft_core", "shell", "cloud"))
    )

  palette <- c(
    core      = "#2C5F7A",
    soft_core = "#5B9CBD",
    shell     = "#F2A766",
    cloud     = "#D06461"
  )

  plot <- ggplot2::ggplot(
    per_genome_cat,
    ggplot2::aes(x = genome_key, y = n_clusters, fill = category)
  ) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(
      title = "Cluster categories per genome",
      x     = NULL,
      y     = "Cluster count",
      fill  = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position  = "bottom",
      panel.grid.major.x = ggplot2::element_blank()
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, plot, width = 10, height = 6, dpi = 300)
    message("category_bar_plot written to: ", output_file)
  }

  invisible(plot)
}
