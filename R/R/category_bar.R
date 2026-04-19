#' Stacked bar plot of cluster categories per genome
#'
#' For each genome, counts how many core / accessory / unique clusters
#' contain at least one CDS from that genome, and draws a stacked bar.
#' Complements the flower plot with an ordered, quantitative view.
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
    # ggplot2's position_stack() draws in the REVERSE order of factor
    # levels, so listing unique → accessory → core here puts core at the
    # bottom of each bar and unique at the top, which is the pan-genome
    # convention (foundation layer first, strain-specific tip last).
    dplyr::mutate(
      category = factor(category, levels = c("unique", "accessory", "core"))
    ) %>%
    dplyr::group_by(genome_key) %>%
    dplyr::mutate(
      pct = n_clusters / sum(n_clusters),
      label = sprintf("%s\n(%.1f%%)", format(n_clusters, big.mark = ","), pct * 100)
    ) %>%
    dplyr::ungroup()

  palette <- .dnmb_presence_pal()

  plot <- ggplot2::ggplot(
    per_genome_cat,
    ggplot2::aes(x = genome_key, y = n_clusters, fill = category)
  ) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      position = ggplot2::position_stack(vjust = 0.5),
      color    = "white",
      size     = 3,
      lineheight = 0.9
    ) +
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
