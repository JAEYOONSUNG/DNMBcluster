#' Cluster presence/absence heatmap
#'
#' Roary/BPGA-style matrix plot — one row per cluster, one column per
#' genome, filled cell wherever a cluster has at least one CDS from
#' that genome. Clusters are sorted by `n_genomes` descending, so the
#' core block sits at the top and the strain-specific tail at the
#' bottom, giving an immediate visual readout of the pan-genome shape.
#'
#' For very large pan-genomes (tens of thousands of clusters) row
#' labels and ticks are hidden to keep the image legible.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return A ggplot object.
#' @export
presence_absence_heatmap <- function(dnmb, output_file = NULL) {
  cat_lookup <- dnmb$cluster_summary %>%
    dplyr::select(cluster_id, category, n_genomes)

  # cluster_id x genome_uid presence pairs from the unique CDS rows.
  long <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(
      dnmb$genome_meta %>% dplyr::select(genome_uid, genome_key),
      by = "genome_uid"
    ) %>%
    dplyr::left_join(cat_lookup, by = "cluster_id")

  # Row order: most-shared clusters at the top, rarest at the bottom.
  cluster_order <- dnmb$cluster_summary %>%
    dplyr::arrange(dplyr::desc(n_genomes), cluster_id) %>%
    dplyr::pull(cluster_id)

  long <- long %>%
    dplyr::mutate(
      cluster_id = factor(cluster_id, levels = cluster_order),
      category   = factor(category, levels = c("core", "accessory", "unique"))
    )

  palette <- c(
    core      = "#2C5F7A",
    accessory = "#F2A766",
    unique    = "#D06461"
  )

  plot <- ggplot2::ggplot(long, ggplot2::aes(x = genome_key, y = cluster_id)) +
    ggplot2::geom_tile(ggplot2::aes(fill = category)) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_discrete(breaks = NULL) +
    ggplot2::labs(
      title = "Cluster presence/absence matrix",
      subtitle = sprintf(
        "%s clusters x %s genomes, sorted by shared-genome count",
        format(length(unique(long$cluster_id)), big.mark = ","),
        length(unique(long$genome_key))
      ),
      x = NULL,
      y = "Cluster",
      fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      axis.ticks.y     = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank(),
      legend.position  = "bottom"
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    # Keep the figure tall and narrow — each genome column only needs
    # ~0.4 in, and the cluster axis benefits from vertical room. Width
    # scales modestly with genome count so runs with more strains still
    # breathe a bit.
    n_genomes <- length(unique(long$genome_key))
    fig_width <- min(max(3, 1.5 + 0.35 * n_genomes), 8)
    ggplot2::ggsave(output_file, plot, width = fig_width, height = 10, dpi = 200)
    message("presence_absence_heatmap written to: ", output_file)
  }

  invisible(plot)
}
