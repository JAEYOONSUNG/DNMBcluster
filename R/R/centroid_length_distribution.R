#' Centroid protein length distribution by pan-genome category
#'
#' Violin + boxplot showing the `aa_length` of each cluster
#' representative split into core / accessory / unique. Core
#' proteins routinely skew longer — universally conserved
#' housekeeping enzymes are structurally constrained — whereas the
#' unique bucket is dominated by short hypotheticals and truncated
#' ORFs. A clear right-shift of the core violin vs the unique violin
#' is the expected biological signal.
#'
#' Length is plotted on a log10 scale because pan-genome length
#' distributions span two or three orders of magnitude.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return A ggplot object.
#' @export
centroid_length_distribution <- function(dnmb, output_file = NULL) {
  cat_lookup <- dnmb$cluster_summary %>%
    dplyr::select(cluster_id, category, representative_uid) %>%
    dplyr::rename(protein_uid = representative_uid)

  df <- cat_lookup %>%
    dplyr::left_join(
      dnmb$id_map %>% dplyr::select(protein_uid, aa_length),
      by = "protein_uid"
    ) %>%
    dplyr::filter(!is.na(aa_length), aa_length > 0L) %>%
    dplyr::mutate(
      category = factor(category, levels = c("core", "accessory", "unique"))
    )

  palette <- .dnmb_presence_pal()

  summary_df <- df %>%
    dplyr::group_by(category) %>%
    dplyr::summarise(
      n      = dplyr::n(),
      median = stats::median(aa_length),
      top    = stats::quantile(aa_length, 0.99),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = sprintf(
        "n=%s\nmed=%d aa",
        format(n, big.mark = ","), as.integer(median)
      )
    )

  plot <- ggplot2::ggplot(
    df, ggplot2::aes(x = category, y = aa_length, fill = category)
  ) +
    ggplot2::geom_violin(alpha = 0.75, color = NA, trim = TRUE) +
    ggplot2::geom_boxplot(
      width = 0.14, fill = "white", outlier.shape = NA,
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      data = summary_df,
      ggplot2::aes(x = category, y = top, label = label),
      inherit.aes = FALSE,
      vjust = -0.2, size = 3.2
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_log10(labels = scales::comma) +
    ggplot2::labs(
      title    = "Centroid protein length by pan-genome category",
      subtitle = "One point per cluster representative; log10 length axis",
      x    = NULL,
      y    = "aa length (log10)",
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position  = "none",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, plot, width = 8, height = 6, dpi = 300)
    message("centroid_length_distribution written to: ", output_file)
  }

  invisible(plot)
}
