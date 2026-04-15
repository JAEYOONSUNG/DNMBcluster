#' Pan-core genome curve plot
#'
#' Draws the classic pan vs core genome curves with a ribbon showing
#' min/max across bootstrap permutations. Data comes from
#' `pan_core_curve.parquet` computed by the Python stage — no curve
#' fitting required at plot time.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return A ggplot object.
#' @export
pan_core_plot <- function(dnmb, output_file = NULL) {
  curve <- dnmb$pan_core_curve %>%
    dplyr::group_by(k) %>%
    dplyr::summarise(
      pan_mean  = mean(pan),
      pan_min   = min(pan),
      pan_max   = max(pan),
      core_mean = mean(core),
      core_min  = min(core),
      core_max  = max(core),
      .groups   = "drop"
    )

  plot <- ggplot2::ggplot(curve, ggplot2::aes(x = k)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = pan_min, ymax = pan_max, fill = "Pan genome"),
      alpha = 0.25
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = pan_mean, color = "Pan genome"),
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = pan_mean, color = "Pan genome"),
      size = 2
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = core_min, ymax = core_max, fill = "Core genome"),
      alpha = 0.25
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = core_mean, color = "Core genome"),
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = core_mean, color = "Core genome"),
      size = 2
    ) +
    ggplot2::scale_color_manual(values = c("Pan genome" = "#3B81A1", "Core genome" = "#D06461")) +
    ggplot2::scale_fill_manual(values  = c("Pan genome" = "#3B81A1", "Core genome" = "#D06461")) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::scale_x_continuous(breaks = seq(min(curve$k), max(curve$k), by = 1)) +
    ggplot2::labs(
      title = "Pan-core curve",
      x     = "Number of genomes",
      y     = "Gene count",
      color = NULL,
      fill  = NULL
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, plot, width = 10, height = 6, dpi = 300)
    message("pan_core_plot written to: ", output_file)
  }

  invisible(plot)
}
