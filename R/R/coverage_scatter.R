#' Member vs representative coverage scatter
#'
#' Every non-centroid cluster member is plotted once, with its
#' MMseqs2 `member_coverage` (query-side) on the x-axis and
#' `rep_coverage` (target-side) on the y-axis. Points on the
#' diagonal are alignments between length-matched sequences; points
#' far from the diagonal flag length-disparate pairs where one
#' sequence aligns over a small fraction of the other (a common
#' sign of short-domain matches against long multi-domain proteins,
#' or truncated/fragmented CDSs).
#'
#' Colored by pan-genome category so core vs accessory vs unique
#' disparate-pair patterns can be compared at a glance.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return A ggplot object.
#' @export
coverage_scatter <- function(dnmb, output_file = NULL) {
  cat_lookup <- dnmb$cluster_summary %>%
    dplyr::select(cluster_id, category)

  df <- dnmb$clusters %>%
    dplyr::filter(
      !is_centroid,
      !is.na(member_coverage),
      !is.na(rep_coverage)
    ) %>%
    dplyr::select(cluster_id, member_coverage, rep_coverage) %>%
    dplyr::left_join(cat_lookup, by = "cluster_id") %>%
    dplyr::mutate(
      category = factor(category, levels = c("core", "accessory", "unique"))
    )

  palette <- .dnmb_presence_pal()

  plot <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = member_coverage, y = rep_coverage, color = category)
  ) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", color = "grey40"
    ) +
    ggplot2::geom_point(alpha = 0.35, size = 1) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::coord_fixed() +
    ggplot2::facet_wrap(~ category, ncol = 3) +
    ggplot2::labs(
      title = "Member vs representative coverage",
      subtitle = "Points off the diagonal indicate length-disparate alignments",
      x = "member_coverage  (query side)",
      y = "rep_coverage  (target side)"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position   = "none",
      strip.background  = ggplot2::element_rect(fill = "#F2F2F2"),
      panel.grid.minor  = ggplot2::element_blank()
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, plot, width = 12, height = 5, dpi = 300)
    message("coverage_scatter written to: ", output_file)
  }

  invisible(plot)
}
