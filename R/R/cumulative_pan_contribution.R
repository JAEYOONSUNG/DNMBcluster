#' Per-genome cumulative pan-genome contribution
#'
#' Walks the input genomes in `genome_uid` order (= deterministic
#' parse order) and, for each step, records how many *new* clusters
#' that genome introduces to the pan-genome that weren't already
#' contributed by an earlier genome. Complements the pan/core curve
#' — where `pan_core_plot` averages over random permutations, this
#' plot attributes each delta to a specific strain, which is what
#' you actually want when asking "which genome adds the most novel
#' content to this collection?".
#'
#' The bar height is the number of new clusters a genome adds; the
#' overlaid step line is the cumulative pan-genome size after that
#' genome is included. A near-flat bar late in the walk means the
#' last strain barely contributed anything new (redundant in the
#' current set); a tall bar late in the walk means a divergent
#' strain was added late and kicked the pan-genome open again.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return A ggplot object.
#' @export
cumulative_pan_contribution <- function(dnmb, output_file = NULL) {
  meta <- dnmb$genome_meta %>%
    dplyr::select(genome_uid, genome_key) %>%
    dplyr::arrange(genome_uid)

  presence <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct()

  per_genome <- split(presence$cluster_id, presence$genome_uid)

  seen <- integer(0)
  new_counts <- integer(nrow(meta))
  for (idx in seq_len(nrow(meta))) {
    gid  <- meta$genome_uid[idx]
    here <- per_genome[[as.character(gid)]]
    if (is.null(here)) here <- integer(0)
    new_clusters <- setdiff(here, seen)
    new_counts[idx] <- length(new_clusters)
    seen <- union(seen, new_clusters)
  }

  df <- meta
  df$new_clusters <- new_counts
  df$cumulative   <- cumsum(new_counts)
  df$order_label  <- factor(df$genome_key, levels = df$genome_key)

  max_cum <- max(df$cumulative)

  plot <- ggplot2::ggplot(df, ggplot2::aes(x = order_label, group = 1)) +
    ggplot2::geom_col(
      ggplot2::aes(y = new_clusters),
      fill = "#2C5F7A", width = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = new_clusters, label = format(new_clusters, big.mark = ",")),
      vjust = -0.4, size = 3, color = "#303030"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = cumulative),
      color = "#C65A11", linewidth = 1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = cumulative),
      color = "#C65A11", size = 2
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      expand = ggplot2::expansion(mult = c(0, 0.12)),
      sec.axis = ggplot2::sec_axis(
        ~ .,
        name   = "Cumulative pan-genome size",
        labels = scales::comma
      )
    ) +
    ggplot2::labs(
      title    = "Per-genome contribution to the pan-genome",
      subtitle = "Blue bars = new clusters added at each step    Orange line = cumulative pan size",
      x = NULL,
      y = "New clusters added"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.text.x         = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y.right  = ggplot2::element_text(color = "#C65A11"),
      axis.text.y.right   = ggplot2::element_text(color = "#C65A11"),
      panel.grid.minor    = ggplot2::element_blank(),
      panel.grid.major.x  = ggplot2::element_blank()
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, plot, width = 11, height = 6, dpi = 300)
    message("cumulative_pan_contribution written to: ", output_file)
  }

  invisible(plot)
}
