#' Gene-content MDS scatter
#'
#' Projects each input genome into 2D via classical multidimensional
#' scaling on the `1 - Jaccard` distance of cluster membership.
#' Strains with near-identical orthologous content collapse onto the
#' same point; strain-level subclades form tight groupings; outliers
#' (divergent assemblies or reclassifications) pop out immediately.
#'
#' Binary presence/absence data makes classical MDS the natural
#' choice — unlike PCA on raw 0/1 columns, it respects set-overlap
#' distance directly. The `Dim 1 / Dim 2` variance explanations in
#' the axis labels are a qualitative readout of how faithful the 2D
#' embedding is.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return A ggplot object.
#' @export
gene_content_mds <- function(dnmb, output_file = NULL) {
  presence <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct()

  meta <- dnmb$genome_meta %>%
    dplyr::select(genome_uid, genome_key) %>%
    dplyr::arrange(genome_uid)

  uids <- meta$genome_uid
  keys <- meta$genome_key
  n    <- length(uids)

  if (n < 3L) {
    warning("gene_content_mds: need at least 3 genomes, got ", n)
    return(invisible(NULL))
  }

  cluster_sets <- split(presence$cluster_id, presence$genome_uid)

  jacc <- matrix(0, nrow = n, ncol = n, dimnames = list(keys, keys))
  for (i in seq_len(n)) {
    a <- cluster_sets[[as.character(uids[i])]]
    for (j in seq_len(n)) {
      b <- cluster_sets[[as.character(uids[j])]]
      un <- length(union(a, b))
      jacc[i, j] <- if (un == 0L) 0 else length(intersect(a, b)) / un
    }
  }

  mds <- stats::cmdscale(stats::as.dist(1 - jacc), k = 2, eig = TRUE)
  coords <- data.frame(
    genome_key = keys,
    Dim1 = mds$points[, 1],
    Dim2 = mds$points[, 2],
    stringsAsFactors = FALSE
  )

  # Variance-explained proxy — cmdscale eigenvalues can be negative
  # for non-Euclidean distances, so normalize against the absolute
  # sum to get an honest qualitative percentage.
  ev       <- mds$eig
  total_ev <- sum(abs(ev))
  pct1     <- if (total_ev > 0) 100 * ev[1] / total_ev else NA_real_
  pct2     <- if (total_ev > 0) 100 * ev[2] / total_ev else NA_real_

  plot <- ggplot2::ggplot(
    coords, ggplot2::aes(x = Dim1, y = Dim2, label = genome_key)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
    ggplot2::geom_point(size = 3.5, color = "#2C5F7A", alpha = 0.85) +
    ggplot2::geom_text(
      hjust = -0.12, vjust = -0.3,
      size  = 3.2, color = "#303030"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title    = "Gene content MDS",
      subtitle = "Classical MDS on (1 - Jaccard) distance of cluster membership",
      x = sprintf("Dim 1 (%.1f%%)", pct1),
      y = sprintf("Dim 2 (%.1f%%)", pct2)
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, plot, width = 8, height = 7, dpi = 300)
    message("gene_content_mds written to: ", output_file)
  }

  invisible(plot)
}
