#' Faceted grid of per-orthogroup gene trees
#'
#' Plots the first `top_n` trees from `per_og_trees()` output as a small-
#' multiple ggtree grid. Useful as a QC thumbnail and as a figure in
#' comparative-genomics write-ups.
#'
#' @param og_result Tibble returned by `per_og_trees()`.
#' @param top_n How many orthogroups to show (default 24, a 6×4 grid).
#' @param ncol Grid column count.
#' @param prefer Selection order: `"largest"` (most members), `"single_copy"`,
#'   or `"random"`.
#' @return A ggplot / patchwork object. Caller saves via `ggsave()`.
#' @export
per_og_tree_grid <- function(og_result,
                             top_n = 24L,
                             ncol  = 4L,
                             prefer = c("largest", "single_copy", "random")) {
  prefer <- match.arg(prefer)
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    stop("ggtree not installed. Install from Bioconductor: BiocManager::install('ggtree')")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork not installed.")
  }

  ok <- og_result[!is.na(og_result$tree_path), , drop = FALSE]
  if (!nrow(ok)) {
    stop("per_og_tree_grid: no successful trees in og_result.")
  }

  pick <- switch(prefer,
    largest     = ok[order(-ok$n_members), , drop = FALSE],
    single_copy = ok[ok$single_copy, , drop = FALSE],
    random      = ok[sample.int(nrow(ok)), , drop = FALSE]
  )
  pick <- pick[seq_len(min(top_n, nrow(pick))), , drop = FALSE]

  # NJ can emit tiny negative branches on near-identical sequences; ggtree
  # warns on each one, which floods the log. Silence by flipping to zero.
  plots <- lapply(seq_len(nrow(pick)), function(i) {
    tr <- ape::read.tree(pick$tree_path[i])
    if (!is.null(tr$edge.length)) tr$edge.length <- pmax(tr$edge.length, 0)
    ggtree::ggtree(tr) +
      ggplot2::labs(title = sprintf("OG_%07d  n=%d", pick$cluster_id[i], pick$n_members[i])) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 7))
  })

  patchwork::wrap_plots(plots, ncol = ncol)
}
