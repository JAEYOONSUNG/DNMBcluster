#' Heatmap of per-species-pair ortholog count
#'
#' Collapses the `reconcile_relationships()` output to a species × species
#' matrix of ortholog pair counts, and plots it as a ggplot tile map. Rows
#' and columns are ordered by `hclust` on the log counts for readability.
#'
#' @param rel Tibble from `reconcile_relationships()`.
#' @param event One of `"ortholog"`, `"paralog"`, or `"xenolog"`. Selects
#'   the event type to tally. `"xenolog"` uses `is_xenolog_candidate`.
#' @param log1p_fill Logical; if TRUE, map `log1p(n)` to fill for wide dynamic range.
#' @return ggplot object.
#' @export
ortholog_matrix_plot <- function(rel,
                                 event = c("ortholog", "paralog", "xenolog"),
                                 log1p_fill = TRUE) {
  event <- match.arg(event)

  tally <- if (event == "xenolog") {
    rel[rel$is_xenolog_candidate %in% TRUE, , drop = FALSE]
  } else {
    rel[rel$event == event, , drop = FALSE]
  }
  if (!nrow(tally)) stop("ortholog_matrix_plot: no rows for event=", event)

  m <- tally %>%
    dplyr::count(species_a, species_b, name = "n") %>%
    tidyr::pivot_wider(names_from = species_b,
                       values_from = n,
                       values_fill = 0L)

  mat <- as.matrix(m[, -1])
  rownames(mat) <- m$species_a

  # Symmetrize (A vs B and B vs A should merge)
  sp <- union(rownames(mat), colnames(mat))
  full <- matrix(0L, length(sp), length(sp), dimnames = list(sp, sp))
  full[rownames(mat), colnames(mat)] <- mat
  full <- full + t(full)
  diag(full) <- 0

  # hclust ordering
  if (nrow(full) >= 3) {
    ord <- stats::hclust(stats::dist(log1p(full)))$order
    full <- full[ord, ord]
  }

  plot_df <- tibble::as_tibble(as.data.frame.table(full, responseName = "n"))
  names(plot_df) <- c("species_a", "species_b", "n")
  plot_df$species_a <- factor(plot_df$species_a, levels = rownames(full))
  plot_df$species_b <- factor(plot_df$species_b, levels = colnames(full))
  plot_df$fill <- if (log1p_fill) log1p(plot_df$n) else plot_df$n

  ggplot2::ggplot(plot_df,
                  ggplot2::aes(species_a, species_b, fill = fill)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(
      name = if (log1p_fill) sprintf("log1p(%s pairs)", event) else sprintf("%s pairs", event)
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = sprintf("%s count, species × species", event)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
