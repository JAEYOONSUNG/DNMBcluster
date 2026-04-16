#' Genome-genome ortholog Jaccard heatmap
#'
#' For every ordered pair of input genomes, computes the Jaccard
#' similarity of their cluster membership sets — i.e. what fraction
#' of orthologous clusters is shared between them — and draws a
#' symmetric heatmap. Rows and columns are reordered via average-
#' linkage hierarchical clustering on `1 - Jaccard`, so closely
#' related strains sit next to each other and any obvious subclade
#' structure jumps out visually.
#'
#' The matrix is symmetric, so text labels pack two additional
#' metrics per cell — no information is wasted on redundant halves:
#'
#' - **Upper triangle** shows the Jaccard value (`0.87`).
#' - **Lower triangle** shows the raw shared-cluster count — the
#'   numerator `|A ∩ B|` — which reveals the absolute scale of the
#'   overlap that Jaccard normalizes away.
#' - **Diagonal** shows each genome's own total cluster count.
#'
#' Cell fill stays keyed to Jaccard (the primary metric) across all
#' cells, so color remains interpretable as a single gradient.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return A ggplot object.
#' @export
genome_jaccard_heatmap <- function(dnmb, output_file = NULL) {
  # Build cluster membership set per genome
  presence <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct()

  genome_keys_df <- dnmb$genome_meta %>%
    dplyr::select(genome_uid, genome_key) %>%
    dplyr::arrange(genome_uid)

  cluster_sets <- split(presence$cluster_id, presence$genome_uid)

  uids <- genome_keys_df$genome_uid
  keys <- genome_keys_df$genome_key
  n    <- length(uids)

  # Pairwise metrics: Jaccard and raw intersection size. Diagonal
  # jaccard is 1 and diagonal intersection is the genome's own
  # cluster count.
  jacc  <- matrix(0, nrow = n, ncol = n, dimnames = list(keys, keys))
  inter <- matrix(0L, nrow = n, ncol = n, dimnames = list(keys, keys))
  own   <- integer(n)
  for (i in seq_len(n)) {
    a <- cluster_sets[[as.character(uids[i])]]
    own[i] <- length(a)
    for (j in seq_len(n)) {
      b <- cluster_sets[[as.character(uids[j])]]
      inter_count <- length(intersect(a, b))
      un_count <- length(union(a, b))
      jacc[i, j]  <- if (un_count == 0L) 0 else inter_count / un_count
      inter[i, j] <- inter_count
    }
  }

  # Reorder via hclust so neighboring rows are related strains.
  ord <- if (n >= 3L) {
    stats::hclust(stats::as.dist(1 - jacc), method = "average")$order
  } else {
    seq_len(n)
  }
  jacc_ord  <- jacc[ord, ord, drop = FALSE]
  inter_ord <- inter[ord, ord, drop = FALSE]
  own_ord   <- own[ord]
  keys_ord  <- keys[ord]

  # Build long-form data with per-cell triangle assignment. Indices
  # are relative to the reordered matrix so upper/lower correspond
  # to the visual layout after rev()-ing the y-axis factor below.
  grid <- expand.grid(i = seq_len(n), j = seq_len(n))
  grid$genome_a <- keys_ord[grid$i]
  grid$genome_b <- keys_ord[grid$j]
  grid$jaccard  <- jacc_ord[cbind(grid$i, grid$j)]
  grid$inter    <- inter_ord[cbind(grid$i, grid$j)]
  grid$tri <- ifelse(
    grid$i == grid$j, "diag",
    ifelse(grid$j > grid$i, "upper", "lower")
  )
  grid$label <- ifelse(
    grid$tri == "upper",
    sprintf("%.2f", grid$jaccard),
    ifelse(
      grid$tri == "lower",
      format(grid$inter, big.mark = ","),
      format(own_ord[grid$i], big.mark = ",")
    )
  )

  grid$genome_a <- factor(grid$genome_a, levels = keys_ord)
  grid$genome_b <- factor(grid$genome_b, levels = rev(keys_ord))

  # Split into three layers so each can carry its own fill scale
  # via ggnewscale. Upper triangle = Jaccard (blue), lower triangle
  # = raw shared-cluster count (orange), diagonal = neutral grey
  # plate labeled with the genome's own cluster count.
  upper_df <- grid[grid$tri == "upper", ]
  lower_df <- grid[grid$tri == "lower", ]
  diag_df  <- grid[grid$tri == "diag", ]

  max_inter <- max(lower_df$inter)
  if (max_inter == 0L) max_inter <- 1L

  plot <- ggplot2::ggplot(mapping = ggplot2::aes(x = genome_a, y = genome_b)) +
    # --- diagonal ------------------------------------------------
    ggplot2::geom_tile(
      data = diag_df, fill = "#D0D0D0",
      color = "white", linewidth = 0.3
    ) +
    # --- upper triangle: Jaccard (blue) --------------------------
    ggplot2::geom_tile(
      data = upper_df, ggplot2::aes(fill = jaccard),
      color = "white", linewidth = 0.3
    ) +
    ggplot2::scale_fill_gradient(
      low = "#F2F2F2", high = "#2C5F7A",
      limits = c(0, 1),
      name   = "Jaccard\n(upper)",
      guide  = ggplot2::guide_colorbar(order = 2)
    ) +
    # --- lower triangle: shared cluster count (orange) -----------
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(
      data = lower_df, ggplot2::aes(fill = inter),
      color = "white", linewidth = 0.3
    ) +
    ggplot2::scale_fill_gradient(
      low = "#F2F2F2", high = "#C65A11",
      limits = c(0, max_inter),
      labels = scales::comma,
      name   = "Shared\nclusters\n(lower)",
      guide  = ggplot2::guide_colorbar(order = 1)
    ) +
    # --- cell labels ---------------------------------------------
    ggplot2::geom_text(
      data = upper_df,
      ggplot2::aes(label = label),
      color = ifelse(upper_df$jaccard > 0.35, "white", "#303030"),
      size = 3
    ) +
    ggplot2::geom_text(
      data = lower_df,
      ggplot2::aes(label = label),
      color = ifelse(lower_df$inter / max_inter > 0.35, "white", "#303030"),
      size = 3
    ) +
    ggplot2::geom_text(
      data = diag_df,
      ggplot2::aes(label = label),
      color = "#303030", size = 3, fontface = "bold"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title    = "Genome-genome ortholog Jaccard",
      subtitle = "Diagonal = genome's own cluster count",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "right"
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    side <- max(5, 1.5 + 0.55 * n)
    ggplot2::ggsave(output_file, plot, width = side, height = side, dpi = 300)
    message("genome_jaccard_heatmap written to: ", output_file)
  }

  invisible(plot)
}
