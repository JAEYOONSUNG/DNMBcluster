#' ANI + POCP genome-genome similarity heatmaps
#'
#' Reads `dnmb/processed/ani_matrix.parquet` and
#' `dnmb/processed/pocp_matrix.parquet` (written by the Python
#' similarity stage — skani for ANI, cluster-based approximation for
#' POCP) and draws each as a pheatmap with automatic hierarchical
#' clustering on `1 - metric/100`.
#'
#' Species-boundary (ANI >= 95%) and genus-boundary (POCP >= 60%)
#' lines are annotated in the title so users don't have to eyeball
#' thresholds. Row and column ordering comes from a single hclust
#' so both matrices are comparable visually when viewed side-by-side.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level pipeline results directory.
#' @param output_file_ani PDF path for the ANI heatmap.
#' @param output_file_pocp PDF path for the POCP heatmap.
#' @return Named list with the two pheatmap objects (or NULL for any
#'   metric whose parquet file is absent).
#' @export
ani_pocp_heatmap <- function(dnmb, results_dir, output_file_ani = NULL,
                             output_file_pocp = NULL) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    warning("pheatmap not installed — skipping ani_pocp_heatmap")
    return(invisible(NULL))
  }

  processed_dir <- file.path(results_dir, "dnmb", "processed")

  load_matrix <- function(parquet_name, value_col) {
    path <- file.path(processed_dir, parquet_name)
    if (!file.exists(path)) return(NULL)
    df <- tibble::as_tibble(arrow::read_parquet(path))
    keys <- sort(unique(c(df$genome_a, df$genome_b)))
    n <- length(keys)
    mat <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(keys, keys))
    for (i in seq_len(nrow(df))) {
      mat[df$genome_a[i], df$genome_b[i]] <- df[[value_col]][i]
    }
    mat
  }

  ani_mat  <- load_matrix("ani_matrix.parquet",  "ani_percent")
  pocp_mat <- load_matrix("pocp_matrix.parquet", "pocp_percent")

  result <- list(ani = NULL, pocp = NULL)

  render_one <- function(mat, title, low, mid, high, threshold_line, output_file) {
    if (is.null(mat)) return(NULL)
    dist_obj <- stats::as.dist(1 - mat / 100)
    col_pal <- grDevices::colorRampPalette(c(low, mid, high))(100)

    p <- pheatmap::pheatmap(
      mat,
      color              = col_pal,
      clustering_distance_rows = dist_obj,
      clustering_distance_cols = dist_obj,
      clustering_method  = "average",
      display_numbers    = TRUE,
      number_format      = "%.1f",
      number_color       = "#202020",
      fontsize_number    = 7,
      fontsize           = 10,
      border_color       = "#FFFFFF",
      main               = title,
      silent             = is.null(output_file)
    )

    if (!is.null(output_file)) {
      dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
      # pheatmap returns a grob; save via its built-in .gtable.
      grDevices::pdf(output_file, width = 9, height = 8)
      grid::grid.newpage()
      grid::grid.draw(p$gtable)
      grDevices::dev.off()
      message("ani_pocp_heatmap written to: ", output_file)
    }
    p
  }

  if (!is.null(ani_mat)) {
    result$ani <- render_one(
      ani_mat,
      title = sprintf(
        "ANI (skani)  — species boundary at 95%%  |  %d genomes",
        nrow(ani_mat)
      ),
      low  = "#F5F5F5", mid = "#8FB4D4", high = "#2C5F7A",
      threshold_line = 95,
      output_file = output_file_ani
    )
  }

  if (!is.null(pocp_mat)) {
    result$pocp <- render_one(
      pocp_mat,
      title = sprintf(
        "POCP  — genus boundary at 60%%  |  %d genomes",
        nrow(pocp_mat)
      ),
      low  = "#F5F5F5", mid = "#F2B672", high = "#C65A11",
      threshold_line = 60,
      output_file = output_file_pocp
    )
  }

  invisible(result)
}
