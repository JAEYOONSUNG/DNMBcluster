#' Run the full DNMBcluster plotting suite on a results directory
#'
#' Entry point invoked by the Python CLI. Loads the DNMB Parquet
#' outputs, then generates every plot defined in the package into
#' `<results_dir>/plots/`.
#'
#' @param results_dir Path produced by `dnmbcluster run`.
#' @param output_dir Optional override; defaults to `<results_dir>/plots`.
#' @return A named list of ggplot objects.
#' @export
run_dnmb_plot <- function(results_dir, output_dir = NULL) {
  if (is.null(output_dir)) {
    output_dir <- file.path(results_dir, "plots")
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  dnmb <- load_dnmb(results_dir)

  plots <- list()

  plots$flower <- tryCatch(
    flower_plot(dnmb, output_file = file.path(output_dir, "flower_plot.pdf")),
    error = function(e) {
      warning("flower_plot failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$pan_core <- tryCatch(
    pan_core_plot(dnmb, output_file = file.path(output_dir, "pan_core_plot.pdf")),
    error = function(e) {
      warning("pan_core_plot failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$category_bar <- tryCatch(
    category_bar_plot(dnmb, output_file = file.path(output_dir, "category_bar.pdf")),
    error = function(e) {
      warning("category_bar_plot failed: ", conditionMessage(e))
      NULL
    }
  )

  message("DNMBcluster plots saved to: ", output_dir)
  invisible(plots)
}
