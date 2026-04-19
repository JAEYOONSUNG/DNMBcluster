#' Render the OrthoFinder-parity HTML report
#'
#' Wraps the bundled Rmd template in `inst/rmd/orthofinder_report.Rmd`.
#' The report summarizes OG count, single-copy OG count, the rooted
#' species tree, per-node HOG counts, and (if available) the pan/core
#' curve. Output is a single self-contained HTML file.
#'
#' @param out_dir Directory produced by `run_orthofinder_like()`
#'   (contains `orthogroup_trees.tsv`, `species_tree_rooted.nwk`,
#'   `HOGs/`, …).
#' @param output_file Path for the rendered HTML. Defaults to
#'   `<out_dir>/orthofinder_report.html`.
#' @param run_title Title shown at the top of the report.
#' @param dnmb_dir Optional path to the DNMB processed directory that
#'   holds `pan_core_curve.parquet`. When NULL the template searches
#'   plausible locations around `out_dir`.
#' @return Path to the rendered HTML (invisibly).
#' @export
render_orthofinder_report <- function(out_dir,
                                      output_file = NULL,
                                      run_title = "DNMBcluster run",
                                      dnmb_dir = NULL) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("render_orthofinder_report requires the 'rmarkdown' Suggests package.")
  }
  tmpl <- system.file("rmd", "orthofinder_report.Rmd", package = "DNMBcluster")
  if (!nzchar(tmpl)) stop("orthofinder_report.Rmd template not found in installed package.")
  if (is.null(output_file)) {
    output_file <- file.path(out_dir, "orthofinder_report.html")
  }
  rmarkdown::render(
    tmpl,
    output_file = output_file,
    params = list(out_dir = normalizePath(out_dir, mustWork = TRUE),
                  dnmb_dir = if (is.null(dnmb_dir)) NULL
                              else normalizePath(dnmb_dir, mustWork = TRUE),
                  run_title = run_title),
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )
  invisible(output_file)
}
