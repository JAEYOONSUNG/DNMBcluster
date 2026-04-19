#' Stitch every DNMBcluster composite into one editorial report PDF
#'
#' Convenience wrapper that renders the standard family of composites
#' produced by `run_orthofinder_like(..., dtl = TRUE)` and concatenates
#' them into a single multi-page PDF suitable for handing to a reviewer
#' or a co-author who would rather not flip through six sidecar files.
#'
#' Page order:
#'
#' \itemize{
#'   \item p.1 `orthofinder_overview.pdf` -- run-level snapshot.
#'   \item p.2 `pangenome_landscape.pdf` -- HOG occupancy + species pairs.
#'   \item p.3 `dtl_summary.pdf` -- DTL reconciliation overview.
#'   \item p.4 `dtl_top_ogs.pdf` -- top D / T / loss OGs.
#'   \item p.5 `duplication_burden.pdf` -- where duplications fall.
#'   \item p.6 `per_og_tree_grid.pdf` -- per-OG tree thumbnails.
#' }
#'
#' Missing composites are skipped silently (with a `verbose` note), so
#' the report still assembles for partial runs.
#'
#' @param out_dir Directory produced by `run_orthofinder_like()`.
#' @param out_pdf Destination PDF. Defaults to `<out_dir>/dnmb_report.pdf`.
#' @param refresh If `TRUE` (default), each composite is re-rendered
#'   before assembly so the report always reflects the current data.
#'   If `FALSE`, only existing PDFs in `out_dir` are stitched.
#' @param verbose Echo progress.
#' @return Path to the written report (invisibly), or `NULL` on failure.
#' @export
dnmb_report <- function(out_dir,
                         out_pdf = NULL,
                         refresh = TRUE,
                         verbose = TRUE) {
  if (!requireNamespace("qpdf", quietly = TRUE)) {
    warning("[dnmb_report] requires the 'qpdf' package; skipping")
    return(NULL)
  }

  steps <- list(
    list(file = "orthofinder_overview.pdf", fn = plot_orthofinder_overview),
    list(file = "pangenome_landscape.pdf",  fn = plot_pangenome_landscape),
    list(file = "dtl_summary.pdf",          fn = plot_dtl_summary),
    list(file = "dtl_top_ogs.pdf",          fn = plot_dtl_top_ogs),
    list(file = "duplication_burden.pdf",   fn = plot_duplication_burden),
    list(file = "per_og_tree_grid.pdf",     fn = NULL)  # rendered by caller
  )

  pages <- character()
  for (st in steps) {
    p <- file.path(out_dir, st$file)
    if (refresh && !is.null(st$fn)) {
      tryCatch(st$fn(out_dir, verbose = FALSE),
                error = function(e) {
                  if (verbose) message(sprintf(
                    "[dnmb_report] %s render failed: %s",
                    st$file, conditionMessage(e)))
                })
    }
    if (file.exists(p)) {
      pages <- c(pages, p)
      if (verbose) message(sprintf("[dnmb_report] + %s", basename(p)))
    } else if (verbose) {
      message(sprintf("[dnmb_report] - %s (missing, skipped)",
                       basename(p)))
    }
  }

  if (!length(pages)) {
    warning("[dnmb_report] no composite PDFs found -- nothing to stitch")
    return(NULL)
  }

  if (is.null(out_pdf)) out_pdf <- file.path(out_dir, "dnmb_report.pdf")
  qpdf::pdf_combine(input = pages, output = out_pdf)
  if (verbose) message(sprintf(
    "[dnmb_report] wrote %s (%d pages)", out_pdf, length(pages)))
  invisible(out_pdf)
}
