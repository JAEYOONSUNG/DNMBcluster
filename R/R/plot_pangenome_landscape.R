#' Pan-genome landscape composite (HOG occupancy + species-pair relationships)
#'
#' Stacks the two N-by-many heatmaps that describe the pan-genome from
#' complementary angles into one editorial figure:
#'
#' \itemize{
#'   \item **Top** root-HOG copy-number heatmap with species tree.
#'         Reads how many copies of each HOG sit in each genome.
#'   \item **Bottom** species-pair triptych (orthologs / paralogs /
#'         xenolog candidates). Reads how genomes relate pairwise.
#' }
#'
#' Both panels are produced via their respective standalone functions
#' with `out_pdf = FALSE` so this composite remains a thin orchestrator.
#'
#' @param out_dir Directory produced by `run_orthofinder_like()`.
#' @param out_pdf Destination PDF.
#' @param width,height PDF dimensions (inches). Auto-scaled if NULL.
#' @param verbose Echo progress.
#' @return Path to the written PDF (invisibly), or `NULL` on failure.
#' @export
plot_pangenome_landscape <- function(out_dir,
                                      out_pdf = NULL,
                                      width   = NULL,
                                      height  = NULL,
                                      verbose = TRUE) {
  for (pkg in c("ggplot2", "patchwork", "ape")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(sprintf("[pangenome_landscape] %s required; skipping", pkg))
      return(NULL)
    }
  }
  pal <- .dnmb_anvio_pal

  hog_p <- tryCatch(
    plot_hog_occupancy_heatmap(out_dir, out_pdf = FALSE, verbose = FALSE),
    error = function(e) NULL
  )
  rel_p <- tryCatch(
    plot_relationships_heatmap(out_dir, out_pdf = FALSE, verbose = FALSE),
    error = function(e) NULL
  )
  if (!is.null(hog_p)) hog_p <- hog_p & .dnmb_panel_card(pal, c(8, 10, 8, 10))
  if (!is.null(rel_p)) rel_p <- rel_p & .dnmb_panel_card(pal, c(8, 10, 8, 10))
  if (is.null(hog_p) && is.null(rel_p)) {
    warning("[pangenome_landscape] both panels unavailable -- skipping")
    return(NULL)
  }

  header <- .pg_header(out_dir, pal)

  parts <- list(hog_p, rel_p)
  parts <- parts[!vapply(parts, is.null, logical(1))]

  body <- patchwork::wrap_plots(
    lapply(parts, patchwork::wrap_elements),
    ncol = 1,
    heights = if (length(parts) == 2L) c(1.05, 1) else 1
  )

  combined <- patchwork::wrap_elements(full = header) /
              patchwork::wrap_elements(full = body) +
    patchwork::plot_layout(heights = c(0.12, 2.1)) +
    patchwork::plot_annotation(
      title    = "Pan-genome landscape",
      caption  = sprintf(
        "header: at-a-glance totals | top: HOG occupancy & copy number | bottom: species-pair relationships\n%s",
        .dnmb_footer()),
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
        plot.title      = ggplot2::element_text(face = "bold", size = 18,
                                                  colour = pal$text,
                                                  family = "serif"),
        plot.caption    = ggplot2::element_text(size = 8,
                                                  colour = pal$subtitle,
                                                  hjust = 1, lineheight = 1.4)
      )
    ) &
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA)
    )

  if (is.null(out_pdf)) out_pdf <- file.path(out_dir, "pangenome_landscape.pdf")

  sp_nwk <- file.path(out_dir, "species_tree_rooted.nwk")
  n_sp <- if (file.exists(sp_nwk))
    length(ape::read.tree(sp_nwk)$tip.label) else 8L
  if (is.null(width))  width  <- max(13, n_sp * 1.0 + 5)
  if (is.null(height)) height <- max(12, n_sp * 0.55 + 9)
  width  <- min(width, 40)
  height <- min(height, 36)

  ggplot2::ggsave(out_pdf, combined, width = width, height = height,
                   limitsize = FALSE, bg = pal$bg)
  if (verbose) message(sprintf("[pangenome_landscape] wrote %s", out_pdf))
  invisible(out_pdf)
}


.pg_header <- function(out_dir, pal) {
  sp_nwk  <- file.path(out_dir, "species_tree_rooted.nwk")
  hog_dir <- file.path(out_dir, "HOGs")
  rel_pq  <- file.path(out_dir, "relationships.parquet")

  n_sp <- if (file.exists(sp_nwk))
    length(ape::read.tree(sp_nwk)$tip.label) else NA_integer_

  n_hogs <- 0L
  hog_files <- list.files(hog_dir, pattern = "\\.tsv$", full.names = TRUE)
  for (p in hog_files) n_hogs <- n_hogs + max(0L, length(readLines(p)) - 3L)

  n_ortho <- n_para <- n_xeno <- NA_integer_
  if (file.exists(rel_pq) && requireNamespace("arrow", quietly = TRUE)) {
    rel <- tryCatch(arrow::read_parquet(rel_pq), error = function(e) NULL)
    if (!is.null(rel) && nrow(rel)) {
      if ("event" %in% names(rel)) {
        n_ortho <- sum(rel$event == "ortholog", na.rm = TRUE)
        n_para  <- sum(rel$event == "paralog",  na.rm = TRUE)
      }
      if ("is_xenolog_candidate" %in% names(rel))
        n_xeno <- sum(rel$is_xenolog_candidate %in% TRUE)
    }
  }

  keys <- c("G", "H")
  vals <- c(if (is.na(n_sp)) 0L else n_sp, n_hogs)
  labs <- c("Genomes", "Root HOGs")
  fills <- c(pal$navy, pal$blue_cool)
  if (!is.na(n_ortho)) {
    keys <- c(keys, "O"); vals <- c(vals, n_ortho)
    labs <- c(labs, "Ortholog pairs"); fills <- c(fills, pal$orange)
  }
  if (!is.na(n_para)) {
    keys <- c(keys, "P"); vals <- c(vals, n_para)
    labs <- c(labs, "Paralog pairs"); fills <- c(fills, pal$red)
  }
  if (!is.na(n_xeno)) {
    keys <- c(keys, "X"); vals <- c(vals, n_xeno)
    labs <- c(labs, "Xenolog cand."); fills <- c(fills, pal$red_strong)
  }

  stat_df <- data.frame(key = keys, value = vals, label = labs,
                         fill = fills, stringsAsFactors = FALSE)
  .dnmb_stat_band(stat_df, pal)
}
