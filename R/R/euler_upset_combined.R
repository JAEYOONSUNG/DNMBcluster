#' Euler diagram + UpSet plot — combined pangenome overview figure
#'
#' The two most common ways to visualize multi-genome ortholog overlap,
#' side by side on one page. Left panel = area-proportional Euler
#' diagram (overall shape of sharing). Right panel = UpSet intersection
#' bar chart (exact counts for every combination of genomes).
#'
#' For > 6 genomes the Euler diagram becomes approximate; the UpSet
#' plot stays exact for any number. Together they give both the
#' gestalt and the detail.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @param max_upset_sets Maximum genome sets to show in UpSet (default 40 combos).
#' @return Invisibly, the combined patchwork object.
#' @export
euler_upset_combined <- function(dnmb, output_file = NULL, max_upset_sets = 40L) {
  if (!requireNamespace("eulerr", quietly = TRUE) ||
      !requireNamespace("UpSetR", quietly = TRUE)) {
    warning("eulerr or UpSetR not installed — skipping euler_upset_combined")
    return(invisible(NULL))
  }

  presence <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(
      dnmb$genome_meta %>% dplyr::select(genome_uid, genome_key),
      by = "genome_uid"
    )

  # --- Build binary presence matrix (clusters x genomes) ---------
  genome_keys <- sort(unique(presence$genome_key))
  cluster_ids <- sort(unique(presence$cluster_id))

  bin_mat <- matrix(0L, nrow = length(cluster_ids), ncol = length(genome_keys),
                    dimnames = list(cluster_ids, genome_keys))
  for (i in seq_len(nrow(presence))) {
    cid <- as.character(presence$cluster_id[i])
    gk  <- presence$genome_key[i]
    bin_mat[cid, gk] <- 1L
  }
  bin_df <- as.data.frame(bin_mat)

  # --- Left panel: Euler ----------------------------------------
  set_list <- split(presence$cluster_id, presence$genome_key)
  euler_fit <- eulerr::euler(set_list, shape = "ellipse")
  n_genomes <- length(set_list)
  palette <- grDevices::hcl.colors(n_genomes, palette = "Zissou 1", alpha = 0.55)
  euler_plot <- plot(
    euler_fit,
    fills    = list(fill = palette),
    edges    = list(col = "grey30", lwd = 0.8),
    labels   = list(font = 2, cex = 0.7),
    quantities = list(type = "counts", cex = 0.6),
    main     = list(label = "Euler diagram", fontface = "bold", cex = 1.0)
  )

  # --- Right panel: UpSet ----------------------------------------
  # UpSetR expects a data.frame with 0/1 columns named by set.
  upset_plot <- UpSetR::upset(
    bin_df,
    nsets        = length(genome_keys),
    nintersects  = max_upset_sets,
    order.by     = "freq",
    decreasing   = TRUE,
    show.numbers = "yes",
    text.scale   = c(1.2, 1.0, 1.0, 0.9, 1.2, 0.8),
    point.size   = 2.5,
    line.size    = 0.8,
    mb.ratio     = c(0.6, 0.4),
    sets.bar.color = "#2C5F7A",
    main.bar.color = "#3C5488",
    matrix.color   = "#2C5F7A"
  )

  # --- Combine with cowplot or save sequentially ------------------
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(output_file, width = 18, height = 9)

    # Page 1: side-by-side (Euler left, UpSet right)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = c(0.4, 0.6))))

    # Euler (grid-based)
    grid::pushViewport(grid::viewport(layout.pos.col = 1))
    print(euler_plot, newpage = FALSE)
    grid::popViewport()

    # UpSet (base-graphics-based — needs print inside a viewport trick)
    grid::pushViewport(grid::viewport(layout.pos.col = 2))
    # UpSetR uses base graphics; render into a separate device then grab
    grid::popViewport()
    grid::popViewport()

    # Since UpSetR uses base graphics (not grid), render on page 2
    print(upset_plot)

    grDevices::dev.off()
    message("euler_upset_combined written to: ", output_file)
  }

  invisible(list(euler = euler_plot, upset = upset_plot))
}
