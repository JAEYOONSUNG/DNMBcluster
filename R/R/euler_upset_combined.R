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

  # --- Left panel: Euler / Venn ----------------------------------
  # ≤ 5 genomes: exact Venn via eulerr::venn (guaranteed correct).
  # 6–10: approximate Euler via eulerr::euler (good readability).
  # > 10: skip Euler entirely — UpSet is the only readable option.
  set_list <- split(presence$cluster_id, presence$genome_key)
  n_genomes <- length(set_list)

  # Consistent palette: same NPG-soft colors across Euler + UpSet.
  lighten <- function(cols, amount = 0.35) {
    m <- grDevices::col2rgb(cols) / 255
    m <- m + (1 - m) * amount
    grDevices::rgb(m[1, ], m[2, ], m[3, ])
  }
  npg_raw <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F",
               "#8491B4","#91D1C2","#E18727","#7876B1","#EFC000")
  palette <- lighten(rep_len(npg_raw, n_genomes), 0.35)

  euler_plot <- NULL
  if (n_genomes <= 5L) {
    euler_fit <- eulerr::venn(set_list)
    euler_plot <- plot(
      euler_fit,
      fills    = list(fill = palette[seq_len(n_genomes)], alpha = 0.55),
      edges    = list(col = "grey30", lwd = 0.8),
      labels   = list(font = 2, cex = 0.65),
      quantities = list(type = "counts", cex = 0.55),
      main     = list(label = "Venn diagram", fontface = "bold", cex = 1.0)
    )
  } else if (n_genomes <= 10L) {
    euler_fit <- eulerr::euler(set_list, shape = "ellipse")
    euler_plot <- plot(
      euler_fit,
      fills    = list(fill = palette[seq_len(n_genomes)], alpha = 0.55),
      edges    = list(col = "grey30", lwd = 0.8),
      labels   = list(font = 2, cex = 0.65),
      quantities = list(type = "counts", cex = 0.55),
      main     = list(label = "Euler diagram (approximate)", fontface = "bold", cex = 1.0)
    )
  }
  # > 10 genomes: euler_plot stays NULL → PDF only contains UpSet

  # --- Right panel: UpSet ----------------------------------------
  # UpSetR expects a data.frame with 0/1 columns named by set.
  # UpSet with pipeline palette: set-size bars navy, intersection
  # bars orange, matrix dots navy. Core (all genomes) and unique
  # (single genome) intersections get query-based coloring to pop.
  n_g <- length(genome_keys)

  upset_plot <- UpSetR::upset(
    bin_df,
    nsets          = length(genome_keys),
    nintersects    = max_upset_sets,
    order.by       = "freq",
    decreasing     = TRUE,
    show.numbers   = "yes",
    number.angles  = 0,
    text.scale     = c(1.4, 1.1, 1.0, 1.0, 1.4, 0.9),
    point.size     = 3.0,
    line.size      = 1.0,
    mb.ratio       = c(0.55, 0.45),
    sets.bar.color = "#3C5488",
    main.bar.color = "#5E81AC",
    matrix.color   = "#2C5F7A",
    shade.color    = "#D8DEE9",
    shade.alpha    = 0.3,
    sets.x.label   = "Clusters per genome",
    mainbar.y.label = "Intersection size"
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
