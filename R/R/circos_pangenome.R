#' Anvi'o-style circular presence/absence pangenome display
#'
#' The iconic pangenome figure: each concentric ring is one genome,
#' each radial slice is one gene cluster. Dark = present, light =
#' absent. Clusters are sorted by the number of genomes they appear
#' in (core at the top/left of the circle, unique at the bottom/right)
#' so the characteristic "core plateau + accessory fade + unique tail"
#' shape is immediately visible.
#'
#' Built with the ``circlize`` package for proper circular coordinate
#' geometry. The outermost track shows genome labels; inner tracks
#' show the binary matrix as colored rectangles.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return Invisible NULL (circlize draws directly to the device).
#' @export
circos_pangenome <- function(dnmb, results_dir = NULL, output_file = NULL) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    warning("circlize not installed — skipping circos_pangenome")
    return(invisible(NULL))
  }

  # --- Build presence/absence binary matrix ----------------------
  # Rows = clusters (sorted by n_genomes desc → core first),
  # Columns = genomes (sorted by genome_uid).
  presence <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(
      dnmb$genome_meta %>% dplyr::select(genome_uid, genome_key),
      by = "genome_uid"
    )

  genome_keys <- dnmb$genome_meta %>%
    dplyr::arrange(genome_uid) %>%
    dplyr::pull(genome_key)

  cluster_ids <- sort(unique(presence$cluster_id))

  # Compute n_genomes per cluster for sorting
  n_per_cluster <- presence %>%
    dplyr::count(cluster_id, name = "n_genomes") %>%
    dplyr::arrange(dplyr::desc(n_genomes), cluster_id)
  sorted_cids <- n_per_cluster$cluster_id

  # Binary matrix: rows = sorted clusters, cols = genomes
  bin_mat <- matrix(0L, nrow = length(sorted_cids), ncol = length(genome_keys),
                    dimnames = list(sorted_cids, genome_keys))
  for (i in seq_len(nrow(presence))) {
    cid <- as.character(presence$cluster_id[i])
    gk  <- presence$genome_key[i]
    if (cid %in% rownames(bin_mat)) {
      bin_mat[cid, gk] <- 1L
    }
  }

  n_clusters <- nrow(bin_mat)
  n_genomes  <- ncol(bin_mat)

  # Category colors
  cat_lookup <- n_per_cluster$n_genomes
  names(cat_lookup) <- as.character(n_per_cluster$cluster_id)
  n_total <- length(genome_keys)

  # Strain-label abbreviation (organism + strain, shortened)
  strain_labels <- dnmb$genome_meta %>%
    dplyr::arrange(genome_uid) %>%
    dplyr::mutate(
      short = dplyr::case_when(
        !is.na(strain) & nzchar(strain) ~ strain,
        !is.na(organism) ~ sub("^(\\S+\\s+\\S+).*", "\\1", organism),
        TRUE ~ genome_key
      )
    ) %>%
    dplyr::pull(short)

  # --- Color scheme: present = dark, absent = light ---------------
  col_present <- "#2C5F7A"
  col_absent  <- "#F0F0F0"

  # --- Draw with circlize ----------------------------------------
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(output_file, width = 14, height = 14)
  }
  par(mar = c(1, 1, 2, 1))

  circlize::circos.clear()
  # Open a 90° gap at the bottom for summary annotations (ANI box,
  # genome size bars, GC content). The gap spans the 6 o'clock
  # position; the pangenome ring fills the remaining 270°.
  # Sector starts at 0° (3 o'clock / right) and runs clockwise
  # for 270°, leaving the upper-right quadrant (0°→90°) empty for
  # the ANI heatmap + GC + genome-size annotation overlays.
  circlize::circos.par(
    gap.degree     = 90,
    cell.padding   = c(0, 0, 0, 0),
    track.margin   = c(0.005, 0.005),
    start.degree   = 0,
    clock.wise     = TRUE
  )

  # Single sector spanning all clusters
  circlize::circos.initialize(
    factors = "pangenome",
    xlim    = c(0, n_clusters)
  )

  # One track per genome (inside-out = first track = outermost)
  for (g in seq_len(n_genomes)) {
    gk <- genome_keys[g]
    circlize::circos.track(
      factors = "pangenome",
      ylim    = c(0, 1),
      bg.border = NA,
      track.height = 0.04,
      panel.fun = function(x, y) {
        for (i in seq_len(n_clusters)) {
          col <- if (bin_mat[i, g] == 1L) col_present else col_absent
          circlize::circos.rect(
            xleft   = i - 1,
            ybottom = 0,
            xright  = i,
            ytop    = 1,
            col     = col,
            border  = NA
          )
        }
        # Genome label on the right edge
        circlize::circos.text(
          x         = n_clusters + n_clusters * 0.01,
          y         = 0.5,
          labels    = strain_labels[g],
          facing    = "inside",
          niceFacing = TRUE,
          cex       = 0.55,
          adj       = c(0, 0.5)
        )
      }
    )
  }

  # Category annotation track (outermost = innermost in draw order)
  circlize::circos.track(
    factors = "pangenome",
    ylim = c(0, 1),
    bg.border = NA,
    track.height = 0.03,
    panel.fun = function(x, y) {
      for (i in seq_len(n_clusters)) {
        cid <- as.character(sorted_cids[i])
        ng <- cat_lookup[cid]
        col <- if (ng >= n_total) "#2C5F7A"
               else if (ng == 1)  "#D06461"
               else               "#F2A766"
        circlize::circos.rect(
          xleft = i - 1, ybottom = 0, xright = i, ytop = 1,
          col = col, border = NA
        )
      }
    }
  )

  # Title + legend for circos panel
  title(
    main = sprintf(
      "Pangenome  (%s clusters x %s genomes)",
      format(n_clusters, big.mark = ","), n_genomes
    ),
    cex.main = 1.3, font.main = 2
  )
  legend(
    "bottomleft",
    legend = c("Core", "Accessory", "Unique", "Present", "Absent"),
    fill   = c("#2C5F7A", "#F2A766", "#D06461", col_present, col_absent),
    border = "grey50", cex = 0.7, bty = "n"
  )

  circlize::circos.clear()

  # --- Gap-area overlay: ANI heatmap + GC + genome size -----------
  # The 90° gap sits in the upper-right quadrant (start.degree=90,
  # clockwise, gap at 0°→90° absolute = 3 o'clock → 12 o'clock).
  # Overlay base R graphics into that region using par(fig=...).

  meta_sorted <- dnmb$genome_meta %>% dplyr::arrange(genome_uid)
  gc_vals  <- meta_sorted$gc_percent
  size_vals <- meta_sorted$total_length
  n_g <- nrow(meta_sorted)

  # --- ANI heatmap in the gap ------------------------------------
  ani_path <- file.path(results_dir, "dnmb", "processed", "ani_matrix.parquet")
  if (file.exists(ani_path)) {
    ani_df <- tibble::as_tibble(arrow::read_parquet(ani_path))
    ani_keys <- sort(unique(c(ani_df$genome_a, ani_df$genome_b)))
    ani_mat <- matrix(NA_real_, nrow = length(ani_keys), ncol = length(ani_keys),
                      dimnames = list(ani_keys, ani_keys))
    for (r in seq_len(nrow(ani_df))) {
      ani_mat[ani_df$genome_a[r], ani_df$genome_b[r]] <- ani_df$ani_percent[r]
    }
    # Cluster
    hc <- stats::hclust(stats::as.dist(100 - ani_mat), method = "complete")
    ord <- hc$order

    # Upper-right quadrant: x = right half, y = upper half
    par(fig = c(0.55, 0.97, 0.55, 0.95), new = TRUE, mar = c(2, 0, 2, 3))
    image(
      ani_mat[ord, ord],
      col = grDevices::colorRampPalette(c("#2C5F7A","#88C0D0","#FFDC91","#E18727","#CA0020"))(100),
      zlim = c(80, 100),
      axes = FALSE
    )
    # Row labels (right side)
    n_a <- length(ord)
    axis(4, at = seq(0, 1, length.out = n_a),
         labels = strain_labels[ord], las = 2, cex.axis = 0.5, tick = FALSE, line = -0.5)
    title(main = "ANI (%)", cex.main = 0.9, font.main = 2)
    box(lwd = 0.5)
  }

  # --- GC% strip below ANI --------------------------------------
  par(fig = c(0.55, 0.97, 0.48, 0.55), new = TRUE, mar = c(0, 0, 0, 3))
  gc_range <- range(gc_vals, na.rm = TRUE)
  gc_col_fun <- grDevices::colorRampPalette(c("#FFFFCC", "#006837"))
  gc_idx <- pmax(1, pmin(100, round((gc_vals - gc_range[1]) / max(diff(gc_range), 0.1) * 99) + 1))
  gc_colors <- gc_col_fun(100)[gc_idx]

  plot(NULL, xlim = c(0, n_g), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  for (g in seq_len(n_g)) {
    rect(g - 1, 0, g, 1, col = gc_colors[g], border = "white", lwd = 0.3)
    text(g - 0.5, 0.5, sprintf("%.1f", gc_vals[g]), cex = 0.4)
  }
  mtext("GC%", side = 4, line = 0.3, cex = 0.6, las = 0)

  # --- Genome size bars below GC ---------------------------------
  par(fig = c(0.55, 0.97, 0.40, 0.48), new = TRUE, mar = c(1, 0, 0, 3))
  max_s <- max(size_vals, na.rm = TRUE)

  barplot(
    size_vals / 1e6,
    col     = "#5E81AC",
    border  = "white",
    space   = 0,
    axes    = FALSE,
    names.arg = rep("", n_g)
  )
  axis(2, cex.axis = 0.5, las = 1)
  mtext("Mb", side = 4, line = 0.3, cex = 0.6, las = 0)

  if (!is.null(output_file)) {
    grDevices::dev.off()
    message("circos_pangenome written to: ", output_file)
  }

  invisible(NULL)
}
