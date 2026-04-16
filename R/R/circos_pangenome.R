#' Circular pangenome with a flat annotation table in the gap
#'
#' Draws a 270-degree presence/absence circos plot and fills the remaining
#' 90-degree gap with a cartesian annotation table whose rows are aligned to
#' the exact ring heights reported by `circlize`.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level results directory.
#' @param output_file Optional PDF path.
#' @export
circos_pangenome <- function(dnmb, results_dir = NULL, output_file = NULL) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    warning("circlize not installed")
    return(invisible(NULL))
  }

  polar_to_xy <- getFromNamespace("polar2Cartesian", "circlize")

  range_or_default <- function(x, default = c(0, 1)) {
    x <- x[is.finite(x)]
    if (!length(x)) default else range(x)
  }

  palette_value <- function(x, palette, limits) {
    if (!is.finite(x)) {
      return("#F0F0F0")
    }
    span <- diff(limits)
    if (!is.finite(span) || span <= 0) {
      return(palette[ceiling(length(palette) / 2)])
    }
    idx <- round((x - limits[1]) / span * (length(palette) - 1)) + 1L
    palette[pmax(1L, pmin(length(palette), idx))]
  }

  fmt_num <- function(x, digits = 1, suffix = "") {
    ifelse(is.finite(x), paste0(formatC(x, format = "f", digits = digits), suffix), "NA")
  }

  fmt_int <- function(x) {
    ifelse(is.finite(x), prettyNum(round(x), big.mark = ",", preserve.width = "none"), "NA")
  }

  presence <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(
      dnmb$genome_meta %>% dplyr::select(genome_uid, genome_key),
      by = "genome_uid"
    )

  tree_path <- if (!is.null(results_dir)) {
    file.path(results_dir, "dnmb", "processed", "phylo_tree.nwk")
  } else {
    ""
  }

  if (file.exists(tree_path) && requireNamespace("ape", quietly = TRUE)) {
    tree <- ape::ladderize(ape::read.tree(tree_path))
    genome_keys <- tree$tip.label[tree$tip.label %in% dnmb$genome_meta$genome_key]
    genome_keys <- c(genome_keys, setdiff(dnmb$genome_meta$genome_key, genome_keys))
  } else {
    genome_keys <- dnmb$genome_meta %>%
      dplyr::arrange(genome_uid) %>%
      dplyr::pull(genome_key)
  }

  n_per <- presence %>%
    dplyr::count(cluster_id, name = "ng") %>%
    dplyr::arrange(dplyr::desc(ng), cluster_id)

  sorted_cids <- n_per$cluster_id
  n_clusters <- length(sorted_cids)
  n_genomes <- length(genome_keys)

  if (!n_clusters || !n_genomes) {
    warning("No clusters or genomes available for circos plotting")
    return(invisible(NULL))
  }

  bin_mat <- matrix(
    0L,
    nrow = n_clusters,
    ncol = n_genomes,
    dimnames = list(sorted_cids, genome_keys)
  )

  for (i in seq_len(nrow(presence))) {
    cid <- as.character(presence$cluster_id[i])
    gk <- presence$genome_key[i]
    if (!is.na(gk) && cid %in% rownames(bin_mat) && gk %in% colnames(bin_mat)) {
      bin_mat[cid, gk] <- 1L
    }
  }

  meta <- dnmb$genome_meta[match(genome_keys, dnmb$genome_meta$genome_key), ]
  organism_short <- ifelse(
    !is.na(meta$organism) & nzchar(meta$organism),
    sub("^(\\S+\\s+\\S+).*", "\\1", meta$organism),
    meta$genome_key
  )
  strain_labels <- ifelse(!is.na(meta$strain) & nzchar(meta$strain), meta$strain, organism_short)
  strain_labels[is.na(strain_labels) | !nzchar(strain_labels)] <- meta$genome_key[is.na(strain_labels) | !nzchar(strain_labels)]

  gc_vals <- meta$gc_percent
  size_vals <- meta$total_length / 1e6
  cds_vals <- meta$n_cds

  col_present <- "#2F6F6F"
  col_absent <- "#F3F3F3"
  col_core <- "#1F4E5F"
  col_accessory <- "#F0A35B"
  col_unique <- "#D06461"

  cluster_track_h <- 0.026
  genome_track_h <- max(0.008, min(0.04, (0.78 - cluster_track_h) / n_genomes))
  genome_track_idx <- seq_len(n_genomes) + 1L

  gc_pal <- grDevices::colorRampPalette(c("#E8F5E9", "#2E7D32"))(100)
  size_pal <- grDevices::colorRampPalette(c("#E3F2FD", "#1565C0"))(100)
  cds_pal <- grDevices::colorRampPalette(c("#FFF3E0", "#EF6C00"))(100)
  gc_limits <- range_or_default(gc_vals)
  size_limits <- range_or_default(size_vals)
  cds_limits <- range_or_default(cds_vals)

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(output_file, width = 16, height = 14)
  } else {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
  }

  on.exit(try(circlize::circos.clear(), silent = TRUE), add = TRUE)
  if (!is.null(output_file)) {
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  }

  graphics::par(mar = c(1.6, 1.2, 3.2, 1.2), xpd = NA)

  circlize::circos.clear()
  circlize::circos.par(
    start.degree = 90,
    clock.wise = FALSE,
    gap.after = 90,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0, 0),
    points.overflow.warning = FALSE,
    canvas.xlim = c(-1.08, 1.85),
    canvas.ylim = c(-1.08, 1.12)
  )
  circlize::circos.initialize("pan", xlim = c(0, n_clusters))

  cluster_cols <- ifelse(
    n_per$ng >= n_genomes, col_core,
    ifelse(n_per$ng == 1L, col_unique, col_accessory)
  )

  circlize::circos.track(
    "pan",
    ylim = c(0, 1),
    bg.border = NA,
    track.height = cluster_track_h,
    panel.fun = function(x, y) {
      circlize::circos.rect(
        xleft = seq_len(n_clusters) - 1,
        ybottom = 0,
        xright = seq_len(n_clusters),
        ytop = 1,
        col = cluster_cols,
        border = NA
      )
    }
  )

  track_top_y <- numeric(n_genomes)
  track_bottom_y <- numeric(n_genomes)

  for (g in seq_len(n_genomes)) {
    circlize::circos.track(
      "pan",
      ylim = c(0, 1),
      bg.border = NA,
      track.height = genome_track_h,
      panel.fun = function(x, y) {
        circlize::circos.rect(
          xleft = seq_len(n_clusters) - 1,
          ybottom = 0,
          xright = seq_len(n_clusters),
          ytop = 1,
          col = ifelse(bin_mat[, g] == 1L, col_present, col_absent),
          border = NA
        )
      }
    )

    top_xy <- polar_to_xy(circlize::circlize(
      x = n_clusters,
      y = 1,
      sector.index = "pan",
      track.index = genome_track_idx[g]
    ))
    bottom_xy <- polar_to_xy(circlize::circlize(
      x = n_clusters,
      y = 0,
      sector.index = "pan",
      track.index = genome_track_idx[g]
    ))

    track_top_y[g] <- top_xy[1, "y"]
    track_bottom_y[g] <- bottom_xy[1, "y"]
  }

  usr <- graphics::par("usr")
  table_x0 <- 0.10
  table_x1 <- usr[2] - 0.06
  col_fracs <- c(0.16, 0.18, 0.20, 0.46)
  col_breaks <- c(0, cumsum(col_fracs) / sum(col_fracs))
  col_left <- table_x0 + diff(range(c(0, table_x1 - table_x0))) * col_breaks[-length(col_breaks)]
  col_right <- table_x0 + diff(range(c(0, table_x1 - table_x0))) * col_breaks[-1]

  header_y0 <- max(track_top_y) + 0.018
  header_y1 <- min(usr[4] - 0.02, header_y0 + 0.06)

  value_cex <- if (n_genomes > 40) {
    0.38
  } else if (n_genomes > 28) {
    0.46
  } else {
    0.56
  }
  label_cex <- value_cex
  header_cex <- min(0.62, value_cex + 0.08)
  left_pad <- (table_x1 - table_x0) * 0.012

  graphics::rect(
    xleft = table_x0,
    ybottom = min(track_bottom_y),
    xright = table_x1,
    ytop = max(track_top_y),
    col = "#FAFAFA",
    border = NA
  )

  for (j in seq_along(col_left)) {
    graphics::rect(
      col_left[j],
      header_y0,
      col_right[j],
      header_y1,
      col = "#E9ECEF",
      border = "#FFFFFF",
      lwd = 1
    )
  }

  for (g in seq_len(n_genomes)) {
    y0 <- track_bottom_y[g]
    y1 <- track_top_y[g]
    y_mid <- (y0 + y1) / 2

    graphics::rect(col_left[1], y0, col_right[1], y1,
                   col = palette_value(gc_vals[g], gc_pal, gc_limits), border = "#FFFFFF", lwd = 1)
    graphics::rect(col_left[2], y0, col_right[2], y1,
                   col = palette_value(size_vals[g], size_pal, size_limits), border = "#FFFFFF", lwd = 1)
    graphics::rect(col_left[3], y0, col_right[3], y1,
                   col = palette_value(cds_vals[g], cds_pal, cds_limits), border = "#FFFFFF", lwd = 1)
    graphics::rect(col_left[4], y0, col_right[4], y1,
                   col = "#FFFFFF", border = "#FFFFFF", lwd = 1)

    graphics::text((col_left[1] + col_right[1]) / 2, y_mid, fmt_num(gc_vals[g], digits = 1), cex = value_cex)
    graphics::text((col_left[2] + col_right[2]) / 2, y_mid, fmt_num(size_vals[g], digits = 2), cex = value_cex)
    graphics::text((col_left[3] + col_right[3]) / 2, y_mid, fmt_int(cds_vals[g]), cex = value_cex)
    graphics::text(col_left[4] + left_pad, y_mid, strain_labels[g], adj = c(0, 0.5), cex = label_cex)
  }

  for (j in seq_along(col_left)) {
    graphics::text(
      x = (col_left[j] + col_right[j]) / 2,
      y = (header_y0 + header_y1) / 2,
      labels = c("GC%", "Mb", "CDS", "Strain")[j],
      cex = header_cex,
      font = 2
    )
  }

  graphics::rect(
    xleft = table_x0,
    ybottom = min(track_bottom_y),
    xright = table_x1,
    ytop = max(track_top_y),
    border = "#DADADA",
    lwd = 1,
    col = NA
  )

  title(
    main = sprintf(
      "Pangenome (%s clusters x %s genomes)",
      format(n_clusters, big.mark = ","),
      n_genomes
    ),
    cex.main = 1.2,
    font.main = 2
  )

  legend(
    "bottomleft",
    legend = c("Core", "Accessory", "Unique", "Present", "Absent"),
    fill = c(col_core, col_accessory, col_unique, col_present, col_absent),
    border = "grey50",
    cex = 0.78,
    bty = "n"
  )

  if (!is.null(output_file)) {
    message("circos_pangenome written to: ", output_file)
  }

  invisible(NULL)
}
