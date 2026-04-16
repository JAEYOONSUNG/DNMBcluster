#' Anvi'o-style circular pangenome + rectangular gap annotations
#'
#' Single-sector circos (270°) for presence/absence rings, with
#' rectangular base-R annotation panels overlaid in the 90° gap
#' (upper-right). Panels are vertically aligned with the circos
#' rings by computing track positions from circlize internals.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level results directory.
#' @param output_file Optional PDF path.
#' @return Invisible NULL.
#' @export
circos_pangenome <- function(dnmb, results_dir = NULL, output_file = NULL) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    warning("circlize not installed — skipping circos_pangenome")
    return(invisible(NULL))
  }

  # --- Data prep --------------------------------------------------
  presence <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(
      dnmb$genome_meta %>% dplyr::select(genome_uid, genome_key),
      by = "genome_uid"
    )

  # Order genomes by core phylogeny if available, else by genome_uid.
  # This ensures tracks match the tree topology so related strains
  # sit next to each other in the circular display.
  tree_path <- if (!is.null(results_dir))
    file.path(results_dir, "dnmb", "processed", "phylo_tree.nwk") else ""

  if (file.exists(tree_path) && requireNamespace("ape", quietly = TRUE)) {
    tree <- ape::read.tree(tree_path)
    tree <- ape::ladderize(tree)
    # Tip labels are genome_keys (before any pretty-name replacement)
    tip_order <- tree$tip.label
    # Only keep tips that are in genome_meta
    all_keys <- dnmb$genome_meta$genome_key
    genome_keys <- tip_order[tip_order %in% all_keys]
    # Append any genomes not in the tree (shouldn't happen but defensive)
    genome_keys <- c(genome_keys, setdiff(all_keys, genome_keys))
  } else {
    genome_keys <- dnmb$genome_meta %>%
      dplyr::arrange(genome_uid) %>% dplyr::pull(genome_key)
  }

  n_per <- presence %>%
    dplyr::count(cluster_id, name = "ng") %>%
    dplyr::arrange(dplyr::desc(ng), cluster_id)
  sorted_cids <- n_per$cluster_id
  n_clusters <- length(sorted_cids)
  n_genomes  <- length(genome_keys)
  n_total    <- n_genomes

  bin_mat <- matrix(0L, nrow = n_clusters, ncol = n_genomes,
                    dimnames = list(sorted_cids, genome_keys))
  for (i in seq_len(nrow(presence))) {
    cid <- as.character(presence$cluster_id[i])
    gk  <- presence$genome_key[i]
    if (cid %in% rownames(bin_mat)) bin_mat[cid, gk] <- 1L
  }

  cat_ng <- n_per$ng
  names(cat_ng) <- as.character(n_per$cluster_id)

  meta <- dnmb$genome_meta[match(genome_keys, dnmb$genome_meta$genome_key), ]
  strain_labels <- ifelse(!is.na(meta$strain) & nzchar(meta$strain),
                          meta$strain,
                          sub("^(\\S+\\s+\\S+).*", "\\1", meta$organism))
  gc_vals   <- meta$gc_percent
  size_vals <- meta$total_length / 1e6
  cds_vals  <- meta$n_cds

  # Present color must differ from the core-category color (#2C5F7A)
  # so the presence/absence rings and the category annotation ring
  # are visually distinguishable.
  col_present <- "#3B8686"  # teal — clearly distinct from core navy (#2C5F7A)
  col_absent  <- "#F5F5F5"
  track_h <- 0.035
  cat_track_h <- 0.025

  # --- PDF --------------------------------------------------------
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(output_file, width = 14, height = 14)
  }
  par(mar = c(1, 1, 2, 1))

  # --- Single-sector circos (270°) --------------------------------
  circlize::circos.clear()
  circlize::circos.par(
    gap.degree   = 90,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.003, 0.003),
    start.degree = 0,
    clock.wise   = TRUE
  )
  circlize::circos.initialize(factors = "pan", xlim = c(0, n_clusters))

  for (g in seq_len(n_genomes)) {
    circlize::circos.track(
      factors = "pan", ylim = c(0, 1), bg.border = NA,
      track.height = track_h,
      panel.fun = function(x, y) {
        for (i in seq_len(n_clusters)) {
          col <- if (bin_mat[i, g] == 1L) col_present else col_absent
          circlize::circos.rect(i - 1, 0, i, 1, col = col, border = NA)
        }
        circlize::circos.text(
          n_clusters + n_clusters * 0.005, 0.5,
          labels = strain_labels[g],
          facing = "inside", niceFacing = TRUE, cex = 0.45,
          adj = c(0, 0.5)
        )
      }
    )
  }

  # Category ring (innermost)
  circlize::circos.track(
    factors = "pan", ylim = c(0, 1), bg.border = NA,
    track.height = cat_track_h,
    panel.fun = function(x, y) {
      for (i in seq_len(n_clusters)) {
        ng <- cat_ng[as.character(sorted_cids[i])]
        col <- if (ng >= n_total) "#2C5F7A"
               else if (ng == 1)  "#D06461"
               else               "#F2A766"
        circlize::circos.rect(i - 1, 0, i, 1, col = col, border = NA)
      }
    }
  )

  title(
    main = sprintf("Pangenome  (%s clusters x %s genomes)",
                   format(n_clusters, big.mark = ","), n_genomes),
    cex.main = 1.2, font.main = 2
  )
  legend("bottomleft",
         legend = c("Core","Accessory","Unique","Present","Absent"),
         fill = c("#2C5F7A","#F2A766","#D06461",col_present,col_absent),
         border = "grey50", cex = 0.7, bty = "n")

  circlize::circos.clear()

  # --- Rectangular gap-area overlays (upper-right) ----------------
  # Rotated 90°: genomes along x-axis (left=outer → right=inner,
  # matching radial order at the 0° gap boundary), annotation rows
  # stacked vertically (GC / Mb / CDS / Strain top→bottom).

  # Compute exact NDC positions matching the circos ring boundaries
  # at the 0° (right-side) gap edge. The circos center is at NDC
  # (0.5, 0.5); its radius in NDC depends on margins but is ~0.45.
  r_ndc <- 0.45
  margin_per_track <- 2 * 0.003  # track.margin = c(0.003, 0.003)
  total_radial <- n_genomes * (track_h + margin_per_track) +
                  (cat_track_h + margin_per_track)

  # At 0° the rings extend rightward from center.
  # Outermost ring outer edge = center + r_ndc (normalized radius 1.0)
  # Innermost ring inner edge = center + r_ndc * (1 - total_radial)
  genome_x_right <- 0.50 + r_ndc * 0.97   # slight inset from outer edge
  genome_x_left  <- 0.50 + r_ndc * (1 - total_radial) * 0.97

  # Annotation rows stacked in the upper part of the gap
  row_top <- 0.96
  n_rows  <- 4
  row_h   <- 0.085
  row_gap <- 0.005

  draw_row <- function(row_idx, values, bar_col, labels, title_text,
                       is_label_only = FALSE) {
    y_top <- row_top - (row_idx - 1) * (row_h + row_gap)
    y_bot <- y_top - row_h

    par(fig = c(genome_x_left, genome_x_right, y_bot, y_top),
        new = TRUE, mar = c(0, 2.5, 0, 0))
    plot(NULL, xlim = c(0, n_genomes), ylim = c(0, 1),
         axes = FALSE, xlab = "", ylab = "")

    if (is_label_only) {
      for (g in seq_len(n_genomes)) {
        rect(g - 1, 0, g, 1, col = "#FAFAFA", border = "#E8E8E8", lwd = 0.3)
        text(g - 0.5, 0.5, labels[g], cex = 0.32, srt = 45)
      }
    } else {
      vmin <- min(values, na.rm = TRUE)
      vmax <- max(values, na.rm = TRUE)
      if (vmax == vmin) vmax <- vmin + 1
      fracs <- (values - vmin) / (vmax - vmin)
      fracs <- pmax(0.05, fracs)

      for (g in seq_len(n_genomes)) {
        rect(g - 1, 0, g, 1, col = "#F5F5F5", border = "#E8E8E8", lwd = 0.3)
        rect(g - 1 + 0.05, 0, g - 0.05, fracs[g],
             col = bar_col[g], border = NA)
        text(g - 0.5, fracs[g] + 0.08, labels[g],
             cex = 0.28, adj = c(0.5, 0))
      }
    }
    mtext(title_text, side = 2, line = 0.3, cex = 0.6, font = 2, las = 2)
  }

  # Row 1: GC%
  gc_pal <- grDevices::colorRampPalette(c("#A8D5A2", "#006837"))(100)
  gc_idx <- pmax(1, pmin(100, round((gc_vals - min(gc_vals, na.rm=TRUE)) /
    max(diff(range(gc_vals, na.rm=TRUE)), 0.1) * 99) + 1))
  draw_row(1, gc_vals, gc_pal[gc_idx], sprintf("%.1f", gc_vals), "GC%")

  # Row 2: Genome size (Mb)
  sz_pal <- grDevices::colorRampPalette(c("#A3C4DC", "#2C5F7A"))(100)
  sz_idx <- pmax(1, pmin(100, round(size_vals / max(size_vals, na.rm=TRUE) * 99) + 1))
  draw_row(2, size_vals, sz_pal[sz_idx], sprintf("%.1f", size_vals), "Mb")

  # Row 3: CDS count
  cds_pal <- grDevices::colorRampPalette(c("#D4C5E2", "#7B68A0"))(100)
  cds_idx <- pmax(1, pmin(100, round(cds_vals / max(cds_vals, na.rm=TRUE) * 99) + 1))
  draw_row(3, cds_vals, cds_pal[cds_idx], format(cds_vals, big.mark=","), "CDS")

  # Row 4: Strain labels
  draw_row(4, NULL, NULL, strain_labels, "Strain", is_label_only = TRUE)

  if (!is.null(output_file)) {
    grDevices::dev.off()
    message("circos_pangenome written to: ", output_file)
  }
  invisible(NULL)
}
