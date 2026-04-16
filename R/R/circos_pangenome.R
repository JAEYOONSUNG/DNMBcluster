#' Anvi'o-style circular pangenome — rings + gap-area annotations
#'
#' Two-sector circlize layout:
#' - **Pangenome sector (270°)**: one concentric ring per genome
#'   showing presence/absence. Sorted by n_genomes (core → unique).
#' - **Summary sector (90°)**: per-genome annotation columns that
#'   are PERFECTLY aligned with the rings because they share the
#'   same radial track structure. Columns: GC%, genome size (Mb),
#'   gene count, and ANI mini-heatmap.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level results directory (for ANI parquet).
#' @param output_file Optional PDF path.
#' @return Invisible NULL.
#' @export
circos_pangenome <- function(dnmb, results_dir = NULL, output_file = NULL) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    warning("circlize not installed — skipping circos_pangenome")
    return(invisible(NULL))
  }

  # --- Build presence/absence binary matrix ----------------------
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

  n_per_cluster <- presence %>%
    dplyr::count(cluster_id, name = "n_genomes") %>%
    dplyr::arrange(dplyr::desc(n_genomes), cluster_id)
  sorted_cids <- n_per_cluster$cluster_id

  n_clusters <- length(sorted_cids)
  n_genomes  <- length(genome_keys)

  bin_mat <- matrix(0L, nrow = n_clusters, ncol = n_genomes,
                    dimnames = list(sorted_cids, genome_keys))
  for (i in seq_len(nrow(presence))) {
    cid <- as.character(presence$cluster_id[i])
    gk  <- presence$genome_key[i]
    if (cid %in% rownames(bin_mat)) bin_mat[cid, gk] <- 1L
  }

  # Category for each cluster
  n_total <- n_genomes
  cat_lookup <- n_per_cluster$n_genomes
  names(cat_lookup) <- as.character(n_per_cluster$cluster_id)

  # Metadata
  meta <- dnmb$genome_meta %>% dplyr::arrange(genome_uid)
  strain_labels <- ifelse(
    !is.na(meta$strain) & nzchar(meta$strain), meta$strain,
    sub("^(\\S+\\s+\\S+).*", "\\1", meta$organism)
  )
  gc_vals   <- meta$gc_percent
  size_vals <- meta$total_length / 1e6
  cds_vals  <- meta$n_cds

  col_present <- "#2C5F7A"
  col_absent  <- "#F0F0F0"

  # --- PDF -------------------------------------------------------
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(output_file, width = 14, height = 14)
  }

  # --- Two-sector circos -----------------------------------------
  circlize::circos.clear()
  circlize::circos.par(
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.005, 0.005),
    start.degree = 0,
    clock.wise   = TRUE,
    gap.after    = c(5, 85)  # 5° gap after pangenome, 85° gap after summary = 90° total gap for summary
  )

  # The summary sector needs ~25% of the arc so its boxes are wide
  # enough to read. circlize sizes sectors proportional to their
  # xlim ranges, so we inflate the summary range to n_clusters / 3.
  summary_range <- n_clusters / 3
  circlize::circos.initialize(
    factors = c("pangenome", "summary"),
    xlim    = matrix(c(0, n_clusters, 0, summary_range), nrow = 2, byrow = TRUE)
  )

  # --- Per-genome tracks (one per genome) -------------------------
  for (g in seq_len(n_genomes)) {
    gk <- genome_keys[g]

    circlize::circos.track(
      factors = c("pangenome", "summary"),
      ylim    = c(0, 1),
      bg.border = NA,
      track.height = 0.035,
      panel.fun = function(x, y) {
        sector <- circlize::CELL_META$sector.index

        if (sector == "pangenome") {
          # Presence/absence rectangles
          for (i in seq_len(n_clusters)) {
            col <- if (bin_mat[i, g] == 1L) col_present else col_absent
            circlize::circos.rect(i - 1, 0, i, 1, col = col, border = NA)
          }
        } else {
          # Summary sector: x range is 0..summary_range. Divide into
          # 4 equal sub-columns scaled to that range.
          sr <- summary_range
          w  <- sr / 4

          # GC% color box
          gc_range <- range(gc_vals, na.rm = TRUE)
          gc_pal <- grDevices::colorRampPalette(c("#FFFFCC", "#006837"))(100)
          gc_i <- max(1, min(100, round((gc_vals[g] - gc_range[1]) / max(diff(gc_range), 0.1) * 99) + 1))
          circlize::circos.rect(0, 0, w, 1, col = gc_pal[gc_i], border = "white")
          circlize::circos.text(w / 2, 0.5, sprintf("%.1f", gc_vals[g]),
                                cex = 0.35, facing = "inside", niceFacing = TRUE)

          # Genome size bar
          h <- size_vals[g] / max(size_vals, na.rm = TRUE)
          circlize::circos.rect(w, 0, 2 * w, h, col = "#5E81AC", border = "white")
          circlize::circos.text(1.5 * w, 0.5, sprintf("%.1f", size_vals[g]),
                                cex = 0.30, facing = "inside", niceFacing = TRUE)

          # CDS count bar
          h2 <- cds_vals[g] / max(cds_vals, na.rm = TRUE)
          circlize::circos.rect(2 * w, 0, 3 * w, h2, col = "#B48EAD", border = "white")
          circlize::circos.text(2.5 * w, 0.5, format(cds_vals[g], big.mark = ","),
                                cex = 0.25, facing = "inside", niceFacing = TRUE)

          # Strain label
          circlize::circos.text(3.5 * w, 0.5, strain_labels[g],
                                cex = 0.40, facing = "inside", niceFacing = TRUE)
        }
      }
    )
  }

  # Category annotation track (innermost)
  circlize::circos.track(
    factors = c("pangenome", "summary"),
    ylim = c(0, 1),
    bg.border = NA,
    track.height = 0.025,
    panel.fun = function(x, y) {
      sector <- circlize::CELL_META$sector.index
      if (sector == "pangenome") {
        for (i in seq_len(n_clusters)) {
          cid <- as.character(sorted_cids[i])
          ng <- cat_lookup[cid]
          col <- if (ng >= n_total) "#2C5F7A"
                 else if (ng == 1)  "#D06461"
                 else               "#F2A766"
          circlize::circos.rect(i - 1, 0, i, 1, col = col, border = NA)
        }
      } else {
        # Summary header labels
        w <- summary_range / 4
        labels <- c("GC%", "Mb", "CDS", "Strain")
        for (j in seq_along(labels)) {
          circlize::circos.text((j - 0.5) * w, 0.5, labels[j],
                                cex = 0.45, facing = "inside", niceFacing = TRUE,
                                font = 2)
        }
      }
    }
  )

  # Title + legend
  title(
    main = sprintf("Pangenome  (%s clusters x %s genomes)",
                   format(n_clusters, big.mark = ","), n_genomes),
    cex.main = 1.2, font.main = 2
  )
  legend(
    "bottomleft",
    legend = c("Core", "Accessory", "Unique", "Present", "Absent"),
    fill   = c("#2C5F7A", "#F2A766", "#D06461", col_present, col_absent),
    border = "grey50", cex = 0.7, bty = "n"
  )

  circlize::circos.clear()

  if (!is.null(output_file)) {
    grDevices::dev.off()
    message("circos_pangenome written to: ", output_file)
  }

  invisible(NULL)
}
