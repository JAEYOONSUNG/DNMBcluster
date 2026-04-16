#' Circular pangenome — 2-sector with flat annotation cells
#'
#' Pangenome sector (75%): presence/absence rings.
#' Info sector (25%): flat colored cells (NO protruding bars) with
#' values as text. Track widths are IDENTICAL because both sectors
#' share the same circos.track() calls.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level results directory.
#' @param output_file Optional PDF path.
#' @export
circos_pangenome <- function(dnmb, results_dir = NULL, output_file = NULL) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    warning("circlize not installed"); return(invisible(NULL))
  }

  presence <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>% dplyr::distinct() %>%
    dplyr::left_join(dnmb$genome_meta %>% dplyr::select(genome_uid, genome_key), by = "genome_uid")

  tree_path <- if (!is.null(results_dir)) file.path(results_dir, "dnmb", "processed", "phylo_tree.nwk") else ""
  if (file.exists(tree_path) && requireNamespace("ape", quietly = TRUE)) {
    tree <- ape::ladderize(ape::read.tree(tree_path))
    genome_keys <- tree$tip.label[tree$tip.label %in% dnmb$genome_meta$genome_key]
    genome_keys <- c(genome_keys, setdiff(dnmb$genome_meta$genome_key, genome_keys))
  } else {
    genome_keys <- dnmb$genome_meta %>% dplyr::arrange(genome_uid) %>% dplyr::pull(genome_key)
  }

  n_per <- presence %>% dplyr::count(cluster_id, name = "ng") %>% dplyr::arrange(dplyr::desc(ng), cluster_id)
  sorted_cids <- n_per$cluster_id; n_clusters <- length(sorted_cids); n_genomes <- length(genome_keys)
  cat_ng <- setNames(n_per$ng, as.character(n_per$cluster_id))

  bin_mat <- matrix(0L, nrow = n_clusters, ncol = n_genomes, dimnames = list(sorted_cids, genome_keys))
  for (i in seq_len(nrow(presence))) {
    cid <- as.character(presence$cluster_id[i]); gk <- presence$genome_key[i]
    if (cid %in% rownames(bin_mat)) bin_mat[cid, gk] <- 1L
  }

  meta <- dnmb$genome_meta[match(genome_keys, dnmb$genome_meta$genome_key), ]
  strain_labels <- ifelse(!is.na(meta$strain) & nzchar(meta$strain), meta$strain,
                          sub("^(\\S+\\s+\\S+).*", "\\1", meta$organism))
  gc_vals <- meta$gc_percent; size_vals <- meta$total_length / 1e6; cds_vals <- meta$n_cds

  col_present <- "#3B8686"; col_absent <- "#F5F5F5"
  track_h <- 0.035; cat_h <- 0.025

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(output_file, width = 16, height = 14)
  }

  circlize::circos.clear()
  sum_range <- n_clusters / 3
  circlize::circos.par(
    cell.padding = c(0,0,0,0), track.margin = c(0.003, 0.003),
    start.degree = 0, clock.wise = TRUE, gap.after = c(3, 87)
  )
  circlize::circos.initialize(
    factors = c("pan", "info"),
    xlim = matrix(c(0, n_clusters, 0, sum_range), nrow = 2, byrow = TRUE)
  )

  # Color palettes for annotations
  gc_range <- range(gc_vals, na.rm = TRUE)
  gc_pal <- grDevices::colorRampPalette(c("#C8E6C9", "#1B5E20"))(100)
  sz_pal <- grDevices::colorRampPalette(c("#BBDEFB", "#1565C0"))(100)
  cd_pal <- grDevices::colorRampPalette(c("#E1BEE7", "#6A1B9A"))(100)

  for (g in seq_len(n_genomes)) {
    circlize::circos.track(
      factors = c("pan", "info"), ylim = c(0, 1), bg.border = NA, track.height = track_h,
      panel.fun = function(x, y) {
        sec <- circlize::CELL_META$sector.index
        if (sec == "pan") {
          for (i in seq_len(n_clusters)) {
            col <- if (bin_mat[i, g] == 1L) col_present else col_absent
            circlize::circos.rect(i-1, 0, i, 1, col = col, border = NA)
          }
        } else {
          w <- sum_range / 4
          # GC% — flat colored cell
          gc_i <- max(1, min(100, round((gc_vals[g] - gc_range[1]) / max(diff(gc_range), 0.1) * 99) + 1))
          circlize::circos.rect(0, 0, w, 1, col = gc_pal[gc_i], border = "white")
          circlize::circos.text(w/2, 0.5, sprintf("%.1f", gc_vals[g]),
                                cex = 0.30, facing = "inside", niceFacing = TRUE)
          # Mb — flat colored cell
          sz_i <- max(1, min(100, round(size_vals[g] / max(size_vals, na.rm=TRUE) * 99) + 1))
          circlize::circos.rect(w, 0, 2*w, 1, col = sz_pal[sz_i], border = "white")
          circlize::circos.text(1.5*w, 0.5, sprintf("%.1f", size_vals[g]),
                                cex = 0.28, facing = "inside", niceFacing = TRUE)
          # CDS — flat colored cell
          cd_i <- max(1, min(100, round(cds_vals[g] / max(cds_vals, na.rm=TRUE) * 99) + 1))
          circlize::circos.rect(2*w, 0, 3*w, 1, col = cd_pal[cd_i], border = "white")
          circlize::circos.text(2.5*w, 0.5, format(cds_vals[g], big.mark=","),
                                cex = 0.24, facing = "inside", niceFacing = TRUE)
          # Strain label — plain
          circlize::circos.rect(3*w, 0, 4*w, 1, col = "#FAFAFA", border = "white")
          circlize::circos.text(3.5*w, 0.5, strain_labels[g],
                                cex = 0.35, facing = "inside", niceFacing = TRUE)
        }
      }
    )
  }

  # Category + headers
  circlize::circos.track(
    factors = c("pan", "info"), ylim = c(0,1), bg.border = NA, track.height = cat_h,
    panel.fun = function(x, y) {
      sec <- circlize::CELL_META$sector.index
      if (sec == "pan") {
        for (i in seq_len(n_clusters)) {
          ng <- cat_ng[as.character(sorted_cids[i])]
          col <- if (ng >= n_genomes) "#2C5F7A" else if (ng == 1) "#D06461" else "#F2A766"
          circlize::circos.rect(i-1, 0, i, 1, col = col, border = NA)
        }
      } else {
        w <- sum_range / 4
        hdrs <- c("GC%", "Mb", "CDS", "Strain")
        for (j in seq_along(hdrs))
          circlize::circos.text((j-0.5)*w, 0.5, hdrs[j], cex = 0.40, font = 2,
                                facing = "inside", niceFacing = TRUE)
      }
    }
  )

  title(main = sprintf("Pangenome  (%s clusters x %s genomes)",
                       format(n_clusters, big.mark=","), n_genomes),
        cex.main = 1.2, font.main = 2)
  legend("bottomleft",
         legend = c("Core","Accessory","Unique","Present","Absent"),
         fill = c("#2C5F7A","#F2A766","#D06461", col_present, col_absent),
         border = "grey50", cex = 0.7, bty = "n")

  circlize::circos.clear()
  if (!is.null(output_file)) { grDevices::dev.off(); message("circos_pangenome written to: ", output_file) }
  invisible(NULL)
}
