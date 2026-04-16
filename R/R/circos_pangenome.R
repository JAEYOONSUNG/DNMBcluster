#' Anvi'o-style circular pangenome — rings + gap annotations
#'
#' Two-sector circlize: pangenome (75%) + summary (25%). All
#' annotations share the same radial tracks as the genome rings
#' so alignment is GUARANTEED. Summary cells show vertical bars
#' for GC%, genome size, CDS count + strain labels.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level results directory.
#' @param output_file Optional PDF path.
#' @export
circos_pangenome <- function(dnmb, results_dir = NULL, output_file = NULL) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    warning("circlize not installed — skipping"); return(invisible(NULL))
  }

  presence <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>% dplyr::distinct() %>%
    dplyr::left_join(dnmb$genome_meta %>% dplyr::select(genome_uid, genome_key), by = "genome_uid")

  # Phylo-ordered genomes if tree available
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
    grDevices::pdf(output_file, width = 14, height = 14)
  }

  circlize::circos.clear()
  # Summary sector: inflated xlim so it takes ~25% of the arc
  sum_range <- n_clusters / 3
  circlize::circos.par(
    cell.padding = c(0, 0, 0, 0), track.margin = c(0.003, 0.003),
    start.degree = 0, clock.wise = TRUE, gap.after = c(3, 87)
  )
  circlize::circos.initialize(
    factors = c("pan", "info"),
    xlim = matrix(c(0, n_clusters, 0, sum_range), nrow = 2, byrow = TRUE)
  )

  for (g in seq_len(n_genomes)) {
    circlize::circos.track(
      factors = c("pan", "info"), ylim = c(0, 1), bg.border = NA, track.height = track_h,
      panel.fun = function(x, y) {
        sec <- circlize::CELL_META$sector.index
        if (sec == "pan") {
          for (i in seq_len(n_clusters)) {
            col <- if (bin_mat[i, g] == 1L) col_present else col_absent
            circlize::circos.rect(i - 1, 0, i, 1, col = col, border = NA)
          }
        } else {
          # 4 sub-columns in the info sector
          w <- sum_range / 4
          # GC%
          gc_r <- range(gc_vals, na.rm = TRUE)
          gc_pal <- grDevices::colorRampPalette(c("#C8E6C9", "#1B5E20"))(100)
          gc_i <- max(1, min(100, round((gc_vals[g] - gc_r[1]) / max(diff(gc_r), 0.1) * 99) + 1))
          frac_gc <- (gc_vals[g] - gc_r[1]) / max(diff(gc_r), 0.1)
          # Bars grow INWARD (from outer edge y=1 toward center y=0)
          # so they visually hang down from the track's outer boundary.
          circlize::circos.rect(0, 0, w, 1, col = "#F5F5F5", border = "white")
          circlize::circos.rect(0, 1 - max(0.05, frac_gc), w, 1, col = gc_pal[gc_i], border = NA)
          circlize::circos.text(w / 2, 0.5, sprintf("%.1f", gc_vals[g]), cex = 0.30,
                                facing = "inside", niceFacing = TRUE)
          # Size (Mb)
          frac_s <- size_vals[g] / max(size_vals, na.rm = TRUE)
          circlize::circos.rect(w, 0, 2*w, 1, col = "#F5F5F5", border = "white")
          circlize::circos.rect(w, 1 - max(0.05, frac_s), 2*w, 1, col = "#5E81AC", border = NA)
          circlize::circos.text(1.5*w, 0.5, sprintf("%.1f", size_vals[g]), cex = 0.28,
                                facing = "inside", niceFacing = TRUE)
          # CDS
          frac_c <- cds_vals[g] / max(cds_vals, na.rm = TRUE)
          circlize::circos.rect(2*w, 0, 3*w, 1, col = "#F5F5F5", border = "white")
          circlize::circos.rect(2*w, 1 - max(0.05, frac_c), 3*w, 1, col = "#B48EAD", border = NA)
          circlize::circos.text(2.5*w, 0.5, format(cds_vals[g], big.mark = ","), cex = 0.24,
                                facing = "inside", niceFacing = TRUE)
          # Strain
          circlize::circos.text(3.5*w, 0.5, strain_labels[g], cex = 0.35,
                                facing = "inside", niceFacing = TRUE)
        }
      }
    )
  }

  # Category + header track
  circlize::circos.track(
    factors = c("pan", "info"), ylim = c(0, 1), bg.border = NA, track.height = cat_h,
    panel.fun = function(x, y) {
      sec <- circlize::CELL_META$sector.index
      if (sec == "pan") {
        for (i in seq_len(n_clusters)) {
          ng <- cat_ng[as.character(sorted_cids[i])]
          col <- if (ng >= n_genomes) "#2C5F7A" else if (ng == 1) "#D06461" else "#F2A766"
          circlize::circos.rect(i - 1, 0, i, 1, col = col, border = NA)
        }
      } else {
        w <- sum_range / 4
        for (j in seq_along(c("GC%", "Mb", "CDS", "Strain")))
          circlize::circos.text((j - 0.5) * w, 0.5, c("GC%", "Mb", "CDS", "Strain")[j],
                                cex = 0.4, font = 2, facing = "inside", niceFacing = TRUE)
      }
    }
  )

  title(main = sprintf("Pangenome  (%s clusters x %s genomes)",
                       format(n_clusters, big.mark = ","), n_genomes),
        cex.main = 1.2, font.main = 2)
  legend("bottomleft",
         legend = c("Core","Accessory","Unique","Present","Absent"),
         fill = c("#2C5F7A","#F2A766","#D06461",col_present,col_absent),
         border = "grey50", cex = 0.7, bty = "n")

  circlize::circos.clear()
  if (!is.null(output_file)) { grDevices::dev.off(); message("circos_pangenome written to: ", output_file) }
  invisible(NULL)
}
