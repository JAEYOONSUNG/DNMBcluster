#' Circular pangenome + annotation table (single page)
#'
#' Top panel: single-sector circos (270°) with presence/absence rings.
#' Bottom panel: ggplot bar chart table showing GC%, genome size, CDS
#' count per genome — same genome order, separate coordinate system,
#' no alignment headaches.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level results directory.
#' @param output_file Optional PDF path.
#' @export
circos_pangenome <- function(dnmb, results_dir = NULL, output_file = NULL) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    warning("circlize not installed"); return(invisible(NULL))
  }

  # --- Data -------------------------------------------------------
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

  # Singleton counts per genome
  singleton_cids <- n_per$cluster_id[n_per$ng == 1]
  singleton_counts <- vapply(genome_keys, function(gk) {
    sum(bin_mat[as.character(singleton_cids), gk])
  }, integer(1))

  col_present <- "#3B8686"; col_absent <- "#F5F5F5"
  track_h <- 0.035; cat_h <- 0.025

  # --- PDF --------------------------------------------------------
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(output_file, width = 14, height = 18)
  }

  # === LEFT: Circos (270°) | RIGHT: Annotation table ===============
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1.2))
  par(mar = c(1, 1, 2, 0))

  circlize::circos.clear()
  circlize::circos.par(
    gap.degree = 90, cell.padding = c(0,0,0,0),
    track.margin = c(0.003, 0.003),
    start.degree = 0, clock.wise = TRUE
  )
  circlize::circos.initialize(factors = "pan", xlim = c(0, n_clusters))

  for (g in seq_len(n_genomes)) {
    circlize::circos.track(
      factors = "pan", ylim = c(0,1), bg.border = NA, track.height = track_h,
      panel.fun = function(x, y) {
        for (i in seq_len(n_clusters)) {
          col <- if (bin_mat[i, g] == 1L) col_present else col_absent
          circlize::circos.rect(i-1, 0, i, 1, col = col, border = NA)
        }
        circlize::circos.text(n_clusters * 1.005, 0.5, strain_labels[g],
                              facing = "inside", niceFacing = TRUE, cex = 0.45, adj = c(0, 0.5))
      }
    )
  }
  circlize::circos.track(
    factors = "pan", ylim = c(0,1), bg.border = NA, track.height = cat_h,
    panel.fun = function(x, y) {
      for (i in seq_len(n_clusters)) {
        ng <- cat_ng[as.character(sorted_cids[i])]
        col <- if (ng >= n_genomes) "#2C5F7A" else if (ng == 1) "#D06461" else "#F2A766"
        circlize::circos.rect(i-1, 0, i, 1, col = col, border = NA)
      }
    }
  )

  title(main = sprintf("Pangenome  (%s clusters x %s genomes)",
                       format(n_clusters, big.mark = ","), n_genomes),
        cex.main = 1.3, font.main = 2)
  legend("bottomleft",
         legend = c("Core","Accessory","Unique","Present","Absent"),
         fill = c("#2C5F7A","#F2A766","#D06461", col_present, col_absent),
         border = "grey50", cex = 0.7, bty = "n")
  circlize::circos.clear()

  # === RIGHT PANEL: Annotation table (genomes as rows) ============
  par(mar = c(1, 0, 2, 1))

  anno_mat <- cbind(
    GC     = gc_vals,
    Mb     = size_vals,
    CDS    = cds_vals / 1000,
    Unique = singleton_counts
  )
  rownames(anno_mat) <- strain_labels

  anno_norm <- apply(anno_mat, 2, function(c) {
    rng <- range(c, na.rm = TRUE)
    if (diff(rng) == 0) rep(0.5, length(c)) else (c - rng[1]) / diff(rng)
  })

  metric_cols <- c(GC = "#4CAF50", Mb = "#5E81AC", CDS = "#B48EAD", Unique = "#D06461")
  n_metrics <- ncol(anno_mat)

  plot(NULL, xlim = c(0, n_metrics + 0.5), ylim = c(0, n_genomes),
       axes = FALSE, xlab = "", ylab = "")

  for (g in seq_len(n_genomes)) {
    y_base <- n_genomes - g  # genome 1 at top
    for (m in seq_len(n_metrics)) {
      x_base <- m - 1
      rect(x_base, y_base, x_base + 0.95, y_base + 0.95,
           col = "#F8F8F8", border = "#E8E8E8", lwd = 0.3)
      bar_w <- anno_norm[g, m] * 0.85 + 0.05
      rect(x_base + 0.02, y_base + 0.02, x_base + 0.02 + bar_w * 0.91, y_base + 0.93,
           col = metric_cols[m], border = NA)

      val_txt <- switch(colnames(anno_mat)[m],
                        GC = sprintf("%.1f", anno_mat[g, m]),
                        Mb = sprintf("%.1f", anno_mat[g, m]),
                        CDS = format(as.integer(anno_mat[g, m] * 1000), big.mark = ","),
                        Unique = format(as.integer(anno_mat[g, m]), big.mark = ","))
      text(x_base + 0.48, y_base + 0.48, val_txt, cex = 0.33)
    }
    # Strain label on right
    text(n_metrics + 0.15, y_base + 0.48, strain_labels[g], cex = 0.38, adj = c(0, 0.5))
  }
  # Column headers at top
  for (m in seq_len(n_metrics)) {
    text(m - 0.52, n_genomes + 0.1, colnames(anno_mat)[m],
         cex = 0.6, font = 2, adj = c(0.5, 0))
  }

  if (!is.null(output_file)) {
    grDevices::dev.off()
    message("circos_pangenome written to: ", output_file)
  }
  invisible(NULL)
}
