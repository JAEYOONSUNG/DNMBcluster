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

  col_present <- "#2C5F7A"
  col_absent  <- "#F0F0F0"
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
  # Compute NDC y-positions matching each track ring.
  # circos center = NDC (0.5, 0.5). The plot radius in NDC ≈ 0.45.
  # Outermost track starts at radius ~0.97 (in circos units), which
  # maps to NDC 0.5 + 0.45*0.97 ≈ 0.94.
  # Each track is track_h=0.035 of the unit radius.
  # In NDC: track_h_ndc = 0.45 * 0.035 ≈ 0.016.
  # Track margins eat ~0.006 per track, so effective ndc per track ≈ 0.018.

  r_ndc <- 0.44  # approximate plot radius in NDC
  track_ndc <- r_ndc * (track_h + 0.006)  # include margins

  # Gap quadrant top: outermost ring position in NDC
  # Since gap is at 0°–90° (upper-right), tracks align along the
  # y-axis at x ≈ 0.5 (the right edge of the circle). The y-position
  # of each genome is: top of outermost ring → bottom of innermost.
  outer_y <- 0.50 + r_ndc * 0.96  # top of genome 1 ring
  n_cols <- 4  # GC, Size, CDS, label

  # Annotation panel: 4 columns side by side, each containing
  # n_genomes rows aligned with the tracks.
  panel_left  <- 0.56
  panel_right <- 0.96
  col_w <- (panel_right - panel_left) / n_cols

  # Draw one annotation column as horizontal bar graphs.
  # Each row = one genome, bar length ∝ value within [vmin, vmax].
  # Bar fills from left, colored by gradient; number label at bar tip.
  draw_col <- function(col_idx, values, bar_col, labels, title_text,
                       is_label_only = FALSE) {
    x1 <- panel_left + (col_idx - 1) * col_w
    x2 <- x1 + col_w
    top_y <- outer_y
    bot_y <- top_y - n_genomes * track_ndc

    par(fig = c(x1, x2, bot_y, top_y), new = TRUE, mar = c(0, 0, 1.2, 0))
    plot(NULL, xlim = c(0, 1), ylim = c(0, n_genomes),
         axes = FALSE, xlab = "", ylab = "")

    if (is_label_only) {
      for (g in seq_len(n_genomes)) {
        y_bot <- n_genomes - g
        rect(0, y_bot, 1, y_bot + 1, col = "#FAFAFA", border = "#E8E8E8", lwd = 0.3)
        text(0.5, y_bot + 0.5, labels[g], cex = 0.38)
      }
    } else {
      # Normalize values to [0, 1] for bar length
      vmin <- min(values, na.rm = TRUE)
      vmax <- max(values, na.rm = TRUE)
      if (vmax == vmin) vmax <- vmin + 1
      fracs <- (values - vmin) / (vmax - vmin)
      fracs <- pmax(0.05, fracs)  # minimum 5% so the bar is visible

      for (g in seq_len(n_genomes)) {
        y_bot <- n_genomes - g
        # Background
        rect(0, y_bot, 1, y_bot + 1, col = "#F5F5F5", border = "#E8E8E8", lwd = 0.3)
        # Proportional bar
        rect(0, y_bot + 0.05, fracs[g], y_bot + 0.95,
             col = bar_col[g], border = NA)
        # Value label
        text(fracs[g] + 0.03, y_bot + 0.5, labels[g],
             cex = 0.35, adj = c(0, 0.5))
      }
    }
    title(main = title_text, cex.main = 0.7, font.main = 2, line = 0.1)
  }

  # GC% — green gradient bars
  gc_pal <- grDevices::colorRampPalette(c("#A8D5A2", "#006837"))(100)
  gc_idx <- pmax(1, pmin(100, round((gc_vals - min(gc_vals, na.rm=TRUE)) /
    max(diff(range(gc_vals, na.rm=TRUE)), 0.1) * 99) + 1))
  draw_col(1, gc_vals, gc_pal[gc_idx],
           sprintf("%.1f", gc_vals), "GC%")

  # Genome size (Mb) — blue gradient bars
  sz_pal <- grDevices::colorRampPalette(c("#A3C4DC", "#2C5F7A"))(100)
  sz_idx <- pmax(1, pmin(100, round(size_vals / max(size_vals, na.rm=TRUE) * 99) + 1))
  draw_col(2, size_vals, sz_pal[sz_idx],
           sprintf("%.1f", size_vals), "Mb")

  # CDS count — purple gradient bars
  cds_pal <- grDevices::colorRampPalette(c("#D4C5E2", "#7B68A0"))(100)
  cds_idx <- pmax(1, pmin(100, round(cds_vals / max(cds_vals, na.rm=TRUE) * 99) + 1))
  draw_col(3, cds_vals, cds_pal[cds_idx],
           format(cds_vals, big.mark = ","), "CDS")

  # Strain label — text only
  draw_col(4, NULL, NULL, strain_labels, "Strain", is_label_only = TRUE)

  if (!is.null(output_file)) {
    grDevices::dev.off()
    message("circos_pangenome written to: ", output_file)
  }
  invisible(NULL)
}
