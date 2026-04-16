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
circos_pangenome <- function(dnmb, output_file = NULL) {
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
    grDevices::pdf(output_file, width = 12, height = 12)
  }

  circlize::circos.clear()
  circlize::circos.par(
    gap.degree     = 2,
    cell.padding   = c(0, 0, 0, 0),
    track.margin   = c(0.005, 0.005),
    start.degree   = 90,
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

  # Title
  title(
    main = sprintf(
      "Pangenome presence/absence  (%s clusters x %s genomes)",
      format(n_clusters, big.mark = ","), n_genomes
    ),
    cex.main = 1.2, font.main = 2
  )

  # Legend
  legend(
    "bottomleft",
    legend = c("Core", "Accessory", "Unique", "Present", "Absent"),
    fill   = c("#2C5F7A", "#F2A766", "#D06461", col_present, col_absent),
    border = "grey50",
    cex    = 0.75,
    bty    = "n"
  )

  circlize::circos.clear()

  if (!is.null(output_file)) {
    grDevices::dev.off()
    message("circos_pangenome written to: ", output_file)
  }

  invisible(NULL)
}
