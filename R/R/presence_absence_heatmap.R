#' Cluster presence/absence heatmap via ComplexHeatmap
#'
#' Publication-grade binary matrix: rows = clusters (sorted by
#' n_genomes, core at top), columns = genomes (hierarchically
#' clustered by Jaccard). Each cell is dark (present) or light
#' (absent). Side annotations show cluster category (core/accessory/
#' unique); top annotations show genus + species coloring.
#'
#' For large pan-genomes (> 5k clusters) the plot compresses
#' gracefully — individual rows are not labeled but the overall
#' pattern (core plateau → accessory gradient → unique tail) is
#' immediately visible.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return A ComplexHeatmap object, or a ggplot if ComplexHeatmap is unavailable.
#' @export
presence_absence_heatmap <- function(dnmb, output_file = NULL) {

  # --- Build binary matrix ----------------------------------------
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

  # Cluster order: core first (most genomes) → unique last
  n_per <- presence %>%
    dplyr::count(cluster_id, name = "ng") %>%
    dplyr::arrange(dplyr::desc(ng), cluster_id)
  sorted_cids <- n_per$cluster_id

  bin_mat <- matrix(0L, nrow = length(sorted_cids), ncol = length(genome_keys),
                    dimnames = list(sorted_cids, genome_keys))
  for (i in seq_len(nrow(presence))) {
    cid <- as.character(presence$cluster_id[i])
    gk  <- presence$genome_key[i]
    if (cid %in% rownames(bin_mat)) bin_mat[cid, gk] <- 1L
  }

  n_total <- length(genome_keys)
  cat_vec <- ifelse(n_per$ng >= n_total, "core",
                    ifelse(n_per$ng == 1, "unique", "accessory"))

  # --- Try ComplexHeatmap, fall back to ggplot --------------------
  if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
      requireNamespace("circlize", quietly = TRUE)) {

    # Row annotation: category
    cat_colors <- .dnmb_presence_pal()
    row_ha <- ComplexHeatmap::rowAnnotation(
      Category = cat_vec,
      col = list(Category = cat_colors),
      show_annotation_name = FALSE,
      width = grid::unit(4, "mm")
    )

    # Column annotation: genus + species from genome_meta
    meta <- dnmb$genome_meta %>% dplyr::arrange(genome_uid)
    genus <- sub("^(\\S+).*", "\\1", meta$organism)
    genus <- ifelse(is.na(genus), "Unknown", genus)
    species <- sub("^\\S+\\s+(\\S+).*", "\\1", meta$organism)
    species <- ifelse(is.na(species) | species == genus, "sp.", species)

    unique_genus <- sort(unique(genus))
    unique_species <- sort(unique(species))
    g_cols <- stats::setNames(
      grDevices::hcl.colors(length(unique_genus), "Dark 3"),
      unique_genus
    )
    s_cols <- stats::setNames(
      grDevices::hcl.colors(length(unique_species), "Set 3"),
      unique_species
    )

    col_ha <- ComplexHeatmap::HeatmapAnnotation(
      Genus   = genus,
      Species = species,
      col = list(Genus = g_cols, Species = s_cols),
      annotation_name_side = "left",
      simple_anno_size = grid::unit(3, "mm")
    )

    # Heatmap color
    col_fun <- circlize::colorRamp2(c(0, 1), c("#F5F5F5", "#2C5F7A"))

    ht <- ComplexHeatmap::Heatmap(
      bin_mat,
      name = "Present",
      col = col_fun,
      top_annotation = col_ha,
      left_annotation = row_ha,
      cluster_rows = FALSE,  # already sorted by n_genomes
      cluster_columns = TRUE,
      clustering_distance_columns = "binary",
      clustering_method_columns = "average",
      show_row_names = FALSE,
      show_column_names = TRUE,
      column_names_gp = grid::gpar(fontsize = 7),
      column_names_rot = 45,
      border = TRUE,
      rect_gp = grid::gpar(col = NA),
      use_raster = nrow(bin_mat) > 2000,
      column_title = sprintf(
        "Presence/absence matrix  (%s clusters x %s genomes)",
        format(nrow(bin_mat), big.mark = ","), ncol(bin_mat)
      ),
      column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
      heatmap_legend_param = list(
        title = "Gene",
        at = c(0, 1),
        labels = c("Absent", "Present")
      )
    )

    if (!is.null(output_file)) {
      dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
      n <- ncol(bin_mat)
      fig_w <- max(8, 4 + n * 0.4)
      fig_h <- max(8, min(16, nrow(bin_mat) * 0.003 + 3))
      grDevices::pdf(output_file, width = fig_w, height = fig_h)
      ComplexHeatmap::draw(ht, merge_legend = TRUE)
      grDevices::dev.off()
      message("presence_absence_heatmap written to: ", output_file)
    }
    return(invisible(ht))

  } else {
    # Fallback: simple ggplot geom_tile (original implementation)
    long <- presence %>%
      dplyr::mutate(
        cluster_id = factor(cluster_id, levels = sorted_cids),
        present = 1L
      )
    cat_df <- data.frame(
      cluster_id = sorted_cids,
      category = factor(cat_vec, levels = c("core","accessory","unique"))
    )
    long <- long %>% dplyr::left_join(cat_df, by = "cluster_id")

    p <- ggplot2::ggplot(long, ggplot2::aes(x = genome_key, y = cluster_id)) +
      ggplot2::geom_tile(ggplot2::aes(fill = category)) +
      ggplot2::scale_fill_manual(values = .dnmb_presence_pal()) +
      ggplot2::scale_y_discrete(breaks = NULL) +
      ggplot2::labs(title = "Presence/absence matrix", x = NULL, y = "Cluster", fill = NULL) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank(),
        legend.position = "bottom"
      )

    if (!is.null(output_file)) {
      dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
      n <- length(unique(long$genome_key))
      fig_w <- min(max(3, 1.5 + 0.35 * n), 8)
      ggplot2::ggsave(output_file, p, width = fig_w, height = 10, dpi = 200)
      message("presence_absence_heatmap written to: ", output_file)
    }
    return(invisible(p))
  }
}
