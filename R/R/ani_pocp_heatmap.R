#' ANI + POCP genome-genome heatmaps via ComplexHeatmap
#'
#' Publication-grade heatmaps with:
#' - Genus + Species annotation color strips (auto-extracted from
#'   ``genome_meta.organism``)
#' - Hierarchical clustering with dendrograms (euclidean / complete)
#' - Italic ``Genus species`` + plain ``strain`` row/col labels
#' - POCP gradient 60-100, ANI gradient 80-100
#' - ``colorspace::diverge_hcl("Vik")`` palette (matching the user's
#'   reference Geobacillus POCP figure)
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level pipeline results directory.
#' @param output_file_ani PDF path for ANI heatmap.
#' @param output_file_pocp PDF path for POCP heatmap.
#' @return Named list with the two ComplexHeatmap objects.
#' @export
ani_pocp_heatmap <- function(dnmb, results_dir,
                             output_file_ani = NULL,
                             output_file_pocp = NULL) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE) ||
      !requireNamespace("colorspace", quietly = TRUE)) {
    warning("ComplexHeatmap/circlize/colorspace not installed -- skipping")
    return(invisible(NULL))
  }

  processed_dir <- file.path(results_dir, "dnmb", "processed")
  result <- list(ani = NULL, pocp = NULL)

  # --- Helper: build square matrix from long-form parquet ---------
  load_matrix <- function(parquet_name, value_col) {
    path <- file.path(processed_dir, parquet_name)
    if (!file.exists(path)) return(NULL)
    df <- tibble::as_tibble(arrow::read_parquet(path))
    keys <- sort(unique(c(df$genome_a, df$genome_b)))
    n <- length(keys)
    mat <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(keys, keys))
    mat[cbind(df$genome_a, df$genome_b)] <- df[[value_col]]
    mat
  }

  # --- Helper: extract Genus / Species from genome_meta -----------
  build_annotation <- function(mat_keys) {
    meta <- dnmb$genome_meta %>%
      dplyr::filter(genome_key %in% mat_keys)

    anno_df <- data.frame(row.names = meta$genome_key)

    # Extract genus (first word of organism) and species (second word)
    org <- meta$organism
    genus   <- ifelse(is.na(org), "Unknown",
                      sub("^(\\S+).*", "\\1", org))
    species <- ifelse(is.na(org), "sp.",
                      sub("^\\S+\\s+(\\S+).*", "\\1", org))
    species <- ifelse(species == genus, "sp.", species)

    anno_df$Genus   <- genus
    anno_df$Species <- species
    anno_df <- anno_df[mat_keys, , drop = FALSE]

    # Auto-generate colors for each unique genus and species
    unique_genus <- sort(unique(anno_df$Genus))
    unique_species <- sort(unique(anno_df$Species))

    lighten <- function(cols, amount = 0.20) {
      m <- grDevices::col2rgb(cols) / 255
      m <- m + (1 - m) * amount
      grDevices::rgb(m[1, ], m[2, ], m[3, ])
    }
    genus_base <- c("#E66101", "#FDB863", "#5E3C99", "#018571",
                    "#CA0020", "#0571B0", "#404040", "#B2ABD2")
    species_base <- c("#CA0020", "#F4A582", "#4DAC26", "#B8E186",
                      "#0571B0", "#92C5DE", "#B2ABD2", "#5E3C99",
                      "#80CDC1", "#018571", "#DFC27D", "#545863",
                      "#BABABA", "#404040", "#E66101", "#FDB863")

    genus_cols <- stats::setNames(
      lighten(rep_len(genus_base, length(unique_genus)), 0.15),
      unique_genus
    )
    species_cols <- stats::setNames(
      lighten(rep_len(species_base, length(unique_species)), 0.15),
      unique_species
    )

    list(
      df = anno_df,
      colors = list(Genus = genus_cols, Species = species_cols)
    )
  }

  # --- Helper: build row labels (organism strain) + col labels (accession) ---
  make_row_labels <- function(mat_keys) {
    meta <- dnmb$genome_meta %>%
      dplyr::filter(genome_key %in% mat_keys) %>%
      dplyr::arrange(match(genome_key, mat_keys))
    vapply(seq_along(mat_keys), function(i) {
      org <- meta$organism[i]
      st  <- meta$strain[i]
      if (is.na(org) || !nzchar(org)) return(mat_keys[i])
      tokens <- strsplit(trimws(org), "\\s+")[[1]]
      gs <- paste(tokens[1:min(2, length(tokens))], collapse = " ")
      rest <- if (!is.na(st) && nzchar(st)) st else ""
      if (nzchar(rest)) paste0(gs, " ", rest) else gs
    }, character(1))
  }

  make_col_labels <- function(mat_keys) {
    # Column labels = GenBank accession (genome_key) -- provides
    # non-redundant information vs the row labels (organism name).
    mat_keys
  }

  # --- Render one heatmap -----------------------------------------
  render_heatmap <- function(mat, title, val_range, palette_name, output_file) {
    if (is.null(mat)) return(NULL)

    anno_info <- build_annotation(rownames(mat))
    row_labels <- make_row_labels(rownames(mat))
    col_labels <- make_col_labels(colnames(mat))

    # Color function: each metric gets its own palette.
    # ANI uses the classic pyANI blue->yellow->red heat gradient;
    # POCP uses the colorspace "Vik" diverging scheme.
    if (palette_name == "ani_custom") {
      col_fun <- circlize::colorRamp2(
        seq(val_range[1], val_range[2], length.out = 5),
        c("#2C5F7A", "#88C0D0", "#FFDC91", "#E18727", "#CA0020")
      )
    } else {
      col_fun <- circlize::colorRamp2(
        seq(val_range[1], val_range[2], length.out = 100),
        colorspace::diverge_hcl(100, palette_name)
      )
    }

    # Clamp values below range floor to the floor so the gradient
    # starts clean (e.g. 55% POCP -> shows as 60% color, not white).
    mat_clamped <- pmax(mat, val_range[1])

    ha <- ComplexHeatmap::HeatmapAnnotation(
      Species = anno_info$df$Species,
      Genus   = anno_info$df$Genus,
      col     = anno_info$colors,
      annotation_name_side = "left",
      annotation_legend_param = list(
        Species = list(title = "Species", ncol = 1),
        Genus   = list(title = "Genus", ncol = 1)
      )
    )

    # Force square cells by setting equal width and height.
    n <- nrow(mat_clamped)
    cell_size <- grid::unit(5, "mm")

    ht <- ComplexHeatmap::Heatmap(
      mat_clamped,
      name = title,
      col  = col_fun,
      top_annotation = ha,
      clustering_distance_rows = "euclidean",
      clustering_distance_columns = "euclidean",
      clustering_method_rows = "complete",
      clustering_method_columns = "complete",
      show_row_dend = TRUE,
      show_column_dend = TRUE,
      row_labels    = row_labels,
      column_labels = col_labels,
      row_names_gp    = grid::gpar(fontsize = 8, fontface = "italic"),
      column_names_gp = grid::gpar(fontsize = 8),
      column_names_rot = 45,
      width  = cell_size * n,
      height = cell_size * n,
      border = TRUE,
      rect_gp = grid::gpar(col = "white", lwd = 0.5),
      heatmap_legend_param = list(
        title = paste0(title, " value"),
        at = seq(val_range[1], val_range[2], length.out = 5),
        labels = sprintf("%.0f", seq(val_range[1], val_range[2], length.out = 5)),
        legend_direction = "vertical"
      )
    )

    if (!is.null(output_file)) {
      dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
      n <- nrow(mat)
      # Square figure so heatmap cells are perfect squares.
      fig_size <- max(8, 4 + n * 0.40)
      grDevices::pdf(output_file, width = fig_size, height = fig_size)
      ComplexHeatmap::draw(
        ht,
        heatmap_legend_side = "right",
        annotation_legend_side = "right",
        merge_legend = TRUE
      )
      grDevices::dev.off()
      message("heatmap written to: ", output_file)
    }
    ht
  }

  # --- ANI heatmap (80-100, classic pyANI blue->red gradient) -------
  ani_mat <- load_matrix("ani_matrix.parquet", "ani_percent")
  result$ani <- render_heatmap(ani_mat, "ANI", c(80, 100), "ani_custom", output_file_ani)

  # --- POCP heatmap (60-100, "Vik" brown/blue diverging) ----------
  pocp_mat <- load_matrix("pocp_matrix.parquet", "pocp_percent")
  result$pocp <- render_heatmap(pocp_mat, "POCP", c(60, 100), "Vik", output_file_pocp)

  invisible(result)
}
