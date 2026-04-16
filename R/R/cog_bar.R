#' COG functional category distribution from EggNOG annotations
#'
#' Reads ``dnmb/processed/eggnog_annotations.parquet`` (from the
#' ``--annotate`` stage) and draws a stacked bar of single-letter COG
#' categories faceted by pan-genome category (core / accessory /
#' unique). Uses the standard NCBI COG color scheme.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level results directory.
#' @param output_file Optional PDF path.
#' @return A ggplot object, or NULL if eggnog_annotations.parquet is absent.
#' @export
cog_bar <- function(dnmb, results_dir, output_file = NULL) {
  egg_path <- file.path(
    results_dir, "dnmb", "processed", "eggnog_annotations.parquet"
  )
  if (!file.exists(egg_path)) {
    message("cog_bar: eggnog_annotations.parquet not found — skipping")
    return(invisible(NULL))
  }

  egg <- tibble::as_tibble(arrow::read_parquet(egg_path))
  cat_lookup <- dnmb$cluster_summary %>%
    dplyr::select(cluster_id, category)

  df <- egg %>%
    dplyr::filter(cog_category != "") %>%
    dplyr::left_join(cat_lookup, by = "cluster_id") %>%
    dplyr::mutate(
      # Some entries have multi-letter COG (e.g. "EG"); take the first
      cog_letter = substr(cog_category, 1, 1),
      category = factor(category, levels = c("core", "accessory", "unique"))
    ) %>%
    dplyr::count(category, cog_letter, name = "n") %>%
    dplyr::group_by(category) %>%
    dplyr::mutate(pct = n / sum(n) * 100) %>%
    dplyr::ungroup()

  # Standard COG single-letter descriptions
  cog_desc <- c(
    J = "Translation", A = "RNA processing", K = "Transcription",
    L = "Replication", B = "Chromatin", D = "Cell cycle",
    Y = "Nuclear", V = "Defense", T = "Signal transduction",
    M = "Cell wall", N = "Motility", Z = "Cytoskeleton",
    W = "Extracellular", U = "Secretion", O = "Chaperones",
    X = "Mobilome", C = "Energy", G = "Carbohydrate",
    E = "Amino acid", F = "Nucleotide", H = "Coenzyme",
    I = "Lipid", P = "Inorganic ion", Q = "Secondary metabolite",
    R = "General prediction", S = "Unknown function"
  )

  # Reorder by total count descending
  letter_order <- df %>%
    dplyr::group_by(cog_letter) %>%
    dplyr::summarise(total = sum(n), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total)) %>%
    dplyr::pull(cog_letter)
  df$cog_letter <- factor(df$cog_letter, levels = rev(letter_order))

  # Labels for x axis: "S - Unknown function"
  x_labels <- paste0(letter_order, " - ", cog_desc[letter_order])
  x_labels[is.na(x_labels)] <- letter_order[is.na(x_labels)]

  # Color palette: use journal-soft 20 colors recycled
  n_cats <- length(letter_order)
  lighten <- function(cols, amount = 0.30) {
    m <- grDevices::col2rgb(cols) / 255
    m <- m + (1 - m) * amount
    grDevices::rgb(m[1, ], m[2, ], m[3, ])
  }
  base <- c(
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#E18727", "#7876B1", "#EFC000",
    "#20854E", "#925E9F", "#CD534C", "#6F99AD", "#EE4C97",
    "#0073C2", "#B09C85", "#7AA6DC", "#42B540", "#FFDC91"
  )
  pal <- lighten(rep_len(base, n_cats), 0.30)
  names(pal) <- rev(letter_order)

  p <- ggplot2::ggplot(
    df, ggplot2::aes(x = cog_letter, y = pct, fill = cog_letter)
  ) +
    ggplot2::geom_col(width = 0.8, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::scale_x_discrete(labels = setNames(
      paste0(letter_order, " - ", cog_desc[letter_order]),
      rev(letter_order)
    )) +
    ggplot2::scale_y_continuous(labels = function(v) paste0(v, "%")) +
    ggplot2::facet_wrap(~ category, ncol = 1, scales = "free_y") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "COG functional category distribution (EggNOG)",
      subtitle = "Per-CDS annotation, faceted by pan-genome category",
      x = NULL, y = "% of annotated CDS"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "#F2F2F2"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, p, width = 10, height = 10, dpi = 300)
    message("cog_bar written to: ", output_file)
  }

  invisible(p)
}
