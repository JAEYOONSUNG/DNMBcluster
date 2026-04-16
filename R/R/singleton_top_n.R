#' Per-genome singleton (unique-gene) counts and top annotations
#'
#' For each input genome, counts how many `unique`-category clusters
#' are owned by that genome and, alongside the count, surfaces the
#' top-N annotation strings by representative aa_length (the longest
#' unique proteins are usually the most biologically interesting —
#' they exclude short hypotheticals and tiny peptides that dominate
#' unique rows in most bacterial assemblies).
#'
#' Bar length = unique cluster count. Text next to the bar =
#' top-N product lines (`•` separated). Use this plot to eyeball
#' which strains carry the most distinctive genetic cargo and
#' preview what that cargo looks like.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @param top_n How many annotations to show per genome (default 3).
#' @return A ggplot object.
#' @export
singleton_top_n <- function(dnmb, output_file = NULL, top_n = 3L) {
  unique_clusters <- dnmb$cluster_summary %>%
    dplyr::filter(category == "unique") %>%
    dplyr::select(cluster_id, representative_uid, representative_product)

  # Each unique cluster has exactly one genome_uid; pull it via the
  # centroid row in the clusters table.
  cluster_genome <- dnmb$clusters %>%
    dplyr::filter(is_centroid) %>%
    dplyr::select(cluster_id, genome_uid, protein_uid)

  aa_len <- dnmb$id_map %>%
    dplyr::select(protein_uid, aa_length)

  merged <- unique_clusters %>%
    dplyr::inner_join(cluster_genome, by = "cluster_id") %>%
    dplyr::left_join(aa_len, by = "protein_uid") %>%
    dplyr::left_join(
      dnmb$genome_meta %>% dplyr::select(genome_uid, genome_key),
      by = "genome_uid"
    )

  per_genome_counts <- merged %>%
    dplyr::group_by(genome_key) %>%
    dplyr::summarise(n_unique = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_unique))

  top_annotations <- merged %>%
    dplyr::filter(!is.na(representative_product), representative_product != "") %>%
    dplyr::arrange(dplyr::desc(aa_length)) %>%
    dplyr::group_by(genome_key) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::summarise(
      top_products = paste0("• ", representative_product, collapse = "\n"),
      .groups = "drop"
    )

  df <- per_genome_counts %>%
    dplyr::left_join(top_annotations, by = "genome_key") %>%
    dplyr::mutate(
      top_products = tidyr::replace_na(top_products, "(no annotations)"),
      genome_key   = factor(genome_key, levels = rev(genome_key))
    )

  max_n <- max(df$n_unique)

  plot <- ggplot2::ggplot(
    df, ggplot2::aes(x = genome_key, y = n_unique)
  ) +
    ggplot2::geom_col(fill = "#D06461", width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = format(n_unique, big.mark = ",")),
      hjust = -0.15, size = 3.3, color = "#303030"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = max_n * 1.05, label = top_products),
      hjust = 0, vjust = 0.5,
      size = 2.8, color = "#505050",
      lineheight = 0.95
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      expand = ggplot2::expansion(mult = c(0, 0.55))
    ) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::labs(
      title    = "Unique (strain-specific) gene contributions",
      subtitle = sprintf(
        "Top %d longest unique representatives shown per genome",
        top_n
      ),
      x = NULL,
      y = "Unique cluster count"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor    = ggplot2::element_blank(),
      panel.grid.major.y  = ggplot2::element_blank(),
      plot.margin         = ggplot2::margin(12, 12, 12, 12)
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, plot, width = 13, height = 7, dpi = 300)
    message("singleton_top_n written to: ", output_file)
  }

  invisible(plot)
}
