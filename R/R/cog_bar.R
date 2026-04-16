#' COG enrichment dot plot — core vs accessory odds ratio
#'
#' For each COG functional category, computes the Fisher's exact test
#' odds ratio of being in the core genome vs the accessory+unique pool.
#' Plots as a horizontal dot plot: x = log2(OR), dot size =
#' -log10(p-value), color = enriched in core (blue) or accessory (orange).
#'
#' Directly answers "what functions define the core genome?" — the
#' standard question reviewers expect in a pan-genome paper.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level results directory (for eggnog_annotations).
#' @param output_file Optional PDF path.
#' @return A ggplot object, or NULL if eggnog data is absent.
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
  cat_lookup <- dnmb$cluster_summary %>% dplyr::select(cluster_id, category)

  df <- egg %>%
    dplyr::filter(cog_category != "") %>%
    dplyr::left_join(cat_lookup, by = "cluster_id") %>%
    dplyr::mutate(
      cog_letter = substr(cog_category, 1, 1),
      is_core = category == "core"
    )

  # COG descriptions
  cog_desc <- c(
    J = "Translation", A = "RNA processing", K = "Transcription",
    L = "Replication/repair", B = "Chromatin", D = "Cell cycle",
    V = "Defense", T = "Signal transduction",
    M = "Cell wall/membrane", N = "Motility",
    U = "Secretion", O = "Chaperones/PTM",
    X = "Mobilome", C = "Energy production",
    G = "Carbohydrate metabolism", E = "Amino acid metabolism",
    F = "Nucleotide metabolism", H = "Coenzyme metabolism",
    I = "Lipid metabolism", P = "Inorganic ion transport",
    Q = "Secondary metabolism",
    R = "General prediction", S = "Unknown function"
  )

  # Fisher's exact test per COG letter: core vs non-core
  letters <- sort(unique(df$cog_letter))
  total_core <- sum(df$is_core)
  total_noncore <- sum(!df$is_core)

  results <- lapply(letters, function(l) {
    in_l  <- df$cog_letter == l
    a <- sum(in_l & df$is_core)      # core + this COG
    b <- sum(in_l & !df$is_core)     # non-core + this COG
    c <- total_core - a              # core + other COG
    d <- total_noncore - b           # non-core + other COG
    ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2))
    data.frame(
      cog_letter = l,
      description = if (l %in% names(cog_desc)) paste0(l, " - ", cog_desc[l]) else l,
      log2OR = log2(max(ft$estimate, 1e-6)),
      pval = ft$p.value,
      n = a + b,
      core_n = a,
      noncore_n = b,
      stringsAsFactors = FALSE
    )
  })
  or_df <- do.call(rbind, results)
  or_df$neg_log_p <- -log10(pmax(or_df$pval, 1e-300))
  or_df$enriched <- ifelse(or_df$log2OR > 0, "Core-enriched", "Accessory-enriched")
  or_df$sig <- ifelse(or_df$pval < 0.001, "***",
                       ifelse(or_df$pval < 0.01, "**",
                              ifelse(or_df$pval < 0.05, "*", "")))

  p <- ggplot2::ggplot(
    or_df,
    ggplot2::aes(
      x = stats::reorder(description, log2OR),
      y = log2OR,
      size = neg_log_p,
      color = enriched
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::geom_text(
      ggplot2::aes(label = sig),
      color = "#303030", size = 3, vjust = -0.8,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      name = "Enrichment",
      values = c("Core-enriched" = "#2C5F7A", "Accessory-enriched" = "#F2A766")
    ) +
    ggplot2::scale_size_continuous(
      name = expression(-log[10](p)),
      range = c(2, 8)
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "COG functional enrichment: Core vs Accessory genome",
      subtitle = "Fisher's exact test odds ratio per COG category  |  * p<0.05  ** p<0.01  *** p<0.001",
      x = NULL,
      y = expression(log[2](Odds~Ratio))
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, p, width = 10, height = 8, dpi = 300)
    message("cog_bar written to: ", output_file)
  }

  invisible(p)
}
