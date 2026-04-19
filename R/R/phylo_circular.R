#' Anvi'o-style circular pangenome display
#'
#' Circular tree with concentric data rings -- the classic anvi'o
#' pangenome display recreated in ggtree + ggnewscale. Separate
#' from the rectangular ``phylo_tree_plot`` which emphasizes the
#' heatmap; this one emphasizes the radial visual density.
#'
#' Rings (inside-out):
#' 1. Core-gene phylogeny (circular, ladderized)
#' 2. Category dots (core navy / accessory orange / unique red)
#' 3. GC% gradient
#' 4. Genome size
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level results directory.
#' @param output_file Optional PDF path.
#' @return A ggtree plot object, or NULL if tree is absent.
#' @export
phylo_circular <- function(dnmb, results_dir, output_file = NULL) {
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    warning("ggtree not installed -- skipping phylo_circular")
    return(invisible(NULL))
  }

  tree_path <- file.path(results_dir, "dnmb", "processed", "phylo_tree.nwk")
  if (!file.exists(tree_path)) {
    message("phylo_circular: tree not found -- skipping")
    return(invisible(NULL))
  }

  tree <- treeio::read.newick(tree_path, node.label = "support")

  # Tip label decoration
  tips_df <- tibble::tibble(label = tree@phylo$tip.label) %>%
    dplyr::left_join(
      dnmb$genome_meta %>%
        dplyr::select(genome_key, organism, strain) %>%
        dplyr::mutate(
          pretty = dplyr::case_when(
            !is.na(organism) & !is.na(strain) & strain != "" &
              !mapply(grepl, strain, organism, fixed = TRUE) ~
              paste0(organism, " ", strain),
            !is.na(organism) ~ organism,
            TRUE ~ genome_key
          )
        ) %>%
        dplyr::select(label = genome_key, pretty),
      by = "label"
    ) %>%
    dplyr::mutate(pretty = dplyr::coalesce(pretty, label))

  tree@phylo$tip.label <- tips_df$pretty[match(
    tree@phylo$tip.label, tips_df$label
  )]

  # Per-genome metadata
  cat_lookup <- dnmb$cluster_summary %>% dplyr::select(cluster_id, category)
  cat_wide <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(cat_lookup, by = "cluster_id") %>%
    dplyr::count(genome_uid, category, name = "n") %>%
    tidyr::pivot_wider(names_from = category, values_from = n, values_fill = 0)

  meta <- dnmb$genome_meta %>%
    dplyr::mutate(label = tips_df$pretty[match(genome_key, tips_df$label)]) %>%
    dplyr::left_join(cat_wide, by = "genome_uid") %>%
    dplyr::select(label, gc_percent, total_length, n_cds,
                  dplyr::any_of(c("core","accessory","unique"))) %>%
    dplyr::mutate(
      core = dplyr::coalesce(core, 0L),
      accessory = dplyr::coalesce(accessory, 0L),
      unique = dplyr::coalesce(unique, 0L)
    ) %>%
    as.data.frame()

  max_x <- max(ggtree::fortify(tree@phylo)$x)
  r1 <- max_x * 1.3
  r2 <- max_x * 1.45
  r3 <- max_x * 1.6
  r4 <- max_x * 1.85
  r5 <- max_x * 2.1

  p <- ggtree::ggtree(tree@phylo, layout = "circular", size = 0.6, ladderize = TRUE) +
    ggplot2::xlim(-0.2, max_x * 3.5)

  p <- ggtree::`%<+%`(p, meta)

  p <- p +
    ggtree::geom_tiplab2(
      ggplot2::aes(label = label),
      fontface = "italic", size = 2.8, color = "#1F2E4A",
      align = TRUE, linesize = 0.12, linetype = 3,
      offset = max_x * 1.2, show.legend = FALSE
    ) +
    # Category dots
    ggplot2::geom_point(ggplot2::aes(x = r1, size = core),
                        color = "#2C5F7A", alpha = 0.7, shape = 16) +
    ggplot2::geom_point(ggplot2::aes(x = r2, size = accessory),
                        color = "#F2A766", alpha = 0.7, shape = 16) +
    ggplot2::geom_point(ggplot2::aes(x = r3, size = unique),
                        color = "#D06461", alpha = 0.7, shape = 16) +
    ggplot2::scale_size_continuous(name = "Cluster\ncount", range = c(1, 5))

  p <- p +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(ggplot2::aes(x = r4, color = gc_percent),
                        size = 3, shape = 15) +
    ggplot2::scale_color_distiller(name = "GC%", palette = "YlGnBu", direction = 1,
                                   na.value = "grey80")

  p <- p +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(ggplot2::aes(x = r5, color = total_length / 1e6),
                        size = 3, shape = 15) +
    ggplot2::scale_color_distiller(name = "Size\n(Mb)", palette = "YlOrRd", direction = 1,
                                   na.value = "grey80")

  p <- p +
    ggplot2::labs(
      title = "Pangenome circular display",
      subtitle = "Core (navy) / Accessory (orange) / Unique (red) | GC% | Genome Mb"
    ) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "#F4F1EA", color = NA),
      panel.background = ggplot2::element_rect(fill = "#F4F1EA", color = NA),
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.35, "cm"),
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, color = "#505050")
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, p, width = 14, height = 12, dpi = 300)
    message("phylo_circular written to: ", output_file)
  }
  invisible(p)
}
