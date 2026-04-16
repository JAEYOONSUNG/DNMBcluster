#' Core-gene IQ-TREE best tree visualization via ggtree
#'
#' Reads `dnmb/processed/phylo_tree.nwk` (written by the Python
#' phylogenomics stage — core concat + MAFFT + IQ-TREE fast mode)
#' and draws it as a publication-quality tree with:
#'
#' - Tips labeled with genome key / organism + strain
#' - Internal nodes colored by SH-aLRT + UFBoot support
#' - Tidy rectangular layout with a scale bar
#' - Softer-than-anvio palette (navy accents, off-white background)
#'
#' Falls back to a short-circuit warning when the tree file doesn't
#' exist — phylogenomics is opt-in via the `--phylo` CLI flag, so
#' this plot only fires when the file is actually there.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level pipeline results directory (needed
#'   to locate `dnmb/processed/phylo_tree.nwk`).
#' @param output_file Optional PDF path.
#' @return A ggtree plot object, or NULL if the tree file is absent.
#' @export
phylo_tree_plot <- function(dnmb, results_dir, output_file = NULL) {
  if (!requireNamespace("ggtree", quietly = TRUE) ||
      !requireNamespace("treeio", quietly = TRUE) ||
      !requireNamespace("ape", quietly = TRUE)) {
    warning("ggtree/treeio/ape not installed — skipping phylo_tree_plot")
    return(invisible(NULL))
  }

  tree_path <- file.path(results_dir, "dnmb", "processed", "phylo_tree.nwk")
  if (!file.exists(tree_path)) {
    message("phylo_tree_plot: ", tree_path, " not found — skipping")
    return(invisible(NULL))
  }

  tree <- treeio::read.newick(tree_path, node.label = "support")
  # IQ-TREE writes node labels as "<SH_aLRT>/<UFBoot>" strings; parse
  # the UFBoot half (the one most downstream users treat as "the"
  # support) and store it as a numeric tip/node attribute.
  node_labels <- tree@phylo$node.label
  if (!is.null(node_labels)) {
    ufboot <- vapply(
      strsplit(node_labels, "/", fixed = TRUE),
      function(parts) {
        if (length(parts) >= 2L) suppressWarnings(as.numeric(parts[[2]]))
        else if (length(parts) == 1L) suppressWarnings(as.numeric(parts[[1]]))
        else NA_real_
      },
      numeric(1)
    )
  } else {
    ufboot <- rep(NA_real_, tree@phylo$Nnode)
  }

  # Tip label decoration — prefer organism+strain from genome_meta
  # when available, fall back to genome_key. RefSeq ``organism``
  # strings already include the strain suffix for many entries
  # (e.g. ``"Geobacillus kaustophilus HTA426"`` with strain
  # ``"HTA426"``), so naively concatenating produces a duplicate tail.
  # Skip the strain append when the strain string is already a
  # substring of organism.
  tips_df <- tibble::tibble(
    label = tree@phylo$tip.label
  ) %>%
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

  tree@phylo$tip.label <- tips_df$pretty[match(tree@phylo$tip.label, tips_df$label)]

  node_support_df <- tibble::tibble(
    node = seq_len(tree@phylo$Nnode) + length(tree@phylo$tip.label),
    ufboot = ufboot
  )

  p <- ggtree::ggtree(
    tree@phylo, size = 0.7, ladderize = TRUE
  ) +
    ggtree::geom_tiplab(
      size = 3.4, color = "#1F2E4A", align = TRUE, linesize = 0.2,
      offset = 0.002
    ) +
    ggtree::geom_nodepoint(
      ggplot2::aes(subset = !isTip),
      size = 2.2, shape = 21, fill = "#F8F8F8", color = "#2C5F7A"
    ) +
    ggplot2::labs(
      title = "Core-gene phylogeny  (IQ-TREE best tree, LG+G4 --fast)",
      subtitle = "Tips = input genomes    Node dots = internal branches with SH-aLRT/UFBoot support"
    ) +
    ggtree::theme_tree2(base_size = 13) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "#FCFCFC", color = NA),
      panel.background = ggplot2::element_rect(fill = "#FCFCFC", color = NA)
    ) +
    ggplot2::xlim(NA, max(tree@phylo$edge.length * 5))

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, p, width = 10, height = 7, dpi = 300)
    message("phylo_tree_plot written to: ", output_file)
  }

  invisible(p)
}
