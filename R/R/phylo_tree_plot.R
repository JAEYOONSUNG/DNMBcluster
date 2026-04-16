#' Core-gene phylogeny — circular multi-ring visualization
#'
#' Publication-quality phylogenomic figure built with ggtree +
#' ggnewscale (no ggtreeExtra dependency — avoids version conflicts).
#' Uses the same manual-ring pattern as the aprE/araA gene-mining
#' scripts: ``geom_tippoint`` at offset x positions with
#' ``new_scale("color")`` between layers.
#'
#' Rings (inside-out):
#' 1. Tree + italic tip labels
#' 2. Category breakdown (3 dots: core / accessory / unique count)
#' 3. GC% gradient
#' 4. Genome size (Mb)
#' 5. Gain/loss branch labels (if available)
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level pipeline results directory.
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

  # --- Tip label decoration -------------------------------------
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

  original_tip_order <- tips_df$label
  tree@phylo$tip.label <- tips_df$pretty[match(
    tree@phylo$tip.label, tips_df$label
  )]

  # --- Per-genome metadata for ring layers ----------------------
  meta <- dnmb$genome_meta %>%
    dplyr::mutate(
      label = tips_df$pretty[match(genome_key, tips_df$label)]
    )

  # Category counts per genome
  cat_lookup <- dnmb$cluster_summary %>% dplyr::select(cluster_id, category)
  cat_wide <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(cat_lookup, by = "cluster_id") %>%
    dplyr::count(genome_uid, category, name = "n") %>%
    tidyr::pivot_wider(names_from = category, values_from = n, values_fill = 0) %>%
    dplyr::left_join(
      meta %>% dplyr::select(genome_uid, label),
      by = "genome_uid"
    )

  # Merge into a single metadata frame keyed by 'label' (= pretty
  # tip name) with no duplicate or conflicting column names.
  ring_data <- meta %>%
    dplyr::select(genome_uid, label, gc_percent, total_length, n_cds) %>%
    dplyr::left_join(
      cat_wide %>% dplyr::select(genome_uid, dplyr::any_of(c("core","accessory","unique"))),
      by = "genome_uid"
    ) %>%
    dplyr::select(-genome_uid) %>%
    dplyr::mutate(
      core      = dplyr::coalesce(core, 0L),
      accessory = dplyr::coalesce(accessory, 0L),
      unique    = dplyr::coalesce(unique, 0L)
    ) %>%
    as.data.frame()

  # --- Gain / loss -----------------------------------------------
  gl_path <- file.path(results_dir, "dnmb", "processed", "gain_loss.parquet")
  gl_tips <- NULL
  if (file.exists(gl_path)) {
    gl_raw <- tibble::as_tibble(arrow::read_parquet(gl_path))
    tree_fort <- ggtree::fortify(tree@phylo)
    tip_node_map <- data.frame(
      child_node = original_tip_order,
      node = seq_len(length(original_tip_order)),
      stringsAsFactors = FALSE
    )
    gl_tips <- gl_raw %>%
      dplyr::inner_join(tip_node_map, by = "child_node") %>%
      dplyr::filter(n_gained > 0 | n_lost > 0) %>%
      dplyr::mutate(
        gl_label = dplyr::case_when(
          n_gained > 0 & n_lost > 0 ~ sprintf("+%d/-%d", n_gained, n_lost),
          n_gained > 0              ~ sprintf("+%d", n_gained),
          TRUE                      ~ sprintf("-%d", n_lost)
        )
      ) %>%
      dplyr::left_join(tree_fort[, c("node", "x", "y")], by = "node")
  }

  # --- Max branch length for ring x-offsets ----------------------
  max_x <- max(ggtree::fortify(tree@phylo)$x) * 1.0
  r1 <- max_x * 1.3   # core ring
  r2 <- max_x * 1.45  # accessory ring
  r3 <- max_x * 1.6   # unique ring
  r4 <- max_x * 1.85  # GC ring
  r5 <- max_x * 2.1   # genome size ring

  # --- Build tree ------------------------------------------------
  p <- ggtree::ggtree(
    tree@phylo, layout = "circular", size = 0.6, ladderize = TRUE
  ) +
    ggplot2::xlim(-0.2, max_x * 3.5)

  # Attach metadata
  p <- ggtree::`%<+%`(p, ring_data)

  # Tip labels
  p <- p +
    ggtree::geom_tiplab2(
      ggplot2::aes(label = label),
      fontface = "italic", size = 2.8, color = "#1F2E4A",
      align = TRUE, linesize = 0.12, linetype = 3,
      offset = max_x * 1.2, show.legend = FALSE
    )

  # Gain/loss branch text
  if (!is.null(gl_tips) && nrow(gl_tips) > 0L) {
    p <- p +
      ggplot2::geom_text(
        data = gl_tips,
        ggplot2::aes(x = x, y = y, label = gl_label),
        inherit.aes = FALSE,
        size = 1.8, color = "#8B4513", fontface = "bold",
        hjust = 1.1, vjust = -0.3
      )
  }

  # --- Ring 1-3: Category dots (core / accessory / unique) -------
  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(x = r1, size = core),
      color = "#2C5F7A", alpha = 0.7, shape = 16
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = r2, size = accessory),
      color = "#F2A766", alpha = 0.7, shape = 16
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = r3, size = unique),
      color = "#D06461", alpha = 0.7, shape = 16
    ) +
    ggplot2::scale_size_continuous(
      name = "Cluster\ncount", range = c(1, 5), guide = "legend"
    )

  # --- Ring 4: GC% gradient dot ---------------------------------
  p <- p +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(
      ggplot2::aes(x = r4, color = gc_percent),
      size = 3, shape = 15
    ) +
    ggplot2::scale_color_distiller(
      name = "GC%", palette = "YlGnBu", direction = 1,
      na.value = "grey80"
    )

  # --- Ring 5: Genome size (Mb) scaled dot ----------------------
  p <- p +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(
      ggplot2::aes(x = r5, color = total_length / 1e6),
      size = 3, shape = 15
    ) +
    ggplot2::scale_color_distiller(
      name = "Genome\nsize (Mb)", palette = "YlOrRd", direction = 1,
      na.value = "grey80"
    )

  # --- Theme + labels -------------------------------------------
  has_gl <- !is.null(gl_tips) && nrow(gl_tips) > 0L
  p <- p +
    ggplot2::labs(
      title = "Core-gene phylogeny",
      subtitle = paste0(
        "IQ-TREE best tree  |  Dots: core (navy) / accessory (orange) / unique (red)  |  ",
        "GC%  |  Genome Mb",
        if (has_gl) "  |  Brown = +gained / -lost" else ""
      )
    ) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = "#FAFAFA", color = NA),
      panel.background  = ggplot2::element_rect(fill = "#FAFAFA", color = NA),
      legend.position   = "right",
      legend.background = ggplot2::element_rect(fill = NA),
      legend.key.size   = ggplot2::unit(0.35, "cm"),
      legend.text       = ggplot2::element_text(size = 7),
      legend.title      = ggplot2::element_text(size = 8, face = "bold"),
      plot.title        = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle     = ggplot2::element_text(size = 8, color = "#505050")
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, p, width = 14, height = 12, dpi = 300)
    message("phylo_tree_plot written to: ", output_file)
  }

  invisible(p)
}
