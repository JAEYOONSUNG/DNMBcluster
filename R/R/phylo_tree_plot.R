#' Core-gene phylogeny — rectangular tree + side annotation columns
#'
#' Rectangular layout with right-aligned heatmap columns showing
#' per-genome pan-genome composition, GC%, genome size, and CDS
#' count via ggtree's ``gheatmap()``. Gain/loss labels annotate
#' tip branches. Clean, paper-friendly layout.
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

  # --- Per-genome metadata for heatmap columns ------------------
  cat_lookup <- dnmb$cluster_summary %>% dplyr::select(cluster_id, category)
  cat_wide <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(cat_lookup, by = "cluster_id") %>%
    dplyr::count(genome_uid, category, name = "n") %>%
    tidyr::pivot_wider(names_from = category, values_from = n, values_fill = 0)

  meta <- dnmb$genome_meta %>%
    dplyr::mutate(
      label = tips_df$pretty[match(genome_key, tips_df$label)]
    ) %>%
    dplyr::left_join(
      cat_wide, by = "genome_uid"
    )

  # Build the heatmap data frame: rows = tip labels, cols = metrics
  hm_df <- data.frame(
    row.names    = meta$label,
    Core         = meta$core,
    Accessory    = meta$accessory,
    Unique       = meta$unique,
    `GC (%)`     = round(meta$gc_percent, 1),
    `Size (Mb)`  = round(meta$total_length / 1e6, 2),
    CDS          = meta$n_cds,
    check.names  = FALSE,
    stringsAsFactors = FALSE
  )

  # --- Gain / loss -----------------------------------------------
  gl_path <- file.path(results_dir, "dnmb", "processed", "gain_loss.parquet")
  gl_df <- NULL
  if (file.exists(gl_path)) {
    gl_raw <- tibble::as_tibble(arrow::read_parquet(gl_path))
    tree_fort <- ggtree::fortify(tree@phylo)
    n_tips <- length(tree@phylo$tip.label)

    # Map child_node names to ggtree node numbers.
    # Tips: original_tip_order[i] → node i.
    # Internal: N0, N1, ... → node (n_tips + 1 + id). But ggtree's
    # internal numbering doesn't follow our N# scheme — it uses ape's
    # node numbering. The fortify data has 'label' for internal nodes
    # = whatever was in tree@phylo$node.label BEFORE we renamed them.
    # Since our gain_loss.py assigns N0..Nk, and tree@phylo$node.label
    # still holds the IQ-TREE support strings, we can't match directly.
    # Instead: match tips by name, and for internal branches match by
    # topology — map (parent_node, child_node) where child is a tip.
    # For full coverage, build a node name → ggtree node# map from
    # the fortify data for tips, then use tree topology for internals.
    tip_map <- data.frame(
      child_node = original_tip_order,
      node = seq_len(n_tips),
      stringsAsFactors = FALSE
    )

    # For internal nodes: gain_loss.py uses N0, N1... from a level-order
    # traversal. ggtree's internal nodes are numbered (n_tips+1)..
    # in ape's postorder. We match by the CHILD field: if child_node
    # is a genome_key (tip), we match directly. If child_node is N#
    # (internal), we need to identify it. For now, annotate only
    # branches whose child_node is a recognized tip — this covers
    # the terminal branches where most visible events occur.
    # Internal branches annotated as sum labels on node points later.
    gl_df <- gl_raw %>%
      dplyr::left_join(tip_map, by = "child_node") %>%
      dplyr::filter(n_gained > 0 | n_lost > 0) %>%
      dplyr::mutate(
        gl_label = dplyr::case_when(
          n_gained > 0 & n_lost > 0 ~ sprintf("+%d/-%d", n_gained, n_lost),
          n_gained > 0              ~ sprintf("+%d", n_gained),
          TRUE                      ~ sprintf("-%d", n_lost)
        ),
        is_tip = !is.na(node)
      )

    # For tip branches: get x, y from fortify
    gl_tip_df <- gl_df %>%
      dplyr::filter(is_tip) %>%
      dplyr::left_join(tree_fort[, c("node", "x", "y")], by = "node")

    # For internal branches: aggregate total gain/loss and show as
    # a summary annotation on the root or as node labels.
    gl_internal <- gl_df %>% dplyr::filter(!is_tip)
  }

  # --- Build rectangular tree -----------------------------------
  p <- ggtree::ggtree(
    tree@phylo, size = 0.7, ladderize = TRUE
  ) +
    ggtree::geom_tiplab(
      fontface = "italic", size = 3.2, color = "#1F2E4A",
      offset = 0.002
    ) +
    ggtree::geom_nodepoint(
      ggplot2::aes(subset = !isTip),
      size = 2.0, shape = 21, fill = "#F8F8F8", color = "#2C5F7A"
    ) +
    ggtree::geom_treescale(x = 0, y = -0.5, fontsize = 2.5, linesize = 0.5)

  # --- Heatmap columns via gheatmap -----------------------------
  p2 <- ggtree::gheatmap(
    p, hm_df,
    offset       = 0.05,
    width        = 0.6,
    colnames_angle   = 45,
    colnames_offset_y = 0.2,
    font.size    = 2.8,
    hjust        = 0,
    color        = "white",
    low          = "#F5F5F5",
    high         = "#2C5F7A"
  ) +
    ggplot2::scale_fill_viridis_c(
      name = "Value", option = "C", na.value = "grey90"
    )

  # --- Gain/loss paired dots (after gheatmap) ---------------------
  # Two dots per tip branch — green (gained) + red (lost), sizes
  # scaled by sqrt(count). Positioned at the branch midpoint with a
  # tiny y offset so they sit above/below the branch line and don't
  # overlap each other.
  if (!is.null(gl_df) && exists("gl_tip_df") && nrow(gl_tip_df) > 0L) {
    tree_fort_fresh <- ggtree::fortify(tree@phylo)
    parent_x <- tree_fort_fresh$x[match(
      tree_fort_fresh$parent[gl_tip_df$node], tree_fort_fresh$node
    )]
    gl_tip_df$mid_x <- (gl_tip_df$x + parent_x) / 2
    gl_tip_df$total <- gl_tip_df$n_gained + gl_tip_df$n_lost

    dot_df <- gl_tip_df %>%
      dplyr::filter(total > 0) %>%
      as.data.frame()

    if (nrow(dot_df) > 0L) {
      max_evt <- max(c(dot_df$n_gained, dot_df$n_lost), na.rm = TRUE)

      p2 <- p2 +
        # Gained (green) — nudged slightly above
        ggnewscale::new_scale("size") +
        ggplot2::geom_point(
          data = dot_df,
          ggplot2::aes(x = mid_x, y = y + 0.18,
                       size = n_gained),
          color = "#4CAF50", alpha = 0.80, shape = 16,
          inherit.aes = FALSE
        ) +
        # Lost (red) — nudged slightly below
        ggplot2::geom_point(
          data = dot_df,
          ggplot2::aes(x = mid_x, y = y - 0.18,
                       size = n_lost),
          color = "#E53935", alpha = 0.80, shape = 16,
          inherit.aes = FALSE
        ) +
        ggplot2::scale_size_continuous(
          name = "Gene\nevents",
          range = c(1, 6),
          limits = c(0, max_evt),
          guide = ggplot2::guide_legend(
            override.aes = list(color = "#555555")
          )
        ) +
        # Text label: +N/-N
        ggplot2::geom_text(
          data = dot_df,
          ggplot2::aes(x = mid_x, y = y, label = gl_label),
          inherit.aes = FALSE,
          size = 1.6, color = "#303030", fontface = "bold",
          hjust = 0.5, vjust = 0.5
        )
    }
  }

  p2 <- p2 +
    ggplot2::labs(
      title = "Core-gene phylogeny  (IQ-TREE)",
      subtitle = paste0(
        "Heatmap: Core / Accessory / Unique counts, GC%, Genome size, CDS count",
        if (!is.null(gl_df)) "  |  Brown = +gained / -lost (Fitch parsimony)" else ""
      )
    ) +
    ggtree::theme_tree2(base_size = 12) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = "#FCFCFC", color = NA),
      panel.background = ggplot2::element_rect(fill = "#FCFCFC", color = NA),
      legend.position  = "right",
      legend.key.size  = ggplot2::unit(0.4, "cm"),
      legend.text      = ggplot2::element_text(size = 7),
      legend.title     = ggplot2::element_text(size = 8, face = "bold"),
      plot.title       = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 8, color = "#505050")
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, p2, width = 14, height = 8, dpi = 300)
    message("phylo_tree_plot written to: ", output_file)
  }

  invisible(p2)
}
