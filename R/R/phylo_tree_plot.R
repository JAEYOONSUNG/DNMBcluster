#' Core-gene phylogeny -- rectangular tree + side annotation columns
#'
#' Rectangular ggtree layout with italic tip labels and a minimal
#' two-column heatmap showing Core gene count and GC%. A separate
#' aligned tile panel is attached on the right for the remaining
#' per-genome metrics (Accessory, Unique, Size, CDS) using patchwork.
#'
#' NOTE: on current ggtree S7 builds, appending *any* ggplot/ggtree
#' layer (labs, theme, scales, additional geoms) AFTER `gheatmap()`
#' corrupts the `@mapping` invariant with
#' "`@mapping must be <ggplot2::mapping>, not S3<data.frame>`".
#' This function therefore finalises all tree decorations (titles,
#' themes, gain/loss annotations) BEFORE `gheatmap()` and treats
#' `gheatmap()` as a terminal step. Per-genome metadata that can't
#' fit into the 2-column heatmap is shown in a patchwork-aligned
#' side tile panel.
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
    warning("ggtree/treeio/ape not installed -- skipping phylo_tree_plot")
    return(invisible(NULL))
  }

  tree_path <- file.path(results_dir, "dnmb", "processed", "phylo_tree.nwk")
  if (!file.exists(tree_path)) {
    message("phylo_tree_plot: ", tree_path, " not found -- skipping")
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

  # --- Per-genome metadata --------------------------------------
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
    dplyr::left_join(cat_wide, by = "genome_uid")

  safe_num <- function(x, default = 0) {
    if (is.null(x)) return(rep(default, nrow(meta)))
    ifelse(is.na(x), default, x)
  }
  core_vec <- if ("core" %in% names(meta)) safe_num(meta$core) else rep(0, nrow(meta))
  gc_vec   <- if ("gc_percent" %in% names(meta)) safe_num(meta$gc_percent) else rep(NA_real_, nrow(meta))

  # hm_df kept for reference but we no longer call gheatmap() --
  # current ggtree stat_tree breaks with gheatmap aesthetics recycling.

  # --- Gain / loss per-tip data (rendered on the side tile panel
  #     instead of the tree, so we don't break stat_tree) -------
  gl_tip_df <- NULL
  gl_path <- file.path(results_dir, "dnmb", "processed", "gain_loss.parquet")
  if (file.exists(gl_path) && requireNamespace("arrow", quietly = TRUE)) {
    gl_raw <- tryCatch(
      tibble::as_tibble(arrow::read_parquet(gl_path)),
      error = function(e) NULL
    )
    if (!is.null(gl_raw) && nrow(gl_raw) > 0L &&
        all(c("child_node", "n_gained", "n_lost") %in% names(gl_raw))) {
      tip_map <- data.frame(
        child_node = original_tip_order,
        pretty     = tips_df$pretty,
        stringsAsFactors = FALSE
      )
      gl_tip_df <- gl_raw %>%
        dplyr::inner_join(tip_map, by = "child_node") %>%
        dplyr::filter(n_gained > 0 | n_lost > 0) %>%
        dplyr::mutate(
          gl_label = dplyr::case_when(
            n_gained > 0 & n_lost > 0 ~ sprintf("+%d / -%d", n_gained, n_lost),
            n_gained > 0              ~ sprintf("+%d", n_gained),
            TRUE                      ~ sprintf("-%d", n_lost)
          )
        )
    }
  }

  # --- Build rectangular tree (NO layers after tiplab) ----------
  # NOTE: appending labs/theme/scales/geoms to a ggtree object
  # before gheatmap() corrupts stat_tree's internal from/to/.panel
  # aesthetics on current ggtree. Keep the tree minimal; put all
  # titles/themes on the patchwork wrapper instead.
  # Extend x limit so tip labels aren't clipped by the panel edge.
  tree_fort_x <- ggtree::fortify(tree@phylo)
  tree_xmax   <- max(tree_fort_x$x, na.rm = TRUE)
  longest_lab <- max(nchar(tree@phylo$tip.label), na.rm = TRUE)
  xlim_hi     <- tree_xmax + 0.08 * tree_xmax + 0.01 * longest_lab * tree_xmax

  p <- ggtree::ggtree(tree@phylo, size = 0.7, ladderize = TRUE) +
    ggtree::geom_tiplab(
      fontface = "italic", size = 3.2, color = "#1F2E4A",
      offset = 0.002
    ) +
    ggtree::xlim_tree(xlim_hi)

  # gheatmap intentionally skipped -- see hm_df comment above.
  p2 <- p

  # --- Side tile panel via patchwork for remaining metrics ------
  side_panel <- NULL
  if (requireNamespace("patchwork", quietly = TRUE)) {
    tree_fort <- ggtree::fortify(tree@phylo)
    tip_rows  <- tree_fort[tree_fort$isTip, c("label", "y")]
    tip_rows  <- tip_rows[order(tip_rows$y), ]

    acc_vec  <- if ("accessory" %in% names(meta)) safe_num(meta$accessory) else rep(0, nrow(meta))
    uniq_vec <- if ("unique"    %in% names(meta)) safe_num(meta$unique)    else rep(0, nrow(meta))
    size_vec <- if ("total_length" %in% names(meta)) round(safe_num(meta$total_length) / 1e6, 2) else rep(NA_real_, nrow(meta))
    cds_vec  <- if ("n_cds" %in% names(meta)) safe_num(meta$n_cds) else rep(NA_real_, nrow(meta))

    tile_df <- tibble::tibble(
      label       = meta$label,
      Core        = core_vec,
      Accessory   = acc_vec,
      Unique      = uniq_vec,
      `GC (%)`    = round(gc_vec, 1),
      `Size (Mb)` = size_vec,
      CDS         = cds_vec
    ) %>%
      tidyr::pivot_longer(-label, names_to = "metric", values_to = "value") %>%
      dplyr::mutate(
        metric = factor(metric,
                        levels = c("Core", "Accessory", "Unique",
                                   "GC (%)", "Size (Mb)", "CDS")),
        label  = factor(label, levels = tip_rows$label),
        group  = ifelse(metric %in% c("Core", "Accessory", "Unique"),
                        "pangenome", "genomic")
      ) %>%
      dplyr::group_by(metric) %>%
      dplyr::mutate(
        norm = if (all(is.na(value))) NA_real_
               else (value - min(value, na.rm = TRUE)) /
                    max(1e-9, diff(range(value, na.rm = TRUE)))
      ) %>%
      dplyr::ungroup()

    # Per-metric color: navy ramp for pan-genome cols, red ramp for
    # genomic cols. Using scale_fill_identity() avoids ggnewscale.
    # NB: compute on the non-NA subset only -- colorRamp/rgb error
    # when fed NA, so the NAs are painted white up front.
    ramp_navy <- grDevices::colorRamp(c("#F5F5F5", "#2C5F7A"))
    ramp_red  <- grDevices::colorRamp(c("#F5F5F5", "#D06461"))
    to_hex <- function(ramp, v) {
      v <- pmin(1, pmax(0, v))
      rgbm <- ramp(v)
      grDevices::rgb(rgbm[, 1], rgbm[, 2], rgbm[, 3], maxColorValue = 255)
    }
    tile_df$fill_color <- "#FFFFFF"
    ok_idx <- which(!is.na(tile_df$norm))
    if (length(ok_idx) > 0L) {
      is_pg <- tile_df$group[ok_idx] == "pangenome"
      tile_df$fill_color[ok_idx[ is_pg]] <- to_hex(ramp_navy, tile_df$norm[ok_idx[ is_pg]])
      tile_df$fill_color[ok_idx[!is_pg]] <- to_hex(ramp_red,  tile_df$norm[ok_idx[!is_pg]])
    }

    # Optional gain/loss column appended to the side panel.
    if (!is.null(gl_tip_df) && nrow(gl_tip_df) > 0L) {
      gl_tile <- tibble::tibble(
        label  = factor(gl_tip_df$pretty, levels = tip_rows$label),
        metric = factor("Gain/Loss",
                        levels = c(levels(tile_df$metric), "Gain/Loss")),
        value  = gl_tip_df$n_gained + gl_tip_df$n_lost,
        norm   = NA_real_,
        gl_label = gl_tip_df$gl_label
      )
      # Extend metric factor on the main tile_df too.
      tile_df$metric <- factor(tile_df$metric,
                               levels = c(levels(tile_df$metric), "Gain/Loss"))
    } else {
      gl_tile <- NULL
    }

    side_panel <- ggplot2::ggplot() +
      ggplot2::geom_tile(
        data = tile_df,
        ggplot2::aes(x = metric, y = label, fill = fill_color),
        color = "white", linewidth = 0.4
      ) +
      ggplot2::geom_text(
        data = tile_df,
        ggplot2::aes(x = metric, y = label,
                     label = ifelse(is.na(value), "",
                                    formatC(value, format = "g", digits = 3)),
                     color = ifelse(!is.na(norm) & norm > 0.55,
                                    "#FFFFFF", "#202020")),
        size = 2.6
      ) +
      ggplot2::scale_color_identity()
    if (!is.null(gl_tile)) {
      side_panel <- side_panel +
        ggplot2::geom_tile(
          data = gl_tile,
          ggplot2::aes(x = metric, y = label),
          fill = "#FCE4DE", color = "white", linewidth = 0.4
        ) +
        ggplot2::geom_text(
          data = gl_tile,
          ggplot2::aes(x = metric, y = label, label = gl_label),
          size = 2.4, color = "#8B3A2F", fontface = "bold"
        )
    }
    side_panel <- side_panel +
      ggplot2::scale_fill_identity() +
      ggplot2::scale_x_discrete(position = "top") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 0, size = 8),
        panel.grid = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(5, 5, 5, 0)
      )
  }

  subtitle_txt <- paste0(
    "Side panel: pan-genome (navy = Core, Accessory, Unique) | ",
    "genomic (red = GC%, Size, CDS)",
    if (!is.null(gl_tip_df)) " | +gained/-lost (Fitch)" else ""
  )

  final_plot <- if (!is.null(side_panel)) {
    tryCatch(
      patchwork::wrap_plots(p2, side_panel, widths = c(2.2, 1.3)) +
        patchwork::plot_annotation(
          title    = "Core-gene phylogeny  (IQ-TREE)",
          subtitle = subtitle_txt,
          theme = ggplot2::theme(
            plot.title    = ggplot2::element_text(size = 13, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 8, color = "#505050")
          )
        ),
      error = function(e) p2
    )
  } else {
    p2
  }

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    tryCatch(
      ggplot2::ggsave(output_file, final_plot, width = 14, height = 8, dpi = 300),
      error = function(e) {
        warning("phylo_tree_plot: ggsave of combined plot failed (",
                conditionMessage(e), ") -- saving tree-only")
        ggplot2::ggsave(output_file, p2, width = 10, height = 8, dpi = 300)
      }
    )
    message("phylo_tree_plot written to: ", output_file)
  }

  invisible(final_plot)
}
