#' Circular pangenome with a right-hand genome metadata panel
#'
#' Draws a 270-degree presence/absence circos plot and overlays a dedicated
#' per-genome metadata panel in the open quadrant. The circular body keeps the
#' dense presence/absence view, while the panel "unwraps" each genome track
#' into a readable row with stacked cluster counts plus GC%, genome size, and
#' CDS summaries.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level results directory.
#' @param output_file Optional PDF path.
#' @export

circos_pangenome_range_or_default <- function(x, default = c(0, 1)) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    default
  } else {
    range(x)
  }
}

circos_pangenome_palette_value <- function(x, palette, limits, na_col = "#F0F0F0") {
  if (!is.finite(x)) {
    return(na_col)
  }

  span <- diff(limits)
  if (!is.finite(span) || span <= 0) {
    return(palette[ceiling(length(palette) / 2)])
  }

  idx <- round((x - limits[1]) / span * (length(palette) - 1)) + 1L
  palette[pmax(1L, pmin(length(palette), idx))]
}

circos_pangenome_fmt_int <- function(x) {
  ifelse(
    is.finite(x),
    prettyNum(round(x), big.mark = ",", preserve.width = "none"),
    "NA"
  )
}

circos_pangenome_fmt_num <- function(x, digits = 1, suffix = "") {
  ifelse(
    is.finite(x),
    paste0(formatC(x, format = "f", digits = digits), suffix),
    "NA"
  )
}

circos_pangenome_pick_labels <- function(meta) {
  organism_short <- ifelse(
    !is.na(meta$organism) & nzchar(meta$organism),
    sub("^(\\S+\\s+\\S+).*", "\\1", meta$organism),
    meta$genome_key
  )

  strain_labels <- ifelse(
    !is.na(meta$strain) & nzchar(meta$strain),
    meta$strain,
    organism_short
  )

  missing_label <- is.na(strain_labels) | !nzchar(strain_labels)
  strain_labels[missing_label] <- meta$genome_key[missing_label]
  strain_labels
}

circos_pangenome_load_pairwise <- function(results_dir, genome_keys,
                                           filename, value_col) {
  if (is.null(results_dir)) return(NULL)
  path <- file.path(results_dir, "dnmb", "processed", filename)
  if (!file.exists(path) || !requireNamespace("arrow", quietly = TRUE)) return(NULL)
  df <- tibble::as_tibble(arrow::read_parquet(path))
  if (!all(c("genome_a", "genome_b", value_col) %in% names(df))) return(NULL)

  n <- length(genome_keys)
  mat <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(genome_keys, genome_keys))
  keep <- df$genome_a %in% genome_keys & df$genome_b %in% genome_keys
  df <- df[keep, , drop = FALSE]
  if (!nrow(df)) return(NULL)
  idx_a <- match(df$genome_a, genome_keys)
  idx_b <- match(df$genome_b, genome_keys)
  mat[cbind(idx_a, idx_b)] <- df[[value_col]]
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (is.na(mat[i, j]) && !is.na(mat[j, i])) {
        mat[i, j] <- mat[j, i]
      }
    }
    if (is.na(mat[i, i])) mat[i, i] <- 100
  }
  mat
}

circos_pangenome_load_ani <- function(results_dir, genome_keys) {
  circos_pangenome_load_pairwise(results_dir, genome_keys,
                                 "ani_matrix.parquet",  "ani_percent")
}
circos_pangenome_load_pocp <- function(results_dir, genome_keys) {
  circos_pangenome_load_pairwise(results_dir, genome_keys,
                                 "pocp_matrix.parquet", "pocp_percent")
}

#' Authoritative genome ordering for the circos + metadata panel
#'
#' One place that decides the left-to-right order of panel columns and the
#' inner-to-outer order of circos ring tracks. Every downstream routine —
#' `prepare_data`, `build_linkage`, ring drawing, panel drawing — must read
#' from this ordering so that the heatmap, dendrogram, ring tracks, and
#' metadata bars all reference the same genome sequence.
#'
#' Strategy: prefer ANI hclust order (gives the dendrogram for free and
#' clusters similar genomes together), fall back to the core-gene newick,
#' finally to input order.
#' @noRd
circos_pangenome_build_ordering <- function(dnmb, results_dir = NULL) {
  all_keys <- dnmb$genome_meta %>%
    dplyr::arrange(genome_uid) %>%
    dplyr::pull(genome_key)

  ani_mat_raw <- circos_pangenome_load_ani(results_dir, all_keys)
  hclust_obj <- NULL
  order_source <- "input"
  genome_keys <- all_keys

  if (!is.null(ani_mat_raw) && nrow(ani_mat_raw) >= 2) {
    d <- stats::as.dist(100 - ani_mat_raw)
    hc <- tryCatch(stats::hclust(d, method = "average"),
                   error = function(e) NULL)
    if (!is.null(hc)) {
      hclust_obj <- hc
      genome_keys <- hc$labels[hc$order]
      order_source <- "ani"
    }
  }

  if (is.null(hclust_obj)) {
    tree_path <- if (!is.null(results_dir)) {
      file.path(results_dir, "dnmb", "processed", "phylo_tree.nwk")
    } else {
      ""
    }
    if (file.exists(tree_path) && requireNamespace("ape", quietly = TRUE)) {
      tree <- ape::ladderize(ape::read.tree(tree_path))
      ordered <- tree$tip.label[tree$tip.label %in% all_keys]
      genome_keys <- c(ordered, setdiff(all_keys, ordered))
      order_source <- "tree"
    }
  }

  ani_mat <- if (!is.null(ani_mat_raw)) {
    ani_mat_raw[genome_keys, genome_keys, drop = FALSE]
  } else {
    NULL
  }

  list(
    genome_keys = genome_keys,
    n = length(genome_keys),
    hclust_obj = hclust_obj,
    order_source = order_source,
    ani_mat = ani_mat
  )
}

#' Shared geometry spec coupling panel columns ↔ ring tracks
#'
#' Forces the panel's per-column width to equal the ring's radial thickness,
#' and anchors the panel's left edge at the innermost ring radius so that
#' each panel column sits at the exact same x-coordinate as the matching
#' ring's end-point on the positive x-axis (angle 0° — the track terminus
#' that faces the open quadrant).
#'
#' Ring draw order is also spelled out here: leftmost panel column (ordering
#' position 1) maps to the INNERMOST genome ring. The main circos loop iterates
#' `track_draw_order` which reverses ordering position so tracks stack
#' outside-in while still respecting the panel ↔ ring alignment.
#' @noRd
circos_pangenome_build_linkage <- function(ordering,
                                           ring_budget = 0.55,
                                           cluster_track_h = 0.022,
                                           track_gap_frac = 0.14,
                                           cell_cap = 0.060,
                                           cell_floor = 0.010,
                                           summary_n = 0L,
                                           summary_n_outer = NULL,
                                           summary_n_inner = 0L,
                                           summary_frac = 0.6) {
  n <- ordering$n
  summary_n <- as.integer(summary_n)
  # Back-compat: if only `summary_n` was supplied, assume everything goes
  # OUTSIDE the genome rings (the original behaviour).
  if (is.null(summary_n_outer)) {
    summary_n_outer <- summary_n - as.integer(summary_n_inner)
  }
  summary_n_outer <- as.integer(summary_n_outer)
  summary_n_inner <- as.integer(summary_n_inner)
  # Per-ring "slot" = drawn ring thickness + two margins (one above, one
  # below). Fixing the column pitch = slot width keeps panel columns
  # aligned with ring midpoints even though there is now a visible gap
  # between tracks.
  #
  # Summary tracks (outside the cluster header) share the same margin but
  # are drawn thinner: summary_h = cell_w * summary_frac. Radial budget:
  #   ring_budget = cluster_track_h + n*cell_w + summary_n*summary_h
  #                 + 2m*(summary_n + 1)
  # with 2m = cell_w * track_gap_frac. Solving for cell_w:
  cell_w <- max(
    cell_floor,
    min(cell_cap,
        (ring_budget - cluster_track_h) /
          (n + summary_n * summary_frac + track_gap_frac * (1 + summary_n)))
  )
  genome_track_h  <- cell_w * (1 - track_gap_frac)       # drawn thickness
  track_margin    <- cell_w * track_gap_frac / 2         # per-side gap
  summary_track_h <- cell_w * summary_frac               # outer-track thickness

  # Anchor: innermost ring's radial bottom at angle 0° is
  #   1 - summary_n*summary_h - 2m*(summary_n + 1) - cluster_track_h - n*cell_w
  # Setting x_block_start to this value keeps col_center[i] ≡ ring i's
  # r_mid at angle 0°, regardless of how many outer summary tracks exist.
  x_block_start <- 1 - summary_n * summary_track_h -
    2 * track_margin * (summary_n + 1) -
    cluster_track_h - n * cell_w

  # Track ordering in circlize is outermost-first (index 1 = outermost):
  #   1..summary_n_outer                        : OUTER summary tracks
  #   summary_n_outer + 1..summary_n_outer + n  : per-genome rings
  #                                                 outermost → innermost
  #   summary_n_outer + n + 1 ..
  #     summary_n_outer + n + summary_n_inner   : INNER summary tracks
  #                                                 (radially adjacent to
  #                                                  the cluster header)
  #   summary_n_outer + n + summary_n_inner + 1 : cluster header (INNERMOST)
  # Leftmost panel column (position 1) → innermost genome ring
  #   (track summary_n_outer + n); rightmost → outermost
  #   (track summary_n_outer + 1).
  track_index_for_panel <- (summary_n_outer + n + 1L) - seq_len(n)
  track_draw_order      <- rev(seq_len(n))               # draws outside-in

  list(
    n                = n,
    cell_w           = cell_w,
    genome_track_h   = genome_track_h,
    track_margin     = track_margin,
    track_gap_frac   = track_gap_frac,
    cluster_track_h  = cluster_track_h,
    summary_n        = summary_n,
    summary_n_outer  = summary_n_outer,
    summary_n_inner  = summary_n_inner,
    summary_track_h  = summary_track_h,
    summary_frac     = summary_frac,
    ring_budget      = ring_budget,
    x_block_start    = x_block_start,
    x_block_end      = x_block_start + n * cell_w,
    col_left         = x_block_start + (seq_len(n) - 1) * cell_w,
    col_right        = x_block_start + seq_len(n) * cell_w,
    col_center       = x_block_start + (seq_len(n) - 0.5) * cell_w,
    track_index      = track_index_for_panel,
    track_draw_order = track_draw_order
  )
}

#' Query ring radial bounds for each panel column (after circos.initialize)
#' @noRd
circos_pangenome_linkage_radii <- function(linkage, sector = "pan") {
  n <- linkage$n
  r_bot <- r_top <- numeric(n)
  for (i in seq_len(n)) {
    ti <- linkage$track_index[i]
    r_bot[i] <- circlize::get.cell.meta.data(
      "cell.bottom.radius", sector.index = sector, track.index = ti)
    r_top[i] <- circlize::get.cell.meta.data(
      "cell.top.radius",    sector.index = sector, track.index = ti)
  }
  list(r_bot = r_bot, r_top = r_top, r_mid = (r_bot + r_top) / 2)
}

circos_pangenome_prepare_data <- function(dnmb, results_dir = NULL,
                                          ordering = NULL) {
  presence <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(
      dnmb$genome_meta %>% dplyr::select(genome_uid, genome_key),
      by = "genome_uid"
    ) %>%
    dplyr::filter(!is.na(genome_key))

  if (is.null(ordering)) {
    ordering <- circos_pangenome_build_ordering(dnmb, results_dir)
  }
  genome_keys  <- ordering$genome_keys
  hclust_obj   <- ordering$hclust_obj
  order_source <- ordering$order_source
  ani_mat      <- ordering$ani_mat
  pocp_mat     <- circos_pangenome_load_pocp(results_dir, genome_keys)
  if (!is.null(pocp_mat)) {
    pocp_mat <- pocp_mat[genome_keys, genome_keys, drop = FALSE]
  }

  n_per <- presence %>%
    dplyr::count(cluster_id, name = "ng") %>%
    dplyr::arrange(dplyr::desc(ng), cluster_id)

  sorted_cids <- n_per$cluster_id
  n_clusters <- length(sorted_cids)
  n_genomes <- length(genome_keys)

  if (!n_clusters || !n_genomes) {
    return(NULL)
  }

  # Per-cluster outer-track stats. `dnmb$clusters` is one row per gene, so
  # count(cluster_id, genome_uid) gives paralog multiplicity; summing rows
  # per cluster gives n_genes; rows per cluster gives n_contrib; the max of
  # the per-(cluster, genome) counts gives max_paralogs.
  cluster_stats_raw <- dnmb$clusters %>%
    dplyr::count(cluster_id, genome_uid, name = "paralogs") %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(
      n_genes      = sum(paralogs),
      n_contrib    = dplyr::n(),
      max_paralogs = max(paralogs),
      .groups      = "drop"
    )
  cluster_stats <- tibble::tibble(cluster_id = sorted_cids) %>%
    dplyr::left_join(cluster_stats_raw, by = "cluster_id") %>%
    dplyr::mutate(
      n_genes      = dplyr::coalesce(n_genes, 0L),
      n_contrib    = dplyr::coalesce(n_contrib, 0L),
      max_paralogs = dplyr::coalesce(max_paralogs, 0L)
    )

  presence$cluster_id <- factor(presence$cluster_id, levels = sorted_cids)
  presence$genome_key <- factor(presence$genome_key, levels = genome_keys)

  bin_mat <- stats::xtabs(~ cluster_id + genome_key, data = presence)
  bin_mat[bin_mat > 0] <- 1L
  bin_mat <- unclass(bin_mat)
  storage.mode(bin_mat) <- "integer"

  meta <- dnmb$genome_meta[match(genome_keys, dnmb$genome_meta$genome_key), , drop = FALSE]

  per_genome <- compute_per_genome_stats(dnmb) %>%
    dplyr::select(genome_key, core, accessory, unique) %>%
    dplyr::mutate(
      core = dplyr::coalesce(core, 0L),
      accessory = dplyr::coalesce(accessory, 0L),
      unique = dplyr::coalesce(unique, 0L)
    )

  panel_df <- meta %>%
    dplyr::left_join(per_genome, by = "genome_key") %>%
    dplyr::mutate(
      label = circos_pangenome_pick_labels(.),
      size_mb = total_length / 1e6,
      core = dplyr::coalesce(core, 0L),
      accessory = dplyr::coalesce(accessory, 0L),
      unique = dplyr::coalesce(unique, 0L),
      total_clusters = core + accessory + unique
    )

  list(
    bin_mat       = bin_mat,
    n_per         = n_per,
    n_clusters    = n_clusters,
    n_genomes     = n_genomes,
    genome_keys   = genome_keys,
    panel_df      = panel_df,
    cluster_stats = cluster_stats,
    ani_mat       = ani_mat,
    pocp_mat      = pocp_mat,
    hclust_obj    = hclust_obj,
    order_source  = order_source,
    ordering      = ordering
  )
}

#' Shared genus/species colour builder
#'
#' Reuses the exact rule from [ani_pocp_heatmap()]: genus = first token of
#' organism, species = second token (fallback "sp."). The base palettes are
#' identical so the circos track colouring, the right-panel strips, and the
#' ANI/POCP heatmap annotations line up across plots.
#' @noRd
circos_pangenome_species_palette <- function(panel_df) {
  org <- panel_df$organism
  genus <- ifelse(is.na(org) | !nzchar(org), "Unknown",
                  sub("^(\\S+).*", "\\1", org))
  species <- ifelse(is.na(org) | !nzchar(org), "sp.",
                    sub("^\\S+\\s+(\\S+).*", "\\1", org))
  species <- ifelse(species == genus, "sp.", species)

  # Calm earth/sage/slate family chosen to harmonise with the ANI gradient
  # (teal → pastel blue → amber → orange → red) and the cluster/metric
  # fills (#1F4E5F core, #F0A35B accessory, #D06461 unique, #EF6C00 CDS,
  # #1565C0 Mb, #2E7D32 GC%). Species tones are medium-light so each
  # ring's "present" cell stays legible against the grey "absent"
  # background; genus tones share the same hue family at lower lightness
  # so the two strips read as a hierarchy, not two unrelated palettes.
  species_base <- c(
    "#8FB2AD",  # sage teal
    "#D2B188",  # warm sand
    "#C39287",  # terracotta
    "#9AAEC7",  # slate blue
    "#AEA5BD",  # dusty mauve
    "#BBA78A",  # taupe
    "#A1B89E",  # moss
    "#C79E8F",  # clay
    "#97AFB4",  # steel teal
    "#C9BD90",  # honey
    "#B3A1B8",  # pale heather
    "#A0B49F",  # pale olive
    "#D6C39D",  # linen
    "#A6B3AB",  # cool stone
    "#C9A59A",  # dusty coral
    "#B5A99A"   # warm grey
  )
  genus_base <- c(
    "#4E7C77",  # forest teal
    "#9C7A4F",  # bronze
    "#8F554A",  # burnt sienna
    "#556F92",  # twilight blue
    "#755A7B",  # plum
    "#6D7549",  # olive
    "#9A604A",  # rust
    "#54606B"   # graphite
  )

  uniq_sp <- sort(unique(species))
  uniq_gn <- sort(unique(genus))
  species_cols <- setNames(rep_len(species_base, length(uniq_sp)), uniq_sp)
  genus_cols   <- setNames(rep_len(genus_base,   length(uniq_gn)), uniq_gn)

  list(
    species = species,
    genus = genus,
    species_cols = species_cols,
    genus_cols = genus_cols,
    # Present-cell colour for each genome (ordered like panel_df rows)
    present_cols = unname(species_cols[species])
  )
}

circos_pangenome_track_radii <- function(n_genomes, sector = "pan") {
  # Track 1 is the cluster header ring; tracks 2..(n_genomes+1) are per-genome.
  radii <- vector("list", n_genomes)
  for (g in seq_len(n_genomes)) {
    ti <- g + 1L
    r_top <- circlize::get.cell.meta.data("cell.top.radius", sector.index = sector, track.index = ti)
    r_bot <- circlize::get.cell.meta.data("cell.bottom.radius", sector.index = sector, track.index = ti)
    radii[[g]] <- c(r_bot = r_bot, r_top = r_top)
  }
  radii
}

circos_pangenome_draw_dendrogram <- function(hc, x_start, cell_w, y_bot, y_top,
                                             col = "#3E4954", lwd = 0.7) {
  if (is.null(hc)) return(invisible(NULL))
  n <- length(hc$order)
  if (n < 2) return(invisible(NULL))
  heights <- hc$height
  max_h <- max(heights, na.rm = TRUE)
  if (!is.finite(max_h) || max_h <= 0) max_h <- 1

  # pos[k] = leaf k's left-to-right position (1..n) in the dendrogram.
  pos <- integer(n)
  pos[hc$order] <- seq_len(n)
  leaf_x <- function(k) x_start + (pos[k] - 0.5) * cell_w
  h_to_y <- function(h) y_bot + (y_top - y_bot) * (h / max_h)

  node_x <- numeric(n - 1)
  for (i in seq_len(n - 1)) {
    a <- hc$merge[i, 1]
    b <- hc$merge[i, 2]
    xa <- if (a < 0) leaf_x(-a) else node_x[a]
    xb <- if (b < 0) leaf_x(-b) else node_x[b]
    ya <- if (a < 0) y_bot      else h_to_y(heights[a])
    yb <- if (b < 0) y_bot      else h_to_y(heights[b])
    ym <- h_to_y(heights[i])

    graphics::segments(xa, ya, xa, ym, col = col, lwd = lwd)
    graphics::segments(xb, yb, xb, ym, col = col, lwd = lwd)
    graphics::segments(xa, ym, xb, ym, col = col, lwd = lwd)
    node_x[i] <- (xa + xb) / 2
  }
  invisible(NULL)
}

#' Default spec for panel metric rows
#'
#' Each entry is a self-describing row that the module layout engine can
#' stack without hardcoding. Callers can extend this list (or pass a
#' custom list via the `panel_rows` argument on `circos_pangenome_*`) to
#' add new bar-type metrics; the module frame, row positions, and row
#' labels all auto-reflow to fit.
#'
#' Row schema:
#'   key        : unique ID (used for logging / debugging)
#'   label      : character vector — one element for a single-line row
#'                label, two elements for a two-line label (e.g.
#'                c("Gene clusters", "(No. of CDS)"))
#'   type       : "bar" (single-value per genome, min-max normalized) or
#'                "stacked" (e.g. core / accessory / unique segments)
#'   fill       : single colour (type="bar") OR named list of colours
#'                (type="stacked")
#'   values     : numeric vector of per-genome values (type="bar") OR
#'                a named list of per-genome vectors (type="stacked",
#'                one entry per segment, names matching `fill`)
#'   total      : per-genome totals used to scale a stacked bar
#'                (type="stacked" only; ignored for "bar")
#'   fmt        : value formatter function (v) -> character
#'   h_frac     : row height as multiple of `cell_w`
#'   axis_type  : "minmax" (bars range from a floor fraction up to vmax;
#'                right-hand axis labels vmin/vmax) or "proportional"
#'                (stacked bars sized as value / profile_max; right-hand
#'                axis labels 0/max)
#' @noRd
circos_pangenome_default_panel_rows <- function(panel_df,
                                                col_core,
                                                col_accessory,
                                                col_unique) {
  list(
    list(
      key       = "gene_clusters",
      label     = c("Gene clusters", "(No. of CDSs)"),
      type      = "stacked",
      fill      = list(core       = col_core,
                       accessory  = col_accessory,
                       unique     = col_unique),
      values    = list(core       = panel_df$core,
                       accessory  = panel_df$accessory,
                       unique     = panel_df$unique),
      total     = panel_df$total_clusters,
      fmt       = function(v) circos_pangenome_fmt_int(v),
      h_frac    = 2.10,
      axis_type = "proportional"
    ),
    list(
      key       = "size_mb",
      label     = "Mb",
      type      = "bar",
      fill      = "#1565C0",
      values    = panel_df$size_mb,
      fmt       = function(v) circos_pangenome_fmt_num(v, 2),
      h_frac    = 1.10,
      axis_type = "minmax"
    ),
    list(
      key       = "gc_percent",
      label     = "GC%",
      type      = "bar",
      fill      = "#2E7D32",
      values    = panel_df$gc_percent,
      fmt       = function(v) circos_pangenome_fmt_num(v, 1),
      h_frac    = 1.10,
      axis_type = "minmax"
    )
  )
}

#' Total module height for a given panel_rows spec
#' @noRd
circos_pangenome_panel_rows_h <- function(panel_rows, cell_w, metric_gap_y,
                                           module_pad) {
  n_rows <- length(panel_rows)
  if (!n_rows) return(0)
  row_h <- vapply(panel_rows, function(r) r$h_frac, numeric(1)) * cell_w
  sum(row_h) + max(0L, n_rows - 1L) * metric_gap_y + 2 * module_pad
}

circos_pangenome_panel_heights <- function(cell_w, n, has_ani,
                                            has_pocp = has_ani,
                                            summary_total_h = NULL) {
  # Legend slot holds one gradient strip + tick labels when only ANI is
  # present, or two stacked strips (ANI + POCP) with ticks when both
  # exist. Each strip is ~0.38 cell_w tall with a 0.08 gap + tick label
  # row (~0.45 cell_w total).
  h_ani_legend <- if (!has_ani) 0 else if (has_pocp) cell_w * 1.35 else cell_w * 0.75
  # Dendrogram height matches the outer summary-track stack's radial
  # height when those tracks exist — visually ties the top of the panel
  # to the outermost rings on the left. Falls back to a cell_w * 4.0
  # default when there are no summary tracks.
  h_dendro <- if (!is.null(summary_total_h) && is.finite(summary_total_h) &&
                  summary_total_h > 0) {
    summary_total_h
  } else {
    cell_w * 4.0
  }
  list(
    h_lab        = cell_w * 3.5,
    h_strip      = cell_w * 0.45,   # 50% of prior 0.9
    h_dendro     = h_dendro,
    # Gene-clusters bar keeps the taller profile because it stacks three
    # segments (core/accessory/unique) and needs room for in-bar counts.
    # GC/Mb bars carry one value each, so they use the shorter
    # h_bar_metric; the heatmap above auto-shifts down accordingly. The
    # ratio (gene clusters : metric) = ~1.9 reads as a clear hierarchy
    # inside the module box (gene clusters is the "anchor" row, Mb + GC
    # are supporting metrics).
    h_bar        = cell_w * 2.10,
    h_bar_metric = cell_w * 1.10,
    h_ani        = if (has_ani) cell_w * n else 0,
    h_ani_legend = h_ani_legend,
    gap_y        = cell_w * 0.3,
    # Tight gap used only between the metric bar rows inside the module
    # (Gene clusters / Mb / GC%). Other section boundaries keep the
    # wider gap_y so the dendro/strip/ANI/label groupings stay visually
    # distinct from the bar-module.
    metric_gap_y = cell_w * 0.45,
    # Padding inside the module-enclosing border (top/bottom + left/right)
    # so the inner bars breathe within the frame.
    module_pad   = cell_w * 0.18
  )
}

#' Total y-extent of the metadata panel (so the main fn can size the canvas)
#' @noRd
circos_pangenome_panel_total_h <- function(cell_w, n, has_ani,
                                            has_pocp = has_ani,
                                            summary_total_h = NULL) {
  h <- circos_pangenome_panel_heights(cell_w, n, has_ani, has_pocp,
                                      summary_total_h = summary_total_h)
  # Row stack (bottom→top): module { gene-clusters (CDS) | Mb | GC% },
  # [ANI legend, ANI matrix], genus, species, dendro, label. The three
  # bar rows share one enclosing border ("module"), so their contribution
  # adds h_bar + 2*h_bar_metric + 2*metric_gap_y + 2*module_pad.
  legend_gap <- if (has_ani) 0.2 * h$gap_y else 0
  module_h <- h$h_bar + 2 * h$h_bar_metric +
    2 * h$metric_gap_y + 2 * h$module_pad
  # Outer gap breakdown (outside the module):
  #   3 * gap_y        : module→ANI_legend, ANI_matrix→genus, dendro→label
  #   0.5 * gap_y      : species→dendro
  #   legend_gap       : ANI_legend→ANI_matrix
  #   0.002            : genus→species
  module_h +
    h$h_ani + h$h_ani_legend + 2 * h$h_strip +
    h$h_dendro + h$h_lab +
    3.5 * h$gap_y + legend_gap + 0.002
}

circos_pangenome_draw_aligned_panel <- function(panel_df,
                                                ani_mat,
                                                pocp_mat = NULL,
                                                hclust_obj = NULL,
                                                order_source = "input",
                                                col_core,
                                                col_accessory,
                                                col_unique,
                                                ani_pal,
                                                ani_limits,
                                                pocp_pal = NULL,
                                                pocp_limits = NULL,
                                                diag_col = "#1B2A35",
                                                species_info = NULL,
                                                linkage,
                                                panel_rows = NULL,
                                                y_bottom = 0) {
  n <- nrow(panel_df)
  if (!n) {
    return(invisible(NULL))
  }

  cell_w        <- linkage$cell_w
  x_block_start <- linkage$x_block_start
  x_block_end   <- linkage$x_block_end
  block_w       <- x_block_end - x_block_start
  track_margin  <- linkage$track_margin
  summary_n     <- linkage$summary_n

  if (is.null(panel_rows)) {
    panel_rows <- circos_pangenome_default_panel_rows(
      panel_df, col_core, col_accessory, col_unique)
  }
  n_module_rows <- length(panel_rows)

  # ── Heights ──────────────────────────────────────────────────────
  summary_total_h_local <- if (!is.null(summary_n) && summary_n > 0) {
    summary_n * (linkage$summary_track_h + 2 * track_margin)
  } else {
    NULL
  }
  heights      <- circos_pangenome_panel_heights(
    cell_w, n, !is.null(ani_mat), !is.null(pocp_mat),
    summary_total_h = summary_total_h_local
  )
  h_lab        <- heights$h_lab
  h_strip      <- heights$h_strip
  h_ani_legend <- heights$h_ani_legend
  gap_y        <- heights$gap_y
  metric_gap_y <- heights$metric_gap_y
  module_pad   <- heights$module_pad
  has_ani      <- !is.null(ani_mat)
  legend_gap   <- if (has_ani) 0.2 * gap_y else 0

  # Module height auto-derives from panel_rows spec; adding or removing
  # rows reflows the module without touching layout math.
  module_h <- circos_pangenome_panel_rows_h(
    panel_rows, cell_w, metric_gap_y, module_pad)

  # ── Upper block Y — anchored to the circos track radii so each panel
  # element's vertical extent EXACTLY matches its corresponding
  # horizontal track at angle 90° apex (top of circle). This is what
  # gives the heatmap/strips/dendrogram their "continuation of the
  # tracks" visual on the panel side.
  ring_radii <- circos_pangenome_linkage_radii(linkage, "pan")
  # Cluster header is the INNERMOST track. Above it sit (radially outward):
  # the inner summary tracks, the genome rings, and finally the outer
  # summary tracks. cluster_ti = outer + n + inner + 1.
  cluster_ti <- linkage$summary_n_outer + linkage$n +
    linkage$summary_n_inner + 1L
  r_cluster_top <- circlize::get.cell.meta.data(
    "cell.top.radius",    sector.index = "pan",
    track.index = cluster_ti)
  r_cluster_bot <- circlize::get.cell.meta.data(
    "cell.bottom.radius", sector.index = "pan",
    track.index = cluster_ti)
  if (!is.null(summary_n) && summary_n > 0) {
    r_summary_top <- circlize::get.cell.meta.data(
      "cell.top.radius",    sector.index = "pan",
      track.index = 1L)
    r_summary_bot <- circlize::get.cell.meta.data(
      "cell.bottom.radius", sector.index = "pan",
      track.index = summary_n)
  } else {
    r_summary_top <- NA_real_
    r_summary_bot <- NA_real_
  }

  # Heatmap Y — each row uses its ring's exact drawn extent (r_bot..r_top)
  # so row h's vertical span matches ring h's visible thickness at apex.
  # Gaps between rows == gaps between rings (no margin extension).
  y_ani_bot <- ring_radii$r_bot[1]
  y_ani_top <- ring_radii$r_top[n]

  # Dendrogram spans the outer summary tracks' radial extent (exact).
  # Species/Genus info is now presented via the 12-o'clock thin strip at
  # the circle apex — no duplicate strip in the panel block.
  y_dendro_bot <- if (!is.na(r_summary_bot)) r_summary_bot else y_ani_top + gap_y
  y_dendro_top <- if (!is.na(r_summary_top)) r_summary_top else y_dendro_bot + cell_w * 4.0

  y_lab_bot <- y_dendro_top + gap_y * 0.6
  y_lab_top <- y_lab_bot + h_lab

  # ── Lower block Y — modules stack independently below the heatmap.
  # Each module (scale, then panel_rows bars) carves out its own slot
  # from `y_cursor`, top→bottom, so adding/removing modules just moves
  # the cursor without re-deriving the rest. Gaps are tightened so the
  # module sits flush against the scale and has maximum headroom to
  # expand into the space above `y_bottom`.
  y_cursor <- y_ani_bot - gap_y * 0.6
  if (has_ani) {
    y_ani_leg_top <- y_cursor
    y_ani_leg_bot <- y_ani_leg_top - h_ani_legend
    y_cursor      <- y_ani_leg_bot - gap_y * 0.35
  } else {
    y_ani_leg_top <- y_ani_leg_bot <- NA_real_
  }
  y_module_top <- y_cursor
  # Interactive proportional-fill: the module EXACTLY fills the
  # vertical budget between the cursor (below the scale) and the
  # panel's lower boundary `y_bottom`. Row heights, inter-row gaps
  # and inner padding scale by the SAME factor so the h_frac ratios
  # declared in `panel_rows` are preserved. Works for any row count
  # without hardcoding — add a 4th row and the three existing rows
  # auto-shrink; drop one and the remaining rows auto-expand.
  available_h <- y_module_top - y_bottom
  fill_scale  <- if (is.finite(available_h) && available_h > 0 &&
                     module_h > 0) {
    available_h / module_h
  } else {
    1
  }
  module_h_eff   <- module_h * fill_scale
  metric_gap_eff <- metric_gap_y * fill_scale
  module_pad_eff <- module_pad   * fill_scale
  y_module_bot   <- y_module_top - module_h_eff

  # Per-row Y inside the module (bottom→top, following panel_rows order).
  # Each row's height = h_frac * cell_w * fill_scale; proportions between
  # rows (e.g. gene_clusters : mb : gc = 2.10 : 1.10 : 1.10) are preserved.
  row_ypos <- vector("list", n_module_rows)
  y_row <- y_module_bot + module_pad_eff
  for (i in seq_len(n_module_rows)) {
    row_h <- panel_rows[[i]]$h_frac * cell_w * fill_scale
    row_ypos[[i]] <- c(bot = y_row, top = y_row + row_h)
    y_row <- y_row + row_h + (if (i < n_module_rows) metric_gap_eff else 0)
  }

  # Left-side row labels
  lab_x <- x_block_start - cell_w * 0.18
  lab_cex <- max(0.45, min(0.80, 14 / (n + 8)))
  draw_row_label <- function(y_mid, text, cex = lab_cex) {
    graphics::text(lab_x, y_mid, labels = text,
                   adj = c(1, 0.5), cex = cex, font = 2, col = "#26303A")
  }
  if (!is.null(hclust_obj)) {
    draw_row_label((y_dendro_top + y_dendro_bot) / 2, "ANI tree",
                   cex = lab_cex * 0.95)
  }
  if (has_ani) {
    ani_mid_y <- (y_ani_top + y_ani_bot) / 2
    # Vertical ANI (left) and POCP (right) labels flanking the heatmap.
    # ANI reads bottom→top (srt=90); POCP reads top→bottom (srt=270) so
    # the user can rotate their head the SAME way for both (head tilts
    # right to read the right-side POCP, left for the left-side ANI).
    graphics::text(x_block_start - cell_w * 0.30, ani_mid_y,
                   "ANI", srt = 90, adj = c(0.5, 0.5),
                   cex = lab_cex * 1.15, font = 2, col = "#26303A")
    if (!is.null(pocp_mat)) {
      graphics::text(x_block_end + cell_w * 0.30, ani_mid_y,
                     "POCP", srt = 270, adj = c(0.5, 0.5),
                     cex = lab_cex * 1.15, font = 2, col = "#26303A")
    }
    # "Scale" label is rendered rotated 90° in the gap between the ANI
    # and POCP gradient strips (see gradient section below) — no
    # horizontal row-label here.
  }

  # Per-row labels inside the module, rotated to read bottom→top
  # (srt = 90). Two-line specs are rendered in a SINGLE text() call
  # with "\n" so R lays them out as two parallel lines — both bold,
  # both sharing the same y anchor but offset perpendicular to the
  # reading direction (native multi-line handling). Font size is
  # UNIFIED across rows and auto-shrunk via strwidth so the longest
  # rotated label fits inside the tightest row's vertical span —
  # adjacent module boxes never get straddled.
  label_strings <- vapply(panel_rows, function(r) {
    if (length(r$label) == 1) r$label[1] else paste(r$label, collapse = "\n")
  }, character(1))
  # strwidth of the LONGEST individual line (pre-rotation width == post-
  # rotation vertical span). Computed at cex=1, scaled later.
  max_line_w <- max(vapply(panel_rows, function(r) {
    max(graphics::strwidth(r$label, cex = 1, units = "user"))
  }, numeric(1)))
  min_row_h <- min(vapply(seq_along(panel_rows), function(i) {
    row_ypos[[i]]["top"] - row_ypos[[i]]["bot"]
  }, numeric(1)))
  fit_cex <- if (max_line_w > 0) (min_row_h * 0.88) / max_line_w else lab_cex
  panel_lab_cex <- max(0.36, min(lab_cex, fit_cex))
  x_lab <- x_block_start - cell_w * 0.42
  for (i in seq_len(n_module_rows)) {
    y_mid <- (row_ypos[[i]]["top"] + row_ypos[[i]]["bot"]) / 2
    graphics::text(x_lab, y_mid, labels = label_strings[i],
                   srt = 90, adj = c(0.5, 0.5),
                   cex = panel_lab_cex, font = 2, col = "#26303A")
  }

  # Dendrogram
  if (!is.null(hclust_obj)) {
    circos_pangenome_draw_dendrogram(
      hclust_obj, x_block_start, cell_w, y_dendro_bot, y_dendro_top
    )
  }

  # Bar horizontal padding ≡ track margin so bar width equals the ring's
  # drawn radial thickness.
  bar_pad_x <- track_margin

  # ── Right-side scale axis helper ─────────────────────────────────
  axis_cex <- max(0.36, min(0.52, 6 / (n + 4)))
  draw_bar_axis <- function(y_bot, y_top, vmin, vmax, label_fmt,
                            axis_origin = 0, axis_top = 1) {
    if (!is.finite(vmin) || !is.finite(vmax)) return(invisible(NULL))
    h_avail  <- y_top - y_bot
    axis_x   <- x_block_end + cell_w * 0.18
    tick_len <- cell_w * 0.28
    y_lo <- y_bot + h_avail * axis_origin
    y_hi <- y_bot + h_avail * axis_top
    graphics::segments(axis_x, y_lo, axis_x, y_hi,
                       col = "#9AA5AF", lwd = 0.6)
    graphics::segments(axis_x, y_lo, axis_x + tick_len, y_lo,
                       col = "#9AA5AF", lwd = 0.6)
    graphics::segments(axis_x, y_hi, axis_x + tick_len, y_hi,
                       col = "#9AA5AF", lwd = 0.6)
    graphics::text(axis_x + tick_len + cell_w * 0.10, y_lo,
                   label_fmt(vmin), adj = c(0, 0.5),
                   cex = axis_cex, col = "#3E4954")
    graphics::text(axis_x + tick_len + cell_w * 0.10, y_hi,
                   label_fmt(vmax), adj = c(0, 0.5),
                   cex = axis_cex, col = "#3E4954")
  }

  # ── Bar row renderers — one per `type` in panel_rows. Each takes the
  # row spec + its Y slot and draws the bars + axis.
  draw_minmax_row <- function(row, y_bot, y_top) {
    vals <- suppressWarnings(as.numeric(row$values))
    finite_vals <- vals[is.finite(vals)]
    if (!length(finite_vals)) return(invisible(NULL))
    vmin <- min(finite_vals)
    vmax <- max(finite_vals)
    h_avail <- y_top - y_bot
    frac_min <- 0.45
    frac_max <- 1.00
    for (g in seq_len(n)) {
      v <- vals[g]
      if (!is.finite(v)) next
      frac <- if (vmax > vmin) {
        frac_min + (frac_max - frac_min) * (v - vmin) / (vmax - vmin)
      } else {
        (frac_min + frac_max) / 2
      }
      x_l <- x_block_start + (g - 1) * cell_w + bar_pad_x
      x_r <- x_block_start + g * cell_w - bar_pad_x
      bar_top_y <- y_bot + h_avail * frac
      graphics::rect(x_l, y_bot, x_r, bar_top_y,
                     col = row$fill, border = "white", lwd = 0.3)
      if (!is.null(row$fmt) && cell_w > 0.028 && n <= 20) {
        graphics::text((x_l + x_r) / 2, bar_top_y - cell_w * 0.08,
                       row$fmt(v), adj = c(1, 0.5), srt = 90,
                       cex = max(0.32, min(0.46, 5 / (n + 4))),
                       col = "white", font = 2)
      }
    }
    if (!is.null(row$fmt) && vmax > vmin) {
      draw_bar_axis(y_bot, y_top, vmin, vmax, row$fmt,
                    axis_origin = frac_min, axis_top = frac_max)
    }
  }

  draw_stacked_row <- function(row, y_bot, y_top) {
    total_vec <- suppressWarnings(as.numeric(row$total))
    profile_max <- max(total_vec, na.rm = TRUE)
    if (!is.finite(profile_max) || profile_max <= 0) profile_max <- 1
    h_avail <- y_top - y_bot
    seg_names <- names(row$fill)
    seg_cex   <- max(0.30, min(0.44, 5 / (n + 4)))
    for (g in seq_len(n)) {
      total_g <- total_vec[g]
      x_l <- x_block_start + (g - 1) * cell_w + bar_pad_x
      x_r <- x_block_start + g * cell_w - bar_pad_x
      if (!is.finite(total_g) || total_g <= 0) next
      bar_h <- h_avail * total_g / profile_max
      seg_vals <- vapply(seg_names,
                         function(nm) row$values[[nm]][g],
                         numeric(1))
      seg_h <- bar_h * seg_vals / total_g
      y_cursor <- y_bot
      for (j in seq_along(seg_h)) {
        if (seg_h[j] <= 0 || !is.finite(seg_h[j])) next
        graphics::rect(x_l, y_cursor, x_r, y_cursor + seg_h[j],
                       col = row$fill[[seg_names[j]]], border = NA)
        if (seg_h[j] > cell_w * 0.45 && cell_w > 0.028 && n <= 20) {
          graphics::text((x_l + x_r) / 2, y_cursor + seg_h[j] / 2,
                         circos_pangenome_fmt_int(seg_vals[j]),
                         adj = c(0.5, 0.5), srt = 90,
                         cex = seg_cex, col = "white", font = 2)
        }
        y_cursor <- y_cursor + seg_h[j]
      }
    }
    if (!is.null(row$fmt)) {
      draw_bar_axis(y_bot, y_top, 0, profile_max, row$fmt,
                    axis_origin = 0.0, axis_top = 1.0)
    }
  }

  # ── Independent per-row module frames ────────────────────────────
  # Each panel_row gets its OWN background + border so modules are
  # visually independent units stacked in the remaining space — not a
  # shared container. Adding/removing rows simply adds/removes slots.
  module_x1 <- x_block_start - cell_w * 0.06
  module_x2 <- x_block_end   + cell_w * 0.06
  row_frame_pad <- min(module_pad_eff * 0.6, metric_gap_eff * 0.35)
  row_frames <- vector("list", n_module_rows)
  for (i in seq_len(n_module_rows)) {
    rb <- row_ypos[[i]]["bot"] - row_frame_pad
    rt <- row_ypos[[i]]["top"] + row_frame_pad
    row_frames[[i]] <- c(bot = rb, top = rt)
    graphics::rect(module_x1, rb, module_x2, rt,
                   col = "#F3F6F8", border = "#8A95A1", lwd = 0.9)
  }

  # ── Per-column rendering (labels, strips, heatmap, bar slots) ────
  col_cex <- max(0.30, min(0.65, 6 / (n + 2)))
  has_pocp <- !is.null(pocp_mat) && !is.null(pocp_pal) && !is.null(pocp_limits)

  for (g in seq_len(n)) {
    x_l <- x_block_start + (g - 1) * cell_w
    x_r <- x_l + cell_w
    x_c <- (x_l + x_r) / 2

    # Rotated genome label
    graphics::text(x_c, y_lab_bot + 0.003, panel_df$label[g],
                   srt = 90, adj = c(0, 0.5), cex = col_cex,
                   col = "#17212B")

    # ANI / POCP split matrix — each row h uses the h-th genome ring's
    # own radii so the row aligns precisely with that ring at angle 90°
    # apex. h=1 sits at the innermost ring (lowest Y); h=n at the
    # outermost (highest Y). The diagonal therefore runs bottom-left →
    # top-right; ANI stays on the LEFT via the upper-left triangle
    # (h > g) and POCP on the RIGHT via the lower-right triangle
    # (h < g), preserving the side labels.
    if (has_ani) {
      cell_cex <- max(0.22, min(0.38, 4.5 / (n + 4)))
      # Square cells: shrink X by track_margin on each side so the drawn
      # cell width matches its vertical drawn height (which is the ring's
      # r_top - r_bot = cell_w - 2 * track_margin).
      x_cl <- x_l + track_margin
      x_cr <- x_r - track_margin
      for (h in seq_len(n)) {
        y_cb <- ring_radii$r_bot[h]
        y_ct <- ring_radii$r_top[h]
        if (h == g) {
          fill <- diag_col
          v_shown <- 100
          pal_ref <- diag_col
        } else if (h > g) {
          v_shown <- ani_mat[h, g]
          fill <- circos_pangenome_palette_value(v_shown, ani_pal, ani_limits)
          pal_ref <- fill
        } else if (has_pocp) {
          v_shown <- pocp_mat[h, g]
          fill <- circos_pangenome_palette_value(v_shown, pocp_pal, pocp_limits)
          pal_ref <- fill
        } else {
          v_shown <- ani_mat[h, g]
          fill <- circos_pangenome_palette_value(v_shown, ani_pal, ani_limits)
          pal_ref <- fill
        }
        graphics::rect(x_cl, y_cb, x_cr, y_ct, col = fill, border = NA)
        if (is.finite(v_shown) && cell_w > 0.024 && n <= 18) {
          lum <- sum(grDevices::col2rgb(pal_ref) / 255 *
                       c(0.299, 0.587, 0.114))
          txt_col <- if (lum < 0.62) "white" else "#1F2A33"
          label <- if (h == g) "100" else sprintf("%.0f", v_shown)
          graphics::text((x_cl + x_cr) / 2, (y_cb + y_ct) / 2,
                         label, adj = c(0.5, 0.5),
                         cex = cell_cex, col = txt_col)
        }
      }
    }
  }

  # Heatmap outline — trimmed by track_margin to match the square cells.
  if (has_ani) {
    graphics::rect(x_block_start + track_margin, y_ani_bot,
                   x_block_end   - track_margin, y_ani_top,
                   border = "#D7DDE5", lwd = 0.6)
  }

  # ── Panel rows (module bars) ─────────────────────────────────────
  for (i in seq_len(n_module_rows)) {
    row <- panel_rows[[i]]
    y_row_bot <- row_ypos[[i]]["bot"]
    y_row_top <- row_ypos[[i]]["top"]
    switch(row$type,
           stacked = draw_stacked_row(row, y_row_bot, y_row_top),
           bar     = draw_minmax_row(row, y_row_bot, y_row_top),
           stop("Unknown panel_row type: ", row$type))
  }

  # ── ANI / POCP gradient scale (below the heatmap) ────────────────
  if (has_ani) {
    strip_h     <- cell_w * 0.38
    gap_between <- cell_w * 0.08
    if (has_pocp) {
      ani_strip_top <- y_ani_leg_top
      ani_strip_bot <- ani_strip_top - strip_h
      pocp_strip_top <- ani_strip_bot - gap_between
      pocp_strip_bot <- pocp_strip_top - strip_h
    } else {
      ani_strip_top <- y_ani_leg_top
      ani_strip_bot <- ani_strip_top - strip_h
      pocp_strip_top <- pocp_strip_bot <- NA_real_
    }

    draw_gradient <- function(pal, limits, y_bot, y_top, label) {
      n_stops <- length(pal)
      stop_w  <- (x_block_end - x_block_start) / n_stops
      for (k in seq_len(n_stops)) {
        graphics::rect(x_block_start + (k - 1) * stop_w, y_bot,
                       x_block_start + k * stop_w,       y_top,
                       col = pal[k], border = NA)
      }
      graphics::rect(x_block_start, y_bot, x_block_end, y_top,
                     border = "#D7DDE5", lwd = 0.5)
      leg_cex <- max(0.38, min(0.52, 6.5 / (n + 4)))
      mid <- (limits[1] + limits[2]) / 2
      y_mid <- (y_bot + y_top) / 2
      pad_x <- cell_w * 0.12
      draw_inline <- function(x, v, anchor_h) {
        fill <- circos_pangenome_palette_value(v, pal, limits)
        lum  <- sum(grDevices::col2rgb(fill) / 255 *
                      c(0.299, 0.587, 0.114))
        txt_col <- if (lum < 0.62) "white" else "#1F2A33"
        graphics::text(x, y_mid, sprintf("%.0f%%", v),
                       adj = c(anchor_h, 0.5), cex = leg_cex,
                       font = 2, col = txt_col)
      }
      draw_inline(x_block_start + pad_x,                 limits[1], 0)
      draw_inline((x_block_start + x_block_end) / 2,     mid,       0.5)
      draw_inline(x_block_end   - pad_x,                 limits[2], 1)
      graphics::text(x_block_end + cell_w * 0.25, y_mid,
                     label, adj = c(0, 0.5), cex = leg_cex,
                     font = 2, col = "#26303A")
    }

    draw_gradient(ani_pal, ani_limits, ani_strip_bot, ani_strip_top, "ANI")
    if (has_pocp) {
      draw_gradient(pocp_pal, pocp_limits,
                    pocp_strip_bot, pocp_strip_top, "POCP")
      # Rotated "Scale" label spanning the gap between ANI and POCP
      # strips, centered on the gap midpoint. Placed on the left gutter
      # (same x column as other rotated row labels) so it reads inline.
      scale_mid_y <- (ani_strip_bot + pocp_strip_top) / 2
      graphics::text(x_block_start - cell_w * 0.42, scale_mid_y,
                     labels = "Scale", srt = 90, adj = c(0.5, 0.5),
                     cex = lab_cex * 0.9, font = 2, col = "#26303A")
    }
  }

  invisible(list(x_right = x_block_end,
                 y_top   = y_lab_top,
                 y_bot   = y_module_bot))
}

circos_pangenome <- function(dnmb, results_dir = NULL, output_file = NULL) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    warning("circlize not installed")
    return(invisible(NULL))
  }

  # Build the authoritative ordering FIRST so prepare_data + linkage stay
  # consistent. Everything downstream (bin_mat rows, panel rows, ring
  # sequence) follows this one source of truth.
  ordering  <- circos_pangenome_build_ordering(dnmb, results_dir)
  plot_data <- circos_pangenome_prepare_data(dnmb, results_dir = results_dir,
                                             ordering = ordering)
  if (is.null(plot_data)) {
    warning("No clusters or genomes available for circos plotting")
    return(invisible(NULL))
  }

  col_absent <- "#F3F3F3"
  col_core <- "#1F4E5F"
  col_accessory <- "#F0A35B"
  col_unique <- "#D06461"

  species_info <- circos_pangenome_species_palette(plot_data$panel_df)
  col_present_per_genome <- species_info$present_cols

  # Summary tracks come in two flavours:
  #   - outer: drawn OUTSIDE the genome rings (outermost in the plot)
  #   - inner: drawn INSIDE the genome rings but just above the cluster
  #            header — placed radially adjacent to the gene-cluster
  #            track so tightly coupled stats (like num-contributing-
  #            genomes) sit next to the cluster category they describe.
  # Each track draws a bar per cluster with height ∝ value/max.
  summary_defs_outer <- list(
    list(key = "n_genes",
         label = "Num genes in gene cluster",
         fill  = "#A96246"),
    list(key = "max_paralogs",
         label = "Max num paralogs",
         fill  = "#6E5486")
  )
  summary_defs_inner <- list(
    list(key = "n_contrib",
         label = "Num contributing genomes",
         fill  = "#4A7B78")
  )
  summary_defs     <- c(summary_defs_outer, summary_defs_inner)
  summary_n_outer  <- length(summary_defs_outer)
  summary_n_inner  <- length(summary_defs_inner)
  summary_n        <- summary_n_outer + summary_n_inner

  # Linkage: forces panel cell_w ≡ ring radial thickness AND anchors panel
  # x-start at the innermost genome ring's radial edge, so col_center[i]
  # is *exactly* the x-coord where ring i terminates at angle 0°. With
  # outer summary tracks, ring_budget is bumped slightly so the genome
  # rings do not shrink too much.
  cluster_track_h <- 0.022
  ring_budget     <- 0.60
  linkage <- circos_pangenome_build_linkage(
    ordering,
    ring_budget     = ring_budget,
    cluster_track_h = cluster_track_h,
    summary_n       = summary_n,
    summary_n_outer = summary_n_outer,
    summary_n_inner = summary_n_inner,
    summary_frac    = 0.6
  )
  cell_w         <- linkage$cell_w
  genome_track_h <- linkage$genome_track_h
  x_block_start  <- linkage$x_block_start
  x_block_end    <- linkage$x_block_end
  n_g <- plot_data$n_genomes
  has_ani  <- !is.null(plot_data$ani_mat)
  has_pocp <- !is.null(plot_data$pocp_mat)

  # Two sequential palettes ramping from a light neutral up to a SHARED
  # deep anchor (`diag_col`). Matching the 100% endpoint across palettes
  # guarantees any 100% cell (diagonal or otherwise) reads with the same
  # colour in both the ANI (lower-left) and POCP (upper-right) triangles,
  # while the mid-range tones stay distinct (cool blues vs warm burgundy)
  # so the two triangles are still visually separable.
  diag_col <- "#1B2A35"
  ani_pal <- grDevices::colorRampPalette(
    c("#F1F6F9", "#C9DCE5", "#7FAFC2", "#3A6F87", diag_col)
  )(120)
  pocp_pal <- grDevices::colorRampPalette(
    c("#FBF3EA", "#E8C9A5", "#C98A66", "#8C4A3A", diag_col)
  )(120)

  # Fixed scale limits so the legend + cell colours stay comparable
  # across datasets (no data-driven auto-stretch): ANI 70–100 %, POCP
  # 60–100 %.
  ani_limits  <- c(70, 100)
  pocp_limits <- c(60, 100)

  # Anchor the panel just above the ring terminus (the +x axis line, y=0)
  # so the gene-cluster bar is flush with the innermost track. Shift it up
  # by one ring slot to leave room for small colour "flare" blocks that
  # cap each ring's 0° endpoint. `y_top` is derived from the shared panel
  # height helper — no duplication of row-count arithmetic.
  # panel_bottom defines the lower boundary of the panel drawing region.
  # It's used by the proportional-fill algorithm as the floor the module
  # expands down to. Pushing this toward 0 (the circle's horizontal
  # diameter at angle 0°) maximises the module's vertical budget, so
  # the bar rows stretch to fill the whole right-side gap without
  # leaving dead space under the last row.
  panel_bottom <- linkage$cell_w * 0.2
  summary_total_h <- if (summary_n > 0) {
    summary_n *
      (linkage$summary_track_h + 2 * linkage$track_margin)
  } else {
    NULL
  }
  panel_height <- circos_pangenome_panel_total_h(
    cell_w, n_g, has_ani, has_pocp,
    summary_total_h = summary_total_h
  )
  # Extra x-margin reserves room for the right-side scale axis on each
  # bar row (line + tick + 2 numeric labels).
  axis_x_budget <- cell_w * 1.9
  canvas_x_max <- max(1.05, x_block_end + axis_x_budget + 0.04)
  # Reserve a single title row above y = 1 (top of the circle). Per-track
  # labels for the outer summary rings are drawn in-place near each
  # track's 90° apex (inside the first-quadrant gap), so no extra
  # legend-stack budget is needed.
  outer_title_budget <- cell_w * 2.8
  canvas_y_max <- max(
    panel_bottom + panel_height + cell_w * 1.6,
    1.0 + outer_title_budget
  )
  canvas_y_max <- min(canvas_y_max, 5.0)

  # Match the output PDF aspect to the canvas so the circle stays round.
  canvas_w <- canvas_x_max - (-1.08)
  canvas_h <- canvas_y_max - (-1.08)
  aspect   <- canvas_h / canvas_w
  pdf_w    <- 14
  pdf_h    <- max(10, pdf_w * aspect)

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(output_file, width = pdf_w, height = pdf_h)
  } else {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
  }

  on.exit(try(circlize::circos.clear(), silent = TRUE), add = TRUE)
  if (!is.null(output_file)) {
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  }

  graphics::par(mar = c(1.4, 1.1, 3.0, 1.2), xpd = NA)
  circlize::circos.clear()
  circlize::circos.par(
    start.degree = 90,
    clock.wise = FALSE,
    gap.after = 90,
    cell.padding = c(0, 0, 0, 0),
    # Symmetric per-ring margin creates a visible gap between adjacent
    # tracks. Column pitch (cell_w) is set to ring slot (drawn height +
    # 2*track_margin) in build_linkage so col_center ≡ ring r_mid still
    # holds despite the gaps.
    track.margin = c(linkage$track_margin, linkage$track_margin),
    points.overflow.warning = FALSE,
    canvas.xlim = c(-1.08, canvas_x_max),
    canvas.ylim = c(-1.08, canvas_y_max)
  )
  circlize::circos.initialize("pan", xlim = c(0, plot_data$n_clusters))

  # Summary-track renderer (one circos.track per summary def). Used for
  # both outer and inner summary groups — the only difference between
  # them is radial position, which is determined by WHEN this function
  # is called (outer: before genome rings; inner: after genome rings,
  # before cluster header).
  draw_summary_track <- function(sd) {
    vals_raw <- plot_data$cluster_stats[[sd$key]]
    vals <- suppressWarnings(as.numeric(vals_raw))
    vals[!is.finite(vals)] <- 0
    vmax <- max(vals, na.rm = TRUE)
    if (!is.finite(vmax) || vmax <= 0) vmax <- 1
    fill <- sd$fill
    local({
      vals_local <- vals
      vmax_local <- vmax
      fill_local <- fill
      circlize::circos.track(
        "pan",
        ylim = c(0, 1),
        bg.border = "#D7DDE5",
        bg.col    = "#FAFCFD",
        track.height = linkage$summary_track_h,
        panel.fun = function(x, y) {
          circlize::circos.rect(
            xleft   = seq_len(plot_data$n_clusters) - 1,
            ybottom = 0,
            xright  = seq_len(plot_data$n_clusters),
            ytop    = pmin(1, vals_local / vmax_local),
            col     = fill_local,
            border  = NA
          )
        }
      )
    })
  }

  # ── Outer summary tracks (drawn OUTERMOST-FIRST, BEFORE genome rings)
  for (sd in summary_defs_outer) draw_summary_track(sd)

  # Draw genome rings OUTSIDE-IN, pulling per-genome data from the
  # REVERSED panel position — linkage$track_draw_order supplies the
  # mapping so that:
  #   track summary_n+1 (outermost genome) ← panel column n (rightmost)
  #   track summary_n+n (innermost genome) ← panel column 1 (leftmost)
  # Every panel column's center x-coordinate equals its ring's radial
  # midpoint on the +x axis — a single straight vertical line ties the
  # column to its ring's terminus at angle 0°.
  for (i_track in seq_len(plot_data$n_genomes)) {
    local({
      g_panel <- linkage$track_draw_order[i_track]
      present_col_g <- col_present_per_genome[g_panel]
      circlize::circos.track(
        "pan",
        ylim = c(0, 1),
        bg.border = NA,
        track.height = genome_track_h,
        panel.fun = function(x, y) {
          circlize::circos.rect(
            xleft = seq_len(plot_data$n_clusters) - 1,
            ybottom = 0,
            xright = seq_len(plot_data$n_clusters),
            ytop = 1,
            col = ifelse(plot_data$bin_mat[, g_panel] == 1L,
                         present_col_g, col_absent),
            border = NA
          )
        }
      )
    })
  }

  # ── Inner summary tracks — drawn AFTER the genome rings but BEFORE
  # the cluster header, so they sit radially adjacent to the cluster
  # track. Used for stats tightly coupled to cluster identity
  # (e.g. num contributing genomes).
  for (sd in summary_defs_inner) draw_summary_track(sd)

  # ── Cluster-header track at the INNERMOST position — drawn LAST so
  # its track.index becomes summary_n_outer + n + summary_n_inner + 1.
  # Carries core/accessory/unique category colours across clusters.
  cluster_category <- ifelse(
    plot_data$n_per$ng >= plot_data$n_genomes, "Core",
    ifelse(plot_data$n_per$ng == 1L, "Unique", "Accessory")
  )
  cluster_cols <- ifelse(
    cluster_category == "Core", col_core,
    ifelse(cluster_category == "Unique", col_unique, col_accessory)
  )
  # Per-category label text colour: white on dark bgs (Core/Unique),
  # near-black on the orange Accessory bg.
  cat_text_col <- c(Core = "#FFFFFF", Accessory = "#1A1A1A", Unique = "#FFFFFF")
  # Compute the longest contiguous run for each category so the curved
  # label sits on a uniform-colour stretch rather than straddling a
  # category boundary.
  cat_runs <- rle(cluster_category)
  run_ends   <- cumsum(cat_runs$lengths)
  run_starts <- c(1L, utils::head(run_ends, -1L) + 1L)
  circlize::circos.track(
    "pan",
    ylim = c(0, 1),
    bg.border = NA,
    track.height = cluster_track_h,
    panel.fun = function(x, y) {
      circlize::circos.rect(
        xleft = seq_len(plot_data$n_clusters) - 1,
        ybottom = 0,
        xright = seq_len(plot_data$n_clusters),
        ytop = 1,
        col = cluster_cols,
        border = NA
      )
      # Curved in-track category labels (Core / Accessory / Unique)
      # placed on the centre of each category's longest contiguous run.
      # Skip runs that are too short to hold the text legibly.
      min_run_len <- max(plot_data$n_clusters * 0.03, 12)
      for (cat_name in c("Core", "Accessory", "Unique")) {
        idx <- which(cat_runs$values == cat_name)
        if (!length(idx)) next
        pick <- idx[which.max(cat_runs$lengths[idx])]
        if (cat_runs$lengths[pick] < min_run_len) next
        mid_x <- (run_starts[pick] + run_ends[pick] - 1) / 2
        circlize::circos.text(
          x = mid_x, y = 0.5, labels = cat_name,
          facing = "bending.inside", niceFacing = TRUE,
          cex = 0.68, font = 2, col = cat_text_col[[cat_name]]
        )
      }
    }
  )

  # ── Outer-track in-place labels near the 90° apex ────────────────────
  # Rings occupy angles [90°, 360°], so the first quadrant near angle 90°
  # (x just above 0, y near each track's r_mid) is the only clear zone
  # where a horizontal label can sit without overlapping an arc. A tiny
  # swatch + text is drawn at each summary track's 90° apex so the reader
  # can read the track label inline with its arc — no separate legend
  # stack is needed.
  cell_w_num <- linkage$cell_w
  outer_top_r <- circlize::get.cell.meta.data(
    "cell.top.radius", sector.index = "pan", track.index = 1)
  label_cex <- max(0.52, min(0.72, 8 / (plot_data$n_genomes + 6)))
  # Summary swatches share the same X slot and thin-bar style as the
  # 12-o'clock species strip: each swatch's Y range == its track's exact
  # r_bot..r_top (width matches track thickness), X width matches the
  # species strip so the apex gap renders as one stacked thin column.
  sw_x1 <- linkage$track_margin * 0.6
  sw_x2 <- sw_x1 + linkage$cell_w * 0.55
  # Track-index map for summary swatches:
  #   outer summary[s] → track index s         (for s in 1..outer)
  #   inner summary[j] → track index outer + n + j (for j in 1..inner)
  # Both groups use the same swatch style so they read as one legend
  # stack at 12 o'clock, regardless of radial position.
  outer_n <- linkage$summary_n_outer
  inner_n <- linkage$summary_n_inner
  swatch_tracks <- c(seq_len(outer_n),
                     outer_n + linkage$n + seq_len(inner_n))
  for (s in seq_along(swatch_tracks)) {
    ti   <- swatch_tracks[s]
    sd   <- summary_defs[[s]]
    vmax <- max(plot_data$cluster_stats[[sd$key]], na.rm = TRUE)
    if (!is.finite(vmax) || vmax <= 0) vmax <- 1
    r_top_s <- circlize::get.cell.meta.data(
      "cell.top.radius", sector.index = "pan", track.index = ti)
    r_bot_s <- circlize::get.cell.meta.data(
      "cell.bottom.radius", sector.index = "pan", track.index = ti)
    graphics::rect(
      sw_x1, r_bot_s,
      sw_x2, r_top_s,
      col = sd$fill, border = "#3E4954", lwd = 0.3
    )
    graphics::text(
      sw_x2 + linkage$cell_w * 0.25, (r_bot_s + r_top_s) / 2,
      sprintf("%s  (max: %s)", sd$label,
              circos_pangenome_fmt_int(vmax)),
      adj = c(0, 0.5), cex = label_cex, font = 2, col = "#26303A"
    )
  }

  # ── Cluster-header swatch + label at 12 o'clock (innermost track).
  # Uses the same thin-bar X slot as the summary swatches above; three
  # stacked category swatches (Core / Accessory / Unique) so the reader
  # can identify the cluster-header colours inline with the track.
  cluster_ti <- outer_n + linkage$n + inner_n + 1L
  r_cluster_top <- circlize::get.cell.meta.data(
    "cell.top.radius", sector.index = "pan", track.index = cluster_ti)
  r_cluster_bot <- circlize::get.cell.meta.data(
    "cell.bottom.radius", sector.index = "pan", track.index = cluster_ti)
  graphics::rect(
    sw_x1, r_cluster_bot,
    sw_x2, r_cluster_top,
    col = col_accessory, border = "#3E4954", lwd = 0.3
  )
  graphics::text(
    sw_x2 + linkage$cell_w * 0.25, (r_cluster_bot + r_cluster_top) / 2,
    "Gene cluster track",
    adj = c(0, 0.5), cex = label_cex, font = 2, col = "#26303A"
  )

  # Centered title above the circle (no legend stack anymore, so the
  # title sits just above the outermost ring's top radius with a small
  # breathing margin).
  title_y <- outer_top_r + cell_w_num * 1.6
  graphics::text(
    0, title_y,
    sprintf("Pangenome  (%s clusters \u00D7 %s genomes)",
            format(plot_data$n_clusters, big.mark = ","),
            plot_data$n_genomes),
    adj = c(0.5, 0.5), cex = 1.05, font = 2, col = "#16202A"
  )

  # Legend: cluster categories on the left, and a species swatch block so
  # the track colours are readable at a glance.
  graphics::legend(
    "bottomleft",
    legend = c("Core cluster", "Accessory cluster", "Unique cluster", "Absent"),
    fill = c(col_core, col_accessory, col_unique, col_absent),
    border = "grey50",
    cex = 0.78,
    bty = "n",
    title = "Gene clusters",
    title.adj = 0
  )
  graphics::legend(
    "bottomright",
    legend = names(species_info$species_cols),
    fill = unname(species_info$species_cols),
    border = "grey50",
    cex = 0.72,
    bty = "n",
    title = "Species (track colour)",
    title.adj = 0,
    ncol = ceiling(length(species_info$species_cols) / 8)
  )

  # Draw the rectangular metadata panel inside the 90° gap, with rows aligned
  # to each track's actual radial bounds so the circular rings and the
  # rectangular rows line up visually. This uses circlize's own user coords
  # (no par(fig=...) sidecar), so track radii and panel rows share one space.
  panel_extents <- circos_pangenome_draw_aligned_panel(
    panel_df = plot_data$panel_df,
    ani_mat = plot_data$ani_mat,
    pocp_mat = plot_data$pocp_mat,
    hclust_obj = plot_data$hclust_obj,
    order_source = plot_data$order_source,
    col_core = col_core,
    col_accessory = col_accessory,
    col_unique = col_unique,
    ani_pal = ani_pal,
    ani_limits = ani_limits,
    pocp_pal = pocp_pal,
    pocp_limits = pocp_limits,
    diag_col = diag_col,
    species_info = species_info,
    linkage = linkage,
    y_bottom = panel_bottom
  )

  # ── Vertical species strip at 12 o'clock: one narrow rectangle per
  # ring, stacked along the y-axis at the apex (angle 90° terminus).
  # Each slot's Y range matches its ring's [r_bot[i], r_top[i]] so the
  # strip shares the ring's exact thickness AND inter-ring gap — it
  # reads as a direct vertical continuation of the rings at 12 o'clock.
  # Placed entirely in the upper-right gap (x > 0, angle < 90°) so it
  # does not overdraw any sector content.
  radii     <- circos_pangenome_linkage_radii(linkage, sector = "pan")
  strip_x_l <- linkage$track_margin * 0.6
  strip_x_r <- strip_x_l + linkage$cell_w * 0.55
  for (i in seq_len(linkage$n)) {
    fill <- species_info$present_cols[i]
    graphics::rect(
      xleft  = strip_x_l, ybottom = radii$r_bot[i],
      xright = strip_x_r, ytop    = radii$r_top[i],
      col    = fill, border = "white", lwd = 0.3
    )
  }

  if (!is.null(output_file)) {
    message("circos_pangenome written to: ", output_file)
  }

  invisible(NULL)
}
