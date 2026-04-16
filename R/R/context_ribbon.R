#' Locus neighborhood context plot with cross-strain homology ribbons
#'
#' Given a user-chosen ``(genome_key, locus_tag)`` anchor, extracts a
#' fixed-bp window around the anchor in every genome that shares the
#' anchor's cluster and draws a stacked gggenes panel with:
#'
#' - **Cross-strain homology ribbons** — a translucent polygon
#'   connects every pair of adjacent strain rows wherever they share
#'   the same DNMBcluster ``cluster_id``. Ribbons are filled in the
#'   cluster's color so the eye follows "this ortholog is here in
#'   every strain" as a continuous gradient stripe across the
#'   stack. (gggenomes-style homology flow without the gggenomes
#'   dependency, which isn't on conda.)
#' - **Anchor-centered + strand-normalized** — coordinates are
#'   shifted so the anchor's midpoint sits at ``x = 0`` on every row;
#'   the ★ marker forms a vertical column down the page. Rows where
#'   the anchor lies on the minus strand get the whole window
#'   reflected so every row reads left-to-right.
#' - **locus_tag labels** printed inside every arrow via
#'   ``gggenes::geom_gene_label`` (clipped automatically on narrow
#'   arrows) plus an italic below-arrow print for the anchor row.
#' - **TSV export** — a ``*.tsv`` sibling of the PDF lists every
#'   CDS in the drawn window with both original bp coordinates and
#'   the shifted/normalized ones, the cluster_id, product annotation,
#'   and an ``is_anchor`` flag. Drop into Excel / pandas / polars
#'   when the visual isn't enough.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param anchor_genome Genome key of the starting locus.
#' @param anchor_locus Locus_tag (or cds_key fallback) of the anchor CDS.
#' @param window_bp Half-width of the visible window in base pairs.
#'   Default 25 000 → 50 kb total view.
#' @param output_file Optional PDF path. A sibling ``.tsv`` is
#'   written next to it with the per-CDS table.
#' @return A ggplot object (invisibly).
#' @export
context_ribbon <- function(dnmb, anchor_genome, anchor_locus,
                           window_bp = 25000L, output_file = NULL) {
  if (!requireNamespace("gggenes", quietly = TRUE)) {
    stop("gggenes not installed — cannot render context ribbon")
  }

  id_map   <- dnmb$id_map
  clusters <- dnmb$clusters
  meta     <- dnmb$genome_meta

  # Resolve anchor → protein_uid + contig + coord + cluster_id.
  anchor_row <- id_map %>%
    dplyr::filter(genome_key == anchor_genome) %>%
    dplyr::filter(locus_tag == anchor_locus | cds_key == anchor_locus)
  if (nrow(anchor_row) == 0L) {
    stop(sprintf(
      "anchor %s::%s not found in id_map", anchor_genome, anchor_locus
    ))
  }
  anchor_row <- anchor_row[1, ]
  anchor_uid <- anchor_row$protein_uid
  anchor_cluster <- clusters$cluster_id[clusters$protein_uid == anchor_uid][1]
  if (is.na(anchor_cluster)) {
    stop("anchor locus has no cluster assignment")
  }

  sharing_uids <- clusters$genome_uid[clusters$cluster_id == anchor_cluster]
  sharing_keys <- unique(meta$genome_key[meta$genome_uid %in% sharing_uids])

  # ---------------------------------------------------------------
  # Per-strain window extraction (anchor-centered, strand-normalized)
  # ---------------------------------------------------------------
  build_panel <- function(gkey) {
    row_ids <- id_map %>% dplyr::filter(genome_key == gkey)
    cluster_members_here <- clusters$protein_uid[
      clusters$cluster_id == anchor_cluster &
      clusters$genome_uid %in% meta$genome_uid[meta$genome_key == gkey]
    ]
    if (length(cluster_members_here) == 0L) return(NULL)
    anchor_here <- row_ids %>%
      dplyr::filter(protein_uid %in% cluster_members_here)
    if (nrow(anchor_here) == 0L) return(NULL)
    anchor_here <- anchor_here[1, ]

    contig_here <- anchor_here$contig
    anchor_mid <- (anchor_here$start + anchor_here$end) / 2
    anchor_strand <- anchor_here$strand
    if (is.na(anchor_strand)) anchor_strand <- 1L

    w <- row_ids %>%
      dplyr::filter(contig == !!contig_here) %>%
      dplyr::filter(end >= anchor_mid - window_bp,
                    start <= anchor_mid + window_bp) %>%
      dplyr::arrange(start) %>%
      dplyr::left_join(
        clusters %>% dplyr::select(protein_uid, pct_identity_fwd),
        by = "protein_uid"
      ) %>%
      dplyr::mutate(
        s_shift = start - anchor_mid,
        e_shift = end   - anchor_mid
      )

    if (anchor_strand < 0) {
      w <- w %>% dplyr::mutate(
        new_start  = -e_shift,
        new_end    = -s_shift,
        new_strand = -strand
      )
      # After reflection the visual left edge maps to a LARGER bp
      # coord than the visual right edge. Keep both so the label
      # layer can print them honestly at the window extremes.
      orig_left  <- anchor_mid + window_bp
      orig_right <- anchor_mid - window_bp
    } else {
      w <- w %>% dplyr::mutate(
        new_start  = s_shift,
        new_end    = e_shift,
        new_strand = strand
      )
      orig_left  <- anchor_mid - window_bp
      orig_right <- anchor_mid + window_bp
    }

    w$genome_key <- gkey
    w$is_anchor  <- w$protein_uid == anchor_here$protein_uid
    w$cluster_id <- clusters$cluster_id[
      match(w$protein_uid, clusters$protein_uid)
    ]
    w$orig_left_bp  <- orig_left
    w$orig_right_bp <- orig_right
    w$contig_name   <- contig_here
    w
  }

  panels <- lapply(sharing_keys, build_panel)
  panels <- panels[!vapply(panels, is.null, logical(1))]
  if (length(panels) == 0L) {
    stop("no genomes share the anchor cluster — nothing to plot")
  }
  df <- dplyr::bind_rows(panels)

  # Fill class per gene: shared clusters get a color, unique genes fade to grey.
  cluster_counts <- table(df$cluster_id)
  shared_clusters <- names(cluster_counts[cluster_counts >= 2L])
  df$fill_class <- ifelse(
    as.character(df$cluster_id) %in% shared_clusters,
    paste0("c", df$cluster_id),
    "unique"
  )
  shared_ids <- sort(unique(df$fill_class[df$fill_class != "unique"]))
  n_shared <- length(shared_ids)

  # Harmonized 20-color palette — curated from NPG, Lancet, NEJM, JCO
  # with the 5 darkest tones removed and the remaining colors
  # reordered warm↔cool so neighboring clusters always contrast in
  # temperature. 45% white blend lifts everything to a uniform
  # pastel register that's easy on the eye even at full opacity.
  journal_raw <- c(
    "#E64B35",  # 1  coral        (warm)
    "#4DBBD5",  # 2  sky blue     (cool)
    "#42B540",  # 3  green        (cool)
    "#F39B7F",  # 4  salmon       (warm)
    "#7876B1",  # 5  lavender     (cool)
    "#E18727",  # 6  orange       (warm)
    "#91D1C2",  # 7  mint         (cool)
    "#EE4C97",  # 8  pink         (warm)
    "#00A087",  # 9  teal         (cool)
    "#EFC000",  # 10 gold         (warm)
    "#6F99AD",  # 11 grey-blue    (cool)
    "#CD534C",  # 12 brick        (warm)
    "#20854E",  # 13 forest       (cool)
    "#FFDC91",  # 14 pale gold    (warm)
    "#8491B4",  # 15 steel        (cool)
    "#925E9F",  # 16 purple       (cool)
    "#7AA6DC",  # 17 light blue   (cool)
    "#B09C85",  # 18 tan          (warm)
    "#0073C2",  # 19 blue         (cool)
    "#3C5488"   # 20 navy         (cool)
  )
  lighten_cols <- function(cols, amount = 0.45) {
    m <- grDevices::col2rgb(cols) / 255
    m <- m + (1 - m) * amount
    grDevices::rgb(m[1, ], m[2, ], m[3, ])
  }
  journal_soft <- lighten_cols(journal_raw, 0.38)

  n_pal <- max(1L, n_shared)
  if (n_pal <= length(journal_soft)) {
    pal_colors <- journal_soft[seq_len(n_pal)]
  } else {
    pal_colors <- grDevices::colorRampPalette(journal_soft)(n_pal)
  }
  palette <- c(unique = "#ECEAE7", setNames(pal_colors, shared_ids))

  # Row order: anchor genome at the top, then the remaining strains
  # sorted by **descending identity to the anchor** (for the anchor
  # cluster specifically). This produces a visual gradient — the most
  # similar strain sits right below the anchor and the most divergent
  # sinks to the bottom, which matches the intuition "fade = further
  # from reference".
  #
  # For each strain in the anchor cluster, compute its identity to
  # the anchor gene (centroid-approx rules apply). Strains NOT in the
  # anchor cluster get identity = 0 so they sink to the bottom.
  anchor_genome_uid_val <- meta$genome_uid[meta$genome_key == anchor_genome][1]
  anchor_pfwd_val <- clusters$pct_identity_fwd[clusters$protein_uid == anchor_uid][1]
  anchor_is_cen <- is.na(anchor_pfwd_val)

  strain_sim <- vapply(unique(df$genome_key), function(gkey) {
    if (gkey == anchor_genome) return(100)
    gkey_uid <- meta$genome_uid[meta$genome_key == gkey][1]
    pfwds <- clusters$pct_identity_fwd[
      clusters$cluster_id == anchor_cluster & clusters$genome_uid == gkey_uid
    ]
    if (length(pfwds) == 0L) return(0)
    pfwd <- pfwds[1]
    if (is.na(pfwd) && anchor_is_cen) return(100)
    if (is.na(pfwd)) return(if (is.na(anchor_pfwd_val)) 100 else anchor_pfwd_val)
    if (anchor_is_cen) return(pfwd)
    return(min(anchor_pfwd_val, pfwd, na.rm = TRUE))
  }, numeric(1))

  sim_order <- names(sort(strain_sim, decreasing = TRUE))
  row_levels <- sim_order
  # Bottom-up factor → row 1 sits at the bottom, highest level at top.
  df$genome_key <- factor(df$genome_key, levels = rev(row_levels))
  df$y_num <- as.numeric(df$genome_key)

  # Per-strain summary: one row per genome with the 4-line display
  # label (Genus / species / strain / accession), window bp bounds,
  # and the numeric y coordinate. Drives the backbone segment geom,
  # the bp-at-edges text, and the y-axis tick labels.
  parse_parts <- function(meta_row) {
    org <- meta_row$organism
    st  <- meta_row$strain
    parts <- list(genus = "", species = "", strain = "")
    if (!is.na(org) && nzchar(org)) {
      tokens <- strsplit(trimws(org), "\\s+")[[1]]
      if (length(tokens) >= 1L) parts$genus   <- tokens[1]
      if (length(tokens) >= 2L) parts$species <- tokens[2]
      if (is.na(st) || !nzchar(st)) {
        if (length(tokens) >= 3L) {
          parts$strain <- paste(tokens[3:length(tokens)], collapse = " ")
        }
      } else {
        parts$strain <- st
      }
    } else if (!is.na(st) && nzchar(st)) {
      parts$strain <- st
    }
    # Strip "strain " prefix routinely stuffed into RefSeq organism tails.
    parts$strain <- sub("^strain\\s+", "", parts$strain, ignore.case = TRUE)
    parts
  }

  per_strain <- df %>%
    dplyr::distinct(genome_key, y_num, orig_left_bp, orig_right_bp, contig_name) %>%
    dplyr::arrange(y_num)

  four_line <- character(nrow(per_strain))
  for (i in seq_len(nrow(per_strain))) {
    gkey <- as.character(per_strain$genome_key[i])
    mr <- meta[meta$genome_key == gkey, ][1, ]
    parts <- parse_parts(mr)
    lines <- c(parts$genus, parts$species, parts$strain, gkey)
    lines <- lines[nzchar(lines)]
    four_line[i] <- paste(lines, collapse = "\n")
  }
  per_strain$row_label <- four_line

  y_breaks <- per_strain$y_num
  y_labels <- per_strain$row_label

  # ---------------------------------------------------------------
  # Cross-strain homology ribbons (polygons between adjacent rows)
  # ---------------------------------------------------------------
  # Vertical row pitch (inches per strain row) — chosen so that
  # genome rows sit about 2/3 as tall as the previous 1.1 in/row
  # setting without making 4-line y-axis labels collide.
  ROW_HEIGHT_IN <- 0.73
  # Ribbon offset in y-units. gggenes draws arrow bodies at a fixed
  # 3.5 mm tall regardless of data scaling, so the half-height in
  # y-units depends on how tall a row is in inches. We inflate the
  # half-height by a tiny buffer (0.015) so the ribbon edge shows
  # through as a thin line hugging the arrow outline.
  ribbon_offset_y <- 1.75 / (ROW_HEIGHT_IN * 25.4) + 0.015

  # Pre-compute the anchor strain's pct_identity_fwd for every cluster
  # shown in the plot so build_ribbons can compute "identity relative
  # to anchor" instead of a generic min(top, bot).
  anchor_genome_uid <- meta$genome_uid[meta$genome_key == anchor_genome][1]
  cluster_anchor_pfwd <- list()
  for (cid in unique(df$cluster_id[!is.na(df$cluster_id)])) {
    pfwd <- clusters$pct_identity_fwd[
      clusters$cluster_id == cid & clusters$genome_uid == anchor_genome_uid
    ]
    # NA when anchor IS the centroid of this cluster; NULL when
    # anchor has no gene in this cluster at all.
    cluster_anchor_pfwd[[as.character(cid)]] <- if (length(pfwd) == 0L) NULL else pfwd[1]
  }

  build_ribbons <- function(df) {
    lvls <- levels(df$genome_key)
    n <- length(lvls)
    if (n < 2L) return(NULL)
    out <- list()
    rib_id <- 0L
    for (i in seq_len(n - 1L)) {
      top_lvl <- lvls[i + 1L]  # visually higher
      bot_lvl <- lvls[i]
      top_df <- df[df$genome_key == top_lvl, ]
      bot_df <- df[df$genome_key == bot_lvl, ]
      shared <- intersect(top_df$cluster_id, bot_df$cluster_id)
      shared <- shared[!is.na(shared)]
      shared <- shared[as.character(shared) %in% shared_clusters]
      for (cid in shared) {
        top_gene <- top_df[top_df$cluster_id == cid, ][1, ]
        bot_gene <- bot_df[bot_df$cluster_id == cid, ][1, ]

        # Anchor-referenced alpha: for each ribbon, compute how
        # similar the BOTTOM row's gene is to the anchor strain's
        # gene in the same cluster. This makes the "fade" relative to
        # a single consistent reference so a cluster that's ~98% in
        # strain B reads differently from one that's ~85% in strain D.
        #
        # When the anchor is the cluster centroid (its pct_identity_fwd
        # is NA → substitute 100), the bottom strain's pct_identity_fwd
        # IS the exact anchor-vs-bottom identity. Otherwise approximate
        # via min(anchor's pfwd, bottom's pfwd).
        bot_pid <- bot_gene$pct_identity_fwd
        if (is.na(bot_pid)) bot_pid <- 100
        ap <- cluster_anchor_pfwd[[as.character(cid)]]
        if (is.null(ap)) {
          # Anchor not a member of this cluster → pairwise fallback
          top_pid <- top_gene$pct_identity_fwd
          if (is.na(top_pid)) top_pid <- 100
          identity_pair <- min(top_pid, bot_pid)
        } else if (is.na(ap)) {
          # Anchor IS the centroid → bot's pfwd is exact
          identity_pair <- bot_pid
        } else {
          identity_pair <- min(ap, bot_pid)
        }
        # Power-law mapping [75, 100] → [0.08, 0.95] with a ^2 curve
        # so subtle identity differences produce visible alpha gaps.
        norm <- max(0, min(1, (identity_pair - 75) / 25))
        alpha_val <- 0.08 + norm^2 * 0.87

        rib_id <- rib_id + 1L
        out[[length(out) + 1L]] <- data.frame(
          x    = c(top_gene$new_start, top_gene$new_end,
                   bot_gene$new_end,   bot_gene$new_start),
          # Ribbon vertices hug the arrow bodies tightly. Offset is
          # derived from the figure's row-pitch so changes to
          # ROW_HEIGHT_IN automatically keep ribbons flush with the
          # gggenes arrow outline.
          y    = c(i + 1L - ribbon_offset_y, i + 1L - ribbon_offset_y,
                   i +     ribbon_offset_y, i +     ribbon_offset_y),
          group         = rib_id,
          fill_class    = top_gene$fill_class,
          identity_pair = identity_pair,
          alpha_val     = alpha_val,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(out) == 0L) return(NULL)
    dplyr::bind_rows(out)
  }
  ribbons_df <- build_ribbons(df)

  # Label text: prefer gene name (e.g. "dnaA") since it's shorter
  # and biologically informative; fall back to locus_tag, then cds_key.
  df$gene_label <- dplyr::case_when(
    !is.na(df$gene) & df$gene != "" ~ df$gene,
    !is.na(df$locus_tag) & df$locus_tag != "" ~ df$locus_tag,
    TRUE ~ df$cds_key
  )

  anchor_df <- dplyr::filter(df, is_anchor)

  # ---------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------
  p <- ggplot2::ggplot()

  if (!is.null(ribbons_df) && nrow(ribbons_df) > 0L) {
    p <- p + ggplot2::geom_polygon(
      data = ribbons_df,
      mapping = ggplot2::aes(
        x = x, y = y, group = group,
        fill = fill_class, alpha = alpha_val
      ),
      color = NA
    ) +
      ggplot2::scale_alpha_identity()
  }

  p <- p +
    ggplot2::geom_vline(
      xintercept = 0, color = "#A0A0A0",
      linetype = "dotted", linewidth = 0.4
    ) +
    # Per-strain backbone: thin black line across the full window
    # behind the ribbons and arrows. Acts as a visual skeleton so
    # gaps in the arrow row read as "nothing there" instead of
    # "undefined space".
    ggplot2::geom_segment(
      data = per_strain,
      mapping = ggplot2::aes(
        x = -window_bp, xend = window_bp,
        y = y_num, yend = y_num
      ),
      color = "#202020", linewidth = 0.35,
      inherit.aes = FALSE
    ) +
    gggenes::geom_gene_arrow(
      data = df,
      mapping = ggplot2::aes(
        xmin = new_start, xmax = new_end,
        y    = y_num,
        fill = fill_class,
        forward = new_strand > 0
      ),
      arrowhead_height  = grid::unit(4, "mm"),
      arrowhead_width   = grid::unit(2.5, "mm"),
      arrow_body_height = grid::unit(3.5, "mm"),
      color             = "#303030",
      size              = 0.25
    ) +
    # Uniform-size gene labels at gene midpoint. check_overlap drops
    # labels that would collide with already-drawn text so narrow
    # gene clusters stay readable instead of becoming an ink blob.
    ggplot2::geom_text(
      data = df,
      mapping = ggplot2::aes(
        x = (new_start + new_end) / 2,
        y = y_num,
        label = gene_label
      ),
      size = 1.8, color = "#303030",
      inherit.aes = FALSE,
      check_overlap = TRUE
    ) +
    ggplot2::geom_text(
      data = anchor_df,
      ggplot2::aes(x = (new_start + new_end) / 2, y = y_num,
                   label = "\u2605"),
      color = "#202020", size = 5.2, inherit.aes = FALSE,
      nudge_y = 0.30
    ) +
    # Per-strain original bp positions at the window extremes.
    # Pinned just inside the panel so they stay visible when
    # coord_cartesian clips the axis. For minus-strand anchors the
    # left label will be the LARGER bp value — that's the honest
    # mapping after window reflection. nudge_y pushes the text well
    # above the backbone line so it never collides with the arrow
    # bodies underneath.
    # Per-strain bp positions at the window extremes, pushed BELOW
    # the backbone line so they never overlap the gene arrows above.
    ggplot2::geom_text(
      data = per_strain,
      mapping = ggplot2::aes(
        x = -window_bp * 0.985, y = y_num,
        label = format(round(orig_left_bp), big.mark = ",")
      ),
      hjust = 0, nudge_y = -0.15,
      color = "#606060", size = 2.5, inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = per_strain,
      mapping = ggplot2::aes(
        x = window_bp * 0.985, y = y_num,
        label = format(round(orig_right_bp), big.mark = ",")
      ),
      hjust = 1, nudge_y = -0.15,
      color = "#606060", size = 2.5, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::scale_x_continuous(
      labels = function(v) paste0(v / 1000, " kb"),
      breaks = scales::pretty_breaks(n = 9)
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      # The bottom row needs deeper headroom so its arrow bodies
      # aren't clipped by the panel edge after the bp-position
      # labels push text up above the backbone. Top needs room for
      # the nudged star + bp labels on the top row.
      expand = ggplot2::expansion(add = c(1.3, 0.8))
    ) +
    ggplot2::coord_cartesian(xlim = c(-window_bp, window_bp), expand = FALSE) +
    gggenes::theme_genes() +
    ggplot2::theme(
      text        = ggplot2::element_text(size = 12),
      # Four-line tick labels (Genus / species / strain / accession)
      # need smaller font + tight lineheight so they fit inside one
      # row slot without overflowing into neighbors.
      axis.text.y = ggplot2::element_text(size = 7.5, lineheight = 0.85),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = sprintf(
        "Locus neighborhood of %s::%s  (cluster %d)",
        anchor_genome, anchor_locus, anchor_cluster
      ),
      subtitle = sprintf(
        paste0(
          "\u00b1%d kb window, anchor-centered, L\u2192R oriented  |  ",
          "%d strains  |  \u2605 = anchor homolog  |  ",
          "ribbon opacity \u221d pairwise pct_identity_fwd"
        ),
        as.integer(window_bp / 1000), length(row_levels)
      ),
      x = "distance from anchor midpoint",
      y = NULL
    )

  # ---------------------------------------------------------------
  # Write the per-CDS table as a sibling TSV
  # ---------------------------------------------------------------
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    # Row pitch (see ROW_HEIGHT_IN above) × strain count + 2 in
    # for the title and top/bottom margins.
    fig_height <- ROW_HEIGHT_IN * length(row_levels) + 2
    ggplot2::ggsave(
      output_file, p, width = 14, height = fig_height, dpi = 300
    )
    message("context_ribbon written to: ", output_file)

    tsv_path <- sub("\\.pdf$", ".tsv", output_file)
    if (tsv_path == output_file) tsv_path <- paste0(output_file, ".tsv")

    # Add product annotation when available via id_map (already joined).
    product_col <- if ("product" %in% names(df)) df$product else NA_character_

    out_df <- data.frame(
      genome_key     = as.character(df$genome_key),
      is_anchor      = df$is_anchor,
      locus_tag      = df$locus_tag,
      cds_key        = df$cds_key,
      cluster_id     = df$cluster_id,
      contig         = df$contig,
      orig_start     = df$start,
      orig_end       = df$end,
      orig_strand    = df$strand,
      window_start   = df$new_start,
      window_end     = df$new_end,
      window_strand  = df$new_strand,
      product        = product_col,
      stringsAsFactors = FALSE
    )
    out_df <- out_df[order(out_df$genome_key, out_df$window_start), ]

    utils::write.table(
      out_df, file = tsv_path,
      sep = "\t", quote = FALSE, row.names = FALSE, na = ""
    )
    message("context_ribbon table: ", tsv_path)
  }

  invisible(p)
}
