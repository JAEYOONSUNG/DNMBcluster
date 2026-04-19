#' Species-pair relationship heatmap (ortholog / paralog / xenolog)
#'
#' Triptych of N x N species-pair heatmaps from `relationships.parquet`.
#' Rows/columns are ordered by the rooted species tree so clades stay
#' contiguous. Color is `log10(count + 1)`; raw counts overlay as text.
#'
#' Anvi'o palette / theme via `.dnmb_anvio_pal` so the figure composes
#' with the rest of the plotting family.
#'
#' @param out_dir `run_orthofinder_like()` output directory.
#' @param out_pdf Destination PDF.
#' @param show_values If TRUE, overlay raw counts on each tile.
#' @param width,height Override auto-scaled PDF dimensions.
#' @param verbose Echo progress.
#' @return Path to the written PDF (invisibly) or `NULL` on failure.
#' @export
plot_relationships_heatmap <- function(out_dir,
                                        out_pdf = NULL,
                                        show_values = TRUE,
                                        width = NULL,
                                        height = NULL,
                                        verbose = TRUE) {
  rel_pq <- file.path(out_dir, "relationships.parquet")
  sp_nwk <- file.path(out_dir, "species_tree_rooted.nwk")
  if (!file.exists(rel_pq)) {
    warning("[relationships_heatmap] missing ", rel_pq, " -- skipping")
    return(NULL)
  }
  for (pkg in c("arrow", "ggplot2", "patchwork")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(sprintf("[relationships_heatmap] requires %s; skipping", pkg))
      return(NULL)
    }
  }

  pal <- .dnmb_anvio_pal

  rel <- as.data.frame(arrow::read_parquet(rel_pq))
  if (!nrow(rel)) {
    warning("[relationships_heatmap] empty relationships -- skipping")
    return(NULL)
  }

  sp_levels <- NULL
  if (file.exists(sp_nwk) && requireNamespace("ape", quietly = TRUE)) {
    tr <- tryCatch(ape::read.tree(sp_nwk), error = function(e) NULL)
    if (!is.null(tr)) {
      tr <- ape::ladderize(tr)
      sp_levels <- tr$tip.label
    }
  }
  if (is.null(sp_levels)) {
    sp_levels <- sort(unique(c(rel$species_a, rel$species_b)))
  }

  rel <- rel[rel$species_a %in% sp_levels & rel$species_b %in% sp_levels, ,
             drop = FALSE]

  build_matrix <- function(df) {
    if (!nrow(df)) return(NULL)
    pairs <- rbind(
      data.frame(a = df$species_a, b = df$species_b, stringsAsFactors = FALSE),
      data.frame(a = df$species_b, b = df$species_a, stringsAsFactors = FALSE)
    )
    tab <- table(factor(pairs$a, levels = sp_levels),
                  factor(pairs$b, levels = sp_levels))
    diag(tab) <- diag(tab) / 2
    long <- as.data.frame(tab, stringsAsFactors = FALSE)
    names(long) <- c("species_a", "species_b", "count")
    long
  }

  orth <- build_matrix(rel[rel$event == "ortholog", , drop = FALSE])
  para <- build_matrix(rel[rel$event == "paralog", , drop = FALSE])
  xeno <- build_matrix(rel[isTRUE_col(rel$is_xenolog_candidate), , drop = FALSE])

  panels <- list(
    orthologs = orth,
    paralogs  = para,
    xenologs  = xeno
  )
  if (!length(panels[!vapply(panels, is.null, logical(1))])) {
    warning("[relationships_heatmap] no relationships to plot -- skipping")
    return(NULL)
  }

  plots <- lapply(names(panels), function(nm) {
    long <- panels[[nm]]
    if (is.null(long)) {
      accent <- .panel_color(nm, pal)
      return(
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0, y = 0.55,
                             label = sprintf("no %s", nm),
                             size = 4.4, colour = pal$subtitle,
                             fontface = "italic") +
          ggplot2::scale_x_continuous(limits = c(-1, 1)) +
          ggplot2::scale_y_continuous(limits = c(0, 1)) +
          ggplot2::labs(title = nm) +
          ggplot2::theme_void(base_size = 11) +
          ggplot2::theme(
            plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
            panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
            plot.title       = ggplot2::element_text(face = "bold", size = 12,
                                                       colour = accent,
                                                       margin = ggplot2::margin(b = 4)),
            plot.margin      = ggplot2::margin(8, 8, 6, 8)
          )
      )
    }
    long$log_count <- log10(long$count + 1)
    long$species_a <- factor(long$species_a, levels = sp_levels)
    long$species_b <- factor(long$species_b, levels = rev(sp_levels))
    accent <- .panel_color(nm, pal)
    g <- ggplot2::ggplot(long,
                         ggplot2::aes(x = .data$species_a,
                                       y = .data$species_b,
                                       fill = .data$log_count)) +
      ggplot2::geom_tile(color = pal$bg, linewidth = 0.35) +
      ggplot2::scale_fill_gradient(
        low = pal$bg, high = accent,
        name = "log10(n+1)",
        guide = ggplot2::guide_colorbar(barwidth = 5, barheight = 0.32,
                                          title.position = "top",
                                          title.hjust = 0.5)
      ) +
      ggplot2::scale_x_discrete(position = "top") +
      ggplot2::coord_equal() +
      ggplot2::labs(title = nm, x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
        panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 0, size = 8,
                                              colour = pal$text),
        axis.text.y = ggplot2::element_text(size = 8, colour = pal$text),
        panel.grid  = ggplot2::element_blank(),
        plot.title  = ggplot2::element_text(face = "bold", size = 12,
                                              colour = accent,
                                              margin = ggplot2::margin(b = 4)),
        legend.position  = "bottom",
        legend.title     = ggplot2::element_text(size = 8.5, face = "bold",
                                                   colour = pal$text),
        legend.text      = ggplot2::element_text(size = 8, colour = pal$text),
        legend.margin    = ggplot2::margin(0, 0, 0, 0),
        plot.margin      = ggplot2::margin(8, 8, 6, 8)
      )
    if (isTRUE(show_values)) {
      g <- g + ggplot2::geom_text(
        ggplot2::aes(label = ifelse(.data$count > 0, .data$count, "")),
        size = 2.5, colour = pal$text
      )
    }
    g
  })

  combined <- patchwork::wrap_plots(plots, nrow = 1, widths = rep(1, length(plots))) +
    patchwork::plot_annotation(
      title = "Species-pair relationships",
      subtitle = sprintf("orthologs=%d  paralogs=%d  xenolog candidates=%d",
                         sum(rel$event == "ortholog"),
                         sum(rel$event == "paralog"),
                         sum(isTRUE_col(rel$is_xenolog_candidate))),
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = pal$bg, colour = NA),
        plot.title    = ggplot2::element_text(face = "bold", size = 14,
                                                colour = pal$text),
        plot.subtitle = ggplot2::element_text(size = 10, colour = pal$subtitle)
      )
    ) &
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal$bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = pal$bg, colour = NA)
    )

  if (isFALSE(out_pdf)) return(invisible(combined))
  if (is.null(out_pdf)) out_pdf <- file.path(out_dir, "relationships_heatmap.pdf")
  n_sp <- length(sp_levels)
  if (is.null(width))  width  <- max(8, length(panels) * (n_sp * 0.4 + 1.6))
  if (is.null(height)) height <- max(5, n_sp * 0.4 + 3.2)
  width  <- min(width, 40)
  height <- min(height, 30)
  ggplot2::ggsave(out_pdf, combined, width = width, height = height,
                   limitsize = FALSE, bg = pal$bg)
  if (verbose) message(sprintf("[relationships_heatmap] wrote %s (%d species, %d panels)",
                                out_pdf, n_sp, length(panels)))
  invisible(out_pdf)
}


.panel_color <- function(name, pal = .dnmb_anvio_pal) {
  switch(name,
    orthologs = pal$navy,
    paralogs  = pal$orange,
    xenologs  = pal$red_strong,
    pal$text
  )
}


isTRUE_col <- function(x) {
  if (is.logical(x)) return(x & !is.na(x))
  v <- tolower(as.character(x))
  !is.na(v) & v %in% c("true", "t", "1")
}
