#' Functional category distribution by pan-genome class
#'
#' Stacked bar showing how CDS products distribute across ~15
#' keyword-based functional categories (see ``functional.py``),
#' faceted by pan-genome category (core / accessory / unique).
#' Core proteins are typically enriched in housekeeping classes
#' (Translation, DNA replication, Energy metabolism) whereas unique
#' proteins are dominated by Hypothetical and Mobile elements.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param results_dir Top-level results directory (needed to locate
#'   ``dnmb/processed/functional_categories.parquet``).
#' @param output_file Optional PDF path.
#' @return A ggplot object, or NULL if the parquet file is absent.
#' @export
functional_bar <- function(dnmb, results_dir, output_file = NULL) {
  func_path <- file.path(
    results_dir, "dnmb", "processed", "functional_categories.parquet"
  )
  if (!file.exists(func_path)) {
    message("functional_bar: ", func_path, " not found — skipping")
    return(invisible(NULL))
  }

  df <- tibble::as_tibble(arrow::read_parquet(func_path))

  # Count per (pan_category, functional_class)
  counts <- df %>%
    dplyr::count(pan_category, functional_class, name = "n") %>%
    dplyr::group_by(pan_category) %>%
    dplyr::mutate(pct = n / sum(n) * 100) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      pan_category = factor(
        pan_category,
        levels = c("core", "accessory", "unique")
      )
    )

  # Order functional classes by total count descending so the biggest
  # category is at the bottom of every bar.
  class_order <- counts %>%
    dplyr::group_by(functional_class) %>%
    dplyr::summarise(total = sum(n), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total)) %>%
    dplyr::pull(functional_class)
  counts$functional_class <- factor(
    counts$functional_class, levels = rev(class_order)
  )

  # Softened journal palette for ~15 classes.
  n_classes <- length(class_order)
  lighten <- function(cols, amount = 0.30) {
    m <- grDevices::col2rgb(cols) / 255
    m <- m + (1 - m) * amount
    grDevices::rgb(m[1, ], m[2, ], m[3, ])
  }
  base_cols <- c(
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#E18727", "#7876B1", "#EFC000",
    "#20854E", "#925E9F", "#CD534C", "#6F99AD", "#EE4C97"
  )
  if (n_classes <= length(base_cols)) {
    pal <- lighten(base_cols[seq_len(n_classes)], 0.30)
  } else {
    pal <- lighten(grDevices::colorRampPalette(base_cols)(n_classes), 0.30)
  }
  names(pal) <- rev(class_order)  # match factor levels

  p <- ggplot2::ggplot(
    counts,
    ggplot2::aes(x = pan_category, y = pct, fill = functional_class)
  ) +
    ggplot2::geom_col(position = "stack", width = 0.7) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::scale_y_continuous(labels = function(v) paste0(v, "%")) +
    ggplot2::labs(
      title    = "Functional category distribution by pan-genome class",
      subtitle = "Product-keyword classification (see functional.py rules)",
      x    = NULL,
      y    = "% of CDS in class",
      fill = NULL
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(ncol = 1, reverse = TRUE)
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position  = "right",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, p, width = 10, height = 7, dpi = 300)
    message("functional_bar written to: ", output_file)
  }

  invisible(p)
}
