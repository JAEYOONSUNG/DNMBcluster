#' OrthoVenn3-style area-proportional Euler diagram of cluster sharing
#'
#' Treats each genome as a set of cluster_ids and draws the
#' area-proportional Euler/Venn layout via the `eulerr` package.
#' `eulerr` handles arbitrary set counts (not just the <= 5 ceiling
#' of classical Venn), so runs with 6+ strains still produce a
#' readable diagram. For very large runs (>= 15 genomes) the output
#' is necessarily approximate — it's meant as an overview companion
#' to the exact genome_jaccard_heatmap and presence_absence_heatmap.
#'
#' Each region is labeled with its absolute cluster count. The color
#' palette is deterministic across reruns (seeded viridis).
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional PDF path.
#' @return An eulerr plot object.
#' @export
ortho_euler <- function(dnmb, output_file = NULL) {
  if (!requireNamespace("eulerr", quietly = TRUE)) {
    warning("eulerr not installed — skipping ortho_euler")
    return(invisible(NULL))
  }

  presence <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(
      dnmb$genome_meta %>% dplyr::select(genome_uid, genome_key),
      by = "genome_uid"
    )

  set_list <- split(presence$cluster_id, presence$genome_key)
  n_genomes <- length(set_list)
  if (n_genomes < 2L) {
    warning("ortho_euler: need at least 2 genomes, got ", n_genomes)
    return(invisible(NULL))
  }

  fit <- eulerr::euler(set_list, shape = "ellipse")

  # Discrete, non-adjacent palette for readability.
  palette <- grDevices::hcl.colors(n_genomes, palette = "Zissou 1", alpha = 0.55)

  plot_obj <- plot(
    fit,
    fills    = list(fill = palette),
    edges    = list(col = "grey30", lwd = 1),
    labels   = list(font = 2, cex = 0.85),
    quantities = list(type = "counts", cex = 0.75),
    main     = list(
      label = "Cluster sharing across genomes",
      fontface = "bold"
    )
  )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    # eulerr returns a grid object; use pdf()/print()/dev.off()
    grDevices::pdf(output_file, width = 10, height = 8)
    print(plot_obj)
    grDevices::dev.off()
    message("ortho_euler written to: ", output_file)
  }

  invisible(plot_obj)
}
