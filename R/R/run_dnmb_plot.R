#' Run the full DNMBcluster plotting suite on a results directory
#'
#' Entry point invoked by the Python CLI. Loads the DNMB Parquet
#' outputs, then generates every plot defined in the package into
#' `<results_dir>/plots/`.
#'
#' @param results_dir Path produced by `dnmbcluster run`.
#' @param output_dir Optional override; defaults to `<results_dir>/plots`.
#' @return A named list of ggplot objects.
#' @export
run_dnmb_plot <- function(results_dir, output_dir = NULL) {
  if (is.null(output_dir)) {
    output_dir <- file.path(results_dir, "plots")
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  dnmb <- load_dnmb(results_dir)

  plots <- list()

  plots$flower <- tryCatch(
    flower_plot(dnmb, output_file = file.path(output_dir, "flower_plot.pdf")),
    error = function(e) {
      warning("flower_plot failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$pan_core <- tryCatch(
    pan_core_plot(dnmb, output_file = file.path(output_dir, "pan_core_plot.pdf")),
    error = function(e) {
      warning("pan_core_plot failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$category_bar <- tryCatch(
    category_bar_plot(dnmb, output_file = file.path(output_dir, "category_bar.pdf")),
    error = function(e) {
      warning("category_bar_plot failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$presence_absence <- tryCatch(
    presence_absence_heatmap(
      dnmb,
      output_file = file.path(output_dir, "presence_absence_heatmap.pdf")
    ),
    error = function(e) {
      warning("presence_absence_heatmap failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$cluster_size <- tryCatch(
    cluster_size_dist(
      dnmb,
      output_file = file.path(output_dir, "cluster_size_dist.pdf")
    ),
    error = function(e) {
      warning("cluster_size_dist failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$identity_dist <- tryCatch(
    identity_distribution(
      dnmb,
      output_file = file.path(output_dir, "identity_distribution.pdf")
    ),
    error = function(e) {
      warning("identity_distribution failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$genome_jaccard <- tryCatch(
    genome_jaccard_heatmap(
      dnmb,
      output_file = file.path(output_dir, "genome_jaccard_heatmap.pdf")
    ),
    error = function(e) {
      warning("genome_jaccard_heatmap failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$coverage_scatter <- tryCatch(
    coverage_scatter(
      dnmb,
      output_file = file.path(output_dir, "coverage_scatter.pdf")
    ),
    error = function(e) {
      warning("coverage_scatter failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$gene_content_mds <- tryCatch(
    gene_content_mds(
      dnmb,
      output_file = file.path(output_dir, "gene_content_mds.pdf")
    ),
    error = function(e) {
      warning("gene_content_mds failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$centroid_length <- tryCatch(
    centroid_length_distribution(
      dnmb,
      output_file = file.path(output_dir, "centroid_length_distribution.pdf")
    ),
    error = function(e) {
      warning("centroid_length_distribution failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$core_gene_conservation <- tryCatch(
    core_gene_conservation(
      dnmb,
      output_file = file.path(output_dir, "core_gene_conservation.pdf")
    ),
    error = function(e) {
      warning("core_gene_conservation failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$cumulative_pan_contribution <- tryCatch(
    cumulative_pan_contribution(
      dnmb,
      output_file = file.path(output_dir, "cumulative_pan_contribution.pdf")
    ),
    error = function(e) {
      warning("cumulative_pan_contribution failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$singleton_top_n <- tryCatch(
    singleton_top_n(
      dnmb,
      output_file = file.path(output_dir, "singleton_top_n.pdf")
    ),
    error = function(e) {
      warning("singleton_top_n failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$ortho_euler <- tryCatch(
    ortho_euler(
      dnmb,
      output_file = file.path(output_dir, "ortho_euler.pdf")
    ),
    error = function(e) {
      warning("ortho_euler failed: ", conditionMessage(e))
      NULL
    }
  )

  plots$functional_bar <- tryCatch(
    functional_bar(
      dnmb,
      results_dir      = results_dir,
      output_file = file.path(output_dir, "functional_bar.pdf")
    ),
    error = function(e) {
      warning("functional_bar failed: ", conditionMessage(e))
      NULL
    }
  )

  # ANI + POCP heatmaps (emit only if the parquet files exist).
  plots$ani_pocp <- tryCatch(
    ani_pocp_heatmap(
      dnmb,
      results_dir      = results_dir,
      output_file_ani  = file.path(output_dir, "ani_heatmap.pdf"),
      output_file_pocp = file.path(output_dir, "pocp_heatmap.pdf")
    ),
    error = function(e) {
      warning("ani_pocp_heatmap failed: ", conditionMessage(e))
      NULL
    }
  )

  # Euler + UpSet combined figure.
  plots$euler_upset <- tryCatch(
    euler_upset_combined(
      dnmb,
      output_file = file.path(output_dir, "euler_upset.pdf")
    ),
    error = function(e) {
      warning("euler_upset_combined failed: ", conditionMessage(e))
      NULL
    }
  )

  # Circos presence/absence (anvi'o-style).
  plots$circos <- tryCatch(
    circos_pangenome(
      dnmb,
      output_file = file.path(output_dir, "circos_pangenome.pdf")
    ),
    error = function(e) {
      warning("circos_pangenome failed: ", conditionMessage(e))
      NULL
    }
  )

  # COG category bar (emits only if --annotate ran eggnog-fast).
  plots$cog_bar <- tryCatch(
    cog_bar(
      dnmb,
      results_dir      = results_dir,
      output_file = file.path(output_dir, "cog_bar.pdf")
    ),
    error = function(e) {
      warning("cog_bar failed: ", conditionMessage(e))
      NULL
    }
  )

  # Circular anvi'o-style display (emits only if --phylo ran).
  plots$phylo_circular <- tryCatch(
    phylo_circular(
      dnmb,
      results_dir = results_dir,
      output_file = file.path(output_dir, "phylo_circular.pdf")
    ),
    error = function(e) {
      warning("phylo_circular failed: ", conditionMessage(e))
      NULL
    }
  )

  # Core-gene IQ-TREE visualization (emits only if --phylo ran).
  plots$phylo_tree <- tryCatch(
    phylo_tree_plot(
      dnmb,
      results_dir = results_dir,
      output_file = file.path(output_dir, "phylo_tree.pdf")
    ),
    error = function(e) {
      warning("phylo_tree_plot failed: ", conditionMessage(e))
      NULL
    }
  )

  message("DNMBcluster plots saved to: ", output_dir)
  invisible(plots)
}
