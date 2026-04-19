test_that("phylo_tree_plot renders without the ggtree @mapping bug", {
  testthat::skip_if_not_installed("ggtree")
  testthat::skip_if_not_installed("treeio")
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("arrow")
  testthat::skip_if_not_installed("patchwork")

  n_genomes <- 4L
  n_ogs     <- 6L

  synth <- make_synthetic_dnmb(n_genomes = n_genomes, n_ogs = n_ogs, seed = 1L)
  base  <- withr::local_tempdir()
  write_synthetic_dnmb_parquets(synth, base)

  # Augment genome_meta with the columns phylo_tree_plot reads.
  gm <- synth$genome_meta
  gm$organism    <- paste("Synth genome", gm$genome_uid)
  gm$strain      <- paste0("strain-", gm$genome_uid)
  gm$gc_percent  <- runif(n_genomes, 38, 48)
  gm$total_length<- runif(n_genomes, 4e6, 5.5e6)
  gm$n_cds       <- sample(3800:4300, n_genomes)
  gm$core        <- rep(n_ogs, n_genomes)
  gm$accessory   <- sample(0:3, n_genomes, replace = TRUE)
  gm$unique      <- sample(0:2, n_genomes, replace = TRUE)
  arrow::write_parquet(gm, file.path(base, "dnmb", "inputs", "genome_meta.parquet"))

  # cluster_summary needs a `category` column for the side panel lookup.
  cs <- synth$clusters %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(n_genomes = dplyr::n_distinct(genome_uid),
                     n_members = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(category = dplyr::case_when(
      n_genomes == !!n_genomes ~ "core",
      n_genomes == 1L          ~ "unique",
      TRUE                     ~ "accessory"
    ))
  arrow::write_parquet(cs, file.path(base, "dnmb", "processed",
                                     "cluster_summary.parquet"))

  # Minimal valid Newick with internal node.label (support).
  tree <- ape::rcoal(n_genomes, tip.label = gm$genome_key)
  tree$node.label <- as.character(sample(60:100, tree$Nnode, replace = TRUE))
  ape::write.tree(tree, file.path(base, "dnmb", "processed", "phylo_tree.nwk"))

  dnmb <- load_dnmb(base)
  out_pdf <- file.path(base, "phylo_tree.pdf")

  # We expect no thrown error AND no fallback warning. The @mapping
  # corruption manifests at ggsave time as "<ggtree> object properties
  # are invalid: @mapping must be <ggplot2::mapping>, not S3<data.frame>".
  # The stat_tree aesthetics bug shows as "Aesthetics must be either
  # length 1 or the same as the data (N). Fix: from, to, .panel".
  #
  # Both failure modes are currently caught by an internal tryCatch in
  # phylo_tree_plot() which warns "ggsave of combined plot failed ... --
  # saving tree-only" and emits the tree alone. We assert that neither
  # an error nor that warning is raised -- suppressWarnings() would hide
  # exactly the signal we care about.
  warnings_seen <- character()
  expect_no_error(
    withCallingHandlers(
      phylo_tree_plot(dnmb, base, output_file = out_pdf),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  )
  fallback_hits <- grep("ggsave of combined plot failed|saving tree-only",
                        warnings_seen, value = TRUE)
  expect_identical(fallback_hits, character(0))

  expect_true(file.exists(out_pdf))
  # Sanity floor only; the fallback_hits check above is what actually
  # distinguishes combined-plot success from tree-only fallback.
  expect_gt(file.info(out_pdf)$size, 1000)
})
