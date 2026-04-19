test_that("run_phylogenomics_sota runs end-to-end on synthetic data", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("igraph")

  dnmb <- make_synthetic_dnmb(n_genomes = 4L, n_ogs = 4L, seed = 3L)
  out <- tempfile("sota_")
  dir.create(out)

  # Keep bootstraps small so the test stays fast but still exercises
  # the support-gating code path.
  res <- run_phylogenomics_sota(
    dnmb, out,
    min_seqs = 2L, max_seqs = 50L,
    threads = 1L, workers = 1L,
    bootstrap = 5L,
    og_bootstrap_reps = 5L,
    min_dup_support = 50,
    plot = TRUE,
    report = FALSE,
    verbose = FALSE
  )

  expect_named(res, c("og_result", "species_tree", "hogs", "relationships",
                      "dtl", "refined", "overlay_pdf"),
               ignore.order = TRUE)
  expect_true(file.exists(file.path(out, "species_tree_rooted.nwk")))
  expect_true(file.exists(file.path(out, "dtl_per_og.tsv")))
  expect_s3_class(res$species_tree, "phylo")
  if (!is.null(res$overlay_pdf)) {
    expect_true(file.exists(res$overlay_pdf))
  }
})
