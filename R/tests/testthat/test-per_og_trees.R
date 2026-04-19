test_that("per_og_trees produces one entry per OG with SC flag set", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()

  dnmb <- make_synthetic_dnmb(n_genomes = 4L, n_ogs = 3L)
  out  <- tempfile("per_og_")
  dir.create(out)

  res <- per_og_trees(dnmb, out, min_seqs = 2L, max_seqs = 50L,
                     method = "nj", threads = 1L, verbose = FALSE)

  expect_equal(nrow(res), 3L)
  expect_true(all(!is.na(res$tree_path)))
  # Each synthetic OG has 1 member per genome → SC-OGs
  expect_true(all(res$single_copy))
  expect_true(file.exists(file.path(out, "orthogroup_trees.tsv")))
})

test_that("per_og_trees with method='iqtree' uses IQ-TREE when available", {
  skip_if_not_installed("ape")
  skip_if_no_aln_engine()
  skip_if_not(nzchar(Sys.which("iqtree2")) || nzchar(Sys.which("iqtree")),
              "iqtree binary not on PATH")

  dnmb <- make_synthetic_dnmb(n_genomes = 4L, n_ogs = 2L)
  out  <- tempfile("per_og_iq_")
  dir.create(out)

  res <- per_og_trees(dnmb, out, min_seqs = 2L, max_seqs = 50L,
                      method = "iqtree", threads = 1L, verbose = FALSE,
                      bootstrap_reps = 1000L)

  expect_true(all(!is.na(res$tree_path)))
  tr <- ape::read.tree(res$tree_path[1])
  expect_s3_class(tr, "phylo")
  # UFBoot labels should be numeric-looking when bootstrap requested.
  if (!is.null(tr$node.label) && any(nzchar(tr$node.label))) {
    nz <- tr$node.label[nzchar(tr$node.label)]
    expect_true(all(grepl("^[0-9.]+$", nz)))
  }
})

test_that("per_og_trees computes SC flag on FULL OG, not subsampled", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()

  # 4 genomes × 3 copies/genome → 12 members/OG. Subsampling to 4 would
  # falsely look like a single-copy OG if the flag were computed after.
  dnmb <- make_synthetic_dnmb(n_genomes = 4L, n_ogs = 1L)
  # duplicate each member twice more to make multi-copy
  extra <- dnmb$clusters
  puid_shift <- max(dnmb$gene_table$protein_uid)
  for (k in 1:2) {
    add_cl <- extra
    add_cl$protein_uid <- add_cl$protein_uid + puid_shift * k
    add_g  <- dnmb$gene_table
    add_g$protein_uid  <- add_g$protein_uid + puid_shift * k
    dnmb$clusters   <- dplyr::bind_rows(dnmb$clusters, add_cl)
    dnmb$gene_table <- dplyr::bind_rows(dnmb$gene_table, add_g)
  }

  out <- tempfile("per_og_mc_")
  dir.create(out)
  res <- per_og_trees(dnmb, out, min_seqs = 2L, max_seqs = 4L,
                     method = "nj", threads = 1L, verbose = FALSE)

  expect_false(res$single_copy[1])
})
