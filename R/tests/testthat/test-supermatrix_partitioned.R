test_that("supermatrix_species_tree returns a valid tree with NJ on SC-OGs", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()

  dnmb <- make_synthetic_dnmb(n_genomes = 4L, n_ogs = 3L, seed = 7L)
  out  <- tempfile("sm_nj_")
  dir.create(out)
  og <- per_og_trees(dnmb, out, min_seqs = 2L, max_seqs = 50L,
                     method = "nj", threads = 1L, verbose = FALSE)

  tr <- supermatrix_species_tree(og, dnmb, method = "nj")
  expect_s3_class(tr, "phylo")
  expect_equal(sort(tr$tip.label), sort(dnmb$genome_meta$genome_key))
})

test_that("supermatrix_species_tree with partitioned=TRUE returns a tree or falls back", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()

  dnmb <- make_synthetic_dnmb(n_genomes = 4L, n_ogs = 3L, seed = 9L)
  out  <- tempfile("sm_part_")
  dir.create(out)
  og <- per_og_trees(dnmb, out, min_seqs = 2L, max_seqs = 50L,
                     method = "nj", threads = 1L, verbose = FALSE)

  tr <- suppressWarnings(
    supermatrix_species_tree(og, dnmb, method = "ml",
                             partitioned = TRUE, gamma_k = 0L)
  )
  expect_s3_class(tr, "phylo")
  expect_equal(length(tr$tip.label), nrow(dnmb$genome_meta))
})
