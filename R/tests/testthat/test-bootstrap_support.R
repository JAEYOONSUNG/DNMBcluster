test_that("per_og_trees populates tree$node.label with bootstrap_reps > 0 (NJ)", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()

  dnmb <- make_synthetic_dnmb(n_genomes = 5L, n_ogs = 2L, seed = 3L)
  out  <- tempfile("per_og_bs_")
  dir.create(out)

  res <- per_og_trees(dnmb, out, min_seqs = 2L, max_seqs = 50L,
                      method = "nj", threads = 1L, verbose = FALSE,
                      bootstrap_reps = 20L)

  expect_true(all(file.exists(res$tree_path)))
  tr <- ape::read.tree(res$tree_path[1])
  # Bootstrap should leave at least one non-empty node label on a 5-tip tree.
  expect_false(is.null(tr$node.label))
  expect_true(any(nzchar(tr$node.label)))
  # Support values, when parseable, must be in [0, 100].
  nums <- suppressWarnings(as.integer(tr$node.label))
  nums <- nums[!is.na(nums)]
  if (length(nums)) {
    expect_true(all(nums >= 0L & nums <= 100L))
  }
})

test_that("per_og_trees omits node.label when bootstrap_reps = 0", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()

  dnmb <- make_synthetic_dnmb(n_genomes = 4L, n_ogs = 1L, seed = 4L)
  out  <- tempfile("per_og_nobs_")
  dir.create(out)

  res <- per_og_trees(dnmb, out, min_seqs = 2L, max_seqs = 50L,
                      method = "nj", threads = 1L, verbose = FALSE,
                      bootstrap_reps = 0L)

  tr <- ape::read.tree(res$tree_path[1])
  # Either NULL or all-empty node labels is acceptable.
  expect_true(is.null(tr$node.label) || all(!nzchar(tr$node.label)))
})
