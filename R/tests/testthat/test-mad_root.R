test_that("mad_root returns a rooted tree on random unrooted input", {
  testthat::skip_if_not_installed("ape")
  set.seed(17)
  tr <- ape::rtree(10, rooted = FALSE)
  tr <- ape::unroot(tr)
  r <- mad_root(tr)
  expect_true(ape::is.rooted(r))
  expect_equal(sort(r$tip.label), sort(tr$tip.label))
})

test_that("mad_root falls back to midpoint when edge lengths are absent", {
  testthat::skip_if_not_installed("ape")
  tr <- ape::read.tree(text = "((A,B),(C,D));")
  r <- mad_root(tr)
  expect_true(ape::is.rooted(r))
})
