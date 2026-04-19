test_that("reconcile_dtl labels a clean duplication as D", {
  testthat::skip_if_not_installed("ape")
  sp <- ape::read.tree(text = "((A,B),(C,D));")
  # Gene tree with two A-lineage copies — classic duplication.
  gt <- ape::read.tree(text = "(((p1_gA,p2_gA),p3_gB),(p4_gC,p5_gD));")
  res <- reconcile_dtl(gt, sp)
  expect_equal(nrow(res$events), gt$Nnode)
  expect_gte(res$summary$n_D, 1)
})

test_that("reconcile_dtl labels a clean speciation history as all-S", {
  testthat::skip_if_not_installed("ape")
  sp <- ape::read.tree(text = "((A,B),(C,D));")
  gt <- ape::read.tree(text = "((p1_gA,p2_gB),(p3_gC,p4_gD));")
  res <- reconcile_dtl(gt, sp)
  expect_equal(res$summary$n_D, 0L)
  expect_equal(res$summary$n_T, 0L)
  expect_gte(res$summary$n_S, gt$Nnode - 1L)
})

test_that("reconcile_dtl tolerates missing species mapping via label parsing", {
  testthat::skip_if_not_installed("ape")
  sp <- ape::read.tree(text = "(A,B);")
  gt <- ape::read.tree(text = "(p1_gA,p2_gB);")
  expect_error(reconcile_dtl(gt, sp), NA)
})

test_that("min_dup_support downgrades low-support D/T to S", {
  testthat::skip_if_not_installed("ape")
  sp <- ape::read.tree(text = "((A,B),(C,D));")
  gt <- ape::read.tree(text = "(((p1_gA,p2_gA)20,p3_gB)80,(p4_gC,p5_gD)90);")
  lo <- reconcile_dtl(gt, sp, min_dup_support = 0)
  hi <- reconcile_dtl(gt, sp, min_dup_support = 50)
  expect_gte(lo$summary$n_D, 1L)
  expect_equal(hi$summary$n_D, 0L)
  expect_gte(hi$summary$n_S, lo$summary$n_S)
})

test_that("long_branch_quantile gate downgrades short-branch T to D", {
  testthat::skip_if_not_installed("ape")
  sp <- ape::read.tree(text = "((A,B),(C,D));")
  gt <- ape::read.tree(
    text = "(((p1_gA:0.1,p2_gC:0.1):0.1,p3_gB:0.1):0.1,(p4_gC:0.1,p5_gD:0.1):0.1);"
  )
  loose <- reconcile_dtl(gt, sp, long_branch_quantile = 1)
  strict <- reconcile_dtl(gt, sp, long_branch_quantile = 0.99)
  expect_true(loose$summary$n_T >= strict$summary$n_T)
})
