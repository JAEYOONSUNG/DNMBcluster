test_that("parse_collinearity reads blocks and pairs correctly", {
  tmp <- tempfile(fileext = ".collinearity")
  writeLines(c(
    "############### Parameters ###############",
    "# MATCH_SCORE: 50",
    "## Alignment 0: score=250.0 e_value=1e-20 N=5 g1&g2 plus",
    "  0-  0: p1_g1   p11_g2   1e-30",
    "  0-  1: p2_g1   p12_g2   1e-28",
    "  0-  2: p3_g1   p13_g2   1e-25",
    "## Alignment 1: score=90.0 e_value=1e-5 N=3 g1&g3 minus",
    "  1-  0: p5_g1   p21_g3   1e-10",
    "  1-  1: p6_g1   p22_g3   1e-09"
  ), tmp)

  out <- parse_collinearity(tmp)
  expect_equal(nrow(out$blocks), 2L)
  expect_equal(out$blocks$sp_a, c("g1", "g1"))
  expect_equal(out$blocks$sp_b, c("g2", "g3"))
  expect_equal(out$blocks$orientation, c("plus", "minus"))
  expect_equal(sum(out$pairs$block_id == 0L), 3L)
  expect_equal(sum(out$pairs$block_id == 1L), 2L)
  expect_true(all(grepl("^p", out$pairs$gene_a)))
})

test_that("run_mcscanx skips with a warning when MCScanX is absent", {
  skip_if(nzchar(Sys.which("MCScanX")), "MCScanX installed; skip absence test")
  td <- tempfile("mcsx_"); dir.create(td)
  expect_warning(res <- run_mcscanx(td), "MCScanX")
  expect_null(res)
})
