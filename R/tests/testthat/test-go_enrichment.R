test_that("go_enrichment Fisher fallback detects over-represented terms", {
  set.seed(10)
  pg <- tibble::tibble(
    protein_uid = c(
      rep(1:30, each = 1),      # GO:0001 annotates 30 proteins
      rep(31:60, each = 1),     # GO:0002 annotates 30 proteins
      rep(1:5,  each = 1)       # GO:0003 annotates 5 only
    ),
    go_id = c(
      rep("GO:0001", 30),
      rep("GO:0002", 30),
      rep("GO:0003", 5)
    )
  )
  universe <- 1:100
  # Target heavily enriched for GO:0001
  target <- c(1:25, 90:92)

  res <- go_enrichment(pg, target, universe,
                       engine = "fisher",
                       p_cutoff = 1,
                       min_annotated = 3L)
  expect_true(nrow(res) >= 1L)
  top <- res$go_id[1]
  expect_equal(top, "GO:0001")
  expect_true(res$p_value[1] < 0.05)
  expect_equal(res$engine[1], "fisher")
})

test_that("go_enrichment handles empty target gracefully", {
  pg <- tibble::tibble(protein_uid = 1:10, go_id = rep("GO:X", 10))
  res <- go_enrichment(pg, target_uids = integer(0),
                       engine = "fisher")
  expect_equal(nrow(res), 0L)
})
