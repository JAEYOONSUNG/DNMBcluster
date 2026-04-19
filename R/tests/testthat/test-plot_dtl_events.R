test_that("plot_dtl_events aggregates counts per species-tree node", {
  skip_if_not_installed("ape")

  sp <- ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
  events <- tibble::tibble(
    cluster_id = c(1L, 1L, 2L, 2L),
    g_node     = c(6L, 7L, 6L, 7L),
    g_n_leaves = c(2L, 2L, 2L, 2L),
    sp_lca     = c(6L, 7L, 6L, 5L),
    event      = c("S", "D", "T", "S")
  )
  losses <- tibble::tibble(sp_parent = 5L, sp_child = 1L, loss_count = 2L)

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  res <- plot_dtl_events(events, sp, losses)

  expect_s3_class(res, "tbl_df")
  expect_equal(sum(res$n_S), 2L)
  expect_equal(sum(res$n_D), 1L)
  expect_equal(sum(res$n_T), 1L)
  expect_equal(sum(res$n_loss_in), 2L)
  expect_equal(res$n_loss_in[1], 2L)
})
