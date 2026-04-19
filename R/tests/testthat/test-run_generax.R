test_that("parse_generax_events tallies D/T/L/S per family", {
  td <- tempfile("grx_"); dir.create(td)
  rec <- file.path(td, "reconciliations"); dir.create(rec)
  # Family 42: 3 speciations, 2 duplications, 1 transfer, 1 loss
  writeLines(c(
    "S node=1 species=g1",
    "D node=2 species=g1",
    "S node=3 species=g2",
    "T node=4 src=g1 dst=g3",
    "S node=5 species=g3",
    "D node=6 species=g4",
    "L node=7 species=g5"
  ), file.path(rec, "OG_0000042_events.txt"))
  # Family 7: all speciations
  writeLines(c("S", "S", "S"), file.path(rec, "OG_0000007_events.txt"))

  out <- parse_generax_events(td)
  expect_equal(sort(unique(out$per_og$cluster_id)), c(7L, 42L))

  row42 <- out$per_og[out$per_og$cluster_id == 42L, ]
  expect_equal(row42$speciations,  3L)
  expect_equal(row42$duplications, 2L)
  expect_equal(row42$transfers,    1L)
  expect_equal(row42$losses,       1L)

  row7 <- out$per_og[out$per_og$cluster_id == 7L, ]
  expect_equal(row7$speciations, 3L)
  expect_equal(row7$duplications, 0L)
})

test_that("run_generax skips with a warning when generax is absent", {
  skip_if(nzchar(Sys.which("generax")), "generax installed; skip absence test")
  td <- tempfile("grx_missing_"); dir.create(td)
  writeLines("", file.path(td, "families.txt"))
  writeLines("(A,B);", file.path(td, "species_tree.nwk"))
  out <- tempfile("grx_out_")
  expect_warning(res <- run_generax(td, out), "generax")
  expect_null(res)
})
