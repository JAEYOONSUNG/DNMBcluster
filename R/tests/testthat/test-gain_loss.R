test_that("Dollo + Fitch match known ground truth on a toy tree", {
  testthat::skip_if_not_installed("ape")
  tr <- ape::read.tree(text = "(((A:1,B:1):1,(C:1,D:1):1):1,E:2);")
  pres <- cbind(
    all_present = c(1, 1, 1, 1, 1),
    single_tip  = c(0, 0, 0, 1, 0),
    clade_ABCD  = c(1, 1, 1, 1, 0),
    patchy      = c(1, 0, 0, 0, 1)  # A and E only
  )
  rownames(pres) <- c("A", "B", "C", "D", "E")

  gl_d <- reconstruct_gain_loss(tr, pres, method = "dollo", return_states = TRUE)
  gl_f <- reconstruct_gain_loss(tr, pres, method = "fitch", return_states = TRUE)

  # Core (all_present) produces 0 gains, 0 losses under both.
  expect_equal(sum(gl_d$states[, "all_present"] == 1L),
               length(tr$tip.label) + tr$Nnode)
  # single_tip: Dollo gives exactly 1 gain on edge leading into D.
  expect_equal(sum(gl_d$edges$gains[gl_d$edges$child_label == "D" & !is.na(gl_d$edges$child_label)]), 1)
  # Patchy (A + E) should trigger HGT candidate.
  hg <- detect_hgt_candidates(tr, pres)
  expect_true(hg$hgt_candidate[hg$hog_id == "patchy"])
})

test_that("Fitch vectorised output equals per-column fold", {
  testthat::skip_if_not_installed("ape")
  set.seed(11)
  tr <- ape::rtree(8, rooted = TRUE)
  pres <- matrix(sample(0:1, 8 * 20, replace = TRUE), nrow = 8)
  rownames(pres) <- tr$tip.label
  gl <- reconstruct_gain_loss(tr, pres, method = "fitch", return_states = TRUE)
  # Every internal state is 0 or 1
  expect_true(all(gl$states %in% c(0L, 1L)))
  # Edge gain + loss never exceeds n_hog
  expect_true(all(gl$edges$gains  <= ncol(pres)))
  expect_true(all(gl$edges$losses <= ncol(pres)))
})

test_that("ancestral_pan_core_sizes returns monotone non-increasing pan down the tree", {
  testthat::skip_if_not_installed("ape")
  set.seed(5)
  tr <- ape::rtree(6, rooted = TRUE)
  pres <- matrix(sample(0:1, 6 * 10, replace = TRUE), nrow = 6)
  rownames(pres) <- tr$tip.label
  gl <- reconstruct_gain_loss(tr, pres, method = "dollo", return_states = TRUE)
  sizes <- ancestral_pan_core_sizes(gl)
  expect_true(all(sizes$core_size <= sizes$pan_size))
})
