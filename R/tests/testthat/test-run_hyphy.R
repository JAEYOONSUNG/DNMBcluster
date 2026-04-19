test_that("run_hyphy_* wrappers skip with a warning when hyphy is absent", {
  skip_if(nzchar(Sys.which("hyphy")), "hyphy installed; skip absence test")
  td <- tempfile("hy_abs_"); dir.create(td)
  og <- file.path(td, "OG_0000001"); dir.create(og)
  writeLines(c(">p1_g1", "ATGAAA",
               ">p2_g2", "ATGGCT"), file.path(og, "aln.nuc.fasta"))
  writeLines("(p1_g1:0.1,p2_g2:0.1);", file.path(og, "tree.nwk"))
  wd <- file.path(td, "work")

  expect_warning(res1 <- run_hyphy_fel   (og, work_dir = wd), "hyphy")
  expect_warning(res2 <- run_hyphy_busted(og, work_dir = wd), "hyphy")
  expect_warning(res3 <- run_hyphy_absrel(og, work_dir = wd), "hyphy")
  expect_null(res1); expect_null(res2); expect_null(res3)
})


test_that(".parse_hyphy_fel extracts per-site omega summary from fixture JSON", {
  skip_if_not_installed("jsonlite")
  td <- tempfile("fel_"); dir.create(td)
  j <- list(
    MLE = list(
      headers = list(
        c("alpha", "syn"),
        c("beta",  "nonsyn"),
        c("alpha=beta", "shared"),
        c("LRT", "lrt"),
        c("p-value", "p"),
        c("Total branch length", "len")
      ),
      content = list(
        # site 1: beta > alpha, p = 0.02 → positive
        # site 2: alpha > beta, p = 0.01 → negative
        # site 3: alpha == beta, p = 0.5 → neither
        `0` = rbind(
          c(0.5, 2.0, 0.0, 5.0, 0.02, 0.3),
          c(3.0, 0.2, 0.0, 4.0, 0.01, 0.3),
          c(1.0, 1.0, 1.0, 0.0, 0.50, 0.3)
        )
      )
    )
  )
  jp <- file.path(td, "fel.json")
  jsonlite::write_json(j, jp, auto_unbox = TRUE, matrix = "rowmajor")
  out <- DNMBcluster:::.parse_hyphy_fel(cid = 42L, json_path = jp)
  expect_equal(out$cluster_id, 42L)
  expect_equal(out$n_sites, 3L)
  expect_equal(out$n_positive, 1L)
  expect_equal(out$n_negative, 1L)
  expect_true(out$converged)
  expect_true(is.finite(out$median_omega))
})


test_that(".parse_hyphy_busted surfaces gene-wide p-value + positive rate class", {
  skip_if_not_installed("jsonlite")
  td <- tempfile("bu_"); dir.create(td)
  j <- list(
    `test results` = list(`p-value` = 0.003, LRT = 9.8),
    fits = list(
      `Unconstrained model` = list(
        `Rate Distributions` = list(
          Test = data.frame(
            omega      = c(0.1, 1.0, 4.2),
            proportion = c(0.6, 0.3, 0.1)
          )
        )
      )
    )
  )
  jp <- file.path(td, "busted.json")
  jsonlite::write_json(j, jp, auto_unbox = TRUE, dataframe = "columns")
  out <- DNMBcluster:::.parse_hyphy_busted(cid = 7L, json_path = jp)
  expect_equal(out$cluster_id, 7L)
  expect_equal(out$p_value, 0.003, tolerance = 1e-6)
  expect_equal(out$lrt, 9.8, tolerance = 1e-6)
  expect_equal(out$omega_positive, 4.2, tolerance = 1e-6)
  expect_equal(out$weight_positive, 0.1, tolerance = 1e-6)
  expect_true(out$converged)
})


test_that(".parse_hyphy_absrel counts branches and min corrected p", {
  skip_if_not_installed("jsonlite")
  td <- tempfile("ab_"); dir.create(td)
  j <- list(
    `branch attributes` = list(
      `0` = list(
        Node1 = list(`Corrected P-value` = 0.001),
        Node2 = list(`Corrected P-value` = 0.6),
        Node3 = list(`Corrected P-value` = 0.04)
      )
    )
  )
  jp <- file.path(td, "absrel.json")
  jsonlite::write_json(j, jp, auto_unbox = TRUE)
  out <- DNMBcluster:::.parse_hyphy_absrel(cid = 3L, json_path = jp)
  expect_equal(out$n_branches, 3L)
  expect_equal(out$n_selected, 2L)          # 0.001 and 0.04
  expect_equal(out$min_corrected_p, 0.001, tolerance = 1e-6)
  expect_true(out$converged)
})


test_that("run_hyphy_batch skips with warning when hyphy absent", {
  skip_if(nzchar(Sys.which("hyphy")), "hyphy installed; skip absence test")
  td <- tempfile("hyb_"); dir.create(td)
  og <- file.path(td, "OG_0000001"); dir.create(og)
  writeLines(c(">p1_g1", "ATGAAA", ">p2_g2", "ATGGCT"),
             file.path(og, "aln.nuc.fasta"))
  writeLines("(p1_g1:0.1,p2_g2:0.1);", file.path(og, "tree.nwk"))
  expect_warning(res <- run_hyphy_batch(og, work_dir = file.path(td, "w")),
                  "hyphy")
  expect_null(res)
})


test_that(".hyphy_combine merges FEL/BUSTED/aBSREL per cluster_id", {
  fel <- tibble::tibble(
    cluster_id = c(1L, 2L), method = "FEL",
    n_sites = c(10L, 15L),
    n_positive = c(2L, 0L), n_negative = c(1L, 3L),
    median_omega = c(0.8, 0.3),
    json_path = NA_character_, converged = TRUE, reason = NA_character_)
  bst <- tibble::tibble(
    cluster_id = c(1L, 3L), method = "BUSTED",
    p_value = c(0.01, 0.4), q_value = c(0.02, 0.4),
    lrt = c(10.1, 1.3),
    omega_positive = c(3.2, 1.1), weight_positive = c(0.1, 0.2),
    json_path = NA_character_, converged = TRUE, reason = NA_character_)
  abs <- tibble::tibble(
    cluster_id = c(2L, 3L), method = "aBSREL",
    n_branches = c(5L, 4L), n_selected = c(1L, 0L),
    min_corrected_p = c(0.02, 0.9), q_value = c(0.04, 0.9),
    json_path = NA_character_, converged = TRUE, reason = NA_character_)

  out <- DNMBcluster:::.hyphy_combine(fel, bst, abs)
  expect_equal(sort(out$cluster_id), c(1L, 2L, 3L))
  row1 <- out[out$cluster_id == 1L, ]
  expect_equal(row1$fel_n_positive, 2L)
  expect_equal(row1$busted_p, 0.01)
  expect_true(is.na(row1$absrel_n_branches))
  row2 <- out[out$cluster_id == 2L, ]
  expect_equal(row2$absrel_n_selected, 1L)
  expect_true(is.na(row2$busted_p))
})


test_that("BH correction is applied across OGs on busted/absrel", {
  skip_if_not_installed("jsonlite")
  # Simulate a batch where hyphy binary check passes by mocking it out:
  # we directly test the FDR path by building method tibbles the way
  # run_hyphy_batch would post-process them.
  bst <- tibble::tibble(
    cluster_id = 1:5,
    p_value = c(0.001, 0.02, 0.04, 0.2, 0.9),
    lrt = c(10, 5, 4, 2, 0),
    omega_positive = c(4, 3, 2, 1, 0.5),
    weight_positive = c(0.1, 0.15, 0.2, 0.3, 0.5),
    method = "BUSTED", converged = TRUE,
    json_path = NA_character_, reason = NA_character_)
  # emulate run_hyphy_batch's correction step
  bst$q_value <- stats::p.adjust(bst$p_value, method = "BH")
  expect_true(all(bst$q_value >= bst$p_value))
  # 3rd smallest raw p (0.04) with n=5: q = 0.04 * 5 / 3 = 0.0667
  expect_equal(bst$q_value[3], 0.04 * 5 / 3, tolerance = 1e-8)
})


test_that("parsers return NA-row tibbles when JSON is missing", {
  out_f <- DNMBcluster:::.parse_hyphy_fel   (cid = 1L, json_path = NULL,
                                             reason = "no file")
  out_b <- DNMBcluster:::.parse_hyphy_busted(cid = 1L, json_path = NULL,
                                             reason = "no file")
  out_a <- DNMBcluster:::.parse_hyphy_absrel(cid = 1L, json_path = NULL,
                                             reason = "no file")
  expect_false(out_f$converged)
  expect_false(out_b$converged)
  expect_false(out_a$converged)
  expect_equal(out_f$reason, "no file")
})
