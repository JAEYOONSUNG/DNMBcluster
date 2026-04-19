test_that("plot_busted_volcano returns a ggplot with expected layers", {
  skip_if_not_installed("ggplot2")
  bst <- tibble::tibble(
    cluster_id = 1:6,
    p_value = c(0.001, 0.01, 0.03, 0.2, 0.6, NA),
    lrt = c(10, 7, 5, 2, 0.3, NA),
    omega_positive = c(4.2, 2.5, 1.5, 0.9, 0.3, NA),
    weight_positive = c(0.1, 0.15, 0.2, 0.3, 0.5, NA),
    method = "BUSTED", converged = c(rep(TRUE, 5), FALSE),
    json_path = NA_character_, reason = NA_character_)
  p <- plot_busted_volcano(bst, q_threshold = 0.05, label_top_n = 2L)
  expect_s3_class(p, "ggplot")
  expect_true(any(grepl("volcano", tolower(p$labels$title))))
  # Aesthetic mapping carries the expected numeric columns.
  expect_true("colour" %in% names(p$mapping))
})


test_that("plot_busted_volcano handles empty/NA-only input", {
  skip_if_not_installed("ggplot2")
  empty <- tibble::tibble(
    cluster_id = integer(0), p_value = numeric(0),
    omega_positive = numeric(0), weight_positive = numeric(0),
    method = character(0), converged = logical(0),
    json_path = character(0), reason = character(0))
  p1 <- plot_busted_volcano(empty)
  expect_s3_class(p1, "ggplot")
  p2 <- plot_busted_volcano(NULL)
  expect_s3_class(p2, "ggplot")
})


test_that("plot_fel_sites builds a ggplot from JSON path or parsed list", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("jsonlite")
  td <- tempfile("fel_plot_"); dir.create(td)
  j <- list(
    MLE = list(
      headers = list(c("alpha","s"), c("beta","ns"), c("eq","sh"),
                     c("LRT","l"), c("p-value","p"), c("len","tbl")),
      content = list(
        `0` = rbind(
          c(0.5, 2.0, 0.0, 5.0, 0.02, 0.3),
          c(3.0, 0.2, 0.0, 4.0, 0.01, 0.3),
          c(1.0, 1.0, 1.0, 0.0, 0.50, 0.3),
          c(0.2, 0.4, 0.0, 1.0, 0.30, 0.3)
        )
      )
    )
  )
  jp <- file.path(td, "fel.json")
  jsonlite::write_json(j, jp, auto_unbox = TRUE, matrix = "rowmajor")
  p_path <- plot_fel_sites(jp)
  expect_s3_class(p_path, "ggplot")
  # Also accept a parsed list.
  parsed <- jsonlite::fromJSON(jp, simplifyVector = TRUE)
  p_list <- plot_fel_sites(parsed)
  expect_s3_class(p_list, "ggplot")
})


test_that("plot_codeml_vs_hyphy restricts to overlapping cluster_ids", {
  skip_if_not_installed("ggplot2")
  cml <- tibble::tibble(cluster_id = c(1L, 2L, 3L),
                         omega = c(0.3, 0.9, 2.1))
  bst <- tibble::tibble(cluster_id = c(2L, 3L, 4L),
                         omega_positive = c(1.0, 3.0, 0.5),
                         p_value = c(0.5, 0.01, 0.8))
  p <- plot_codeml_vs_hyphy(cml, bst)
  expect_s3_class(p, "ggplot")
  # inner_join of cluster_id 1/2/3 and 2/3/4 → 2/3 only.
  expect_equal(nrow(p$data), 2L)
  # Empty-overlap case returns an informative empty plot, not an error.
  cml2 <- tibble::tibble(cluster_id = integer(0), omega = numeric(0))
  bst2 <- tibble::tibble(cluster_id = integer(0), omega_positive = numeric(0),
                         p_value = numeric(0))
  expect_s3_class(plot_codeml_vs_hyphy(cml2, bst2), "ggplot")
})


test_that("plot_selection_overview tolerates missing FEL or BUSTED", {
  skip_if_not_installed("ggplot2")
  sel <- list(
    busted = tibble::tibble(
      cluster_id = 1:3, p_value = c(0.01, 0.3, 0.9),
      lrt = c(5, 1, 0.1), omega_positive = c(2.5, 1.1, 0.4),
      weight_positive = c(0.2, 0.3, 0.4),
      method = "BUSTED", converged = TRUE,
      json_path = NA_character_, reason = NA_character_),
    fel    = NULL,
    absrel = NULL
  )
  p <- plot_selection_overview(sel)
  # Either a single ggplot (no patchwork) or a patchwork object.
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))

  # Fully empty input returns an empty-placeholder plot.
  expect_s3_class(plot_selection_overview(NULL), "ggplot")
})
