test_that("score_foldseek_hog uses unordered-pair denominator and keeps HOGs with no hits", {
  td <- tempfile("fsk_"); dir.create(td)
  tsv <- file.path(td, "foldseek_aln.tsv")
  writeLines(c(
    paste("query", "target", "prob", "fident", "alnlen", "evalue", "bits",
          "alntmscore", sep = "\t"),
    # HOG A (2 members): reciprocal high-prob hit → one supported pair.
    paste("p1", "p2", "0.95", "0.80", "200", "1e-20", "300", "0.90", sep = "\t"),
    paste("p2", "p1", "0.94", "0.80", "200", "1e-20", "300", "0.90", sep = "\t"),
    # HOG B (3 members → 3 expected pairs):
    #   - p3/p4 reciprocal, both > 0.9 → supported
    #   - p4/p5 one-way-only with prob 0.5 → observed but below cutoff
    #   - p3/p5 never shows up → missing pair (still counts in denom)
    paste("p3", "p4", "0.92", "0.75", "180", "1e-18", "260", "0.82", sep = "\t"),
    paste("p4", "p3", "0.95", "0.78", "180", "1e-18", "280", "0.84", sep = "\t"),
    paste("p4", "p5", "0.50", "0.60", "150", "1e-10", "180", "0.55", sep = "\t"),
    # cross-HOG noise (should be filtered out).
    paste("p1", "p3", "0.60", "0.55", "120", "1e-08", "150", "0.40", sep = "\t")
  ), tsv)

  hog_table <- tibble::tibble(
    hog_id      = c("A", "A", "B", "B", "B", "C", "C"),
    protein_uid = c(1L, 2L, 3L, 4L, 5L, 10L, 11L)
  )
  res <- score_foldseek_hog(tsv, hog_table, prob_cutoff = 0.9)

  expect_equal(sort(res$per_hog$hog_id), c("A", "B", "C"))

  a <- res$per_hog[res$per_hog$hog_id == "A", ]
  expect_equal(a$n_expected_pairs, 1L)
  expect_equal(a$n_pairs,          1L)   # reciprocal hit collapsed
  expect_equal(a$n_supported,      1L)
  expect_equal(a$support_frac,     1)

  b <- res$per_hog[res$per_hog$hog_id == "B", ]
  expect_equal(b$n_expected_pairs, 3L)
  expect_equal(b$n_pairs,          2L)   # (p3,p4) + (p4,p5); (p3,p5) missing
  expect_equal(b$n_supported,      1L)
  expect_equal(b$support_frac,     1 / 3)

  # HOG C never shows up in the hits table. Under the old directed-hit
  # counter it was dropped; the new denominator must keep it at frac 0.
  c_ <- res$per_hog[res$per_hog$hog_id == "C", ]
  expect_equal(c_$n_expected_pairs, 1L)
  expect_equal(c_$n_pairs,          0L)
  expect_equal(c_$n_supported,      0L)
  expect_equal(c_$support_frac,     0)
})

test_that("score_foldseek_hog returns honest zeros when the alignment TSV is empty", {
  td <- tempfile("fsk_empty_"); dir.create(td)
  tsv <- file.path(td, "foldseek_aln.tsv")
  writeLines(paste("query", "target", "prob", "fident", "alnlen",
                   "evalue", "bits", "alntmscore", sep = "\t"), tsv)

  hog_table <- tibble::tibble(hog_id = c("A", "A"), protein_uid = c(1L, 2L))
  res <- score_foldseek_hog(tsv, hog_table, prob_cutoff = 0.9)
  expect_equal(nrow(res$per_hog), 1L)
  expect_equal(res$per_hog$n_expected_pairs, 1L)
  expect_equal(res$per_hog$n_supported,       0L)
  expect_equal(res$per_hog$support_frac,      0)
})

test_that("annotate_hogs_with_foldseek rewrites HOG TSVs with support cols", {
  out <- tempfile("fsk_annot_"); dir.create(out)
  dir.create(file.path(out, "HOGs"))
  hogs_path <- file.path(out, "HOGs", "N001_clade.tsv")
  writeLines(c(
    "# node_id=5",
    "# clade_leaves=g1,g2,g3",
    "hog_id\tcluster_id\tn_members\tmember_leaves\tmember_genome_keys",
    "N001.HOG00001\t42\t3\tp1_g1,p2_g2,p3_g3\tg1,g2,g3",
    "N001.HOG00002\t43\t2\tp4_g1,p5_g2\tg1,g2"
  ), hogs_path)

  per_hog <- tibble::tibble(
    hog_id           = c("N001.HOG00001", "N001.HOG00002"),
    n_pairs          = c(3L, 1L),
    n_supported      = c(3L, 0L),
    n_expected_pairs = c(3L, 1L),
    support_frac     = c(1.0, 0.0),
    mean_tmscore     = c(0.88, 0.42)
  )
  rewritten <- annotate_hogs_with_foldseek(out, per_hog)
  expect_length(rewritten, 1L)

  new_lines <- readLines(hogs_path)
  header <- new_lines[grep("^hog_id", new_lines)]
  expect_true(grepl("fsk_support_frac", header))
  expect_true(grepl("n_fsk_expected",   header))

  first_row <- strsplit(new_lines[length(new_lines) - 1L], "\t")[[1]]
  expect_equal(first_row[1], "N001.HOG00001")
  # Original 5 cols + n_fsk_pairs + n_fsk_expected + fsk_support_frac + fsk_mean_tmscore
  expect_equal(as.integer(first_row[6]), 3L)
  expect_equal(as.integer(first_row[7]), 3L)
  expect_equal(as.numeric(first_row[8]), 1.0)
  expect_equal(as.numeric(first_row[9]), 0.88, tolerance = 1e-2)
})

test_that("run_foldseek_hog skips when foldseek is absent", {
  skip_if(nzchar(Sys.which("foldseek")), "foldseek installed; skip absence test")
  td <- tempfile("fsk_in_"); dir.create(td)
  out <- tempfile("fsk_out_")
  dummy <- tibble::tibble(hog_id = character(), protein_uid = integer())
  expect_warning(
    res <- run_foldseek_hog(dnmb = list(), hog_table = dummy,
                             input_dir = td, out_dir = out),
    "foldseek"
  )
  expect_null(res)
})
