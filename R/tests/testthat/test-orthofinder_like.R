test_that("run_orthofinder_like round-trips a tiny synthetic dataset", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()

  dnmb <- make_synthetic_dnmb(n_genomes = 4L, n_ogs = 5L)
  out  <- tempfile("of_like_")
  dir.create(out)

  res <- run_orthofinder_like(
    dnmb, out,
    method = "nj", min_seqs = 2L, max_seqs = 50L,
    threads = 1L, verbose = FALSE
  )

  expect_named(res, c("og_result", "species_tree", "hogs",
                      "relationships", "dtl", "refined"))
  expect_s3_class(res$species_tree, "phylo")
  expect_true(file.exists(file.path(out, "species_tree_rooted.nwk")))
  expect_true(nrow(res$og_result) >= 1L)
  expect_null(res$dtl)
})

test_that("run_orthofinder_like with dtl=TRUE writes per-OG DTL summary", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()

  dnmb <- make_synthetic_dnmb(n_genomes = 4L, n_ogs = 3L, seed = 2L)
  out  <- tempfile("of_like_dtl_")
  dir.create(out)

  res <- run_orthofinder_like(
    dnmb, out,
    method = "nj", min_seqs = 2L, max_seqs = 50L,
    threads = 1L, verbose = FALSE, dtl = TRUE
  )

  expect_false(is.null(res$dtl))
  expect_true(file.exists(file.path(out, "dtl_per_og.tsv")))
  expect_true(file.exists(file.path(out, "dtl_events.tsv")))
  expect_true(all(c("n_S", "n_D", "n_T", "n_loss") %in% names(res$dtl$per_og)))
})

test_that("incremental_add assigns new proteins and preserves label scheme", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("arrow")

  dnmb <- make_synthetic_dnmb(n_genomes = 3L, n_ogs = 2L)
  base <- tempfile("base_")
  dir.create(base)
  write_synthetic_dnmb_parquets(dnmb, base)
  og <- per_og_trees(dnmb, base, min_seqs = 2L, max_seqs = 50L,
                    method = "nj", threads = 1L, verbose = FALSE)

  new_seq <- dnmb$gene_table$translation[1]
  new_fa <- tempfile(fileext = ".fa")
  writeLines(c(">query_A", new_seq), new_fa)

  out_dir <- tempfile("inc_")
  dir.create(out_dir)
  res <- incremental_add(base, new_fa, out_dir,
                         min_bitscore = 10,
                         new_genome_key = "NEWK",
                         threads = 1L, method = "nj")

  expect_true(file.exists(file.path(out_dir, "incremental", "assignments.tsv")))
  expect_true(file.exists(file.path(out_dir, "incremental", "new_genome.tsv")))
  assign_tbl <- utils::read.table(
    file.path(out_dir, "incremental", "assignments.tsv"),
    sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  expect_equal(nrow(assign_tbl), 1L)
  expect_false(is.na(assign_tbl$cluster_id[1]))

  if (nrow(res) >= 1L && !is.na(res$tree_path[1])) {
    tr <- ape::read.tree(res$tree_path[1])
    new_tip <- grep(paste0("_g", max(dnmb$genome_meta$genome_uid) + 1L),
                    tr$tip.label, value = TRUE)
    expect_true(length(new_tip) >= 1L)
    expect_match(new_tip[1], "^p\\d+_g\\d+$")
  }
})
