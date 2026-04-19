mk_fa <- function(path, names, seqs) {
  ln <- character(2L * length(names))
  ln[seq(1L, length(ln), by = 2L)] <- paste0(">", names)
  ln[seq(2L, length(ln), by = 2L)] <- seqs
  writeLines(ln, path)
}


test_that("back_translate threads codons under an AA MSA with gaps", {
  td <- tempfile("bt_"); dir.create(td)
  og <- file.path(td, "orthogroups", "OG_0000001")
  dir.create(og, recursive = TRUE)

  # AA alignment with a gap column:
  #   seq1  M K - P
  #   seq2  M - A P
  mk_fa(file.path(og, "aln.fasta"),
        c("p1_g1", "p2_g2"),
        c("MK-P", "M-AP"))
  # CDS: one stop codon appended to exercise the trim path.
  mk_fa(file.path(td, "cds.fa"),
        c("p1_g1", "p2_g2"),
        c("ATGAAACCGTAA",       # M K P + stop
          "ATGGCTCCGTAA"))      # M A P + stop

  og_result <- tibble::tibble(
    cluster_id = 1L,
    aln_path   = file.path(og, "aln.fasta")
  )
  res <- back_translate_og_codons(og_result,
                                  cds_source = file.path(td, "cds.fa"),
                                  write_phylip = TRUE)
  expect_true(res$ok)
  expect_equal(res$n_seqs, 2L)
  expect_equal(res$n_codons, 4L)

  nuc <- readLines(res$nuc_path)
  # Row order follows input. Expect "ATGAAA---CCG" and "ATG---GCTCCG".
  seqs <- nuc[!startsWith(nuc, ">")]
  expect_equal(seqs[1], "ATGAAA---CCG")
  expect_equal(seqs[2], "ATG---GCTCCG")

  phy <- readLines(res$phy_path)
  expect_match(phy[1], "^ *2 +12$")   # 2 seqs × 12 bp
  expect_true(any(grepl("ATGAAA---CCG", phy)))
})


test_that("back_translate accepts bacterial alt-start codons at column 1", {
  td <- tempfile("bt_alt_"); dir.create(td)
  og <- file.path(td, "OG"); dir.create(og)
  mk_fa(file.path(og, "aln.fasta"),
        c("p1_g1", "p2_g2"),
        c("MK", "MK"))
  # p1_g1 starts with TTG (alt-start → M at position 1 under code 11);
  # p2_g2 uses canonical ATG. Both encode K at column 2.
  mk_fa(file.path(td, "cds.fa"),
        c("p1_g1", "p2_g2"),
        c("TTGAAA", "ATGAAA"))

  og_result <- tibble::tibble(
    cluster_id = 1L,
    aln_path   = file.path(og, "aln.fasta")
  )
  res <- back_translate_og_codons(og_result,
                                  cds_source = file.path(td, "cds.fa"),
                                  genetic_code = 11L,
                                  write_phylip = FALSE)
  expect_true(res$ok)
  nuc <- readLines(res$nuc_path)
  seqs <- nuc[!startsWith(nuc, ">")]
  expect_equal(seqs, c("TTGAAA", "ATGAAA"))
})


test_that("back_translate flags CDS/AA length mismatch rather than silently padding", {
  td <- tempfile("bt_mis_"); dir.create(td)
  og <- file.path(td, "OG"); dir.create(og)
  mk_fa(file.path(og, "aln.fasta"),
        c("p1_g1", "p2_g2"),
        c("MKP", "MAP"))
  # p1_g1 has one codon missing → 6 nt for 3 AA.
  mk_fa(file.path(td, "cds.fa"),
        c("p1_g1", "p2_g2"),
        c("ATGAAA",       # only 2 codons
          "ATGGCTCCG"))

  og_result <- tibble::tibble(
    cluster_id = 1L,
    aln_path   = file.path(og, "aln.fasta")
  )
  res <- back_translate_og_codons(og_result,
                                  cds_source = file.path(td, "cds.fa"),
                                  write_phylip = FALSE)
  expect_false(res$ok)
  expect_match(res$reason, "length mismatch|divisible")
  expect_true(is.na(res$nuc_path))
})


test_that("back_translate catches mid-alignment codon/AA disagreement", {
  td <- tempfile("bt_bad_"); dir.create(td)
  og <- file.path(td, "OG"); dir.create(og)
  mk_fa(file.path(og, "aln.fasta"),
        c("p1_g1"),
        c("MKP"))
  # CDS encodes M K S (AGT at pos 3), but AA says P — mismatch at col 3.
  mk_fa(file.path(td, "cds.fa"),
        c("p1_g1"),
        c("ATGAAAAGT"))

  og_result <- tibble::tibble(
    cluster_id = 1L,
    aln_path   = file.path(og, "aln.fasta")
  )
  res <- back_translate_og_codons(og_result,
                                  cds_source = file.path(td, "cds.fa"),
                                  write_phylip = FALSE)
  expect_false(res$ok)
  expect_match(res$reason, "mismatch")
})
