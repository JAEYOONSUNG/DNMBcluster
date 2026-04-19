test_that(".export_generax emits a parser-valid families.txt layout", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()

  dnmb <- make_synthetic_dnmb(n_genomes = 4L, n_ogs = 3L, seed = 11L)
  out  <- tempfile("exp_grx_"); dir.create(out)
  per_res <- per_og_trees(dnmb, out, min_seqs = 2L, max_seqs = 50L,
                          method = "nj", threads = 1L, verbose = FALSE)
  # Synthesize a species tree the exporter can copy alongside.
  sp_nwk <- file.path(out, "sp.nwk")
  ape::write.tree(ape::rtree(nrow(dnmb$genome_meta),
                             tip.label = dnmb$genome_meta$genome_key),
                  sp_nwk)

  bundle <- export_og_bundles(dnmb, per_res, out,
                              tools = "generax",
                              species_tree = sp_nwk,
                              aa_model = "WAG")
  fams_path <- file.path(bundle$generax, "families.txt")
  expect_true(file.exists(fams_path))
  ln <- readLines(fams_path, warn = FALSE)

  # Header must be literally [FAMILIES] on its own line.
  expect_equal(trimws(ln[1]), "[FAMILIES]")
  # Family entries must be "- <name>", not "- FAMILY: <name>".
  fam_entries <- grep("^\\s*-\\s", ln, value = TRUE)
  expect_true(length(fam_entries) > 0L)
  expect_true(all(!grepl("FAMILY:", fam_entries)))
  expect_true(all(grepl("^- OG_\\d{7}$", fam_entries)))

  # alignment / mapping / starting_gene_tree must be absolute paths that
  # exist — GeneRax resolves them relative to its own cwd.
  kv <- function(key) {
    rows <- grep(paste0("^", key, "\\s*="), ln, value = TRUE)
    trimws(sub(paste0("^", key, "\\s*=\\s*"), "", rows))
  }
  for (key in c("alignment", "mapping", "starting_gene_tree")) {
    vals <- kv(key)
    expect_true(length(vals) > 0L)
    expect_true(all(startsWith(vals, "/")), info = paste(key, "must be absolute"))
    expect_true(all(file.exists(vals)), info = paste(key, "must exist"))
  }

  # aa_model is plumbed through.
  expect_true(any(grepl("^subst_model\\s*=\\s*WAG\\+G4$", ln)))

  # species_tree.nwk copied next to families.txt for convenience.
  expect_true(file.exists(file.path(bundle$generax, "species_tree.nwk")))

  # mapping.link uses phyldog "species:gene1;gene2" direction.
  fam_dirs <- list.dirs(bundle$generax, recursive = FALSE)
  fam_dirs <- fam_dirs[grepl("/OG_\\d{7}$", fam_dirs)]
  expect_true(length(fam_dirs) > 0L)
  link <- readLines(file.path(fam_dirs[1], "mapping.link"), warn = FALSE)
  expect_true(all(grepl("^[^:]+:.+$", link)))
  # Right-hand side carries at least one gene id of the p<uid>_g<genome> shape.
  rhs <- sub("^[^:]+:", "", link)
  expect_true(all(grepl("p\\d+_g\\d+", rhs)))
})

test_that(".export_paml emits per-OG dirs with aln.aa.fasta + tree + codon ctl", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()

  dnmb <- make_synthetic_dnmb(n_genomes = 3L, n_ogs = 2L, seed = 12L)
  out  <- tempfile("exp_paml_"); dir.create(out)
  per_res <- per_og_trees(dnmb, out, min_seqs = 2L, max_seqs = 50L,
                          method = "nj", threads = 1L, verbose = FALSE)
  bundle <- export_og_bundles(dnmb, per_res, out, tools = "paml")
  root <- bundle$paml

  # Per-OG dirs (not flat .fasta files — that layout is not consumable
  # by run_codeml_m0()).
  fam_dirs <- list.dirs(root, recursive = FALSE)
  fam_dirs <- fam_dirs[grepl("/OG_\\d{7}$", fam_dirs)]
  expect_true(length(fam_dirs) > 0L)
  for (fd in fam_dirs) {
    expect_true(file.exists(file.path(fd, "aln.aa.fasta")))
    expect_true(file.exists(file.path(fd, "tree.nwk")))
  }

  # Control template must be codon-oriented (seqtype = 1) so users who
  # produce a back-translated aln.nuc.phy can run codeml M0 directly.
  ctl <- readLines(file.path(root, "codeml.ctl.template"), warn = FALSE)
  expect_true(any(grepl("^\\s*seqtype\\s*=\\s*1", ctl)))
  expect_true(any(grepl("^\\s*seqfile\\s*=\\s*aln\\.nuc\\.phy", ctl)))
  expect_true(file.exists(file.path(root, "README.txt")))
})

test_that(".export_hyphy emits per-OG dirs with AA MSA + tree + HyPhy README", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()

  dnmb <- make_synthetic_dnmb(n_genomes = 3L, n_ogs = 2L, seed = 14L)
  out  <- tempfile("exp_hyphy_"); dir.create(out)
  per_res <- per_og_trees(dnmb, out, min_seqs = 2L, max_seqs = 50L,
                          method = "nj", threads = 1L, verbose = FALSE)
  bundle <- export_og_bundles(dnmb, per_res, out, tools = "hyphy")
  root <- bundle$hyphy

  fam_dirs <- list.dirs(root, recursive = FALSE)
  fam_dirs <- fam_dirs[grepl("/OG_\\d{7}$", fam_dirs)]
  expect_true(length(fam_dirs) > 0L)
  for (fd in fam_dirs) {
    expect_true(file.exists(file.path(fd, "aln.aa.fasta")))
    expect_true(file.exists(file.path(fd, "tree.nwk")))
  }

  readme <- readLines(file.path(root, "README.txt"), warn = FALSE)
  expect_true(any(grepl("hyphy fel",    readme)))
  expect_true(any(grepl("hyphy busted", readme)))
  expect_true(any(grepl("hyphy absrel", readme)))
  expect_true(any(grepl("back_translate_og_codons", readme)))
})


test_that(".export_mcscanx emits real contig IDs and FASTA-derived pairs", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  skip_if_no_aln_engine()

  dnmb <- make_synthetic_dnmb(n_genomes = 3L, n_ogs = 2L, seed = 13L)
  # Fake two contigs per genome so we can verify the exporter does not
  # collapse them onto a single "g<uid>" pseudo-chromosome.
  dnmb$id_map$contig <- paste0("contig",
                               (dnmb$id_map$protein_uid %% 2L) + 1L)
  out  <- tempfile("exp_mcs_"); dir.create(out)
  per_res <- per_og_trees(dnmb, out, min_seqs = 2L, max_seqs = 50L,
                          method = "nj", threads = 1L, verbose = FALSE)
  bundle <- export_og_bundles(dnmb, per_res, out, tools = "mcscanx")
  root <- bundle$mcscanx
  gff <- utils::read.table(file.path(root, "genomes.gff"), sep = "\t",
                           header = FALSE, stringsAsFactors = FALSE)
  # Column 1 is chrom (per-genome-tag + safe contig), and >= 2 distinct
  # chrom values per genome prefix would be impossible under the old
  # "g<uid>" collapse.
  chroms <- unique(gff$V1)
  expect_true(any(grepl("_contig1$", chroms)))
  expect_true(any(grepl("_contig2$", chroms)))

  # Blast file exists and every gene ID is of the p<uid>_g<gid> shape.
  blast <- utils::read.table(file.path(root, "genomes.blast"), sep = "\t",
                             header = FALSE, stringsAsFactors = FALSE)
  expect_true(all(grepl("^p\\d+_g\\d+$", blast$V1)))
  expect_true(all(grepl("^p\\d+_g\\d+$", blast$V2)))
})
