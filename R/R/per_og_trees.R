#' Build per-orthogroup alignments and gene trees
#'
#' For every cluster with at least `min_seqs` members, writes an amino-acid
#' FASTA and a Newick tree to `out_dir/orthogroups/<cluster_id>/`.
#'
#' Engines are chosen in this order of preference:
#'   1. `DECIPHER::AlignSeqs` + `phangorn` ML/NJ   (pure R — no external bin)
#'   2. MAFFT + IQ-TREE via `system2`               (if binaries on PATH)
#'
#' Both paths produce identical output layout, so downstream consumers
#' (Phase B HOG cutter, Phase C reconciler, GeneRax export) do not need
#' to know which engine ran.
#'
#' @param dnmb Result of `load_dnmb()`. Must have been run on results
#'   produced with `--level protein` so `gene_table$translation` is set.
#' @param out_dir Directory under which `orthogroups/<cluster_id>/` is
#'   created. The directory is created if missing.
#' @param min_seqs Skip clusters with fewer than this many members.
#'   Default 3 (two-leaf trees carry no topological information).
#' @param max_seqs Cap large clusters by random subsampling. NULL disables.
#' @param method Tree-building strategy: `"nj"` (distance NJ, fast) or
#'   `"ml"` (phangorn `pml` + NNI, accurate).
#' @param threads Thread count for alignment and tree engines.
#' @param verbose Logical; print a progress line per OG.
#' @return A tibble with columns
#'   `cluster_id, n_members, n_genomes, aln_path, tree_path, single_copy`.
#' @export
per_og_trees <- function(dnmb,
                         out_dir,
                         min_seqs = 3L,
                         max_seqs = 500L,
                         method = c("nj", "ml"),
                         threads = 2L,
                         verbose = TRUE) {
  method <- match.arg(method)
  stopifnot(is.list(dnmb), "translation" %in% names(dnmb$gene_table))

  og_root <- file.path(out_dir, "orthogroups")
  dir.create(og_root, recursive = TRUE, showWarnings = FALSE)

  engine <- .pick_aln_engine()
  if (verbose) message(sprintf("[per_og_trees] engine=%s method=%s", engine, method))

  n_total_genomes <- nrow(dnmb$genome_meta)

  members <- dnmb$clusters %>%
    dplyr::select(cluster_id, protein_uid, genome_uid) %>%
    dplyr::inner_join(
      dnmb$gene_table %>% dplyr::select(protein_uid, translation),
      by = "protein_uid"
    )

  by_cluster <- split(members, members$cluster_id)
  keep_mask <- vapply(by_cluster, function(df) nrow(df) >= min_seqs, logical(1))
  by_cluster <- by_cluster[keep_mask]

  if (verbose) {
    message(sprintf("[per_og_trees] building %d / %d orthogroups (min_seqs=%d)",
                    length(by_cluster), length(keep_mask), min_seqs))
  }

  rows <- vector("list", length(by_cluster))
  for (i in seq_along(by_cluster)) {
    df <- by_cluster[[i]]
    cid <- df$cluster_id[1]

    if (!is.null(max_seqs) && nrow(df) > max_seqs) {
      df <- df[sample.int(nrow(df), max_seqs), , drop = FALSE]
    }

    og_dir <- file.path(og_root, sprintf("OG_%07d", cid))
    dir.create(og_dir, showWarnings = FALSE)
    aln_path <- file.path(og_dir, "aln.fasta")
    tree_path <- file.path(og_dir, "tree.nwk")

    ok <- tryCatch({
      .align_and_tree(df, aln_path, tree_path, engine, method, threads)
      TRUE
    }, error = function(e) {
      if (verbose) message(sprintf("  OG %d failed: %s", cid, conditionMessage(e)))
      FALSE
    })

    n_genomes_in_og <- length(unique(df$genome_uid))
    single_copy <- ok && n_genomes_in_og == n_total_genomes &&
      nrow(df) == n_total_genomes

    rows[[i]] <- tibble::tibble(
      cluster_id  = cid,
      n_members   = nrow(df),
      n_genomes   = n_genomes_in_og,
      aln_path    = if (ok) aln_path else NA_character_,
      tree_path   = if (ok) tree_path else NA_character_,
      single_copy = single_copy
    )

    if (verbose && i %% 50 == 0) {
      message(sprintf("  [%d/%d] OG_%07d n=%d", i, length(by_cluster), cid, nrow(df)))
    }
  }

  result <- dplyr::bind_rows(rows)

  summary_path <- file.path(out_dir, "orthogroup_trees.tsv")
  utils::write.table(result, summary_path, sep = "\t",
                     quote = FALSE, row.names = FALSE)

  sc <- result[result$single_copy & !is.na(result$aln_path), ]
  sc_path <- file.path(out_dir, "single_copy_OGs.tsv")
  utils::write.table(sc[, c("cluster_id", "n_members", "aln_path", "tree_path")],
                     sc_path, sep = "\t", quote = FALSE, row.names = FALSE)

  if (verbose) {
    message(sprintf("[per_og_trees] done. built=%d single_copy=%d summary=%s",
                    sum(!is.na(result$aln_path)), sum(result$single_copy),
                    summary_path))
  }
  result
}


#' Identify single-copy orthogroups
#'
#' An OG is single-copy if it is present in exactly one copy per genome
#' across all `n_total` genomes. These are the rows that can be
#' concatenated into a phylogenomic supermatrix.
#'
#' @param dnmb Result of `load_dnmb()`.
#' @return Tibble: `cluster_id, n_members, n_genomes`.
#' @export
single_copy_ogs <- function(dnmb) {
  n_total <- nrow(dnmb$genome_meta)
  dnmb$clusters %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(
      n_members = dplyr::n(),
      n_genomes = dplyr::n_distinct(genome_uid),
      .groups   = "drop"
    ) %>%
    dplyr::filter(n_members == n_total, n_genomes == n_total)
}


# ---- internal helpers ------------------------------------------------------

.pick_aln_engine <- function() {
  has_decipher <- requireNamespace("DECIPHER", quietly = TRUE) &&
    requireNamespace("Biostrings", quietly = TRUE)
  if (has_decipher) return("decipher")
  mafft <- Sys.which("mafft")
  if (nzchar(mafft)) return("mafft")
  stop("No alignment engine available. Install Bioconductor 'DECIPHER' ",
       "(and 'Biostrings'), or put MAFFT on PATH.")
}

.pick_tree_engine <- function(method) {
  if (method == "ml" && requireNamespace("phangorn", quietly = TRUE)) {
    return("phangorn_ml")
  }
  if (requireNamespace("phangorn", quietly = TRUE)) return("phangorn_nj")
  if (requireNamespace("ape", quietly = TRUE)) return("ape_nj")
  stop("Need either 'phangorn' or 'ape' installed for tree building.")
}

.align_and_tree <- function(df, aln_path, tree_path, engine, method, threads) {
  # protein_uid is uint64 → rendered via bit64/integer64 in R. `%d`
  # chokes on those; paste0 goes through format.integer64 cleanly.
  labels <- paste0("p", df$protein_uid, "_g", df$genome_uid)
  seqs   <- df$translation

  if (engine == "decipher") {
    aa <- Biostrings::AAStringSet(seqs)
    names(aa) <- labels
    aligned <- DECIPHER::AlignSeqs(aa, verbose = FALSE, processors = threads)
    Biostrings::writeXStringSet(aligned, aln_path)
  } else if (engine == "mafft") {
    in_tmp <- tempfile(fileext = ".fa")
    on.exit(unlink(in_tmp), add = TRUE)
    .write_fasta(labels, seqs, in_tmp)
    res <- system2("mafft",
                   args = c("--auto", "--thread", threads, shQuote(in_tmp)),
                   stdout = aln_path, stderr = FALSE)
    if (res != 0) stop("mafft failed (exit=", res, ")")
  }

  tree_eng <- .pick_tree_engine(method)
  if (tree_eng == "phangorn_ml") {
    aln <- phangorn::read.phyDat(aln_path, format = "fasta", type = "AA")
    d <- phangorn::dist.ml(aln)
    start <- ape::nj(d)
    fit <- phangorn::pml(start, data = aln)
    fit <- phangorn::optim.pml(fit, optNni = TRUE,
                               control = phangorn::pml.control(trace = 0))
    tree <- fit$tree
  } else if (tree_eng == "phangorn_nj") {
    aln <- phangorn::read.phyDat(aln_path, format = "fasta", type = "AA")
    d <- phangorn::dist.ml(aln)
    tree <- ape::nj(d)
  } else {  # ape_nj
    aln <- ape::read.FASTA(aln_path, type = "AA")
    d <- ape::dist.aa(aln)
    tree <- ape::nj(d)
  }
  ape::write.tree(tree, tree_path)
}

.write_fasta <- function(labels, seqs, path) {
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  for (i in seq_along(labels)) {
    cat(">", labels[i], "\n", seqs[i], "\n", sep = "", file = con)
  }
}
