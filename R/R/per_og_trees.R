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
#' @param method Tree-building strategy: `"nj"` (distance NJ, fast),
#'   `"ml"` (phangorn `pml` + NNI, accurate), or `"iqtree"` (IQ-TREE
#'   `--fast` with UFBoot — fastest ML for large OGs; falls back to
#'   phangorn `"ml"` if the `iqtree2`/`iqtree` binary is not on PATH).
#' @param threads Thread count for alignment and tree engines.
#' @param verbose Logical; print a progress line per OG.
#' @param force If FALSE (default), OGs whose `aln.fasta` and `tree.nwk`
#'   already exist on disk are skipped (checkpointed resume). Pass TRUE
#'   to rebuild everything. Partial outputs (only the alignment, not the
#'   tree) are treated as incomplete and re-run.
#' @param workers If > 1 and the `future.apply` Suggests package is
#'   installed, OGs are processed in parallel using the user-supplied
#'   `future::plan()` (or a default multisession plan). If unset, a
#'   single-threaded loop is used. Mutually beneficial with
#'   `threads = 1L` inside each worker.
#' @param trim_gaps If TRUE (default), drop alignment columns with gap
#'   fraction exceeding `max_gap_frac` before tree building. Tracks
#'   OrthoFinder / standard phylogenomic practice (e.g. trimAl -gt 0.5).
#' @param max_gap_frac Column-gap fraction threshold for `trim_gaps`.
#' @param model Amino-acid substitution model for ML branch. Default
#'   `"LG"` (better than WAG/JTT for prokaryotic proteins).
#' @param gamma_k Number of Γ rate categories for ML (0 disables). Default 4.
#' @param bootstrap_reps Per-OG bootstrap replicates (0 disables). When > 0,
#'   `tree$node.label` is populated with integer support percentages. NJ
#'   uses `ape::boot.phylo`; ML uses `phangorn::bootstrap.pml`.
#' @return A tibble with columns
#'   `cluster_id, n_members, n_genomes, aln_path, tree_path, single_copy`.
#' @export
per_og_trees <- function(dnmb,
                         out_dir,
                         min_seqs = 3L,
                         max_seqs = 500L,
                         method = c("nj", "ml", "iqtree"),
                         threads = 2L,
                         verbose = TRUE,
                         force = FALSE,
                         workers = 1L,
                         trim_gaps = TRUE,
                         max_gap_frac = 0.5,
                         model = "LG",
                         gamma_k = 4L,
                         bootstrap_reps = 0L) {
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

  process_one <- function(df) {
    cid <- df$cluster_id[1]
    full_n_members <- nrow(df)
    full_n_genomes <- length(unique(df$genome_uid))
    is_sc <- full_n_members == n_total_genomes &&
      full_n_genomes == n_total_genomes

    if (!is.null(max_seqs) && nrow(df) > max_seqs) {
      df <- df[sample.int(nrow(df), max_seqs), , drop = FALSE]
    }

    og_dir <- file.path(og_root, sprintf("OG_%07d", cid))
    dir.create(og_dir, showWarnings = FALSE)
    aln_path <- file.path(og_dir, "aln.fasta")
    tree_path <- file.path(og_dir, "tree.nwk")

    cached <- !force &&
      file.exists(aln_path) && file.size(aln_path) > 0 &&
      file.exists(tree_path) && file.size(tree_path) > 0

    ok <- if (cached) TRUE else tryCatch({
      .align_and_tree(df, aln_path, tree_path, engine, method, threads,
                      trim_gaps = trim_gaps, max_gap_frac = max_gap_frac,
                      model = model, gamma_k = gamma_k,
                      bootstrap_reps = bootstrap_reps)
      TRUE
    }, error = function(e) {
      if (verbose) message(sprintf("  OG %d failed: %s", cid, conditionMessage(e)))
      FALSE
    })

    tibble::tibble(
      cluster_id  = cid,
      n_members   = full_n_members,
      n_genomes   = full_n_genomes,
      aln_path    = if (ok) aln_path else NA_character_,
      tree_path   = if (ok) tree_path else NA_character_,
      single_copy = ok && is_sc
    )
  }

  use_parallel <- workers > 1L && requireNamespace("future.apply", quietly = TRUE)
  if (workers > 1L && !use_parallel && verbose) {
    message("  workers>1 requested but 'future.apply' not installed; running serially.")
  }

  if (use_parallel) {
    if (!requireNamespace("future", quietly = TRUE)) {
      use_parallel <- FALSE
    } else {
      prev_plan <- future::plan()
      # Only set a default plan if the caller hasn't already chosen one.
      if (inherits(prev_plan, "sequential")) {
        future::plan(future::multisession, workers = workers)
        on.exit(future::plan(prev_plan), add = TRUE)
      }
    }
  }

  if (use_parallel) {
    rows <- future.apply::future_lapply(
      by_cluster, process_one,
      future.seed = TRUE,
      future.globals = TRUE,
      future.packages = c("dplyr", "tibble", "ape", "phangorn")
    )
  } else {
    rows <- vector("list", length(by_cluster))
    for (i in seq_along(by_cluster)) {
      rows[[i]] <- process_one(by_cluster[[i]])
      if (verbose && i %% 50 == 0) {
        message(sprintf("  [%d/%d] cluster_id=%d",
                        i, length(by_cluster), rows[[i]]$cluster_id))
      }
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
  if (method == "iqtree") {
    bin <- .iqtree_binary()
    if (nzchar(bin)) return("iqtree")
    if (requireNamespace("phangorn", quietly = TRUE)) {
      message("  iqtree2/iqtree not found on PATH; falling back to phangorn ML.")
      return("phangorn_ml")
    }
  }
  if (method == "ml" && requireNamespace("phangorn", quietly = TRUE)) {
    return("phangorn_ml")
  }
  if (requireNamespace("phangorn", quietly = TRUE)) return("phangorn_nj")
  if (requireNamespace("ape", quietly = TRUE)) return("ape_nj")
  stop("Need either 'phangorn' or 'ape' installed for tree building.")
}

.iqtree_binary <- function() {
  for (cand in c("iqtree2", "iqtree")) {
    p <- Sys.which(cand)
    if (nzchar(p)) return(unname(p))
  }
  ""
}

.align_and_tree <- function(df, aln_path, tree_path, engine, method, threads,
                            trim_gaps = TRUE, max_gap_frac = 0.5,
                            model = "LG", gamma_k = 4L,
                            bootstrap_reps = 0L) {
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

  if (isTRUE(trim_gaps)) {
    .trim_alignment(aln_path, max_gap_frac = max_gap_frac)
  }

  tree_eng <- .pick_tree_engine(method)
  if (tree_eng == "iqtree") {
    tree <- .iqtree_build(aln_path, threads = threads, model = model,
                          bootstrap_reps = bootstrap_reps)
  } else if (tree_eng == "phangorn_ml") {
    aln <- phangorn::read.phyDat(aln_path, format = "fasta", type = "AA")
    d <- phangorn::dist.ml(aln)
    start <- ape::nj(d)
    use_gamma <- isTRUE(gamma_k > 0L)
    fit <- if (use_gamma) {
      phangorn::pml(start, data = aln, model = model, k = as.integer(gamma_k))
    } else {
      phangorn::pml(start, data = aln, model = model)
    }
    fit <- phangorn::optim.pml(
      fit,
      model   = model,
      optNni  = TRUE,
      optGamma = use_gamma,
      control = phangorn::pml.control(trace = 0)
    )
    tree <- fit$tree
    if (bootstrap_reps > 0L) {
      bs <- phangorn::bootstrap.pml(fit, bs = bootstrap_reps,
                                     optNni = TRUE,
                                     control = phangorn::pml.control(trace = 0))
      supp <- ape::prop.clades(tree, bs)
      tree$node.label <- .format_support(supp, bootstrap_reps)
    }
  } else if (tree_eng == "phangorn_nj") {
    aln <- phangorn::read.phyDat(aln_path, format = "fasta", type = "AA")
    d <- phangorn::dist.ml(aln)
    tree <- ape::nj(d)
    if (bootstrap_reps > 0L) {
      tree <- .bootstrap_nj_phyDat(tree, aln, reps = bootstrap_reps)
    }
  } else {  # ape_nj
    aln <- ape::read.FASTA(aln_path, type = "AA")
    d <- ape::dist.aa(aln)
    tree <- ape::nj(d)
    if (bootstrap_reps > 0L) {
      tree <- .bootstrap_nj_ape(tree, aln, reps = bootstrap_reps)
    }
  }
  ape::write.tree(tree, tree_path)
}

# IQ-TREE --fast ML with optional UFBoot. Writes treefile to a temp prefix
# and returns the parsed ape::phylo. UFBoot percentages end up in node.label.
.iqtree_build <- function(aln_path, threads, model = "LG",
                          bootstrap_reps = 0L) {
  bin <- .iqtree_binary()
  if (!nzchar(bin)) stop("iqtree2/iqtree binary not on PATH.")
  tmp <- tempfile("iqtree_")
  on.exit(unlink(paste0(tmp, "*")), add = TRUE)
  iq_model <- if (nzchar(model)) paste0(model, "+G4") else "LG+G4"
  # IQ-TREE refuses UFBoot with <4 sequences ("It makes no sense to
  # perform bootstrap with less than 4 sequences."). Count *unique*
  # sequences — identical sequences get dedup'd internally by iqtree
  # before this check, so 5 header records with 2 unique seqs will
  # still trip the error (observed on geobacillus runs).
  n_seqs <- tryCatch(
    {
      lines <- readLines(aln_path, warn = FALSE)
      hdr_idx <- which(startsWith(lines, ">"))
      if (!length(hdr_idx)) {
        0L
      } else {
        seqs <- vapply(seq_along(hdr_idx), function(i) {
          start <- hdr_idx[i] + 1L
          end <- if (i < length(hdr_idx)) hdr_idx[i + 1L] - 1L else length(lines)
          if (start > end) "" else paste0(lines[start:end], collapse = "")
        }, character(1))
        length(unique(seqs[nzchar(seqs)]))
      }
    },
    error = function(e) NA_integer_
  )
  do_ufboot <- bootstrap_reps > 0L && !is.na(n_seqs) && n_seqs >= 4L
  # IQ-TREE constraint: --fast (FastTree-style hill climb) is incompatible
  # with -bb (UFBoot). Pick the combination the user actually asked for.
  args <- c("-s", shQuote(aln_path),
            "-m", iq_model,
            "-T", as.character(max(1L, as.integer(threads))),
            "--prefix", shQuote(tmp),
            "-redo", "-quiet")
  if (do_ufboot) {
    args <- c(args,
              "-bb", as.character(max(1000L, as.integer(bootstrap_reps))),
              "-alrt", "1000")
  } else {
    args <- c(args, "--fast")
  }
  log_path <- paste0(tmp, ".log")
  rc <- system2(bin, args = args, stdout = FALSE, stderr = log_path)
  # Fallback: if UFBoot was requested but iqtree rejected it (rare edge
  # case: gap-normalized seqs dedupe below 4 even though our raw-string
  # counter said otherwise), retry with --fast and no bootstrap rather
  # than bubbling the OG failure up.
  if (rc != 0 && do_ufboot) {
    err0 <- if (file.exists(log_path)) paste(readLines(log_path, warn = FALSE),
                                              collapse = "\n") else ""
    if (grepl("less than 4 sequences", err0, fixed = TRUE)) {
      args <- c("-s", shQuote(aln_path),
                "-m", iq_model,
                "-T", as.character(max(1L, as.integer(threads))),
                "--prefix", shQuote(tmp),
                "-redo", "-quiet", "--fast")
      rc <- system2(bin, args = args, stdout = FALSE, stderr = log_path)
      do_ufboot <- FALSE
    }
  }
  if (rc != 0) {
    err <- if (file.exists(log_path)) paste(readLines(log_path, warn = FALSE),
                                              collapse = "\n") else ""
    stop("iqtree failed (exit=", rc, "): ", substr(err, 1, 500))
  }
  tree_file <- paste0(tmp, ".treefile")
  if (!file.exists(tree_file))
    stop("iqtree produced no treefile: ", tree_file)
  tr <- ape::read.tree(tree_file)
  if (do_ufboot && !is.null(tr$node.label)) {
    tr$node.label <- sub("/.*$", "", tr$node.label)  # keep UFBoot, drop aLRT
  }
  tr
}

# Drop alignment columns whose gap fraction exceeds max_gap_frac.
# Operates in-place on a FASTA file. Preserves sequence order and names.
.trim_alignment <- function(path, max_gap_frac = 0.5) {
  lines <- readLines(path, warn = FALSE)
  if (!length(lines)) return(invisible())
  is_hdr <- startsWith(lines, ">")
  # Rebuild (name, seq) pairs — each record may span multiple lines.
  hdr_idx <- which(is_hdr)
  if (!length(hdr_idx)) return(invisible())
  seqs <- character(length(hdr_idx))
  for (i in seq_along(hdr_idx)) {
    start <- hdr_idx[i] + 1L
    end   <- if (i < length(hdr_idx)) hdr_idx[i + 1L] - 1L else length(lines)
    seqs[i] <- paste0(lines[start:end], collapse = "")
  }
  names(seqs) <- sub("^>", "", lines[hdr_idx])
  if (!length(seqs)) return(invisible())

  ncol <- unique(nchar(seqs))
  if (length(ncol) != 1L) return(invisible())  # not aligned; skip silently
  if (ncol == 0L) return(invisible())

  mat <- do.call(rbind, strsplit(seqs, "", fixed = TRUE))
  gap_frac <- colMeans(mat == "-" | mat == ".")
  keep <- gap_frac <= max_gap_frac
  if (!any(keep)) return(invisible())  # refuse to empty the alignment
  if (all(keep)) return(invisible())   # nothing to trim

  mat <- mat[, keep, drop = FALSE]
  trimmed <- apply(mat, 1L, paste0, collapse = "")

  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  for (i in seq_along(trimmed)) {
    cat(">", names(seqs)[i], "\n", trimmed[i], "\n",
        sep = "", file = con)
  }
  invisible()
}

.format_support <- function(supp, reps) {
  # prop.clades returns NA for unsupported or external nodes. Map NA→0.
  pct <- round(100 * ifelse(is.na(supp), 0, supp) / reps)
  as.character(pct)
}

.bootstrap_nj_phyDat <- function(tree, aln, reps) {
  # True per-site bootstrap: expand phyDat to a full-site char matrix,
  # resample columns with replacement, rebuild phyDat, redo NJ.
  # (Site-pattern resampling in phangorn::subset.phyDat is a different
  # statistic — sites with the same pattern get weight-multiplied, not
  # independently resampled — so we bypass it here to match the
  # semantics of ape::boot.phylo.)
  m <- as.character(aln)   # taxa × sites char matrix
  n_site <- ncol(m)
  bs_trees <- vector("list", reps)
  for (i in seq_len(reps)) {
    idx <- sample.int(n_site, replace = TRUE)
    aln_bs <- phangorn::phyDat(m[, idx, drop = FALSE], type = "AA")
    d <- phangorn::dist.ml(aln_bs)
    bs_trees[[i]] <- ape::nj(d)
  }
  supp <- ape::prop.clades(tree, bs_trees)
  tree$node.label <- .format_support(supp, reps)
  tree
}

.bootstrap_nj_ape <- function(tree, aln, reps) {
  m <- as.matrix(aln)  # DNAbin / AAbin
  bs_trees <- vector("list", reps)
  for (i in seq_len(reps)) {
    idx <- sample.int(ncol(m), replace = TRUE)
    sub <- m[, idx, drop = FALSE]
    class(sub) <- class(m)
    d <- ape::dist.aa(sub)
    bs_trees[[i]] <- ape::nj(d)
  }
  supp <- ape::prop.clades(tree, bs_trees)
  tree$node.label <- .format_support(supp, reps)
  tree
}

.write_fasta <- function(labels, seqs, path) {
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  for (i in seq_along(labels)) {
    cat(">", labels[i], "\n", seqs[i], "\n", sep = "", file = con)
  }
}
