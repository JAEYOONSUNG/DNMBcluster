#' Incrementally add a new genome (or proteome) to a prior DNMB run
#'
#' Implements the SHOOT-style add-species workflow in pure R + optional
#' DIAMOND. Instead of recomputing the whole pan-genome, only the
#' orthogroups touched by the new proteins are re-aligned and re-treed.
#' Unaffected OGs and the rooted species tree are preserved.
#'
#' Steps:
#' \enumerate{
#'   \item Load base run via `load_dnmb(base_dir)`.
#'   \item Extract per-OG representative proteins (is_centroid == TRUE)
#'         into a single FASTA.
#'   \item Score new proteins against representatives -- DIAMOND if on PATH,
#'         else Biostrings::pairwiseAlignment (slow but pure R).
#'   \item Assign each new protein to the best-hit OG above
#'         `min_bitscore`, else flag as orphan.
#'   \item Re-align and re-tree only the OGs that received >= 1 new member.
#'   \item Write `incremental/` with assignments, touched-OG trees,
#'         orphan FASTA, and a summary TSV.
#' }
#'
#' Grafting new tips onto the rooted species tree is out of scope for
#' this MVP -- the user runs `stag_species_tree()` + `stride_root()`
#' again on the updated per_og_trees output to refresh the species tree.
#'
#' @param base_dir Path to a prior DNMB run directory (same input as
#'   `load_dnmb()`), which also contains the `per_og_trees()` output
#'   under `<base_dir>/orthogroups/` and `orthogroup_trees.tsv`.
#' @param new_fasta Path to a FASTA of amino-acid sequences to add.
#'   IDs may be anything; each new sequence is assigned a synthetic
#'   `protein_uid` and `genome_uid` so its gene-tree leaf label
#'   (`p<protein_uid>_g<genome_uid>`) is parseable by the same
#'   regex (`sub(".*_g", "", lab)`) used across the package.
#' @param out_dir Output directory. Creates `<out_dir>/incremental/`.
#' @param min_bitscore Minimum DIAMOND (or pairwise) bitscore for an
#'   OG assignment. Below this the protein is considered orphan.
#' @param new_genome_key Genome label for the new proteome (single string).
#'   Recorded in the assignment table and written to `new_genome.tsv`
#'   alongside the synthetic `new_genome_uid`.
#' @param new_genome_uid Optional integer genome_uid for the new proteome.
#'   Defaults to `max(existing genome_uid) + 1`, guaranteeing it does not
#'   collide with any genome already registered in `dnmb$genome_meta`.
#' @param threads Thread count passed to DIAMOND / alignment.
#' @param method Tree-building method for touched OGs: `"nj"` or `"ml"`.
#' @return Tibble summary of touched OGs with re-alignment and re-tree paths.
#' @export
incremental_add <- function(base_dir,
                            new_fasta,
                            out_dir,
                            min_bitscore = 50,
                            new_genome_key = "NEW",
                            new_genome_uid = NULL,
                            threads = 2L,
                            method = c("nj", "ml")) {
  method <- match.arg(method)

  dnmb <- load_dnmb(base_dir)
  base_og_tsv <- file.path(base_dir, "orthogroup_trees.tsv")
  if (!file.exists(base_og_tsv)) {
    stop("base_dir is missing orthogroup_trees.tsv; run per_og_trees() first.")
  }
  base_ogs <- utils::read.table(base_og_tsv, sep = "\t", header = TRUE,
                                 stringsAsFactors = FALSE)

  if (is.null(new_genome_uid)) {
    new_genome_uid <- as.integer(max(dnmb$genome_meta$genome_uid, 0L)) + 1L
  }
  if (new_genome_uid %in% dnmb$genome_meta$genome_uid) {
    stop("new_genome_uid ", new_genome_uid,
         " already exists in base run genome_meta.")
  }

  out <- file.path(out_dir, "incremental")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  # 1. Build representatives FASTA (one centroid per OG) ----------------
  reps_path <- file.path(out, "representatives.fasta")
  .write_representatives(dnmb, reps_path)

  # 2. Score new proteins vs reps ---------------------------------------
  new_aa <- .read_fasta_aa(new_fasta)
  hits <- .score_new_vs_reps(new_fasta, reps_path, dnmb, threads = threads)

  # Mint a synthetic protein_uid for every input sequence so downstream
  # labels match the package-wide `p<protein_uid>_g<genome_uid>` convention.
  # protein_uid may be integer64 (bit64::integer64); stay in that domain
  # so large DNMB runs do not silently truncate at .Machine$integer.max.
  base_puid <- if (requireNamespace("bit64", quietly = TRUE)) {
    as_int64 <- tryCatch(bit64::as.integer64(0L), error = function(e) NULL)
    if (!is.null(as_int64)) bit64::as.integer64(max(dnmb$gene_table$protein_uid, 0))
    else max(as.numeric(dnmb$gene_table$protein_uid), 0)
  } else {
    max(as.numeric(dnmb$gene_table$protein_uid), 0)
  }
  new_puids <- base_puid + seq_along(new_aa)
  id_book <- tibble::tibble(
    new_id       = names(new_aa),
    protein_uid  = new_puids,
    genome_uid   = new_genome_uid
  )

  # 3. Assign each new protein -----------------------------------------
  assign <- .assign_best_hit(hits, min_bitscore)
  # DIAMOND drops queries with no hits; add them back as orphans so the
  # assignment accounting matches the input FASTA row count.
  missing_ids <- setdiff(names(new_aa), assign$new_id)
  if (length(missing_ids)) {
    assign <- dplyr::bind_rows(
      assign,
      tibble::tibble(new_id = missing_ids,
                     cluster_id = NA_integer_,
                     bitscore = NA_real_)
    )
  }
  assign <- dplyr::left_join(assign, id_book, by = "new_id")
  assign_path <- file.path(out, "assignments.tsv")
  utils::write.table(assign, assign_path, sep = "\t",
                     quote = FALSE, row.names = FALSE)

  # Register the synthetic genome so the user can join back to the label scheme.
  utils::write.table(
    tibble::tibble(genome_uid = new_genome_uid,
                   genome_key = new_genome_key,
                   source     = normalizePath(new_fasta, mustWork = FALSE)),
    file.path(out, "new_genome.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # 4. Re-align / re-tree touched OGs ----------------------------------
  touched <- assign[!is.na(assign$cluster_id), ]
  touched_cids <- unique(touched$cluster_id)
  message(sprintf("[incremental] %d new proteins assigned across %d OGs (%d orphan)",
                  nrow(touched), length(touched_cids),
                  sum(is.na(assign$cluster_id))))

  rows <- list()
  for (cid in touched_cids) {
    base_row <- base_ogs[base_ogs$cluster_id == cid, ]
    if (!nrow(base_row) || is.na(base_row$aln_path)) next
    new_ids <- touched$new_id[touched$cluster_id == cid]
    new_meta <- id_book[match(new_ids, id_book$new_id), ]

    dest <- file.path(out, "orthogroups", sprintf("OG_%07d", cid))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    aln_dest <- file.path(dest, "aln.fasta")
    tree_dest <- file.path(dest, "tree.nwk")

    ok <- tryCatch({
      .realign_with_new(base_row$aln_path, new_aa, new_ids,
                        new_meta$protein_uid, new_meta$genome_uid,
                        aln_dest, tree_dest, method, threads)
      TRUE
    }, error = function(e) { message("  OG ", cid, " failed: ", conditionMessage(e)); FALSE })

    rows[[length(rows) + 1]] <- tibble::tibble(
      cluster_id    = cid,
      n_new_members = length(new_ids),
      aln_path      = if (ok) aln_dest else NA_character_,
      tree_path     = if (ok) tree_dest else NA_character_
    )
  }
  result <- dplyr::bind_rows(rows)

  # 5. Orphan FASTA -----------------------------------------------------
  orphan_ids <- assign$new_id[is.na(assign$cluster_id)]
  if (length(orphan_ids)) {
    orphan_path <- file.path(out, "orphans.fasta")
    .write_subset_fasta(new_aa, orphan_ids, orphan_path)
    message(sprintf("[incremental] wrote %d orphan sequences to %s",
                    length(orphan_ids), orphan_path))
  }

  summary_path <- file.path(out, "touched_OGs.tsv")
  utils::write.table(result, summary_path, sep = "\t",
                     quote = FALSE, row.names = FALSE)
  message(sprintf("[incremental] summary -> %s", summary_path))
  result
}


# ---- helpers ---------------------------------------------------------------

.write_representatives <- function(dnmb, out_path) {
  reps <- dnmb$clusters[dnmb$clusters$is_centroid, ]
  reps <- reps[!duplicated(reps$cluster_id), ]
  reps <- dplyr::inner_join(
    reps[, c("cluster_id", "protein_uid", "genome_uid")],
    dnmb$gene_table[, c("protein_uid", "translation")],
    by = "protein_uid"
  )
  con <- file(out_path, "w"); on.exit(close(con), add = TRUE)
  for (i in seq_len(nrow(reps))) {
    cat(">OG_", sprintf("%07d", reps$cluster_id[i]),
        "|p", format(reps$protein_uid[i], scientific = FALSE),
        "\n", reps$translation[i], "\n", sep = "", file = con)
  }
}

.read_fasta_aa <- function(path) {
  if (requireNamespace("Biostrings", quietly = TRUE)) {
    aa <- Biostrings::readAAStringSet(path)
    return(stats::setNames(as.character(aa), names(aa)))
  }
  lines <- readLines(path)
  headers <- grep("^>", lines)
  out <- character(length(headers))
  nm  <- character(length(headers))
  ends <- c(utils::tail(headers - 1, -1), length(lines))
  for (i in seq_along(headers)) {
    nm[i] <- sub("^>\\s*(\\S+).*", "\\1", lines[headers[i]])
    out[i] <- paste(lines[(headers[i] + 1):ends[i]], collapse = "")
  }
  stats::setNames(out, nm)
}

.score_new_vs_reps <- function(new_fa, reps_fa, dnmb, threads) {
  diamond <- Sys.which("diamond")
  if (nzchar(diamond)) {
    db <- tempfile(fileext = ".dmnd")
    out <- tempfile(fileext = ".m8")
    on.exit({ unlink(db); unlink(out) }, add = TRUE)
    res <- system2("diamond",
                   c("makedb", "--in", shQuote(reps_fa), "-d", shQuote(sub("\\.dmnd$", "", db))),
                   stdout = FALSE, stderr = FALSE)
    if (res != 0) stop("diamond makedb failed (exit=", res, ")")
    res <- system2("diamond",
                   c("blastp", "-q", shQuote(new_fa), "-d", shQuote(sub("\\.dmnd$", "", db)),
                     "-o", shQuote(out), "-k", "5", "--threads", threads,
                     "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                     "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                     "evalue", "bitscore"),
                   stdout = FALSE, stderr = FALSE)
    if (res != 0) stop("diamond blastp failed (exit=", res, ")")
    df <- utils::read.table(out, sep = "\t", header = FALSE,
                            stringsAsFactors = FALSE,
                            col.names = c("qseqid","sseqid","pident","length",
                                           "mismatch","gapopen","qstart","qend",
                                           "sstart","send","evalue","bitscore"))
    return(df)
  }
  # Pure-R fallback: Biostrings pairwise vs every representative
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Neither DIAMOND nor Biostrings available -- cannot score new vs reps.")
  }
  message("[incremental] DIAMOND not on PATH -- using Biostrings pairwise (slow)")
  reps <- .read_fasta_aa(reps_fa)
  new_aa <- .read_fasta_aa(new_fa)
  reps_set <- Biostrings::AAStringSet(reps)
  new_set  <- Biostrings::AAStringSet(new_aa)
  rows <- list()
  n_reps <- length(reps_set)
  for (q in names(new_set)) {
    # Biostrings requires length(pattern) == length(subject) (or both == 1).
    # Replicate the single query to match reps_set so we get one score per rep.
    q_rep <- new_set[q][rep(1L, n_reps)]
    scores <- Biostrings::pairwiseAlignment(
      pattern = q_rep, subject = reps_set,
      substitutionMatrix = "BLOSUM62",
      scoreOnly = TRUE, type = "local")
    k <- min(5, length(scores))
    top <- order(scores, decreasing = TRUE)[seq_len(k)]
    for (idx in top) {
      rows[[length(rows) + 1]] <- data.frame(
        qseqid = q, sseqid = names(reps_set)[idx],
        bitscore = scores[idx], evalue = NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

.assign_best_hit <- function(hits, min_bitscore) {
  # One row per query with best-bitscore OG, or NA if below cutoff
  if (!nrow(hits)) return(tibble::tibble(new_id = character(), cluster_id = integer(), bitscore = numeric()))
  best <- hits %>%
    dplyr::group_by(qseqid) %>%
    dplyr::slice_max(bitscore, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  cid <- sub("^OG_0*", "", sub("\\|.*$", "", best$sseqid))
  cid <- suppressWarnings(as.integer(cid))
  keep <- best$bitscore >= min_bitscore & !is.na(cid)
  tibble::tibble(
    new_id     = best$qseqid,
    cluster_id = ifelse(keep, cid, NA_integer_),
    bitscore   = best$bitscore
  )
}

.realign_with_new <- function(prior_aln, new_aa, new_ids,
                              new_protein_uids, new_genome_uids,
                              aln_dest, tree_dest, method, threads) {
  prior <- .read_fasta_aa(prior_aln)
  addition <- new_aa[new_ids]
  # Tree-leaf labels must match the `p<protein_uid>_g<genome_uid>` scheme
  # every downstream parser uses (regex: sub(".*_g", "", lab)).
  names(addition) <- paste0("p", new_protein_uids, "_g", new_genome_uids)
  # Concatenate unaligned prior + new; re-align the whole thing. For
  # very large alignments a MAFFT --add path would be faster, but
  # keeping a single code path is cleaner.
  prior_stripped <- gsub("-", "", prior)
  combined <- c(prior_stripped, addition)

  tmp_in <- tempfile(fileext = ".fa")
  on.exit(unlink(tmp_in), add = TRUE)
  con <- file(tmp_in, "w")
  for (nm in names(combined)) {
    cat(">", nm, "\n", combined[nm], "\n", sep = "", file = con)
  }
  close(con)

  if (requireNamespace("DECIPHER", quietly = TRUE) &&
      requireNamespace("Biostrings", quietly = TRUE)) {
    aa <- Biostrings::AAStringSet(combined)
    aligned <- DECIPHER::AlignSeqs(aa, verbose = FALSE, processors = threads)
    Biostrings::writeXStringSet(aligned, aln_dest)
  } else if (nzchar(Sys.which("mafft"))) {
    res <- system2("mafft",
                   c("--auto", "--thread", threads, shQuote(tmp_in)),
                   stdout = aln_dest, stderr = FALSE)
    if (res != 0) stop("mafft failed (exit=", res, ")")
  } else {
    stop("Need DECIPHER or MAFFT to realign touched OG.")
  }

  # Tree
  if (method == "ml" && requireNamespace("phangorn", quietly = TRUE)) {
    aln <- phangorn::read.phyDat(aln_dest, format = "fasta", type = "AA")
    d <- phangorn::dist.ml(aln)
    start <- ape::nj(d)
    fit <- phangorn::pml(start, data = aln)
    fit <- phangorn::optim.pml(fit, optNni = TRUE,
                               control = phangorn::pml.control(trace = 0))
    tree <- fit$tree
  } else if (requireNamespace("phangorn", quietly = TRUE)) {
    aln <- phangorn::read.phyDat(aln_dest, format = "fasta", type = "AA")
    tree <- ape::nj(phangorn::dist.ml(aln))
  } else {
    aln <- ape::read.FASTA(aln_dest, type = "AA")
    tree <- ape::nj(ape::dist.aa(aln))
  }
  ape::write.tree(tree, tree_dest)
}

.write_subset_fasta <- function(seq_vec, ids, path) {
  con <- file(path, "w"); on.exit(close(con), add = TRUE)
  for (id in ids) {
    cat(">", id, "\n", seq_vec[[id]], "\n", sep = "", file = con)
  }
}
