#' Run Foldseek structural orthology and score per-HOG support
#'
#' Foldseek is the modern linear-time structural aligner: it encodes
#' 3-D tertiary contacts into a 3Di alphabet and runs MMseqs2-style
#' k-mer prefilter + gapped SW. That makes it 10^4–10^6× faster than
#' classical DALI/TM-align pairwise searches, which is why it is now
#' the default choice for whole-proteome structural similarity.
#'
#' This wrapper supports two modes:
#' \itemize{
#'   \item \strong{pdb}: user supplies a directory of predicted or
#'         experimental PDB/mmCIF files (e.g. AlphaFold DB, ESM-Fold).
#'         We call \code{foldseek easy-search} against itself.
#'   \item \strong{prosT5}: no structures available. Foldseek's
#'         \code{--prostt5-model} option lets it predict 3Di from AA
#'         sequence alone, so we can get structural orthology at
#'         sequence-search cost. Requires the ProstT5 weights (see
#'         \code{prosT5_weights}).
#' }
#'
#' Per-HOG support scores are computed over \emph{unordered} within-HOG
#' pairs: the denominator is \code{choose(n_members, 2)} for each HOG,
#' and a pair counts as supported if either direction's Foldseek hit
#' has \code{prob >= 0.9}. Missing pairs (no Foldseek hit) count as
#' unsupported rather than being excluded, which is the honest measure
#' of structural corroboration — a HOG that Foldseek never found hits
#' for should score low, not be dropped from the denominator.
#'
#' @param dnmb Result of `load_dnmb()`.
#' @param hog_table Result of `hog_cut()` (must have `hog_id, protein_uid`).
#' @param input_dir For `mode = "pdb"`, directory of per-protein PDBs
#'   named `p<protein_uid>.pdb`. For `mode = "prosT5"`, path to a
#'   combined AA FASTA with headers `p<protein_uid>`.
#' @param out_dir Output directory; created if missing.
#' @param mode `"pdb"` (default, structures on disk) or `"prosT5"`
#'   (AA-only, 3Di predicted).
#' @param prosT5_weights Path to ProstT5 weights directory. Required
#'   when \code{mode = "prosT5"}.
#' @param sensitivity Foldseek \code{-s} parameter. Default 9.5 —
#'   the high-sensitivity preset (>7.5 is recommended for remote
#'   homologs). Lower values (4–6) run noticeably faster.
#' @param threads Threads for Foldseek.
#' @param prob_cutoff Posterior probability cut-off for counting a
#'   pair as structurally supported. Default 0.9.
#' @param binary Explicit path to `foldseek`. NULL → `Sys.which`.
#' @return Named list: `per_hog` (`hog_id, n_pairs, n_supported,
#'   support_frac, mean_tmscore`), `pairs` (full Foldseek alignment
#'   table tibble), `out_dir`. NULL with a warning when the binary is
#'   missing.
#' @export
run_foldseek_hog <- function(dnmb,
                             hog_table,
                             input_dir,
                             out_dir,
                             mode           = c("pdb", "prosT5"),
                             prosT5_weights = NULL,
                             sensitivity    = 9.5,
                             threads        = 4L,
                             prob_cutoff    = 0.9,
                             binary         = NULL) {
  mode <- match.arg(mode)
  bin <- if (!is.null(binary) && nzchar(binary)) binary else unname(Sys.which("foldseek"))
  if (!nzchar(bin)) {
    warning("[run_foldseek_hog] 'foldseek' binary not found on PATH; skipping.")
    return(NULL)
  }
  stopifnot("hog_id" %in% names(hog_table), "protein_uid" %in% names(hog_table))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  hits_tsv <- file.path(out_dir, "foldseek_aln.tsv")
  tmp_dir  <- file.path(out_dir, "tmp"); dir.create(tmp_dir, showWarnings = FALSE)

  # Filter `input_dir` down to HOG-member proteins before searching.
  # Without this, foldseek does full-proteome all-vs-all over every file
  # under input_dir — O(N^2) in total protein count — but we only need
  # hits between proteins that share a HOG. Staging costs one hardlink
  # (or FASTA filter) per member and cuts runtime by orders of magnitude
  # on real pangenomes.
  member_uids <- unique(as.integer(hog_table$protein_uid))
  scoped_dir <- file.path(out_dir, "search_input")
  unlink(scoped_dir, recursive = TRUE)
  dir.create(scoped_dir, showWarnings = FALSE, recursive = TRUE)
  scoped_input <- .foldseek_scope_input(input_dir, scoped_dir,
                                        member_uids, mode)

  fmt <- "query,target,prob,fident,alnlen,evalue,bits,alntmscore"
  args <- c("easy-search", shQuote(scoped_input), shQuote(scoped_input),
            shQuote(hits_tsv), shQuote(tmp_dir),
            "--format-mode", "4",
            "--format-output", fmt,
            "-s", as.character(sensitivity),
            "--threads", as.character(threads),
            "-e", "0.001")
  if (mode == "prosT5") {
    if (is.null(prosT5_weights) || !nzchar(prosT5_weights))
      stop("mode='prosT5' requires prosT5_weights = <path to weights>")
    args <- c(args, "--prostt5-model", shQuote(prosT5_weights))
  }
  rc <- system2(bin, args = args, stdout = FALSE, stderr = FALSE)
  if (rc != 0) stop("foldseek failed (exit=", rc, ")")
  if (!file.exists(hits_tsv)) stop("foldseek produced no hits: ", hits_tsv)

  score_foldseek_hog(hits_tsv, hog_table, prob_cutoff)
}


# Stage the HOG-member subset of `input_dir` under `scoped_dir`. In pdb
# mode we link/copy each p<uid>.pdb/.cif/.ent that exists; in prosT5
# mode we filter the concatenated AA FASTA to the headers we need.
# Returns the path foldseek should search over (a directory for pdb
# mode, a file for prosT5 mode).
.foldseek_scope_input <- function(input_dir, scoped_dir, member_uids, mode) {
  if (mode == "pdb") {
    exts <- c(".pdb", ".pdb.gz", ".cif", ".cif.gz", ".ent", ".ent.gz")
    for (uid in member_uids) {
      base <- paste0("p", uid)
      for (ext in exts) {
        src <- file.path(input_dir, paste0(base, ext))
        if (file.exists(src)) {
          dst <- file.path(scoped_dir, paste0(base, ext))
          ok <- suppressWarnings(file.link(src, dst))
          if (!isTRUE(ok)) file.copy(src, dst, overwrite = TRUE)
          break
        }
      }
    }
    scoped_dir
  } else {
    # prosT5: input_dir is a combined AA FASTA with p<uid> headers.
    keep_set <- paste0("p", as.character(member_uids))
    out_fa <- file.path(scoped_dir, "members.faa")
    con_in  <- file(input_dir, "r"); on.exit(close(con_in), add = TRUE)
    con_out <- file(out_fa, "w");    on.exit(close(con_out), add = TRUE)
    keep <- FALSE
    while (length(ln <- readLines(con_in, n = 1L, warn = FALSE))) {
      if (startsWith(ln, ">")) {
        hdr <- sub("^>", "", ln)
        # header may carry extra tokens after the id; accept first token.
        first <- strsplit(hdr, "\\s+", perl = TRUE)[[1]][1]
        keep <- first %in% keep_set
      }
      if (keep) writeLines(ln, con_out)
    }
    out_fa
  }
}


#' Score per-HOG structural support from a Foldseek alignment table
#'
#' Separated so pre-computed Foldseek outputs (e.g. from a cluster
#' run) can be fed in without re-running the binary.
#'
#' @param hits_tsv Path to Foldseek alignment TSV with columns
#'   `query, target, prob, fident, alnlen, evalue, bits, alntmscore`
#'   (format-mode=4 header preserved).
#' @param hog_table Tibble with `hog_id`, `protein_uid`.
#' @param prob_cutoff Posterior probability threshold.
#' @return Named list: `per_hog`, `pairs`, `hits_tsv`.
#' @export
score_foldseek_hog <- function(hits_tsv, hog_table, prob_cutoff = 0.9) {
  stopifnot(file.exists(hits_tsv))
  pairs <- utils::read.table(hits_tsv, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, check.names = FALSE)
  pairs <- tibble::as_tibble(pairs)

  # Expected denominator: choose(n_members, 2) unordered pairs per HOG,
  # computed from hog_table — independent of what foldseek returned so a
  # HOG with zero recovered pairs still scores support_frac = 0 instead
  # of being dropped.
  ht <- tibble::as_tibble(hog_table) %>%
    dplyr::distinct(hog_id = as.character(.data$hog_id),
                    protein_uid = as.integer(.data$protein_uid))
  hog_sizes <- ht %>%
    dplyr::group_by(.data$hog_id) %>%
    dplyr::summarise(n_members = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(n_expected_pairs = as.integer(.data$n_members *
                                                  (.data$n_members - 1L) / 2L))

  if (!nrow(pairs)) {
    per_hog <- hog_sizes %>%
      dplyr::transmute(hog_id = .data$hog_id,
                       n_pairs = 0L,
                       n_supported = 0L,
                       n_expected_pairs = .data$n_expected_pairs,
                       support_frac = 0,
                       mean_tmscore = NA_real_)
    return(list(per_hog = per_hog, pairs = pairs, hits_tsv = hits_tsv))
  }
  pairs$q_uid <- suppressWarnings(as.integer(sub("^p", "",
                                                  sub("\\..*$", "", pairs$query))))
  pairs$t_uid <- suppressWarnings(as.integer(sub("^p", "",
                                                  sub("\\..*$", "", pairs$target))))
  uid2hog <- stats::setNames(ht$hog_id, as.character(ht$protein_uid))
  pairs$q_hog <- unname(uid2hog[as.character(pairs$q_uid)])
  pairs$t_hog <- unname(uid2hog[as.character(pairs$t_uid)])
  intra <- pairs[!is.na(pairs$q_hog) &
                   !is.na(pairs$t_hog) &
                   pairs$q_hog == pairs$t_hog &
                   pairs$q_uid != pairs$t_uid, , drop = FALSE]

  if (!nrow(intra)) {
    per_hog <- hog_sizes %>%
      dplyr::transmute(hog_id = .data$hog_id,
                       n_pairs = 0L,
                       n_supported = 0L,
                       n_expected_pairs = .data$n_expected_pairs,
                       support_frac = 0,
                       mean_tmscore = NA_real_)
    return(list(per_hog = per_hog, pairs = pairs, hits_tsv = hits_tsv))
  }

  # Canonicalize to unordered pairs so reciprocal (A->B, B->A) hits
  # collapse to one row. Supported = either direction clears prob_cutoff.
  intra$pair_key <- ifelse(intra$q_uid < intra$t_uid,
                           paste(intra$q_uid, intra$t_uid, sep = "_"),
                           paste(intra$t_uid, intra$q_uid, sep = "_"))
  intra$alntmscore_num <- suppressWarnings(as.numeric(intra$alntmscore))
  observed <- intra %>%
    dplyr::group_by(.data$q_hog, .data$pair_key) %>%
    dplyr::summarise(
      best_prob  = max(.data$prob, na.rm = TRUE),
      best_tm    = max(.data$alntmscore_num, na.rm = TRUE),
      .groups    = "drop"
    )
  observed$best_tm[!is.finite(observed$best_tm)] <- NA_real_

  observed_per_hog <- observed %>%
    dplyr::group_by(hog_id = .data$q_hog) %>%
    dplyr::summarise(
      n_pairs      = dplyr::n(),
      n_supported  = sum(.data$best_prob >= prob_cutoff, na.rm = TRUE),
      mean_tmscore = mean(.data$best_tm, na.rm = TRUE),
      .groups      = "drop"
    )

  per_hog <- hog_sizes %>%
    dplyr::left_join(observed_per_hog, by = "hog_id") %>%
    dplyr::mutate(
      n_pairs      = dplyr::coalesce(as.integer(.data$n_pairs), 0L),
      n_supported  = dplyr::coalesce(as.integer(.data$n_supported), 0L),
      support_frac = ifelse(.data$n_expected_pairs > 0L,
                            .data$n_supported / .data$n_expected_pairs,
                            NA_real_),
      mean_tmscore = .data$mean_tmscore
    ) %>%
    dplyr::select("hog_id", "n_pairs", "n_supported",
                  "n_expected_pairs", "support_frac", "mean_tmscore")
  list(per_hog = per_hog, pairs = pairs, hits_tsv = hits_tsv)
}

.foldseek_empty_perhog <- function() tibble::tibble(
  hog_id = character(0), n_pairs = integer(0), n_supported = integer(0),
  n_expected_pairs = integer(0),
  support_frac = numeric(0), mean_tmscore = numeric(0)
)


#' Annotate per-node HOG TSVs with Foldseek structural support
#'
#' Reads every `HOGs/N*_*.tsv` produced by `cut_hogs()` and appends three
#' columns — `n_fsk_pairs`, `fsk_support_frac`, `fsk_mean_tmscore` — so
#' structural corroboration is visible next to the sequence-based HOG
#' membership. Files are rewritten in place. HOGs absent from the
#' Foldseek scoring table get NA in the new columns (typical cause:
#' singleton HOG with no intra-pairs to align).
#'
#' @param out_dir Directory containing `HOGs/`.
#' @param per_hog Tibble from `score_foldseek_hog()`/`run_foldseek_hog()$per_hog`
#'   (`hog_id`, `n_pairs`, `support_frac`, `mean_tmscore`).
#' @return Invisibly: character vector of rewritten paths.
#' @export
annotate_hogs_with_foldseek <- function(out_dir, per_hog) {
  stopifnot(dir.exists(file.path(out_dir, "HOGs")))
  stopifnot("hog_id" %in% names(per_hog))
  n_expected <- if ("n_expected_pairs" %in% names(per_hog))
    per_hog$n_expected_pairs else rep(NA_integer_, nrow(per_hog))
  lookup <- stats::setNames(
    list(per_hog$n_pairs, n_expected, per_hog$support_frac,
         per_hog$mean_tmscore),
    c("n_fsk_pairs", "n_fsk_expected", "fsk_support_frac",
      "fsk_mean_tmscore")
  )
  lookup <- lapply(lookup, function(v) stats::setNames(v, per_hog$hog_id))

  files <- list.files(file.path(out_dir, "HOGs"), pattern = "^N\\d+_.*\\.tsv$",
                      full.names = TRUE)
  rewritten <- character()
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    if (!length(lines)) next
    header_idx <- which(!startsWith(lines, "#"))[1]
    if (is.na(header_idx)) next
    preamble <- lines[seq_len(header_idx - 1L)]
    header   <- lines[header_idx]
    body     <- lines[-seq_len(header_idx)]
    if (grepl("\\tfsk_support_frac\\b", header)) next  # already annotated
    header <- paste(header, "n_fsk_pairs", "n_fsk_expected",
                    "fsk_support_frac", "fsk_mean_tmscore", sep = "\t")
    if (length(body)) {
      hog_ids <- vapply(strsplit(body, "\t", fixed = TRUE), `[`,
                        character(1), 1L)
      n_p  <- unname(lookup$n_fsk_pairs[hog_ids])
      n_e  <- unname(lookup$n_fsk_expected[hog_ids])
      sf   <- unname(lookup$fsk_support_frac[hog_ids])
      tm   <- unname(lookup$fsk_mean_tmscore[hog_ids])
      body <- paste(body,
                    ifelse(is.na(n_p), "", format(n_p)),
                    ifelse(is.na(n_e), "", format(n_e)),
                    ifelse(is.na(sf),  "", formatC(sf, format = "g", digits = 4)),
                    ifelse(is.na(tm),  "", formatC(tm, format = "g", digits = 4)),
                    sep = "\t")
    }
    writeLines(c(preamble, header, body), f)
    rewritten <- c(rewritten, f)
  }
  invisible(rewritten)
}
