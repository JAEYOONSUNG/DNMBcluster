#' Back-translate AA alignments to codon alignments (pure-R pal2nal)
#'
#' DNMBcluster's main pipeline keeps only amino-acid translations, so any
#' codon-level analysis (codeml M0 dN/dS, HyPhy FEL/BUSTED, RELAX) needs
#' a separate back-translation step. This helper takes a directory of
#' per-OG AA alignments written by `per_og_trees()` and a nucleotide CDS
#' FASTA (or directory of per-genome FASTAs) supplied by the user, then
#' writes `aln.nuc.fasta` next to each `aln.fasta`, threading codons
#' under the AA alignment column by column.
#'
#' Features:
#' \itemize{
#'   \item Codon-aware gap handling — one AA gap becomes one "---" codon.
#'   \item Stop-codon trim: if the CDS ends in a stop codon, it is dropped
#'         before back-translation so alignment columns correspond 1:1.
#'   \item AA mismatch check: each back-translated codon is decoded with
#'         the standard (or user-supplied) genetic code and compared to
#'         the AA column; mismatched OGs are reported to
#'         `back_translate_summary.tsv` and skipped rather than silently
#'         producing junk codons.
#'   \item No external binary required — works inside the R package.
#' }
#'
#' @param og_result Output of `per_og_trees()`. Must contain `aln_path`.
#' @param cds_source Either (a) a path to a single nucleotide FASTA whose
#'   headers match the alignment headers (`p<uid>_g<genome_uid>`), or
#'   (b) a directory containing one FASTA per genome (`.fasta`/`.fna`),
#'   where headers carry the protein_uid in the first whitespace-delimited
#'   token (`p<uid>` or `<uid>`).
#' @param out_suffix File name written inside each OG dir. Default
#'   `"aln.nuc.fasta"`. Codeml's PHYLIP form is produced by
#'   `write_codon_phylip()` if requested.
#' @param genetic_code Integer PAML genetic-code ID (0 = standard, 4 =
#'   vertebrate mt, 11 = bacterial). Only 0 and 11 are implemented in
#'   pure R; 11 differs only in a handful of stop/start codons and is
#'   the right choice for DNMBcluster's bacterial pangenomes. Default 11.
#' @param write_phylip If TRUE, also write `aln.nuc.phy` in PAML PHYLIP
#'   sequential form so `run_codeml_m0()` / HyPhy can consume it directly.
#' @param verbose If TRUE, prints one line per OG with the outcome.
#' @return Tibble: `cluster_id, nuc_path, phy_path, n_seqs, n_codons,
#'   ok, reason`.
#' @export
back_translate_og_codons <- function(og_result,
                                      cds_source,
                                      out_suffix = "aln.nuc.fasta",
                                      genetic_code = 11L,
                                      write_phylip = TRUE,
                                      verbose = FALSE) {
  stopifnot("aln_path" %in% names(og_result))
  table <- .read_cds_index(cds_source)
  tt    <- .codon_table(genetic_code)

  ok <- og_result[!is.na(og_result$aln_path) &
                    file.exists(og_result$aln_path), , drop = FALSE]
  rows <- vector("list", nrow(ok))
  for (i in seq_len(nrow(ok))) {
    cid      <- ok$cluster_id[i]
    aln_path <- ok$aln_path[i]
    og_dir   <- dirname(aln_path)
    nuc_path <- file.path(og_dir, out_suffix)
    phy_path <- file.path(og_dir, "aln.nuc.phy")
    res <- tryCatch(
      .back_translate_one(aln_path, table, tt),
      error = function(e) list(ok = FALSE, reason = conditionMessage(e))
    )
    if (isTRUE(res$ok)) {
      .write_codon_fasta(res$names, res$codons, nuc_path)
      if (isTRUE(write_phylip))
        .write_codon_phylip(res$names, res$codons, phy_path)
    }
    rows[[i]] <- tibble::tibble(
      cluster_id = cid,
      nuc_path   = if (isTRUE(res$ok)) nuc_path else NA_character_,
      phy_path   = if (isTRUE(res$ok) && isTRUE(write_phylip))
                     phy_path else NA_character_,
      n_seqs    = if (isTRUE(res$ok)) length(res$names) else NA_integer_,
      n_codons  = if (isTRUE(res$ok)) res$n_codons       else NA_integer_,
      ok        = isTRUE(res$ok),
      reason    = if (isTRUE(res$ok)) NA_character_ else res$reason
    )
    if (verbose) message(sprintf("[back_translate] OG %d: %s",
                                  cid,
                                  if (isTRUE(res$ok))
                                    sprintf("%d seqs x %d codons",
                                            length(res$names), res$n_codons)
                                  else paste("skipped —", res$reason)))
  }
  dplyr::bind_rows(rows)
}


# Build name -> nucleotide sequence lookup from either a single FASTA or
# a directory of FASTAs. Accepts headers of the form "p<uid>_g<gid>",
# "p<uid>", or plain "<uid>". Returns a character vector keyed by the
# fully-qualified "p<uid>_g<gid>" header if present, else by the uid-only
# form. The alignment header is tried in both forms at lookup time.
.read_cds_index <- function(cds_source) {
  files <- if (dir.exists(cds_source))
    list.files(cds_source, pattern = "\\.(fa|fna|fasta)(\\.gz)?$",
               full.names = TRUE, ignore.case = TRUE)
  else cds_source
  if (!length(files) || !all(file.exists(files)))
    stop("cds_source is neither a readable FASTA nor a dir with FASTAs: ",
         cds_source)

  lookup <- new.env(parent = emptyenv())
  for (f in files) {
    con <- if (grepl("\\.gz$", f)) gzfile(f, "r") else file(f, "r")
    header <- NULL; buf <- character()
    flush <- function() {
      if (is.null(header)) return(invisible())
      seq <- paste(buf, collapse = "")
      # strip whitespace, keep case for stop detection; downstream .toupper
      seq <- gsub("\\s", "", seq)
      # Accept header forms: "p<uid>_g<gid>", "p<uid>", "<uid>".
      hdr <- strsplit(header, "\\s+")[[1]][1]
      # Full form wins if present.
      if (grepl("^p\\d+_g\\d+$", hdr)) {
        lookup[[hdr]] <- seq
      } else if (grepl("^p\\d+$", hdr)) {
        lookup[[hdr]] <- seq
      } else if (grepl("^\\d+$", hdr)) {
        lookup[[paste0("p", hdr)]] <- seq
      } else {
        lookup[[hdr]] <- seq
      }
      header <<- NULL; buf <<- character()
    }
    while (length(ln <- readLines(con, n = 1L, warn = FALSE))) {
      if (startsWith(ln, ">")) {
        flush()
        header <- sub("^>", "", ln)
      } else if (!is.null(header)) {
        buf <- c(buf, ln)
      }
    }
    flush(); close(con)
  }
  lookup
}


# Return a named character vector mapping each 64-codon string to its
# single-letter AA under the requested genetic code. Currently implements
# NCBI code 1 (standard) and 11 (bacterial/archaeal/plant plastid).
.codon_table <- function(genetic_code = 11L) {
  if (!genetic_code %in% c(0L, 1L, 11L))
    stop("back_translate: only NCBI codes 1 and 11 are implemented (got ",
         genetic_code, ")")
  bases <- c("T", "C", "A", "G")
  codons <- character(64)
  k <- 0L
  for (a in bases) for (b in bases) for (c in bases) {
    k <- k + 1L
    codons[k] <- paste0(a, b, c)
  }
  # Standard code 1.
  aa <- strsplit(
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ""
  )[[1]]
  if (identical(genetic_code, 11L) || identical(genetic_code, 0L)) {
    # NCBI 11 differs from 1 only in alternate *start* codons (TTG, CTG,
    # GTG, ATA, ATT, ATC all treated as M when at position 1). The stop
    # set is identical to standard (TAA, TAG, TGA). For back-translation
    # purposes (AA verification mid-alignment) the standard table is
    # correct; start-codon M detection is handled separately on the
    # leftmost column.
  }
  stats::setNames(aa, codons)
}


.back_translate_one <- function(aln_path, cds_lookup, codon_tbl) {
  fa <- .read_fasta(aln_path)
  if (!length(fa$names)) stop("empty alignment")
  aa_len <- unique(nchar(fa$seqs))
  if (length(aa_len) != 1L)
    stop("alignment rows have inconsistent lengths (not aligned)")

  codon_rows <- character(length(fa$names))
  for (i in seq_along(fa$names)) {
    nm <- fa$names[i]
    nuc <- .lookup_cds(cds_lookup, nm)
    if (is.null(nuc)) stop("CDS missing for ", nm)
    nuc <- toupper(gsub("U", "T", nuc))
    nuc <- gsub("[^ACGTN]", "N", nuc)
    # Trim trailing stop.
    if (nchar(nuc) >= 3L) {
      tail_cd <- substr(nuc, nchar(nuc) - 2L, nchar(nuc))
      if (tail_cd %in% c("TAA", "TAG", "TGA"))
        nuc <- substr(nuc, 1L, nchar(nuc) - 3L)
    }
    if (nchar(nuc) %% 3L != 0L)
      stop("CDS for ", nm, " not divisible by 3 after stop trim")
    n_expected_aa <- sum(strsplit(fa$seqs[i], "")[[1]] != "-")
    if (nchar(nuc) / 3L != n_expected_aa)
      stop(sprintf("CDS/AA length mismatch for %s: %d codons vs %d non-gap AA",
                    nm, as.integer(nchar(nuc) / 3L), n_expected_aa))
    aa_vec <- strsplit(fa$seqs[i], "")[[1]]
    out <- character(length(aa_vec))
    p <- 1L
    for (j in seq_along(aa_vec)) {
      if (aa_vec[j] == "-") {
        out[j] <- "---"
      } else {
        cd <- substr(nuc, p, p + 2L)
        aa_decoded <- codon_tbl[cd]
        # Mid-alignment mismatch => user gave mismatched AA/CDS pair.
        # Allow first-column M vs alt-start codons under code 11.
        first_col <- j == 1L
        alt_starts <- c("TTG", "CTG", "GTG", "ATA", "ATT", "ATC")
        if (!is.na(aa_decoded) && aa_decoded != aa_vec[j]) {
          if (!(first_col && aa_vec[j] == "M" && cd %in% alt_starts)) {
            stop(sprintf(
              "CDS/AA mismatch for %s at col %d: codon %s decodes to %s but AA is %s",
              nm, j, cd, aa_decoded, aa_vec[j]))
          }
        }
        out[j] <- cd
        p <- p + 3L
      }
    }
    codon_rows[i] <- paste(out, collapse = "")
  }
  list(ok = TRUE, names = fa$names, codons = codon_rows,
       n_codons = aa_len)
}


.lookup_cds <- function(env, header_name) {
  hdr <- strsplit(header_name, "\\s+")[[1]][1]
  if (!is.null(env[[hdr]])) return(env[[hdr]])
  # try uid-only form from "p<uid>_g<gid>".
  puid <- sub("_g.*$", "", hdr)
  if (!is.null(env[[puid]])) return(env[[puid]])
  # strip leading p
  bare <- sub("^p", "", puid)
  if (!is.null(env[[bare]])) return(env[[bare]])
  NULL
}


.read_fasta <- function(path) {
  ln <- readLines(path, warn = FALSE)
  idx <- which(startsWith(ln, ">"))
  if (!length(idx)) return(list(names = character(0), seqs = character(0)))
  ends <- c(idx[-1] - 1L, length(ln))
  seqs <- vapply(seq_along(idx), function(i)
    paste(ln[(idx[i] + 1L):ends[i]], collapse = ""), character(1))
  list(names = sub("^>", "", ln[idx]), seqs = seqs)
}


.write_codon_fasta <- function(names, codons, path) {
  out <- character(2L * length(names))
  out[seq(1L, length(out), by = 2L)] <- paste0(">", names)
  out[seq(2L, length(out), by = 2L)] <- codons
  writeLines(out, path)
}


.write_codon_phylip <- function(names, codons, path) {
  n <- length(names); L <- nchar(codons[1])
  # PAML sequential PHYLIP: two spaces, N L; then "name  seq".
  header <- sprintf(" %d  %d", n, L)
  body <- paste0(format(names, width = max(10L, max(nchar(names)))),
                 "  ", codons)
  writeLines(c(header, body), path)
}
