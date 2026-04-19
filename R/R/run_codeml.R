#' Run PAML codeml site-model M0 per single-copy OG and parse omega
#'
#' Given a directory of single-copy-core OGs with codon-aware nucleotide
#' alignments (PHYLIP or FASTA) and rooted gene trees, runs codeml with
#' model M0 (one-ratio across sites and lineages) and returns the
#' maximum-likelihood dN/dS estimate per OG. This is the classical
#' "background omega" comparative-genomics summary and the cheapest
#' PAML flavour; site/branch models are a separate tool.
#'
#' The function is deliberately defensive: missing binary, missing
#' alignments, or unparseable output return NA for that OG rather than
#' aborting the batch, so a single broken family does not tank the
#' pipeline.
#'
#' @param og_dirs Character vector of per-OG directories, each holding
#'   a nucleotide alignment and a gene tree (see `aln_name` / `tree_name`).
#'   Typically one of the `orthogroups/OG_*` subdirs, filtered to the
#'   single-copy-core rows. Multi-copy families are silently skipped
#'   (codeml M0 assumes one tip per species).
#' @param aln_name File name inside each OG dir of the codon alignment
#'   (PAML PHYLIP interleaved or sequential). Default `"aln.nuc.phy"`.
#' @param tree_name File name of the rooted gene tree. Default `"tree.nwk"`.
#' @param work_dir Working directory for per-OG control files and
#'   outputs. Created if missing.
#' @param threads Reserved for future parallelism. Currently serial —
#'   codeml itself is single-threaded.
#' @param binary Explicit path to `codeml`. NULL → `Sys.which`.
#' @return Tibble: `cluster_id, omega, dN, dS, log_lik, converged`.
#'   NULL with a warning when codeml is not on PATH.
#' @export
run_codeml_m0 <- function(og_dirs,
                          aln_name  = "aln.nuc.phy",
                          tree_name = "tree.nwk",
                          work_dir,
                          threads   = 1L,
                          binary    = NULL) {
  bin <- if (!is.null(binary) && nzchar(binary)) binary else unname(Sys.which("codeml"))
  if (!nzchar(bin)) {
    warning("[run_codeml_m0] 'codeml' binary not found on PATH; skipping.")
    return(NULL)
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- vector("list", length(og_dirs))
  for (i in seq_along(og_dirs)) {
    ogd <- og_dirs[i]
    cid <- .codeml_cluster_id(ogd)
    aln <- file.path(ogd, aln_name)
    tre <- file.path(ogd, tree_name)
    if (!file.exists(aln) || !file.exists(tre)) {
      rows[[i]] <- .codeml_na(cid, reason = "missing alignment or tree file")
      next
    }
    sub_wd <- file.path(work_dir, sprintf("OG_%07d", cid))
    dir.create(sub_wd, showWarnings = FALSE)
    # Stage inputs under safe basenames inside sub_wd so (a) the ctl can
    # reference short relative paths regardless of caller cwd, and (b)
    # spaces in the caller's path (PAML's ctl parser splits on whitespace)
    # do not corrupt the control file.
    staged_aln <- file.path(sub_wd, "in.aln")
    staged_tre <- file.path(sub_wd, "in.tree")
    file.copy(aln, staged_aln, overwrite = TRUE)
    file.copy(tre, staged_tre, overwrite = TRUE)
    check <- .codeml_validate_input(staged_aln, staged_tre)
    if (!check$ok) {
      writeLines(check$reason, file.path(sub_wd, "skip_reason.txt"))
      rows[[i]] <- .codeml_na(cid, reason = check$reason)
      next
    }
    ctl <- .write_codeml_m0_ctl(sub_wd, "in.aln", "in.tree")
    rc <- .run_codeml(bin, ctl, sub_wd)
    rows[[i]] <- if (rc == 0L) .parse_codeml_m0(sub_wd, cid)
                 else .codeml_na(cid, reason = sprintf("codeml exit=%d", rc))
  }
  dplyr::bind_rows(rows)
}


.codeml_na <- function(cid, reason = NA_character_) tibble::tibble(
  cluster_id = cid, omega = NA_real_, dN = NA_real_, dS = NA_real_,
  log_lik = NA_real_, converged = FALSE, reason = reason
)

# Validate that the staged alignment is codon-shaped (length divisible by 3,
# DNA alphabet) and that the tree has one tip per sequence. Returns
# list(ok, reason).
.codeml_validate_input <- function(aln_path, tree_path) {
  ln <- readLines(aln_path, warn = FALSE)
  is_fasta <- any(startsWith(ln, ">"))
  seqs <- if (is_fasta) {
    idx <- which(startsWith(ln, ">"))
    ends <- c(idx[-1] - 1L, length(ln))
    vapply(seq_along(idx), function(i)
      paste(ln[(idx[i] + 1L):ends[i]], collapse = ""), character(1))
  } else {
    # PHYLIP: first token per block line is the id, rest is seq. Skip header.
    body <- ln[-1L]
    body <- body[nzchar(trimws(body))]
    toks <- strsplit(trimws(body), "\\s+", perl = TRUE)
    vapply(toks, function(t) paste(t[-1L], collapse = ""), character(1))
  }
  seqs <- gsub("[\\s-]", "", seqs, perl = TRUE)
  if (!length(seqs) || !any(nzchar(seqs)))
    return(list(ok = FALSE, reason = "empty alignment"))
  # Codon sanity: length divisible by 3 and restricted to ACGTN (case-insensitive).
  if (any(nchar(seqs) %% 3L != 0L))
    return(list(ok = FALSE, reason = "alignment length not divisible by 3 (not a codon alignment)"))
  non_dna_frac <- mean(grepl("[^ACGTNacgtn]", seqs))
  if (non_dna_frac > 0)
    return(list(ok = FALSE, reason = "alignment contains non-DNA characters (AA input?)"))
  tr <- tryCatch(ape::read.tree(tree_path), error = function(e) NULL)
  if (is.null(tr))
    return(list(ok = FALSE, reason = "unreadable tree"))
  if (length(unique(tr$tip.label)) != length(tr$tip.label))
    return(list(ok = FALSE, reason = "tree has duplicated tips (multi-copy family?)"))
  list(ok = TRUE, reason = NA_character_)
}

.codeml_cluster_id <- function(og_dir) {
  m <- regmatches(basename(og_dir), regexpr("[0-9]+$", basename(og_dir)))
  if (length(m) && nzchar(m)) as.integer(m) else NA_integer_
}

.write_codeml_m0_ctl <- function(wd, aln, tre) {
  ctl <- file.path(wd, "codeml.ctl")
  # seqfile/treefile/outfile are written as basenames only: .run_codeml()
  # chdirs into `wd` before launching, and PAML's ctl parser splits on
  # whitespace so absolute paths containing spaces corrupt the control
  # file. Callers stage inputs inside `wd` under safe basenames.
  writeLines(c(
    paste("seqfile =", aln),
    paste("treefile =", tre),
    "outfile = out.txt",
    "noisy = 0", "verbose = 0", "runmode = 0",
    "seqtype = 1",           # codons
    "CodonFreq = 2",         # F3x4
    "model = 0",             # one omega across branches
    "NSsites = 0",           # M0
    "icode = 0",
    "fix_kappa = 0", "kappa = 2",
    "fix_omega = 0", "omega = 0.4",
    "cleandata = 1"
  ), ctl)
  ctl
}

.run_codeml <- function(bin, ctl, wd) {
  # codeml is positional: argv[1] = control file. Run in wd so relative
  # "rst/rub" intermediate files land next to the outputs, and so the
  # ctl's basename-only seqfile/treefile entries resolve correctly.
  old <- getwd(); on.exit(setwd(old), add = TRUE)
  setwd(wd)
  log_path <- "codeml.log"
  system2(bin, args = shQuote(basename(ctl)),
          stdout = log_path, stderr = log_path)
}

.parse_codeml_m0 <- function(wd, cid) {
  out_path <- file.path(wd, "out.txt")
  if (!file.exists(out_path)) return(.codeml_na(cid, "no out.txt produced"))
  lines <- readLines(out_path, warn = FALSE)

  omega <- .grep_num(lines, "omega\\s*\\(dN/dS\\)\\s*=\\s*([0-9.eE+-]+)")
  dn    <- .grep_num(lines, "dN\\s*=\\s*([0-9.eE+-]+)")
  ds    <- .grep_num(lines, "dS\\s*=\\s*([0-9.eE+-]+)")
  ll    <- .grep_num(lines, "lnL.*:\\s*([0-9.eE+-]+)")
  tibble::tibble(
    cluster_id = cid,
    omega      = omega,
    dN         = dn,
    dS         = ds,
    log_lik    = ll,
    converged  = !is.na(omega) && !is.na(ll),
    reason     = if (!is.na(omega) && !is.na(ll)) NA_character_
                 else "unparseable omega/lnL"
  )
}

.grep_num <- function(lines, pattern) {
  hits <- regmatches(lines, regexec(pattern, lines))
  hits <- hits[vapply(hits, length, integer(1)) >= 2L]
  if (!length(hits)) return(NA_real_)
  suppressWarnings(as.numeric(hits[[length(hits)]][2]))
}
