#' Run HyPhy FEL / BUSTED / aBSREL per OG and parse the JSON results
#'
#' Modern SOTA selection tests for bacterial pan-genomes. Post-2023 high-
#' impact microbial-evolution papers treat HyPhy 2.5 as the default and
#' cite codeml M0 only as a sensitivity check (Murrell & Kosakovsky Pond;
#' Wolf 2021 PNAS dN/dS-inflation; Jeffares MBE 2023). DNMBcluster keeps
#' codeml M0 for fast baseline omega and adds these wrappers for
#'
#' \itemize{
#'   \item `run_hyphy_fel()`   — per-site pervasive selection (alpha vs beta
#'         with LRT p-value), replaces codeml M7/M8.
#'   \item `run_hyphy_busted()` — gene-wide episodic selection test with
#'         synonymous rate variation (BUSTED[S]); power/type-I exceeds
#'         codeml branch-site Model A on the short-branch alignments
#'         typical of bacterial pan-genomes.
#'   \item `run_hyphy_absrel()` — adaptive branch-site REL, per-branch
#'         selection without a priori foreground. Current default for
#'         branch-level tests in microbial papers.
#' }
#'
#' Every wrapper follows DNMBcluster's binary-absent convention: if
#' `hyphy` is not on PATH, it warns once and returns NULL so a pipeline
#' run continues. Per-OG directories receive their own JSON + log, and
#' the returned tibble mirrors `run_codeml_m0()`'s column discipline
#' (`cluster_id`, a small set of scalar summaries, `converged`, `reason`).
#'
#' Input expectation: each `og_dir` contains a codon alignment
#' (`aln.nuc.fasta` by default, FASTA — HyPhy does not take PHYLIP
#' directly) and a rooted gene tree (`tree.nwk`). The codon MSA is what
#' `back_translate_og_codons()` writes.
#'
#' @name run_hyphy
NULL


#' @rdname run_hyphy
#' @param og_dirs Character vector of per-OG directories.
#' @param aln_name File name of the codon FASTA alignment inside each OG
#'   dir. Default `"aln.nuc.fasta"`.
#' @param tree_name File name of the rooted gene tree. Default `"tree.nwk"`.
#' @param work_dir Output root; each OG lands under
#'   `work_dir/OG_<id>/<method>.json`.
#' @param srv If TRUE (default), BUSTED uses synonymous-rate variation
#'   (BUSTED[S]) — the recommended mode for bacteria.
#' @param threads Passed to `hyphy` via `CPU=<n>`; HyPhy natively threads
#'   a single analysis. Default 1.
#' @param binary Explicit path to `hyphy`. NULL → `Sys.which`.
#' @return Tibble with per-OG summaries; NULL if binary missing.
#' @export
run_hyphy_fel <- function(og_dirs,
                          aln_name  = "aln.nuc.fasta",
                          tree_name = "tree.nwk",
                          work_dir,
                          threads   = 1L,
                          binary    = NULL) {
  .hyphy_batch("fel", og_dirs, aln_name, tree_name, work_dir,
               extra_args = character(0), threads = threads, binary = binary,
               parser = .parse_hyphy_fel)
}


#' @rdname run_hyphy
#' @export
run_hyphy_busted <- function(og_dirs,
                             aln_name  = "aln.nuc.fasta",
                             tree_name = "tree.nwk",
                             work_dir,
                             srv       = TRUE,
                             threads   = 1L,
                             binary    = NULL) {
  # BUSTED[S] is controlled by --srv Yes/No in HyPhy >=2.5.52.
  extra <- c("--srv", if (isTRUE(srv)) "Yes" else "No")
  .hyphy_batch("busted", og_dirs, aln_name, tree_name, work_dir,
               extra_args = extra, threads = threads, binary = binary,
               parser = .parse_hyphy_busted)
}


#' @rdname run_hyphy
#' @export
run_hyphy_absrel <- function(og_dirs,
                             aln_name  = "aln.nuc.fasta",
                             tree_name = "tree.nwk",
                             work_dir,
                             threads   = 1L,
                             binary    = NULL) {
  .hyphy_batch("absrel", og_dirs, aln_name, tree_name, work_dir,
               extra_args = character(0), threads = threads, binary = binary,
               parser = .parse_hyphy_absrel)
}


#' Run FEL + BUSTED + aBSREL in a single pass, in parallel, with FDR
#'
#' Wraps `run_hyphy_fel()`, `run_hyphy_busted()`, and `run_hyphy_absrel()`
#' into one call that:
#' \itemize{
#'   \item processes OGs in parallel via `future.apply::future_lapply`
#'         (one worker per OG; HyPhy is threaded internally too but
#'         `workers * per-analysis threads` dominates at scale);
#'   \item returns a list with one tibble per method plus a `combined`
#'         wide-format summary keyed by `cluster_id`;
#'   \item applies Benjamini-Hochberg correction across OGs for the
#'         BUSTED gene-wide p-value and the aBSREL minimum-corrected p —
#'         the two multiple-testing corrections a 2025 reviewer expects
#'         to see on a pan-genome scan.
#' }
#'
#' @param og_dirs Character vector of per-OG directories (each with
#'   `aln.nuc.fasta` + `tree.nwk`, typically from
#'   `back_translate_og_codons()` + `per_og_trees()`).
#' @param work_dir Output root. Per-OG subdirs land under
#'   `work_dir/<method>/OG_<id>/`.
#' @param methods Subset of `c("fel","busted","absrel")`. Default all.
#' @param aln_name,tree_name Alignment / tree file names inside each OG dir.
#' @param srv BUSTED synonymous-rate-variation flag. Default TRUE
#'   (BUSTED[S] — recommended for bacteria).
#' @param workers Parallel OG workers. 1 = serial. > 1 uses
#'   `future.apply::future_lapply` when the package is installed and sets
#'   `future::plan(multisession, workers = workers)` unless the caller
#'   has a plan installed already (`"multisession"` check).
#' @param threads HyPhy `CPU=<n>` per analysis. Default 1. With
#'   `workers = k` and `threads = 1` you get k-wide parallelism across OGs;
#'   inverting is also fine for very large families.
#' @param binary Explicit path to `hyphy`. NULL → `Sys.which`.
#' @return List with elements `fel`, `busted`, `absrel`, `combined`, and
#'   `work_dir`. Each method tibble includes a `q_value` column (BH-
#'   corrected across OGs) where applicable. NULL with warning if hyphy
#'   is absent.
#' @export
run_hyphy_batch <- function(og_dirs,
                            work_dir,
                            methods   = c("fel", "busted", "absrel"),
                            aln_name  = "aln.nuc.fasta",
                            tree_name = "tree.nwk",
                            srv       = TRUE,
                            workers   = 1L,
                            threads   = 1L,
                            binary    = NULL) {
  methods <- match.arg(methods, several.ok = TRUE)
  bin <- if (!is.null(binary) && nzchar(binary)) binary else
           unname(Sys.which("hyphy"))
  if (!nzchar(bin)) {
    warning("[run_hyphy_batch] 'hyphy' binary not found on PATH; skipping.")
    return(NULL)
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  # Spawn a multisession plan only if workers>1 and the user has not set
  # one already — avoids stomping on caller-configured plans (SLURM, etc).
  use_parallel <- isTRUE(workers > 1L) &&
                  requireNamespace("future.apply", quietly = TRUE)
  if (use_parallel) {
    if (!inherits(future::plan(), "multisession") &&
        !inherits(future::plan(), "cluster") &&
        !inherits(future::plan(), "multicore")) {
      old_plan <- future::plan(future::multisession, workers = workers)
      on.exit(future::plan(old_plan), add = TRUE)
    }
  }

  run_one_method <- function(method) {
    fn <- switch(method,
                 fel    = run_hyphy_fel,
                 busted = run_hyphy_busted,
                 absrel = run_hyphy_absrel)
    method_wd <- file.path(work_dir, method)
    dir.create(method_wd, showWarnings = FALSE, recursive = TRUE)
    if (use_parallel) {
      rows <- future.apply::future_lapply(og_dirs, function(d) {
        args <- list(og_dirs = d, aln_name = aln_name,
                     tree_name = tree_name, work_dir = method_wd,
                     threads = threads, binary = bin)
        if (identical(method, "busted")) args$srv <- srv
        do.call(fn, args)
      }, future.seed = TRUE)
      dplyr::bind_rows(rows)
    } else {
      args <- list(og_dirs = og_dirs, aln_name = aln_name,
                   tree_name = tree_name, work_dir = method_wd,
                   threads = threads, binary = bin)
      if ("busted" %in% method) args$srv <- srv
      do.call(fn, args)
    }
  }

  results <- list(fel = NULL, busted = NULL, absrel = NULL)
  for (m in methods) results[[m]] <- run_one_method(m)

  # BH correction across OGs, per method. FEL is per-site (no OG-level
  # p), so we correct its n_positive/n_negative as-is and instead flag
  # OGs with any positive site at FDR<=0.1 by BH-correcting the vector
  # of min site p-values if available; keeping the wrapper simple, we
  # only apply BH to the two OG-scalar tests.
  if (!is.null(results$busted) && "p_value" %in% names(results$busted)) {
    results$busted$q_value <-
      stats::p.adjust(results$busted$p_value, method = "BH")
  }
  if (!is.null(results$absrel) && "min_corrected_p" %in% names(results$absrel)) {
    results$absrel$q_value <-
      stats::p.adjust(results$absrel$min_corrected_p, method = "BH")
  }

  # Wide-format combined summary, one row per cluster_id.
  combined <- .hyphy_combine(results$fel, results$busted, results$absrel)
  results$combined <- combined
  results$work_dir <- work_dir
  results
}


.hyphy_combine <- function(fel, busted, absrel) {
  cids <- sort(unique(c(
    if (!is.null(fel))    fel$cluster_id    else integer(0),
    if (!is.null(busted)) busted$cluster_id else integer(0),
    if (!is.null(absrel)) absrel$cluster_id else integer(0)
  )))
  if (!length(cids))
    return(tibble::tibble(cluster_id = integer(0)))
  out <- tibble::tibble(cluster_id = cids)
  if (!is.null(fel)) {
    f <- dplyr::select(fel,
                       cluster_id,
                       fel_n_sites    = n_sites,
                       fel_n_positive = n_positive,
                       fel_n_negative = n_negative,
                       fel_median_omega = median_omega)
    out <- dplyr::left_join(out, f, by = "cluster_id")
  }
  if (!is.null(busted)) {
    b <- dplyr::select(busted,
                       cluster_id,
                       busted_p      = p_value,
                       busted_q      = dplyr::any_of("q_value"),
                       busted_lrt    = lrt,
                       busted_omega_pos = omega_positive,
                       busted_w_pos  = weight_positive)
    out <- dplyr::left_join(out, b, by = "cluster_id")
  }
  if (!is.null(absrel)) {
    a <- dplyr::select(absrel,
                       cluster_id,
                       absrel_n_branches  = n_branches,
                       absrel_n_selected  = n_selected,
                       absrel_min_p       = min_corrected_p,
                       absrel_q           = dplyr::any_of("q_value"))
    out <- dplyr::left_join(out, a, by = "cluster_id")
  }
  out
}


# ----------------- internals ---------------------------------------------

.hyphy_batch <- function(method, og_dirs, aln_name, tree_name, work_dir,
                         extra_args, threads, binary, parser) {
  bin <- if (!is.null(binary) && nzchar(binary)) binary else
           unname(Sys.which("hyphy"))
  if (!nzchar(bin)) {
    warning(sprintf("[run_hyphy_%s] 'hyphy' binary not found on PATH; skipping.",
                    method))
    return(NULL)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning(sprintf("[run_hyphy_%s] package 'jsonlite' required to parse HyPhy output; skipping.",
                    method))
    return(NULL)
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- vector("list", length(og_dirs))
  for (i in seq_along(og_dirs)) {
    ogd <- og_dirs[i]
    cid <- .hyphy_cluster_id(ogd)
    aln <- file.path(ogd, aln_name)
    tre <- file.path(ogd, tree_name)
    sub_wd <- file.path(work_dir, sprintf("OG_%07d", cid))
    dir.create(sub_wd, showWarnings = FALSE)
    json <- file.path(sub_wd, sprintf("%s.json", method))
    log_path <- file.path(sub_wd, sprintf("%s.log", method))

    if (!file.exists(aln) || !file.exists(tre)) {
      rows[[i]] <- parser(cid, NULL,
                         reason = "missing codon alignment or tree")
      next
    }

    # Stage under safe basenames to dodge whitespace-in-path issues.
    staged_aln <- file.path(sub_wd, "in.aln.fasta")
    staged_tre <- file.path(sub_wd, "in.tree.nwk")
    file.copy(aln, staged_aln, overwrite = TRUE)
    file.copy(tre, staged_tre, overwrite = TRUE)

    args <- c(
      sprintf("CPU=%d", as.integer(threads)),
      method,
      "--alignment", staged_aln,
      "--tree",      staged_tre,
      "--output",    json,
      extra_args
    )
    rc <- system2(bin, args = args, stdout = log_path, stderr = log_path)
    if (rc != 0L || !file.exists(json)) {
      rows[[i]] <- parser(cid, NULL,
                         reason = sprintf("hyphy %s exit=%d", method, rc))
      next
    }
    parsed <- tryCatch(parser(cid, json, reason = NA_character_),
                       error = function(e)
                         parser(cid, NULL,
                                reason = paste("parse error:", conditionMessage(e))))
    rows[[i]] <- parsed
  }
  dplyr::bind_rows(rows)
}


.hyphy_cluster_id <- function(og_dir) {
  m <- regmatches(basename(og_dir), regexpr("[0-9]+$", basename(og_dir)))
  if (length(m) && nzchar(m)) as.integer(m) else NA_integer_
}


# FEL JSON: MLE$content[["0"]] is a list of per-site rows; column order
# per MLE$headers is [alpha, beta, alpha=beta (shared), LRT, p-value,
# total branch length, ...]. We summarise n_sites, n_positive (beta>alpha
# with p<=0.1), n_negative (alpha>beta with p<=0.1), median omega over
# sites with alpha>0. The 0.1 default is HyPhy's own site-level cutoff.
.parse_hyphy_fel <- function(cid, json_path, reason = NA_character_) {
  empty <- tibble::tibble(
    cluster_id = cid, method = "FEL",
    n_sites = NA_integer_, n_positive = NA_integer_,
    n_negative = NA_integer_, median_omega = NA_real_,
    json_path = NA_character_, converged = FALSE,
    reason = reason
  )
  if (is.null(json_path) || !file.exists(json_path)) return(empty)
  j <- jsonlite::fromJSON(json_path, simplifyVector = TRUE)
  mle <- j[["MLE"]]
  if (is.null(mle) || is.null(mle$content)) return(empty)
  mat <- mle$content[["0"]]
  if (is.null(mat) || !length(mat)) return(empty)
  if (!is.matrix(mat)) mat <- do.call(rbind, mat)
  # Column indices are stable across HyPhy 2.5: alpha=1, beta=2, p=5.
  alpha <- as.numeric(mat[, 1L])
  beta  <- as.numeric(mat[, 2L])
  pval  <- as.numeric(mat[, 5L])
  n <- length(alpha)
  pos <- sum(beta > alpha & pval <= 0.1, na.rm = TRUE)
  neg <- sum(alpha > beta & pval <= 0.1, na.rm = TRUE)
  omega <- ifelse(alpha > 0, beta / alpha, NA_real_)
  tibble::tibble(
    cluster_id = cid, method = "FEL",
    n_sites = n, n_positive = pos, n_negative = neg,
    median_omega = if (any(!is.na(omega))) stats::median(omega, na.rm = TRUE) else NA_real_,
    json_path = json_path, converged = TRUE,
    reason = NA_character_
  )
}


# BUSTED JSON: episodic-selection test result lives under
# `test results` -> `p-value` (plus LRT); per-branch / per-rate-class fits
# are under `fits`. We pull the gene-wide p-value and omega of the
# positive-rate class if present.
.parse_hyphy_busted <- function(cid, json_path, reason = NA_character_) {
  empty <- tibble::tibble(
    cluster_id = cid, method = "BUSTED",
    p_value = NA_real_, lrt = NA_real_,
    omega_positive = NA_real_, weight_positive = NA_real_,
    json_path = NA_character_, converged = FALSE,
    reason = reason
  )
  if (is.null(json_path) || !file.exists(json_path)) return(empty)
  j <- jsonlite::fromJSON(json_path, simplifyVector = TRUE)
  tr <- j[["test results"]]
  p  <- if (!is.null(tr) && "p-value" %in% names(tr))
          as.numeric(tr[["p-value"]]) else NA_real_
  lrt <- if (!is.null(tr) && "LRT" %in% names(tr))
           as.numeric(tr[["LRT"]]) else NA_real_
  # Positive-rate class of the unconstrained model.
  fits <- j[["fits"]]
  unc  <- fits[["Unconstrained model"]]
  om_pos <- NA_real_; w_pos <- NA_real_
  if (!is.null(unc)) {
    rates <- unc[["Rate Distributions"]][["Test"]]
    if (!is.null(rates)) {
      # rates is a data.frame with columns "omega" and "proportion"
      # (column names historically bracketed "omega" vs "0"/"1"/"2").
      omega_col <- if ("omega" %in% names(rates)) rates[["omega"]] else rates[[1L]]
      prop_col  <- if ("proportion" %in% names(rates)) rates[["proportion"]] else rates[[2L]]
      omega_col <- as.numeric(omega_col); prop_col <- as.numeric(prop_col)
      hi <- which.max(omega_col)
      if (length(hi)) {
        om_pos <- omega_col[hi]; w_pos <- prop_col[hi]
      }
    }
  }
  tibble::tibble(
    cluster_id = cid, method = "BUSTED",
    p_value = p, lrt = lrt,
    omega_positive = om_pos, weight_positive = w_pos,
    json_path = json_path,
    converged = !is.na(p),
    reason = if (!is.na(p)) NA_character_ else "no test result in JSON"
  )
}


# aBSREL JSON: `branch attributes$0` is a named list per branch, each
# with fields `Corrected P-value` / `Uncorrected P-value` / `Rate classes`.
# We count branches significant at 0.05 after holm correction (aBSREL's
# own default) and report the min corrected p-value.
.parse_hyphy_absrel <- function(cid, json_path, reason = NA_character_) {
  empty <- tibble::tibble(
    cluster_id = cid, method = "aBSREL",
    n_branches = NA_integer_, n_selected = NA_integer_,
    min_corrected_p = NA_real_,
    json_path = NA_character_, converged = FALSE,
    reason = reason
  )
  if (is.null(json_path) || !file.exists(json_path)) return(empty)
  j <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)
  ba <- j[["branch attributes"]]
  if (is.null(ba)) return(empty)
  first <- ba[[1L]]
  if (is.null(first)) return(empty)
  pvals <- vapply(first, function(b) {
    v <- b[["Corrected P-value"]]
    if (is.null(v)) NA_real_ else as.numeric(v)
  }, numeric(1))
  n_br <- length(pvals)
  n_sel <- sum(pvals <= 0.05, na.rm = TRUE)
  min_p <- if (any(!is.na(pvals))) min(pvals, na.rm = TRUE) else NA_real_
  tibble::tibble(
    cluster_id = cid, method = "aBSREL",
    n_branches = n_br, n_selected = n_sel,
    min_corrected_p = min_p,
    json_path = json_path,
    converged = n_br > 0L,
    reason = if (n_br > 0L) NA_character_ else "no branch attributes in JSON"
  )
}
