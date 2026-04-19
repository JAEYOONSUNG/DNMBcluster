#' Run GeneRax on a GeneRax-format export bundle and aggregate DTL events
#'
#' Uses the bundle produced by `export_og_bundles(..., tools="generax",
#' species_tree=...)` — i.e. a directory containing `families.txt`, a
#' `species_tree.nwk`, and one subdirectory per family with
#' `aln.fasta`, `mapping.link`, and `tree.nwk`.
#'
#' GeneRax runs a maximum-likelihood reconciliation under the UndatedDTL
#' model, which is the most rigorous DTL inference available short of
#' ALE. The returned tibble mirrors `reconcile_dtl()` so it plugs into
#' `plot_dtl_events()` unchanged.
#'
#' @param bundle_dir Path produced by the GeneRax exporter. Must contain
#'   `families.txt` and `species_tree.nwk`.
#' @param out_dir Output directory (GeneRax's `--prefix`). Created if
#'   missing. Must not pre-exist with stale runs unless `force = TRUE`.
#' @param rec_model `"UndatedDTL"` (default; duplication + transfer +
#'   loss) or `"UndatedDL"` (no transfer, safe for eukaryotes).
#' @param strategy `"SPR"` (default; rearrange gene trees) or `"EVAL"`
#'   (just score the starting trees, faster).
#' @param threads Number of MPI ranks. GeneRax parallelizes via MPI
#'   (not OpenMP / a `--threads` flag), so when `threads > 1` we launch
#'   via `mpiexec -np <threads>` if `mpiexec` is on PATH; otherwise a
#'   one-line warning is emitted and the run proceeds serially.
#' @param force If TRUE, deletes any existing `out_dir` before running.
#' @param binary Explicit path to `generax`. If NULL, looks on PATH.
#' @return Named list: `events` (tibble `cluster_id, event, count`),
#'   `per_og` (tibble `cluster_id, duplications, transfers, losses,
#'   speciations, likelihood`), `out_dir`. NULL with a warning when the
#'   binary is not available.
#' @export
run_generax <- function(bundle_dir,
                        out_dir,
                        rec_model = c("UndatedDTL", "UndatedDL"),
                        strategy  = c("SPR", "EVAL"),
                        threads   = 4L,
                        force     = FALSE,
                        binary    = NULL) {
  rec_model <- match.arg(rec_model)
  strategy  <- match.arg(strategy)
  bin <- if (!is.null(binary) && nzchar(binary)) binary else unname(Sys.which("generax"))
  if (!nzchar(bin)) {
    warning("[run_generax] 'generax' binary not found on PATH; skipping.")
    return(NULL)
  }
  stopifnot(dir.exists(bundle_dir))
  families <- file.path(bundle_dir, "families.txt")
  species  <- file.path(bundle_dir, "species_tree.nwk")
  if (!file.exists(families) || !file.exists(species)) {
    stop("GeneRax bundle missing families.txt/species_tree.nwk in ", bundle_dir)
  }
  if (dir.exists(out_dir)) {
    if (isTRUE(force)) unlink(out_dir, recursive = TRUE)
    else if (length(list.files(out_dir))) {
      stop("out_dir exists and is non-empty: ", out_dir,
           " (pass force=TRUE to overwrite)")
    }
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # NOTE: `--strategy` does not exist in GeneRax 2.1.3; the flag is
  # `--geneSearchStrategy`. We accept the short/legacy form via the R
  # argument but emit the actual CLI flag here.
  args <- c("--families", families,
            "--species-tree", species,
            "--rec-model", rec_model,
            "--geneSearchStrategy", strategy,
            "--per-family-rates",
            "--prefix", out_dir)

  # GeneRax parallelizes via MPI, not via a `--threads` CLI flag. When
  # the user asks for >1 rank and mpiexec is available, launch through
  # it; otherwise warn once and run serially so we don't fail a whole
  # SOTA run just because MPI isn't installed.
  run_bin <- bin
  run_args <- args
  n_rank <- if (is.null(threads)) 1L else as.integer(threads)
  if (!is.na(n_rank) && n_rank > 1L) {
    mpi <- unname(Sys.which("mpiexec"))
    if (nzchar(mpi)) {
      run_bin <- mpi
      run_args <- c("-np", as.character(n_rank), bin, args)
    } else {
      warning("[run_generax] threads > 1 requested but 'mpiexec' is not ",
              "on PATH; running GeneRax serially. Install openmpi to ",
              "enable MPI parallelism.")
    }
  }

  log_path <- file.path(out_dir, "generax.log")
  rc <- system2(run_bin, args = run_args,
                stdout = log_path, stderr = log_path)
  if (rc != 0) {
    tail_msg <- if (file.exists(log_path)) {
      ln <- readLines(log_path, warn = FALSE)
      paste(utils::tail(ln, 40), collapse = "\n")
    } else ""
    stop("generax failed (exit=", rc, "):\n", tail_msg)
  }

  parse_generax_events(out_dir)
}


#' Parse GeneRax reconciliation output into tidy DTL tibbles
#'
#' Walks a GeneRax `--prefix` directory, finds every `reconciliations/
#' <family>_events.txt`, and sums D/T/L/S counts per family. The event
#' table has the same shape as `reconcile_dtl()$events`, so existing
#' plotting helpers work unchanged.
#'
#' @param out_dir GeneRax output directory.
#' @return Named list with `events`, `per_og`, `out_dir`.
#' @export
parse_generax_events <- function(out_dir) {
  stopifnot(dir.exists(out_dir))
  rec_dir <- file.path(out_dir, "reconciliations")
  if (!dir.exists(rec_dir)) rec_dir <- out_dir   # older GeneRax layout
  files <- list.files(rec_dir, pattern = "_events\\.txt$",
                      recursive = TRUE, full.names = TRUE)
  if (!length(files)) {
    warning("[parse_generax_events] no *_events.txt found under ", rec_dir)
    return(list(events = .empty_events(), per_og = .empty_per_og(),
                out_dir = out_dir))
  }
  per_og_rows <- vector("list", length(files))
  evt_rows    <- vector("list", length(files))
  for (i in seq_along(files)) {
    lines <- readLines(files[i], warn = FALSE)
    cid <- .generax_cluster_id(files[i])
    counts <- c(S = 0L, D = 0L, T = 0L, L = 0L)
    for (ln in lines) {
      # Example line: "event: D node=... species=..."
      tok <- regmatches(ln, regexpr("^\\s*(S|D|T|L|SL|TL)\\b", ln))
      if (!length(tok) || !nzchar(tok)) {
        tok <- regmatches(ln, regexpr("event:\\s*(S|D|T|L)", ln))
        tok <- sub("event:\\s*", "", tok)
      }
      tok <- trimws(tok)
      if (tok %in% c("SL", "TL")) tok <- substr(tok, 2L, 2L)  # count the L
      if (nzchar(tok) && tok %in% names(counts)) counts[tok] <- counts[tok] + 1L
    }
    for (e in names(counts)) {
      evt_rows[[length(evt_rows) + 1L]] <- tibble::tibble(
        cluster_id = cid, event = e, count = counts[[e]]
      )
    }
    per_og_rows[[i]] <- tibble::tibble(
      cluster_id   = cid,
      speciations  = counts[["S"]],
      duplications = counts[["D"]],
      transfers    = counts[["T"]],
      losses       = counts[["L"]],
      likelihood   = .generax_likelihood(files[i])
    )
  }
  list(
    events  = dplyr::bind_rows(evt_rows),
    per_og  = dplyr::bind_rows(per_og_rows),
    out_dir = out_dir
  )
}

.empty_events <- function() tibble::tibble(
  cluster_id = integer(0), event = character(0), count = integer(0)
)
.empty_per_og <- function() tibble::tibble(
  cluster_id = integer(0), speciations = integer(0), duplications = integer(0),
  transfers = integer(0), losses = integer(0), likelihood = numeric(0)
)
.generax_cluster_id <- function(path) {
  nm <- sub("_events\\.txt$", "", basename(path))
  # OG_0000042 → 42; also accept plain integer filenames.
  m <- regmatches(nm, regexpr("[0-9]+$", nm))
  if (length(m) && nzchar(m)) as.integer(m) else NA_integer_
}
.generax_likelihood <- function(events_path) {
  # GeneRax writes a sibling <family>.ll or the likelihood in a log;
  # both are optional. Return NA silently if neither is present.
  ll_path <- sub("_events\\.txt$", ".ll", events_path)
  if (file.exists(ll_path)) {
    v <- suppressWarnings(as.numeric(trimws(readLines(ll_path, warn = FALSE)[1])))
    if (length(v) && !is.na(v)) return(v)
  }
  NA_real_
}
