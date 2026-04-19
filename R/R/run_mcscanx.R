#' Run MCScanX on an exported OG bundle and parse collinear blocks
#'
#' Takes an `export/mcscanx` directory produced by `export_og_bundles()`
#' and invokes the external `MCScanX` binary. Returns the parsed
#' `.collinearity` file as a tidy tibble of synteny blocks plus their
#' anchor pairs. Skips gracefully with a warning when the binary is
#' missing so batch pipelines keep running.
#'
#' @param bundle_dir Path to the `mcscanx` export directory. Must
#'   contain `genomes.gff` and `genomes.blast` (produced by
#'   `export_og_bundles(..., tools = "mcscanx")`).
#' @param prefix Input prefix MCScanX reads. Defaults to `"genomes"`,
#'   matching the exporter.
#' @param match_score,gap_penalty,e_value,max_gaps,match_size MCScanX
#'   tuning flags. Defaults mirror OrthoFinder's synteny plug-in:
#'   match_score=50, gap_penalty=-1, e_value=1e-5, max_gaps=25,
#'   match_size=5.
#' @param binary Explicit path to `MCScanX`. If NULL, uses `Sys.which`.
#' @return Named list: `blocks` tibble (`block_id, score, e_value,
#'   n_pairs, sp_a, sp_b, orientation`), `pairs` tibble (`block_id,
#'   gene_a, gene_b`), `collinearity_path`. NULL when binary not on
#'   PATH (with a `warning()`).
#' @export
run_mcscanx <- function(bundle_dir,
                        prefix       = "genomes",
                        match_score  = 50L,
                        gap_penalty  = -1L,
                        e_value      = 1e-5,
                        max_gaps     = 25L,
                        match_size   = 5L,
                        binary       = NULL) {
  bin <- if (!is.null(binary) && nzchar(binary)) binary else unname(Sys.which("MCScanX"))
  if (!nzchar(bin)) {
    warning("[run_mcscanx] 'MCScanX' binary not found on PATH; skipping.")
    return(NULL)
  }
  stopifnot(dir.exists(bundle_dir))
  gff   <- file.path(bundle_dir, paste0(prefix, ".gff"))
  blast <- file.path(bundle_dir, paste0(prefix, ".blast"))
  if (!file.exists(gff) || !file.exists(blast)) {
    stop("MCScanX inputs missing in ", bundle_dir,
         " (expected ", gff, " and ", blast, ")")
  }

  args <- c(file.path(bundle_dir, prefix),
            "-k", as.character(match_score),
            "-g", as.character(gap_penalty),
            "-e", format(e_value, scientific = TRUE),
            "-m", as.character(max_gaps),
            "-s", as.character(match_size))
  rc <- system2(bin, args = args, stdout = FALSE, stderr = FALSE)
  if (rc != 0) stop("MCScanX failed (exit=", rc, ")")
  coll <- file.path(bundle_dir, paste0(prefix, ".collinearity"))
  if (!file.exists(coll)) stop("MCScanX produced no .collinearity file: ", coll)

  parse_collinearity(coll)
}


#' Parse an MCScanX `.collinearity` file into tibbles
#'
#' Separated out so it can also be called on `.collinearity` outputs
#' produced elsewhere (e.g., by a prior MCScanX run on a shared cluster).
#'
#' @param path Path to a `.collinearity` file.
#' @return Named list with `blocks`, `pairs`, `collinearity_path`.
#' @export
parse_collinearity <- function(path) {
  stopifnot(file.exists(path))
  txt <- readLines(path, warn = FALSE)
  header_re <- "^##\\s*Alignment\\s+(\\d+):\\s*score=([0-9.eE+-]+)\\s+e_value=([0-9.eE+-]+)\\s+N=(\\d+)\\s+(\\S+)&(\\S+)\\s+(plus|minus)"

  hdr_idx <- grep("^## Alignment", txt)
  blk_rows <- vector("list", length(hdr_idx))
  pair_rows <- list()
  for (i in seq_along(hdr_idx)) {
    line <- txt[hdr_idx[i]]
    m <- regmatches(line, regexec(header_re, line))[[1]]
    if (length(m) != 8L) next
    bid <- as.integer(m[2])
    blk_rows[[i]] <- tibble::tibble(
      block_id    = bid,
      score       = suppressWarnings(as.numeric(m[3])),
      e_value     = suppressWarnings(as.numeric(m[4])),
      n_pairs     = suppressWarnings(as.integer(m[5])),
      sp_a        = m[6],
      sp_b        = m[7],
      orientation = m[8]
    )
    start <- hdr_idx[i] + 1L
    end   <- if (i < length(hdr_idx)) hdr_idx[i + 1L] - 1L else length(txt)
    body  <- txt[start:end]
    body  <- body[nzchar(body) & !startsWith(body, "#")]
    if (!length(body)) next
    # Pair lines: "<blk>-<rank>:  <gene_a>  <gene_b>  <evalue>"
    # Pair lines look like "  0-  0:  gene_a  gene_b  evalue". After
    # trimming leading whitespace, splitting on \s+ yields tokens
    # ["0-", "0:", gene_a, gene_b, evalue] — so genes sit at 3/4.
    body  <- trimws(body)
    parts <- strsplit(body, "\\s+")
    ok <- vapply(parts, function(p) length(p) >= 5L, logical(1))
    parts <- parts[ok]
    if (!length(parts)) next
    gene_a <- vapply(parts, `[`, character(1), 3L)
    gene_b <- vapply(parts, `[`, character(1), 4L)
    pair_rows[[length(pair_rows) + 1L]] <- tibble::tibble(
      block_id = bid,
      gene_a   = gene_a,
      gene_b   = gene_b
    )
  }

  blocks <- dplyr::bind_rows(blk_rows)
  pairs  <- if (length(pair_rows)) dplyr::bind_rows(pair_rows) else
    tibble::tibble(block_id = integer(), gene_a = character(), gene_b = character())
  list(blocks = blocks, pairs = pairs, collinearity_path = path)
}
