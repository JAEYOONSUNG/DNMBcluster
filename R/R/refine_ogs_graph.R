#' Refine orthogroups via length-normalized edge graph + community split
#'
#' DNMBcluster's upstream clustering is identity-greedy (MMseqs2/CD-HIT/
#' DIAMOND/USEARCH), which is fast but can over-merge paralog families
#' when member length varies (fusion genes, domain-shuffled families,
#' aberrant-length isoforms). OrthoFinder avoids this by scoring edges
#' with length-normalized bitscores before MCL clustering.
#'
#' This function post-processes an existing DNMB OG table with the same
#' idea, without rerunning upstream clustering:
#'
#' \enumerate{
#'   \item For each OG with \code{>= min_size} members, compute all-vs-all
#'         bitscores within that OG (DIAMOND if available, else
#'         \code{Biostrings::pairwiseAlignment} local BLOSUM62).
#'   \item Normalize each edge by the geometric mean of the two members'
#'         lengths: \code{score_norm = bitscore / sqrt(len_a * len_b)}.
#'   \item Threshold edges at \code{edge_quantile} of per-OG distribution
#'         (default 0.25 -- keep top 75%). This cuts weak cross-subfamily
#'         links before clustering.
#'   \item Community-detect the trimmed graph:
#'         \itemize{
#'           \item \code{"mcl"} -- the \code{MCL} CRAN package (if installed)
#'           \item \code{"louvain"} -- \code{igraph::cluster_louvain} (fallback)
#'           \item \code{"walktrap"} -- \code{igraph::cluster_walktrap}
#'         }
#'   \item If the community algorithm returns > 1 cluster the OG is
#'         \emph{split} into sub-OGs, each getting a new
#'         \code{refined_id = <orig>_<subidx>}.
#' }
#'
#' @param dnmb Result of \code{load_dnmb()}.
#' @param out_dir Destination directory; writes \code{refined_OGs.tsv}.
#' @param min_size Skip OGs smaller than this (can't meaningfully split).
#' @param method Community detection algorithm: mcl / louvain / walktrap.
#' @param inflation MCL inflation parameter (only used when \code{method = "mcl"}).
#' @param edge_quantile Quantile for edge trimming; 0 disables trimming.
#' @param threads Thread count for DIAMOND.
#' @return Tibble with columns \code{cluster_id, refined_id, protein_uid,
#'   genome_uid, n_subgroups}. Original OGs that weren't split get
#'   \code{refined_id == cluster_id} and \code{n_subgroups == 1}.
#' @export
refine_ogs_graph <- function(dnmb,
                             out_dir,
                             min_size = 5L,
                             method   = c("mcl", "louvain", "walktrap"),
                             inflation = 1.5,
                             edge_quantile = 0.25,
                             threads = 2L) {
  method <- match.arg(method)
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("refine_ogs_graph requires the 'igraph' package.")
  }
  have_diamond <- nzchar(Sys.which("diamond"))
  have_biostr  <- requireNamespace("Biostrings", quietly = TRUE)
  if (!have_diamond && !have_biostr) {
    stop("refine_ogs_graph needs DIAMOND on PATH or the Biostrings package.")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Join cluster_id <-> protein sequence
  members <- dnmb$clusters %>%
    dplyr::select(cluster_id, protein_uid, genome_uid) %>%
    dplyr::inner_join(
      dnmb$gene_table %>% dplyr::select(protein_uid, translation, length),
      by = "protein_uid"
    )
  by_og <- split(members, members$cluster_id)

  out_rows <- list()
  for (cid_chr in names(by_og)) {
    df <- by_og[[cid_chr]]
    cid <- as.integer(cid_chr)
    if (nrow(df) < min_size) {
      out_rows[[length(out_rows) + 1]] <- tibble::tibble(
        cluster_id   = cid,
        refined_id   = as.character(cid),
        protein_uid  = df$protein_uid,
        genome_uid   = df$genome_uid,
        n_subgroups  = 1L
      )
      next
    }

    edges <- .intra_og_scores(df, have_diamond, threads)
    if (!nrow(edges)) {
      out_rows[[length(out_rows) + 1]] <- tibble::tibble(
        cluster_id   = cid,
        refined_id   = as.character(cid),
        protein_uid  = df$protein_uid,
        genome_uid   = df$genome_uid,
        n_subgroups  = 1L
      )
      next
    }

    # Length-normalize
    len_map <- stats::setNames(df$length, as.character(df$protein_uid))
    edges$len_a <- unname(len_map[as.character(edges$a)])
    edges$len_b <- unname(len_map[as.character(edges$b)])
    edges$score_norm <- edges$bitscore / sqrt(edges$len_a * edges$len_b)

    # Trim weak edges
    if (edge_quantile > 0 && nrow(edges) >= 10) {
      thr <- stats::quantile(edges$score_norm, edge_quantile, na.rm = TRUE)
      edges <- edges[edges$score_norm >= thr, , drop = FALSE]
    }

    # Community-detect
    labels <- .run_community(edges, df$protein_uid, method, inflation)

    n_sub <- length(unique(labels))
    sub_ids <- if (n_sub == 1L) rep(as.character(cid), nrow(df)) else paste0(cid, "_", labels[as.character(df$protein_uid)])

    out_rows[[length(out_rows) + 1]] <- tibble::tibble(
      cluster_id   = cid,
      refined_id   = sub_ids,
      protein_uid  = df$protein_uid,
      genome_uid   = df$genome_uid,
      n_subgroups  = n_sub
    )
  }

  result <- dplyr::bind_rows(out_rows)
  path <- file.path(out_dir, "refined_OGs.tsv")
  utils::write.table(result, path, sep = "\t", quote = FALSE, row.names = FALSE)

  n_split <- sum(result$n_subgroups > 1) / pmax(nrow(result), 1)
  message(sprintf("[refine_ogs_graph] %d proteins -> %d refined OGs (%d original). Splits: %d OGs affected. -> %s",
                  nrow(result),
                  length(unique(result$refined_id)),
                  length(unique(result$cluster_id)),
                  length(unique(result$cluster_id[result$n_subgroups > 1])),
                  path))
  result
}


# ---- helpers ---------------------------------------------------------------

.intra_og_scores <- function(df, have_diamond, threads) {
  if (have_diamond) return(.intra_diamond(df, threads))
  .intra_biostrings(df)
}

.intra_diamond <- function(df, threads) {
  fa <- tempfile(fileext = ".fa")
  on.exit(unlink(fa), add = TRUE)
  con <- file(fa, "w")
  for (i in seq_len(nrow(df))) {
    cat(">", format(df$protein_uid[i], scientific = FALSE), "\n",
        df$translation[i], "\n", sep = "", file = con)
  }
  close(con)
  db <- tempfile(fileext = ".dmnd")
  out <- tempfile(fileext = ".m8")
  on.exit({ unlink(db); unlink(out) }, add = TRUE)
  res <- system2("diamond", c("makedb", "--in", shQuote(fa),
                                "-d", shQuote(sub("\\.dmnd$", "", db))),
                  stdout = FALSE, stderr = FALSE)
  if (res != 0) return(data.frame())
  res <- system2("diamond", c("blastp", "-q", shQuote(fa),
                                "-d", shQuote(sub("\\.dmnd$", "", db)),
                                "-o", shQuote(out), "-k", as.character(nrow(df)),
                                "--threads", threads,
                                "--outfmt", "6", "qseqid", "sseqid", "bitscore"),
                  stdout = FALSE, stderr = FALSE)
  if (res != 0 || !file.exists(out) || file.info(out)$size == 0) return(data.frame())
  hits <- utils::read.table(out, sep = "\t", header = FALSE,
                             stringsAsFactors = FALSE,
                             col.names = c("a", "b", "bitscore"))
  hits <- hits[hits$a != hits$b, , drop = FALSE]
  if (!nrow(hits)) return(hits)
  # DIAMOND emits both a->b and b->a with (typically unequal) bitscores.
  # Fold to an undirected edge list: canonicalize the pair and average
  # the two bitscores. Leaving duplicates in inflates every edge twice
  # and lets Louvain/MCL mistake directional asymmetry for community
  # structure.
  ord <- pmin(hits$a, hits$b) != hits$a
  tmp <- hits$a[ord]; hits$a[ord] <- hits$b[ord]; hits$b[ord] <- tmp
  stats::aggregate(bitscore ~ a + b, data = hits, FUN = mean)
}

.intra_biostrings <- function(df) {
  n <- nrow(df)
  if (n < 2) return(data.frame())
  ss <- Biostrings::AAStringSet(df$translation)
  names(ss) <- as.character(df$protein_uid)
  pw_align <- .dnmb_pairwise_fn()
  # Upper triangle pairwise — replicate the singleton subject to match
  # pattern length (Biostrings >= 2.77 / pwalign require equal lengths).
  rows <- list()
  for (i in 1:(n - 1)) {
    j <- (i + 1):n
    subj <- ss[rep(i, length(j))]
    scores <- suppressWarnings(pw_align(
      pattern = ss[j], subject = subj,
      substitutionMatrix = "BLOSUM62",
      scoreOnly = TRUE, type = "local"))
    rows[[i]] <- data.frame(a = names(ss)[i], b = names(ss)[j],
                             bitscore = scores, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

.run_community <- function(edges, all_uids, method, inflation) {
  uids <- as.character(all_uids)
  if (!nrow(edges)) return(stats::setNames(rep("1", length(uids)), uids))

  g <- igraph::graph_from_data_frame(
    edges[, c("a", "b")], directed = FALSE,
    vertices = data.frame(name = uids, stringsAsFactors = FALSE))
  igraph::E(g)$weight <- edges$score_norm

  if (method == "mcl" && requireNamespace("MCL", quietly = TRUE)) {
    adj <- as.matrix(igraph::as_adjacency_matrix(g, attr = "weight"))
    diag(adj) <- 0
    res <- MCL::mcl(adj, addLoops = TRUE, inflation = inflation, allow1 = TRUE)
    if (is.list(res) && !is.null(res$Cluster)) {
      return(stats::setNames(as.character(res$Cluster), uids))
    }
    method <- "louvain"  # fall through
  }
  if (method == "walktrap") {
    cl <- igraph::cluster_walktrap(g, weights = igraph::E(g)$weight)
  } else {
    cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  }
  stats::setNames(as.character(igraph::membership(cl)), uids)
}
