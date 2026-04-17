#' Label leaf pairs in gene trees as orthologs, paralogs, or xenologs
#'
#' For each gene tree, walks every internal node and classifies its
#' event using the species-overlap heuristic against a rooted species
#' tree:
#'
#' \itemize{
#'   \item \strong{Speciation} — child subtrees have disjoint species sets.
#'         All leaf pairs spanning the children are labeled \emph{orthologs}.
#'   \item \strong{Duplication} — child subtrees share >=1 species.
#'         All leaf pairs spanning the children are labeled \emph{paralogs}.
#'   \item \strong{Candidate xenolog} — the gene-tree branch length to the
#'         split is an outlier (> `xeno_z` stdev above the mean) compared
#'         to the species-tree distance between the same species pair.
#'         These are flagged separately in an \code{is_xenolog_candidate}
#'         column.
#' }
#'
#' Output is written as `relationships.parquet` if arrow is available,
#' else as `relationships.tsv`.
#'
#' @param og_result Tibble from `per_og_trees()`.
#' @param rooted_species_tree Output of `stride_root()`.
#' @param dnmb `load_dnmb()` result.
#' @param out_dir Directory to write `relationships.{parquet,tsv}`.
#' @param xeno_z Z-score cutoff for xenolog-candidate flag.
#' @return Tibble with columns
#'   `cluster_id, gene_a, gene_b, species_a, species_b, event,
#'    is_xenolog_candidate`.
#' @export
reconcile_relationships <- function(og_result,
                                    rooted_species_tree,
                                    dnmb,
                                    out_dir,
                                    xeno_z = 3.0) {
  if (!requireNamespace("ape", quietly = TRUE)) stop("ape required.")
  ok <- og_result[!is.na(og_result$tree_path), , drop = FALSE]
  if (!nrow(ok)) stop("reconcile_relationships: no trees to reconcile.")

  uid2key <- setNames(dnmb$genome_meta$genome_key, as.character(dnmb$genome_meta$genome_uid))
  sp_cophen <- ape::cophenetic.phylo(rooted_species_tree)

  all_rows <- vector("list", nrow(ok))
  for (i in seq_len(nrow(ok))) {
    tr <- ape::read.tree(ok$tree_path[i])
    if (!is.null(tr$edge.length)) tr$edge.length <- pmax(tr$edge.length, 0)
    rows <- .classify_tree_events(tr, ok$cluster_id[i], uid2key, sp_cophen, xeno_z)
    all_rows[[i]] <- rows
  }
  df <- dplyr::bind_rows(all_rows)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("arrow", quietly = TRUE)) {
    parq <- file.path(out_dir, "relationships.parquet")
    arrow::write_parquet(df, parq)
    message("[reconcile] wrote ", parq, "  (", nrow(df), " pairs)")
  } else {
    tsv <- file.path(out_dir, "relationships.tsv")
    utils::write.table(df, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[reconcile] wrote ", tsv)
  }
  df
}


.classify_tree_events <- function(tr, cluster_id, uid2key, sp_cophen, xeno_z) {
  tips <- tr$tip.label
  n_tip <- length(tips)
  if (n_tip < 2) return(NULL)

  gidx <- sub(".*_g", "", tips)
  species_key <- unname(uid2key[gidx])

  gene_cophen <- ape::cophenetic.phylo(tr)

  rows <- list()
  n_int <- tr$Nnode
  for (nd in (n_tip + 1):(n_tip + n_int)) {
    children <- tr$edge[tr$edge[, 1] == nd, 2]
    if (length(children) != 2) next
    side <- lapply(children, function(ch) {
      if (ch <= n_tip) ch else .tips_under(tr, ch)
    })
    sp_a <- unique(species_key[side[[1]]])
    sp_b <- unique(species_key[side[[2]]])
    shared <- intersect(sp_a, sp_b)
    event <- if (length(shared) >= 1) "duplication" else "speciation"

    # Emit one row per cross-side leaf pair
    for (ia in side[[1]]) for (ib in side[[2]]) {
      sa <- species_key[ia]; sb <- species_key[ib]
      # Xenolog candidate: gene distance >> expected species distance
      expected <- if (!is.na(sa) && !is.na(sb) && sa %in% rownames(sp_cophen)
                       && sb %in% colnames(sp_cophen)) sp_cophen[sa, sb] else NA_real_
      observed <- gene_cophen[ia, ib]
      rows[[length(rows) + 1]] <- list(
        cluster_id           = cluster_id,
        gene_a               = tips[ia],
        gene_b               = tips[ib],
        species_a            = sa,
        species_b            = sb,
        event                = if (event == "speciation") "ortholog" else "paralog",
        observed_distance    = observed,
        expected_distance    = expected,
        is_xenolog_candidate = !is.na(expected) && expected > 0 &&
                                  (observed - expected) / max(expected, 1e-6) > xeno_z
      )
    }
  }
  dplyr::bind_rows(rows)
}

.tips_under <- function(tr, node) {
  n_tip <- length(tr$tip.label)
  if (node <= n_tip) return(node)
  stack <- node
  tips <- integer()
  while (length(stack)) {
    cur <- stack[1]; stack <- stack[-1]
    ch <- tr$edge[tr$edge[, 1] == cur, 2]
    leaf <- ch[ch <= n_tip]
    tips <- c(tips, leaf)
    stack <- c(stack, ch[ch > n_tip])
  }
  tips
}
