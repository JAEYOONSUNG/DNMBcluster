#' Cut Hierarchical Orthogroups (HOGs) from a rooted species tree
#'
#' For each internal node N of the rooted species tree, emits a HOG
#' table: for every OG, restrict the member list to genes from genomes
#' whose genome_key is in N's descendant leaf set. HOGs at the root
#' recapitulate the flat OG table; HOGs at lower internal nodes give
#' lineage-specific orthogroups that are robust to paralog inflation.
#'
#' Output layout under `out_dir`:
#' \preformatted{
#'   HOGs/
#'     N0_root.tsv        # root — full OGs
#'     N1_<clade>.tsv     # each internal node, leaf set in file header
#'     ...
#' }
#'
#' @param rooted_species_tree `ape::phylo` rooted tree (from `stride_root()`).
#' @param dnmb `load_dnmb()` result.
#' @param out_dir Directory under which `HOGs/` is written.
#' @param min_members Drop HOG rows with fewer than this many members.
#' @return Tibble with columns `node_id, clade_label, n_leaves,
#'   n_ogs_retained, n_members, path`.
#' @export
cut_hogs <- function(rooted_species_tree, dnmb, out_dir, min_members = 2L) {
  if (!requireNamespace("ape", quietly = TRUE)) stop("ape required.")
  dir.create(file.path(out_dir, "HOGs"), recursive = TRUE, showWarnings = FALSE)

  tree <- rooted_species_tree
  n_tip <- length(tree$tip.label)
  n_int <- tree$Nnode

  # key→uid map so we can filter by genome_uid in $clusters
  key2uid <- setNames(dnmb$genome_meta$genome_uid, dnmb$genome_meta$genome_key)

  node_ids <- seq_len(n_tip + n_int)
  internal <- node_ids[node_ids > n_tip]
  # Include root as node_id = n_tip + 1
  rows <- list()

  for (nd in internal) {
    desc <- .descendants(tree, nd)
    keys <- tree$tip.label[desc]
    uids <- unname(key2uid[keys])
    uids <- uids[!is.na(uids)]
    if (!length(uids)) next

    # Subset clusters to members whose genome is in this clade, then
    # drop clusters that fall below min_members within the clade.
    sub <- dnmb$clusters[dnmb$clusters$genome_uid %in% uids, ]
    counts <- dplyr::count(sub, cluster_id, name = "n_clade_members")
    keep_cids <- counts$cluster_id[counts$n_clade_members >= min_members]
    sub <- sub[sub$cluster_id %in% keep_cids, ]
    if (!nrow(sub)) next

    clade_label <- if (nd == n_tip + 1L) "root"
                   else paste(head(keys, 3), collapse = "+")
    safe_label <- gsub("[^A-Za-z0-9._-]+", "_", clade_label)
    path <- file.path(out_dir, "HOGs", sprintf("N%03d_%s.tsv", nd - n_tip, safe_label))

    con <- file(path, "w")
    cat("# node_id=", nd, "\n", sep = "", file = con)
    cat("# clade_leaves=", paste(keys, collapse = ","), "\n", sep = "", file = con)
    cat("cluster_id\tn_members\tmember_protein_uids\tmember_genome_uids\n", file = con)
    by_og <- split(sub, sub$cluster_id)
    for (og in by_og) {
      cat(og$cluster_id[1], "\t",
          nrow(og), "\t",
          paste(og$protein_uid, collapse = ","), "\t",
          paste(og$genome_uid,  collapse = ","), "\n",
          sep = "", file = con)
    }
    close(con)

    rows[[length(rows) + 1]] <- tibble::tibble(
      node_id        = nd - n_tip,
      clade_label    = clade_label,
      n_leaves       = length(keys),
      n_ogs_retained = length(by_og),
      n_members      = nrow(sub),
      path           = path
    )
  }
  dplyr::bind_rows(rows)
}


.descendants <- function(tr, node) {
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
  sort(unique(tips))
}
