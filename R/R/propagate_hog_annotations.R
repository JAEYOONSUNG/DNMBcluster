#' Propagate per-protein annotations to per-HOG consensus
#'
#' Given a per-protein annotation table and the HOG tables written by
#' `cut_hogs()`, compute a consensus annotation for each HOG at each
#' internal species-tree node. "Consensus" is a plurality vote with
#' confidence = (top-count / total-annotated-members).
#'
#' Typical use: after `cut_hogs()`, feed `dnmb$gene_table` (joined with
#' whatever annotation source you have — COG, KO, GO, Pfam, product —
#' on `protein_uid`). Columns are propagated independently so you can
#' pass multiple annotation schemes at once.
#'
#' @param hogs Tibble returned by `cut_hogs()`.
#' @param annotations Tibble with a `protein_uid` key and one or more
#'   annotation columns. Multi-value columns should be split upstream
#'   (one row per protein_uid per annotation).
#' @param ann_cols Character vector of column names in `annotations` to
#'   propagate. Defaults to all columns except `protein_uid`.
#' @param min_confidence Drop HOG-level calls whose plurality share is
#'   below this threshold (set to 0 to always emit). Default 0.5.
#' @param out_path If non-NULL, write the consensus tibble to this TSV.
#' @return Tibble with columns `node_id, hog_id, cluster_id,
#'   annotation_column, value, count, total, confidence`.
#' @export
propagate_hog_annotations <- function(hogs,
                                      annotations,
                                      ann_cols = NULL,
                                      min_confidence = 0.5,
                                      out_path = NULL) {
  stopifnot("protein_uid" %in% names(annotations))
  if (is.null(ann_cols)) {
    ann_cols <- setdiff(names(annotations), "protein_uid")
  }
  if (!length(ann_cols)) stop("propagate_hog_annotations: no annotation columns.")
  if (!nrow(hogs)) {
    return(tibble::tibble(
      node_id = integer(), hog_id = character(), cluster_id = integer(),
      annotation_column = character(), value = character(),
      count = integer(), total = integer(), confidence = numeric()
    ))
  }

  # Read each per-node HOG table once, resolve leaves → protein_uid.
  # Leaves use `p<protein_uid>_g<genome_uid>` scheme; the protein_uid
  # is everything between "p" and "_g".
  all_rows <- list()
  for (k in seq_len(nrow(hogs))) {
    node_id <- hogs$node_id[k]
    path    <- hogs$path[k]
    if (!file.exists(path)) next
    tbl <- utils::read.table(path, sep = "\t", header = TRUE, comment.char = "#",
                             stringsAsFactors = FALSE)
    if (!nrow(tbl)) next
    for (r in seq_len(nrow(tbl))) {
      leaves <- strsplit(tbl$member_leaves[r], ",", fixed = TRUE)[[1]]
      puid_chr <- sub("_g.*$", "", sub("^p", "", leaves))
      # Match on character to stay safe with integer64 / large protein_uid.
      ann_uid_chr <- as.character(annotations$protein_uid)
      hit <- ann_uid_chr %in% puid_chr
      if (!any(hit)) next
      sub_ann <- annotations[hit, , drop = FALSE]
      if (!nrow(sub_ann)) next

      for (col in ann_cols) {
        vals <- sub_ann[[col]]
        vals <- vals[!is.na(vals) & nzchar(as.character(vals))]
        if (!length(vals)) next
        counts <- sort(table(vals), decreasing = TRUE)
        top    <- names(counts)[1]
        total  <- length(vals)
        conf   <- unname(counts[1]) / total
        if (conf < min_confidence) next
        all_rows[[length(all_rows) + 1]] <- tibble::tibble(
          node_id           = node_id,
          hog_id            = tbl$hog_id[r],
          cluster_id        = tbl$cluster_id[r],
          annotation_column = col,
          value             = as.character(top),
          count             = as.integer(counts[1]),
          total             = as.integer(total),
          confidence        = conf
        )
      }
    }
  }
  result <- dplyr::bind_rows(all_rows)
  if (!is.null(out_path) && nrow(result)) {
    utils::write.table(result, out_path, sep = "\t",
                       quote = FALSE, row.names = FALSE)
  }
  result
}
