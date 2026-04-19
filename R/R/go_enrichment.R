#' GO enrichment on HOGs / gene sets via topGO or Fisher fallback
#'
#' Given a universe of protein→GO term mappings (typically from the
#' DNMB EggNOG-fast stage) and a target gene set (e.g. HOG members,
#' a clade-specific HOG node, or any cluster_id subset), tests each
#' GO term for over-representation.
#'
#' Two engines:
#' \itemize{
#'   \item `topGO` (preferred) — Fisher with `elim` decorrelation,
#'         ontology graph aware. Requires the Bioconductor `topGO`
#'         package.
#'   \item Fisher fallback — plain `stats::fisher.test` per term.
#'         No ontology decorrelation, but no Bioconductor dep.
#' }
#'
#' @param protein_go Tibble with columns `protein_uid`, `go_id`
#'   (long format, one row per term; multiple rows per protein).
#'   Usually comes from the EggNOG-fast output or a custom join.
#' @param target_uids Integer vector of `protein_uid`s in the target
#'   set (e.g. `dnmb$clusters$protein_uid` filtered by HOG membership).
#' @param universe_uids Integer vector of the background universe.
#'   If NULL, uses `unique(protein_go$protein_uid)`.
#' @param ontology GO ontology to test: `"BP"`, `"MF"`, or `"CC"`.
#'   When `engine = "topGO"`, required. Ignored for the Fisher
#'   fallback (all terms pooled).
#' @param engine `"auto"` picks topGO if installed, else Fisher.
#' @param p_cutoff Term-level p-value threshold for the returned table.
#' @param min_annotated Drop terms with fewer than this many annotated
#'   proteins in the universe (default 5).
#' @return Tibble with `go_id, term, n_target, n_universe, n_annotated,
#'   expected, p_value, fdr, engine`. Sorted by p_value ascending.
#' @export
go_enrichment <- function(protein_go,
                          target_uids,
                          universe_uids = NULL,
                          ontology = c("BP", "MF", "CC"),
                          engine = c("auto", "topGO", "fisher"),
                          p_cutoff = 0.05,
                          min_annotated = 5L) {
  ontology <- match.arg(ontology)
  engine   <- match.arg(engine)
  stopifnot(all(c("protein_uid", "go_id") %in% names(protein_go)))

  if (is.null(universe_uids)) universe_uids <- unique(protein_go$protein_uid)
  target_uids <- intersect(unique(target_uids), universe_uids)
  if (!length(target_uids))
    return(.go_empty_result())

  use_topgo <- engine == "topGO" ||
    (engine == "auto" && requireNamespace("topGO", quietly = TRUE))

  if (use_topgo) {
    res <- tryCatch(.go_topgo(protein_go, target_uids, universe_uids,
                               ontology, min_annotated),
                    error = function(e) {
                      message("[go_enrichment] topGO failed (", conditionMessage(e),
                              "), falling back to Fisher.")
                      NULL
                    })
    if (!is.null(res)) return(.go_filter(res, p_cutoff))
  }
  .go_filter(.go_fisher(protein_go, target_uids, universe_uids,
                        min_annotated), p_cutoff)
}

.go_empty_result <- function() {
  tibble::tibble(
    go_id        = character(0),
    term         = character(0),
    n_target     = integer(0),
    n_universe   = integer(0),
    n_annotated  = integer(0),
    expected     = numeric(0),
    p_value      = numeric(0),
    fdr          = numeric(0),
    engine       = character(0)
  )
}

.go_filter <- function(tbl, p_cutoff) {
  if (!nrow(tbl)) return(tbl)
  tbl$fdr <- stats::p.adjust(tbl$p_value, method = "BH")
  tbl <- tbl[tbl$p_value <= p_cutoff, , drop = FALSE]
  tbl[order(tbl$p_value), , drop = FALSE]
}

.go_fisher <- function(pg, target, universe, min_ann) {
  pg <- pg[pg$protein_uid %in% universe, , drop = FALSE]
  n_universe <- length(universe)
  n_target   <- length(target)

  by_term <- split(pg$protein_uid, pg$go_id)
  rows <- list()
  for (gid in names(by_term)) {
    ann <- unique(by_term[[gid]])
    if (length(ann) < min_ann) next
    hit <- sum(target %in% ann)
    if (hit == 0L) next
    mat <- matrix(
      c(hit,
        length(ann) - hit,
        n_target - hit,
        n_universe - n_target - (length(ann) - hit)),
      nrow = 2
    )
    mat[mat < 0] <- 0
    p <- tryCatch(stats::fisher.test(mat, alternative = "greater")$p.value,
                  error = function(e) NA_real_)
    rows[[length(rows) + 1]] <- tibble::tibble(
      go_id        = gid,
      term         = NA_character_,
      n_target     = hit,
      n_universe   = n_target,
      n_annotated  = length(ann),
      expected     = n_target * length(ann) / n_universe,
      p_value      = p,
      fdr          = NA_real_,
      engine       = "fisher"
    )
  }
  if (!length(rows)) return(.go_empty_result())
  dplyr::bind_rows(rows)
}

.go_topgo <- function(pg, target, universe, ontology, min_ann) {
  pg <- pg[pg$protein_uid %in% universe, , drop = FALSE]
  uid2go <- split(pg$go_id, pg$protein_uid)
  uid2go <- lapply(uid2go, unique)

  all_uids <- as.character(universe)
  gene_list <- factor(as.integer(all_uids %in% as.character(target)),
                      levels = c(0, 1))
  names(gene_list) <- all_uids

  # topGOdata requires annotation in list form keyed by gene ID
  godata <- methods::new(
    "topGOdata",
    ontology      = ontology,
    allGenes      = gene_list,
    annot         = topGO::annFUN.gene2GO,
    gene2GO       = stats::setNames(uid2go, as.character(names(uid2go))),
    nodeSize      = min_ann
  )
  test  <- methods::new("elimCount", testStatistic = topGO::GOFisherTest,
                        name = "Fisher elim")
  result <- topGO::getSigGroups(godata, test)
  scores <- topGO::score(result)
  tab <- topGO::GenTable(godata, elimFisher = result,
                         orderBy = "elimFisher",
                         topNodes = length(scores))
  tibble::tibble(
    go_id        = tab$GO.ID,
    term         = tab$Term,
    n_target     = suppressWarnings(as.integer(tab$Significant)),
    n_universe   = length(target),
    n_annotated  = suppressWarnings(as.integer(tab$Annotated)),
    expected     = suppressWarnings(as.numeric(tab$Expected)),
    p_value      = suppressWarnings(as.numeric(tab$elimFisher)),
    fdr          = NA_real_,
    engine       = "topGO"
  )
}
