#' Export per-OG alignments + trees in tool-specific layouts
#'
#' Once `per_og_trees()` has written `orthogroups/<id>/{aln.fasta,tree.nwk}`,
#' downstream comparative-genomics tools each want the same files in a
#' slightly different shape. This helper copies/rewrites them once so a
#' user can hand the bundle to any of:
#'
#' \itemize{
#'   \item \strong{generax}: one directory per OG, `families.txt` index
#'   \item \strong{ale}: `.ale` files (tree posterior — here we emit the ML
#'         tree as a single-sample placeholder for downstream reconciliation)
#'   \item \strong{paml}: `*.fasta` codeml-style alignments + control stub
#'   \item \strong{mcscanx}: `.blast` pairs + `.gff` positions (requires
#'         \code{dnmb$id_map} with contig/start/end — already present)
#' }
#'
#' @param dnmb Result of `load_dnmb()`.
#' @param og_result Output of `per_og_trees()`.
#' @param out_dir Directory under which `export/<tool>/` is written.
#' @param tools Character vector from `c("generax","ale","paml","mcscanx")`.
#' @param species_tree Optional path to rooted species tree (required by
#'   generax/ale; if NULL the caller must add it before running the tool).
#' @return Named list of output paths (one entry per tool).
#' @export
export_og_bundles <- function(dnmb,
                              og_result,
                              out_dir,
                              tools = c("generax", "ale", "paml", "mcscanx"),
                              species_tree = NULL) {
  tools <- match.arg(tools, several.ok = TRUE)
  ok <- og_result[!is.na(og_result$tree_path), , drop = FALSE]
  if (!nrow(ok)) stop("export_og_bundles: og_result has no successful trees.")

  out <- list()
  if ("generax" %in% tools) out$generax <- .export_generax(dnmb, ok, out_dir, species_tree)
  if ("ale"     %in% tools) out$ale     <- .export_ale(ok, out_dir)
  if ("paml"    %in% tools) out$paml    <- .export_paml(ok, out_dir)
  if ("mcscanx" %in% tools) out$mcscanx <- .export_mcscanx(dnmb, ok, out_dir)
  out
}


.export_generax <- function(dnmb, ok, out_dir, species_tree) {
  root <- file.path(out_dir, "export", "generax")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  families_path <- file.path(root, "families.txt")
  con <- file(families_path, "w"); on.exit(close(con), add = TRUE)

  gid2name <- setNames(dnmb$genome_meta$genome_key, as.character(dnmb$genome_meta$genome_uid))

  for (i in seq_len(nrow(ok))) {
    cid <- ok$cluster_id[i]
    fam <- sprintf("OG_%07d", cid)
    fam_dir <- file.path(root, fam)
    dir.create(fam_dir, showWarnings = FALSE)
    file.copy(ok$aln_path[i],  file.path(fam_dir, "aln.fasta"),  overwrite = TRUE)
    file.copy(ok$tree_path[i], file.path(fam_dir, "tree.nwk"),   overwrite = TRUE)

    map_path <- file.path(fam_dir, "mapping.link")
    members <- dnmb$clusters[dnmb$clusters$cluster_id == cid, ]
    lines <- paste0("p", members$protein_uid, "_g", members$genome_uid,
                    ":", gid2name[as.character(members$genome_uid)])
    writeLines(lines, map_path)

    cat(sprintf("- FAMILY: %s\nalignment = %s\nmapping = %s\nstarting_gene_tree = %s\nsubst_model = LG+G4\n\n",
                fam, file.path(fam, "aln.fasta"),
                file.path(fam, "mapping.link"),
                file.path(fam, "tree.nwk")),
        file = con)
  }
  if (!is.null(species_tree)) {
    file.copy(species_tree, file.path(root, "species_tree.nwk"), overwrite = TRUE)
  }
  root
}


.export_ale <- function(ok, out_dir) {
  root <- file.path(out_dir, "export", "ale")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(ok))) {
    dest <- file.path(root, sprintf("OG_%07d.tree", ok$cluster_id[i]))
    file.copy(ok$tree_path[i], dest, overwrite = TRUE)
  }
  root
}


.export_paml <- function(ok, out_dir) {
  root <- file.path(out_dir, "export", "paml")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(ok))) {
    dest <- file.path(root, sprintf("OG_%07d.fasta", ok$cluster_id[i]))
    file.copy(ok$aln_path[i], dest, overwrite = TRUE)
  }
  writeLines(c(
    "seqfile  = OG_XXXXXXX.fasta",
    "treefile = OG_XXXXXXX.tree",
    "outfile  = OG_XXXXXXX.out",
    "noisy = 3",
    "verbose = 1",
    "seqtype = 2",   # AA
    "model = 2",     # Empirical+F
    "aaRatefile = lg.dat"
  ), file.path(root, "codeml.ctl.template"))
  root
}


.export_mcscanx <- function(dnmb, ok, out_dir) {
  root <- file.path(out_dir, "export", "mcscanx")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)

  # .gff: species_tag\tgene_id\tstart\tend
  gff_path <- file.path(root, "genomes.gff")
  gff <- dnmb$id_map %>%
    dplyr::transmute(
      sp    = paste0("g", genome_uid),
      id    = paste0("p", protein_uid, "_g", genome_uid),
      start = as.integer(start),
      end   = as.integer(end)
    )
  utils::write.table(gff[, c("sp", "id", "start", "end")],
                     gff_path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)

  # .blast: pairs within each OG (all-vs-all minus self)
  blast_path <- file.path(root, "genomes.blast")
  con <- file(blast_path, "w"); on.exit(close(con), add = TRUE)
  for (i in seq_len(nrow(ok))) {
    cid <- ok$cluster_id[i]
    mem <- dnmb$clusters[dnmb$clusters$cluster_id == cid, ]
    ids <- paste0("p", mem$protein_uid, "_g", mem$genome_uid)
    if (length(ids) < 2) next
    pairs <- utils::combn(ids, 2)
    # crude placeholder: bitscore=100, evalue=1e-50
    lines <- sprintf("%s\t%s\t100.0\t200\t0\t0\t1\t200\t1\t200\t1e-50\t100",
                     pairs[1,], pairs[2,])
    cat(paste(lines, collapse = "\n"), "\n", file = con)
  }
  root
}
