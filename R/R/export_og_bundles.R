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
#'   \item \strong{paml}: per-OG directories `OG_XXXXXXX/{aln.aa.fasta,
#'         tree.nwk}` plus a codon-ready `codeml.ctl.template`. DNMBcluster
#'         only has amino-acid sequences, so the emitted alignment is AA
#'         and is NOT directly consumable by `run_codeml_m0()` (which
#'         requires `aln.nuc.phy`). The README explains the required
#'         pal2nal/TranslatorX back-translation step.
#'   \item \strong{mcscanx}: `.blast` pairs + `.gff` positions using real
#'         contig IDs and genomic coordinates from \code{dnmb$id_map}.
#'   \item \strong{hyphy}: per-OG directories holding `aln.aa.fasta` +
#'         `tree.nwk` plus a README with the one-line `hyphy fel/busted/
#'         absrel` invocations. Like the PAML bundle this leaves codon
#'         back-translation to the caller; `back_translate_og_codons()`
#'         writes the `aln.nuc.fasta` HyPhy actually consumes.
#' }
#'
#' @param dnmb Result of `load_dnmb()`.
#' @param og_result Output of `per_og_trees()`.
#' @param out_dir Directory under which `export/<tool>/` is written.
#' @param tools Character vector from `c("generax","ale","paml","hyphy","mcscanx","dlcpar")`.
#' @param species_tree Optional path to rooted species tree (required by
#'   generax / ale / dlcpar; if NULL the caller must add it before running
#'   the tool).
#' @param aa_model Amino-acid substitution model used when building the
#'   per-OG alignments/trees upstream. Written into the GeneRax family
#'   file as `subst_model = <aa_model>+G4` so reconciliation scores the
#'   same model the trees were built under. Default `"LG"`.
#' @return Named list of output paths (one entry per tool).
#' @export
export_og_bundles <- function(dnmb,
                              og_result,
                              out_dir,
                              tools = c("generax", "ale", "paml", "hyphy",
                                        "mcscanx", "dlcpar"),
                              species_tree = NULL,
                              aa_model = "LG") {
  tools <- match.arg(tools, several.ok = TRUE)
  ok <- og_result[!is.na(og_result$tree_path), , drop = FALSE]
  if (!nrow(ok)) stop("export_og_bundles: og_result has no successful trees.")

  out <- list()
  if ("generax" %in% tools) out$generax <- .export_generax(dnmb, ok, out_dir,
                                                            species_tree,
                                                            aa_model = aa_model)
  if ("ale"     %in% tools) out$ale     <- .export_ale(ok, out_dir)
  if ("paml"    %in% tools) out$paml    <- .export_paml(ok, out_dir)
  if ("hyphy"   %in% tools) out$hyphy   <- .export_hyphy(ok, out_dir)
  if ("mcscanx" %in% tools) out$mcscanx <- .export_mcscanx(dnmb, ok, out_dir)
  if ("dlcpar"  %in% tools) out$dlcpar  <- .export_dlcpar(dnmb, ok, out_dir, species_tree)
  out
}


.export_generax <- function(dnmb, ok, out_dir, species_tree, aa_model = "LG") {
  # GeneRax's subst_model accepts standard RAxML-style names with a
  # gamma suffix, e.g. LG+G4, WAG+G4. Accept plain AA model names from
  # upstream and append "+G4" if the caller did not already include a
  # rate-heterogeneity suffix.
  subst_model <- if (grepl("\\+", aa_model)) aa_model
                 else paste0(aa_model, "+G4")
  root <- file.path(out_dir, "export", "generax")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  families_path <- file.path(root, "families.txt")
  con <- file(families_path, "w"); on.exit(close(con), add = TRUE)
  # GeneRax family file format requires the [FAMILIES] header and
  # `- <family_name>` on its own line (not `- FAMILY: <name>`, which the
  # parser silently drops as an invalid family).
  cat("[FAMILIES]\n", file = con)

  gid2name <- stats::setNames(dnmb$genome_meta$genome_key, as.character(dnmb$genome_meta$genome_uid))

  for (i in seq_len(nrow(ok))) {
    cid <- ok$cluster_id[i]
    fam <- sprintf("OG_%07d", cid)
    fam_dir <- file.path(root, fam)
    dir.create(fam_dir, showWarnings = FALSE)
    file.copy(ok$aln_path[i],  file.path(fam_dir, "aln.fasta"),  overwrite = TRUE)
    file.copy(ok$tree_path[i], file.path(fam_dir, "tree.nwk"),   overwrite = TRUE)

    map_path <- file.path(fam_dir, "mapping.link")
    # Read membership straight from the alignment FASTA so we match
    # whatever cluster_id space is in play (original vs. refined). The
    # header format written by per_og_trees is "p<uid>_g<genome_uid>".
    hdrs <- grep("^>", readLines(ok$aln_path[i], warn = FALSE), value = TRUE)
    gene_ids <- sub("^>", "", hdrs)
    gids <- sub(".*_g", "", gene_ids)
    sp_names <- gid2name[gids]
    keep <- !is.na(sp_names) & nzchar(sp_names)
    gene_ids <- gene_ids[keep]; sp_names <- sp_names[keep]
    # GeneRax "link" format: one line per species as
    #   speciesID:gene1;gene2;...
    # (not gene:species — that shape makes GeneRax reject the family.)
    by_sp <- split(gene_ids, sp_names)
    lines <- vapply(names(by_sp),
                    function(sp) paste0(sp, ":", paste(by_sp[[sp]], collapse = ";")),
                    character(1))
    writeLines(unname(lines), map_path)

    # Emit absolute paths — GeneRax resolves entries relative to its
    # working directory, not the families.txt parent, so relative paths
    # produce "alignment file does not exist" and the family is dropped.
    cat(sprintf("- %s\nalignment = %s\nmapping = %s\nstarting_gene_tree = %s\nsubst_model = %s\n\n",
                fam,
                normalizePath(file.path(fam_dir, "aln.fasta"),      mustWork = TRUE),
                normalizePath(file.path(fam_dir, "mapping.link"),   mustWork = TRUE),
                normalizePath(file.path(fam_dir, "tree.nwk"),       mustWork = TRUE),
                subst_model),
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

  # Emit per-OG subdirs holding both the AA alignment and the gene tree,
  # matching the layout `run_codeml_m0()` expects via `og_dirs`. DNMB has
  # no nucleotide sequences, so we cannot write `aln.nuc.phy` directly —
  # the README documents the back-translation step.
  for (i in seq_len(nrow(ok))) {
    cid <- ok$cluster_id[i]
    fam_dir <- file.path(root, sprintf("OG_%07d", cid))
    dir.create(fam_dir, showWarnings = FALSE, recursive = TRUE)
    file.copy(ok$aln_path[i],  file.path(fam_dir, "aln.aa.fasta"),
              overwrite = TRUE)
    file.copy(ok$tree_path[i], file.path(fam_dir, "tree.nwk"),
              overwrite = TRUE)
  }
  writeLines(c(
    "# Codon-level codeml control (site model M0) — expects the codon",
    "# alignment at aln.nuc.phy next to tree.nwk. Produce aln.nuc.phy",
    "# by back-translating aln.aa.fasta with pal2nal/TranslatorX using",
    "# the original CDS nucleotide FASTAs; DNMBcluster does not keep",
    "# nucleotide sequences in the main pipeline.",
    "seqfile  = aln.nuc.phy",
    "treefile = tree.nwk",
    "outfile  = out.txt",
    "noisy    = 0",
    "verbose  = 0",
    "runmode  = 0",
    "seqtype  = 1",   # codons
    "CodonFreq = 2",  # F3x4
    "model    = 0",   # one omega across branches (M0)
    "NSsites  = 0",
    "icode    = 0",
    "fix_kappa = 0", "kappa = 2",
    "fix_omega = 0", "omega = 0.4",
    "cleandata = 1"
  ), file.path(root, "codeml.ctl.template"))
  writeLines(c(
    "DNMBcluster PAML bundle",
    "========================",
    "",
    "Each OG_XXXXXXX/ directory contains:",
    "  aln.aa.fasta   — amino-acid MSA (from per_og_trees())",
    "  tree.nwk       — gene tree (ML/NJ, rooted by DNMBcluster)",
    "",
    "To run codeml M0 (run_codeml_m0()):",
    "  1. Back-translate aln.aa.fasta -> aln.nuc.phy with the original",
    "     CDS nucleotide FASTAs (pal2nal.pl -output paml ...). DNMB does",
    "     NOT ship nucleotide sequences, so you must supply them.",
    "  2. Copy codeml.ctl.template into each OG dir as codeml.ctl.",
    "  3. run_codeml_m0(og_dirs = list.dirs(\".\", recursive=FALSE),",
    "                   aln_name = \"aln.nuc.phy\",",
    "                   tree_name = \"tree.nwk\",",
    "                   work_dir = \"paml_work\")",
    "",
    "run_codeml_m0 validates codon alphabet + length-%%-3 and writes a",
    "skip_reason.txt when an input fails, so AA-only families are flagged",
    "instead of silently producing NA omega."
  ), file.path(root, "README.txt"))
  root
}


.export_hyphy <- function(ok, out_dir) {
  root <- file.path(out_dir, "export", "hyphy")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  # HyPhy's fel/busted/absrel drivers all take `--alignment <codon.fa>
  # --tree <tree.nwk> --output <out.json>`. We stage the AA alignment
  # and gene tree; the caller back-translates aln.aa.fasta into
  # aln.nuc.fasta before invoking hyphy. This mirrors the PAML bundle —
  # DNMB has no CDS sequences, so we cannot emit the codon FASTA here.
  for (i in seq_len(nrow(ok))) {
    cid <- ok$cluster_id[i]
    fam_dir <- file.path(root, sprintf("OG_%07d", cid))
    dir.create(fam_dir, showWarnings = FALSE, recursive = TRUE)
    file.copy(ok$aln_path[i],  file.path(fam_dir, "aln.aa.fasta"),
              overwrite = TRUE)
    file.copy(ok$tree_path[i], file.path(fam_dir, "tree.nwk"),
              overwrite = TRUE)
  }
  writeLines(c(
    "DNMBcluster HyPhy bundle",
    "=========================",
    "",
    "Each OG_XXXXXXX/ directory contains:",
    "  aln.aa.fasta   — amino-acid MSA (from per_og_trees())",
    "  tree.nwk       — rooted gene tree",
    "",
    "Prepare the codon alignment HyPhy expects:",
    "  # inside R:",
    "  back_translate_og_codons(og_result, cds_source = \"<cds.fasta-or-dir>\")",
    "  # → writes aln.nuc.fasta next to every aln.aa.fasta.",
    "",
    "Run HyPhy (binary from bioconda `hyphy`, 2.5.52+):",
    "  hyphy fel    --alignment aln.nuc.fasta --tree tree.nwk --output fel.json",
    "  hyphy busted --alignment aln.nuc.fasta --tree tree.nwk --output busted.json \\",
    "               --srv Yes        # BUSTED[S] — recommended for bacteria",
    "  hyphy absrel --alignment aln.nuc.fasta --tree tree.nwk --output absrel.json",
    "",
    "Or batch via R:",
    "  og_dirs <- list.dirs(\".\", recursive = FALSE)",
    "  run_hyphy_fel   (og_dirs, work_dir = \"hyphy_work\")",
    "  run_hyphy_busted(og_dirs, work_dir = \"hyphy_work\", srv = TRUE)",
    "  run_hyphy_absrel(og_dirs, work_dir = \"hyphy_work\")",
    "",
    "Reviewer notes (2025 bacterial pan-genome papers):",
    "  • HyPhy FEL/BUSTED has replaced codeml M7/M8 and branch-site Model A",
    "    as the default selection test; codeml M0 is kept as a fast baseline.",
    "  • BUSTED[S] (srv = TRUE) corrects for synonymous-rate variation and",
    "    is the recommended mode for closely-related bacterial alignments.",
    "  • Wolf et al. 2021 PNAS: bacterial dN/dS is inflated at short",
    "    branches. Consider HyPhy MSS or filtering to sufficiently",
    "    diverged pairs before interpreting omega as selection strength."
  ), file.path(root, "README.txt"))
  root
}


.export_dlcpar <- function(dnmb, ok, out_dir, species_tree) {
  root <- file.path(out_dir, "export", "dlcpar")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  gid2name <- stats::setNames(dnmb$genome_meta$genome_key, as.character(dnmb$genome_meta$genome_uid))

  # DLCpar expects: (1) a species tree .stree, (2) per-family gene tree
  # with leaf names of the form "species__geneid", (3) a .smap file
  # mapping species regex → species name. We emit all three.
  if (!is.null(species_tree) && file.exists(species_tree)) {
    file.copy(species_tree, file.path(root, "species.stree"), overwrite = TRUE)
  }

  smap_path <- file.path(root, "genes.smap")
  writeLines(sprintf("*_g%s\t%s", names(gid2name), gid2name), smap_path)

  for (i in seq_len(nrow(ok))) {
    cid <- ok$cluster_id[i]
    tr <- ape::read.tree(ok$tree_path[i])
    if (!is.null(tr$edge.length)) tr$edge.length <- pmax(tr$edge.length, 0)
    # Rewrite leaves from "p<uid>_g<gid>" → "<species>__p<uid>_g<gid>"
    tr$tip.label <- vapply(tr$tip.label, function(lab) {
      gid <- sub(".*_g", "", lab)
      paste0(gid2name[gid], "__", lab)
    }, character(1))
    dest <- file.path(root, sprintf("OG_%07d.tree", cid))
    ape::write.tree(tr, dest)
  }
  writeLines(c(
    "# Example DLCpar invocation:",
    "#   dlcpar search -s species.stree -S genes.smap -I .tree -O .dlcpar OG_*.tree",
    "# Leaf format used here: <species_key>__p<uid>_g<genome_uid>",
    "# smap regex matches \"*_g<genome_uid>\" for each input genome."
  ), file.path(root, "README.txt"))
  root
}

.export_mcscanx <- function(dnmb, ok, out_dir) {
  root <- file.path(out_dir, "export", "mcscanx")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)

  # .gff columns for MCScanX: chromosome_id\tgene_id\tstart\tend.
  # MCScanX derives collinear blocks from adjacency on each chromosome,
  # so the real contig/replicon ID must go in column 1 — collapsing all
  # contigs onto "g<genome_uid>" would fuse unrelated contigs into one
  # pseudo-chromosome and invalidate every block. We prefix the contig
  # with a short per-genome tag so MCScanX's species-code parsing still
  # works across genomes while keeping contigs distinct within a genome.
  gid_tag <- stats::setNames(
    sprintf("g%03d", seq_along(dnmb$genome_meta$genome_uid)),
    as.character(dnmb$genome_meta$genome_uid)
  )
  safe_contig <- function(x) gsub("[^A-Za-z0-9_.-]", "_", as.character(x))
  gff_path <- file.path(root, "genomes.gff")
  gff <- dnmb$id_map %>%
    dplyr::transmute(
      chrom = paste0(unname(gid_tag[as.character(genome_uid)]),
                     "_", safe_contig(contig)),
      id    = paste0("p", protein_uid, "_g", genome_uid),
      start = as.integer(start),
      end   = as.integer(end)
    ) %>%
    dplyr::arrange(chrom, start)
  utils::write.table(gff[, c("chrom", "id", "start", "end")],
                     gff_path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)

  # .blast: pairs within each OG. Derive membership from the exported
  # alignment FASTA (headers of the form "p<uid>_g<genome_uid>") so the
  # export is robust to refined cluster IDs (refine_graph = TRUE rewrites
  # dnmb$clusters inside run_orthofinder_like(), but the caller's dnmb
  # may still hold the original space). Streams pairs line-by-line so
  # large families do not materialize combn() up front.
  blast_path <- file.path(root, "genomes.blast")
  con <- file(blast_path, "w"); on.exit(close(con), add = TRUE)
  fmt <- "%s\t%s\t100.0\t200\t0\t0\t1\t200\t1\t200\t1e-50\t100\n"
  for (i in seq_len(nrow(ok))) {
    hdrs <- grep("^>", readLines(ok$aln_path[i], warn = FALSE), value = TRUE)
    ids <- unique(sub("^>", "", hdrs))
    n <- length(ids)
    if (n < 2L) next
    for (a in seq_len(n - 1L)) {
      ida <- ids[a]
      for (b in (a + 1L):n) {
        cat(sprintf(fmt, ida, ids[b]), file = con)
      }
    }
  }
  root
}
