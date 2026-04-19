#' Build a concatenation-based (supermatrix) species tree
#'
#' Phylogenomic supermatrix approach: for every single-copy orthogroup
#' (SC-OG), concatenate its amino-acid alignment column-wise across
#' species. The resulting supermatrix is one long alignment with one
#' row per species. Build the species tree from that matrix via
#' `phangorn::pml + optim.pml` (ML) or NJ (distance).
#'
#' This is the classic alternative to STAG's distance-summarization
#' approach. It is more robust when species coverage is high and
#' conflict among gene trees is low. It is unusable when there are no
#' SC-OGs or when paralog decay has replaced single-copy signal with
#' many-copy noise — in those cases prefer `stag_species_tree()`.
#'
#' Species that lack a sequence in an SC-OG receive gaps (`-`) for that
#' block. Species with zero total coverage are dropped with a warning.
#'
#' @param og_result Tibble from `per_og_trees()`.
#' @param dnmb `load_dnmb()` result.
#' @param method `"nj"` (fast, distance) or `"ml"` (phangorn pml + NNI).
#' @param model Substitution model passed to `phangorn::pml(..., model=...)`
#'   (ML only; default `"LG"`).
#' @param max_ogs Cap the number of SC-OGs used (picked by member count,
#'   largest first). NULL = use all.
#' @param bootstrap Non-parametric bootstrap replicates. 0 disables.
#' @param bootstrap_seed RNG seed for reproducibility.
#' @param partitioned If TRUE and `method = "ml"`, perform a partitioned
#'   ML fit via `phangorn::pmlPart`, letting each SC-OG partition carry
#'   its own Γ rate / edge-scaler while sharing topology. Falls back to
#'   concatenated ML if `pmlPart` is unavailable or errors.
#' @param gamma_k Γ rate categories (ML only). 0 disables. Default 4.
#' @param out_path If non-NULL, write the tree to this Newick file.
#' @return An `ape::phylo` tree with tip labels = `genome_key`. When
#'   `bootstrap > 0`, `$node.label` carries 0–100 bipartition support.
#' @export
supermatrix_species_tree <- function(og_result,
                                     dnmb,
                                     method = c("nj", "ml"),
                                     model = "LG",
                                     max_ogs = NULL,
                                     bootstrap = 0L,
                                     bootstrap_seed = 42L,
                                     partitioned = FALSE,
                                     gamma_k = 4L,
                                     out_path = NULL) {
  method <- match.arg(method)
  if (!requireNamespace("ape", quietly = TRUE)) stop("ape required.")

  ok <- og_result[!is.na(og_result$aln_path) & og_result$single_copy, , drop = FALSE]
  if (!nrow(ok)) stop("supermatrix_species_tree: no SC-OGs available.")

  if (!is.null(max_ogs) && nrow(ok) > max_ogs) {
    ok <- ok[order(-ok$n_members), , drop = FALSE][seq_len(max_ogs), , drop = FALSE]
  }

  uid2key <- stats::setNames(dnmb$genome_meta$genome_key,
                      as.character(dnmb$genome_meta$genome_uid))
  species <- unname(uid2key)

  blocks <- .concatenate_sc_ogs(ok, species, uid2key)
  super  <- blocks$matrix
  block_widths <- blocks$widths
  keep_sp <- rowSums(super != "-") > 0L
  if (!all(keep_sp)) {
    warning("supermatrix_species_tree: dropping species with zero SC-OG coverage: ",
            paste(species[!keep_sp], collapse = ", "))
    super <- super[keep_sp, , drop = FALSE]
  }
  if (nrow(super) < 3L) stop("supermatrix: < 3 species have any coverage.")

  main_tree <- if (partitioned && method == "ml" &&
                   length(block_widths) >= 2L &&
                   requireNamespace("phangorn", quietly = TRUE)) {
    tryCatch(
      .supermatrix_partitioned_ml(super, block_widths, model, gamma_k),
      error = function(e) {
        warning("supermatrix: pmlPart failed (", conditionMessage(e),
                "); falling back to concatenated ML.")
        .supermatrix_tree_from_char_matrix(super, method, model, gamma_k)
      }
    )
  } else {
    .supermatrix_tree_from_char_matrix(super, method, model, gamma_k)
  }

  if (bootstrap > 0L) {
    set.seed(bootstrap_seed)
    n_col <- ncol(super)
    reps <- vector("list", bootstrap)
    for (b in seq_len(bootstrap)) {
      idx <- sample.int(n_col, n_col, replace = TRUE)
      reps[[b]] <- .supermatrix_tree_from_char_matrix(
        super[, idx, drop = FALSE], method, model, gamma_k)
    }
    support <- ape::prop.clades(main_tree, reps, rooted = FALSE)
    main_tree$node.label <- round(100 * support / bootstrap)
  }

  if (!is.null(out_path)) ape::write.tree(main_tree, out_path)
  main_tree
}


# ---- helpers ---------------------------------------------------------------

.concatenate_sc_ogs <- function(ok, species, uid2key) {
  mats <- list()
  widths <- integer()
  for (i in seq_len(nrow(ok))) {
    aln_path <- ok$aln_path[i]
    if (!file.exists(aln_path)) next
    aa <- if (requireNamespace("Biostrings", quietly = TRUE)) {
      bs <- Biostrings::readAAStringSet(aln_path)
      stats::setNames(as.character(bs), names(bs))
    } else {
      .read_fasta_chr(aln_path)
    }
    gidx <- sub(".*_g", "", names(aa))
    keys <- unname(uid2key[gidx])
    if (anyDuplicated(keys)) next  # not truly single-copy after all; skip
    keep <- !is.na(keys)
    aa <- aa[keep]; keys <- keys[keep]
    if (!length(aa)) next

    width <- unique(nchar(aa))
    if (length(width) != 1L) next
    block <- matrix("-", nrow = length(species), ncol = width,
                    dimnames = list(species, NULL))
    for (j in seq_along(aa)) {
      chars <- strsplit(aa[j], "", fixed = TRUE)[[1]]
      block[keys[j], ] <- chars
    }
    mats[[length(mats) + 1]] <- block
    widths <- c(widths, width)
  }
  if (!length(mats)) stop(".concatenate_sc_ogs: no usable SC-OG alignments.")
  list(matrix = do.call(cbind, mats), widths = widths)
}

.read_fasta_chr <- function(path) {
  lines <- readLines(path)
  headers <- grep("^>", lines)
  out <- character(length(headers))
  nm  <- character(length(headers))
  ends <- c(utils::tail(headers - 1, -1), length(lines))
  for (i in seq_along(headers)) {
    nm[i] <- sub("^>\\s*(\\S+).*", "\\1", lines[headers[i]])
    out[i] <- paste(lines[(headers[i] + 1):ends[i]], collapse = "")
  }
  stats::setNames(out, nm)
}

.supermatrix_tree_from_char_matrix <- function(mat, method, model, gamma_k = 0L) {
  concat <- apply(mat, 1L, paste, collapse = "")
  tmp <- tempfile(fileext = ".fa")
  on.exit(unlink(tmp), add = TRUE)
  con <- file(tmp, "w")
  for (sp in names(concat)) {
    cat(">", sp, "\n", concat[sp], "\n", sep = "", file = con)
  }
  close(con)

  use_gamma <- isTRUE(gamma_k > 0L)
  if (method == "ml" && requireNamespace("phangorn", quietly = TRUE)) {
    aln   <- phangorn::read.phyDat(tmp, format = "fasta", type = "AA")
    dm    <- phangorn::dist.ml(aln, model = model)
    start <- ape::nj(dm)
    fit   <- if (use_gamma) {
      phangorn::pml(start, data = aln, model = model, k = as.integer(gamma_k))
    } else {
      phangorn::pml(start, data = aln, model = model)
    }
    fit   <- phangorn::optim.pml(
      fit, model = model, optNni = TRUE, optGamma = use_gamma,
      control = phangorn::pml.control(trace = 0))
    return(fit$tree)
  }
  if (requireNamespace("phangorn", quietly = TRUE)) {
    aln <- phangorn::read.phyDat(tmp, format = "fasta", type = "AA")
    return(ape::nj(phangorn::dist.ml(aln, model = model)))
  }
  aln <- ape::read.FASTA(tmp, type = "AA")
  ape::nj(ape::dist.aa(aln))
}

# Per-partition ML via phangorn::pmlPart. Builds one phyDat per block
# directly in memory (phangorn::phyDat accepts a character matrix of
# AA letters — avoids 1-tempfile-per-partition for runs with 1000+ OGs).
# Joint-optimizes tree + edge lengths across partitions with
# per-partition rate / Γ shape.
.supermatrix_partitioned_ml <- function(mat, block_widths, model, gamma_k) {
  stopifnot(sum(block_widths) == ncol(mat))
  use_gamma <- isTRUE(gamma_k > 0L)

  aln_full <- phangorn::phyDat(mat, type = "AA")
  start <- ape::nj(phangorn::dist.ml(aln_full, model = model))

  boundaries <- c(0L, cumsum(block_widths))
  fits <- vector("list", length(block_widths))
  for (k in seq_along(block_widths)) {
    cols <- (boundaries[k] + 1L):boundaries[k + 1L]
    aln_k <- phangorn::phyDat(mat[, cols, drop = FALSE], type = "AA")
    fits[[k]] <- if (use_gamma) {
      phangorn::pml(start, data = aln_k, model = model, k = as.integer(gamma_k))
    } else {
      phangorn::pml(start, data = aln_k, model = model)
    }
  }

  joint <- phangorn::pmlPart(
    if (use_gamma) edge + nni + shape ~ . else edge + nni ~ .,
    fits,
    control = phangorn::pml.control(trace = 0)
  )
  joint$fits[[1]]$tree
}
