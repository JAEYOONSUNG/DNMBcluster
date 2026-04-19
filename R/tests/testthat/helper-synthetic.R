has_aln_engine <- function() {
  (requireNamespace("DECIPHER", quietly = TRUE) &&
     requireNamespace("Biostrings", quietly = TRUE)) ||
    nzchar(Sys.which("mafft"))
}

skip_if_no_aln_engine <- function() {
  testthat::skip_if_not(has_aln_engine(),
                        "No alignment engine (DECIPHER/Biostrings or MAFFT) available.")
}

make_synthetic_dnmb <- function(n_genomes = 4L, n_ogs = 5L, seed = 1L) {
  set.seed(seed)
  aa <- c("A","C","D","E","F","G","H","I","K","L",
          "M","N","P","Q","R","S","T","V","W","Y")

  genome_meta <- tibble::tibble(
    genome_uid = seq_len(n_genomes),
    genome_key = paste0("g", seq_len(n_genomes))
  )

  gene_rows    <- list()
  cluster_rows <- list()
  id_rows      <- list()
  puid <- 0L

  for (og in seq_len(n_ogs)) {
    root_seq <- paste(sample(aa, 120, replace = TRUE), collapse = "")
    for (g in seq_len(n_genomes)) {
      puid <- puid + 1L
      # tiny per-genome mutation
      vec <- strsplit(root_seq, "")[[1]]
      mutate_at <- sample.int(length(vec), size = 3 + og)
      vec[mutate_at] <- sample(aa, length(mutate_at), replace = TRUE)
      seq <- paste(vec, collapse = "")
      gene_rows[[length(gene_rows) + 1]] <- tibble::tibble(
        protein_uid = puid,
        translation = seq,
        length      = nchar(seq)
      )
      cluster_rows[[length(cluster_rows) + 1]] <- tibble::tibble(
        cluster_id  = og,
        protein_uid = puid,
        genome_uid  = g,
        is_centroid = g == 1L
      )
      id_rows[[length(id_rows) + 1]] <- tibble::tibble(
        protein_uid = puid,
        genome_uid  = g,
        contig      = paste0("c", g),
        start       = puid * 1000L,
        end         = puid * 1000L + 300L
      )
    }
  }

  list(
    genome_meta = genome_meta,
    gene_table  = dplyr::bind_rows(gene_rows),
    clusters    = dplyr::bind_rows(cluster_rows),
    id_map      = dplyr::bind_rows(id_rows)
  )
}

write_synthetic_dnmb_parquets <- function(dnmb, base_dir) {
  if (!requireNamespace("arrow", quietly = TRUE))
    testthat::skip("arrow package not available")
  inputs_dir <- file.path(base_dir, "dnmb", "inputs")
  processed_dir <- file.path(base_dir, "dnmb", "processed")
  dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

  gm <- dnmb$genome_meta
  if (!"organism" %in% names(gm)) gm$organism <- gm$genome_key
  if (!"strain"   %in% names(gm)) gm$strain   <- gm$genome_key

  arrow::write_parquet(dnmb$id_map,      file.path(inputs_dir, "id_map.parquet"))
  arrow::write_parquet(dnmb$gene_table,  file.path(inputs_dir, "gene_table.parquet"))
  arrow::write_parquet(gm,               file.path(inputs_dir, "genome_meta.parquet"))
  arrow::write_parquet(dnmb$clusters,    file.path(processed_dir, "clusters.parquet"))

  pa <- dnmb$clusters %>%
    dplyr::distinct(cluster_id, genome_uid) %>%
    dplyr::mutate(present = 1L)
  arrow::write_parquet(pa, file.path(processed_dir, "presence_absence.parquet"))

  pcc <- tibble::tibble(n_genomes = seq_len(nrow(gm)),
                        pan = seq_len(nrow(gm)),
                        core = rev(seq_len(nrow(gm))))
  arrow::write_parquet(pcc, file.path(processed_dir, "pan_core_curve.parquet"))

  cs <- dnmb$clusters %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(n_genomes = dplyr::n_distinct(genome_uid),
                     n_members = dplyr::n(), .groups = "drop")
  arrow::write_parquet(cs, file.path(processed_dir, "cluster_summary.parquet"))
  invisible(base_dir)
}
