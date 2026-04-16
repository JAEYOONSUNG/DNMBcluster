#' Load all DNMBcluster Parquet outputs into a named list
#'
#' Reads the seven core dataframes produced by the Python pipeline:
#' `id_map`, `gene_table`, `genome_meta`, `clusters`, `presence_absence`,
#' `pan_core_curve`, and `cluster_summary`. All are loaded via
#' `arrow::read_parquet` (fast, zero-copy where possible) and returned as
#' a list suitable for downstream plot functions.
#'
#' Layout (post-reorg): `results_dir/dnmb/{inputs, raw, processed}`.
#' Engine-agnostic input artifacts (id_map, gene_table, genome_meta) live
#' in `inputs/`; canonical clustering parquet (clusters,
#' presence_absence, pan_core_curve, cluster_summary, cluster_long) lives
#' in `processed/`.
#'
#' @param results_dir Path to a DNMBcluster results directory.
#' @return Named list with seven tibbles.
#' @export
load_dnmb <- function(results_dir) {
  dnmb_dir <- file.path(results_dir, "dnmb")
  inputs_dir <- file.path(dnmb_dir, "inputs")
  processed_dir <- file.path(dnmb_dir, "processed")

  if (!dir.exists(dnmb_dir)) {
    stop("DNMB results directory not found: ", dnmb_dir)
  }

  read_at <- function(dir, name) {
    path <- file.path(dir, paste0(name, ".parquet"))
    if (!file.exists(path)) {
      stop("Missing DNMB Parquet file: ", path)
    }
    tibble::as_tibble(arrow::read_parquet(path))
  }

  list(
    id_map           = read_at(inputs_dir, "id_map"),
    gene_table       = read_at(inputs_dir, "gene_table"),
    genome_meta      = read_at(inputs_dir, "genome_meta"),
    clusters         = read_at(processed_dir, "clusters"),
    presence_absence = read_at(processed_dir, "presence_absence"),
    pan_core_curve   = read_at(processed_dir, "pan_core_curve"),
    cluster_summary  = read_at(processed_dir, "cluster_summary")
  )
}


#' Compute per-genome core / unique / accessory / absent cluster counts
#'
#' Replicates the four BPGA flower-plot categories without needing a
#' pre-computed `stats.xls` — the numbers come directly from the
#' `clusters` dataframe.
#'
#' @param dnmb Output of `load_dnmb()`.
#' @return Tibble with one row per genome and columns
#'   `genome_uid, genome_key, organism, strain, core, unique, accessory, absent`.
#' @export
compute_per_genome_stats <- function(dnmb) {
  n_total <- nrow(dnmb$genome_meta)

  # cluster_id x genome_uid presence bool
  presence_pairs <- dnmb$clusters %>%
    dplyr::select(cluster_id, genome_uid) %>%
    dplyr::distinct()

  n_per_cluster <- presence_pairs %>%
    dplyr::count(cluster_id, name = "n_genomes")

  # For each genome, walk its cluster memberships and tally the four
  # categories. Vectorized via a single join + group_by; no row loops.
  scored <- presence_pairs %>%
    dplyr::left_join(n_per_cluster, by = "cluster_id") %>%
    dplyr::group_by(genome_uid) %>%
    dplyr::summarise(
      core      = sum(n_genomes == n_total),
      unique    = sum(n_genomes == 1L),
      accessory = sum(n_genomes > 1L & n_genomes < n_total),
      .groups   = "drop"
    )

  # Exclusively absent: clusters where n_genomes == (n_total - 1) and
  # this genome is NOT one of the members. Compute via anti-join.
  exclusively_absent <- lapply(dnmb$genome_meta$genome_uid, function(gid) {
    candidate_clusters <- n_per_cluster$cluster_id[n_per_cluster$n_genomes == n_total - 1]
    present_clusters <- presence_pairs$cluster_id[presence_pairs$genome_uid == gid]
    data.frame(
      genome_uid = gid,
      absent     = length(setdiff(candidate_clusters, present_clusters))
    )
  })
  exclusively_absent <- do.call(rbind, exclusively_absent)

  dnmb$genome_meta %>%
    dplyr::select(genome_uid, genome_key, organism, strain) %>%
    dplyr::left_join(scored, by = "genome_uid") %>%
    dplyr::left_join(exclusively_absent, by = "genome_uid") %>%
    dplyr::mutate(
      core      = dplyr::coalesce(core, 0L),
      unique    = dplyr::coalesce(unique, 0L),
      accessory = dplyr::coalesce(accessory, 0L),
      absent    = dplyr::coalesce(absent, 0L)
    )
}
