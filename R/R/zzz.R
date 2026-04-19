utils::globalVariables(c(
  ".",
  "aa_length", "absent", "accessory", "alpha_val", "angle",
  "bitscore",
  "category", "cds_key", "cluster_id", "cog_category", "contig",
  "contig_name", "core", "core_max", "core_mean", "core_min",
  "cumulative",
  "description", "Dim1", "Dim2",
  "e_shift", "end", "enriched",
  "fill", "fill_class", "functional_class",
  "gc_percent", "gene_label", "genome_a", "genome_b", "genome_key",
  "genome_uid", "gl_label", "group",
  "hjust",
  "is_anchor", "is_centroid", "is_tip", "isTip",
  "jaccard",
  "k",
  "label", "locus_tag", "log2OR",
  "max_paralogs", "median", "median_identity", "member_coverage",
  "mid_x",
  "n", "n_cds", "n_clusters", "n_contrib", "n_gained", "n_genes",
  "n_genomes", "n_lost", "n_members", "n_unique",
  "neg_log_p", "new_cluster", "new_end", "new_start", "new_strand",
  "ng", "node",
  "order_label", "organism", "orig_left_bp", "orig_right_bp",
  "pan", "pan_category", "pan_max", "pan_mean", "pan_min",
  "paralogs", "pct", "pct_identity_fwd", "protein_uid",
  "qseqid",
  "rep_coverage", "representative_product", "representative_uid",
  "s_shift", "sig", "species_a", "species_b", "start", "strain",
  "strand",
  "text_angle", "text_x", "text_y", "top", "top_products",
  "total", "total_length", "translation",
  "x", "y", "y_num"
))

# Resolve pairwiseAlignment across Biostrings / pwalign generations.
# Biostrings >= 2.77.1 (Bioconductor 3.19+) moved pairwiseAlignment into
# its own pwalign package and defuncted the old entry point. Older
# Biostrings still carries it. Returning a closure means the dispatch
# cost is paid once per caller.
.dnmb_pairwise_fn <- function() {
  if (requireNamespace("pwalign", quietly = TRUE))
    return(pwalign::pairwiseAlignment)
  if (requireNamespace("Biostrings", quietly = TRUE))
    return(Biostrings::pairwiseAlignment)
  stop("Neither pwalign nor Biostrings installed; cannot run pairwise alignment.")
}
