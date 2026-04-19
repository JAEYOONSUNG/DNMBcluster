suppressPackageStartupMessages({
  # When DNMB_PKG_ROOT is set, always load from source (development
  # use). Otherwise fall back to the installed package.
  pkg_root <- Sys.getenv("DNMB_PKG_ROOT", unset = "")
  if (nzchar(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(pkg_root, quiet = TRUE)
  } else if (requireNamespace("DNMBcluster", quietly = TRUE)) {
    library(DNMBcluster)
  } else {
    stop("Install DNMBcluster, or set DNMB_PKG_ROOT and install devtools.")
  }
  library(dplyr)
  library(tibble)
  library(tidyr)
})

# =============================================================== #
# Full DNMBcluster figure gallery demo
# =============================================================== #
# Builds a realistic synthetic pangenome (8 genomes x 200 OGs) with
# core / accessory / unique categories, ANI/POCP matrices, a fake
# core-gene newick tree, COG/eggNOG/functional annotations, then
# calls run_dnmb_plot() + the 5 selection plots in one shot.
#
# Output directory is controlled by the DNMB_DEMO_OUT env var;
# defaults to a fresh tempfile("dnmb_demo_") path.
# Example:
#   DNMB_DEMO_OUT=/tmp/my_demo Rscript full_demo.R

set.seed(7)

# 8 genomes is the sweet spot -- readable Euler diagram + fast
# eulerr::euler() convergence (12 sets takes >10 min to fit).
N_GENOMES <- 8L
N_OGS     <- 200L

out_root <- Sys.getenv("DNMB_DEMO_OUT", unset = tempfile("dnmb_demo_"))

# Safety guard: the script auto-wipes `out_root` before every run. A
# typo like DNMB_DEMO_OUT=$HOME or DNMB_DEMO_OUT=/ would then nuke
# the user's home. Refuse anything that isn't clearly a dnmb_demo
# scratch directory (tempdir-rooted, or explicitly prefixed).
.dnmb_demo_safe_path <- function(p) {
  if (!nzchar(p)) return(FALSE)
  norm <- normalizePath(p, winslash = "/", mustWork = FALSE)
  if (!nzchar(norm) || norm %in% c("", "/")) return(FALSE)
  home <- normalizePath(Sys.getenv("HOME", unset = "~"),
                        winslash = "/", mustWork = FALSE)
  cwd  <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if (norm %in% c(home, cwd)) return(FALSE)
  tmp  <- normalizePath(tempdir(), winslash = "/", mustWork = FALSE)
  under_tmp <- nzchar(tmp) && startsWith(paste0(norm, "/"), paste0(tmp, "/"))
  has_prefix <- grepl("dnmb[_-]?demo", basename(norm), ignore.case = TRUE)
  under_tmp || has_prefix
}
if (dir.exists(out_root)) {
  if (!.dnmb_demo_safe_path(out_root)) {
    stop("Refusing to wipe DNMB_DEMO_OUT='", out_root,
         "'. Use a path under tempdir() or one whose basename ",
         "contains 'dnmb_demo'.")
  }
  unlink(out_root, recursive = TRUE)
}
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

inputs_dir    <- file.path(out_root, "dnmb", "inputs")
processed_dir <- file.path(out_root, "dnmb", "processed")
annot_dir     <- file.path(out_root, "dnmb", "annotations")
dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Genome meta --------------------------------------------
organisms <- c("Bacillus subtilis", "Bacillus cereus", "Bacillus anthracis",
                "Bacillus thuringiensis", "Bacillus licheniformis",
                "Bacillus pumilus", "Bacillus atrophaeus", "Bacillus megaterium",
                "Bacillus mycoides", "Bacillus weihenstephanensis",
                "Bacillus amyloliquefaciens", "Bacillus velezensis")
strains <- sprintf("strain-%02d", seq_len(N_GENOMES))

genome_meta <- tibble::tibble(
  genome_uid   = seq_len(N_GENOMES),
  genome_key   = sprintf("g%02d", seq_len(N_GENOMES)),
  organism     = organisms[seq_len(N_GENOMES)],
  strain       = strains,
  gc_percent   = round(rnorm(N_GENOMES, 42, 3), 2),
  total_length = round(rnorm(N_GENOMES, 4.2e6, 3e5)),
  n_cds        = round(rnorm(N_GENOMES, 4100, 150))
)

# ---- 2. Clusters with proper core/accessory/unique ----------------
n_core      <- 120L
n_accessory <- 60L
n_unique    <- N_OGS - n_core - n_accessory

# Pool of product annotations to sample from for each cluster rep.
products_pool <- c(
  "DNA-directed RNA polymerase subunit beta",
  "50S ribosomal protein L2", "30S ribosomal protein S3",
  "ATP synthase subunit alpha", "elongation factor Tu",
  "chaperonin GroEL", "DNA gyrase subunit A",
  "cold-shock protein CspA", "two-component sensor histidine kinase",
  "glycosyltransferase family 2", "ABC transporter permease",
  "transposase IS256 family", "phage tail fiber protein",
  "acyl-CoA dehydrogenase family member", "LysR family transcriptional regulator",
  "sporulation protein SpoIIID", "polyketide synthase module",
  "nonribosomal peptide synthetase", "cytochrome c oxidase subunit II",
  "methyltransferase-like protein", "lipoprotein signal peptidase",
  "hypothetical protein", "restriction endonuclease family",
  "MFS transporter superfamily", "flagellar hook-associated protein"
)

cluster_rows <- list()
id_rows      <- list()
gene_rows    <- list()
rep_uid_per_og      <- integer(N_OGS)
rep_product_per_og  <- character(N_OGS)
puid <- 0L
aa <- c("A","C","D","E","F","G","H","I","K","L",
        "M","N","P","Q","R","S","T","V","W","Y")

for (og in seq_len(N_OGS)) {
  if (og <= n_core) {
    members <- seq_len(N_GENOMES)
  } else if (og <= n_core + n_accessory) {
    k <- sample(2:(N_GENOMES - 1L), 1)
    members <- sort(sample(seq_len(N_GENOMES), k))
  } else {
    members <- sample(seq_len(N_GENOMES), 1L)
  }
  # Longer reps for core OGs, short reps for unique (matches biology).
  rep_len <- switch(
    if (og <= n_core) "core" else if (og <= n_core + n_accessory) "acc" else "uniq",
    core = round(rnorm(1, 300, 80)),
    acc  = round(rnorm(1, 220, 90)),
    uniq = round(rnorm(1, 120, 60))
  )
  rep_len <- max(60L, rep_len)

  root <- paste(sample(aa, rep_len, replace = TRUE), collapse = "")
  product <- sample(products_pool, 1L)
  rep_product_per_og[og] <- product

  for (g in members) {
    puid <- puid + 1L
    vec <- strsplit(root, "")[[1]]
    mut <- sample.int(length(vec), size = min(length(vec), sample(3:10, 1)))
    vec[mut] <- sample(aa, length(mut), replace = TRUE)
    seq <- paste(vec, collapse = "")
    aa_len <- nchar(seq)
    is_cent <- g == members[1L]
    if (is_cent) rep_uid_per_og[og] <- puid

    gene_rows[[length(gene_rows) + 1L]] <- tibble::tibble(
      protein_uid = puid, translation = seq, length = aa_len,
      aa_length   = aa_len)
    cluster_rows[[length(cluster_rows) + 1L]] <- tibble::tibble(
      cluster_id = og, protein_uid = puid, genome_uid = g,
      is_centroid = is_cent,
      pct_identity_fwd = if (is_cent) 100.0 else
        pmin(100, pmax(45, rnorm(1, 88, 6))),
      member_coverage  = if (is_cent) 1.0 else
        pmin(1, pmax(0.4, rnorm(1, 0.93, 0.05))),
      rep_coverage     = if (is_cent) 1.0 else
        pmin(1, pmax(0.4, rnorm(1, 0.96, 0.04))))
    id_rows[[length(id_rows) + 1L]] <- tibble::tibble(
      protein_uid = puid, genome_uid = g,
      contig = paste0("c", g),
      start = puid * 1200L, end = puid * 1200L + 600L,
      aa_length = aa_len)
  }
}
clusters <- dplyr::bind_rows(cluster_rows)
id_map   <- dplyr::bind_rows(id_rows)
genes    <- dplyr::bind_rows(gene_rows)

# ---- 3. DNMB processed Parquets --------------------------------
cs <- clusters %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::summarise(n_genomes = dplyr::n_distinct(genome_uid),
                   n_members = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(category = dplyr::case_when(
    n_genomes == N_GENOMES ~ "core",
    n_genomes == 1L        ~ "unique",
    TRUE                   ~ "accessory"),
    representative_uid     = rep_uid_per_og[cluster_id],
    representative_product = rep_product_per_og[cluster_id])

pa <- clusters %>%
  dplyr::distinct(cluster_id, genome_uid) %>%
  dplyr::mutate(present = 1L)

# Pan/core accumulation curve -- Heaps' law shape, 10 bootstrap reps per k.
# load_dnmb / pan_core_plot expects the column `k`.
n_boot <- 10L
boot_rows <- list()
for (k in seq_len(N_GENOMES)) {
  for (b in seq_len(n_boot)) {
    pan  <- round(60 + 130 * k^0.6  + rnorm(1, 0, 4))
    core <- round(200 * k^(-0.25)   + rnorm(1, 0, 2))
    boot_rows[[length(boot_rows) + 1L]] <- tibble::tibble(
      k = k, bootstrap = b, pan = pan, core = core)
  }
}
pcc <- dplyr::bind_rows(boot_rows)

arrow::write_parquet(id_map,     file.path(inputs_dir, "id_map.parquet"))
arrow::write_parquet(genes,      file.path(inputs_dir, "gene_table.parquet"))
arrow::write_parquet(genome_meta, file.path(inputs_dir, "genome_meta.parquet"))
arrow::write_parquet(clusters,   file.path(processed_dir, "clusters.parquet"))
arrow::write_parquet(pa,         file.path(processed_dir, "presence_absence.parquet"))
arrow::write_parquet(pcc,        file.path(processed_dir, "pan_core_curve.parquet"))
arrow::write_parquet(cs,         file.path(processed_dir, "cluster_summary.parquet"))

# ---- 4. ANI + POCP matrices (pairwise, percent scale) ----------
# Plot code expects columns `ani_percent` and `pocp_percent`.
ani_raw <- matrix(85 + 13 * runif(N_GENOMES^2), N_GENOMES, N_GENOMES)
ani_raw <- (ani_raw + t(ani_raw)) / 2
diag(ani_raw) <- 100
rownames(ani_raw) <- genome_meta$genome_key
colnames(ani_raw) <- genome_meta$genome_key

ani_long <- tibble::as_tibble(as.data.frame(as.table(ani_raw))) %>%
  dplyr::rename(genome_a = Var1, genome_b = Var2, ani_percent = Freq) %>%
  dplyr::mutate(genome_a = as.character(genome_a),
                genome_b = as.character(genome_b))
arrow::write_parquet(ani_long, file.path(processed_dir, "ani_matrix.parquet"))

pocp_raw <- ani_raw - 5 + rnorm(N_GENOMES^2, 0, 2)
pocp_raw <- (pocp_raw + t(pocp_raw)) / 2
pocp_raw[] <- pmin(100, pmax(40, pocp_raw))
diag(pocp_raw) <- 100
pocp_long <- tibble::as_tibble(as.data.frame(as.table(pocp_raw))) %>%
  dplyr::rename(genome_a = Var1, genome_b = Var2, pocp_percent = Freq) %>%
  dplyr::mutate(genome_a = as.character(genome_a),
                genome_b = as.character(genome_b))
arrow::write_parquet(pocp_long, file.path(processed_dir, "pocp_matrix.parquet"))

# ---- 5. Fake eggNOG / COG annotation tables --------------------
cog_letters <- c("J","A","K","L","B","D","Y","V","T","M","N","Z","W","U","O",
                 "C","G","E","F","H","I","P","Q","R","S")
# cog_annotations.parquet: per-protein COG letters (used nowhere directly,
# but keep for completeness).
cog_annot <- tibble::tibble(
  protein_uid = genes$protein_uid,
  cog_category = sample(cog_letters, nrow(genes), replace = TRUE,
                         prob = c(rep(2, 8), rep(4, 8), rep(6, 8), 8))
)
arrow::write_parquet(cog_annot, file.path(annot_dir, "cog_annotations.parquet"))

# eggnog_annotations.parquet: per-cluster COG letter + description.
# Required columns: cluster_id, cog_category. Let core clusters skew
# toward housekeeping letters (J,L,C,E), accessory toward regulators/
# transporters (K,T,P,G), and unique toward mobile/unknown (X,S).
draw_cog <- function(category) {
  if (category == "core") {
    sample(c("J","L","C","E","F","H","G","M","O"), 1,
           prob = c(6,5,5,5,4,3,3,3,3))
  } else if (category == "accessory") {
    sample(c("K","T","P","G","V","N","U","R"), 1,
           prob = c(5,5,4,4,3,3,3,2))
  } else {
    sample(c("S","R","X","L"), 1, prob = c(6,3,3,2))
  }
}
eggnog <- cs %>%
  dplyr::rowwise() %>%
  dplyr::mutate(cog_category = draw_cog(category)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cluster_id, cog_category) %>%
  dplyr::mutate(description = "auto-synthetic")
arrow::write_parquet(eggnog, file.path(processed_dir, "eggnog_annotations.parquet"))

# functional_categories.parquet: cluster_id, pan_category, functional_class
functional_classes <- c(
  "Translation", "DNA replication", "Energy metabolism", "Amino acid biosynthesis",
  "Transcription regulation", "Membrane transport", "Carbohydrate metabolism",
  "Signal transduction", "Mobile elements", "Hypothetical",
  "Secondary metabolism", "Stress response", "Cell wall / envelope",
  "Sporulation", "Chaperones / PTM"
)
draw_func <- function(category) {
  if (category == "core") {
    sample(functional_classes[1:4], 1, prob = c(4, 4, 3, 3))
  } else if (category == "accessory") {
    sample(functional_classes[5:8], 1, prob = c(4, 3, 3, 2))
  } else {
    sample(functional_classes[9:15], 1, prob = c(3, 5, 2, 2, 2, 2, 2))
  }
}
func <- cs %>%
  dplyr::rowwise() %>%
  dplyr::mutate(functional_class = draw_func(category)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cluster_id, functional_class) %>%
  dplyr::mutate(pan_category = cs$category[match(cluster_id, cs$cluster_id)])
arrow::write_parquet(func, file.path(processed_dir, "functional_categories.parquet"))

# ---- 6. Fake core-gene newick tree -----------------------------
tree <- ape::rtree(N_GENOMES, tip.label = genome_meta$genome_key)
tree$edge.length <- pmax(tree$edge.length, 0.01)
# treeio::read.newick(node.label = "support") requires numeric node labels.
tree$node.label <- as.character(sample(60:100, tree$Nnode, replace = TRUE))
ape::write.tree(tree, file.path(processed_dir, "phylo_tree.nwk"))

cat("[demo] synthetic DNMB assets written\n")
cat("  genomes     = ", N_GENOMES, "\n")
cat("  OGs         = ", N_OGS, " (core=", n_core, "/acc=", n_accessory,
     "/unique=", n_unique, ")\n", sep = "")

# ---- 7. Load + run the full plot suite -------------------------
dnmb <- load_dnmb(out_root)

# Build a synthetic selection result to exercise all 5 HyPhy figures.
set.seed(42)
K <- 40L
busted <- tibble::tibble(
  cluster_id = sort(sample(cs$cluster_id, K)),
  method = "BUSTED",
  p_value = c(runif(5, 1e-4, 5e-3), runif(5, 5e-3, 5e-2),
               runif(K - 10L, 0.05, 0.99)),
  lrt = NA_real_,
  omega_positive = c(rlnorm(5, log(4), 0.4), rlnorm(5, log(2), 0.4),
                      rlnorm(K - 10L, log(0.6), 0.4)),
  weight_positive = c(runif(5, 0.05, 0.25), runif(5, 0.05, 0.3),
                       runif(K - 10L, 0.01, 0.15)),
  json_path = NA_character_, converged = TRUE, reason = NA_character_
)
busted$q_value <- stats::p.adjust(busted$p_value, method = "BH")

fel <- tibble::tibble(
  cluster_id = busted$cluster_id, method = "FEL",
  n_sites = sample(50:400, K, replace = TRUE),
  n_positive = c(sample(3:12, 5, replace = TRUE),
                  sample(0:3, K - 5L, replace = TRUE)),
  n_negative = sample(5:60, K, replace = TRUE),
  median_omega = c(rlnorm(5, log(2.5), 0.3),
                    rlnorm(K - 5L, log(0.5), 0.4)),
  json_path = NA_character_, converged = TRUE, reason = NA_character_
)

absrel <- tibble::tibble(
  cluster_id = busted$cluster_id, method = "aBSREL",
  n_branches = sample(4:12, K, replace = TRUE),
  n_selected = c(sample(1:3, 5, replace = TRUE),
                 sample(0:1, K - 5L, replace = TRUE)),
  min_corrected_p = c(runif(5, 1e-4, 5e-3), runif(5, 5e-3, 5e-2),
                       runif(K - 10L, 0.05, 0.99)),
  json_path = NA_character_, converged = TRUE, reason = NA_character_
)
absrel$q_value <- stats::p.adjust(absrel$min_corrected_p, method = "BH")

codeml <- tibble::tibble(
  cluster_id = busted$cluster_id,
  omega = pmin(busted$omega_positive * runif(K, 0.15, 0.35) + 0.05, 3),
  dN = NA_real_, dS = NA_real_, log_lik = NA_real_,
  converged = TRUE, reason = NA_character_
)

# Synthetic FEL JSON for the first OG so run_dnmb_plot() can render
# the per-site FEL scatter automatically.
fel_dir  <- file.path(out_root, "dnmb", "processed")
fel_json <- file.path(fel_dir, "fel_demo_OG0000001.json")
dir.create(fel_dir, showWarnings = FALSE, recursive = TRUE)
{
  n_sites_j <- 120L
  alpha_v <- rexp(n_sites_j, 2); beta_v <- rexp(n_sites_j, 2)
  pos_ix  <- sample(n_sites_j, 6); neg_ix <- sample(setdiff(seq_len(n_sites_j), pos_ix), 8)
  beta_v[pos_ix]  <- alpha_v[pos_ix]  + rexp(6, 1) + 1.5
  alpha_v[neg_ix] <- beta_v[neg_ix]   + rexp(8, 1) + 1.5
  pval <- runif(n_sites_j, 0.01, 0.9)
  pval[pos_ix] <- runif(6, 0.001, 0.05)
  pval[neg_ix] <- runif(8, 0.001, 0.05)
  mat <- cbind(alpha_v, beta_v, abs(beta_v - alpha_v),
               rchisq(n_sites_j, 1), pval, runif(n_sites_j, 0.1, 2))
  jsonlite::write_json(
    list(analysis = list(info = "HyPhy FEL (synthetic demo)"),
         input = list(`number of sites` = n_sites_j,
                      `number of sequences` = N_GENOMES),
         MLE = list(
           headers = list(list("alpha","syn"), list("beta","non-syn"),
                          list("alpha=beta","LRT vs null"),
                          list("LRT","LRT"), list("p-value","p"),
                          list("Total branch length","t"))
           , content = list(`0` = mat))),
    fel_json, auto_unbox = TRUE, matrix = "rowmajor")
}
fel$json_path[1L] <- fel_json

selection <- list(busted = busted, fel = fel, absrel = absrel,
                  codeml = codeml)

cat("[demo] running run_dnmb_plot() + selection plots...\n")
plots_dir <- file.path(out_root, "plots")
run_dnmb_plot(out_root, output_dir = plots_dir, selection = selection)

cat("\n[demo] Output directory: ", plots_dir, "\n", sep = "")
pdf_files <- list.files(plots_dir, pattern = "\\.pdf$", full.names = FALSE)
cat("  PDFs generated (", length(pdf_files), "):\n", sep = "")
for (f in sort(pdf_files)) cat("    - ", f, "\n", sep = "")
