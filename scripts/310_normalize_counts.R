# 310_normalize_counts.R
# DEGES normalization of the high-purity paired samples, per two-group
# comparison and tissue. Cases are the outlier-removed (210), pooled-purity
# (220) cases whose relative purity clears PURITY_THRESHOLD. The DEGES core
# (MUREN normalization + permutation Brunner-Munzel screening, TCC X-Y-X) lives
# in lib/norm_deges.R; this script prepares each unit's count matrix.
# Input : processed/thyr_case_purity.rds    (from 220; case, group, tumor_purity)
#         processed/thyr_se_raw.rds         (from 120; single count assay)
#         lib/norm_muren_helpers.R, lib/norm_muren.R   (muren_norm)
#         lib/stat_brunnermunzel.R          (brunnermunzel_mc_test)
#         lib/norm_deges.R  (deges_muren_bm, muren_to_norm_factors)
# Output: processed/thyr_normalized_counts.rds
#           list(date, config, comparisons, clean_cases, units)
#           units = { R_Tumor, R_Normal, B_Tumor, B_Normal }; each holds a
#           DGEList with DEGES-MUREN norm.factors and a diagnostics list.
#
# Units: R_Sporadic vs R_High and B_Sporadic vs B_High, each x {Tumor, Normal}.
# Per unit: protein_coding -> filterByExpr -> the DEGES core. The returned
# scaling coefficients are applied to the filterByExpr DGEList.

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
})

source(file.path(paths$root, "lib", "norm_muren_helpers.R"))
source(file.path(paths$root, "lib", "norm_muren.R"))
source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))
source(file.path(paths$root, "lib", "norm_deges.R"))

# --- Configuration ---------------------------------------------------------
PURITY_THRESHOLD <- 0.6 # min pooled (common-scale) relative purity to include
ITERATION <- 3L # DEGES iterations (exact count; iDEGES)
FDR <- 0.10 # BH adjusted-p cutoff for potential DEGs
FLOOR_PDEG <- 0.05 # floorPDEG fraction forced as potential DEGs (TCC)
MUREN_METHOD <- "lts" # MUREN pairwise regression
WORKERS <- 16L # MUREN parallel workers

# Brunner-Munzel screening (lib/stat_brunnermunzel.R).
# "auto" enumerates all C(n, nx) allocations exactly when that is affordable
# and samples otherwise. Both units here are well inside the budget
# (C(28,12) = 3.0e7 and C(36,9) = 9.4e7), so the screen's p-values are exact:
# no 1/(B + 1) floor, no sampling error, and no dependence on BM_MC_SEED.
# BM_MC_B / BM_MC_SEED apply only if a unit ever falls back to sampling.
BM_METHOD <- "auto"
BM_EXACT_MAX <- 1e8 # largest C(n, nx) still enumerated exactly
BM_EXACT_THREADS <- 16L # threads for the enumeration
BM_MC_B <- 999999L # Monte Carlo permutations (fallback only)
BM_MC_SEED <- 19860426L # fixed seed (Chernobyl accident date 1986-04-26)
BM_MC_ALTERNATIVE <- "two.sided"
BM_MC_EST <- "original"

options(
  brunnermunzel.exact.max.allocations = BM_EXACT_MAX,
  brunnermunzel.exact.threads = BM_EXACT_THREADS
)

# BLAS/OMP single-threaded (MUREN parallelizes the outer loop).
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}

# --- Load inputs -----------------------------------------------------------
purity_path <- file.path(paths$processed, "thyr_case_purity.rds")
if (!file.exists(purity_path)) {
  stop("thyr_case_purity.rds not found (run 220 first)")
}
purity <- readRDS(purity_path)
message("Purity table: ", nrow(purity), " cases (outlier-removed by 210)")

se_path <- file.path(paths$processed, "thyr_se_raw.rds")
if (!file.exists(se_path)) stop("thyr_se_raw.rds not found (run 120 first)")
se <- readRDS(se_path)
message(
  "SE: ", nrow(se), " genes x ", ncol(se), " samples ; assay: ",
  assayNames(se)
)

# --- Select high-purity cases ----------------------------------------------
# Outliers were already removed by 210 before purity estimation; here we keep
# cases whose pooled (common-scale) relative purity clears the threshold. The
# threshold lives here so a sensitivity check re-runs only 310, not the 220
# ContamDE step.
clean <- purity[purity$tumor_purity >= PURITY_THRESHOLD, , drop = FALSE]
message(sprintf(
  "High-purity cases (>= %.2f): %d / %d", PURITY_THRESHOLD, nrow(clean), nrow(purity)
))
message("Selected cases by group:")
print(table(clean$group))

# --- Resolve _merged Tumor / Normal sample per case ------------------------
cd <- as.data.frame(colData(se))
is_merged <- grepl("_merged", cd$sample_submitter_id)
m <- cd[is_merged, , drop = FALSE]
t_rows <- m[m$sample_type == "Primary Tumor", , drop = FALSE]
n_rows <- m[m$sample_type == "Solid Tissue Normal", , drop = FALSE]
tumor_of <- setNames(t_rows$sample_submitter_id, t_rows$case_submitter_id)
normal_of <- setNames(n_rows$sample_submitter_id, n_rows$case_submitter_id)

sample_ids_for <- function(cases, tissue) {
  ids <- if (tissue == "Tumor") tumor_of[cases] else normal_of[cases]
  if (any(is.na(ids))) {
    stop("A clean case is missing its _merged ", tissue, " sample.")
  }
  unname(ids)
}

# --- Gene annotation and counts --------------------------------------------
gene_info <- as.data.frame(rowData(se))
is_protein_coding <- gene_info$gene_type == "protein_coding"
counts_all <- assay(se)

# --- Normalize one unit (comparison x tissue) ------------------------------
process_unit <- function(samples1, samples2, group_labels) {
  sample_groups <- c(
    rep(group_labels[1], length(samples1)),
    rep(group_labels[2], length(samples2))
  )
  all_samples <- c(samples1, samples2)

  # protein_coding -> filterByExpr (reference gene set for this unit)
  count_matrix <- counts_all[is_protein_coding, all_samples, drop = FALSE]
  gene_info_pc <- gene_info[is_protein_coding, , drop = FALSE]
  dge_initial <- DGEList(counts = count_matrix, group = factor(sample_groups))
  keep <- filterByExpr(dge_initial, group = dge_initial$samples$group)
  count_matrix_filtered <- count_matrix[keep, , drop = FALSE]
  gene_info_filtered <- gene_info_pc[keep, , drop = FALSE]
  message(
    "    protein_coding: ", nrow(count_matrix),
    " ; after filterByExpr: ", nrow(count_matrix_filtered)
  )

  # TCC-faithful DEGES core (MUREN + Brunner-Munzel) on the filterByExpr set.
  deges <- deges_muren_bm(
    counts = count_matrix_filtered, group = sample_groups,
    iteration = ITERATION, fdr = FDR, floor_pdeg = FLOOR_PDEG,
    n_perm = BM_MC_B, seed = BM_MC_SEED,
    alternative = BM_MC_ALTERNATIVE, est = BM_MC_EST,
    bm_method = BM_METHOD,
    muren_method = MUREN_METHOD, workers = WORKERS
  )

  # Apply the final scaling coefficients to the filterByExpr DGEList; carry
  # norm.factors and the MUREN coefficients on samples.
  dge_final <- DGEList(
    counts = count_matrix_filtered,
    group = factor(sample_groups),
    genes = gene_info_filtered
  )
  nf_final <- muren_to_norm_factors(
    deges$scaling_coeff[colnames(dge_final)], dge_final$samples$lib.size
  )
  dge_final$samples$norm.factors <- nf_final
  dge_final$samples$scaling_coeff <-
    deges$scaling_coeff[colnames(dge_final)]
  # The saved norm.factors are the audit trail for the MUREN coefficients, so
  # they must stay exactly what those coefficients imply. Setting one without
  # the other would otherwise be silent.
  stopifnot(isTRUE(all.equal(
    nf_final,
    muren_to_norm_factors(
      dge_final$samples$scaling_coeff, dge_final$samples$lib.size
    )
  )))
  message(sprintf(
    "    norm.factors range: [%.4f, %.4f] ; scaling_coeff range: [%.4f, %.4f]",
    min(nf_final), max(nf_final),
    min(dge_final$samples$scaling_coeff), max(dge_final$samples$scaling_coeff)
  ))

  list(
    dgelist = dge_final,
    diagnostics = list(
      groups = group_labels,
      n_samples = c(length(samples1), length(samples2)),
      n_protein_coding = nrow(count_matrix),
      n_after_filter = nrow(count_matrix_filtered),
      n_iterations = deges$iterations$n,
      iter_n_genes = deges$iterations$n_non_deg,
      iter_deg_count = deges$iterations$deg_count,
      iter_method = deges$iterations$exclusion_method,
      final_deg_count = deges$iterations$deg_count[deges$iterations$n]
    )
  )
}

# --- Process each comparison x tissue --------------------------------------
comparisons <- list(
  R = c("R_Sporadic", "R_High"),
  B = c("B_Sporadic", "B_High")
)
tissues <- c("Tumor", "Normal")

units <- list()
for (comp_key in names(comparisons)) {
  grp <- comparisons[[comp_key]]
  cases1 <- clean$case_submitter_id[clean$group == grp[1]]
  cases2 <- clean$case_submitter_id[clean$group == grp[2]]
  message(
    "Comparison ", comp_key, " (", grp[1], " vs ", grp[2], "): ",
    length(cases1), " + ", length(cases2), " cases"
  )
  if (length(cases1) == 0 || length(cases2) == 0) {
    stop("Comparison ", comp_key, " has an empty group after cleaning.")
  }
  for (tissue in tissues) {
    unit <- paste0(comp_key, "_", tissue)
    message("  Unit ", unit)
    s1 <- sample_ids_for(cases1, tissue)
    s2 <- sample_ids_for(cases2, tissue)
    units[[unit]] <- process_unit(s1, s2, grp)
  }
}

# --- Assemble and save -----------------------------------------------------
thyr_normalized_counts <- list(
  date = Sys.Date(),
  config = list(
    iteration = ITERATION,
    fdr = FDR,
    floor_pdeg = FLOOR_PDEG,
    muren_method = MUREN_METHOD,
    workers = WORKERS,
    bm_test = "permutation_brunnermunzel",
    bm_method = BM_METHOD,
    bm_exact_max = BM_EXACT_MAX,
    bm_B = BM_MC_B,
    bm_seed = BM_MC_SEED,
    bm_alternative = BM_MC_ALTERNATIVE,
    bm_est = BM_MC_EST
  ),
  comparisons = comparisons,
  clean_cases = clean,
  units = units
)

out <- file.path(paths$processed, "thyr_normalized_counts.rds")
saveRDS(thyr_normalized_counts, out)
message(
  "Saved: ", out, " (", length(units), " units: ",
  paste(names(units), collapse = ", "), ")"
)
