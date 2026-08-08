# 310_normalize_counts.R
# DEGES normalization of the high-purity paired samples, per two-group
# comparison and tissue. Case adoption (outlier-removed, purity-passing) is
# fixed upstream in 230 (include_main_bm). The DEGES core
# (MUREN normalization + permutation Brunner-Munzel screening, TCC X-Y-X) lives
# in lib/norm_deges.R; this script prepares each unit's count matrix.
# Input : processed/thyr_analysis_cohorts.rds (from 230; include_main_bm)
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
source(file.path(paths$root, "lib", "stat_storey.R"))
source(file.path(paths$root, "lib", "norm_deges.R"))

# --- Configuration ---------------------------------------------------------
# Shared constants (PURITY_THRESHOLD, DEGES_FDR, WORKERS, SEED, EXACT_THREADS,
# BM_EXACT_MAX) come from config.R via setup.R.
ITERATION <- 3L # DEGES iterations (exact count; iDEGES)
FLOOR_PDEG <- 0.05 # floorPDEG fraction forced as potential DEGs (TCC)
MUREN_METHOD <- "lts" # MUREN pairwise regression

# Brunner-Munzel screening (lib/stat_brunnermunzel.R). Declared exact: both
# units are well inside the enumeration budget (C(28,12) = 3.0e7 and
# C(36,9) = 9.4e7), so the screen's p-values carry no 1/(B + 1) floor, no
# sampling error, and no seed dependence. The declaration matches the path
# actually taken (reorg plan v2 s3.2): if a future cohort ever made
# enumeration impossible the run stops rather than silently sampling.
BM_METHOD <- "exact"
BM_ALTERNATIVE <- "two.sided"

options(
  brunnermunzel.exact.max.allocations = BM_EXACT_MAX,
  brunnermunzel.exact.threads = EXACT_THREADS
)

# BLAS/OMP single-threaded (MUREN parallelizes the outer loop).
pin_blas_threads()

# --- Load inputs -----------------------------------------------------------
cohorts_path <- file.path(paths$processed, "thyr_analysis_cohorts.rds")
if (!file.exists(cohorts_path)) {
  stop("thyr_analysis_cohorts.rds not found (run 230 first)")
}
cohorts <- readRDS(cohorts_path)

se_path <- file.path(paths$processed, "thyr_se_raw.rds")
if (!file.exists(se_path)) stop("thyr_se_raw.rds not found (run 120 first)")
se <- readRDS(se_path)
message(
  "SE: ", nrow(se), " genes x ", ncol(se), " samples ; assay: ",
  assayNames(se)
)

# --- Select the main BM cohort (adoption fixed in 230) ----------------------
# Purity-threshold application, PC-OD flags and pairing all live in 230; a
# purity-threshold sensitivity check re-runs 230 + 310.
included <- cohorts[cohorts$include_main_bm, , drop = FALSE]
clean <- data.frame(
  case_submitter_id = included$case_submitter_id,
  group = paste0(
    ifelse(included$driver == "RET", "R", "B"), "_", included$band
  ),
  cohort = ifelse(included$driver == "RET", "R", "B"),
  tumor_purity = included$tumor_purity,
  stringsAsFactors = FALSE
)
clean <- clean[order(clean$group, clean$case_submitter_id), , drop = FALSE]
rownames(clean) <- NULL
message("Main BM cohort (include_main_bm): ", nrow(clean), " cases")
message("Selected cases by group:")
print(table(clean$group))

tumor_of <- setNames(included$tumor_id, included$case_submitter_id)
normal_of <- setNames(included$normal_id, included$case_submitter_id)

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
    iteration = ITERATION, fdr = DEGES_FDR, floor_pdeg = FLOOR_PDEG,
    n_perm = 0L, seed = SEED, # both unused on the exact path
    alternative = BM_ALTERNATIVE,
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
      iter_pi0 = deges$iterations$pi0,
      iter_removed_frac = deges$iterations$removed_frac,
      iter_jaccard = deges$iterations$jaccard,
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
    fdr = DEGES_FDR,
    fdr_method = "storey_plugin_lambda0.5",
    floor_pdeg = FLOOR_PDEG,
    muren_method = MUREN_METHOD,
    workers = WORKERS,
    bm_test = "permutation_brunnermunzel",
    bm_method = BM_METHOD,
    bm_exact_max = BM_EXACT_MAX,
    bm_alternative = BM_ALTERNATIVE
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
