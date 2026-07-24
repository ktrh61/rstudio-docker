# 070_normalize_counts.R
# DEGES normalization of the clean paired samples from 060, per two-group
# comparison and tissue. The DEGES core (MUREN normalization + permutation
# Brunner-Munzel screening, TCC X-Y-X) lives in utils/deges_muren_bm.R; this
# script prepares each unit's count matrix and runs it.
# Input : processed/thyr_case_outliers.rds  (from 060)
#         processed/thyr_se_raw.rds         (from 021; single count assay)
#         utils/utils_improved.R, utils/norm_improved.R   (muren_norm)
#         utils/brunnermunzel_mc.R          (brunnermunzel_mc_test)
#         utils/deges_muren_bm.R  (deges_muren_bm, muren_to_norm_factors)
# Output: processed/thyr_normalized_counts.rds
#           list(date, config, comparisons, clean_cases, units)
#           units = { R_Tumor, R_Normal, B_Tumor, B_Normal }; each holds a
#           DGEList with DEGES-MUREN norm.factors and a diagnostics list.
#
# Units: R_Sporadic vs R_High and B_Sporadic vs B_High, each x {Tumor, Normal}.
# Per unit: protein_coding -> filterByExpr -> Cook's transient removal -> the
# DEGES core. The returned scaling coefficients are applied to the filterByExpr
# DGEList (before Cook's removal).

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(DESeq2)
})

source(file.path(paths$root, "utils", "utils_improved.R"))
source(file.path(paths$root, "utils", "norm_improved.R"))
source(file.path(paths$root, "utils", "brunnermunzel_mc.R"))
source(file.path(paths$root, "utils", "deges_muren_bm.R"))

# --- Configuration ---------------------------------------------------------
ITERATION <- 3L # DEGES iterations (exact count; iDEGES)
FDR <- 0.10 # BH adjusted-p cutoff for potential DEGs
FLOOR_PDEG <- 0.05 # floorPDEG fraction forced as potential DEGs (TCC)
COOKS_QUANTILE <- 0.99 # Cook's distance F-quantile cutoff
MUREN_METHOD <- "lts" # MUREN pairwise regression
WORKERS <- 16L # MUREN parallel workers

# Brunner-Munzel screening (utils/brunnermunzel_mc.R).
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
outliers_path <- file.path(paths$processed, "thyr_case_outliers.rds")
if (!file.exists(outliers_path)) {
  stop("thyr_case_outliers.rds not found (run 060 first)")
}
outliers <- readRDS(outliers_path)
message("Outlier table: ", nrow(outliers), " cases")

se_path <- file.path(paths$processed, "thyr_se_raw.rds")
if (!file.exists(se_path)) stop("thyr_se_raw.rds not found (run 021 first)")
se <- readRDS(se_path)
message(
  "SE: ", nrow(se), " genes x ", ncol(se), " samples ; assay: ",
  assayNames(se)
)

# --- Select clean cases ----------------------------------------------------
# A clean case has neither a tumor nor a normal outlier flag from 060.
clean <- outliers[
  outliers$has_outlier_tumor == 0 & outliers$has_outlier_normal == 0, ,
  drop = FALSE
]
message("Clean cases: ", nrow(clean), " / ", nrow(outliers))
message("Clean cases by group:")
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

# --- Cook's distance outlier genes (DESeq2) --------------------------------
# Transient removal for normalization only; the reported gene set is unchanged.
detect_cook_outliers <- function(count_matrix, sample_groups, quantile_cutoff) {
  coldata <- data.frame(
    group = factor(sample_groups),
    row.names = colnames(count_matrix)
  )
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix, colData = coldata, design = ~group
  )
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  cooks <- assays(dds)[["cooks"]]

  n_samples <- ncol(count_matrix)
  n_params <- ncol(model.matrix(~ factor(sample_groups)))
  f_cutoff <- qf(quantile_cutoff, df1 = n_params, df2 = n_samples - n_params)

  max_cooks <- apply(cooks, 1, function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  })
  outliers <- !is.na(max_cooks) & max_cooks > f_cutoff
  list(outliers = outliers, threshold = f_cutoff, outlier_count = sum(outliers))
}

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

  # Cook's distance transient removal (normalization only)
  cook <- detect_cook_outliers(
    count_matrix_filtered, sample_groups, COOKS_QUANTILE
  )
  count_for_norm <- if (cook$outlier_count > 0) {
    count_matrix_filtered[!cook$outliers, , drop = FALSE]
  } else {
    count_matrix_filtered
  }
  message("    Cook's outlier genes: ", cook$outlier_count)

  # TCC-faithful DEGES core (MUREN + MC Brunner-Munzel) on the Cook-removed set.
  deges <- deges_muren_bm(
    counts = count_for_norm, group = sample_groups,
    iteration = ITERATION, fdr = FDR, floor_pdeg = FLOOR_PDEG,
    n_perm = BM_MC_B, seed = BM_MC_SEED,
    alternative = BM_MC_ALTERNATIVE, est = BM_MC_EST,
    bm_method = BM_METHOD,
    muren_method = MUREN_METHOD, workers = WORKERS
  )

  # Apply the final scaling coefficients to the filterByExpr DGEList (before
  # Cook's removal); carry norm.factors and the MUREN coefficients on samples.
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
  message(sprintf(
    "    norm.factors range: [%.4f, %.4f]", min(nf_final), max(nf_final)
  ))

  list(
    dgelist = dge_final,
    diagnostics = list(
      groups = group_labels,
      n_samples = c(length(samples1), length(samples2)),
      n_protein_coding = nrow(count_matrix),
      n_after_filter = nrow(count_matrix_filtered),
      cook_outliers = cook$outlier_count,
      cook_threshold = cook$threshold,
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
    cooks_quantile = COOKS_QUANTILE,
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
