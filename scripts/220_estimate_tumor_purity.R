# 060_estimate_tumor_purity.R
# Estimate relative tumor purity per case with ContamDE (MUREN normalization),
# AFTER outlier detection (050) and with the two exposure arms of each cohort
# POOLED. Pooling puts High and Sporadic on one common relative scale (ContamDE
# is comparable only within the set of pairs passed together) and adds pairs for
# a steadier reference profile. It is licensed here because the two arms share
# the driver's tumour transcriptome (RET/PTC or BRAF) and the exposure effect is
# weak, so a common tumour reference holds; verified on this data (within-group
# purity rank preserved 0.93-0.99 vs per-arm estimation; High vs Sporadic tumour
# profile correlation 0.99). The purity axis (tumour vs normal contamination) is
# orthogonal to the exposure sample-mixture axis, so pooling does not dilute the
# exposure signal.
# Input : processed/thyr_case_outliers.rds  (from 050; targets + outlier flags)
#         processed/thyr_se_raw.rds          (from 021; single count assay)
#         utils/utils_improved.R, utils/norm_improved.R,
#         utils/contamde_purity_functions.R
# Output: processed/thyr_case_purity.rds
#           columns: case_submitter_id, group, cohort, tumor_purity
#
# Outlier cases (either tissue flagged by 050) are dropped before estimation and
# do not appear in the output. The relative-purity score (max-one within a
# cohort) is emitted for every retained case at every purity; the >= threshold
# selection is applied downstream in 070, so a threshold sensitivity check needs
# only 070 re-run, not this ContamDE step.

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
})

source(file.path(paths$root, "utils", "utils_improved.R"))
source(file.path(paths$root, "utils", "norm_improved.R"))
source(file.path(paths$root, "utils", "contamde_purity_functions.R"))

# MUREN worker count. The published/reproduction run on the Xeon host uses 4L
# (canonical); raised here only to speed up development iterations. Worker count
# changes speed, not the MUREN result.
workers <- 16L
message("MUREN workers: ", workers, " (development; canonical is 4L)")

# BLAS/OMP single-threaded (MUREN parallelizes the outer loop).
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}

# --- Load inputs -----------------------------------------------------------
outliers_path <- file.path(paths$processed, "thyr_case_outliers.rds")
if (!file.exists(outliers_path)) stop("thyr_case_outliers.rds not found (run 050 first)")
targets <- readRDS(outliers_path)
message("Target cases: ", nrow(targets))

se_path <- file.path(paths$processed, "thyr_se_raw.rds")
if (!file.exists(se_path)) stop("thyr_se_raw.rds not found (run 021 first)")
se <- readRDS(se_path)

# --- Drop outlier cases (order correction: clean before purity) ------------
is_clean <- targets$has_outlier_tumor == 0 & targets$has_outlier_normal == 0
message(
  "Dropping ", sum(!is_clean), " outlier case(s); ",
  sum(is_clean), " retained for purity estimation"
)
targets <- targets[is_clean, , drop = FALSE]
targets$cohort <- ifelse(targets$driver == "RET", "R", "B")

# --- Gene filter (protein_coding + filterByExpr) ---------------------------
gene_info <- as.data.frame(rowData(se))
is_protein_coding <- gene_info$gene_type == "protein_coding"
counts_all <- assay(se)

filter_genes <- function(normal_ids, tumor_ids) {
  cg <- counts_all[, c(normal_ids, tumor_ids), drop = FALSE]
  grp <- factor(c(
    rep("Normal", length(normal_ids)),
    rep("Tumor", length(tumor_ids))
  ))
  keep_expr <- edgeR::filterByExpr(cg, group = grp)
  is_protein_coding & keep_expr
}

# --- Per-cohort pooled purity estimation -----------------------------------
purity_by_cohort <- list()

for (coh in c("R", "B")) {
  gi <- targets[targets$cohort == coh, , drop = FALSE]
  if (nrow(gi) == 0) {
    message("Cohort ", coh, ": no cases; skipped")
    next
  }
  message(
    "Cohort ", coh, ": ", nrow(gi), " pairs (",
    paste(names(table(gi$group)), table(gi$group), sep = "=", collapse = " "), ")"
  )

  keep <- filter_genes(gi$normal_id, gi$tumor_id)
  np <- nrow(gi)
  counts <- cbind(
    counts_all[keep, gi$normal_id, drop = FALSE],
    counts_all[keep, gi$tumor_id, drop = FALSE]
  )
  colnames(counts) <- c(
    paste0("Normal_", seq_len(np)),
    paste0("Tumor_", seq_len(np))
  )
  message(
    "  genes after filter: ", sum(keep), " ; matrix ",
    nrow(counts), " x ", ncol(counts)
  )

  res <- contamde_purity(
    counts = counts,
    subtype = NULL,
    covariate = NULL,
    contaminated = TRUE,
    pairwise_method = "lts",
    workers = workers,
    verbose = FALSE
  )

  purity_by_cohort[[coh]] <- data.frame(
    case_submitter_id = gi$case_submitter_id,
    group = gi$group,
    cohort = coh,
    tumor_purity = as.numeric(res$proportion),
    stringsAsFactors = FALSE
  )
  message(sprintf(
    "  purity (common scale): mean=%.3f sd=%.3f range=[%.3f, %.3f]",
    mean(res$proportion), sd(res$proportion),
    min(res$proportion), max(res$proportion)
  ))
}

# --- Assemble and save -----------------------------------------------------
thyr_case_purity <- do.call(rbind, purity_by_cohort)
rownames(thyr_case_purity) <- NULL
thyr_case_purity <- thyr_case_purity[
  order(thyr_case_purity$group, thyr_case_purity$case_submitter_id), ,
  drop = FALSE
]

out <- file.path(paths$processed, "thyr_case_purity.rds")
saveRDS(thyr_case_purity, out)
message(
  "Saved: ", out, " (", nrow(thyr_case_purity), " cases x ",
  ncol(thyr_case_purity), " columns)"
)
