# reo_lowmid_purity.R  (PROVISIONAL / 仮置き)
# Layer-1 QC for the REO out-of-sample validation: estimate R_Low / R_Mid tumour
# purity on the SAME common scale as the training arms by pooling the whole RET
# cohort (Sporadic + Low + Mid + High) in one ContamDE run, then check whether
# the REO reversal score is driven by purity within a band. If score is
# independent of purity, purity filtering the validation set would not change the
# conclusion; if correlated, a full re-filter (deferred layer 2) is warranted.
# Provisional: not wired into the pipeline; mirrors 220 pooled ContamDE.
# Input : processed/thyr_clinical.rds, thyr_case_assigned_share.rds,
#         thyr_se_raw.rds, thyr_reo_evaluation.rds
#         lib/norm_muren_helpers.R, lib/norm_muren.R,
#         lib/purity_contamde.R, lib/reo_lowmid_cases.R
# Output: processed/thyr_reo_lowmid_purity_provisional.rds

source("setup.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
})
source(file.path(paths$root, "lib", "norm_muren_helpers.R"))
source(file.path(paths$root, "lib", "norm_muren.R"))
source(file.path(paths$root, "lib", "purity_contamde.R"))
source(file.path(paths$root, "lib", "reo_lowmid_cases.R"))

workers <- 16L # development; canonical is 4L
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}

# --- Resolve the whole RET cohort (all four bands) with normal pairs --------
clinical <- readRDS(file.path(paths$processed, "thyr_clinical.rds"))
assigned_share <- as.data.frame(readRDS(file.path(paths$processed, "thyr_case_assigned_share.rds")))
se <- readRDS(file.path(paths$processed, "thyr_se_raw.rds"))
eval <- readRDS(file.path(paths$processed, "thyr_reo_evaluation.rds"))

cases <- resolve_ret_exposure_cases(clinical, assigned_share, se)
cases <- cases[!is.na(cases$tumor_id) & !is.na(cases$normal_id), , drop = FALSE]
message("RET cohort with pairs: ",
  paste(names(table(cases$band)), table(cases$band), sep = "=", collapse = " "))

# --- Pooled ContamDE on the whole RET cohort (common scale) ----------------
gene_info <- as.data.frame(rowData(se))
is_pc <- gene_info$gene_type == "protein_coding"
counts_all <- assay(se, "stranded_second")
grp <- factor(c(rep("Normal", nrow(cases)), rep("Tumor", nrow(cases))))
keep <- is_pc & edgeR::filterByExpr(
  counts_all[, c(cases$normal_id, cases$tumor_id), drop = FALSE], group = grp
)
np <- nrow(cases)
counts <- cbind(counts_all[keep, cases$normal_id, drop = FALSE],
  counts_all[keep, cases$tumor_id, drop = FALSE])
colnames(counts) <- c(paste0("Normal_", seq_len(np)), paste0("Tumor_", seq_len(np)))
res <- contamde_purity(counts = counts, subtype = NULL, covariate = NULL,
  contaminated = TRUE, pairwise_method = "lts", workers = workers, verbose = FALSE)
cases$tumor_purity <- as.numeric(res$proportion)

message("Common-scale purity by band (median):")
for (b in c("R_Sporadic", "R_Low", "R_Mid", "R_High")) {
  w <- cases$tumor_purity[cases$band == b]
  if (length(w)) message(sprintf("  %-11s n=%2d median %.3f range [%.3f, %.3f]",
    b, length(w), stats::median(w), min(w), max(w)))
}

# --- Diagnostic: does purity drive the REO reversal score within Low/Mid? --
score_of <- setNames(eval$samples$score, eval$samples$case_submitter_id)
cases$score <- unname(score_of[cases$case_submitter_id])
lm <- cases[cases$band %in% c("R_Low", "R_Mid") & !is.na(cases$score), , drop = FALSE]
message("\nScore vs purity within R_Low/Mid (if ~0, purity does not drive the validation):")
for (b in c("R_Low", "R_Mid")) {
  d <- lm[lm$band == b, ]
  if (nrow(d) >= 4) {
    r <- suppressWarnings(cor(d$tumor_purity, d$score, method = "spearman"))
    message(sprintf("  %-6s Spearman(purity, score) = %+.3f (n=%d)", b, r, nrow(d)))
  }
}
r_all <- suppressWarnings(cor(lm$tumor_purity, lm$score, method = "spearman"))
message(sprintf("  pooled Low+Mid Spearman(purity, score) = %+.3f (n=%d)", r_all, nrow(lm)))

out_rds <- file.path(paths$processed, "thyr_reo_lowmid_purity_provisional.rds")
saveRDS(cases[, c("case_submitter_id", "band", "assigned_share", "tumor_purity", "score")], out_rds)
message("Saved: ", out_rds)
