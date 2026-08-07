# reo_lowmid_purity.R
# Layer-1 QC for the REO out-of-sample validation: estimate R_Low / R_Mid tumour
# purity on the SAME common scale as the training arms by pooling the whole RET
# cohort (Sporadic + Low + Mid + High) in one ContamDE run, then check whether
# the REO reversal score is driven by purity within a band. If score is
# independent of purity, purity filtering the validation set would not change the
# conclusion; if correlated, a full re-filter (deferred layer 2) is warranted.
# A diagnostic outside the numbered stream (reorg plan v2 s2.6); mirrors 220
# pooled ContamDE. Run after the main chain.
# Input : processed/thyr_case_design.rds (from 140),
#         thyr_se_raw.rds, thyr_reo_evaluation.rds
#         lib/norm_muren_helpers.R, lib/norm_muren.R,
#         lib/purity_contamde.R
# Output: diagnostics/output/reo_lowmid_purity.rds

source("setup.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
})
source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "norm_muren_helpers.R"))
source(file.path(paths$root, "lib", "norm_muren.R"))
source(file.path(paths$root, "lib", "purity_contamde.R"))
source(file.path(paths$root, "lib", "gene_filter.R"))

# WORKERS comes from config.R via setup.R (development 16L; canonical is 4L).
pin_blas_threads()

# --- Resolve the whole RET cohort (all four bands) from the design table ----
design <- readRDS(file.path(paths$processed, "thyr_case_design.rds"))
se <- readRDS(file.path(paths$processed, "thyr_se_raw.rds"))
eval <- readRDS(file.path(paths$processed, "thyr_reo_evaluation.rds"))

ret <- design[design$driver %in% "RET" & !is.na(design$band), , drop = FALSE]
cases <- data.frame(
  case_submitter_id = ret$case_submitter_id,
  dose_mgy = ret$dose_mgy,
  assigned_share = ret$assigned_share_approx,
  band = paste0("R_", ret$band),
  tumor_id = ret$tumor_id,
  normal_id = ret$normal_id,
  stringsAsFactors = FALSE
)
cases <- cases[!is.na(cases$tumor_id) & !is.na(cases$normal_id), , drop = FALSE]
message("RET cohort with pairs: ",
  paste(names(table(cases$band)), table(cases$band), sep = "=", collapse = " "))

# --- Pooled ContamDE on the whole RET cohort (common scale) ----------------
gene_info <- as.data.frame(rowData(se))
is_pc <- gene_info$gene_type == "protein_coding"
counts_all <- se_single_assay(se)
keep <- filter_pc_expr_mask(counts_all, is_pc, cases$normal_id, cases$tumor_id)
np <- nrow(cases)
counts <- cbind(counts_all[keep, cases$normal_id, drop = FALSE],
  counts_all[keep, cases$tumor_id, drop = FALSE])
colnames(counts) <- c(paste0("Normal_", seq_len(np)), paste0("Tumor_", seq_len(np)))
res <- contamde_purity(counts = counts, subtype = NULL, covariate = NULL,
  contaminated = TRUE, pairwise_method = "lts", workers = WORKERS, verbose = FALSE)
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

out_dir <- file.path(paths$root, "diagnostics", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_rds <- file.path(out_dir, "reo_lowmid_purity.rds")
saveRDS(cases[, c("case_submitter_id", "band", "assigned_share", "tumor_purity", "score")], out_rds)
message("Saved: ", out_rds)
