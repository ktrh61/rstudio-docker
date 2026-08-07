# reo_lowmid_outliers.R  (PROVISIONAL / 仮置き)
# Layer-1 QC for the REO out-of-sample validation: run PC-OD outlier detection
# on the R_Low / R_Mid tumours (which 530 currently uses unfiltered), then check
# whether the graded reversal-score relationship and the Mid > Low test survive
# removing any flagged sample. Provisional: not wired into the pipeline; mirrors
# 210_detect_outliers.R's PC-OD but on the intermediate-exposure bands.
# Input : processed/thyr_clinical.rds, thyr_case_assigned_share.rds,
#         thyr_se_raw.rds, thyr_reo_evaluation.rds
#         lib/qc_pc_od.R, lib/reo_lowmid_cases.R, lib/stat_brunnermunzel.R
# Output: processed/thyr_reo_lowmid_outliers_provisional.rds

source("setup.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
})
source(file.path(paths$root, "lib", "qc_pc_od.R"))
source(file.path(paths$root, "lib", "reo_lowmid_cases.R"))
source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))

# --- Load and resolve R_Low / R_Mid tumours --------------------------------
clinical <- readRDS(file.path(paths$processed, "thyr_clinical.rds"))
assigned_share <- as.data.frame(readRDS(file.path(paths$processed, "thyr_case_assigned_share.rds")))
se <- readRDS(file.path(paths$processed, "thyr_se_raw.rds"))
eval <- readRDS(file.path(paths$processed, "thyr_reo_evaluation.rds"))

cases <- resolve_ret_exposure_cases(clinical, assigned_share, se)
cases <- cases[cases$band %in% c("R_Low", "R_Mid") & !is.na(cases$tumor_id), , drop = FALSE]
message("R_Low/Mid tumours: ",
  paste(names(table(cases$band)), table(cases$band), sep = "=", collapse = " "))

# --- PC-OD per band on tumour log-CPM (same scale as 210) ------------------
counts_all <- assay(se, "stranded_second")
run_pcod <- function(ids) {
  y <- edgeR::DGEList(counts_all[, ids, drop = FALSE])
  keep <- edgeR::filterByExpr(y)
  logCPM <- edgeR::cpm(y[keep, ], log = TRUE, normalized.lib.sizes = FALSE, prior.count = 2)
  PC_OD(logCPM)
}
cases$is_outlier <- 0L
for (b in c("R_Low", "R_Mid")) {
  idx <- which(cases$band == b)
  out <- run_pcod(cases$tumor_id[idx])
  cases$is_outlier[idx[out]] <- 1L
  message("  ", b, ": ", length(out), " outlier(s)",
    if (length(out)) paste0(" (", paste(cases$case_submitter_id[idx[out]], collapse = ", "), ")") else "")
}

# --- Robustness: gradient and Mid > Low with vs without outliers -----------
score_of <- setNames(eval$samples$score, eval$samples$case_submitter_id)
cases$score <- unname(score_of[cases$case_submitter_id])
report <- function(keep, label) {
  d <- cases[keep, , drop = FALSE]
  lo <- d$score[d$band == "R_Low"]; mi <- d$score[d$band == "R_Mid"]
  bm <- brunnermunzel_mc_test(lo, mi, alternative = "less", method = "auto", seed = 19860426L)
  message(sprintf(
    "  %-16s Low median %.1f (n=%d) | Mid median %.1f (n=%d) | Mid>Low BM p=%.4f, P(Low<Mid)=%.3f",
    label, stats::median(lo), length(lo), stats::median(mi), length(mi),
    bm$p.value, unname(bm$estimate)))
}
message("Gradient robustness (out-of-sample Low vs Mid):")
report(rep(TRUE, nrow(cases)), "all samples")
report(cases$is_outlier == 0, "outliers removed")

# --- Save ------------------------------------------------------------------
out_rds <- file.path(paths$processed, "thyr_reo_lowmid_outliers_provisional.rds")
saveRDS(cases[, c("case_submitter_id", "band", "assigned_share", "tumor_id", "score", "is_outlier")], out_rds)
message("Saved: ", out_rds, " (", sum(cases$is_outlier), " outlier(s) flagged)")
