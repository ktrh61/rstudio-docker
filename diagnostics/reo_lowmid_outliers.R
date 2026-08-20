# reo_lowmid_outliers.R
# Ancillary QC for the REO intermediate-band application: apply the PC-OD
# procedure to R_Low / R_Mid tumours, record any flags without excluding cases,
# and descriptively recompute the graded comparison after removing flagged
# samples. This diagnostic does not redefine the application cohort. It runs
# outside the numbered stream and mirrors 210_detect_outliers.R's PC-OD on the
# intermediate-exposure bands.
# Input : processed/thyr_case_design.rds (from 140),
#         thyr_se_raw.rds, thyr_reo_evaluation.rds
#         lib/qc_pc_od.R, lib/stat_brunnermunzel.R
# Output: diagnostics/output/reo_lowmid_outliers.rds

source("setup.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
})
source(file.path(paths$root, "lib", "qc_pc_od.R"))
source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "gene_filter.R"))
source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))

# --- Load and resolve R_Low / R_Mid tumours (from the design table) --------
design <- readRDS(file.path(paths$processed, "thyr_case_design.rds"))
se <- readRDS(file.path(paths$processed, "thyr_se_raw.rds"))
eval <- readRDS(file.path(paths$processed, "thyr_reo_evaluation.rds"))

ret <- design[design$driver %in% "RET" & !is.na(design$band), , drop = FALSE]
cases <- data.frame(
  case_submitter_id = ret$case_submitter_id,
  dose_mgy = ret$dose_mgy,
  assigned_share = ret$assigned_share,
  band = paste0("R_", ret$band),
  tumor_id = ret$tumor_id,
  normal_id = ret$normal_id,
  stringsAsFactors = FALSE
)
cases <- cases[cases$band %in% c("R_Low", "R_Mid") & !is.na(cases$tumor_id), , drop = FALSE]
message("R_Low/Mid tumours: ",
  paste(names(table(cases$band)), table(cases$band), sep = "=", collapse = " "))

# --- PC-OD per band on tumour log-CPM (same scale as 210) ------------------
counts_all <- se_single_assay(se)
run_pcod <- function(ids) {
  # keep_lib_sizes = FALSE is this diagnostic's historical behaviour (210 uses
  # TRUE); the divergence is preserved, unification is a recorded open item.
  PC_OD(logcpm_for_qc(counts_all, ids, keep_lib_sizes = FALSE))
}
cases$is_outlier <- 0L
for (b in c("R_Low", "R_Mid")) {
  idx <- which(cases$band == b)
  out <- run_pcod(cases$tumor_id[idx])
  cases$is_outlier[idx[out]] <- 1L
  message("  ", b, ": ", length(out), " outlier(s)",
    if (length(out)) paste0(" (", paste(cases$case_submitter_id[idx[out]], collapse = ", "), ")") else "")
}

# --- Descriptive comparison with and without flagged samples ---------------
score_of <- setNames(eval$samples$score, eval$samples$case_submitter_id)
cases$score <- unname(score_of[cases$case_submitter_id])
report <- function(keep, label) {
  d <- cases[keep, , drop = FALSE]
  lo <- d$score[d$band == "R_Low"]; mi <- d$score[d$band == "R_Mid"]
  bm <- brunnermunzel_mc_test(lo, mi, alternative = "less", method = "auto", seed = SEED)
  message(sprintf(
    "  %-16s Low median %.1f (n=%d) | Mid median %.1f (n=%d) | Mid>Low BM p=%.4f, Pr(Low<Mid)=%.3f",
    label, stats::median(lo), length(lo), stats::median(mi), length(mi),
    bm$p.value, unname(bm$estimate)))
}
message("Diagnostic Low vs Mid comparison with and without flagged samples:")
report(rep(TRUE, nrow(cases)), "all samples")
report(cases$is_outlier == 0, "outliers removed")

# --- Save ------------------------------------------------------------------
out_dir <- file.path(paths$root, "diagnostics", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_rds <- file.path(out_dir, "reo_lowmid_outliers.rds")
saveRDS(cases[, c("case_submitter_id", "band", "assigned_share", "tumor_id", "score", "is_outlier")], out_rds)
message("Saved: ", out_rds, " (", sum(cases$is_outlier), " outlier(s) flagged)")
