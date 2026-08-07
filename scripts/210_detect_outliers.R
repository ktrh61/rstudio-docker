# 210_detect_outliers.R
# Flag sample outliers with PC-OD on the analysis target cases (driver +
# Sporadic/High band + paired, taken from the case-design table), BEFORE
# tumor-purity estimation. Detecting outliers first keeps the downstream
# ContamDE purity estimate (220) from being contaminated by an anomalous
# sample.
# Input : processed/thyr_case_design.rds  (from 140; driver, band, pairs)
#         processed/thyr_se_raw.rds        (from 120; single count assay)
#         lib/qc_pc_od.R                   (PC_OD)
# Output: processed/thyr_case_outliers.rds
#           columns: case_submitter_id, group, driver, exposure,
#                    tumor_id, normal_id, has_outlier_tumor, has_outlier_normal
#
# PC-OD runs on eight group x tissue sub-matrices over ALL target cases (no
# purity filter). Each sub-matrix is filterByExpr-reduced and turned into
# log-CPM without normalization (normalized.lib.sizes = FALSE, prior.count = 2,
# the edgeR default), the natural pre-normalization scale for sample outlier
# detection: raw library sizes preserve composition outliers that normalization
# would mask, and the log scale keeps high-count genes from dominating.

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
})

source(file.path(paths$root, "lib", "qc_pc_od.R"))

# --- Load inputs -----------------------------------------------------------
se_path <- file.path(paths$processed, "thyr_se_raw.rds")
if (!file.exists(se_path)) stop("thyr_se_raw.rds not found (run 120 first)")
se <- readRDS(se_path)
message(
  "SE: ", nrow(se), " genes x ", ncol(se), " samples ; assay: ",
  assayNames(se)
)

design_path <- file.path(paths$processed, "thyr_case_design.rds")
if (!file.exists(design_path)) stop("thyr_case_design.rds not found (run 140 first)")
design <- readRDS(design_path)
message("Case design: ", nrow(design), " cases")

# --- Assemble target cases (adoption from the design table) ----------------
# Driver-classified AND band Sporadic/High AND paired.
targets <- design[
  !is.na(design$driver) &
    design$band %in% c("Sporadic", "High") &
    design$is_paired, ,
  drop = FALSE
]
targets <- data.frame(
  case_submitter_id = targets$case_submitter_id,
  driver = targets$driver,
  exposure = targets$band,
  tumor_id = targets$tumor_id,
  normal_id = targets$normal_id,
  stringsAsFactors = FALSE
)
prefix <- ifelse(targets$driver == "RET", "R", "B")
targets$group <- paste0(prefix, "_", targets$exposure)
message("Target cases by group:")
print(table(targets$group))

# --- PC-OD per (group x tissue) --------------------------------------------
counts_all <- assay(se)

run_pcod <- function(sample_ids) {
  cts <- counts_all[, sample_ids, drop = FALSE]
  y <- edgeR::DGEList(counts = cts)
  keep <- edgeR::filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = TRUE]
  logCPM <- edgeR::cpm(
    y,
    log = TRUE,
    normalized.lib.sizes = FALSE,
    prior.count = 2
  )
  PC_OD(logCPM)
}

targets$has_outlier_tumor <- 0L
targets$has_outlier_normal <- 0L

groups <- c("R_Sporadic", "R_High", "B_Sporadic", "B_High")
for (g in groups) {
  gi_idx <- which(targets$group == g)
  if (length(gi_idx) == 0) {
    message("Group ", g, ": no cases; skipped")
    next
  }
  gi <- targets[gi_idx, , drop = FALSE]
  message("Group ", g, ": ", nrow(gi), " cases")

  out_t <- run_pcod(gi$tumor_id)
  targets$has_outlier_tumor[gi_idx[out_t]] <- 1L
  message(
    "  Tumor : ", length(out_t), " outlier(s)",
    if (length(out_t)) paste0(" (", paste(gi$case_submitter_id[out_t], collapse = ", "), ")") else ""
  )

  out_n <- run_pcod(gi$normal_id)
  targets$has_outlier_normal[gi_idx[out_n]] <- 1L
  message(
    "  Normal: ", length(out_n), " outlier(s)",
    if (length(out_n)) paste0(" (", paste(gi$case_submitter_id[out_n], collapse = ", "), ")") else ""
  )
}

# --- Assemble and save -----------------------------------------------------
out_tbl <- targets[, c(
  "case_submitter_id", "group", "driver", "exposure",
  "tumor_id", "normal_id", "has_outlier_tumor", "has_outlier_normal"
), drop = FALSE]
out_tbl <- out_tbl[order(out_tbl$group, out_tbl$case_submitter_id), , drop = FALSE]
rownames(out_tbl) <- NULL

message("Cases with a tumor outlier : ", sum(out_tbl$has_outlier_tumor))
message("Cases with a normal outlier: ", sum(out_tbl$has_outlier_normal))

out <- file.path(paths$processed, "thyr_case_outliers.rds")
saveRDS(out_tbl, out)
message(
  "Saved: ", out, " (", nrow(out_tbl), " cases x ",
  ncol(out_tbl), " columns)"
)
