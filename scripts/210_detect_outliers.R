# 210_detect_outliers.R
# Assemble the analysis target cases (driver + exposure + paired) and flag
# sample outliers with PC-OD, BEFORE tumor-purity estimation. Detecting
# outliers first keeps the downstream ContamDE purity estimate (220) from being
# contaminated by an anomalous sample.
# Input : processed/thyr_clinical.rds             (from 030; driver columns)
#         processed/thyr_case_assigned_share.rds  (from 130; dose_mgy, AS)
#         processed/thyr_se_raw.rds                (from 120; single count assay)
#         lib/qc_pc_od.R                   (PC_OD)
# Output: processed/thyr_case_outliers.rds
#           columns: case_submitter_id, group, driver, exposure,
#                    tumor_id, normal_id, has_outlier_tumor, has_outlier_normal
#
# Group definition (this dataset):
#   driver  : BRAF  = Designated_Driver "BRAF.MutV600E" with no co-mutation
#             RET   = Designated_Driver in {CCDC6-RET, NCOA4-RET, RET-OTHER}
#   exposure: dose_mgy == 0 -> Sporadic ; dose_mgy > 0 with AS >= 66.6 -> High.
#   group   = {R,B}_{Sporadic,High} ; R = RET driver, B = BRAF driver.
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

# AS threshold (percent) separating High from Low/Mid among exposed cases.
as_high_threshold <- 66.6

# --- Load inputs -----------------------------------------------------------
se_path <- file.path(paths$processed, "thyr_se_raw.rds")
if (!file.exists(se_path)) stop("thyr_se_raw.rds not found (run 120 first)")
se <- readRDS(se_path)
message(
  "SE: ", nrow(se), " genes x ", ncol(se), " samples ; assay: ",
  assayNames(se)
)

clinical_path <- file.path(paths$processed, "thyr_clinical.rds")
if (!file.exists(clinical_path)) stop("thyr_clinical.rds not found (run 030 first)")
clinical <- readRDS(clinical_path)
message("Clinical: ", nrow(clinical), " cases")

as_path <- file.path(paths$processed, "thyr_case_assigned_share.rds")
if (!file.exists(as_path)) stop("thyr_case_assigned_share.rds not found (run 130 first)")
assigned_share <- as.data.frame(readRDS(as_path))
message("Assigned Share: ", nrow(assigned_share), " cases")

# --- Driver classification (030) -------------------------------------------
# RET: exact match on the three RET fusion values of Designated_Driver.
# BRAF: BRAF.MutV600E with no co-mutation (CandidateDriverMutation BRAF-only or
# empty/NA in both WGS and RNA). The two masks are exclusive for this dataset.
ret_values <- c("CCDC6-RET", "NCOA4-RET", "RET-OTHER")
ret_mask <- clinical$Designated_Driver %in% ret_values

mut_ok <- function(x) is.na(x) | x == "" | x == "BRAF"
braf_mask <- !is.na(clinical$Designated_Driver) &
  clinical$Designated_Driver == "BRAF.MutV600E" &
  mut_ok(clinical$WGS_CandidateDriverMutation) &
  mut_ok(clinical$RNA_CandidateDriverMutation)

if (any(braf_mask & ret_mask)) {
  stop("A case matched both BRAF and RET driver masks (unexpected).")
}

driver <- rep(NA_character_, nrow(clinical))
driver[braf_mask] <- "BRAF"
driver[ret_mask] <- "RET"
message("Driver: BRAF ", sum(braf_mask), " ; RET ", sum(ret_mask))

driver_tbl <- data.frame(
  case_submitter_id = as.character(clinical$REBC_ID),
  driver = driver,
  stringsAsFactors = FALSE
)
driver_tbl <- driver_tbl[!is.na(driver_tbl$driver), , drop = FALSE]

# --- Exposure classification (130) -----------------------------------------
# dose_mgy == 0 -> Sporadic ; dose_mgy > 0 with AS >= threshold -> High.
# Other cases (dose_mgy NA, or exposed with AS < threshold) drop out of scope.
as_tbl <- data.frame(
  case_submitter_id = as.character(assigned_share$REBC_ID),
  dose_mgy = as.numeric(assigned_share$dose_mgy),
  assigned_share = as.numeric(assigned_share$assigned_share_approx),
  stringsAsFactors = FALSE
)

exposure <- rep(NA_character_, nrow(as_tbl))
is_sporadic <- is.finite(as_tbl$dose_mgy) & as_tbl$dose_mgy == 0
is_high <- is.finite(as_tbl$dose_mgy) & as_tbl$dose_mgy > 0 &
  is.finite(as_tbl$assigned_share) & as_tbl$assigned_share >= as_high_threshold
exposure[is_sporadic] <- "Sporadic"
exposure[is_high] <- "High"
as_tbl$exposure <- exposure

# --- Pair resolution (120 colData) -----------------------------------------
# One _merged Primary Tumor and one _merged Solid Tissue Normal per case.
cd <- as.data.frame(colData(se))
is_merged <- grepl("_merged", cd$sample_submitter_id)
m <- cd[is_merged, , drop = FALSE]
t_rows <- m[m$sample_type == "Primary Tumor", , drop = FALSE]
n_rows <- m[m$sample_type == "Solid Tissue Normal", , drop = FALSE]

dup_t <- names(which(table(t_rows$case_submitter_id) > 1))
dup_n <- names(which(table(n_rows$case_submitter_id) > 1))
if (length(dup_t) || length(dup_n)) {
  stop(
    "A case has >1 _merged sample of the same type (pairing not unique): ",
    paste(c(dup_t, dup_n), collapse = ", ")
  )
}

tumor_of <- setNames(t_rows$sample_submitter_id, t_rows$case_submitter_id)
normal_of <- setNames(n_rows$sample_submitter_id, n_rows$case_submitter_id)
paired_cases <- intersect(names(tumor_of), names(normal_of))
message("Paired cases (_merged tumor & normal): ", length(paired_cases))

pair_tbl <- data.frame(
  case_submitter_id = paired_cases,
  tumor_id = unname(tumor_of[paired_cases]),
  normal_id = unname(normal_of[paired_cases]),
  stringsAsFactors = FALSE
)

# --- Assemble target cases -------------------------------------------------
# Driver-classified AND exposure-banded AND paired.
targets <- merge(driver_tbl, as_tbl[, c("case_submitter_id", "exposure")],
  by = "case_submitter_id"
)
targets <- targets[!is.na(targets$exposure), , drop = FALSE]
targets <- merge(targets, pair_tbl, by = "case_submitter_id")

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
