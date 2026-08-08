# 140_build_case_design.R
# Build the single case-design table: every clinical case with its driver
# classification, exposure band, and _merged tumor/normal sample resolution.
# No case is dropped here; ineligibility is recorded in status columns so that
# downstream scripts (230 and the diagnostics) filter one shared table instead
# of re-deriving drivers, bands, or pairs.
# Input : processed/thyr_clinical.rds             (from 030)
#         processed/thyr_case_assigned_share.rds  (from 130)
#         processed/thyr_se_raw.rds               (from 120; colData only)
#         lib/cohort_design.R
# Output: processed/thyr_case_design.rds
#           columns: case_submitter_id, designated_driver_raw, driver,
#             driver_status, dose_mgy, assigned_share,
#             assigned_share_status, band, band_status,
#             tumor_id, normal_id, has_tumor, has_normal, is_paired
#
# Guard: every driver-classified (RET/BRAF) exposed case must carry an IREP
# assigned_share -- a silent gap here would silently shrink the analysis
# cohorts, so it stops. (The guard lives here, not in 130, because driver
# classification first exists in this script.)

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

source(file.path(paths$root, "lib", "cohort_design.R"))

# --- Load inputs -----------------------------------------------------------
clinical_path <- file.path(paths$processed, "thyr_clinical.rds")
if (!file.exists(clinical_path)) stop("thyr_clinical.rds not found (run 030 first)")
clinical <- readRDS(clinical_path)
message("Clinical: ", nrow(clinical), " cases")

as_path <- file.path(paths$processed, "thyr_case_assigned_share.rds")
if (!file.exists(as_path)) stop("thyr_case_assigned_share.rds not found (run 130 first)")
assigned_share <- as.data.frame(readRDS(as_path))
message("Assigned Share: ", nrow(assigned_share), " cases")

se_path <- file.path(paths$processed, "thyr_se_raw.rds")
if (!file.exists(se_path)) stop("thyr_se_raw.rds not found (run 120 first)")
se <- readRDS(se_path)

# --- Key uniqueness --------------------------------------------------------
if (anyDuplicated(clinical$REBC_ID)) stop("Duplicate REBC_ID in clinical table.")
if (anyDuplicated(assigned_share$REBC_ID)) stop("Duplicate REBC_ID in assigned-share table.")

# --- Assemble the design table ---------------------------------------------
design <- classify_driver(clinical)

as_tbl <- data.frame(
  case_submitter_id = as.character(assigned_share$REBC_ID),
  dose_mgy = as.numeric(assigned_share$dose_mgy),
  assigned_share = as.numeric(assigned_share$assigned_share),
  assigned_share_status = as.character(assigned_share$assigned_share_status),
  stringsAsFactors = FALSE
)
design <- merge(design, as_tbl, by = "case_submitter_id", all = TRUE, sort = TRUE)

classified_gap <- design$driver %in% c("RET", "BRAF") &
  design$dose_mgy > 0 & !is.finite(design$assigned_share)
if (any(classified_gap)) {
  stop(
    "Driver-classified exposed case(s) lack an IREP assigned_share: ",
    paste(design$case_submitter_id[classified_gap], collapse = ", ")
  )
}

band_tbl <- assign_band(design$dose_mgy, design$assigned_share)
design$band <- band_tbl$band
design$band_status <- band_tbl$band_status

pairs <- resolve_merged_pairs(colData(se))
design <- merge(design, pairs, by = "case_submitter_id", all.x = TRUE, sort = TRUE)
design$has_tumor <- !is.na(design$tumor_id)
design$has_normal <- !is.na(design$normal_id)
design$is_paired <- design$has_tumor & design$has_normal

rownames(design) <- NULL

message("Design table: ", nrow(design), " cases")
message("Driver:")
print(table(design$driver, useNA = "ifany"))
message("Band:")
print(table(design$band, useNA = "ifany"))
message("Driver x band (paired only):")
print(table(design$driver[design$is_paired], design$band[design$is_paired]))

out <- file.path(paths$processed, "thyr_case_design.rds")
saveRDS(design, out)
message("Saved: ", out, " (", nrow(design), " x ", ncol(design), ")")
