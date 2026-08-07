# 230_finalize_analysis_cohorts.R
# Finalize per-analysis cohort membership from the case-design table plus the
# QC results: every adoption condition becomes an explicit boolean, and each
# analysis gets one include_* flag. Downstream scripts read these flags and
# never re-derive eligibility. The case-flow table for the paper is generated
# from this object, not recounted per script.
# Input : processed/thyr_case_design.rds   (from 140)
#         processed/thyr_case_outliers.rds (from 210; PC-OD flags)
#         processed/thyr_case_purity.rds   (from 220; pooled relative purity)
# Output: processed/thyr_analysis_cohorts.rds
#           design columns + pcod_tumor_ok, pcod_normal_ok, tumor_purity,
#           relative_purity_ok, include_main_bm, include_reo_training,
#           include_reo_evaluation, first_exclusion_reason_main_bm
#         processed/thyr_cohort_flow.rds
#           stepwise case-flow counts for the main BM cohort
#
# Adoption conditions:
#   main BM        : driver RET|BRAF, band Sporadic|High, is_paired,
#                    PC-OD clean (both tissues), relative purity >= PURITY_THRESHOLD
#   REO training   : include_main_bm & driver RET (same filters as main BM)
#   REO evaluation : driver RET, band Low|Mid, has_tumor (no pairing, PC-OD,
#                    or purity condition -- deliberate, see 530 header)

source("setup.R")

# --- Load inputs -----------------------------------------------------------
design_path <- file.path(paths$processed, "thyr_case_design.rds")
if (!file.exists(design_path)) stop("thyr_case_design.rds not found (run 140 first)")
design <- readRDS(design_path)

outliers_path <- file.path(paths$processed, "thyr_case_outliers.rds")
if (!file.exists(outliers_path)) stop("thyr_case_outliers.rds not found (run 210 first)")
outliers <- readRDS(outliers_path)

purity_path <- file.path(paths$processed, "thyr_case_purity.rds")
if (!file.exists(purity_path)) stop("thyr_case_purity.rds not found (run 220 first)")
purity <- readRDS(purity_path)

# --- Eligibility booleans --------------------------------------------------
cohorts <- design
cohorts$eligible_driver_main <- !is.na(cohorts$driver)
cohorts$eligible_band_main <- cohorts$band %in% c("Sporadic", "High")

flags <- outliers[, c("case_submitter_id", "has_outlier_tumor", "has_outlier_normal")]
cohorts <- merge(cohorts, flags, by = "case_submitter_id", all.x = TRUE, sort = TRUE)
# NA = PC-OD did not evaluate this case (out of the 210 target set).
cohorts$pcod_tumor_ok <- !is.na(cohorts$has_outlier_tumor) & cohorts$has_outlier_tumor == 0
cohorts$pcod_normal_ok <- !is.na(cohorts$has_outlier_normal) & cohorts$has_outlier_normal == 0

pur <- purity[, c("case_submitter_id", "tumor_purity")]
cohorts <- merge(cohorts, pur, by = "case_submitter_id", all.x = TRUE, sort = TRUE)
cohorts$relative_purity_ok <- !is.na(cohorts$tumor_purity) &
  cohorts$tumor_purity >= PURITY_THRESHOLD

cohorts$has_tumor_reo <- cohorts$has_tumor
cohorts$eligible_band_reo <- cohorts$band %in% c("Low", "Mid")

# --- Analysis membership ---------------------------------------------------
cohorts$include_main_bm <- cohorts$eligible_driver_main &
  cohorts$eligible_band_main &
  cohorts$is_paired &
  cohorts$pcod_tumor_ok &
  cohorts$pcod_normal_ok &
  cohorts$relative_purity_ok

cohorts$include_reo_training <- cohorts$include_main_bm & cohorts$driver %in% "RET"

cohorts$include_reo_evaluation <- cohorts$driver %in% "RET" &
  cohorts$eligible_band_reo &
  cohorts$has_tumor_reo

# --- First exclusion reason (fixed order, main BM) -------------------------
reason <- rep(NA_character_, nrow(cohorts))
reason[!cohorts$relative_purity_ok] <- "low_relative_purity"
reason[!(cohorts$pcod_tumor_ok & cohorts$pcod_normal_ok)] <- "pcod_outlier"
reason[!cohorts$is_paired] <- "unpaired"
reason[!cohorts$eligible_band_main] <- "band_not_sporadic_or_high"
reason[!cohorts$eligible_driver_main] <- "driver_unclassified"
reason[cohorts$include_main_bm] <- NA_character_
cohorts$first_exclusion_reason_main_bm <- reason

rownames(cohorts) <- NULL

# --- Case-flow table (main BM) ---------------------------------------------
step_counts <- function(mask) {
  data.frame(
    n_total = sum(mask),
    n_RET = sum(mask & cohorts$driver %in% "RET"),
    n_BRAF = sum(mask & cohorts$driver %in% "BRAF")
  )
}
s0 <- rep(TRUE, nrow(cohorts))
s1 <- cohorts$eligible_driver_main
s2 <- s1 & cohorts$eligible_band_main
s3 <- s2 & cohorts$is_paired
s4 <- s3 & cohorts$pcod_tumor_ok & cohorts$pcod_normal_ok
s5 <- s4 & cohorts$relative_purity_ok
flow <- cbind(
  step = c(
    "all_cases", "driver_classified", "band_sporadic_or_high",
    "paired", "pcod_clean", "purity_pass"
  ),
  do.call(rbind, lapply(list(s0, s1, s2, s3, s4, s5), step_counts))
)

message("Cohort flow (main BM):")
print(flow)
message("include_main_bm by group:")
print(table(
  paste0(ifelse(cohorts$driver == "RET", "R", "B"), "_", cohorts$band)[cohorts$include_main_bm]
))
message(
  "REO training: ", sum(cohorts$include_reo_training),
  " ; REO evaluation: ", sum(cohorts$include_reo_evaluation),
  " (Low ", sum(cohorts$include_reo_evaluation & cohorts$band == "Low"),
  " / Mid ", sum(cohorts$include_reo_evaluation & cohorts$band == "Mid"), ")"
)

out <- file.path(paths$processed, "thyr_analysis_cohorts.rds")
saveRDS(cohorts, out)
message("Saved: ", out, " (", nrow(cohorts), " x ", ncol(cohorts), ")")

flow_out <- file.path(paths$processed, "thyr_cohort_flow.rds")
saveRDS(flow, flow_out)
message("Saved: ", flow_out)
