# Case-design primitives: driver classification, exposure-band assignment, and
# _merged tumor/normal pair resolution. Single source of truth consumed by
# scripts/140_build_case_design.R; downstream scripts filter the resulting
# design table instead of re-deriving any of this.
#
# Band rule (single convention; boundary-valued cases do not occur in this
# cohort, so unifying the historical closure variants changes no assignment):
#   dose_mgy == 0                                -> Sporadic (AS not required)
#   dose_mgy > 0 & 0 < AS <= AS_LOW_MAX          -> Low
#   dose_mgy > 0 & AS_LOW_MAX < AS < AS_HIGH_MIN -> Mid
#   dose_mgy > 0 & AS >= AS_HIGH_MIN             -> High
# AS_LOW_MAX / AS_HIGH_MIN come from config.R.

RET_DRIVER_VALUES <- c("CCDC6-RET", "NCOA4-RET", "RET-OTHER")
TUMOR_SAMPLE_TYPE <- "Primary Tumor"
NORMAL_SAMPLE_TYPE <- "Solid Tissue Normal"

# Driver classification. RET: exact match on the three RET fusion values of
# Designated_Driver. BRAF: BRAF.MutV600E with no co-mutation
# (CandidateDriverMutation BRAF-only or empty/NA in both WGS and RNA).
# Everything else keeps driver = NA with the reason in driver_status.
classify_driver <- function(clinical) {
  ret_mask <- clinical$Designated_Driver %in% RET_DRIVER_VALUES

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

  driver_status <- rep("other_or_unclassified", nrow(clinical))
  driver_status[ret_mask] <- "ret_fusion"
  driver_status[braf_mask] <- "braf_v600e_no_comutation"

  data.frame(
    case_submitter_id = as.character(clinical$REBC_ID),
    designated_driver_raw = as.character(clinical$Designated_Driver),
    driver = driver,
    driver_status = driver_status,
    stringsAsFactors = FALSE
  )
}

# Exposure-band assignment. Judgment order: dose validity first, dose == 0 is
# Sporadic regardless of AS, only then is AS consulted.
assign_band <- function(dose_mgy, assigned_share) {
  n <- length(dose_mgy)
  band <- rep(NA_character_, n)
  band_status <- rep(NA_character_, n)

  dose_bad <- !is.finite(dose_mgy) | dose_mgy < 0
  band_status[dose_bad] <- "dose_missing_or_invalid"

  sporadic <- !dose_bad & dose_mgy == 0
  band[sporadic] <- "Sporadic"
  band_status[sporadic] <- "ok"

  exposed <- !dose_bad & dose_mgy > 0
  as_bad <- exposed & !is.finite(assigned_share)
  band_status[as_bad] <- "as_missing"
  as_ok <- exposed & is.finite(assigned_share)

  low <- as_ok & assigned_share > 0 & assigned_share <= AS_LOW_MAX
  mid <- as_ok & assigned_share > AS_LOW_MAX & assigned_share < AS_HIGH_MIN
  high <- as_ok & assigned_share >= AS_HIGH_MIN
  band[low] <- "Low"
  band[mid] <- "Mid"
  band[high] <- "High"
  band_status[low | mid | high] <- "ok"
  band_status[as_ok & !(low | mid | high)] <- "as_out_of_range"

  data.frame(band = band, band_status = band_status, stringsAsFactors = FALSE)
}

# _merged tumor/normal pair resolution from SE colData. One _merged sample per
# case and sample type (stop otherwise); a missing tumor or normal is NA, and
# the case row is kept (has_tumor / has_normal / is_paired record the state).
resolve_merged_pairs <- function(coldata) {
  cd <- as.data.frame(coldata)
  m <- cd[grepl("_merged", cd$sample_submitter_id), , drop = FALSE]
  t_rows <- m[m$sample_type == TUMOR_SAMPLE_TYPE, , drop = FALSE]
  n_rows <- m[m$sample_type == NORMAL_SAMPLE_TYPE, , drop = FALSE]

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
  cases <- sort(union(names(tumor_of), names(normal_of)))

  data.frame(
    case_submitter_id = cases,
    tumor_id = unname(tumor_of[cases]),
    normal_id = unname(normal_of[cases]),
    stringsAsFactors = FALSE
  )
}
