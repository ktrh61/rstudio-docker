# 130_compute_assigned_share.R
# Attach the per-case Assigned Share to the prepared clinical table. The values
# are IREP results computed by the researcher (NIH IREP thyroid model,
# cited as version 5.7.3,
# "Assigned Share associated with the expected value of ERR"; radiation type
# electrons E>15keV; sex as recorded; dose entered in cSv = mGy/10; birth year = 1986 -
# AGE_EXPOSURE, exposure year 1986, surgery year = birth year + AGE_SURGERY;
# all IREP assumptions and settings at their defaults: user-defined
# uncertainty distribution Lognormal(1,1), 10,000 iterations, random number
# seed 99),
# read from the canonical file raw/clinical/thyr_irep_assigned_share.csv
# (2 columns, exposed driver-classified cases only; provenance, the input-
# identity verification against this clinical table, and the REBC-AC8W
# adjudication are recorded in reorg plan v2 B.11).
#
# The earlier local Monte Carlo approximation (lib/dose_thyroid_as.R) is
# retired from the pipeline and archived with annotation in _archive/lib/
# (B.11; it tracks the IREP values at median |diff| ~0.1 AS points).
#
# Input : processed/thyr_clinical.rds               (from 030; key REBC_ID)
#         raw/clinical/thyr_irep_assigned_share.csv (researcher-computed IREP)
# Output: processed/thyr_case_assigned_share.rds  (REBC_ID-keyed, 440 rows)
#
# Output columns:
#   REBC_ID                case key (S1 form, no YQ)
#   assigned_share         IREP AS (percent, 0-100), or NA
#   assigned_share_status  "irep"                  value present
#                          "not_required_sporadic" dose == 0; AS undefined
#                          "not_in_reference"      dose > 0 but not in the
#                                                  reference (driver-
#                                                  unclassified cases; they
#                                                  enter no analysis cohort)
#   dose_mgy, age_exposure, age_surgery   inputs copied from the clinical table
#
# 140 stops if any driver-classified exposed case lacks a value (the guard
# lives there because driver classification does not exist yet at this step).

source("setup.R")

library(data.table)

# --- Load inputs -----------------------------------------------------------
irep_path <- file.path(paths$root, "raw", "clinical", "thyr_irep_assigned_share.csv")
if (!file.exists(irep_path)) {
  stop("IREP Assigned Share file not found: ", irep_path)
}
irep <- fread(irep_path, colClasses = list(character = "case_submitter_id"))
if (!identical(names(irep), c("case_submitter_id", "assigned_share"))) {
  stop("Unexpected columns in IREP file: ", paste(names(irep), collapse = ", "))
}
if (anyDuplicated(irep$case_submitter_id)) {
  stop("Duplicate case_submitter_id in IREP file.")
}
if (any(!is.finite(irep$assigned_share) |
        irep$assigned_share <= 0 | irep$assigned_share >= 100)) {
  stop("IREP assigned_share values must be finite and in (0, 100).")
}
message("IREP Assigned Share: ", nrow(irep), " cases")

clinical_path <- file.path(paths$processed, "thyr_clinical.rds")
if (!file.exists(clinical_path)) {
  stop("Clinical table not found: ", clinical_path, " (run 030 first)")
}
clinical <- readRDS(clinical_path)
message("Clinical loaded: ", nrow(clinical), " cases")

needed <- c("REBC_ID", "DOSE", "AGE_EXPOSURE", "AGE_SURGERY")
missing_cols <- setdiff(needed, names(clinical))
if (length(missing_cols) > 0) {
  stop("Clinical missing expected columns: ", paste(missing_cols, collapse = ", "))
}

res <- data.table(
  REBC_ID      = as.character(clinical$REBC_ID),
  dose_mgy     = as.numeric(clinical$DOSE),
  age_exposure = as.numeric(clinical$AGE_EXPOSURE),
  age_surgery  = as.numeric(clinical$AGE_SURGERY)
)

# --- Join and status -------------------------------------------------------
unknown <- setdiff(irep$case_submitter_id, res$REBC_ID)
if (length(unknown) > 0) {
  stop("IREP file has cases absent from the clinical table: ",
       paste(unknown, collapse = ", "))
}

res[, assigned_share := irep$assigned_share[match(REBC_ID, irep$case_submitter_id)]]

if (any(res$dose_mgy == 0 & !is.na(res$assigned_share))) {
  stop("IREP file carries a value for a dose == 0 case (AS is undefined there).")
}

res[, assigned_share_status := fifelse(
  !is.na(assigned_share), "irep",
  fifelse(dose_mgy == 0, "not_required_sporadic", "not_in_reference")
)]

setcolorder(res, c(
  "REBC_ID",
  "assigned_share",
  "assigned_share_status",
  "dose_mgy", "age_exposure", "age_surgery"
))
setkey(res, REBC_ID)

message("Assigned Share attached: ", nrow(res), " cases")
print(table(res$assigned_share_status))

# --- Save ------------------------------------------------------------------
saveRDS(res, file.path(paths$processed, "thyr_case_assigned_share.rds"))
message(
  "Saved: processed/thyr_case_assigned_share.rds (",
  nrow(res), " cases x ", ncol(res), " columns)"
)
