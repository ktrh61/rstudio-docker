# 041_compute_assigned_share.R
# Attach a per-case approximate Assigned Share to the clinical master. The
# quantity is the "Assigned Share associated with the expected value of ERR",
# computed by compute_thyroid_as() as a Monte Carlo approximation of IREP; it is
# a hand estimate, not an IREP-produced value. Attached to all 440 cases; no
# target selection or grouping here (deferred to 042).
#
# Input : processed/thyr_clinical_master.rds   (from 040; master key case_id)
#         utils/thyroid_as_functions.R         (compute_thyroid_as)
# Output: processed/thyr_case_assigned_share.rds  (case_id-keyed, 440 rows)
#         output/thyr_case_assigned_share.csv     (human-readable)
#
# Output columns:
#   case_id                       master key (S1 form, no YQ)
#   assigned_share_approx         approximate AS (percent, 0-100), or NA
#   assigned_share_approx_status  "computed", or a reason string when NA
#   dose_mgy, age_exposure,       inputs copied from the master
#     age_surgery
#
# A case is computed only when DOSE, AGE_EXPOSURE and AGE_SURGERY are all finite
# and DOSE is non-negative; otherwise assigned_share_approx is NA with the
# reason recorded. n_iter and seed use the compute_thyroid_as() defaults.

source("setup.R")

suppressPackageStartupMessages({
  library(data.table)
})

# --- Load estimator and master ---------------------------------------------
as_fun_path <- file.path(paths$root, "utils", "thyroid_as_functions.R")
if (!file.exists(as_fun_path)) {
  stop("Assigned Share functions not found: ", as_fun_path)
}
source(as_fun_path)

master_path <- file.path(paths$processed, "thyr_clinical_master.rds")
if (!file.exists(master_path)) {
  stop("Clinical master not found: ", master_path, " (run 040 first)")
}
master <- readRDS(master_path)
message("Master loaded: ", nrow(master), " cases")

needed <- c("case_id", "DOSE", "AGE_EXPOSURE", "AGE_SURGERY")
missing_cols <- setdiff(needed, names(master))
if (length(missing_cols) > 0) {
  stop("Master missing expected columns: ", paste(missing_cols, collapse = ", "))
}

res <- data.table(
  case_id      = as.character(master$case_id),
  dose_mgy     = as.numeric(master$DOSE),
  age_exposure = as.numeric(master$AGE_EXPOSURE),
  age_surgery  = as.numeric(master$AGE_SURGERY)
)

# --- Computability check ---------------------------------------------------
finite_all <- is.finite(res$dose_mgy) &
  is.finite(res$age_exposure) &
  is.finite(res$age_surgery)
dose_ok <- is.finite(res$dose_mgy) & res$dose_mgy >= 0

status <- rep("computed", nrow(res))
status[!finite_all] <- "missing_input"
status[finite_all & !dose_ok] <- "negative_dose"
computable <- status == "computed"

# --- Compute ---------------------------------------------------------------
# compute_thyroid_as() restores the caller's RNG state on exit, so seeds do not
# accumulate across cases.
approx <- rep(NA_real_, nrow(res))
for (i in which(computable)) {
  approx[i] <- compute_thyroid_as(
    dose_mGy     = res$dose_mgy[i],
    age_exposure = res$age_exposure[i],
    age_surgery  = res$age_surgery[i]
  )
}

res[, assigned_share_approx := approx]
res[, assigned_share_approx_status := status]

setcolorder(res, c(
  "case_id",
  "assigned_share_approx",
  "assigned_share_approx_status",
  "dose_mgy", "age_exposure", "age_surgery"
))
setkey(res, case_id)

message(
  "Assigned Share (approx) computed: ", sum(status == "computed"),
  " / ", nrow(res)
)
message("  NA (missing_input): ", sum(status == "missing_input"))
message("  NA (negative_dose): ", sum(status == "negative_dose"))

# --- Save ------------------------------------------------------------------
saveRDS(res, file.path(paths$processed, "thyr_case_assigned_share.rds"))
message(
  "Saved: processed/thyr_case_assigned_share.rds (",
  nrow(res), " cases x ", ncol(res), " columns)"
)

fwrite(res, file.path(paths$output, "thyr_case_assigned_share.csv"), na = "NA")
message("Saved readable: output/thyr_case_assigned_share.csv")
