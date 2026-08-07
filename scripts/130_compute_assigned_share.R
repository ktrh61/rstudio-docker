# 130_compute_assigned_share.R
# Attach a per-case approximate Assigned Share to the prepared clinical table.
# The quantity is the "Assigned Share associated with the expected value of ERR",
# computed by compute_thyroid_as() as a Monte Carlo approximation of IREP; it is
# a hand estimate, not an IREP-produced value. Attached to all 440 cases; no
# target selection or grouping here.
#
# Input : processed/thyr_clinical.rds          (from 030; key REBC_ID)
#         lib/dose_thyroid_as.R         (compute_thyroid_as)
# Output: processed/thyr_case_assigned_share.rds  (REBC_ID-keyed, 440 rows)
#
# Output columns:
#   REBC_ID                       case key (S1 form, no YQ)
#   assigned_share_approx         approximate AS (percent, 0-100), or NA
#   assigned_share_approx_status  "computed", or a reason string when NA
#   dose_mgy, age_exposure,       inputs copied from the clinical table
#     age_surgery
#
# A case is computed only when DOSE, AGE_EXPOSURE and AGE_SURGERY are all
# finite; otherwise assigned_share_approx is NA with the reason recorded.
# n_iter and seed use the compute_thyroid_as() defaults.

source("setup.R")

library(data.table)

# --- Load estimator and clinical table -------------------------------------
as_fun_path <- file.path(paths$root, "lib", "dose_thyroid_as.R")
if (!file.exists(as_fun_path)) {
  stop("Assigned Share functions not found: ", as_fun_path)
}
source(as_fun_path)

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

# --- Computability check ---------------------------------------------------
finite_all <- is.finite(res$dose_mgy) &
  is.finite(res$age_exposure) &
  is.finite(res$age_surgery)

status <- rep("computed", nrow(res))
status[!finite_all] <- "missing_input"
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
  "REBC_ID",
  "assigned_share_approx",
  "assigned_share_approx_status",
  "dose_mgy", "age_exposure", "age_surgery"
))
setkey(res, REBC_ID)

message(
  "Assigned Share (approx) computed: ", sum(status == "computed"),
  " / ", nrow(res)
)
message("  NA (missing_input): ", sum(status == "missing_input"))

# --- Save ------------------------------------------------------------------
saveRDS(res, file.path(paths$processed, "thyr_case_assigned_share.rds"))
message(
  "Saved: processed/thyr_case_assigned_share.rds (",
  nrow(res), " cases x ", ncol(res), " columns)"
)
