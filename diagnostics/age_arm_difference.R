# age_arm_difference.R
# Between-arm age-at-surgery difference, estimated with CI (claim map C-15;
# reading fixed there before this script's first run).
#
# Purpose: satisfy the ratified presentation style for age structure (plan
# v2 section 0.5, round 2: description + estimation with CI) by disclosing
# the magnitude and uncertainty of the High-vs-Sporadic age difference in
# each driver stratum of the main BM cohort. This is a disclosure, not a
# confounding test: no p-value is computed or printed, and no claim of the
# paper moves with these estimates (C-15). The primary answer to the age
# objection stays on the literature side (objections ledger Q-13 (i)).
#
# Estimators, per stratum (R: 12 Sporadic vs 15 High; B: 27 vs 9):
#   - Hodges-Lehmann shift: median over all pairwise differences
#     High_i - Sporadic_j, in years.
#   - Rank effect P(Sporadic < High) + 0.5 P(=): the same Brunner-Munzel
#     effect estimator as 410 (.bm_effect from lib/stat_brunnermunzel.R;
#     only the estimator is called -- no test is invoked, so nothing here
#     produces a p-value).
#   - 95% percentile bootstrap CI for both, resampling within arm with
#     replacement, B = 9999 replicates, seed 19450809 (diagnostic seed base).
#
# Age at exposure is structurally not comparable between arms (the Sporadic
# arm is unexposed by definition, AGE_EXPOSURE all NA); the per-arm NA
# counts are printed to document this and nothing is estimated.
#
# Consistency check: per-arm n and median [range] printed below must match
# N-09 / N-12 / N-13 of the numbers ledger.
#
# Input : processed/thyr_analysis_cohorts.rds (include_main_bm, driver, band)
#         processed/thyr_clinical.rds         (REBC_ID, AGE_SURGERY, AGE_EXPOSURE)
# Output: diagnostics/output/age_arm_difference.rds

source("setup.R")
source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))

SEED <- 19450809L
B_BOOT <- 9999L

cohorts <- readRDS(file.path(paths$processed, "thyr_analysis_cohorts.rds"))
clinical <- readRDS(file.path(paths$processed, "thyr_clinical.rds"))

m <- merge(
  cohorts[cohorts$include_main_bm, c("case_submitter_id", "driver", "band")],
  clinical[, c("REBC_ID", "AGE_SURGERY", "AGE_EXPOSURE")],
  by.x = "case_submitter_id", by.y = "REBC_ID", all.x = TRUE
)
stopifnot(
  nrow(m) == 63,
  !any(is.na(m$AGE_SURGERY)),
  m$band %in% c("Sporadic", "High")
)

hl_shift <- function(x, y) median(outer(y, x, "-")) # y = High, x = Sporadic
fmt_med <- function(v) {
  sprintf("%g [%g-%g]", median(v), min(v), max(v))
}

cat(
  "Age-at-surgery arm difference (High - Sporadic), main BM cohort;",
  "HL shift + BM effect, percentile bootstrap B =", B_BOOT,
  ", seed", SEED, "\n"
)

set.seed(SEED)
result_rows <- list()
boot_draws <- list()
for (drv in c("RET", "BRAF")) {
  x <- m$AGE_SURGERY[m$driver == drv & m$band == "Sporadic"]
  y <- m$AGE_SURGERY[m$driver == drv & m$band == "High"]

  obs_hl <- hl_shift(x, y)
  obs_eff <- .bm_effect(x, y)
  boot_hl <- numeric(B_BOOT)
  boot_eff <- numeric(B_BOOT)
  for (b in seq_len(B_BOOT)) {
    xb <- sample(x, replace = TRUE)
    yb <- sample(y, replace = TRUE)
    boot_hl[b] <- hl_shift(xb, yb)
    boot_eff[b] <- .bm_effect(xb, yb)
  }
  ci_hl <- unname(quantile(boot_hl, c(0.025, 0.975)))
  ci_eff <- unname(quantile(boot_eff, c(0.025, 0.975)))

  stratum <- if (drv == "RET") "R" else "B"
  cat(sprintf(
    "\n== %s (%s) ==  Sporadic n=%d %s | High n=%d %s\n",
    stratum, drv, length(x), fmt_med(x), length(y), fmt_med(y)
  ))
  cat(sprintf(
    "  HL shift (years)      %+.1f  [%.1f, %.1f]\n",
    obs_hl, ci_hl[1], ci_hl[2]
  ))
  cat(sprintf(
    "  P(Sporadic<High) (BM)  %.3f  [%.3f, %.3f]\n",
    obs_eff, ci_eff[1], ci_eff[2]
  ))

  result_rows[[stratum]] <- data.frame(
    stratum = stratum, driver = drv,
    n_sporadic = length(x), n_high = length(y),
    median_sporadic = median(x), median_high = median(y),
    hl_shift = obs_hl, hl_ci_lo = ci_hl[1], hl_ci_hi = ci_hl[2],
    bm_effect = obs_eff, eff_ci_lo = ci_eff[1], eff_ci_hi = ci_eff[2]
  )
  boot_draws[[stratum]] <- data.frame(hl = boot_hl, effect = boot_eff)
}

cat("\nAge at exposure (not estimable between arms -- Sporadic unexposed):\n")
for (drv in c("RET", "BRAF")) {
  for (bnd in c("Sporadic", "High")) {
    v <- m$AGE_EXPOSURE[m$driver == drv & m$band == bnd]
    cat(sprintf(
      "  %s_%-9s n=%2d  NA %2d / non-NA %2d\n",
      if (drv == "RET") "R" else "B", bnd,
      length(v), sum(is.na(v)), sum(!is.na(v))
    ))
  }
}

result <- list(
  date = Sys.time(),
  config = list(seed = SEED, b_boot = B_BOOT, ci = "percentile 2.5/97.5"),
  summary = do.call(rbind, result_rows),
  bootstrap = boot_draws
)
out_dir <- file.path(paths$root, "diagnostics", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(result, file.path(out_dir, "age_arm_difference.rds"))
cat("\nSaved:", file.path(out_dir, "age_arm_difference.rds"), "\n")
