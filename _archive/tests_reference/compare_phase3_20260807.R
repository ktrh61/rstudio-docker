# compare_phase3_20260807.R
# Phase-3 regression check (dev QA): the refactored pipeline must reproduce the
# case-ID invariants captured in tests/reference/baseline_20260807/ exactly.
# Statistical values (norm factors etc.) are compared at numerical tolerance;
# any mismatch stops with a message. Run from the repository root AFTER
# re-running 140 -> 210 -> 220 -> 230 -> 310 -> 530 on the refactored scripts.

source("setup.R")

base_dir <- file.path(paths$root, "tests", "reference", "baseline_20260807")
TOL <- 1e-9

fail <- function(...) stop("PHASE3 MISMATCH: ", ..., call. = FALSE)

# --- 1. Main-BM ledger (outliers + purity) ---------------------------------
ref <- read.csv(file.path(base_dir, "cases_main_bm.csv"), stringsAsFactors = FALSE)
outliers <- readRDS(file.path(paths$processed, "thyr_case_outliers.rds"))
purity <- readRDS(file.path(paths$processed, "thyr_case_purity.rds"))
new <- merge(
  outliers,
  purity[, c("case_submitter_id", "tumor_purity")],
  by = "case_submitter_id", all.x = TRUE, sort = TRUE
)
if (!identical(nrow(new), nrow(ref))) fail("ledger row count ", nrow(new), " vs ", nrow(ref))
for (col in c("case_submitter_id", "group", "driver", "exposure", "tumor_id",
              "normal_id", "has_outlier_tumor", "has_outlier_normal")) {
  if (!identical(as.character(new[[col]]), as.character(ref[[col]]))) {
    fail("ledger column ", col)
  }
}
dp <- abs(new$tumor_purity - ref$tumor_purity)
if (any(stats::na.omit(dp) > TOL) || !identical(is.na(new$tumor_purity), is.na(ref$tumor_purity))) {
  fail("tumor_purity (max diff ", max(dp, na.rm = TRUE), ")")
}
message("1. cases_main_bm: OK (", nrow(new), " rows)")

# --- 2. Units (sample membership + factors) --------------------------------
ref_units <- read.csv(file.path(base_dir, "unit_samples.csv"), stringsAsFactors = FALSE)
normalized <- readRDS(file.path(paths$processed, "thyr_normalized_counts.rds"))
new_units <- do.call(rbind, lapply(names(normalized$units), function(unit) {
  s <- normalized$units[[unit]]$dgelist$samples
  data.frame(
    unit = unit, sample_id = rownames(s), group = as.character(s$group),
    lib.size = s$lib.size, norm.factors = s$norm.factors,
    scaling_coeff = s$scaling_coeff, stringsAsFactors = FALSE
  )
}))
key <- function(d) paste(d$unit, d$sample_id)
new_units <- new_units[order(key(new_units)), ]
ref_units <- ref_units[order(key(ref_units)), ]
if (!identical(key(new_units), key(ref_units))) fail("unit sample membership")
if (!identical(new_units$group, ref_units$group)) fail("unit group labels")
if (!identical(as.numeric(new_units$lib.size), as.numeric(ref_units$lib.size))) {
  fail("lib.size")
}
for (col in c("norm.factors", "scaling_coeff")) {
  d <- max(abs(new_units[[col]] - ref_units[[col]]))
  if (d > TOL) fail(col, " (max diff ", d, ")")
  message("2. ", col, ": OK (max diff ", format(d, digits = 3), ")")
}
message("2. unit_samples: OK (", nrow(new_units), " rows)")

# --- 3. REO evaluation ------------------------------------------------------
ref_reo <- read.csv(file.path(base_dir, "reo_eval_cases.csv"), stringsAsFactors = FALSE)
reo_eval <- readRDS(file.path(paths$processed, "thyr_reo_evaluation.rds"))
new_reo <- reo_eval$samples
new_reo <- new_reo[order(new_reo$case_submitter_id), ]
ref_reo <- ref_reo[order(ref_reo$case_submitter_id), ]
if (!identical(new_reo$case_submitter_id, ref_reo$case_submitter_id)) fail("REO case set")
for (col in c("band", "tumor_id", "class")) {
  if (!identical(as.character(new_reo[[col]]), as.character(ref_reo[[col]]))) fail("REO ", col)
}
if (!identical(as.integer(new_reo$score), as.integer(ref_reo$score))) fail("REO score")
for (col in c("assigned_share", "dose_mgy")) {
  d <- max(abs(new_reo[[col]] - ref_reo[[col]]))
  if (d > TOL) fail("REO ", col, " (max diff ", d, ")")
}
message("3. reo_eval_cases: OK (", nrow(new_reo), " rows)")

# --- 4. REO read-B / config -------------------------------------------------
ref_sum <- read.csv(file.path(base_dir, "reo_summary.csv"), stringsAsFactors = FALSE)
ref_kv <- setNames(ref_sum$value, ref_sum$key)
stopifnot(abs(reo_eval$read_B$p_value - as.numeric(ref_kv[["read_B.p_value"]])) < TOL)
stopifnot(identical(reo_eval$read_B$method, ref_kv[["read_B.method"]]))
stopifnot(abs(reo_eval$read_B$effect_P_low_lt_mid -
  as.numeric(ref_kv[["read_B.effect_P_low_lt_mid"]])) < TOL)
stopifnot(reo_eval$read_A_threshold == as.numeric(ref_kv[["read_A_threshold"]]))
message("4. reo read_A/read_B: OK (p = ", reo_eval$read_B$p_value, ")")

message("\nPHASE 3 REGRESSION: ALL CHECKS PASSED")
