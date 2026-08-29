# tab_case_characteristics.R  (Table 2)
# Case characteristics by analysis group: main BM (R_Sporadic/R_High/
# B_Sporadic/B_High) and REO evaluation (R_Low/R_Mid). Formatting only --
# no computation beyond counts/medians of frozen inputs (plus the stratum-by-
# band pool sizes, N-95); canonical sources
# for the ledger remain the primary rds (N-09, N-12, N-13). CSV is a
# typesetting convenience. Age-difference footnote values (N-64, N-65) come
# from diagnostics/output/age_arm_difference.rds and are appended as notes.
# Input : processed/thyr_analysis_cohorts.rds, processed/thyr_clinical.rds,
#         diagnostics/output/age_arm_difference.rds
# Output: output/tables/tab_case_characteristics.csv,
#         output/tables/tab_case_characteristics_notes.txt (+ printed table)

source("setup.R")

co <- readRDS(file.path(paths$processed, "thyr_analysis_cohorts.rds"))
cl <- readRDS(file.path(paths$processed, "thyr_clinical.rds"))
age_path <- file.path(paths$root, "diagnostics", "output", "age_arm_difference.rds")
if (!file.exists(age_path)) {
  stop("age_arm_difference.rds not found (run diagnostics/age_arm_difference.R first)")
}
age <- readRDS(age_path)
required_age_columns <- c(
  "stratum", "hl_shift", "hl_ci_lo", "hl_ci_hi",
  "bm_effect", "eff_ci_lo", "eff_ci_hi"
)
stopifnot(
  is.data.frame(age$summary),
  all(required_age_columns %in% names(age$summary)),
  setequal(age$summary$stratum, c("R", "B")),
  identical(age$config$b_boot, 9999L),
  identical(age$config$seed, 19450809L)
)
m <- merge(
  co, cl[, c("REBC_ID", "SEX", "AGE_SURGERY", "AGE_EXPOSURE", "Designated_Driver")],
  by.x = "case_submitter_id", by.y = "REBC_ID", all.x = TRUE
)
sel <- m[m$include_main_bm | m$include_reo_evaluation, ]
sel$group <- paste0(ifelse(sel$driver == "RET", "R", "B"), "_", sel$band)

fmt_med <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0) return("NA")
  sprintf("%g [%g-%g]", stats::median(v), min(v), max(v))
}
# Stratum-by-band pool each analysed group was drawn from: every case of the
# stratum in that band before the sample-based steps (pair, outlier screen,
# purity). Cell = pool size (cases with a tumor-normal pair). Replaces the
# former Table S1 (researcher decision 2026-08-29).
pool_of <- function(g) {
  d <- if (startsWith(g, "R_")) "RET" else "BRAF"
  b <- sub("^[RB]_", "", g)
  pool <- co[co$driver %in% d & co$band %in% b, ]
  sprintf("%d (%d)", nrow(pool), sum(pool$is_paired))
}
rows <- lapply(split(sel, sel$group), function(s) {
  drv <- sort(table(s$Designated_Driver), decreasing = TRUE)
  data.frame(
    group = s$group[1], n = nrow(s), pool = pool_of(s$group[1]),
    female = sum(s$SEX %in% c("F", "Female", "female")),
    male = sum(s$SEX %in% c("M", "Male", "male")),
    age_surgery = fmt_med(s$AGE_SURGERY),
    age_exposure = fmt_med(s$AGE_EXPOSURE),
    driver_detail = paste(sprintf("%s %d", names(drv), drv), collapse = "; ")
  )
})
ord <- c("R_Sporadic", "R_Low", "R_Mid", "R_High", "B_Sporadic", "B_High")
tab <- do.call(rbind, rows[ord[ord %in% names(rows)]])
rownames(tab) <- NULL
print(tab)

age_row <- function(stratum) {
  row <- age$summary[age$summary$stratum == stratum, , drop = FALSE]
  stopifnot(nrow(row) == 1L)
  row
}
format_age_estimates <- function(label, row) {
  sprintf(
    paste0(
      "%s: Hodges–Lehmann difference %+.1f years ",
      "(95%% percentile bootstrap CI, %.1f to %.1f) and ",
      "Brunner–Munzel relative effect θ=%.3f (%.3f to %.3f)"
    ),
    label, row$hl_shift, row$hl_ci_lo, row$hl_ci_hi,
    row$bm_effect, row$eff_ci_lo, row$eff_ci_hi
  )
}
age_note <- paste0(
  "Age-at-surgery differences compare High-AS with dose-zero cases. ",
  format_age_estimates("RET", age_row("R")), "; ",
  format_age_estimates("BRAF", age_row("B")), ". ",
  "Here, θ=Pr(X<Y)+0.5Pr(X=Y), with X denoting dose-zero and Y High-AS. ",
  "Intervals were obtained from 9,999 resamples drawn separately within each group ",
  "(seed 19450809); no p-values were calculated."
)
cat("\nTable 1 footnote:\n", age_note, "\n", sep = "")

if (any(is.na(sel$SEX))) cat("NOTE: SEX has NA values; recheck coding\n")
cat("sex coding seen:", paste(unique(sel$SEX), collapse = ", "), "\n")

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, file.path(out_dir, "tab_case_characteristics.csv"),
                 row.names = FALSE)
writeLines(age_note, file.path(out_dir, "tab_case_characteristics_notes.txt"),
           useBytes = TRUE)
cat("Saved:", file.path(out_dir, "tab_case_characteristics.csv"), "\n")
cat("Saved:", file.path(out_dir, "tab_case_characteristics_notes.txt"), "\n")
