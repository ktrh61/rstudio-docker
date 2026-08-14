# tab_case_characteristics.R  (Tab.2(仮))
# Case characteristics by analysis group: main BM (R_Sporadic/R_High/
# B_Sporadic/B_High) and REO evaluation (R_Low/R_Mid). Formatting only --
# no computation beyond counts/medians of frozen inputs; canonical sources
# for the ledger remain the primary rds (N-09, N-12, N-13). CSV is a
# typesetting convenience. Age-difference footnote values (N-64, N-65) come
# from diagnostics/output/age_arm_difference.rds and are appended as notes.
# Input : processed/thyr_analysis_cohorts.rds, processed/thyr_clinical.rds
# Output: output/tables/tab_case_characteristics.csv (+ printed table)

source("setup.R")

co <- readRDS(file.path(paths$processed, "thyr_analysis_cohorts.rds"))
cl <- readRDS(file.path(paths$processed, "thyr_clinical.rds"))
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
rows <- lapply(split(sel, sel$group), function(s) {
  drv <- sort(table(s$Designated_Driver), decreasing = TRUE)
  data.frame(
    group = s$group[1], n = nrow(s),
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

if (any(is.na(sel$SEX))) cat("NOTE: SEX has NA values; recheck coding\n")
cat("sex coding seen:", paste(unique(sel$SEX), collapse = ", "), "\n")

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, file.path(out_dir, "tab_case_characteristics.csv"),
                 row.names = FALSE)
cat("Saved:", file.path(out_dir, "tab_case_characteristics.csv"), "\n")
