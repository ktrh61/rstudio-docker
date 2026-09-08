# supplementary_data_2_set_level_results.R  (Supplementary Data 2)
# Complete disclosure of the 420 set-level results: every unit x family x set
# with size, ES, NES, p, q_bh (and redundancy flag). Formatting only -- a
# concatenation of the frozen enrichment rds; no computation. List columns
# (e.g. leading_edge) are collapsed to ";"-joined strings.
# Input : processed/thyr_enrichment_test.rds (from 420)
# Output: output/tables/supplementary_data_2.csv

source("setup.R")

en <- readRDS(file.path(paths$processed, "thyr_enrichment_test.rds"))
flatten <- function(df) {
  for (nm in names(df)) {
    if (is.list(df[[nm]])) {
      df[[nm]] <- vapply(df[[nm]], function(v) paste(v, collapse = ";"), "")
    }
  }
  df
}
rows <- lapply(names(en$units), function(u) {
  sets <- en$units[[u]]
  cbind(unit = u, flatten(as.data.frame(sets)))
})
tab <- do.call(rbind, rows)
cat("rows:", nrow(tab), " columns:", paste(names(tab), collapse = ", "), "\n")
cat("q_bh < 0.10 rows (expect 0):", sum(tab$q_bh < 0.10), "\n")

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, file.path(out_dir, "supplementary_data_2.csv"),
                 row.names = FALSE)
cat("Saved:", file.path(out_dir, "supplementary_data_2.csv"), "\n")
