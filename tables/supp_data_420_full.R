# supp_data_420_full.R  (Supp.Data.1(仮))
# Complete disclosure of the 420 set-level results: every unit x family x set
# with size, ES, NES, p, q_bh (and redundancy flag). Formatting only -- a
# concatenation of the frozen enrichment rds; no computation. List columns
# (e.g. leading_edge) are collapsed to ";"-joined strings.
# Input : processed/thyr_enrichment_test.rds
# Output: output/tables/supp_data_420_full.csv

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
  sets <- en$units[[u]]$sets
  if (is.null(sets)) sets <- en$units[[u]]
  cbind(unit = u, flatten(as.data.frame(sets)))
})
tab <- do.call(rbind, rows)
cat("rows:", nrow(tab), " columns:", paste(names(tab), collapse = ", "), "\n")
cat("q_bh < 0.10 rows (expect 0):", sum(tab$q_bh < 0.10), "\n")

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, file.path(out_dir, "supp_data_420_full.csv"),
                 row.names = FALSE)
cat("Saved:", file.path(out_dir, "supp_data_420_full.csv"), "\n")
