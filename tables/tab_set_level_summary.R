# tab_set_level_summary.R  (Tab.4(仮))
# Set-level (420) summary per unit x collection -- sets tested and min q_bh --
# joined with the D6 held-out calibration per-cell P(>=1 discovery) and CI.
# Formatting only; values must match N-24, N-27, N-28 on first run.
# The calibration cell table is located defensively (first data.frame in the
# rds carrying a p_any column) because the rds layout predates this script.
# Input : processed/thyr_enrichment_test.rds
#         diagnostics/output/gsea_null_calibration.rds
# Output: output/tables/tab_set_level_summary.csv (+ printed table)

source("setup.R")

en <- readRDS(file.path(paths$processed, "thyr_enrichment_test.rds"))
cal <- readRDS(file.path(paths$root, "diagnostics", "output",
                         "gsea_null_calibration.rds"))

find_cells <- function(obj) {
  if (is.data.frame(obj) && "p_any" %in% names(obj)) return(obj)
  if (is.list(obj)) {
    for (el in obj) {
      r <- find_cells(el)
      if (!is.null(r)) return(r)
    }
  }
  NULL
}
cells <- find_cells(cal)
if (is.null(cells)) stop("calibration cell table (column p_any) not found")
cat("calibration table columns:", paste(names(cells), collapse = ", "), "\n")

rows <- list()
for (u in names(en$units)) {
  sets <- en$units[[u]]$sets
  if (is.null(sets)) sets <- en$units[[u]]
  for (col in unique(sets$collection)) {
    blk <- sets[sets$collection == col, ]
    rows[[length(rows) + 1]] <- data.frame(
      unit = u, collection = col,
      n_sets = nrow(blk),
      min_q_bh = round(min(blk$q_bh), 3)
    )
  }
}
tab <- do.call(rbind, rows)
print(tab)
cat("\nD6 calibration cells (join by unit/collection at assembly):\n")
print(cells)

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, file.path(out_dir, "tab_set_level_summary.csv"),
                 row.names = FALSE)
utils::write.csv(cells, file.path(out_dir, "tab_set_level_calibration.csv"),
                 row.names = FALSE)
cat("Saved: 2 csv under", out_dir, "\n")
