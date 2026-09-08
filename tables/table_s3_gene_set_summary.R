# table_s3_gene_set_summary.R  (Table S3)
# Set-level (420) summary per contrast x collection: sets tested after the
# size filter and the minimum within-collection q_bh. Formatting only; values
# must match N-27, N-28 on first run.
# Input : processed/thyr_enrichment_test.rds (from 420)
# Output: output/tables/table_s3.csv (+ printed table)

source("setup.R")

en <- readRDS(file.path(paths$processed, "thyr_enrichment_test.rds"))

rows <- list()
for (u in names(en$units)) {
  sets <- en$units[[u]]
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

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "table_s3.csv")
utils::write.csv(tab, out, row.names = FALSE)
cat("Saved:", out, "\n")
