# table_s4_null_calibration.R  (Table S4)
# Held-out complete-null calibration of the set-level procedure, one row per
# contrast x collection cell: sets tested, replicates, replicates with at
# least one discovery, P(>=1 discovery) with exact binomial 95% CI, and the
# mean and maximum discoveries per replicate. Formatting only; values must
# match N-24 on first run.
# Input : processed/thyr_gsea_null_calibration.rds (from 415)
# Output: output/tables/table_s4.csv (+ printed table)

source("setup.R")

cal <- readRDS(file.path(paths$processed, "thyr_gsea_null_calibration.rds"))
tab <- cal$summary
print(tab, digits = 3)

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "table_s4.csv")
utils::write.csv(tab, out, row.names = FALSE)
cat("Saved:", out, "\n")
