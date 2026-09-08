# table_s7_reo_panel.R  (Table S7)
# Ten-pair relative-expression-ordering panel: pair identifiers, higher- and
# lower-expressed genes, construction-set median difference, reversal rate and
# r0_q10 (the 10th percentile of |log2 TPM difference| among dose-zero samples),
# in panel order. Formatting only.
# Input : processed/thyr_reo_panel.rds (from 520)
# Output: output/tables/table_s7.csv

source("setup.R")

panel <- readRDS(file.path(paths$processed, "thyr_reo_panel.rds"))$panel
tab <- panel[, c("pair_id", "up", "up_name", "down", "down_name",
                 "median_diff", "reversal_rate", "r0_q10")]
print(tab, row.names = FALSE)

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "table_s7.csv")
utils::write.csv(tab, out, row.names = FALSE)
cat("Saved:", out, "\n")
