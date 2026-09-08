# table_s5_ora_annotation.R  (Table S5)
# Descriptive over-representation annotation of the RET-tumor q<0.10 list:
# every family x list x set (18,576 rows = 6,192 sets x higher/lower/combined),
# ordered by family, list and hypergeometric p. Formatting only; the
# family x list counts at q_bh<0.10 must match N-59 on first run.
# Input : processed/thyr_deg_ora_annotation.rds (from 430)
# Output: output/tables/table_s5.csv (+ printed counts)

source("setup.R")

ora <- readRDS(file.path(paths$processed, "thyr_deg_ora_annotation.rds"))
tb <- ora$table
tb <- tb[order(tb$family, tb$list, tb$p_hyper),
         c("family", "list", "pathway", "set_size", "overlap", "expected", "p_hyper", "q_bh")]
tb$expected <- round(tb$expected, 4)

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "table_s5.csv")
utils::write.csv(tb, out, row.names = FALSE)

n_hit <- with(tb[tb$q_bh < FDR_CUT, ], table(family, list))
cat("rows:", nrow(tb), "\n")
print(n_hit)
cat("Saved:", out, "\n")
