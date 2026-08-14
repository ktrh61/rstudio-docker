# tab_gene_level_summary.R  (Tab.3(仮))
# Gene-level (410) summary per unit: tested genes, pi0, DEG counts and
# directions, primary omnibus HC p. Formatting only -- values must match
# N-15, N-16, N-17, N-18, N-19, N-20 on first run (consistency check).
# Input : processed/thyr_expression_test.rds
# Output: output/tables/tab_gene_level_summary.csv (+ printed table)

source("setup.R")

x <- readRDS(file.path(paths$processed, "thyr_expression_test.rds"))
rows <- lapply(names(x$units), function(u) {
  un <- x$units[[u]]
  g <- un$genes
  deg <- g$q_storey < 0.10
  hc <- un$omnibus[un$omnibus$test == "hc", ]
  data.frame(
    unit = u,
    n_tested = nrow(g),
    pi0 = round(un$pi0$estimate, 3),
    deg_q10 = sum(deg),
    up = sum(deg & g$effect > 0.5),
    down = sum(deg & g$effect < 0.5),
    min_p_exact = format(min(g$p_exact), digits = 3),
    hc_p = hc$p[1]
  )
})
tab <- do.call(rbind, rows)
print(tab)

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, file.path(out_dir, "tab_gene_level_summary.csv"),
                 row.names = FALSE)
cat("Saved:", file.path(out_dir, "tab_gene_level_summary.csv"), "\n")
