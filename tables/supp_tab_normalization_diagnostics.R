# supp_tab_normalization_diagnostics.R
# Per-contrast expression filtering and DEGES-MUREN normalization diagnostics.
# Formatting only: reads the frozen normalization object and does not refit,
# reshuffle, renormalize, or retest any data.
# Input : processed/thyr_normalized_counts.rds
# Output: output/tables/supp_tab_normalization_diagnostics.csv

source("setup.R")

obj <- readRDS(file.path(paths$processed, "thyr_normalized_counts.rds"))

rows <- lapply(names(obj$units), function(unit) {
  item <- obj$units[[unit]]
  d <- item$diagnostics
  nf <- item$dgelist$samples$norm.factors
  final_iteration <- length(d$iter_pi0)

  data.frame(
    contrast = unit,
    reference_group = d$groups[1],
    high_group = d$groups[2],
    n_reference = d$n_samples[1],
    n_high = d$n_samples[2],
    protein_coding_genes = d$n_protein_coding,
    genes_after_filterByExpr = d$n_after_filter,
    deges_iterations = d$n_iterations,
    final_screen_pi0 = round(d$iter_pi0[final_iteration], 3),
    final_jaccard = round(d$iter_jaccard[final_iteration], 3),
    norm_factor_min = round(min(nf), 4),
    norm_factor_max = round(max(nf), 4),
    stringsAsFactors = FALSE
  )
})

tab <- do.call(rbind, rows)
print(tab)

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "supp_tab_normalization_diagnostics.csv")
utils::write.csv(tab, out_path, row.names = FALSE)
message("Saved: ", out_path)
