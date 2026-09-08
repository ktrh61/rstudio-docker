# table_s6_between_stratum_concordance.R  (Table S6)
# Between-stratum concordance of the exposure contrast, one row per tissue
# pair (normal, tumor): Spearman rho over shared genes, the central 95% of the
# shuffle reference (2.5/97.5 percentiles of rho_null), two-sided p and the
# number of shuffles. Formatting only; values must match N-33, N-34 on first
# run.
# Input : processed/thyr_signature_agreement.rds (from 440)
# Output: output/tables/table_s6.csv (+ printed table)

source("setup.R")

sa <- readRDS(file.path(paths$processed, "thyr_signature_agreement.rds"))

rows <- lapply(names(sa$pairs), function(nm) {
  el <- sa$pairs[[nm]]
  q <- quantile(el$rho_null, c(0.025, 0.975))
  data.frame(
    pair = nm,
    units = paste(el$units, collapse = " x "),
    n_shared_genes = el$n_shared,
    rho = el$rho,
    interval_lo = unname(q[1]),
    interval_hi = unname(q[2]),
    p_two_sided = el$p_two_sided,
    n_perm = el$n_perm
  )
})
tab <- do.call(rbind, rows)
print(tab)

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "table_s6.csv")
utils::write.csv(tab, out, row.names = FALSE)
cat("Saved:", out, "\n")
