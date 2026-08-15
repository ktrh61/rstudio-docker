# supp_tab_concordance.R  (Supp.Tab.3(仮))
# Between-arm concordance of the exposure contrast: the two pairs (normal,
# tumor) as two rows -- rho, shuffle reference band, two-sided p, shared
# genes. Formatting only; values must match N-33, N-34 on first run.
# The rds layout predates this script; the first run (2026-08-15) documented
# the actual structure: sa$pairs$normal / sa$pairs$tumor, each with units,
# n_shared, rho, rho_null (9,999 values), p_two_sided, n_perm. The reference
# interval is the central 95% of rho_null (2.5/97.5 percentiles) -- verified
# to reproduce N-33 ([-0.3914, +0.3930]) and N-34 ([-0.4615, +0.4580]).
# Input : diagnostics/output/signature_agreement.rds
# Output: output/tables/supp_tab_concordance.csv (+ printed structure)

source("setup.R")

sa <- readRDS(file.path(paths$root, "diagnostics", "output",
                        "signature_agreement.rds"))
cat("top-level names:", paste(names(sa), collapse = ", "), "\n")
utils::str(sa, max.level = 2)

rows <- list()
for (nm in names(sa$pairs)) {
  el <- sa$pairs[[nm]]
  q <- quantile(el$rho_null, c(0.025, 0.975))
  rows[[length(rows) + 1]] <- data.frame(
    pair = nm,
    units = paste(el$units, collapse = " x "),
    n_shared_genes = el$n_shared,
    rho = el$rho,
    interval_lo = unname(q[1]),
    interval_hi = unname(q[2]),
    p_two_sided = el$p_two_sided,
    n_perm = el$n_perm
  )
}
if (length(rows) == 0) {
  cat("NOTE: field names differ from guesses -- read str() above and adjust\n")
} else {
  tab <- do.call(rbind, rows)
  print(tab)
  out_dir <- file.path(paths$output, "tables")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(tab, file.path(out_dir, "supp_tab_concordance.csv"),
                   row.names = FALSE)
  cat("Saved:", file.path(out_dir, "supp_tab_concordance.csv"), "\n")
}
