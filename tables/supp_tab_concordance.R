# supp_tab_concordance.R  (Supp.Tab.3(仮))
# Between-arm concordance of the exposure contrast: the two pairs (normal,
# tumor) as two rows -- rho, shuffle reference band, two-sided p, shared
# genes. Formatting only; values must match N-33, N-34 on first run.
# The rds layout predates this script; fields are located defensively and
# printed so the first run documents the mapping.
# Input : diagnostics/output/signature_agreement.rds
# Output: output/tables/supp_tab_concordance.csv (+ printed structure)

source("setup.R")

sa <- readRDS(file.path(paths$root, "diagnostics", "output",
                        "signature_agreement.rds"))
cat("top-level names:", paste(names(sa), collapse = ", "), "\n")
utils::str(sa, max.level = 2)

pick <- function(el, keys) {
  out <- lapply(keys, function(k) if (!is.null(el[[k]])) el[[k]] else NA)
  names(out) <- keys
  out
}
rows <- list()
for (nm in names(sa)) {
  el <- sa[[nm]]
  if (!is.list(el)) next
  v <- pick(el, c("rho", "n_genes", "band_lo", "band_hi", "p", "p_value",
                  "null_lo", "null_hi"))
  if (all(is.na(unlist(v)))) next
  rows[[length(rows) + 1]] <- cbind(pair = nm, as.data.frame(v))
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
