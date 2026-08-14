# supp_tab_package_versions.R  (Supp.Tab.4(仮))
# Package-version table for the Methods software statement: the pinned
# versions from the canonical container build (docker/versions.tsv).
# Formatting only -- reads the tracked build input verbatim. The run-time
# session record (run/xeon_provenance/session_info.txt) lives outside the
# repository; assembly can cite it alongside (N-02).
# Input : docker/versions.tsv
# Output: output/tables/supp_tab_package_versions.csv (+ printed table)

source("setup.R")

tsv <- file.path(paths$root, "docker", "versions.tsv")
lines <- readLines(tsv)
lines <- lines[!grepl("^#", lines) & nzchar(lines)]
tab <- utils::read.delim(text = paste(lines, collapse = "\n"),
                         stringsAsFactors = FALSE)
tab <- tab[order(tab$package), ]
print(tab, row.names = FALSE)
cat("packages:", nrow(tab), "\n")

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, file.path(out_dir, "supp_tab_package_versions.csv"),
                 row.names = FALSE)
cat("Saved:", file.path(out_dir, "supp_tab_package_versions.csv"), "\n")
