# supp_tab_package_versions.R  (Table S2)
# Package-version table for the software statement. Lists the R packages the
# reported work actually loads, derived by scanning the code for library(),
# requireNamespace() and pkg:: calls, and joined to the pinned versions of the
# canonical container build (docker/versions.tsv). Roles:
#   pipeline    -- scripts/, lib/, figures/, tables/, config.R, setup.R
#   diagnostics -- diagnostics/ (reported ancillary analyses)
#   tests       -- tests/testthat (reference implementations and framework)
# Packages pinned in the container but not loaded by any of these files (e.g.
# TCC, whose protocol the in-house DEGES follows but which is not called;
# BiocManager; recommended packages pulled in as dependencies) are omitted:
# the container recipe records the complete environment. Stops if a loaded
# package has no pinned version (researcher decision 2026-08-29).
# Input : docker/versions.tsv + the code tree
# Output: output/tables/supp_tab_package_versions.csv (+ printed table)
source("setup.R")

tsv <- file.path(paths$root, "docker", "versions.tsv")
lines <- readLines(tsv)
lines <- lines[!grepl("^#", lines) & nzchar(lines)]
pinned <- utils::read.delim(text = paste(lines, collapse = "\n"),
                            stringsAsFactors = FALSE)

roles <- list(
  pipeline = c(list.files(file.path(paths$root, "scripts"), "\\.R$", full.names = TRUE),
               list.files(file.path(paths$root, "lib"), "\\.R$", full.names = TRUE),
               list.files(file.path(paths$root, "figures"), "\\.R$", full.names = TRUE),
               list.files(file.path(paths$root, "tables"), "\\.R$", full.names = TRUE),
               file.path(paths$root, c("config.R", "setup.R"))),
  diagnostics = list.files(file.path(paths$root, "diagnostics"), "\\.R$", full.names = TRUE),
  tests = list.files(file.path(paths$root, "tests", "testthat"), "\\.R$", full.names = TRUE)
)
used_in <- function(files) {
  txt <- unlist(lapply(files, readLines, warn = FALSE))
  txt <- txt[!grepl("^\\s*#", txt)]  # comments do not count as use
  hits <- c(
    regmatches(txt, gregexpr("library\\(([A-Za-z][A-Za-z0-9.]*)\\)", txt)),
    regmatches(txt, gregexpr("requireNamespace\\(\"([A-Za-z][A-Za-z0-9.]*)\"", txt)),
    regmatches(txt, gregexpr("\\b([A-Za-z][A-Za-z0-9.]*)::", txt))
  )
  hits <- unlist(hits)
  pk <- gsub("^library\\(|^requireNamespace\\(\"|\\)$|\"$|::$", "", hits)
  unique(pk)
}
found <- lapply(roles, used_in)
base_pkgs <- c(rownames(utils::installed.packages(priority = "base")), "R")
rows <- do.call(rbind, lapply(names(found), function(r) {
  p <- setdiff(found[[r]], base_pkgs)
  if (length(p) == 0) return(NULL)
  data.frame(package = p, role = r, stringsAsFactors = FALSE)
}))
# first role wins in the fixed order pipeline > diagnostics > tests
rows$role <- factor(rows$role, levels = names(roles))
rows <- rows[order(rows$role), ]
rows <- rows[!duplicated(rows$package), ]
rows$role <- as.character(rows$role)

missing <- setdiff(rows$package, pinned$package)
if (length(missing) > 0) {
  stop("Loaded packages without a pinned version: ", paste(missing, collapse = ", "))
}
tab <- merge(rows, pinned[, c("package", "version")], by = "package")
tab <- tab[order(tab$role != "pipeline", tab$role != "diagnostics", tolower(tab$package)), ]
tab <- tab[, c("package", "version", "role")]
print(tab, row.names = FALSE)
cat("packages:", nrow(tab), "| pinned but not loaded:",
    paste(setdiff(pinned$package, tab$package), collapse = ", "), "\n")

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, file.path(out_dir, "supp_tab_package_versions.csv"),
                 row.names = FALSE)
cat("Saved:", file.path(out_dir, "supp_tab_package_versions.csv"), "\n")
