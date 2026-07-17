# setup.R — paths definition and working-directory check. No package loading.

# Treat the current working directory as the repository root.
.root <- getwd()
.markers <- c("Dockerfile", "setup.R")
.missing <- .markers[!file.exists(file.path(.root, .markers))]
if (length(.missing) > 0) {
  stop(
    "Not at repository root: ", .root,
    " (missing: ", paste(.missing, collapse = ", "), ")"
  )
}

paths <- list(
  root      = .root,
  raw       = file.path(.root, "raw"),
  processed = file.path(.root, "processed"),
  output    = file.path(.root, "output"),
  meta      = file.path(.root, "meta")
)

for (.p in paths[c("raw", "processed", "output", "meta")]) {
  if (!dir.exists(.p)) dir.create(.p, recursive = TRUE, showWarnings = FALSE)
}

rm(.root, .markers, .missing, .p)
