# verify_environment.R — build fails unless the image reproduces the
# qualified environment (reorg plan v2 B.12 and phase-6 unit):
#   1. R is 4.5.3 and links the REFERENCE BLAS/LAPACK (never OpenBLAS)
#   2. every explicitly installed package matches the proven version
#      (versions.tsv = parity targets measured in the 4-1b dev container)
#   3. every explicit package actually loads
#   4. the Rcpp toolchain compiles (310/410 build stat_brunnermunzel.cpp
#      at run time)

stopifnot(identical(paste(R.version$major, R.version$minor, sep = "."), "4.5.3"))

blas <- extSoftVersion()[["BLAS"]]
lapack <- La_library()
cat("BLAS  :", blas, "\nLAPACK:", lapack, "\n")
if (grepl("openblas", blas, ignore.case = TRUE) ||
    grepl("openblas", lapack, ignore.case = TRUE)) {
  stop("OpenBLAS is linked; this image requires the reference BLAS/LAPACK.")
}
stopifnot(grepl("/blas/", blas), grepl("/lapack/", lapack))

versions <- read.delim("/tmp/docker-build/versions.tsv",
  header = TRUE, stringsAsFactors = FALSE, comment.char = "#"
)
bad <- character(0)
for (i in seq_len(nrow(versions))) {
  p <- versions$package[i]
  # packageVersion() normalizes "-" to "."; compare in that normal form
  want <- gsub("-", ".", versions$version[i], fixed = TRUE)
  have <- as.character(packageVersion(p))
  if (!identical(have, want)) bad <- c(bad, sprintf("%s: %s != %s", p, have, want))
}
if (length(bad) > 0) stop("Version mismatch:\n", paste(bad, collapse = "\n"))

for (p in versions$package) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
cat("Loaded", nrow(versions), "packages at pinned versions.\n")

stopifnot(Rcpp::evalCpp("1 + 1") == 2)
cat("Rcpp toolchain OK.\n")
cat("Environment verification passed.\n")
