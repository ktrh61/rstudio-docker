# install_packages.R — build-time package installation (P3M snapshot pinned
# via Rprofile.site). Explicit set = recommended tier + packages the active
# pipeline code references (extraction recorded in reorg plan v2 phase-6
# unit); dependencies resolve from the same snapshot.

ncpus <- parallel::detectCores()

recommended <- c(
  "boot", "class", "cluster", "codetools", "foreign", "KernSmooth",
  "lattice", "MASS", "Matrix", "mgcv", "nlme", "nnet", "rpart",
  "spatial", "survival"
)
cran <- c(
  "R.utils", "Rcpp", "RhpcBLASctl", "brunnermunzel", "data.table",
  "digest", "doSNOW", "foreach", "iterators", "ggplot2", "ggrepel",
  "matrixStats", "msigdbr", "testthat", "BiocManager"
)
bioc <- c(
  "edgeR", "limma", "fgsea", "TCC", "SummarizedExperiment",
  "GenomicRanges", "GenomicDataCommons", "rtracklayer"
)

install.packages(c(recommended, cran), Ncpus = ncpus)
BiocManager::install(bioc, update = FALSE, ask = FALSE, Ncpus = ncpus)

missing <- setdiff(
  c(recommended, cran, bioc),
  rownames(installed.packages())
)
if (length(missing) > 0) {
  stop("Packages failed to install: ", paste(missing, collapse = ", "))
}
