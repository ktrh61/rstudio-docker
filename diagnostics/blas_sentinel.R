# blas_sentinel.R
# Decision material for reorg plan v2 appendix C (BLAS policy A/B/C; the
# pending half of docs/decisions.md D-020). Runs MUREN's step-1 scaling on
# the real R_Tumor unit counts and fingerprints the coefficient vector, so
# two amd64 machines can be compared bit-for-bit.
#
# Run TWICE per machine, in the project container, from the repo root:
#   Rscript diagnostics/blas_sentinel.R
#   OPENBLAS_CORETYPE=HASWELL Rscript diagnostics/blas_sentinel.R
# Compare the printed digests across machines. input_digest must match first
# (same processed artifact); coeff_digest equality is the measurement.
# Removable after the BLAS decision is recorded.

source("setup.R")

source(file.path(paths$root, "lib", "norm_muren_helpers.R"))
source(file.path(paths$root, "lib", "norm_muren.R"))

pin_blas_threads()

normalized <- readRDS(file.path(paths$processed, "thyr_normalized_counts.rds"))
counts <- normalized$units$R_Tumor$dgelist$counts

sc <- muren_norm(
  counts,
  refs = "saturated", pairwise_method = "lts", single_param = TRUE,
  res_return = "scaling_coeff", workers = 4L
)

cpu <- grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]
cat("---- blas_sentinel ----\n")
cat("cpu           :", sub(".*: ", "", cpu), "\n")
cat("coretype_env  :", Sys.getenv("OPENBLAS_CORETYPE", unset = "(unset)"), "\n")
cat("blas          :", extSoftVersion()[["BLAS"]], "\n")
cat("la_library    :", La_library(), "\n")
cat("input_digest  :", digest::digest(counts), "\n")
cat("coeff_digest  :", digest::digest(sc), "\n")
cat("coeff_head    :", paste(sprintf("%.17g", head(sc, 3)), collapse = " "), "\n")
