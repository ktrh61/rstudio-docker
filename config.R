# config.R — constants that must agree across scripts. Loaded by setup.R, so
# every script that sources setup.R sees these. Stage-specific parameters stay
# in their own script and are recorded in that script's output config.

# Canonical seed (Chernobyl accident date 1986-04-26). Used by 310 (DEGES MC
# fallback), 410 (label shuffles), 530 (REO Low-vs-Mid BM), and diagnostics.
SEED <- 19860426L

# MUREN parallel workers. Development value; the canonical publication run uses
# 4L (worker count does not affect results; switch at the final clean run).
WORKERS <- 16L

# Brunner-Munzel exact-enumeration threads and allocation cap (310/410/420).
EXACT_THREADS <- 16L
BM_EXACT_MAX <- 1e8

# Pin BLAS/OpenMP to one thread so explicit parallelism is the only parallelism.
pin_blas_threads <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
  invisible(NULL)
}

# Assigned-share band boundaries (percent). Band rule (single convention;
# boundary cases do not occur in this cohort, verified 2026-07-28):
#   dose_mgy == 0                     -> Sporadic (AS not required)
#   dose_mgy > 0 & 0 < AS <= 33.3    -> Low
#   dose_mgy > 0 & 33.3 < AS < 66.6  -> Mid
#   dose_mgy > 0 & AS >= 66.6        -> High
AS_LOW_MAX <- 33.3
AS_HIGH_MIN <- 66.6

# Minimum pooled (common-scale) relative tumor purity for the main BM cohort.
PURITY_THRESHOLD <- 0.6

# REO dead zone: |log2 TPM difference| below this does not count as an order.
DEAD_ZONE <- log2(1.2)

# Display FDR cutoff for gene-level figures/messages (permutation FDR).
FDR_CUT <- 0.10

# BH cutoff for the DEGES potential-DEG screen inside 310. Conceptually
# distinct from FDR_CUT despite the shared value.
DEGES_FDR <- 0.10

# Label shuffles for the empirical null (410; 420 inherits the value from
# 410's saved config so it always matches the artifact under analysis).
N_PERM <- 999L
