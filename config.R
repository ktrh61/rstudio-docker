# config.R — constants that must agree across scripts. Loaded by setup.R, so
# every script that sources setup.R sees these. Stage-specific parameters stay
# in their own script and are recorded in that script's output config.

# Canonical seed (Chernobyl accident date 1986-04-26). Used by 410 (label
# shuffles) and 530 (REO Low-vs-Mid Monte Carlo BM); diagnostics draw their
# own documented seeds.
SEED <- 19860426L

# Parallelism: the canonical declaration. The publication run executes on the
# adoption machine (Xeon, 4C/8T), and a distributed verification script must
# not presuppose more. No stage nests these knobs -- the peak concurrency is
# 4 at every step, and the reference BLAS adds no threads of its own.
WORKERS <- 4L

# Brunner-Munzel exact-enumeration threads and allocation cap (310/410).
EXACT_THREADS <- 4L
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

# Protocol-wide inference threshold: Storey q < FDR_CUT at the gene level
# (410) and the set level (420). Also the display cutoff for figures/messages.
FDR_CUT <- 0.10

# BH cutoff for the DEGES potential-DEG screen inside 310. Conceptually
# distinct from FDR_CUT despite the shared value.
DEGES_FDR <- 0.10

# Label shuffles for the empirical null (410; 420 consumes 410's saved
# per-unit perm_index so the shuffles always match the artifact under
# analysis). Canonical resolution per reorg plan v2 D4: simulation noise must
# not decide claims, and the per-set p floor must clear the BH boundary of the
# smallest family.
N_PERM <- 9999L
