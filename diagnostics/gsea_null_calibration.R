# gsea_null_calibration.R
# Held-out null calibration of the spec-B gene-set inference (reorg plan v2
# D6, decided 2026-08-07). A linear diagnostic that reports a measurement; it
# gates nothing by mechanism. The pre-registered constraint is ordering only:
# 420 is not run on real data until this measurement is consistent with the
# nominal level, and any inconsistency is diagnosed and resolved by amendment,
# not by code.
#
# Design (shared-null leave-out): with a seed independent of the inference
# seed, draw B + R label shuffles per unit and compute the block ES for all of
# them once. Shuffles 1..B form one shared null pool; each of the R later
# shuffles is treated as a pseudo-observation against that pool and pushed
# through the exact inference 420 applies (NES against the pool, per-set
# permutation p, BH within collection, q < FDR_CUT). Under the global null
# every shuffle is exchangeable, so the extra cost over one analysis is
# ~(B + R)/B. Because the replicates share the null pool, their discovery
# indicators are weakly positively correlated; the binomial CI understates
# that, which is recorded here rather than corrected.
#
# History (reorg plan v2 appendix B): this measurement is what replaced the
# pooled tail-ratio FDR of the original D2 -- measured miscalibrated here
# (P(>=1) 0.14 pooled, 0.24 worst cell; a restandardized variant 0.22/0.44)
# before any real-data run; per-set p + BH measured 0.045 and was adopted.
# The superseded measurement logs stay in diagnostics/output/.
#
# Under the complete null FDR = P(at least one discovery), so the measurement
# per unit x collection is: the share of replicates with >= 1 discovery
# (with an exact binomial CI) and the mean discovery count. Consistency with
# the nominal 0.10 is judged by a human reading the table (reorg plan v2
# s3.0-3/-4); this doubles as the integration test of the D2 + D3 machinery.
# Input : processed/thyr_normalized_counts.rds  (from 310)
#         lib/stat_brunnermunzel.R, lib/gsea_permutation.R,
#         lib/gsea_collections.R
# Output: diagnostics/output/gsea_null_calibration.rds
#           list(date, config, summary, counts)

source("setup.R")

suppressPackageStartupMessages({
  library(edgeR)
  library(msigdbr)
  library(parallel)
})

source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))
source(file.path(paths$root, "lib", "gsea_permutation.R"))
source(file.path(paths$root, "lib", "gsea_collections.R"))
source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "annotation.R"))

# --- Configuration ---------------------------------------------------------
SPECIES <- "Homo sapiens"
CALIB_SEED <- 20260807L # decision date of D6; independent of SEED
R_REPLICATES <- 100L
# B (the shared null pool) matches the analysis N_PERM from config.R, so the
# calibration measures the machinery at the resolution 420 actually uses.

pin_blas_threads()

# --- Load inputs -----------------------------------------------------------
norm_path <- file.path(paths$processed, "thyr_normalized_counts.rds")
if (!file.exists(norm_path)) stop("thyr_normalized_counts.rds not found (310)")
normalized <- readRDS(norm_path)

collections <- load_gene_set_collections(SPECIES)
gene_sets <- collections$gene_sets

out_dir <- file.path(paths$root, "diagnostics", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Per-unit calibration --------------------------------------------------
calibrate_unit <- function(dgelist, unit) {
  arms <- unit_arms(dgelist$samples$group, unit)
  cpm_matrix <- edgeR::cpm(
    dgelist,
    normalized.lib.sizes = TRUE, prior.count = 0, log = FALSE
  )[, c(arms$sporadic, arms$high), drop = FALSE]
  storage.mode(cpm_matrix) <- "double"
  ensembl <- strip_ensembl_version(rownames(cpm_matrix))
  cpm_matrix <- cpm_matrix[!duplicated(ensembl), , drop = FALSE]
  rownames(cpm_matrix) <- ensembl[!duplicated(ensembl)]
  nx <- length(arms$sporadic)
  n <- ncol(cpm_matrix)

  index <- lapply(gene_sets, gsea_pathway_index,
    gene_ids = rownames(cpm_matrix),
    min_size = GSEA_MIN_SET_SIZE, max_size = GSEA_MAX_SET_SIZE
  )
  collection_of <- rep(names(index), lengths(index))
  flat <- unlist(unname(index), recursive = FALSE)

  total <- N_PERM + R_REPLICATES
  set.seed(CALIB_SEED)
  perm_index <- vapply(seq_len(total), function(i) sample(n), integer(n))

  es_columns <- mclapply(seq_len(total), function(i) {
    gsea_block_scores(
      gsea_normal_scores(brunnermunzel_statistics(
        cpm_matrix[, perm_index[, i], drop = FALSE], nx
      )),
      flat
    )
  }, mc.cores = WORKERS)
  es <- gsea_bind_null_columns(es_columns, length(flat))
  rownames(es) <- names(flat)
  pool <- es[, seq_len(N_PERM), drop = FALSE]
  pseudo <- es[, N_PERM + seq_len(R_REPLICATES), drop = FALSE]

  counts_by_collection <- lapply(names(index), function(collection) {
    keep <- collection_of == collection
    vapply(seq_len(R_REPLICATES), function(r) {
      standardized <- gsea_nes(pseudo[keep, r], pool[keep, , drop = FALSE])
      pval <- gsea_pathway_pvalues(standardized$nes, standardized$nes_null)
      sum(p.adjust(pval, method = "BH") < FDR_CUT, na.rm = TRUE)
    }, integer(1))
  })
  names(counts_by_collection) <- names(index)

  rows <- lapply(names(index), function(collection) {
    counts <- counts_by_collection[[collection]]
    n_any <- sum(counts >= 1L)
    ci <- stats::binom.test(n_any, R_REPLICATES)$conf.int
    data.frame(
      unit = unit, collection = collection,
      m_sets = sum(collection_of == collection),
      replicates = R_REPLICATES, n_any_discovery = n_any,
      p_any = n_any / R_REPLICATES,
      ci_lo = ci[1], ci_hi = ci[2],
      mean_discoveries = mean(counts),
      max_discoveries = max(counts),
      stringsAsFactors = FALSE
    )
  })
  list(summary = do.call(rbind, rows), counts = counts_by_collection)
}

message(sprintf(
  "Calibration: B = %d shared null, R = %d replicates, seed %d, q < %.2f",
  N_PERM, R_REPLICATES, CALIB_SEED, FDR_CUT
))
results <- lapply(names(normalized$units), function(unit) {
  message("  Unit ", unit)
  calibrate_unit(normalized$units[[unit]]$dgelist, unit)
})
names(results) <- names(normalized$units)

calibration_summary <- do.call(rbind, lapply(results, `[[`, "summary"))
rownames(calibration_summary) <- NULL
print(calibration_summary, digits = 3)

# --- Save ------------------------------------------------------------------
gsea_null_calibration <- list(
  date = Sys.Date(),
  config = list(
    calib_seed = CALIB_SEED,
    n_null_pool = N_PERM,
    n_replicates = R_REPLICATES,
    q_threshold = FDR_CUT,
    nominal = FDR_CUT,
    design = "shared-null leave-out (reorg plan v2 D6)",
    min_set_size = GSEA_MIN_SET_SIZE,
    max_set_size = GSEA_MAX_SET_SIZE,
    radiation_curation = collections$curation,
    msigdbr_version = as.character(packageVersion("msigdbr"))
  ),
  summary = calibration_summary,
  counts = lapply(results, `[[`, "counts")
)

out <- file.path(out_dir, "gsea_null_calibration.rds")
saveRDS(gsea_null_calibration, out)
message("Saved: ", out)
