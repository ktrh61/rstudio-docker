# signature_agreement.R
# Between-driver-arm signature agreement, normal and tumor tissue (reorg plan
# v2 s0.5 6th confirmation and appendix B.7; reading rules for both pairs
# fixed there before the corresponding first run).
#
# For each tissue, the exposed-vs-sporadic contrast is summarized per gene by
# the signed Brunner-Munzel statistic, separately in the RET arm and the BRAF
# arm, and the two genome-wide profiles are compared by Spearman rank
# correlation over shared genes -- threshold-free, no dose assumption.
#
# The normal-tissue pair is the hypothesis-bearing comparison (shared
# glandular exposure memory; Abend/Ory). The tumor pair is a descriptive
# symmetric completion only: the design's premise is that tumor expression is
# dominated by driver biology, so no similarity between the two
# driver-conditioned contrasts is predicted, and no outcome tests the core
# hypothesis (appendix B.7).
#
# Each point estimate travels with a permutation-calibrated reference (same
# convention as the pi0 reporting, v2 D1): labels are shuffled independently
# within each unit and the shuffle-b vectors of the two arms are correlated,
# giving the null spread of rho when neither cohort carries label-aligned
# structure. The shuffles are drawn HERE with a distinct seed per unit --
# 410's saved perm_index is deliberately not reused, because 410 seeds every
# unit identically (equal-n units hold bit-identical index matrices; harmless
# for 410/420's per-unit inference, fatal for a null whose premise is
# cross-unit independence). Inter-gene correlation makes the null spread far
# wider than 1/sqrt(n_genes), which is why a bare coefficient is
# uninterpretable without it. An observed rho outside the spread says the two
# contrasts share label-aligned structure; it does NOT by itself say the
# structure is the exposure trace rather than a shared covariate (age
# structure, purity) -- that separation belongs to the covariate diagnostics
# (s0.5 2nd).
# Input : processed/thyr_expression_test.rds    (from 410; effect, statistic,
#         perm_index)
#         processed/thyr_normalized_counts.rds  (from 310; per-unit DGEList)
# Output: diagnostics/output/signature_agreement.rds

source("setup.R")

suppressPackageStartupMessages({
  library(edgeR)
  library(parallel)
})

source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))
source(file.path(paths$root, "lib", "gsea_permutation.R")) # bind helper
source(file.path(paths$root, "lib", "units.R"))

pin_blas_threads()

PAIRS <- list(
  normal = c("R_Normal", "B_Normal"),
  tumor = c("R_Tumor", "B_Tumor")
)
# One distinct shuffle seed per unit (see header). Base = the shared
# diagnostic base seed (v2 B.10); offset = unit position.
AGREEMENT_SEED_BASE <- 19450809L

expression_test <- readRDS(
  file.path(paths$processed, "thyr_expression_test.rds")
)
normalized <- readRDS(
  file.path(paths$processed, "thyr_normalized_counts.rds")
)

signed_vector <- function(unit) {
  g <- expression_test$units[[unit]]$genes
  setNames(g$statistic * sign(g$effect - 0.5), g$gene_id)
}

unit_cpm <- function(unit) {
  dgelist <- normalized$units[[unit]]$dgelist
  arms <- unit_arms(dgelist$samples$group, unit)
  m <- edgeR::cpm(
    dgelist,
    normalized.lib.sizes = TRUE, prior.count = 0, log = FALSE
  )[, c(arms$sporadic, arms$high), drop = FALSE]
  storage.mode(m) <- "double"
  list(cpm = m, nx = length(arms$sporadic))
}

unit_perm_index <- function(unit, n) {
  set.seed(AGREEMENT_SEED_BASE + match(unit, names(normalized$units)))
  vapply(seq_len(N_PERM), function(i) sample(n), integer(n))
}

compare_pair <- function(pair_name, units) {
  vec_1 <- signed_vector(units[1])
  vec_2 <- signed_vector(units[2])
  shared <- intersect(names(vec_1), names(vec_2))
  rho <- stats::cor(vec_1[shared], vec_2[shared], method = "spearman")

  data_1 <- unit_cpm(units[1])
  data_2 <- unit_cpm(units[2])
  perm_1 <- unit_perm_index(units[1], ncol(data_1$cpm))
  perm_2 <- unit_perm_index(units[2], ncol(data_2$cpm))
  n_perm <- N_PERM
  rho_null <- as.numeric(gsea_bind_null_columns(
    mclapply(seq_len(n_perm), function(i) {
      v1 <- brunnermunzel_statistics(
        data_1$cpm[shared, perm_1[, i], drop = FALSE], data_1$nx
      )
      v2 <- brunnermunzel_statistics(
        data_2$cpm[shared, perm_2[, i], drop = FALSE], data_2$nx
      )
      stats::cor(v1, v2, method = "spearman")
    }, mc.cores = WORKERS),
    1L
  ))

  p_two <- (sum(abs(rho_null) >= abs(rho)) + 1) / (n_perm + 1)
  message(sprintf(
    "%-6s %s x %s : rho = %+.4f over %d shared genes",
    pair_name, units[1], units[2], rho, length(shared)
  ))
  message(sprintf(
    "       shuffle null: median %+.4f, 95%% band [%+.4f, %+.4f] ; two-sided p %.4f",
    stats::median(rho_null),
    stats::quantile(rho_null, 0.025), stats::quantile(rho_null, 0.975),
    p_two
  ))
  list(
    units = units, n_shared = length(shared), rho = rho,
    rho_null = rho_null, p_two_sided = p_two, n_perm = n_perm
  )
}

results <- lapply(names(PAIRS), function(nm) compare_pair(nm, PAIRS[[nm]]))
names(results) <- names(PAIRS)

out_dir <- file.path(paths$root, "diagnostics", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out <- file.path(out_dir, "signature_agreement.rds")
saveRDS(
  list(
    date = Sys.Date(),
    config = list(
      metric = "signed BM statistic (statistic * sign(effect - 0.5))",
      method = "spearman",
      null = "paired label shuffles, one distinct seed per unit",
      seed_base = AGREEMENT_SEED_BASE,
      n_perm = N_PERM,
      readings = "fixed in reorg plan v2 s0.5 6th / appendix B.7 before running"
    ),
    pairs = results
  ),
  out
)
message("Saved: ", out)
