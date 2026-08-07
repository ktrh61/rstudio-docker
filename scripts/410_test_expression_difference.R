# 410_test_expression_difference.R
# Test each analysis unit for a High-vs-Sporadic expression difference on the
# DEGES-normalized counts from 310, with the exact permutation Brunner-Munzel
# test. The products are a complete ranked gene table with Storey q-values
# (the protocol-wide inference, q < 0.10), a unit-level omnibus permutation
# test, and the shared permutation index that 420 reuses.
# Input : processed/thyr_normalized_counts.rds  (from 310; per-unit DGEList)
#         lib/stat_brunnermunzel.R
# Output: processed/thyr_expression_test.rds
#           list(date, config, units)
#           units = { R_Tumor, R_Normal, B_Tumor, B_Normal }; each holds
#             n_samples : Sporadic and High column counts
#             genes     : gene_id, effect, statistic, p_exact, q_storey, rank
#             pi0       : estimate, and the shuffle-null pi0 distribution
#                         (vector + quantiles) that calibrates it
#             q_curve   : rejections R(alpha) over a q-threshold grid
#             omnibus   : test, alpha, cut, observed, null_median, null_lo,
#                         null_hi, p
#             null      : n_perm and quantiles of the pooled null statistic
#             perm_index / perm_index_hash : the label shuffles, shared with
#                         420 by object (reorg plan v2 s3.2)
#
# Group x is Sporadic and group y is High, so effect > 0.5 means the gene runs
# higher in the exposed arm. p_exact enumerates every C(n, nx) allocation and
# needs no seed.
#
# Inference is Storey q on the exact p-values with the plug-in pi0 estimate at
# lambda = 0.5 (q = pi0_hat * BH; reorg plan v2 D1), the same correction the
# DEGES screen in 310 uses. Under a weakly diffuse alternative the empirical
# p-value distribution runs nearly parallel to the correction's self-consistency
# line, so R(alpha) jumps by thousands of genes over a short q interval ("the
# cliff"). That is a property of this signal shape, not a defect: q_curve
# records R(alpha) so reports can show the whole curve instead of a single
# threshold count. pi0 is reported with permutation-calibrated uncertainty:
# each label shuffle is scored against its own gene-wise null (row ranks across
# shuffles) and the plug-in estimator is applied per shuffle, giving the
# distribution the estimate would have under the global null.
#
# The historical permutation-calibrated FDR column (fdr_perm, pi0 = 1) was
# retired with the other sensitivity columns (reorg plan v2 appendix B.2);
# its record lives in the baseline and the amendment record.
#
# The omnibus rows share the shuffle null. "count" asks whether the number of genes
# past a fixed null quantile exceeds what the shuffles produce, "max" asks the
# same of the single most extreme gene, and "hc" is Higher Criticism, which
# scans the threshold instead of fixing one. PRIMARY_OMNIBUS names the row
# that carries the inferential claim at this level; the others are descriptive.
# Its role after DEG recovery is fixed in reorg plan v2 s3.1b: the pi0-free
# existence proof and the between-unit quiescence instrument, not a gatekeeper
# for 420.

source("setup.R")

suppressPackageStartupMessages({
  library(edgeR)
})

source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))
source(file.path(paths$root, "lib", "units.R"))

# --- Configuration ---------------------------------------------------------
# Shared constants (N_PERM, SEED, EXACT_THREADS, BM_EXACT_MAX, FDR_CUT) come
# from config.R via setup.R.
OMNIBUS_ALPHA <- c(1e-2, 1e-3, 1e-4) # per-gene null quantiles for the omnibus
HC_ALPHA0 <- 0.1 # fraction of the p-value range Higher Criticism scans
PRIMARY_OMNIBUS <- "hc" # the pre-specified inferential row

options(
  brunnermunzel.exact.max.allocations = BM_EXACT_MAX,
  brunnermunzel.exact.threads = EXACT_THREADS
)

pin_blas_threads()

# --- Load inputs -----------------------------------------------------------
norm_path <- file.path(paths$processed, "thyr_normalized_counts.rds")
if (!file.exists(norm_path)) {
  stop("thyr_normalized_counts.rds not found (run 310 first)")
}
normalized <- readRDS(norm_path)
message("Units: ", paste(names(normalized$units), collapse = ", "))

# --- Storey q (plug-in pi0, lambda = 0.5) ----------------------------------
# Identical to the DEGES screen in lib/norm_deges.R (reorg plan v2 D1):
# pi0_hat = min(1, mean(p > 0.5) / 0.5), q = min(1, pi0_hat * BH). Monotone
# because BH is monotone.
storey_q <- function(p) {
  pi0_hat <- min(1, mean(p > 0.5) / 0.5)
  list(pi0 = pi0_hat, q = pmin(1, pi0_hat * p.adjust(p, method = "BH")))
}

# The plug-in estimator applied to each shuffle, scored against its own
# gene-wise null: within each gene's row of null statistics, the empirical
# p of shuffle j is the fraction of shuffles at least as extreme. This keeps
# per-gene null heterogeneity (tie patterns differ across genes) out of the
# calibration, which pooling would mix in.
storey_pi0_null <- function(null_statistic) {
  n_perm <- ncol(null_statistic)
  null_p <- apply(null_statistic, 1L, function(s) {
    (n_perm - rank(s, ties.method = "min") + 1) / n_perm
  }) # n_perm x genes
  apply(null_p, 1L, function(p) min(1, mean(p > 0.5) / 0.5))
}

# --- Higher Criticism ------------------------------------------------------
# Scanning the threshold instead of fixing one, so no per-gene cut-off has to
# be chosen. The p-values are empirical against the pooled shuffle null, which
# both keeps real and shuffled genes on one scale and makes the statistic
# usable under inter-gene correlation; Higher Criticism's asymptotic null
# assumes independence and must not be used here.
higher_criticism <- function(statistic, pooled_null, n_null, alpha0) {
  n <- length(statistic)
  p <- sort(
    (n_null - findInterval(statistic, pooled_null, left.open = TRUE) + 1) /
      (n_null + 1)
  )
  index <- seq_len(n)
  usable <- index <= max(1L, floor(alpha0 * n)) & p > 1 / n
  if (!any(usable)) {
    return(NA_real_)
  }
  max((sqrt(n) * (index / n - p) / sqrt(p * (1 - p)))[usable])
}

# --- Omnibus test ----------------------------------------------------------
omnibus_table <- function(statistic, null_statistic, alpha, n_perm) {
  pooled <- sort(as.numeric(null_statistic))
  counted <- lapply(alpha, function(a) {
    cut <- pooled[max(1L, round(length(pooled) * (1 - a)))]
    null_count <- colSums(null_statistic >= cut)
    observed <- sum(statistic >= cut)
    data.frame(
      test = "count", alpha = a, cut = cut, observed = observed,
      null_median = stats::median(null_count),
      null_lo = unname(stats::quantile(null_count, 0.025)),
      null_hi = unname(stats::quantile(null_count, 0.975)),
      p = (sum(null_count >= observed) + 1) / (n_perm + 1),
      stringsAsFactors = FALSE
    )
  })
  null_max <- apply(null_statistic, 2L, max)
  observed_max <- max(statistic)
  counted[[length(counted) + 1L]] <- data.frame(
    test = "max", alpha = NA_real_, cut = NA_real_, observed = observed_max,
    null_median = stats::median(null_max),
    null_lo = unname(stats::quantile(null_max, 0.025)),
    null_hi = unname(stats::quantile(null_max, 0.975)),
    p = (sum(null_max >= observed_max) + 1) / (n_perm + 1),
    stringsAsFactors = FALSE
  )

  n_null <- length(pooled)
  observed_hc <- higher_criticism(statistic, pooled, n_null, HC_ALPHA0)
  null_hc <- apply(null_statistic, 2L, function(s) {
    higher_criticism(s, pooled, n_null, HC_ALPHA0)
  })
  counted[[length(counted) + 1L]] <- data.frame(
    test = "hc", alpha = HC_ALPHA0, cut = NA_real_, observed = observed_hc,
    null_median = stats::median(null_hc),
    null_lo = unname(stats::quantile(null_hc, 0.025)),
    null_hi = unname(stats::quantile(null_hc, 0.975)),
    p = (sum(null_hc >= observed_hc) + 1) / (n_perm + 1),
    stringsAsFactors = FALSE
  )
  do.call(rbind, counted)
}

# --- Per-unit test ---------------------------------------------------------
test_unit <- function(dgelist, unit) {
  arms <- unit_arms(dgelist$samples$group, unit)
  sporadic <- arms$sporadic
  high <- arms$high
  cpm_matrix <- edgeR::cpm(
    dgelist,
    normalized.lib.sizes = TRUE, prior.count = 0, log = FALSE
  )[, c(sporadic, high), drop = FALSE]
  storage.mode(cpm_matrix) <- "double"
  nx <- length(sporadic)
  n <- ncol(cpm_matrix)

  statistic <- abs(unname(brunnermunzel_statistics(cpm_matrix, nx)))
  effect <- unname(brunnermunzel_effects(cpm_matrix, nx))
  p_exact <- as.numeric(
    brunnermunzel_pvalues(cpm_matrix, nx, method = "exact")
  )

  # The label shuffles as an explicit index object, so 420 can consume the
  # same permutations by reference rather than by seed coincidence (s3.2).
  set.seed(SEED)
  perm_index <- vapply(seq_len(N_PERM), function(i) sample(n), integer(n))
  null_statistic <- vapply(
    seq_len(N_PERM),
    function(i) {
      abs(brunnermunzel_statistics(
        cpm_matrix[, perm_index[, i], drop = FALSE], nx
      ))
    },
    numeric(nrow(cpm_matrix))
  )

  storey <- storey_q(p_exact)
  pi0_null <- storey_pi0_null(null_statistic)

  genes <- data.frame(
    gene_id = rownames(cpm_matrix),
    effect = effect,
    statistic = statistic,
    p_exact = p_exact,
    q_storey = storey$q,
    rank = rank(-statistic, ties.method = "min"),
    stringsAsFactors = FALSE
  )
  genes <- genes[order(genes$rank, genes$gene_id), , drop = FALSE]
  rownames(genes) <- NULL

  q_thresholds <- seq(0.01, 0.30, by = 0.01)
  q_curve <- data.frame(
    threshold = q_thresholds,
    n_deg = vapply(
      q_thresholds, function(a) sum(storey$q < a), integer(1)
    )
  )

  omnibus <- omnibus_table(statistic, null_statistic, OMNIBUS_ALPHA, N_PERM)
  message(sprintf(
    "  %-9s %d vs %d ; pi0_hat %.3f (null med %.3f) ; q<%.2f %d ; %s p %.3f ; all %s",
    unit, nx, n - nx, storey$pi0, stats::median(pi0_null),
    FDR_CUT, sum(storey$q < FDR_CUT),
    PRIMARY_OMNIBUS, omnibus$p[omnibus$test == PRIMARY_OMNIBUS],
    paste(sprintf("%.3f", omnibus$p), collapse = "/")
  ))

  list(
    n_samples = c(Sporadic = nx, High = n - nx),
    genes = genes,
    pi0 = list(
      estimate = storey$pi0,
      lambda = 0.5,
      null = pi0_null,
      null_quantiles = stats::quantile(
        pi0_null, c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
      )
    ),
    q_curve = q_curve,
    omnibus = omnibus,
    null = list(
      n_perm = N_PERM,
      quantiles = stats::quantile(
        as.numeric(null_statistic),
        c(0.5, 0.9, 0.99, 0.999, 0.9999, 1)
      )
    ),
    perm_index = perm_index,
    perm_index_hash = digest::digest(perm_index)
  )
}

# --- Run every unit --------------------------------------------------------
message("Testing units (all = count 1e-2/1e-3/1e-4, max, hc):")
units <- lapply(names(normalized$units), function(unit) {
  test_unit(normalized$units[[unit]]$dgelist, unit)
})
names(units) <- names(normalized$units)

# --- Assemble and save -----------------------------------------------------
thyr_expression_test <- list(
  date = Sys.Date(),
  config = list(
    bm_method = "exact",
    alternative = "two.sided",
    n_perm = N_PERM,
    perm_seed = SEED,
    omnibus_alpha = OMNIBUS_ALPHA,
    hc_alpha0 = HC_ALPHA0,
    primary_omnibus = PRIMARY_OMNIBUS,
    exact_max = BM_EXACT_MAX,
    inference = "Storey q < 0.10 on exact p (plug-in pi0, lambda = 0.5)",
    pi0_lambda = 0.5,
    q_threshold = FDR_CUT,
    perm_sharing = "per-unit perm_index consumed by 420; hash per unit",
    perm_index_hash = vapply(
      units, function(u) u$perm_index_hash, character(1)
    ),
    reference_group = "Sporadic"
  ),
  units = units
)

out <- file.path(paths$processed, "thyr_expression_test.rds")
saveRDS(thyr_expression_test, out)
message("Saved: ", out, " (", length(units), " units)")
