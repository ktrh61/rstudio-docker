# 080_test_expression_difference.R
# Test each analysis unit for a High-vs-Sporadic expression difference on the
# DEGES-normalized counts from 070, with the exact permutation Brunner-Munzel
# test. No differential-expression call is made: the products are a complete
# ranked gene table and a unit-level omnibus permutation test, which is what
# the downstream gene-set analysis consumes.
# Input : processed/thyr_normalized_counts.rds  (from 070; per-unit DGEList)
#         utils/brunnermunzel_mc.R
# Output: processed/thyr_expression_test.rds
#           list(date, config, units)
#           units = { R_Tumor, R_Normal, B_Tumor, B_Normal }; each holds
#             n_samples : Sporadic and High column counts
#             genes     : gene_id, effect, statistic, p_exact, fdr_perm, rank
#             omnibus   : test, alpha, cut, observed, null_median, null_lo,
#                         null_hi, p
#             null      : n_perm and quantiles of the pooled null statistic
#
# Group x is Sporadic and group y is High, so effect > 0.5 means the gene runs
# higher in the exposed arm. p_exact enumerates every C(n, nx) allocation and
# needs no seed. fdr_perm is a permutation-calibrated FDR on the statistic:
# the null comes from N_PERM whole-column shuffles, which preserve inter-gene
# correlation, and pi0 is held at 1. The omnibus "count" rows ask whether the
# number of genes past a null quantile exceeds what those same shuffles
# produce; the "max" row asks the same of the single most extreme gene.

source("setup.R")

suppressPackageStartupMessages({
  library(edgeR)
})

source(file.path(paths$root, "utils", "brunnermunzel_mc.R"))

# --- Configuration ---------------------------------------------------------
N_PERM <- 999L # label shuffles for the empirical null
PERM_SEED <- 19860426L # fixed seed (Chernobyl accident date 1986-04-26)
OMNIBUS_ALPHA <- c(1e-2, 1e-3, 1e-4) # per-gene null quantiles for the omnibus
EXACT_MAX <- 1e8 # largest C(n, nx) still enumerated exactly
EXACT_THREADS <- 16L # threads for the enumeration

options(
  brunnermunzel.exact.max.allocations = EXACT_MAX,
  brunnermunzel.exact.threads = EXACT_THREADS
)

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}

# --- Load inputs -----------------------------------------------------------
norm_path <- file.path(paths$processed, "thyr_normalized_counts.rds")
if (!file.exists(norm_path)) {
  stop("thyr_normalized_counts.rds not found (run 070 first)")
}
normalized <- readRDS(norm_path)
message("Units: ", paste(names(normalized$units), collapse = ", "))

# --- Permutation-calibrated FDR --------------------------------------------
# Step-up on the statistic itself rather than on p_exact: the two are
# equivalent within a tie pattern, and the shuffled null already carries the
# mixture of tie patterns that the gene set actually has.
permutation_fdr <- function(statistic, null_statistic, n_perm) {
  ordering <- order(statistic, decreasing = TRUE)
  sorted_null <- sort(null_statistic)
  expected <- (length(sorted_null) -
    findInterval(statistic[ordering], sorted_null, left.open = TRUE)) / n_perm
  fdr <- pmin(1, expected / seq_along(ordering))
  out <- numeric(length(statistic))
  out[ordering] <- rev(cummin(rev(fdr)))
  out
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
  do.call(rbind, counted)
}

# --- Per-unit test ---------------------------------------------------------
test_unit <- function(dgelist, unit) {
  group <- as.character(dgelist$samples$group)
  sporadic <- which(grepl("Sporadic", group, fixed = TRUE))
  high <- which(grepl("High", group, fixed = TRUE))
  if (length(sporadic) + length(high) != ncol(dgelist)) {
    stop("Unit ", unit, " has samples outside Sporadic / High.")
  }
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

  set.seed(PERM_SEED)
  null_statistic <- vapply(
    seq_len(N_PERM),
    function(i) {
      abs(brunnermunzel_statistics(
        cpm_matrix[, sample(n), drop = FALSE], nx
      ))
    },
    numeric(nrow(cpm_matrix))
  )

  genes <- data.frame(
    gene_id = rownames(cpm_matrix),
    effect = effect,
    statistic = statistic,
    p_exact = p_exact,
    fdr_perm = permutation_fdr(statistic, null_statistic, N_PERM),
    rank = rank(-statistic, ties.method = "min"),
    stringsAsFactors = FALSE
  )
  genes <- genes[order(genes$rank, genes$gene_id), , drop = FALSE]
  rownames(genes) <- NULL

  omnibus <- omnibus_table(statistic, null_statistic, OMNIBUS_ALPHA, N_PERM)
  message(sprintf(
    "  %-9s %d vs %d ; min p_exact %.3e ; fdr_perm<0.10 %d ; omnibus p %s",
    unit, nx, n - nx, min(p_exact), sum(genes$fdr_perm < 0.10),
    paste(sprintf("%.3f", omnibus$p), collapse = "/")
  ))

  list(
    n_samples = c(Sporadic = nx, High = n - nx),
    genes = genes,
    omnibus = omnibus,
    null = list(
      n_perm = N_PERM,
      quantiles = stats::quantile(
        as.numeric(null_statistic),
        c(0.5, 0.9, 0.99, 0.999, 0.9999, 1)
      )
    )
  )
}

# --- Run every unit --------------------------------------------------------
message("Testing units (omnibus p reported as 1e-2/1e-3/1e-4/max):")
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
    perm_seed = PERM_SEED,
    omnibus_alpha = OMNIBUS_ALPHA,
    exact_max = EXACT_MAX,
    fdr = "permutation-calibrated, pi0 = 1",
    reference_group = "Sporadic"
  ),
  units = units
)

out <- file.path(paths$processed, "thyr_expression_test.rds")
saveRDS(thyr_expression_test, out)
message("Saved: ", out, " (", length(units), " units)")
