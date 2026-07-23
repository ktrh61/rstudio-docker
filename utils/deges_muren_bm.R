# deges_muren_bm.R
# TCC-faithful DEGES normalization (the X-Y-X / iDEGES protocol of TCC; Kadota
# et al.) with MUREN as the normalizer and a Monte Carlo studentized
# Brunner-Munzel permutation test as the DEG screen. MUREN replaces TMM and the
# Monte Carlo Brunner-Munzel test replaces edgeR's exact test; the iteration,
# floorPDEG screen, and degenerate guard follow TCC::calcNormFactors.
#
# Requires (sourced by the caller before use):
#   utils/norm_improved.R      (muren_norm)
#   utils/brunnermunzel_mc.R   (brunnermunzel_mc_test)
#   package edgeR
#
# deges_muren_bm(counts, group, iteration, fdr, floor_pdeg, ...) runs, on a
# fixed gene x sample count matrix:
#   STEP 1 : initial MUREN scaling coefficients on the full matrix.
#   STEP 2 : repeat `iteration` times -- screen DEGs on the FULL matrix with the
#            current coefficients (MC Brunner-Munzel -> BH; normal cutoff set vs
#            a floorPDEG raw-p rank set, the larger adopted), then re-estimate
#            MUREN coefficients from the non-DEG genes. Runs exactly `iteration`
#            times; stops early only when no non-DEG genes remain (TCC guard).
# Returns the final per-sample scaling coefficients and per-iteration
# diagnostics. Convert the coefficients to edgeR norm.factors for a target
# DGEList with muren_to_norm_factors().

# --- MUREN scaling coefficient -> edgeR norm.factors -----------------------
# scaling_coeff is the per-sample size by which raw counts must be divided.
# Divide by the target DGEList's lib.size and centre geometrically, so
# lib.size * norm.factors stays proportional to scaling_coeff. Follows
# utils/contamde_purity_functions.R. scaling_coeff and lib_size share order.
muren_to_norm_factors <- function(scaling_coeff, lib_size) {
  if (any(!is.finite(scaling_coeff)) || any(scaling_coeff <= 0)) {
    stop("MUREN scaling coefficients must be finite and positive.")
  }
  nf <- scaling_coeff / lib_size
  nf <- nf / exp(mean(log(nf)))
  effective <- lib_size * nf
  ratio <- effective / scaling_coeff
  if (max(abs(ratio / ratio[1] - 1)) > 1e-10) {
    stop("effective library sizes not proportional to MUREN coefficients.")
  }
  unname(nf)
}

# --- MUREN scaling coefficients for a count matrix -------------------------
.deges_muren_scaling <- function(counts, muren_method, workers) {
  sc <- muren_norm(
    counts,
    refs = "saturated",
    pairwise_method = muren_method,
    single_param = TRUE,
    res_return = "scaling_coeff",
    workers = workers
  )
  if (any(!is.finite(sc)) || any(sc <= 0)) {
    stop("MUREN scaling coefficients must be finite and positive.")
  }
  sc
}

# --- Per-gene Monte Carlo Brunner-Munzel p-values --------------------------
# The same n_perm/seed/alternative/est go to every gene; a fixed seed
# makes the result reproducible and enables the null-distribution cache.
.deges_bm_pvalues <- function(cpm_matrix, group, group_levels,
                              n_perm, seed, alternative, est) {
  g1 <- which(group == group_levels[1])
  g2 <- which(group == group_levels[2])
  n_genes <- nrow(cpm_matrix)
  pvalues <- numeric(n_genes)
  for (i in seq_len(n_genes)) {
    pvalues[i] <- brunnermunzel_mc_test(
      cpm_matrix[i, g1], cpm_matrix[i, g2],
      alternative = alternative, est = est, B = n_perm, seed = seed
    )$p.value
  }
  pvalues
}

# --- DEGES normalization (MUREN + MC Brunner-Munzel) -----------------------
deges_muren_bm <- function(counts, group, iteration = 1L,
                           fdr, floor_pdeg,
                           n_perm, seed, alternative = "two.sided",
                           est = "original",
                           muren_method = "lts", workers = 3L) {
  if (!exists("muren_norm", mode = "function", inherits = TRUE)) {
    stop("muren_norm() must be loaded (source utils/norm_improved.R).")
  }
  if (!exists("brunnermunzel_mc_test", mode = "function", inherits = TRUE)) {
    stop("brunnermunzel_mc_test() must be loaded (utils/brunnermunzel_mc.R).")
  }
  group <- factor(group)
  group_levels <- levels(group)
  if (length(group_levels) != 2L) {
    stop("group must have exactly two levels.")
  }
  sample_ids <- colnames(counts)
  n_genes <- nrow(counts)

  # STEP 1: initial MUREN scaling on the full matrix.
  scaling_coeff <- .deges_muren_scaling(counts, muren_method, workers)
  full_dge <- edgeR::DGEList(counts = counts, group = group)

  n_non_deg <- integer(iteration)
  deg_count <- integer(iteration)
  exclusion_method <- character(iteration)
  n_done <- 0L

  for (i in seq_len(iteration)) {
    # DE screen on the FULL matrix with the current scaling coefficients.
    nf <- muren_to_norm_factors(
      scaling_coeff[sample_ids], full_dge$samples$lib.size
    )
    full_dge$samples$norm.factors <- nf
    cpm_matrix <- edgeR::cpm(
      full_dge, normalized.lib.sizes = TRUE, prior.count = 0, log = FALSE
    )
    pvalues <- .deges_bm_pvalues(
      cpm_matrix, group, group_levels, n_perm, seed, alternative, est
    )
    p_adj <- p.adjust(pvalues, method = "BH")

    # Normal cutoff set vs floorPDEG raw-p rank set; adopt the larger.
    normal_deg <- which(p_adj < fdr)
    floor_deg <- which(
      rank(pvalues, ties.method = "min") <= n_genes * floor_pdeg
    )
    if (length(normal_deg) < length(floor_deg)) {
      potential_deg <- floor_deg
      method_used <- "floor"
    } else {
      potential_deg <- normal_deg
      method_used <- "p_adj_cutoff"
    }
    non_deg <- setdiff(seq_len(n_genes), potential_deg)

    n_done <- i
    n_non_deg[i] <- length(non_deg)
    deg_count[i] <- length(normal_deg)
    exclusion_method[i] <- method_used
    message(sprintf(
      "      iteration %d: %d DEGs (%s), %d non-DEG genes retained",
      i, length(normal_deg), method_used, length(non_deg)
    ))

    # TCC degenerate guard: stop if removing DEGs leaves no non-DEG genes.
    if (length(non_deg) == 0L) {
      message("      no non-DEG genes remain; stopping DEGES iterations.")
      break
    }

    # MUREN re-estimation from the non-DEG genes.
    scaling_coeff <- .deges_muren_scaling(
      counts[non_deg, , drop = FALSE], muren_method, workers
    )
  }

  keep <- seq_len(n_done)
  list(
    scaling_coeff = scaling_coeff,
    iterations = list(
      n = n_done,
      n_non_deg = n_non_deg[keep],
      deg_count = deg_count[keep],
      exclusion_method = exclusion_method[keep]
    )
  )
}
