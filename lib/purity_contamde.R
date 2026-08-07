# ==============================================================================
# MUREN-normalized, ContamDE-inspired relative tumor-purity estimation
#
# This file implements the initial matched-pair proportion score used by the
# original contamDE development code; it does not implement the iterative
# pseudo-likelihood estimator in Shen et al. (2016), Bioinformatics 32:705-712.
#
# Deliberate adaptations for the present analysis:
# - MUREN single-parameter MASS-LTS estimates the sample scaling coefficients.
# - Paired limma-voom with robust empirical Bayes selects working genes.
# - Up- and down-regulated genes are combined after aligning their signs, as in
#   the authors' development code.
# - Both the paper's mean-one normalized proportion and the development code's
#   maximum-one relative score are returned explicitly.
# ==============================================================================

# Note: Requires prior loading of:
# source("./lib/norm_muren_helpers.R")
# source("./lib/norm_muren.R")

# Required packages
for (pkg in c("edgeR", "limma", "statmod")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required.")
  }
}

.validate_paired_counts <- function(counts) {
  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }
  if (!is.matrix(counts) || !is.numeric(counts)) {
    stop("'counts' must be a numeric matrix or numeric data.frame.")
  }
  if (nrow(counts) < 2L) {
    stop("'counts' must contain at least two genes.")
  }
  if (ncol(counts) < 4L || ncol(counts) %% 2L != 0L) {
    stop("'counts' must contain at least two matched Normal/Tumor pairs.")
  }
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("'counts' must contain finite non-negative values.")
  }
  if (any(colSums(counts) <= 0)) {
    stop("Every sample must have a positive library size.")
  }

  n_pairs <- ncol(counts) %/% 2L
  sample_names <- colnames(counts)
  if (is.null(sample_names)) {
    stop(
      "Column names must identify ordered pairs as ",
      "'Normal_<pair>' followed by 'Tumor_<pair>'."
    )
  }

  normal_names <- sample_names[seq_len(n_pairs)]
  tumor_names <- sample_names[n_pairs + seq_len(n_pairs)]
  if (!all(grepl("^Normal_.+", normal_names)) ||
      !all(grepl("^Tumor_.+", tumor_names))) {
    stop(
      "The first half of columns must be named 'Normal_<pair>' and ",
      "the second half 'Tumor_<pair>'."
    )
  }

  normal_pairs <- sub("^Normal_", "", normal_names)
  tumor_pairs <- sub("^Tumor_", "", tumor_names)
  if (!identical(normal_pairs, tumor_pairs) || anyDuplicated(normal_pairs)) {
    stop("Normal/Tumor column names do not define unique matched pairs.")
  }

  list(counts = counts, n_pairs = n_pairs, pair_ids = normal_pairs)
}

# ==============================================================================
# limma + voom with MUREN scaling factors
# ==============================================================================
limma_voom_purity <- function(counts, workers) {
  if (missing(workers)) {
    stop("'workers' must be explicitly specified.")
  }

  validated <- .validate_paired_counts(counts)
  counts <- validated$counts
  n_pairs <- validated$n_pairs

  if (!is.numeric(workers) || length(workers) != 1L ||
      !is.finite(workers) || workers <= 0 || workers != floor(workers)) {
    stop("'workers' must be a fixed positive integer.")
  }
  workers <- as.integer(workers)

  if (!exists("muren_norm", mode = "function", inherits = TRUE)) {
    stop("muren_norm() must be loaded before limma_voom_purity().")
  }

  d <- edgeR::DGEList(counts = counts)

  # muren_norm(..., scaling_coeff) returns the sample coefficient by which raw
  # counts must be divided. It is an effective sample size, not an edgeR
  # norm.factor.
  muren_size <- muren_norm(
    counts,
    refs = "saturated",
    pairwise_method = "lts",
    single_param = TRUE,
    res_return = "scaling_coeff",
    workers = workers
  )
  muren_size <- as.numeric(muren_size)
  names(muren_size) <- colnames(counts)
  if (length(muren_size) != ncol(counts) ||
      any(!is.finite(muren_size)) || any(muren_size <= 0)) {
    stop("MUREN scaling coefficients contain non-finite or non-positive values.")
  }

  # edgeR defines effective library size as lib.size * norm.factors. Convert the
  # MUREN coefficient to that representation, then center only the edgeR
  # norm.factor. The resulting effective sizes remain proportional to MUREN.
  library_size <- as.numeric(d$samples$lib.size)
  edge_norm_factors <- muren_size / library_size
  edge_norm_factors <- edge_norm_factors /
    exp(mean(log(edge_norm_factors)))
  d$samples$norm.factors <- unname(edge_norm_factors)
  effective_library_size <- library_size * edge_norm_factors

  effective_ratio <- effective_library_size / muren_size
  if (max(abs(effective_ratio / effective_ratio[[1L]] - 1)) > 1e-10) {
    stop("edgeR effective library sizes are not proportional to MUREN coefficients.")
  }

  # voom transformation (normalize.method fixed to "none")
  # The patient-pair fixed effect corresponds to the shared alpha_ij in the
  # matched-sample contamDE model.
  condition <- factor(c(rep("Normal", n_pairs), rep("Tumor", n_pairs)),
    levels = c("Normal", "Tumor")
  )
  pair <- factor(rep(validated$pair_ids, 2))
  design <- stats::model.matrix(~ 0 + condition + pair)
  v <- limma::voom(d, design, normalize.method = "none")

  fit <- limma::lmFit(v, design)
  contrs <- limma::makeContrasts(
    contrasts = "conditionTumor - conditionNormal",
    levels = design
  )
  fit2 <- limma::contrasts.fit(fit, contrs)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

  return(list(
    counts = counts,
    size = muren_size,
    muren_scaling_coeff = muren_size,
    edge_norm_factors = edge_norm_factors,
    effective_library_size = effective_library_size,
    design = design,
    p.limma = fit2$p.value[, 1],
    log2FC.limma = fit2$coefficients[, 1],
    fit = fit2
  ))
}

# ==============================================================================
# Main purity estimation function (lightweight contamDE)
# ==============================================================================
contamde_purity <- function(counts,
                            subtype = NULL,
                            covariate = NULL,
                            contaminated = TRUE,
                            pairwise_method = "lts",
                            workers,
                            verbose = TRUE) {
  if (missing(workers)) {
    stop("'workers' must be explicitly specified.")
  }

  # These arguments are retained temporarily for existing call sites. This
  # purity-only estimator has one fixed model and does not implement subtype,
  # covariate, uncontaminated, or alternative MUREN-method branches.
  if (!is.null(subtype) || !is.null(covariate)) {
    stop("'subtype' and 'covariate' are not part of this purity-only estimator.")
  }
  if (!isTRUE(contaminated)) {
    stop("'contaminated' must be TRUE for purity estimation.")
  }
  if (!identical(pairwise_method, "lts")) {
    stop("This estimator uses the fixed MUREN pairwise method 'lts'.")
  }

  validated <- .validate_paired_counts(counts)
  counts <- validated$counts
  n_pairs <- validated$n_pairs
  n_genes <- nrow(counts)

  if (verbose) cat("Starting tumor purity estimation...\n")

  if (verbose) {
    cat(sprintf("  Input: %d genes, %d matched pairs\n", n_genes, n_pairs))
  }

  if (verbose) cat("  Running limma+voom with MUREN normalization...\n")

  d <- limma_voom_purity(
    counts = counts,
    workers = workers
  )

  p_limma <- d$p.limma
  log2_fc_limma <- d$log2FC.limma
  muren_size <- d$muren_scaling_coeff

  if (length(p_limma) != n_genes || length(log2_fc_limma) != n_genes ||
      any(!is.finite(p_limma)) || any(!is.finite(log2_fc_limma))) {
    stop("limma results are incomplete or non-finite.")
  }

  # Use the MUREN coefficient itself as the count divisor. Re-centering this
  # coefficient would change the +1 pseudocount scale used by contamDE.
  count_norm <- sweep(counts, 2, muren_size, `/`)
  idx_n <- seq_len(n_pairs)
  idx_t <- n_pairs + seq_len(n_pairs)

  # Pseudocount +1 follows the authors' development code.
  y <- log2((count_norm[, idx_t, drop = FALSE] + 1) /
    (count_norm[, idx_n, drop = FALSE] + 1))
  colnames(y) <- validated$pair_ids

  if (verbose) cat("  Estimating relative tumor proportions...\n")

  p_adj <- stats::p.adjust(p_limma, method = "BH")

  # Retain at most the 1000 smallest adjusted P-values among candidates, as in
  # the authors' development code.
  n_sig <- sum(p_adj < 0.1)
  if (n_sig > 1000L) {
    idx_sig <- which(p_adj < 0.1)
    keep <- idx_sig[order(p_adj[idx_sig])][seq_len(1000L)]
    p_adj[-keep] <- 1
    if (verbose) {
      cat(sprintf(
        "  Info: Limited working genes to 1000 (from %d)\n",
        n_sig
      ))
    }
  }

  up <- which(p_adj < 0.1 & log2_fc_limma > log2(1.5))
  down <- which(p_adj < 0.1 & log2_fc_limma < -log2(1.5))

  if (verbose) {
    cat(sprintf(
      "  Working genes - Up: %d, Down: %d\n",
      length(up), length(down)
    ))
  }

  y_up <- if (length(up) > 0L) {
    y[up, , drop = FALSE]
  } else {
    matrix(0, nrow = 0L, ncol = n_pairs)
  }
  y_down <- if (length(down) > 0L) {
    y[down, , drop = FALSE]
  } else {
    matrix(0, nrow = 0L, ncol = n_pairs)
  }

  sumup <- colSums(y_up)
  sumdown <- colSums(y_down)
  raw_score <- sumup - sumdown
  names(sumup) <- names(sumdown) <- names(raw_score) <- validated$pair_ids

  if (length(up) + length(down) == 0L ||
      any(!is.finite(raw_score)) || any(raw_score <= 0)) {
    stop(
      "The working-gene score is not strictly positive for every pair; ",
      "the relative-purity estimator is not established."
    )
  }

  normalized_proportion <- raw_score / mean(raw_score)
  relative_score <- raw_score / max(raw_score)

  if (verbose) {
    cat(sprintf(
      "  Relative score summary: mean=%.3f, sd=%.3f\n",
      mean(relative_score), stats::sd(relative_score)
    ))
    cat("  Relative tumor-proportion estimation completed.\n")
  }

  return(list(
    # Backward-compatible maximum-one score used by the current filtering code.
    proportion = relative_score,
    relative_score = relative_score,
    # Mean-one identifiability constraint used in the contamDE paper.
    normalized_proportion = normalized_proportion,
    raw_score = raw_score,
    score_components = list(sumup = sumup, sumdown = sumdown),
    counts = counts,
    size = muren_size,
    muren_scaling_coeff = muren_size,
    edge_norm_factors = d$edge_norm_factors,
    effective_library_size = d$effective_library_size,
    design = d$design,
    y = y,
    p_limma = p_limma,
    p_adjusted = p_adj,
    log2FC_limma = log2_fc_limma,
    n_pairs = n_pairs,
    n_genes = n_genes,
    informative_genes = list(up = up, down = down),
    informative_gene_names = list(
      up = if (is.null(rownames(counts))) NULL else rownames(counts)[up],
      down = if (is.null(rownames(counts))) NULL else rownames(counts)[down]
    ),
    normalization_method = "MUREN_LTS",
    working_gene_method = "paired_voom_robust_ebayes_signed_bidirectional",
    estimation_date = Sys.time()
  ))
}

# ==============================================================================
# Quality assessment function
# ==============================================================================
assess_purity_quality <- function(purity_result, threshold = 0.6, verbose = TRUE) {
  if (!"proportion" %in% names(purity_result)) {
    stop("Input must be result from contamde_purity function.")
  }

  # 'proportion' is the maximum-one relative score retained for compatibility;
  # the threshold is not an absolute histological tumor-purity threshold.
  proportions <- purity_result$proportion
  n_samples <- length(proportions)

  # Calculate statistics
  purity_stats <- list(
    n_samples = n_samples,
    mean_purity = mean(proportions),
    median_purity = median(proportions),
    sd_purity = sd(proportions),
    min_purity = min(proportions),
    max_purity = max(proportions),
    threshold = threshold
  )

  # Quality assessment
  high_purity <- proportions >= threshold
  n_high_purity <- sum(high_purity)
  retention_rate <- n_high_purity / n_samples

  purity_stats$n_high_purity <- n_high_purity
  purity_stats$retention_rate <- retention_rate
  purity_stats$high_purity_indices <- which(high_purity)

  if (verbose) {
    cat("\n=== Purity Quality Assessment ===\n")
    cat(sprintf("Total samples: %d\n", n_samples))
    cat(sprintf(
      "Mean purity: %.3f +/- %.3f\n",
      purity_stats$mean_purity, purity_stats$sd_purity
    ))
    cat(sprintf(
      "Range: [%.3f, %.3f]\n",
      purity_stats$min_purity, purity_stats$max_purity
    ))
    cat(sprintf(
      "High purity (>=%.1f): %d/%d (%.1f%%)\n",
      threshold, n_high_purity, n_samples, retention_rate * 100
    ))

    if (retention_rate >= 0.7) {
      cat("Retention: Excellent (>=70%)\n")
    } else if (retention_rate >= 0.5) {
      cat("Retention: Acceptable (50-70%)\n")
    } else {
      cat("Retention: Poor (<50%)\n")
    }
  }

  return(purity_stats)
}

# ==============================================================================
# Helper function for filtered sample list creation
# ==============================================================================
create_purity_filtered_lists <- function(original_sample_lists, purity_results,
                                         threshold = 0.6, verbose = TRUE) {
  if (verbose) cat("\nCreating purity-filtered sample lists...\n")

  filtered_lists <- list()

  for (group_name in names(original_sample_lists)) {
    if (!group_name %in% names(purity_results) ||
      is.null(purity_results[[group_name]])) {
      if (verbose) {
        cat(sprintf("  %s: No purity results available\n", group_name))
      }
      filtered_lists[[group_name]] <- list(
        tumor = character(0),
        normal = character(0),
        cases = character(0)
      )
      next
    }

    # Get purity assessment
    purity_result <- purity_results[[group_name]]
    quality_stats <- assess_purity_quality(purity_result, threshold, verbose = FALSE)

    # Get high purity indices
    high_purity_indices <- quality_stats$high_purity_indices

    if (length(high_purity_indices) == 0) {
      if (verbose) {
        cat(sprintf("  %s: No high purity samples\n", group_name))
      }
      filtered_lists[[group_name]] <- list(
        tumor = character(0),
        normal = character(0),
        cases = character(0)
      )
      next
    }

    # Filter samples based on purity
    original_group <- original_sample_lists[[group_name]]
    n_original <- length(original_group$tumor)

    filtered_tumor <- original_group$tumor[high_purity_indices]
    filtered_normal <- original_group$normal[high_purity_indices]
    filtered_cases <- original_group$cases[high_purity_indices]

    filtered_lists[[group_name]] <- list(
      tumor = filtered_tumor,
      normal = filtered_normal,
      cases = filtered_cases
    )

    if (verbose) {
      cat(sprintf(
        "  %s: %d -> %d samples (%.1f%% retention)\n",
        group_name, n_original, length(filtered_tumor),
        length(filtered_tumor) / n_original * 100
      ))
    }
  }

  if (verbose) cat("Purity filtering completed!\n")
  return(filtered_lists)
}

cat("MUREN-normalized ContamDE-inspired purity functions loaded successfully!\n")
cat("Key features:\n")
cat("  - Paired voom with MUREN effective library sizes\n")
cat("  - Robust eBayes and signed bidirectional working genes\n")
cat("  - Mean-one normalized proportions and maximum-one relative scores\n")
cat("Main functions:\n")
cat("  - contamde_purity(): Estimate relative tumor proportions\n")
cat("  - assess_purity_quality(): Relative-score assessment\n")
cat("  - create_purity_filtered_lists(): Sample filtering\n")
