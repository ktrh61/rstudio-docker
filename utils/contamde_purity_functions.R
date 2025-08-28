# ==============================================================================
# ContamDE Purity Estimation Functions v6 (Lightweight Version)
# contamde_purity_functions_v6.R
# ==============================================================================
# 
# Purpose: Tumor purity estimation only (no DEG analysis)
# Modifications for v6:
# - Updated MUREN integration for compatibility with utils_improved.R
# - LTS option support for better stability
# - Simplified verbose output

# Note: Requires prior loading of:
# source("./utils/utils_improved.R")
# source("./utils/norm_improved.R")

# ==============================================================================
# limma + voom with MUREN scaling factors
# ==============================================================================
limma_voom_purity <- function(counts, pairwise_method = "lts",
                              workers = "auto", voom_norm = "quantile") {
  
  if (!is.matrix(counts)) {
    counts <- as.matrix(counts)
  }
  
  ncol_counts <- ncol(counts) / 2
  stopifnot(ncol_counts == floor(ncol(counts) / 2))
  
  # Create DGEList
  d <- edgeR::DGEList(counts = counts)
  
  # MUREN scaling coefficients using improved version
  # pairwise_method: "lts", "mode", "median", "trim10", "huber"
  c <- muren_norm(
    counts,
    refs = "saturated",
    pairwise_method = pairwise_method,
    single_param = TRUE,
    res_return = "scaling_coeff",
    workers = workers
  )
  
  # Validate coefficients
  if (any(!is.finite(c)) || any(c <= 0)) {
    stop("MUREN scaling coefficients contain non-finite or non-positive values.")
  }
  
  # Rescale to geometric mean = 1
  geo_mean <- exp(mean(log(c)))
  c <- c / geo_mean
  
  # Apply to edgeR
  d$samples$norm.factors <- as.numeric(c)
  
  # voom transformation
  condition <- factor(c(rep("Normal", ncol_counts), rep("Tumor", ncol_counts)))
  design <- stats::model.matrix(~ 0 + condition)
  colnames(design) <- levels(condition)
  v <- limma::voom(d, design, normalize.method = voom_norm)
  
  # Effective library sizes
  size <- d$samples$lib.size * d$samples$norm.factors
  size <- size / mean(size)
  
  # limma fit for p-values and log2FC
  fit <- limma::lmFit(v, design)
  contrs <- limma::makeContrasts(contrasts = "Tumor - Normal",
                                 levels = design)
  fit2 <- limma::contrasts.fit(fit, contrs)
  fit2 <- limma::eBayes(fit2, trend = TRUE)
  
  return(list(
    counts = counts,
    size = size,
    p.limma = fit2$p.value[, 1],
    log2FC.limma = fit2$coefficients[, 1]
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
                            workers = "auto",
                            voom_norm = "quantile",
                            prior.count = 0,
                            verbose = TRUE) {
  
  if (verbose) cat("Starting tumor purity estimation...\n")
  
  counts <- as.matrix(counts)
  ncol_counts <- ncol(counts) / 2
  stopifnot(ncol_counts == floor(ncol(counts) / 2))
  nrow_counts <- nrow(counts)
  
  if (verbose) {
    cat(sprintf("  Input: %d genes, %d paired samples\n", nrow_counts, ncol_counts))
  }
  
  # limma+voom with MUREN normalization
  if (verbose) cat("  Running limma+voom with MUREN normalization...\n")
  
  d <- limma_voom_purity(
    counts = counts, 
    pairwise_method = pairwise_method,
    workers = workers, 
    voom_norm = voom_norm
  )
  
  p_limma <- d$p.limma
  log2_fc_limma <- d$log2FC.limma
  size <- d$size
  
  # Normalize counts by MUREN size factors
  count_norm <- t(t(counts) / size)
  
  # Handle zero counts
  if (prior.count > 0) {
    count_norm_valid <- count_norm
    valid_genes <- rep(TRUE, nrow(count_norm))
  } else {
    valid_genes <- apply(count_norm, 1, function(x) all(x > 0))
    count_norm_valid <- count_norm[valid_genes, , drop = FALSE]
  }
  
  counts_valid <- counts[valid_genes, , drop = FALSE]
  
  if (verbose) {
    n_exc <- sum(!valid_genes)
    cat(sprintf("  Excluded %d genes with zero counts (%.1f%%)\n", 
                n_exc, 100 * n_exc / nrow(count_norm)))
  }
  
  # Split normal/tumor and calculate log2 ratios
  idx_n <- seq_len(ncol_counts)
  idx_t <- setdiff(seq_len(ncol(count_norm_valid)), idx_n)
  
  if (prior.count > 0) {
    y <- log2((count_norm_valid[, idx_t, drop=FALSE] + prior.count) /
                (count_norm_valid[, idx_n, drop=FALSE] + prior.count))
  } else {
    y <- log2(count_norm_valid[, idx_t, drop=FALSE]) -
      log2(count_norm_valid[, idx_n, drop=FALSE])
  }
  
  # Design matrix (simplified)
  if (is.null(subtype) || length(unique(subtype)) == 1) {
    subtype <- rep(1, ncol_counts)
    if (is.null(covariate)) {
      design <- matrix(subtype, ncol = 1)
    } else {
      covariate <- matrix(covariate, ncol = ncol_counts)
      design <- stats::model.matrix(~ covariate)
    }
  } else {
    subtype <- factor(subtype)
    if (is.null(covariate)) {
      design <- stats::model.matrix(~ 0 + subtype)
    } else {
      covariate0 <- matrix(covariate, nrow = ncol_counts)
      design <- stats::model.matrix(~ 0 + subtype + covariate0)
      rm(covariate0)
    }
  }
  
  # Initialize informative gene lists
  up <- integer(0)
  down <- integer(0)
  
  # Estimate purity
  w_hat <- rep(1, ncol_counts)
  
  if (contaminated) {
    if (verbose) cat("  Estimating tumor purity proportions...\n")
    
    p_limma <- p_limma[valid_genes]
    log2_fc_limma <- log2_fc_limma[valid_genes]
    
    # Calculate adjusted p-values
    p_adj <- tryCatch(
      qvalue::qvalue(p_limma, pi0.method = "bootstrap")$qvalues,
      error = function(e) {
        stop("qvalue calculation failed: ", e$message)
      }
    )
    
    # Define informative genes
    log2_fc <- log2_fc_limma
    if (sum(p_adj < 0.1) > 1e3) {
      p_adj[-order(p_adj)[1:1e3]] <- 1
    }
    
    up <- which(p_adj < 0.1 & log2_fc > log2(1.5))
    down <- which(p_adj < 0.1 & log2_fc < -log2(1.5))
    
    if (verbose) {
      cat(sprintf("  Informative genes - Up: %d, Down: %d\n", 
                  length(up), length(down)))
    }
    
    # Calculate purity estimates
    if (length(up) > 0) {
      y_up <- y[up, , drop = FALSE]
    } else {
      y_up <- matrix(0, 0, ncol_counts)
    }
    
    if (length(down) > 0) {
      y_down <- y[down, , drop = FALSE]
    } else {
      y_down <- matrix(0, 0, ncol_counts)
    }
    
    sumup <- colSums(y_up)
    sumdown <- colSums(y_down)
    sum.max <- max(sumup - sumdown)
    
    if (!is.finite(sum.max) || sum.max <= 0) {
      warning("No informative genes found. Setting purity to 1.0.")
      w_hat <- rep(1.0, ncol_counts)
    } else {
      w_hat <- (sumup - sumdown) / sum.max
      w_hat <- pmax(0, pmin(1, w_hat))
    }
    
    if (verbose) {
      cat(sprintf("  Purity summary: mean=%.3f, sd=%.3f\n", 
                  mean(w_hat), sd(w_hat)))
    }
  } else {
    if (verbose) {
      cat("  Contamination correction disabled. Setting purity to 1.0.\n")
    }
  }
  
  # Validate results
  stopifnot(length(w_hat) == ncol_counts, all(is.finite(w_hat)))
  
  if (verbose) {
    cat("  Purity estimation completed successfully!\n")
  }
  
  return(list(
    proportion = w_hat,
    counts = counts_valid,
    size = size,
    design = design,
    y = y,
    n_pairs = ncol_counts,
    n_genes = nrow(counts_valid),
    informative_genes = list(up = up, down = down),
    normalization_method = "MUREN_LTS",
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
    cat(sprintf("Mean purity: %.3f ± %.3f\n", 
                purity_stats$mean_purity, purity_stats$sd_purity))
    cat(sprintf("Range: [%.3f, %.3f]\n", 
                purity_stats$min_purity, purity_stats$max_purity))
    cat(sprintf("High purity (≥%.1f): %d/%d (%.1f%%)\n",
                threshold, n_high_purity, n_samples, retention_rate * 100))
    
    if (retention_rate >= 0.7) {
      cat("Retention: Excellent (≥70%)\n")
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
      cat(sprintf("  %s: %d → %d samples (%.1f%% retention)\n",
                  group_name, n_original, length(filtered_tumor),
                  length(filtered_tumor) / n_original * 100))
    }
  }
  
  if (verbose) cat("Purity filtering completed!\n")
  return(filtered_lists)
}

cat("ContamDE purity functions v6 loaded successfully!\n")
cat("Using MUREN normalization with LTS option for stability\n")
cat("Main functions:\n")
cat("  - contamde_purity(): Estimate tumor purity\n")
cat("  - assess_purity_quality(): Quality assessment\n")
cat("  - create_purity_filtered_lists(): Sample filtering\n")