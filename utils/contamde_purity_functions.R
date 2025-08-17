# ============================================================================== 
# ContamDE Purity Estimation Functions (Lightweight Version)
# contamde_purity_functions.R
# ==============================================================================
# 
# Purpose: Tumor purity estimation only (no DEG analysis)
# Based on: contamde_functions_v2.R
# Modifications: Removed DEG analysis, focused on proportion estimation

# ---- limma + voom (quantile維持) with MUREN scaling factors -----------------
limma_voom_purity <- function(counts) {
  if (!is.matrix(counts)) counts <- as.matrix(counts)
  ncol_counts <- ncol(counts) / 2
  
  # edgeR: DGE + TMM（lib.size算出のため）。norm.factorsは後でMURENに置換
  d <- edgeR::DGEList(counts = counts)
  d <- edgeR::calcNormFactors(d)
  
  # MUREN スケーリング係数（幾何平均=1へリスケール）
  r <- d$counts
  c <- MUREN::muren_norm(
    r,
    workers = max(1, parallel::detectCores() - 1),
    res_return = "scaling_coeff"
  )
  if (any(!is.finite(c)) || any(c <= 0)) {
    stop("MUREN scaling coefficients contain non-finite or non-positive values.")
  }
  geo_mean <- exp(mean(log(c)))
  c <- c / geo_mean
  d$samples$norm.factors <- c
  
  # 一貫した size（mean=1）：lib.size * norm.factors を平均1にスケール
  size <- d$samples$lib.size * d$samples$norm.factors
  size <- size / mean(size)
  
  # ペア＋条件
  sample_pair <- as.factor(rep(seq_len(ncol_counts), 2))
  sample_condition <- as.factor(c(rep(0, ncol_counts), rep(1, ncol_counts)))
  
  design <- stats::model.matrix(~ 0 + sample_condition + sample_pair)
  
  # voom は分布補正として quantile を維持
  v <- limma::voom(d, design, normalize.method = "quantile")
  fit <- limma::lmFit(v, design)
  contr <- "sample_condition1 - sample_condition0"
  contrs <- limma::makeContrasts(contrasts = contr, levels = colnames(design))
  fit2 <- limma::contrasts.fit(fit, contrs)
  fit2 <- limma::eBayes(fit2, trend = TRUE)
  
  p_limma <- fit2$p.value[, 1]
  log2_fc_limma <- fit2$coefficients[, 1]
  
  list(
    counts = counts,
    size = size,
    p.limma = p_limma,
    log2FC.limma = log2_fc_limma
  )
}

# ---- Purity estimation only (lightweight contamDE) -------------------------
contamde_purity <- function(counts,
                            subtype = NULL,
                            covariate = NULL,
                            contaminated = TRUE,
                            verbose = TRUE) {
  
  if (verbose) cat("Starting tumor purity estimation...\n")
  
  counts <- as.matrix(counts)
  ncol_counts <- ncol(counts) / 2  # paired samples
  nrow_counts <- nrow(counts)      # genes
  
  if (verbose) cat(sprintf("Input: %d genes, %d paired samples\n", 
                           nrow_counts, ncol_counts))
  
  # voom + limma（p値と log2FC の初期情報）
  if (verbose) cat("Running limma+voom for initial estimates...\n")
  d <- limma_voom_purity(counts = counts)
  p_limma <- d$p.limma
  log2_fc_limma <- d$log2FC.limma
  
  ## design matrix
  if (is.null(subtype) || length(unique(subtype)) == 1) {
    subtype <- rep(1, ncol_counts)
    n_sub <- 1
    if (is.null(covariate)) {
      design <- matrix(subtype, ncol = 1)
    } else {
      covariate <- matrix(covariate, ncol = ncol_counts)
      design <- stats::model.matrix(~ covariate)
    }
  } else {
    subtype <- factor(subtype)
    n_sub <- length(unique(subtype))
    if (is.null(covariate)) {
      design <- stats::model.matrix(~ 0 + subtype)
    } else {
      covariate0 <- matrix(covariate, nrow = ncol_counts)
      design <- stats::model.matrix(~ 0 + subtype + covariate0)
      rm(covariate0)
    }
  }
  ncol_design <- ncol(design)
  
  # sample size factor（limma_voom の出力を使用）
  size <- d$size
  count_norm   <- t(t(counts) / size)
  count_normal <- count_norm[,  seq_len(ncol_counts)] + 1
  count_tumor  <- count_norm[, -seq_len(ncol_counts)] + 1
  
  # log2 ratio
  y <- log2(count_tumor) - log2(count_normal)
  
  if (verbose) cat("Calculating log2 tumor/normal ratios...\n")
  
  # purity proportion estimates
  w_hat <- rep(1, ncol_counts)
  
  if (contaminated) {
    if (verbose) cat("Estimating tumor purity proportions...\n")
    
    # qvalueは必須（ユーザー要件）
    tryCatch({
      p_adj <- qvalue::qvalue(p_limma, pi0.method = "bootstrap")$qvalues
    }, error = function(e) {
      stop("qvalue calculation failed. This is required for purity estimation.")
    })
    
    log2_fc <- log2_fc_limma
    if (sum(p_adj < 0.1) > 1e3) p_adj[-order(p_adj)[1:1e3]] <- 1
    
    up   <- which(p_adj < 0.1 & log2_fc >  log2(1.5))
    down <- which(p_adj < 0.1 & log2_fc < -log2(1.5))
    
    if (verbose) {
      cat(sprintf("Informative genes - Up: %d, Down: %d\n", 
                  length(up), length(down)))
    }
    
    y_up <- if (length(up))   y[up, , drop = FALSE]   else matrix(0, 0, ncol_counts)
    y_down <- if (length(down)) y[down, , drop = FALSE] else matrix(0, 0, ncol_counts)
    
    sumup <- colSums(y_up)
    sumdown <- colSums(y_down)
    sum.max <- max(sumup - sumdown)
    
    if (!is.finite(sum.max) || sum.max <= 0) {
      warning("No informative genes found for purity estimation. Setting purity to 1.0.")
      w_hat <- rep(1.0, ncol_counts)
    } else {
      w_hat <- (sumup - sumdown) / sum.max
      # Ensure purity is between 0 and 1
      w_hat <- pmax(0, pmin(1, w_hat))
    }
    
    if (verbose) {
      cat("Purity estimation summary:\n")
      print(summary(w_hat))
    }
  } else {
    if (verbose) cat("Contamination correction disabled. Setting purity to 1.0.\n")
  }
  
  # Validation checks
  if (length(w_hat) != ncol_counts) {
    stop(sprintf("Purity vector length (%d) != number of pairs (%d).", 
                 length(w_hat), ncol_counts))
  }
  
  if (any(!is.finite(w_hat))) {
    stop("Non-finite values in purity estimates.")
  }
  
  if (verbose) {
    cat(sprintf("Purity estimation completed successfully!\n"))
    cat(sprintf("Mean purity: %.3f, Range: [%.3f, %.3f]\n",
                mean(w_hat), min(w_hat), max(w_hat)))
  }
  
  # Return only essential purity information
  list(
    proportion = w_hat,                    # Main output: purity estimates
    counts = counts,                       # Original count data
    size = size,                          # Size factors used
    design = design,                      # Design matrix
    y = y,                               # Log2 ratios
    n_pairs = ncol_counts,               # Number of pairs
    n_genes = nrow_counts,               # Number of genes
    informative_genes = list(up = up, down = down),  # Genes used for estimation
    normalization_method = "MUREN",       # Normalization used
    estimation_date = Sys.time()         # Timestamp
  )
}

# ---- Quality control function for purity estimates --------------------------
assess_purity_quality <- function(purity_result, threshold = 0.6, verbose = TRUE) {
  
  if (!"proportion" %in% names(purity_result)) {
    stop("Input must be result from contamde_purity function.")
  }
  
  proportions <- purity_result$proportion
  n_samples <- length(proportions)
  
  # Basic statistics
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
    cat("=== Purity Quality Assessment ===\n")
    cat(sprintf("Total samples: %d\n", n_samples))
    cat(sprintf("Mean purity: %.3f ± %.3f\n", 
                purity_stats$mean_purity, purity_stats$sd_purity))
    cat(sprintf("Purity range: [%.3f, %.3f]\n", 
                purity_stats$min_purity, purity_stats$max_purity))
    cat(sprintf("High purity samples (≥%.1f): %d/%d (%.1f%%)\n",
                threshold, n_high_purity, n_samples, retention_rate * 100))
    
    if (retention_rate >= 0.7) {
      cat("✅ Excellent retention rate (≥70%)\n")
    } else if (retention_rate >= 0.5) {
      cat("⚠️ Acceptable retention rate (50-70%)\n")
    } else {
      cat("❌ Poor retention rate (<50%)\n")
    }
  }
  
  return(purity_stats)
}

# ---- Helper function to create filtered sample lists ------------------------
create_purity_filtered_lists <- function(original_sample_lists, purity_results, 
                                         threshold = 0.6, verbose = TRUE) {
  
  if (verbose) cat("Creating purity-filtered sample lists...\n")
  
  filtered_lists <- list()
  
  for (group_name in names(original_sample_lists)) {
    if (verbose) cat(sprintf("Processing group %s...\n", group_name))
    
    if (!group_name %in% names(purity_results) || 
        is.null(purity_results[[group_name]])) {
      if (verbose) cat(sprintf("No purity results for %s, skipping...\n", group_name))
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
      if (verbose) cat(sprintf("No high purity samples in %s\n", group_name))
      filtered_lists[[group_name]] <- list(
        tumor = character(0),
        normal = character(0),
        cases = character(0)
      )
      next
    }
    
    # Extract original sample information
    original_group <- original_sample_lists[[group_name]]
    n_original <- length(original_group$tumor)
    
    # Filter samples based on high purity indices
    filtered_tumor <- original_group$tumor[high_purity_indices]
    filtered_normal <- original_group$normal[high_purity_indices]
    filtered_cases <- original_group$cases[high_purity_indices]
    
    filtered_lists[[group_name]] <- list(
      tumor = filtered_tumor,
      normal = filtered_normal,
      cases = filtered_cases
    )
    
    if (verbose) {
      cat(sprintf("%s: %d → %d samples (%.1f%% retention)\n",
                  group_name, n_original, length(filtered_tumor),
                  length(filtered_tumor) / n_original * 100))
    }
  }
  
  if (verbose) cat("Purity filtering completed!\n")
  return(filtered_lists)
}

cat("ContamDE purity estimation functions loaded successfully!\n")
cat("Main function: contamde_purity()\n")
cat("Quality assessment: assess_purity_quality()\n")
cat("Sample filtering: create_purity_filtered_lists()\n")