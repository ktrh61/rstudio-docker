# ============================================================================== 
# ContamDE Purity Estimation Functions (Lightweight Version)
# contamde_purity_functions.R
# ==============================================================================
# 
# Purpose: Tumor purity estimation only (no DEG analysis)
# Based on: contamde_functions_v2.R
# Modifications: Removed DEG analysis, focused on proportion estimation

# ---- limma + voom (quantile維持) with MUREN scaling factors -----------------
# ---- limma + voom (quantile維持) with MUREN scaling factors -----------------
# counts: genes x samples（前半 = normal, 後半 = tumor のペア順）
# refs:   MURENの参照（デフォルトは全サンプル）
# workers: "auto" か整数
# voom_norm: "quantile" のままにしておく（変更したければ "none" に）
limma_voom_purity <- function(counts, refs = "saturated",
                              workers = "auto", voom_norm = "quantile") {
  if (!is.matrix(counts)) counts <- as.matrix(counts)
  ncol_counts <- ncol(counts) / 2
  stopifnot(ncol_counts == floor(ncol(counts) / 2))  # ペア前提
  
  # edgeR: DGEList（TMMは使わない。lib.size は列和で入る）
  d <- edgeR::DGEList(counts = counts)
  
  # MUREN スケーリング係数（“counts を割る”側）
  c <- muren_norm(
    counts,
    refs        = refs,
    single_param= TRUE,
    res_return  = "scaling_coeff",
    workers     = workers
  )
  if (any(!is.finite(c)) || any(c <= 0)) {
    stop("MUREN scaling coefficients contain non-finite or non-positive values.")
  }
  
  # 定数倍のリスケール（幾何平均=1）は相対関係を壊さず安全
  geo_mean <- exp(mean(log(c)))
  c <- c / geo_mean
  
  # edgeR に適用（有効サイズ = lib.size * norm.factors）
  d$samples$norm.factors <- as.numeric(c)
  
  # voom（分布補正として quantile を維持したいなら voom_norm="quantile"）
  design <- stats::model.matrix(~ 0 + factor(c(rep(0, ncol_counts), rep(1, ncol_counts))))
  v <- limma::voom(d, design, normalize.method = voom_norm)
  
  # 係数ベクトル（のちの手計算用に）: mean 1 に揃えた有効サイズ
  size <- d$samples$lib.size * d$samples$norm.factors
  size <- size / mean(size)
  
  # limma で p と log2FC（情報遺伝子の選抜にだけ使う）
  fit <- limma::lmFit(v, design)
  contrs <- limma::makeContrasts(contrasts = "factor...1 - factor...0",
                                 levels = colnames(design))
  fit2 <- limma::contrasts.fit(fit, contrs)
  fit2 <- limma::eBayes(fit2, trend = TRUE)
  
  list(
    counts         = counts,
    size           = size,
    p.limma        = fit2$p.value[, 1],
    log2FC.limma   = fit2$coefficients[, 1]
  )
}

# ---- Purity estimation only (lightweight contamDE) -------------------------
contamde_purity <- function(counts,
                            subtype = NULL,
                            covariate = NULL,
                            contaminated = TRUE,
                            muren_refs = "saturated",    # ★ 参照を外から指定可
                            workers = "auto",
                            voom_norm = "quantile",
                            prior.count = 0,             # 0のまま=従来どおり全ゼロを除外
                            verbose = TRUE) {
  
  if (verbose) cat("Starting tumor purity estimation...\n")
  
  counts <- as.matrix(counts)
  ncol_counts <- ncol(counts) / 2
  stopifnot(ncol_counts == floor(ncol(counts) / 2))
  nrow_counts <- nrow(counts)
  
  if (verbose) cat(sprintf("Input: %d genes, %d paired samples\n", nrow_counts, ncol_counts))
  
  # limma+voom（p/logFC のみ利用）＋ MURENサイズ（size）を取得
  if (verbose) cat("Running limma+voom (MUREN-scaled) for initial estimates...\n")
  d <- limma_voom_purity(counts = counts, refs = muren_refs,
                         workers = workers, voom_norm = voom_norm)
  p_limma        <- d$p.limma
  log2_fc_limma  <- d$log2FC.limma
  size           <- d$size
  
  # MURENサイズで正規化（mean(size)=1）
  count_norm <- t(t(counts) / size)
  
  # ゼロ処理：従来と同じく「全サンプル>0」を既定に。
  # もし prior.count を入れたいなら prior.count>0 を渡すとゼロ除外を回避できる。
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
    cat(sprintf("Excluded %d genes with zero counts (%.2f%%)\n", n_exc, 100 * n_exc / nrow(count_norm)))
  }
  
  # normal / tumor に分割して log2 比を作る
  idx_n <- seq_len(ncol_counts)
  idx_t <- setdiff(seq_len(ncol(count_norm_valid)), idx_n)
  
  if (prior.count > 0) {
    y <- log2((count_norm_valid[, idx_t, drop=FALSE] + prior.count) /
                (count_norm_valid[, idx_n, drop=FALSE] + prior.count))
  } else {
    y <- log2(count_norm_valid[, idx_t, drop=FALSE]) -
      log2(count_norm_valid[, idx_n, drop=FALSE])
  }
  
  # デザイン（従来どおりの簡易形）
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
  
  # 汚染あり想定の純度推定
  w_hat <- rep(1, ncol_counts)
  if (contaminated) {
    if (verbose) cat("Estimating tumor purity proportions...\n")
    
    p_limma       <- p_limma[valid_genes]
    log2_fc_limma <- log2_fc_limma[valid_genes]
    
    # qvalue 必須
    p_adj <- tryCatch(qvalue::qvalue(p_limma, pi0.method = "bootstrap")$qvalues,
                      error = function(e) stop("qvalue calculation failed: ", e$message))
    
    # 情報遺伝子の定義（従来どおり）
    log2_fc <- log2_fc_limma
    if (sum(p_adj < 0.1) > 1e3) p_adj[-order(p_adj)[1:1e3]] <- 1
    up   <- which(p_adj < 0.1 & log2_fc >  log2(1.5))
    down <- which(p_adj < 0.1 & log2_fc < -log2(1.5))
    
    if (verbose) cat(sprintf("Informative genes - Up: %d, Down: %d\n", length(up), length(down)))
    
    y_up   <- if (length(up))   y[up,   , drop = FALSE] else matrix(0, 0, ncol_counts)
    y_down <- if (length(down)) y[down, , drop = FALSE] else matrix(0, 0, ncol_counts)
    
    sumup   <- colSums(y_up)
    sumdown <- colSums(y_down)
    sum.max <- max(sumup - sumdown)
    
    if (!is.finite(sum.max) || sum.max <= 0) {
      warning("No informative genes found for purity estimation. Setting purity to 1.0.")
      w_hat <- rep(1.0, ncol_counts)
    } else {
      w_hat <- (sumup - sumdown) / sum.max
      w_hat <- pmax(0, pmin(1, w_hat))
    }
    
    if (verbose) {
      cat("Purity estimation summary:\n")
      print(summary(w_hat))
    }
  } else if (verbose) {
    cat("Contamination correction disabled. Setting purity to 1.0.\n")
  }
  
  # 最終チェック
  stopifnot(length(w_hat) == ncol_counts, all(is.finite(w_hat)))
  
  if (verbose) {
    cat("Purity estimation completed successfully!\n")
    cat(sprintf("Mean purity: %.3f, Range: [%.3f, %.3f]\n", mean(w_hat), min(w_hat), max(w_hat)))
  }
  
  list(
    proportion         = w_hat,
    counts             = counts_valid,
    size               = size,
    design             = design,
    y                  = y,
    n_pairs            = ncol_counts,
    n_genes            = nrow(counts_valid),
    informative_genes  = list(up = up, down = down),
    normalization_method = "MUREN",
    estimation_date    = Sys.time()
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