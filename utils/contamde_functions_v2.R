# ---- limma + voom (quantile維持) with MUREN scaling factors -----------------
limma_voom <- function(counts) {
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

# ---- Empirical Bayes for residual variances (robust/non-robust) --------------
ebayes <- function(var_res, dg, winsor_tail_p = c(0.05, 0.1), robust = TRUE) {
  n <- length(var_res)
  df1 <- dg
  prob <- winsor_tail_p
  prob[2] <- 1 - winsor_tail_p[2]
  xq <- stats::quantile(var_res, prob)
  
  if (!robust) {
    zg <- log(var_res)
    z_mean <- mean(zg)
    z_var <- stats::var(zg)
    
    fund0 <- function(x) {
      psigamma(x / 2, deriv = 1) + psigamma(dg / 2, deriv = 1) - z_var
    }
    d0 <- stats::uniroot(f = fund0, lower = 0.1, upper = 500)$root
    log.s02 <- digamma(d0 / 2) - log(d0 / dg) - digamma(dg / 2) + z_mean
    s02 <- exp(log.s02)
    df <- d0 + dg
    sigma_g2_hat <- (d0 * s02 + dg * var_res) / df
    return(list(sigma.g2.hat = sigma_g2_hat, d0 = d0, s02 = s02))
  }
  
  # robust
  win_s_g <- pmax(pmin(var_res, xq[2]), xq[1])
  zg <- log(win_s_g)
  z_mean <- mean(zg)
  z_var <- stats::var(zg)
  
  sub_fun <- function(x) x / (1 + x)
  sub_inv <- function(x) x / (1 - x)
  gauss_quad <- statmod::gauss.quad.prob(128, dist = "uniform")
  
  win_f_moments <- function(df1, df2) {
    fq <- stats::qf(
      p = c(winsor_tail_p[1], 1 - winsor_tail_p[2]),
      df1 = df1, df2 = df2
    )
    zq <- log(fq)
    q <- sub_fun(fq)
    q21 <- q[2] - q[1]
    nodes <- q[1] + (q[2] - q[1]) * gauss_quad$nodes
    inv_nodes <- sub_inv(nodes)
    pdf_nodes <- stats::df(inv_nodes, df1 = df1, df2 = df2)
    log.inv.nodes <- log(inv_nodes)
    h_mean <- log.inv.nodes * pdf_nodes / (1 - nodes)^2
    mean <- q21 * sum(gauss_quad$weight * h_mean) + sum(zq * winsor_tail_p)
    h_var <- (log.inv.nodes - mean)^2 / (1 - nodes)^2 * pdf_nodes
    var <- q21 * sum(gauss_quad$weight * h_var) + sum((zq - mean)^2 * winsor_tail_p)
    list(mean = mean, var = var)
  }
  
  inf_mom <- win_f_moments(df1 = dg, df2 = Inf)
  ub_fun <- log(z_var / inf_mom$var)
  
  if (ub_fun <= 0) {
    df2_corrected <- rep(Inf, n)
    log.s02 <- z_mean - inf_mom$mean
    s02 <- exp(log.s02)
    sigma_g2_hat <- rep(0, n)
    return(list(sigma.g2.hat = sigma_g2_hat, d0 = df2_corrected, s02 = s02))
  }
  
  non_robust <- ebayes(var_res, dg, winsor_tail_p = winsor_tail_p, robust = FALSE)
  fun <- function(x) {
    df2 <- sub_inv(x)
    mom <- win_f_moments(df1 = df1, df2 = df2)
    log(z_var / mom$var)
  }
  lbx <- sub_fun(non_robust$d0)
  lb_fun <- fun(lbx)
  
  if (lb_fun >= 0) {
    df2 <- non_robust$d0
  } else {
    u <- stats::uniroot(fun, interval = c(lbx, 1), tol = 1e-08,
                        f.lower = lb_fun, f.upper = ub_fun)
    df2 <- sub_inv(u$root)
  }
  
  mom <- win_f_moments(df1 = df1, df2 = df2)
  log.s02 <- z_mean - mom$mean
  s02 <- exp(log.s02)
  
  f_stat <- exp(log(var_res) - log.s02)
  log.p.value <- stats::pf(f_stat, df1 = df1, df2 = df2,
                           lower.tail = FALSE, log.p = TRUE)
  
  r <- rank(var_res)
  log.prior.p <- log(n - r + 0.5) - log(n)
  log.prob.not.outlier <- pmin(log.p.value - log.prior.p, 0)
  prob_not_outlier <- exp(log.prob.not.outlier)
  prob_outlier <- -expm1(log.prob.not.outlier)
  
  likelihood_fun <- function(x) {
    stats::df(max(var_res) / s02, df1 = df1, df2 = x, log = TRUE)
  }
  df2_outlier <- stats::optimise(likelihood_fun, c(0, df2), maximum = TRUE)$maximum
  df2_corrected <- prob_not_outlier * df2 + prob_outlier * df2_outlier
  
  if (any(log.prob.not.outlier < 0)) {
    o <- order(log.p.value)
    df2_ordered <- df2_corrected[o]
    cummean_df2 <- cumsum(df2_ordered) / seq_along(df2_ordered)
    min.index <- which.min(cummean_df2)
    df2_ordered[1:min.index] <- cummean_df2[min.index]
    df2_corrected[o] <- cummax(df2_ordered)
  } else {
    df2_corrected <- rep.int(df2, n)
  }
  
  df <- df2_corrected + dg
  sigma_g2_hat <- (df2_corrected * s02 + dg * var_res) / df
  list(sigma.g2.hat = sigma_g2_hat, d0 = df2_corrected, s02 = s02)
}

# ---- contamDE-lm main (qvalue維持、GAM予測一括、回帰ループ最適化) -----------
contamde_lm <- function(counts,
                        subtype = NULL,
                        covariate = NULL,
                        contaminated = TRUE,
                        robust = TRUE) {
  counts <- as.matrix(counts)
  ncol_counts <- ncol(counts) / 2  # paired samples
  nrow_counts <- nrow(counts)      # genes
  
  # voom + limma（p値と log2FC の初期情報）
  d <- limma_voom(counts = counts)
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
  
  # purity proportion estimates (optional)
  w_hat <- rep(1, ncol_counts)
  if (contaminated) {
    # qvalueは必須（ユーザー要件）
    p_adj <- qvalue::qvalue(p_limma, pi0.method = "bootstrap")$qvalues
    log2_fc <- log2_fc_limma
    if (sum(p_adj < 0.1) > 1e3) p_adj[-order(p_adj)[1:1e3]] <- 1
    
    up   <- which(p_adj < 0.1 & log2_fc >  log2(1.5))
    down <- which(p_adj < 0.1 & log2_fc < -log2(1.5))
    y_up <- if (length(up))   y[up, , drop = FALSE]   else matrix(0, 0, ncol_counts)
    y_down <- if (length(down)) y[down, , drop = FALSE] else matrix(0, 0, ncol_counts)
    
    sumup <- colSums(y_up)
    sumdown <- colSums(y_down)
    sum.max <- max(sumup - sumdown)
    if (!is.finite(sum.max) || sum.max <= 0) {
      stop("Invalid sum.max in purity estimation (no informative genes?).")
    }
    w_hat <- (sumup - sumdown) / sum.max
  }
  
  # ... w_hat を計算し終えた直後に挿入 ...
  if (length(w_hat) != nrow(design)) {
    stop(sprintf("w_hat length (%d) != nrow(design) (%d).", length(w_hat), nrow(design)))
  }
  if (ncol(counts) / 2 != nrow(design)) {
    stop("ncol(counts)/2 と design の行数が一致しません。ペア設計を確認してください。")
  }
  
  x_w <- sweep(design, 1, w_hat, `*`)  # ← ここで行方向に重み付け
  dg <- ncol_counts - ncol_design
  if (dg <= 0) stop("Degrees of freedom (dg) must be positive.")
  
  ## residual std (projection)
  # res_ig: (I - X(X'X)^-1 X') * Y  with weighting implied later
  # ここは元設計どおりだが、X'Xの逆行列は下で毎回扱うので最終的な検定では再計算する
  res_ig <- (diag(ncol_counts) - x_w %*% solve(crossprod(x_w)) %*% t(x_w)) %*% t(y)
  var_res <- colSums(res_ig^2) / dg
  if (any(!is.finite(var_res)) || any(var_res <= 0)) {
    stop("Non-finite or non-positive residual variances.")
  }
  
  log_vg <- log(var_res)
  
  m_coverage <- (count_normal + count_tumor) / 2
  log_c <- log(m_coverage)
  mlog_c <- rowMeans(log_c)
  
  ## GAM fit（予測はベクトル一括→行列へ整形）
  gamres <- mgcv::gam(
    y ~ s(x, k = 4),
    data = data.frame(x = mlog_c, y = log_vg)
  )
  pred_vec <- stats::predict(gamres, newdata = data.frame(x = as.vector(log_c)))
  if (length(pred_vec) != length(log_c)) {
    stop("Length mismatch in GAM prediction.")
  }
  predictors <- matrix(pred_vec, nrow = nrow_counts, ncol = ncol_counts, byrow = FALSE)
  
  v_ig <- exp(predictors)
  weight <- 1 / v_ig
  weight <- weight / rowMeans(weight)
  
  ## EB estimation on var_res
  ebayes_result <- ebayes(var_res = var_res, dg = dg, robust = robust)
  hatsg2 <- ebayes_result$sigma.g2.hat
  df <- ebayes_result$d0 + dg
  
  ## Weighted regression per gene（高速化：diag(W)回避＋Cholesky）
  res <- matrix(0, nrow_counts, (ncol_design + ncol_design^2))
  Xt <- t(x_w)
  
  for (g in seq_len(nrow_counts)) {
    wg <- weight[g, ]
    if (any(!is.finite(wg)) || any(wg <= 0)) {
      stop(sprintf("Non-finite or non-positive weights at gene %d.", g))
    }
    y_g <- y[g, ]
    
    # X'WX and X'Wy without forming diag(W)
    X_w <- sweep(x_w, 1, wg, `*`)     # (N x K)
    XWX <- crossprod(x_w, X_w)        # (K x K)
    XWy <- crossprod(x_w, wg * y_g)   # (K)
    
    R <- tryCatch(chol(XWX), error = function(e) NULL)
    if (is.null(R)) {
      stop(sprintf("chol(XWX) failed at gene %d (matrix not SPD).", g))
    }
    thetag_hat <- backsolve(R, forwardsolve(t(R), XWy))
    
    XWX_inv <- chol2inv(R)
    cov_theta <- XWX_inv * hatsg2[g]
    
    res[g, ] <- c(thetag_hat, as.vector(cov_theta))
  }
  
  # t-stats and p-values
  diag_id <- (0:(ncol_design - 1)) * ncol_design + (1:ncol_design)
  log2_fc <- matrix(res[, (1:ncol_design), drop = FALSE], nrow = nrow_counts)
  log2_fc_cov <- matrix(res[, -(1:ncol_design), drop = FALSE], nrow = nrow_counts)
  t_stat <- abs(log2_fc) / sqrt(log2_fc_cov[, diag_id, drop = FALSE])
  
  p_eb <- matrix(0, nrow_counts, ncol_design)
  for (i in seq_len(ncol_design)) {
    p_eb[, i] <- 2 * (1 - stats::pt(t_stat[, i], df = df))
    large_t <- which(t_stat[, i] > 10)
    if (length(large_t) > 0) {
      p_eb[large_t, i] <- 2 * (-stats::pt(t_stat[large_t, i], df = df[large_t], log.p = TRUE))
    }
  }
  
  if (is.null(covariate)) {
    colnames(p_eb) <- paste0("subtype-", seq_len(n_sub), " vs. normal")
  } else {
    colnames(p_eb) <- c(
      paste0("subtype-", seq_len(n_sub), " vs. normal"),
      paste0("covariate-", seq_len(ncol(covariate)))
    )
  }
  
  list(
    counts = d$counts,
    p.contamDE.lm = p_eb,
    log2FC = log2_fc,
    log2_fc_cov = log2_fc_cov,
    proportion = w_hat,
    design = design,
    df = df,
    weight = weight,
    y = y
  )
}
