suppressPackageStartupMessages({
  library(Rcpp)
  library(edgeR)
  library(limma)
  library(mgcv)
  library(qvalue)
  library(statmod)
})

## すでに: Rcpp::sourceCpp("src/contamde_gls.cpp") を実行しておくこと
## gls_per_gene_cpp() が使える状態になっている前提

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
  
  # 一貫した size（mean=1）
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
    d0 <- stats::uniroot(fund0, c(0.1, 500))$root
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
    u <- stats::uniroot(
      fun,
      interval = c(lbx, 1), tol = 1e-08, f.lower = lb_fun, f.upper = ub_fun
    )
    df2 <- sub_inv(u$root)
  }
  
  mom <- win_f_moments(df1 = df1, df2 = df2)
  log.s02 <- z_mean - mom$mean
  s02 <- exp(log.s02)
  
  f_stat <- exp(log(var_res) - log.s02)
  log.p.value <- stats::pf(
    f_stat,
    df1 = df1, df2 = df2, lower.tail = FALSE, log.p = TRUE
  )
  
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

# ---- contamDE-lm main (qvalue維持、GAM予測一括、GLSはRcpp) -------------------
contamde_lm <- function(counts,
                        subtype = NULL,
                        covariate = NULL,
                        contaminated = TRUE,
                        robust = TRUE) {
  counts <- as.matrix(counts)
  ncol_counts <- ncol(counts) / 2
  nrow_counts <- nrow(counts)
  
  d <- limma_voom(counts = counts)
  p_limma <- d$p.limma
  log2_fc_limma <- d$log2FC.limma
  
  ## design
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
  
  ## size 正規化
  size <- d$size
  count_norm   <- t(t(counts) / size)
  count_normal <- count_norm[,  seq_len(ncol_counts)] + 1
  count_tumor  <- count_norm[, -seq_len(ncol_counts)] + 1
  
  ## log2 ratio
  y <- log2(count_tumor) - log2(count_normal)
  
  ## purity proportion
  w_hat <- rep(1, ncol_counts)
  if (contaminated) {
    p_adj <- qvalue::qvalue(p_limma, pi0.method = "bootstrap")$qvalues
    log2_fc <- log2_fc_limma
    if (sum(p_adj < 0.1) > 1e3) p_adj[-order(p_adj)[1:1e3]] <- 1
    
    up   <- which(p_adj < 0.1 & log2_fc >  log2(1.5))
    down <- which(p_adj < 0.1 & log2_fc < -log2(1.5))
    y_up <- if (length(up))     y[up, , drop = FALSE]   else matrix(0, 0, ncol_counts)
    y_down <- if (length(down)) y[down, , drop = FALSE] else matrix(0, 0, ncol_counts)
    
    sumup <- colSums(y_up); sumdown <- colSums(y_down)
    sum.max <- max(sumup - sumdown)
    if (!is.finite(sum.max) || sum.max <= 0) {
      stop("Invalid sum.max in purity estimation (no informative genes?).")
    }
    w_hat <- (sumup - sumdown) / sum.max
  }
  
  # ガード & x_w
  if (length(w_hat) != nrow(design)) {
    stop(sprintf("w_hat length (%d) != nrow(design) (%d).", length(w_hat), nrow(design)))
  }
  if (ncol(counts) / 2 != nrow(design)) {
    stop("ncol(counts)/2 と design の行数が一致しません。ペア設計を確認してください。")
  }
  x_w <- sweep(design, 1, w_hat, `*`)
  
  dg <- ncol_counts - ncol_design
  if (dg <= 0) stop(sprintf("Degrees of freedom <= 0 (N=%d, K=%d).", ncol_counts, ncol_design))
  
  ## EB 用の残差分散
  P <- diag(ncol_counts) - x_w %*% solve(crossprod(x_w)) %*% t(x_w)
  res_ig <- P %*% t(y)
  var_res <- colSums(res_ig^2) / dg
  if (any(!is.finite(var_res)) || any(var_res <= 0)) {
    stop("Non-finite or non-positive residual variances. x_w / weights を確認。")
  }
  
  log_vg <- log(var_res)
  m_coverage <- (count_normal + count_tumor) / 2
  log_c <- log(m_coverage)
  mlog_c <- rowMeans(log_c)
  
  ## GAM（予測は一括）
  gamres <- mgcv::gam(y ~ s(x, k = 4), data = data.frame(x = mlog_c, y = log_vg))
  pred_vec <- stats::predict(gamres, newdata = data.frame(x = as.vector(log_c)))
  if (length(pred_vec) != length(log_c)) stop("Length mismatch in GAM prediction.")
  predictors <- matrix(pred_vec, nrow = nrow_counts, ncol = ncol_counts, byrow = FALSE)
  
  v_ig <- exp(predictors)
  weight <- 1 / v_ig
  weight <- weight / rowMeans(weight)
  
  ## EB 推定
  ebayes_result <- ebayes(var_res = var_res, dg = dg, robust = robust)
  hatsg2 <- ebayes_result$sigma.g2.hat
  df <- ebayes_result$d0 + dg
  
  ## ---- Rcpp: GLS per gene ----
  cpp_out <- gls_per_gene_cpp(
    X = x_w,
    W = weight,
    Y = y,
    sigma2 = hatsg2,
    return_full_cov = TRUE
  )
  
  log2_fc <- cpp_out$coef                    # (G x K)
  K <- ncol_design; G <- nrow_counts
  cov_cube <- cpp_out$cov_full               # (K x K x G)
  log2_fc_cov <- matrix(0, G, K * K)
  for (g in seq_len(G)) log2_fc_cov[g, ] <- as.vector(cov_cube[, , g])
  
  ## t統計 & p
  diag_id <- (0:(K - 1)) * K + (1:K)
  t_stat <- abs(log2_fc) / sqrt(log2_fc_cov[, diag_id, drop = FALSE])
  
  p_eb <- matrix(0, G, K)
  for (i in seq_len(K)) {
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
