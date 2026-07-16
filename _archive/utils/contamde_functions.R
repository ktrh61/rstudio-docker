limma_voom <- function(counts) {
  d <- edgeR::DGEList(counts)
  d <- edgeR::calcNormFactors(d)
  
  r <- d$counts
  c <- MUREN::muren_norm(r,
                         workers = parallel::detectCores() - 1, res_return = "scaling_coeff"
  )
  d$samples$norm.factors <- c
  
  size <- d$samples[, 2] * d$samples[, 3]
  size <- size / mean(size)
  
  ncol_counts <- ncol(counts) / 2
  sample_pair <- as.factor(rep(1:ncol_counts, 2)) # nolint: indentation_linter.
  sample_condition <- as.factor(c(rep(0, ncol_counts), rep(1, ncol_counts)))
  
  design <- model.matrix(~ 0 + sample_condition + sample_pair)
  d <- limma::voom(d, design, normalize.method = "quantile")
  d <- limma::lmFit(d, design)
  contr <- "sample_condition1-sample_condition0"
  contrs <- limma::makeContrasts(contrasts = contr, levels = colnames(design))
  d <- limma::contrasts.fit(d, contrs)
  d <- limma::eBayes(d, trend = TRUE)
  p_limma <- d$p.value
  log2_fc_limma <- d$coefficient
  res <- list(
    counts = counts,
    size = size,
    p.limma = p_limma,
    log2FC.limma = log2_fc_limma
  )
  return(res)
  rm(sample_pair, sample_condition)
}

###  different estimation of deviations and degrees of freedom
ebayes <- function(var_res, dg, winsor_tail_p = c(0.05, 0.1), robust = TRUE) {
  n <- length(var_res)
  df1 <- dg
  prob <- winsor_tail_p
  prob[2] <- 1 - winsor_tail_p[2]
  xq <- quantile(var_res, prob)
  
  ## non robust  (without outliers)
  if (!robust) {
    zg <- log(var_res) ## fisher's z distribution
    z_mean <- mean(zg)
    z_var <- var(zg)
    
    fund0 <- function(x) {
      trigamma(x / 2) + trigamma(dg / 2) - z_var
    }
    
    d0 <- uniroot(fund0, c(0.1, 500))$root
    
    ## the estimation of s0^2
    
    log.s02 <- digamma(d0 / 2) - log(d0 / dg) - digamma(dg / 2) + z_mean
    s02 <- exp(log.s02)
    
    ## residual variance EB estimation of sigma_g^2
    df <- d0 + dg
    sigma_g2_hat <- (d0 * s02 + dg * var_res) / df
    
    result <- list(sigma.g2.hat = sigma_g2_hat, d0 = d0, s02 = s02)
    return(result)
  }
  
  ## robust (with outliers)
  win_s_g <- pmax(pmin(var_res, xq[2]), xq[1])
  zg <- log(win_s_g)
  z_mean <- mean(zg)
  z_var <- var(zg)
  
  sub_fun <- function(x) x / (1 + x)
  sub_inv <- function(x) x / (1 - x)
  gauss_quad <- statmod::gauss.quad.prob(128, dist = "uniform")
  
  win_f_moments <- function(df1 = df1, df2 = df2) {
    fq <- qf(
      p = c(winsor_tail_p[1], 1 - winsor_tail_p[2]),
      df1 = df1, df2 = df2
    )
    zq <- log(fq)
    q <- sub_fun(fq)
    q21 <- q[2] - q[1]
    nodes <- q[1] + (q[2] - q[1]) * gauss_quad$nodes
    inv_nodes <- sub_inv(nodes)
    pdf_nodes <- df(inv_nodes, df1 = df1, df2 = df2)
    log.inv.nodes <- log(inv_nodes)
    h_mean <- log.inv.nodes * pdf_nodes / (1 - nodes)^2
    mean <- q21 * sum(gauss_quad$weight * h_mean) + sum(zq * winsor_tail_p)
    h_var <- (log.inv.nodes - mean)^2 / (1 - nodes)^2 * pdf_nodes
    var <-
      q21 * sum(gauss_quad$weight * h_var) + sum((zq - mean)^2 * winsor_tail_p)
    
    list(mean = mean, var = var)
  }
  
  inf_mom <- win_f_moments(df1 = dg, df2 = Inf)
  ub_fun <- log(z_var / inf_mom$var)
  
  if (ub_fun <= 0) {
    df2 <- rep(Inf, n)
    log.s02 <- z_mean - inf_mom$mean
    s02 <- exp(log.s02)
    sigma_g2_hat <- rep(0, n)
    
    return(list(sigma.g2.hat = sigma_g2_hat, d0 = df2_corrected, s02 = s02))
  }
  
  non_robust <- ebayes(
    var_res, dg,
    winsor_tail_p = c(0.05, 0.1), robust = FALSE
  )
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
    u <- uniroot(
      fun,
      interval = c(lbx, 1), tol = 1e-08, f.lower = lb_fun, f.upper = ub_fun
    )
    df2 <- sub_inv(u$root)
  }
  
  mom <- win_f_moments(df1 = df1, df2 = df2)
  log.s02 <- z_mean - mom$mean
  s02 <- exp(log.s02)
  
  f_stat <- exp(log(var_res) - log.s02)
  
  log.p.value <- pf(
    f_stat,
    df1 = df1, df2 = df2, lower.tail = FALSE, log.p = TRUE
  )
  
  r <- rank(var_res)
  log.prior.p <- log(n - r + 0.5) - log(n)
  log.prob.not.outlier <- pmin(log.p.value - log.prior.p, 0)
  prob_not_outlier <- exp(log.prob.not.outlier)
  prob_outlier <- -expm1(log.prob.not.outlier)
  
  likelihood_fun <- function(x) {
    df(max(var_res) / s02, df1 = df1, df2 = x, log = TRUE)
  }
  
  df2_outlier <- optimise(likelihood_fun, c(0, df2), maximum = TRUE)$maximum
  df2_corrected <- prob_not_outlier * df2 + prob_outlier * df2_outlier
  
  if (any(log.prob.not.outlier < 0)) {
    o <- order(log.p.value)
    df2_ordered <- df2_corrected[o]
    cummean_df2 <- cumsum(df2_ordered) / (1:n)
    min.index <- which.min(cummean_df2)
    
    df2_ordered[1:min.index] <- cummean_df2[min.index]
    df2_corrected[o] <- cummax(df2_ordered)
  } else {
    df2_outlier <- df2
    df2_corrected <- rep.int(df2, n)
  }
  
  df <- df2_corrected + dg
  sigma_g2_hat <- (df2_corrected * s02 + dg * var_res) / df
  
  result <- list(sigma.g2.hat = sigma_g2_hat, d0 = df2_corrected, s02 = s02)
  
  return(result)
}

contamde_lm <- function(counts,
                        subtype = NULL,
                        covariate = NULL,
                        contaminated = TRUE,
                        robust = TRUE) {
  counts <- as.matrix(counts)
  ncol_counts <- ncol(counts) / 2 # the number of paired samples
  nrow_counts <- nrow(counts) # the number of total genes
  
  # use voom strategy in LIMMA package
  # to estimate the initial value of purity proportion
  d <- limma_voom(counts = counts)
  p_limma <- d$p.limma
  log2_fc_limma <- d$log2FC.limma
  
  ## design matrix ##
  # only one tumor subtype (without subtype information)
  if (is.null(subtype) || length(unique(subtype)) == 1) {
    subtype <- rep(1, ncol_counts)
    n_sub <- 1
    if (is.null(covariate)) {
      design <- matrix(subtype, ncol = 1)
    } else {
      covariate <- matrix(covariate, ncol = ncol_counts)
      design <- model.matrix(~covariate)
    }
  } else {
    # with subtype information
    subtype <- factor(subtype)
    n_sub <- length(unique(subtype))
    if (is.null(covariate)) {
      design <- model.matrix(~ 0 + subtype)
    } else {
      covariate0 <- matrix(covariate, nrow = ncol_counts)
      design <- model.matrix(~ 0 + subtype + covariate0)
      rm(covariate0)
    }
  }
  
  ncol_design <- ncol(design)
  
  # sample size factor
  size <- d$size
  count_norm <- t(t(counts) / size)
  count_normal <- count_norm[, (1:ncol_counts)] + 1
  count_tumor <- count_norm[, -(1:ncol_counts)] + 1
  
  # do transformation
  y <- log2(count_tumor) - log2(count_normal)
  
  w_hat <- rep(1, ncol_counts)
  if (contaminated) {
    ##########################################
    # proportion estimation
    ##########################################
    p_adj <- qvalue::qvalue(p_limma, pi0.method = "bootstrap")$qvalues
    log2_fc <- log2_fc_limma
    if (sum(p_adj < 0.1) > 1e3) p_adj[-order(p_adj)[1:1e3]] <- 1
    
    up <- which(p_adj < 0.1 & log2_fc > log2(1.5))
    down <- which(p_adj < 0.1 & log2_fc < -log2(1.5))
    y_up <- y[up, ]
    y_down <- y[down, ]
    
    ## estimate ws
    sumup <- colSums(y_up)
    sumdown <- colSums(y_down)
    sum.max <- max(sumup - sumdown)
    w_hat <- (sumup - sumdown) / sum.max # purity proportion estimates
  }
  
  x_w <- design * w_hat
  dg <- ncol_counts - ncol_design # degree of freedom
  
  ## residual standard deviation
  res_ig <-
    (diag(ncol_counts) - x_w %*% solve((t(x_w) %*% x_w)) %*% t(x_w)) %*% t(y)
  var_res <- colSums(res_ig^2) / dg
  
  log_vg <- log(var_res) # log residual deviations
  
  m_coverage <- (count_normal + count_tumor) / 2
  log_c <- log(m_coverage)
  mlog_c <- rowMeans(log_c)
  
  ## gam fit
  gamres <- mgcv::gam(
    y ~ s(x, k = 4),
    data = data.frame(x = mlog_c, y = log_vg)
  )
  predictors <- apply(log_c, 2, function(x) {
    predict(gamres, data.frame(x = x))
  })
  
  v_ig <- exp(predictors)
  
  weight <- 1 / (v_ig) # the inverse of variances are the precision weights
  weight <- weight / rowMeans(weight)
  
  ## ebayes estimation
  
  ebayes_result <- ebayes(var_res = var_res, dg, robust = robust)
  hatsg2 <- ebayes_result$sigma.g2.hat
  df <- ebayes_result$d0 + dg
  
  ## t-test
  res <- matrix(0, nrow_counts, (ncol_design + ncol_design^2))
  for (g in 1:nrow_counts) {
    y_g <- as.numeric(y[g, ])
    weight_gi <- weight[g, ]
    w_g <- diag(weight_gi)
    xwx_inv <- solve(t(x_w) %*% w_g %*% x_w)
    thetag_hat <- as.numeric(xwx_inv %*% t(x_w) %*% w_g %*% y_g)
    cov_theta <- xwx_inv * hatsg2[g]
    res[g, ] <- c(thetag_hat, cov_theta)
  }
  
  ##########################
  ##########################
  diag_id <- (0:(ncol_design - 1)) * ncol_design + (1:ncol_design)
  log2_fc <- matrix(res[, (1:ncol_design)], nrow = nrow_counts)
  log2_fc_cov <- matrix(res[, -(1:ncol_design)], nrow = nrow_counts)
  t_stat <- abs(log2_fc) / sqrt(log2_fc_cov[, diag_id])
  ##########################
  ##########################
  p_eb <- matrix(0, nrow_counts, ncol_design)
  
  for (i in 1:ncol_design) {
    p_eb[, i] <- 2 * (1 - pt(t_stat[, i], df = df))
    large_t <- which(t_stat[, i] > 10)
    if (length(large_t) > 0) {
      p_eb[large_t, i] <- 2 * (
        -pt(t_stat[large_t, i], df = df[large_t], log.p = TRUE)
      )
    }
  }
  
  if (is.null(covariate)) {
    colnames(p_eb) <- paste0("subtype-", 1:n_sub, " vs. normal")
  } else {
    colnames(p_eb) <- c(
      paste0("subtype-", 1:n_sub, " vs. normal"),
      paste0("covariate-", seq_len(ncol(covariate)))
    )
  }
  d$p.contamDE.robust <- p_eb
  
  ############################################################
  d$proportion <- w_hat
  d$design <- design
  d$log2FC <- log2_fc
  d$log2_fc_cov <- log2_fc_cov
  d$df <- df
  d$y <- y
  d$weight <- weight
  return(list(
    counts = d$counts,
    p.contamDE.lm = d$p.contamDE.robust,
    log2FC = d$log2FC,
    log2_fc_cov = d$log2_fc_cov,
    proportion = d$proportion,
    design = d$design,
    df = d$df,
    weight = d$weight,
    y = d$y
  ))
}