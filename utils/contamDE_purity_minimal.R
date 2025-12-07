# ==============================================================================
# contamDE_purity_minimal.R
# Minimal contamDE purity: TMM→MUREN, BH→qvalue, DEG部分は削除
# 依存: edgeR, limma, qvalue, (utils/norm_improved.R: muren_norm)
# ==============================================================================

# --- voom with MUREN scaling (TMMの置換だけ) ---
limma_voom_muren <- function(counts,
                             pairwise_method = "lts",
                             workers = "auto") {
  counts <- as.matrix(counts)
  N <- ncol(counts) / 2
  stopifnot(N == floor(ncol(counts) / 2))
  
  # DGEList
  d <- edgeR::DGEList(counts = counts)
  
  # MUREN scaling coefficients (single_param, scaling_coeff)
  nf <- muren_norm(
    counts,
    refs = "saturated",
    pairwise_method = pairwise_method,
    single_param = TRUE,
    res_return = "scaling_coeff",
    workers = workers
  )
  
  # 係数の検証 + 幾何平均=1に再スケール
  if (any(!is.finite(nf)) || any(nf <= 0)) {
    stop("MUREN scaling coefficients contain non-finite or non-positive values.")
  }
  nf <- nf / exp(mean(log(nf)))
  d$samples$norm.factors <- as.numeric(nf)
  
  # デザイン（Normal/Tumor）
  condition <- factor(c(rep("Normal", N), rep("Tumor", N)))
  design <- stats::model.matrix(~ 0 + condition)
  colnames(design) <- levels(condition)
  
  # voom（サイズ補正済みなので normalize.method = "none"）
  v <- limma::voom(d, design, normalize.method = "none")
  
  # 対比: Tumor - Normal
  fit  <- limma::lmFit(v, design)
  ctr  <- limma::makeContrasts(contrasts = "Tumor - Normal", levels = design)
  fit2 <- limma::contrasts.fit(fit, ctr)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  
  # size（有効ライブラリサイズの相対スケール）
  size <- d$samples$lib.size * d$samples$norm.factors
  size <- size / mean(size)
  
  list(
    p.limma       = fit2$p.value[, 1],
    log2FC.limma  = fit2$coefficients[, 1],
    size          = size
  )
}

# --- contamDE 純度推定（DEG部は削除、qvalueでinformative選定） ---
contamDE_purity_minimal <- function(counts,
                                    subtype        = NULL,
                                    covariate      = NULL,
                                    is.contaminated = TRUE,
                                    pairwise_method = "lts",
                                    workers         = "auto",
                                    prior.count     = 1,      # 元ネタ(+1)準拠
                                    verbose         = TRUE) {
  if (verbose) cat("Starting tumor purity estimation (minimal)...\n")
  
  counts <- as.matrix(counts)
  N <- ncol(counts) / 2
  G <- nrow(counts)
  stopifnot(N == floor(ncol(counts) / 2))
  
  # 初期ランク付け: limma+voom (MUREN)
  d <- limma_voom_muren(counts, pairwise_method = pairwise_method, workers = workers)
  p.limma      <- d$p.limma
  log2FC.limma <- d$log2FC.limma
  size         <- d$size
  
  # デザイン（元ネタ準拠）
  if (is.null(subtype) || length(unique(subtype)) == 1) {
    subtype <- rep(1, N)
    if (is.null(covariate)) {
      design <- matrix(subtype, ncol = 1)
    } else {
      covariate <- matrix(covariate, ncol = N)
      design <- stats::model.matrix(~ covariate)
    }
  } else {
    subtype <- factor(subtype)
    if (is.null(covariate)) {
      design <- stats::model.matrix(~ 0 + subtype)
    } else {
      covariate0 <- matrix(covariate, nrow = N)
      design <- stats::model.matrix(~ 0 + subtype + covariate0)
    }
  }
  K <- ncol(design)
  
  # MURENサイズで正規化 → N/T に分割（擬似カウントは +prior.count）
  count.norm  <- t(t(counts) / size)
  count.normal <- count.norm[,  (1:N)] + prior.count
  count.tumor  <- count.norm[, -(1:N)] + prior.count
  
  # y = log2(T) - log2(N)
  y <- log2(count.tumor) - log2(count.normal)
  
  # 純度推定
  w_hat <- rep(1, N)
  up <- down <- integer(0)
  
  if (is.contaminated) {
    # 多重補正: BH→qvalue に置換（上限1000はオリジナル踏襲）
    p_adj <- tryCatch(
      qvalue::qvalue(p.limma, pi0.method = "bootstrap")$qvalues,
      error = function(e) stop("qvalue calculation failed: ", e$message)
    )
    if (sum(p_adj < 0.1, na.rm = TRUE) > 1000L) {
      idx_sig    <- which(is.finite(p_adj) & p_adj < 0.1)
      idx_sorted <- idx_sig[order(p_adj[idx_sig])]
      keep       <- idx_sorted[seq_len(1000L)]
      p_adj[-keep] <- 1
    }
    
    # informative 遺伝子（元ネタ閾値）
    up   <- which(p_adj < 0.1 & log2FC.limma >  log2(1.5))
    down <- which(p_adj < 0.1 & log2FC.limma < -log2(1.5))
    
    y_up   <- if (length(up))   y[up,  , drop = FALSE] else matrix(0, 0, N)
    y_down <- if (length(down)) y[down,, drop = FALSE] else matrix(0, 0, N)
    
    sumup   <- colSums(y_up)
    sumdown <- colSums(y_down)
    sum.max <- max(sumup - sumdown)
    
    if (!is.finite(sum.max) || sum.max <= 0) {
      w_hat <- rep(1.0, N)
    } else {
      w_hat <- (sumup - sumdown) / sum.max
      w_hat <- pmax(0, pmin(1, w_hat))
    }
  }
  
  if (verbose) {
    cat(sprintf("Informative genes: up=%d, down=%d\n", length(up), length(down)))
    cat(sprintf("Purity summary: mean=%.3f, sd=%.3f\n", mean(w_hat), stats::sd(w_hat)))
  }
  
  # 返り値（DEG結果は含めない）
  list(
    proportion          = w_hat,
    size                = size,
    design              = design,
    y                   = y,
    n_pairs             = N,
    n_genes             = G,
    informative_genes   = list(up = up, down = down),
    p_init              = p.limma,
    log2FC_init         = log2FC.limma,
    normalization_method = paste0("MUREN_", toupper(pairwise_method)),
    estimation_date     = Sys.time()
  )
}