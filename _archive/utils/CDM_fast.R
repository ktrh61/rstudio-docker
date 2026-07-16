# CDM_fast.R — 評価・改善版
# 高次元低サンプルサイズ (HDLSS) データ向け Cross-Data Matrix 主成分分析
# 
# 主な改善点:
# - エラーハンドリング強化
# - 数値安定性の向上
# - 進捗表示の追加
# - 結果検証機能

# Input:
#   X: 数値行列 (features x samples) もしくは (samples x features)。デフォルトは genes x samples を想定
#   by_sample: TRUE のとき X を (samples x features) とみなす
#   center, scale.: 事前標準化オプション（scale. は列スケール）
#   split: 長さ2の整数ベクトル。前半/後半のサンプル数（例 c(n1, n2)）。NULLなら floor(n/2)
#   k: 返す主成分数（NULLで可能な最大）
#   return_scores: スコア（PC座標）も返す
#   verbose: 進捗表示
#
# Output (list):
#   values: 固有値（CDM 由来の擬似固有値）
#   vectors: 固有ベクトル（ロードings）
#   scores: スコア（必要時）
#   split_idx: 使った分割インデックス
#   centered, scaled: 前処理フラグ
#   variance_explained: 寄与率
#   cumulative_variance: 累積寄与率

CDM_fast <- function(X,
                     by_sample = FALSE,
                     center = TRUE,
                     scale. = FALSE,
                     split = NULL,
                     k = NULL,
                     return_scores = TRUE,
                     verbose = TRUE) {
  
  if (verbose) cat("Starting CDM_fast analysis...\n")
  
  # 入力検証の強化
  if (!is.numeric(X)) stop("X must be numeric matrix or data.frame")
  if (!is.matrix(X)) X <- as.matrix(X)
  
  # NA/Inf チェック
  if (any(!is.finite(X))) {
    na_count <- sum(!is.finite(X))
    warning("Found ", na_count, " non-finite values in X. These will be problematic.")
    if (na_count / length(X) > 0.1) {
      stop("Too many non-finite values (>10%). Please clean the data first.")
    }
  }
  
  # 形状統一: genes x samples にする
  if (by_sample) X <- t(X)  # -> genes x samples
  p <- nrow(X); n <- ncol(X)
  
  if (verbose) cat("Data dimensions: ", p, " features x ", n, " samples\n")
  
  if (n < 2) stop("Need at least 2 samples.")
  if (p < 1) stop("Need at least 1 feature.")
  
  # HDLSS 警告
  if (p > n && verbose) {
    cat("High-dimensional data detected (p=", p, " > n=", n, "). CDM is appropriate.\n")
  }
  
  # 前処理
  original_X <- X  # 元データの保存（必要時）
  
  if (center) {
    if (verbose) cat("Centering features...\n")
    row_means <- rowMeans(X, na.rm = TRUE)
    X <- X - row_means
  }
  
  if (isTRUE(scale.)) {
    if (verbose) cat("Scaling features...\n")
    # 改善: より安全な標準化
    s <- apply(X, 1, function(x) sqrt(sum(x^2, na.rm = TRUE) / (length(x) - 1)))
    zero_var <- s < .Machine$double.eps
    if (any(zero_var)) {
      if (verbose) cat("Warning: ", sum(zero_var), " features have zero variance\n")
      s[zero_var] <- 1  # ゼロ分散の特徴量は標準化しない
    }
    X <- X / s
  }
  
  # サンプル分割の改善
  if (is.null(split)) {
    n1 <- floor(n/2); n2 <- n - n1
    if (verbose) cat("Auto-split: n1=", n1, ", n2=", n2, "\n")
  } else {
    if (length(split) != 2 || !is.numeric(split)) {
      stop("split must be numeric vector of length 2: c(n1, n2)")
    }
    n1 <- as.integer(split[1]); n2 <- as.integer(split[2])
    if (n1 <= 0 || n2 <= 0) stop("Both split values must be positive")
    if (n1 + n2 != n) stop("sum(split) = ", n1 + n2, " must equal number of samples = ", n)
    if (verbose) cat("Manual split: n1=", n1, ", n2=", n2, "\n")
  }
  
  # 最小サンプル数チェック
  if (n1 < 2 || n2 < 2) {
    warning("Very small split sizes (n1=", n1, ", n2=", n2, "). Results may be unstable.")
  }
  
  idx <- seq_len(n)
  idx1 <- idx[1:n1]
  idx2 <- idx[(n1+1):n]
  
  X1 <- X[, idx1, drop = FALSE]
  X2 <- X[, idx2, drop = FALSE]
  
  if (verbose) cat("Computing Cross-Data Matrix (CDM)...\n")
  
  # Cross-Data Matrix (CDM) 計算の改善
  # 数値安定性のための改善
  tryCatch({
    # 正規化係数の計算を安全に
    norm_factor <- sqrt(max(n1-1, 1) * max(n2-1, 1))
    if (norm_factor < .Machine$double.eps) {
      warning("Extremely small normalization factor. Results may be unreliable.")
      norm_factor <- 1
    }
    
    C <- crossprod(X1, X2) / norm_factor
    
    if (verbose) cat("CDM dimensions: ", nrow(C), " x ", ncol(C), "\n")
    
    # CDM の品質チェック
    if (any(!is.finite(C))) {
      stop("Non-finite values in Cross-Data Matrix. Check input data quality.")
    }
    
  }, error = function(e) {
    stop("Error computing Cross-Data Matrix: ", e$message)
  })
  
  # 特異値分解の改善
  if (verbose) cat("Performing SVD...\n")
  
  tryCatch({
    sv <- svd(C)
    d <- sv$d
    u <- sv$u  # in sample space of X1 (n1 x r)
    v <- sv$v  # in sample space of X2 (n2 x r)
    
    if (verbose) {
      cat("SVD completed. ", length(d), " singular values computed\n")
      cat("Singular value range: ", sprintf("%.2e", min(d)), " to ", sprintf("%.2e", max(d)), "\n")
    }
    
    # 数値的ランクの推定
    tol <- max(dim(C)) * max(d) * .Machine$double.eps
    num_rank <- sum(d > tol)
    if (verbose && num_rank < length(d)) {
      cat("Numerical rank: ", num_rank, " (out of ", length(d), ")\n")
    }
    
  }, error = function(e) {
    stop("Error in SVD computation: ", e$message)
  })
  
  # 特徴空間のロードings（近似）を再構成 - 改善版
  if (verbose) cat("Reconstructing loadings in feature space...\n")
  
  tryCatch({
    eps <- .Machine$double.eps
    d_safe <- pmax(d, eps)
    
    # メモリ効率的な計算
    W1 <- X1 %*% sweep(u, 2, d_safe, "/")
    W2 <- X2 %*% sweep(v, 2, d_safe, "/")
    loadings <- (W1 + W2) / 2  # p x r
    
    if (verbose) cat("Loadings dimensions: ", nrow(loadings), " x ", ncol(loadings), "\n")
    
  }, error = function(e) {
    stop("Error reconstructing loadings: ", e$message)
  })
  
  # 直交化（Gram-Schmidt）- 改善版
  if (verbose) cat("Orthogonalizing loadings...\n")
  
  tryCatch({
    qrW <- qr(loadings, tol = .Machine$double.eps)
    vectors <- qr.Q(qrW)  # p x r_ortho
    r_eff <- ncol(vectors)
    
    if (verbose) cat("Effective rank after orthogonalization: ", r_eff, "\n")
    
    # 直交性の確認
    if (r_eff > 1) {
      ortho_check <- max(abs(crossprod(vectors) - diag(r_eff)))
      if (ortho_check > 1e-10 && verbose) {
        cat("Warning: Orthogonality check failed (max deviation: ", sprintf("%.2e", ortho_check), ")\n")
      }
    }
    
  }, error = function(e) {
    stop("Error in orthogonalization: ", e$message)
  })
  
  # 固有値相当の計算と寄与率
  values <- d[seq_len(r_eff)]^2
  
  # 寄与率の計算
  total_var <- sum(values)
  if (total_var > 0) {
    variance_explained <- values / total_var
    cumulative_variance <- cumsum(variance_explained)
  } else {
    variance_explained <- rep(0, length(values))
    cumulative_variance <- rep(0, length(values))
  }
  
  # 返す主成分数を制限
  if (!is.null(k)) {
    if (!is.numeric(k) || k < 1) stop("k must be a positive integer")
    k <- min(as.integer(k), r_eff)
    vectors <- vectors[, seq_len(k), drop = FALSE]
    values  <- values[seq_len(k)]
    variance_explained <- variance_explained[seq_len(k)]
    cumulative_variance <- cumulative_variance[seq_len(k)]
    if (verbose) cat("Returning first ", k, " components\n")
  } else {
    k <- r_eff
    if (verbose) cat("Returning all ", k, " components\n")
  }
  
  # 結果オブジェクトの構築
  out <- list(
    values  = values,
    vectors = vectors,
    split_idx = list(idx1 = idx1, idx2 = idx2),
    centered = center,
    scaled = isTRUE(scale.),
    variance_explained = variance_explained,
    cumulative_variance = cumulative_variance,
    n_features = p,
    n_samples = n,
    n_components = k
  )
  
  # スコアの計算
  if (isTRUE(return_scores)) {
    if (verbose) cat("Computing scores...\n")
    tryCatch({
      # scores = t(vectors) %*% X (全サンプル)
      out$scores <- crossprod(vectors, X)  # k x n
      # 転置してサンプル x 成分にする
      out$scores <- t(out$scores)  # n x k
      rownames(out$scores) <- colnames(original_X)
      colnames(out$scores) <- paste0("PC", seq_len(k))
    }, error = function(e) {
      warning("Error computing scores: ", e$message)
      out$scores <- NULL
    })
  }
  
  # クラス設定
  class(out) <- c("cdm_fast", class(out))
  
  if (verbose) {
    cat("CDM analysis completed successfully!\n")
    if (k <= 5) {
      cat("First", min(k, 5), "components explain", 
          sprintf("%.1f%%", cumulative_variance[min(k, 5)] * 100), "of variance\n")
    }
  }
  
  return(out)
}

# プロット関数の追加
plot.cdm_fast <- function(x, components = c(1, 2), 
                          labels = NULL, 
                          title = "CDM Analysis", ...) {
  if (is.null(x$scores)) {
    stop("Scores not available. Re-run with return_scores = TRUE")
  }
  
  if (max(components) > ncol(x$scores)) {
    stop("Requested components exceed available components")
  }
  
  pc1 <- components[1]
  pc2 <- components[2]
  
  plot(x$scores[, pc1], x$scores[, pc2],
       xlab = paste0("PC", pc1, " (", sprintf("%.1f%%", x$variance_explained[pc1] * 100), ")"),
       ylab = paste0("PC", pc2, " (", sprintf("%.1f%%", x$variance_explained[pc2] * 100), ")"),
       main = title, ...)
  
  if (!is.null(labels)) {
    text(x$scores[, pc1], x$scores[, pc2], labels, pos = 3, cex = 0.8)
  }
  
  invisible(x)
}

# 使用例とテスト関数
test_CDM_fast <- function() {
  cat("Testing CDM_fast with simulated data...\n")
  
  # テストデータ生成
  set.seed(123)
  p <- 1000  # 遺伝子数
  n <- 20    # サンプル数
  X <- matrix(rnorm(p * n), nrow = p, ncol = n)
  
  # テスト実行
  result <- CDM_fast(X, verbose = TRUE)
  
  cat("\nTest Results:\n")
  cat("- Dimensions: ", result$n_features, "x", result$n_samples, "\n")
  cat("- Components: ", result$n_components, "\n")
  cat("- First 3 PCs explain: ", sprintf("%.1f%%", result$cumulative_variance[3] * 100), "\n")
  
  # プロット（簡単な確認）
  if (!is.null(result$scores)) {
    plot(result)
  }
  
  return(result)
}

cat("CDM_fast function loaded successfully.\n")
cat("Run test_CDM_fast() to verify functionality.\n")