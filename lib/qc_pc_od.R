# PC-OD: PC-scores-based outlier detection (Nakayama, Yata & Aoshima 2024,
# Jpn J Stat Data Sci 7:739-766, doi:10.1007/s42081-024-00255-0). Iterated
# Grubbs test (Grubbs 1969, Technometrics 11:1-21) on the first PC's dual
# vector; in-house implementation following the paper's algorithm.
PC_OD <- function(X, alpha = 0.05) {
  n_all <- ncol(X)
  samples <- seq_len(n_all)
  inliers <- samples
  d <- nrow(X) # 反復中は列のみ削るので d は不変

  while (TRUE) {
    n <- length(inliers)
    X_cent <- X - rowMeans(X)

    # 第1右特異ベクトル dualvec のみを、小さい方のグラム行列の
    # 固有分解から求める（フルSVDを避ける）
    if (d <= n) {
      G <- tcrossprod(X_cent) # A A^T (d×d)
      u <- eigen(G, symmetric = TRUE)$vectors[, 1]
      w <- crossprod(X_cent, u) # A^T u = sigma1 * v1
      dualvec <- w / sqrt(sum(w^2)) # 単位ノルムに正規化
    } else {
      G <- crossprod(X_cent) # A^T A (n×n)
      dualvec <- eigen(G, symmetric = TRUE)$vectors[, 1]
    }

    # Step3
    t_value <- qt(alpha / (2 * n), n - 2)
    t0 <- ((n - 1) / sqrt(n)) * sqrt(t_value^2 / (n - 2 + t_value^2))
    T_PC1 <- max(abs(dualvec)) * sqrt(n - 1)

    # Step4
    if (T_PC1 <= t0) {
      break
    } else {
      outlier <- which.max(abs(dualvec))
      X <- X[, -outlier]
      inliers <- inliers[-outlier]
    }
  }

  setdiff(samples, inliers)
}
