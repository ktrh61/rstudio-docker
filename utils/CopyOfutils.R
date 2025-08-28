# utils.R — helpers for MUREN (revised, drop-in)

# ---- deps (軽いチェック) ----
if (!requireNamespace("matrixStats", quietly = TRUE))
  stop("Package 'matrixStats' is required for filter_gene_l().")
if (!requireNamespace("MASS", quietly = TRUE))
  stop("Package 'MASS' is required for robust regressions.")
# robustbase は任意（あれば優先）

# ---- constants ----
TOL <- .Machine$double.eps ^ 0.5

# ---- small helpers ----
is.wholenumber <- function(x, tol = TOL) abs(x - round(x)) < tol
lg <- function(x) log2(1 + x)          # log-space used inside MUREN
ep <- function(x) 2^x - 1              # inverse of lg

# ---- gene filter (light, group-agnostic) ----
# Keep genes with max count >= trim AND present (>0) in at least 2 samples
filter_gene_l <- function(reads, trim) {
  matrixStats::rowMaxs(reads) >= trim & rowSums(reads > 0) >= 2
}

# ---- robust regression backend (C if available) ----
# Prefer robustbase::ltsReg (fast C implementation); fallback to MASS::ltsreg.
.reg_backend <- function(formula, ...) {
  if (requireNamespace("robustbase", quietly = TRUE)) {
    # robustbase::ltsReg: ~1 で切片のみ、~ x で傾き+切片。intercept は式に従うので明示不要
    robustbase::ltsReg(formula, ...)
  } else {
    MASS::ltsreg(formula, ...)
  }
}

# ---- single-parameter regression (location on log space) ----
# Returns a scalar coefficient per (sample, ref) pair
# method は options("muren_pair_method") で切替（muren_norm がセットする）
reg_sp <- function(s_k, s_r, ...) {
  y <- s_r - s_k
  method <- getOption("muren_pair_method", "lts")
  
  if (method == "median") {
    return(stats::median(y, na.rm = TRUE))
  }
  if (method == "trim10") {
    return(mean(y, trim = 0.10, na.rm = TRUE))
  }
  if (method == "huber") {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      warning("MASS::rlm not available; falling back to median()")
      return(stats::median(y, na.rm = TRUE))
    }
    fit <- MASS::rlm(y ~ 1, psi = MASS::psi.huber, maxit = 20)
    return(unname(coef(fit)[1]))
  }
  # default: LTS（従来互換）
  .reg_backend(y ~ 1, ...)$coefficients[1]
}

# ---- mode-based shift (alternative robust location) ----
# Returns a scalar (mode of differences)
mode_sp <- function(s_k, s_r, ...) {
  d <- stats::density(s_r - s_k, ...)
  d$x[which.max(d$y)]
}

# ---- double-parameter regression (non-linear/power) ----
# Returns fitted values (vector length = n_gene)
reg_dp <- function(s_k, s_r, ...) {
  .reg_backend(s_r ~ s_k, ...)$fitted.values
}

# ---- legacy task builder (kept for backward compatibility) ----
# Not used by the revised muren_norm(), but safe to keep.
get_tasks <- function(k, reg_wapper, refs) {
  parse(text = paste(
    reg_wapper,
    "(log_raw_reads_mx[, ", k, "],",
    "log_raw_reads_mx[, ", refs, "], ...)",
    collapse = "\n", sep = ""
  ))
}

# ---- median polish: sample effects from pairwise results (single-param) ----
# fitted_n: numeric vector of length = number of (sample,ref) pairs
# locations maps vector back to a (ref x sample) matrix
polish_coeff <- function(fitted_n, n_exp, locations, unused_refs, maxiter) {
  rs_mx <- matrix(NA_real_, nrow = n_exp, ncol = n_exp)
  rs_mx[locations] <- fitted_n
  if (length(unused_refs) > 0) rs_mx <- rs_mx[-unused_refs, , drop = FALSE]
  m <- stats::medpolish(rs_mx, na.rm = TRUE, trace.iter = FALSE, maxiter = maxiter)
  # column effects + overall (sample effects in log2 space)
  m$overall + m$col
}

# ---- median polish per gene (double-param path) ----
polish_one_gene <- function(fitted_n, n_exp, locations, unused_refs, maxiter) {
  rs_mx <- matrix(NA_real_, nrow = n_exp, ncol = n_exp)
  rs_mx[locations] <- fitted_n
  if (length(unused_refs) > 0) rs_mx <- rs_mx[-unused_refs, , drop = FALSE]
  m <- stats::medpolish(rs_mx, na.rm = TRUE, trace.iter = FALSE, maxiter = maxiter)
  m$overall + m$col
}
