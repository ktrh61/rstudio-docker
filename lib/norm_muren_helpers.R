# ==============================================================================
# MUREN helpers
# Statistical model follows the original MUREN implementation. The LTS engine
# is MASS::ltsreg; no automatic estimator selection or fallback is provided.
# ==============================================================================

if (!requireNamespace("matrixStats", quietly = TRUE)) {
  stop("Package 'matrixStats' is required for filter_gene_l().")
}
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required for LTS regression.")
}

TOL <- .Machine$double.eps^0.5

is.wholenumber <- function(x, tol = TOL) {
  abs(x - round(x)) < tol
}

lg <- function(x) {
  log2(1 + x)
}

ep <- function(x) {
  2^x - 1
}

filter_gene_l <- function(reads, trim) {
  matrixStats::rowMaxs(reads) >= trim
}

# Single-parameter MUREN: robust location shift between a target and reference.
reg_sp <- function(s_k, s_r, ...) {
  MASS::ltsreg(s_r - s_k ~ 1, ...)$coefficients[[1L]]
}

# Original mode-based alternative for the single-parameter form.
mode_sp <- function(s_k, s_r, ...) {
  d <- stats::density(s_r - s_k, ...)
  d$x[which.max(d$y)]
}

# Double-parameter MUREN: fitted reference log-counts for each gene.
reg_dp <- function(s_k, s_r, ...) {
  MASS::ltsreg(s_r ~ s_k, ...)$fitted.values
}

polish_coeff <- function(fitted_n, n_exp, locations, unused_refs, maxiter) {
  rs_mx <- matrix(NA_real_, nrow = n_exp, ncol = n_exp)
  rs_mx[locations] <- fitted_n

  if (length(unused_refs) > 0L) {
    rs_mx <- rs_mx[-unused_refs, , drop = FALSE]
  }

  fit <- stats::medpolish(
    rs_mx,
    na.rm = TRUE,
    trace.iter = FALSE,
    maxiter = maxiter
  )

  fit$overall + fit$col
}

polish_one_gene <- function(fitted_n, n_exp, locations, unused_refs, maxiter) {
  rs_mx <- matrix(NA_real_, nrow = n_exp, ncol = n_exp)
  rs_mx[locations] <- fitted_n

  if (length(unused_refs) > 0L) {
    rs_mx <- rs_mx[-unused_refs, , drop = FALSE]
  }

  fit <- stats::medpolish(
    rs_mx,
    na.rm = TRUE,
    trace.iter = FALSE,
    maxiter = maxiter
  )

  fit$overall + fit$col
}
