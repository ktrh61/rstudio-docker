# ==============================================================================
# utils.R — Improved helpers for MUREN normalization
# Version: Enhanced performance and robustness
# ==============================================================================

# ---- Dependencies Check ----
if (!requireNamespace("matrixStats", quietly = TRUE))
  stop("Package 'matrixStats' is required for filter_gene_l().")
if (!requireNamespace("MASS", quietly = TRUE))
  stop("Package 'MASS' is required for robust regressions.")
# robustbase is optional (preferred if available for speed)

# ---- Constants ----
TOL <- .Machine$double.eps ^ 0.5

# ---- Small Helper Functions ----
is.wholenumber <- function(x, tol = TOL) {
  abs(x - round(x)) < tol
}

lg <- function(x) {
  log2(1 + x)  # log-space transformation used inside MUREN
}

ep <- function(x) {
  2^x - 1      # inverse of lg transformation
}

# ---- Enhanced Gene Filtering ----
# Keep genes with max count >= trim AND present (>0) in at least 2 samples
# This prevents zero-variance genes from entering the normalization
filter_gene_l <- function(reads, trim) {
  max_counts <- matrixStats::rowMaxs(reads)
  min_samples <- rowSums(reads > 0)
  
  # Return logical vector for genes passing both criteria
  (max_counts >= trim) & (min_samples >= 2)
}

# ---- Robust Regression Backend Selection ----
# Prefer robustbase::ltsReg (fast C implementation) over MASS::ltsreg
.reg_backend <- function(formula, ...) {
  if (requireNamespace("robustbase", quietly = TRUE)) {
    # robustbase::ltsReg: fast C implementation
    # Formula handling: ~1 for intercept-only, ~x for slope+intercept
    robustbase::ltsReg(formula, ...)
  } else {
    # Fallback to MASS implementation
    MASS::ltsreg(formula, ...)
  }
}

# ---- Single-Parameter Regression (Location Shift) ----
# Returns a scalar coefficient per (sample, ref) pair
# Method can be controlled via options("muren_pair_method")
reg_sp <- function(s_k, s_r, ...) {
  y <- s_r - s_k
  method <- getOption("muren_pair_method", "lts")
  
  # Fast median-based location shift
  if (method == "median") {
    return(stats::median(y, na.rm = TRUE))
  }
  
  # Trimmed mean (removes 10% extreme values)
  if (method == "trim10") {
    return(mean(y, trim = 0.10, na.rm = TRUE))
  }
  
  # Huber robust regression (M-estimator)
  if (method == "huber") {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      warning("MASS::rlm not available for Huber regression; falling back to median")
      return(stats::median(y, na.rm = TRUE))
    }
    tryCatch({
      fit <- MASS::rlm(y ~ 1, psi = MASS::psi.huber, maxit = 20)
      return(unname(coef(fit)[1]))
    }, error = function(e) {
      warning("Huber regression failed; falling back to median: ", e$message)
      return(stats::median(y, na.rm = TRUE))
    })
  }
  
  # Default: LTS regression (backward compatible)
  tryCatch({
    .reg_backend(y ~ 1, ...)$coefficients[1]
  }, error = function(e) {
    warning("LTS regression failed; falling back to median: ", e$message)
    return(stats::median(y, na.rm = TRUE))
  })
}

# ---- Mode-Based Shift (Alternative Robust Location) ----
# Returns the mode of the difference distribution
mode_sp <- function(s_k, s_r, ...) {
  tryCatch({
    d <- stats::density(s_r - s_k, ...)
    d$x[which.max(d$y)]
  }, error = function(e) {
    warning("Density estimation failed; falling back to median: ", e$message)
    return(stats::median(s_r - s_k, na.rm = TRUE))
  })
}

# ---- Double-Parameter Regression (Non-linear Correction) ----
# Returns fitted values (vector of length = n_genes)
reg_dp <- function(s_k, s_r, ...) {
  tryCatch({
    .reg_backend(s_r ~ s_k, ...)$fitted.values
  }, error = function(e) {
    warning("Double-parameter regression failed; returning input: ", e$message)
    return(s_r)  # Return reference values if regression fails
  })
}

# ---- Legacy Task Builder (Backward Compatibility) ----
# Not used by the revised muren_norm(), but kept for compatibility
get_tasks <- function(k, reg_wapper, refs) {
  parse(text = paste(
    reg_wapper,
    "(log_raw_reads_mx[, ", k, "],",
    "log_raw_reads_mx[, ", refs, "], ...)",
    collapse = "\n", sep = ""
  ))
}

# ---- Median Polish: Sample Effects from Pairwise Results ----
# For single-parameter regression results
# fitted_n: numeric vector of pairwise regression results
# locations: indices mapping vector back to (ref x sample) matrix
polish_coeff <- function(fitted_n, n_exp, locations, unused_refs, maxiter) {
  # Initialize matrix for pairwise results
  rs_mx <- matrix(NA_real_, nrow = n_exp, ncol = n_exp)
  rs_mx[locations] <- fitted_n
  
  # Remove unused references (rows) if any
  if (length(unused_refs) > 0) {
    rs_mx <- rs_mx[-unused_refs, , drop = FALSE]
  }
  
  # Apply median polish to extract sample effects
  tryCatch({
    m <- stats::medpolish(rs_mx, 
                          na.rm = TRUE, 
                          trace.iter = FALSE, 
                          maxiter = maxiter)
    
    # Return column effects + overall (sample effects in log2 space)
    m$overall + m$col
  }, error = function(e) {
    warning("Median polish failed; returning zero effects: ", e$message)
    return(rep(0, ncol(rs_mx)))
  })
}

# ---- Median Polish Per Gene (Double-Parameter Path) ----
# For gene-specific normalization in double-parameter regression
polish_one_gene <- function(fitted_n, n_exp, locations, unused_refs, maxiter) {
  # Initialize matrix for this gene's pairwise results
  rs_mx <- matrix(NA_real_, nrow = n_exp, ncol = n_exp)
  rs_mx[locations] <- fitted_n
  
  # Remove unused references (rows) if any
  if (length(unused_refs) > 0) {
    rs_mx <- rs_mx[-unused_refs, , drop = FALSE]
  }
  
  # Apply median polish for this gene
  tryCatch({
    m <- stats::medpolish(rs_mx, 
                          na.rm = TRUE, 
                          trace.iter = FALSE, 
                          maxiter = maxiter)
    
    # Return column effects + overall (sample effects for this gene)
    m$overall + m$col
  }, error = function(e) {
    warning("Gene-specific median polish failed; returning zero effects: ", e$message)
    return(rep(0, ncol(rs_mx)))
  })
}

# ---- Performance Enhancement Functions ----

# Set optimal MUREN method based on data size and requirements
set_muren_method <- function(n_genes, n_samples, priority = "balanced") {
  total_comparisons <- n_genes * n_samples * (n_samples - 1)
  
  if (priority == "speed") {
    if (total_comparisons > 1e6) {
      options(muren_pair_method = "median")
      message("Large dataset detected: Using median method for maximum speed")
    } else {
      options(muren_pair_method = "trim10")
      message("Using trimmed mean method for good speed/robustness balance")
    }
  } else if (priority == "robustness") {
    options(muren_pair_method = "lts")
    message("Using LTS method for maximum robustness")
  } else {  # balanced
    if (total_comparisons > 5e5) {
      options(muren_pair_method = "trim10")
      message("Using trimmed mean method for balanced performance")
    } else {
      options(muren_pair_method = "lts")
      message("Using LTS method for optimal robustness")
    }
  }
}

# Check and report MUREN configuration
check_muren_config <- function() {
  method <- getOption("muren_pair_method", "lts")
  has_robustbase <- requireNamespace("robustbase", quietly = TRUE)
  has_mass <- requireNamespace("MASS", quietly = TRUE)
  
  cat("MUREN Configuration:\n")
  cat(sprintf("  Method: %s\n", method))
  cat(sprintf("  robustbase available: %s\n", has_robustbase))
  cat(sprintf("  MASS available: %s\n", has_mass))
  
  if (method %in% c("lts", "huber") && !has_mass && !has_robustbase) {
    warning("Neither robustbase nor MASS available for robust regression. Consider installing robustbase for better performance.")
  }
  
  invisible(list(method = method, robustbase = has_robustbase, mass = has_mass))
}

# ---- Initialization Function ----
# Call this to set up optimal MUREN configuration
initialize_muren <- function(n_genes = NULL, n_samples = NULL, priority = "balanced") {
  if (!is.null(n_genes) && !is.null(n_samples)) {
    set_muren_method(n_genes, n_samples, priority)
  }
  check_muren_config()
  invisible(TRUE)
}

# ---- End of utils.R ----
cat("Enhanced MUREN utils loaded successfully.\n")
cat("Use initialize_muren(n_genes, n_samples) to optimize configuration.\n")