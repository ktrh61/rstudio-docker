# analysis_v7/utils/robust_pca_utils.R
# Robust PCA utilities for CDM-based outlier detection with adaptive optimization

# ============================================================================
# Thread Control
# ============================================================================
setup_threads <- function() {
  if(requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    message("BLAS/OMP threads set to 1 via RhpcBLASctl")
  } else {
    Sys.setenv(
      OMP_NUM_THREADS = 1,
      OPENBLAS_NUM_THREADS = 1,
      MKL_NUM_THREADS = 1,
      BLIS_NUM_THREADS = 1,
      VECLIB_MAXIMUM_THREADS = 1,
      ARMA_OPENMP_THREADS = 1
    )
    message("Thread control via environment variables")
  }
}

# ============================================================================
# Preprocessing with robust scaling
# ============================================================================
rprep <- function(X, rel_eps = 1e-8) {
  med <- apply(X, 2, median, na.rm = TRUE)
  madv <- apply(X, 2, mad, constant = 1.4826, na.rm = TRUE)
  iqr <- apply(X, 2, IQR, na.rm = TRUE) / 1.349
  
  scalev <- madv
  small <- scalev < rel_eps
  scalev[small] <- iqr[small]
  
  eps <- rel_eps * max(1, median(abs(X), na.rm = TRUE))
  keep <- scalev > eps
  
  if(sum(!keep) > 0) {
    cols_removed <- paste(which(!keep), collapse = ", ")
    warning(sprintf("Removing %d near-constant genes (columns %s)", 
                    sum(!keep), cols_removed))
    X <- X[, keep, drop = FALSE]
    med <- med[keep]
    scalev <- scalev[keep]
  }
  
  sweep(sweep(X, 2, med, "-"), 2, scalev, "/")
}

# ============================================================================
# Permutation Parallel Analysis with controlled parallelism
# ============================================================================
permPA <- function(X, R = 300, alpha = 0.05, maxK = 8, seed = NULL, mc.cores = 1) {
  # Simple seed setting (RNGkind already set in main script)
  if(!is.null(seed)) set.seed(seed)
  
  # Real eigenvalues
  s <- svd(X, nu = 0, nv = 0)$d^2
  n <- nrow(X)
  p <- ncol(X)
  L <- min(min(n, p), maxK)
  
  # Worker function
  worker <- function(r) {
    X0 <- apply(X, 2, sample)
    sr <- svd(X0, nu = 0, nv = 0)$d^2
    head(c(sr, rep(0, L)), L)
  }
  
  # Parallel or serial execution
  if(mc.cores > 1 && R >= 100) {
    library(parallel)
    chunks <- split(seq_len(R), cut(seq_len(R), mc.cores, labels = FALSE))
    null_eigs <- mclapply(chunks, function(ch) {
      do.call(rbind, lapply(ch, worker))
    }, mc.cores = mc.cores, mc.preschedule = TRUE, mc.set.seed = TRUE)
    null_eigs <- do.call(rbind, null_eigs)
  } else {
    null_eigs <- do.call(rbind, lapply(seq_len(R), worker))
  }
  
  # Calculate threshold
  thr <- apply(null_eigs, 2, quantile, probs = 1 - alpha)
  K0 <- sum(head(s, L) > thr)
  
  list(
    K0 = K0, 
    thr = thr, 
    eigs = head(s, L), 
    null_eigs = null_eigs
  )
}

# ============================================================================
# Orthogonalization via QR
# ============================================================================
orth <- function(M) {
  if(ncol(M) == 0 || nrow(M) == 0) return(M)
  qr.Q(qr(M))
}

# ============================================================================
# Main robust PCA function with adaptive B
# ============================================================================
cdm_robust_pca <- function(logcpm_group,
                           alpha = 0.05,
                           R = 300,
                           B_init = 50,
                           B_step = 25,
                           B_max = 200,
                           frac = 0.8,
                           maxK = 8,
                           seed = 1,
                           verbose = FALSE,
                           mc.cores_pa = 1) {
  
  set.seed(seed)
  
  # Preprocessing
  X <- rprep(logcpm_group)
  
  if(verbose) {
    cat("Data after preprocessing:", nrow(X), "x", ncol(X), "\n")
    cat("Permutation PA: R =", R, ", mc.cores =", mc.cores_pa, "\n")
  }
  
  # Permutation Parallel Analysis
  pa <- permPA(X, R = R, alpha = alpha, maxK = maxK, 
               seed = seed + 1, mc.cores = mc.cores_pa)
  K <- max(1, min(pa$K0, nrow(X) - 1, ncol(X), maxK))
  
  if(verbose) {
    cat("PA selected K =", pa$K0, ", using K =", K, "\n")
  }
  
  # Full CDM fit
  fit_full <- CDM_fast3_arma(X, verbose = FALSE)
  V0 <- fit_full$vectors
  
  # Ensure row names are set for V0
  if(is.null(rownames(V0))) {
    rownames(V0) <- colnames(X)
  }
  
  # Adaptive stability evaluation with gene alignment
  stability_adaptive <- function(K_test) {
    acc <- numeric(0)
    n <- nrow(X)
    
    # Prepare full basis with row names
    V0_subset <- V0[, 1:min(K_test, ncol(V0)), drop = FALSE]
    U0 <- orth(V0_subset)
    
    B_current <- 0
    B_next <- B_init
    
    while(B_current < B_max) {
      B_todo <- min(B_next, B_max - B_current)
      cc <- numeric(B_todo)
      
      for(b in seq_len(B_todo)) {
        idx <- sample.int(n, size = ceiling(frac * n), replace = FALSE)
        X_sub <- X[idx, , drop = FALSE]
        
        # Subsample CDM
        fit_b <- CDM_fast3_arma(X_sub, verbose = FALSE)
        Vb <- fit_b$vectors
        
        # Set row names for Vb
        if(is.null(rownames(Vb))) {
          rownames(Vb) <- colnames(X_sub)
        }
        
        # Find common genes for alignment
        common_genes <- intersect(rownames(U0), rownames(Vb))
        
        # Need minimum number of common genes
        if(length(common_genes) < 10) {
          cc[b] <- 0
          next
        }
        
        K_use <- min(K_test, ncol(Vb))
        if(K_use > 0) {
          # Extract and orthogonalize common subsets
          U0_common <- orth(U0[common_genes, , drop = FALSE])
          Vb_subset <- Vb[common_genes, 1:K_use, drop = FALSE]
          Ub_common <- orth(Vb_subset)
          
          K_fin <- min(ncol(U0_common), ncol(Ub_common))
          if(K_fin > 0) {
            # Compute principal angles via SVD
            S <- svd(t(U0_common[, 1:K_fin, drop = FALSE]) %*% 
                       Ub_common[, 1:K_fin, drop = FALSE])$d
            cc[b] <- min(S)
          } else {
            cc[b] <- 0
          }
        } else {
          cc[b] <- 0
        }
      }
      
      acc <- c(acc, cc)
      B_current <- length(acc)
      
      # Early termination check (after 50 samples)
      if(B_current >= 50) {
        iqr_val <- IQR(acc, na.rm = TRUE)
        if(iqr_val < 0.05) {
          if(verbose) {
            cat("  Early stop at B =", B_current, 
                "(IQR =", sprintf("%.4f", iqr_val), "< 0.05)\n")
          }
          break
        }
      }
      
      # Progress report
      if(verbose && B_current == B_init) {
        cat("  B =", B_current, ", median stability =", 
            sprintf("%.3f", median(acc)), ", IQR =", 
            sprintf("%.4f", IQR(acc, na.rm = TRUE)), "\n")
      }
      
      B_next <- B_step
    }
    
    if(verbose && B_current >= B_max) {
      cat("  Reached B_max =", B_max, "\n")
    }
    
    acc
  }
  
  # K adjustment loop
  cc <- stability_adaptive(K)
  med_cc <- median(cc, na.rm = TRUE)
  
  if(verbose) {
    cat("K =", K, ", median stability =", sprintf("%.3f", med_cc), 
        ", B =", length(cc), "\n")
  }
  
  # Adjust K if stability is insufficient
  while(med_cc < 0.8 && K > 1) {
    K <- K - 1
    cc <- stability_adaptive(K)
    med_cc <- median(cc, na.rm = TRUE)
    if(verbose) {
      cat("K =", K, ", median stability =", sprintf("%.3f", med_cc), 
          ", B =", length(cc), "\n")
    }
  }
  
  # Return comprehensive results
  list(
    K = K,
    cdm = fit_full,
    pa = pa,
    stability_cc = cc,
    final_median_cc = med_cc,
    final_B = length(cc),
    convergence = list(
      stable = med_cc >= 0.8,
      IQR = IQR(cc, na.rm = TRUE),
      B_used = length(cc)
    )
  )
}

# ============================================================================
# Utility function for progress reporting
# ============================================================================
report_progress <- function(group_name, result) {
  cat(sprintf("\n%s: n = %d, K = %d, stability = %.3f, B = %d\n",
              group_name,
              result$cdm$scores %>% nrow(),
              result$K,
              result$final_median_cc,
              result$final_B))
  
  if(result$convergence$IQR < 0.05) {
    cat("  (Early convergence achieved)\n")
  }
}