# 07_pca_outlier_detection.R - Stage 1 PCA Outlier Detection (v7 Enhanced)
# Purpose: Remove technical outliers using CDM with PA-selected K
# Method: Adaptive thresholds with robust safeguards

source("analysis_v7/setup.R")

cat("\n=== Stage 1: PCA Outlier Detection (v7 Enhanced) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# Load packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(Rcpp)
  library(RcppArmadillo)
  library(parallel)
})

# Thread control
lock_single_thread <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
    message("BLAS/OMP set to 1 via RhpcBLASctl")
  } else {
    Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1")
    message("BLAS/OMP set to 1 via env vars")
  }
}
lock_single_thread()
RNGkind("L'Ecuyer-CMRG")

# Compile CDM
cat("\nCompiling CDM implementation...\n")
sourceCpp("utils/CDM_fast3_arma_enhanced.cpp")

# Configuration with adaptive thresholds
CONFIG <- list(
  PA_R = 200,
  PA_ALPHA = 0.05,
  PA_L = 8,
  PA_SEED = 123,
  # Adaptive outlier detection parameters
  MIN_OUTLIERS_PCT = 0.00,    # Minimum outlier percentage (0% = no forced outliers)
  MAX_OUTLIERS_PCT = 0.20,    # Maximum outlier percentage (20%)
  IQR_MULTIPLIER = 3.0,       # IQR multiplier for outlier detection
  MAD_MULTIPLIER = 5.0,       # MAD multiplier for outlier detection
  USE_ADAPTIVE = TRUE,        # Use adaptive thresholds
  USE_OD = TRUE,
  PRIOR_COUNT = 0.5,
  MC_CORES = min(3L, max(1L, parallel::detectCores() - 1L))
)

# Group-specific seed generation
seed_from_name <- function(name, base = 123L) {
  if (requireNamespace("digest", quietly = TRUE)) {
    hx <- digest::digest(paste0(base, "::", name), algo = "xxhash32")
    high <- strtoi(substr(hx, 1, 4), base = 16L)
    low <- strtoi(substr(hx, 5, 8), base = 16L)
    seed_val <- bitwOr(bitwShiftL(high, 16L), low)
    seed_val <- bitwAnd(seed_val, 0x7FFFFFFF)
    if(is.na(seed_val) || seed_val <= 0) {
      seed_val <- abs(sum(utf8ToInt(name)) * base) %% 2147483647L
    }
    as.integer(max(1L, seed_val))
  } else {
    as.integer(max(1L, abs(sum(utf8ToInt(name)) * base) %% 2147483647L))
  }
}

# PA function with column centering
pa_select_K <- function(X, R = 200, alpha = 0.05, L = 8, seed = 123) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  X_t <- t(X)  # samples x genes
  n <- nrow(X_t); p <- ncol(X_t)
  L <- min(L, n - 1L, p)
  if (L < 1L) return(list(K = 1L, real = numeric(0), threshold = numeric(0)))
  
  # Column centering
  mu <- colMeans(X_t)
  center <- function(M) sweep(M, 2, mu, "-")
  
  has_rs <- requireNamespace("RSpectra", quietly = TRUE)
  topL <- function(M) {
    if (has_rs && L < min(n, p) / 2) {
      sv2 <- try(RSpectra::svds(center(M), k = L, nu = 0, nv = 0,
                                opts = list(tol = 1e-6, maxitr = 1000))$d^2,
                 silent = TRUE)
      if (inherits(sv2, "try-error")) {
        sv <- svd(center(M), nu = 0, nv = 0)$d^2
        sv[seq_len(min(L, length(sv)))]
      } else sv2
    } else {
      sv <- svd(center(M), nu = 0, nv = 0)$d^2
      sv[seq_len(min(L, length(sv)))]
    }
  }
  
  real_sv  <- topL(X_t)
  actual_L <- length(real_sv)
  cat("  Running PA (", R, " permutations, L=", actual_L, ")...\n", sep = "")
  
  null_sv <- replicate(R, {
    Xp <- apply(X_t, 2, sample)   # Column-wise shuffle
    topL(Xp)
  })
  if (is.null(dim(null_sv))) null_sv <- matrix(null_sv, nrow = 1L)
  
  thr <- apply(null_sv, 1, quantile, probs = 1 - alpha, na.rm = TRUE)
  K   <- max(1L, sum(real_sv > thr * 1.00001))
  cat("  PA selected K =", K, "\n")
  
  list(K = K, real = real_sv, threshold = thr)
}

# Enhanced outlier detection with robust safeguards
cdm_outlier_detection <- function(cdm_result, X, K, config = CONFIG) {
  # X orientation → samples x genes
  if (ncol(X) == nrow(cdm_result$vectors)) Xuse <- X else Xuse <- t(X)
  n <- nrow(Xuse)
  
  # Define outlier constraints upfront
  min_outliers <- ceiling(n * config$MIN_OUTLIERS_PCT)
  max_outliers <- floor(n * config$MAX_OUTLIERS_PCT)
  
  # Column centering
  mu <- colMeans(Xuse)
  Xc <- sweep(Xuse, 2, mu, "-")
  
  K <- max(1L, min(K, n - 1L, ncol(cdm_result$vectors)))
  V <- qr.Q(qr(cdm_result$vectors[, 1:K, drop = FALSE]))
  S <- Xc %*% V  # Projection scores
  
  # SD calculation with robust normalization
  Z <- apply(S, 2, function(x) {
    m <- median(x); s <- mad(x, constant = 1.4826)
    if (!is.finite(s) || s <= 1e-8) s <- IQR(x)/1.349
    if (!is.finite(s) || s <= 1e-8) return(rep(NA_real_, length(x)))
    (x - m) / s
  })
  
  keep_pc <- which(colSums(!is.na(Z)) == n)
  if (length(keep_pc) == 0L) stop("All PCs dropped due to ~zero MAD/IQR.")
  if (length(keep_pc) < K) message("  Dropped ", K - length(keep_pc), " PC(s) with ~zero spread.")
  
  Z  <- Z[, keep_pc, drop = FALSE]
  SD <- sqrt(rowSums(Z^2))
  
  # Adaptive threshold for SD with robust guards
  if (config$USE_ADAPTIVE) {
    # Robust statistics with guards
    q1 <- as.numeric(quantile(SD, 0.25, na.rm = TRUE))
    q3 <- as.numeric(quantile(SD, 0.75, na.rm = TRUE))
    iqr <- q3 - q1
    med_SD <- median(SD, na.rm = TRUE)
    mad_SD <- mad(SD, constant = 1.4826, na.rm = TRUE)
    
    # Calculate thresholds with guards against zero variance
    thr_SD_iqr <- if (is.finite(iqr) && iqr > 1e-8) {
      q3 + config$IQR_MULTIPLIER * iqr
    } else {
      Inf
    }
    
    thr_SD_mad <- if (is.finite(mad_SD) && mad_SD > 1e-8) {
      med_SD + config$MAD_MULTIPLIER * mad_SD
    } else {
      Inf
    }
    
    # Use the more conservative (higher) threshold
    thr_SD <- max(thr_SD_iqr, thr_SD_mad, na.rm = TRUE)
    
    # Fallback if both methods fail
    if (!is.finite(thr_SD)) {
      thr_SD <- as.numeric(quantile(SD, 0.975, na.rm = TRUE))
      message("  Warning: Using fallback quantile threshold for SD")
    }
    
    # Ensure minimum and maximum outlier constraints
    sorted_SD <- sort(SD, decreasing = TRUE)
    if (min_outliers > 0 && sum(SD > thr_SD) < min_outliers) {
      thr_SD <- sorted_SD[min_outliers] - 1e-10
    }
    if (max_outliers > 0 && sum(SD > thr_SD) > max_outliers) {
      thr_SD <- sorted_SD[max_outliers] + 1e-10
    }
  } else {
    # Fall back to quantile method
    thr_SD <- quantile(SD, 0.975, na.rm = TRUE)
  }
  
  hit_SD <- SD > thr_SD
  
  # OD calculation with adaptive threshold
  hit_OD <- rep(FALSE, n)
  thr_OD <- NA_real_
  OD <- NULL
  
  if (config$USE_OD) {
    OD_raw <- sqrt(pmax(rowSums(Xc^2) - rowSums(S^2), 0))
    OD <- log1p(OD_raw)
    
    if (config$USE_ADAPTIVE) {
      # Robust statistics with guards
      q1_OD <- as.numeric(quantile(OD, 0.25, na.rm = TRUE))
      q3_OD <- as.numeric(quantile(OD, 0.75, na.rm = TRUE))
      iqr_OD <- q3_OD - q1_OD
      med_OD <- median(OD, na.rm = TRUE)
      mad_OD <- mad(OD, constant = 1.4826, na.rm = TRUE)
      
      # Calculate thresholds with guards
      thr_OD_iqr <- if (is.finite(iqr_OD) && iqr_OD > 1e-8) {
        q3_OD + config$IQR_MULTIPLIER * iqr_OD
      } else {
        Inf
      }
      
      thr_OD_mad <- if (is.finite(mad_OD) && mad_OD > 1e-8) {
        med_OD + config$MAD_MULTIPLIER * mad_OD
      } else {
        Inf
      }
      
      thr_OD <- max(thr_OD_iqr, thr_OD_mad, na.rm = TRUE)
      
      # Fallback if both methods fail
      if (!is.finite(thr_OD)) {
        thr_OD <- as.numeric(quantile(OD, 0.975, na.rm = TRUE))
        message("  Warning: Using fallback quantile threshold for OD")
      }
      
      # Apply constraints
      sorted_OD <- sort(OD, decreasing = TRUE)
      if (min_outliers > 0 && sum(OD > thr_OD) < min_outliers) {
        thr_OD <- sorted_OD[min_outliers] - 1e-10
      }
      if (max_outliers > 0 && sum(OD > thr_OD) > max_outliers) {
        thr_OD <- sorted_OD[max_outliers] + 1e-10
      }
    } else {
      thr_OD <- quantile(OD, 0.975, na.rm = TRUE)
    }
    
    hit_OD <- OD > thr_OD
  }
  
  # Combine outliers
  outliers <- hit_SD | hit_OD
  
  # Enforce global maximum outlier constraint
  if (max_outliers > 0L && sum(outliers) > max_outliers) {
    # Calculate robust Z-scores for ranking
    mad_SD_global <- mad(SD, constant = 1.4826, na.rm = TRUE)
    zSD <- if (mad_SD_global > 1e-8) {
      (SD - thr_SD) / mad_SD_global
    } else {
      SD - thr_SD  # Fallback to raw difference
    }
    
    if (config$USE_OD && !is.null(OD)) {
      mad_OD_global <- mad(OD, constant = 1.4826, na.rm = TRUE)
      zOD <- if (mad_OD_global > 1e-8) {
        (OD - thr_OD) / mad_OD_global
      } else {
        OD - thr_OD  # Fallback to raw difference
      }
    } else {
      zOD <- rep(-Inf, length(zSD))
    }
    
    # Take maximum deviation across metrics
    fused <- pmax(zSD, zOD, na.rm = TRUE)
    
    # Keep only top max_outliers
    keep_idx <- order(fused, decreasing = TRUE)[seq_len(max_outliers)]
    new_outliers <- rep(FALSE, length(outliers))
    new_outliers[keep_idx] <- TRUE
    
    n_reduced <- sum(outliers) - sum(new_outliers)
    if (n_reduced > 0) {
      message(sprintf("  Reduced outliers from %d to %d (max constraint)", 
                      sum(outliers), sum(new_outliers)))
    }
    outliers <- new_outliers
  }
  
  # Return detailed results for diagnostics
  list(
    outliers = outliers,
    n_outliers = sum(outliers),
    K = ncol(V),
    SD = SD,
    OD = OD,
    thr_SD = thr_SD,
    thr_OD = thr_OD,
    hit_SD = sum(hit_SD),
    hit_OD = if(config$USE_OD) sum(hit_OD) else 0,
    min_outliers = min_outliers,
    max_outliers = max_outliers
  )
}

# Process one group with diagnostics
process_group <- function(se, group_samples, group_name, keep_genes) {
  # Ensure single thread in worker
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  } else {
    Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1")
  }
  
  cat("\nProcessing", group_name, "(n =", length(group_samples), ")\n")
  
  # Extract counts with unified gene set
  counts <- assay(se, "stranded_second")[keep_genes, group_samples, drop = FALSE]
  if (any(is.na(counts))) counts[is.na(counts)] <- 0
  cat("  Using", nrow(counts), "common genes\n")
  
  # logCPM transformation
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge)
  logcpm <- cpm(dge, log = TRUE, prior.count = CONFIG$PRIOR_COUNT)
  
  # PA with group-specific seed
  pa_result <- pa_select_K(
    logcpm,
    R = CONFIG$PA_R,
    alpha = CONFIG$PA_ALPHA,
    L = CONFIG$PA_L,
    seed = seed_from_name(group_name, CONFIG$PA_SEED)
  )
  
  # CDM
  cdm_result <- CDM_fast3_arma(logcpm, verbose = FALSE)
  
  # Outlier detection with adaptive thresholds
  outlier_result <- cdm_outlier_detection(cdm_result, logcpm, K = pa_result$K, config = CONFIG)
  
  # Diagnostic output
  cat(sprintf("  SD outliers: %d (threshold: %.3f)\n", 
              outlier_result$hit_SD, outlier_result$thr_SD))
  if (CONFIG$USE_OD) {
    cat(sprintf("  OD outliers: %d (threshold: %.3f)\n", 
                outlier_result$hit_OD, outlier_result$thr_OD))
  }
  cat(sprintf("  Total outliers: %d / %d (%.1f%%), constraints: [%d, %d]\n", 
              outlier_result$n_outliers, length(group_samples),
              outlier_result$n_outliers / length(group_samples) * 100,
              outlier_result$min_outliers, outlier_result$max_outliers))
  
  # Return both indices and IDs
  idx_out <- which(outlier_result$outliers)
  list(
    samples_idx = group_samples,
    samples_id  = colnames(se)[group_samples],
    out_idx     = group_samples[idx_out],
    out_id      = colnames(se)[group_samples[idx_out]],
    n_outliers  = outlier_result$n_outliers,
    K           = pa_result$K,
    diagnostics = list(
      SD_range = range(outlier_result$SD),
      OD_range = if(!is.null(outlier_result$OD)) range(outlier_result$OD) else NULL,
      thr_SD = outlier_result$thr_SD,
      thr_OD = outlier_result$thr_OD,
      hit_SD = outlier_result$hit_SD,
      hit_OD = outlier_result$hit_OD
    )
  )
}

# Main analysis
cat("\nLoading data...\n")
se <- readRDS(paste0(paths$processed, "thyr_se_merged_err.rds"))
metadata <- as.data.frame(colData(se))

# Ensure column names
if(is.null(colnames(se))) {
  colnames(se) <- rownames(metadata)
}

# Define groups
groups <- list(
  R0_tumor = which(!is.na(metadata$group) & metadata$group == "R0" & metadata$is_tumor),
  R0_normal = which(!is.na(metadata$group) & metadata$group == "R0" & metadata$is_normal),
  R1_tumor = which(!is.na(metadata$group) & metadata$group == "R1" & metadata$is_tumor),
  R1_normal = which(!is.na(metadata$group) & metadata$group == "R1" & metadata$is_normal),
  B0_tumor = which(!is.na(metadata$group) & metadata$group == "B0" & metadata$is_tumor),
  B0_normal = which(!is.na(metadata$group) & metadata$group == "B0" & metadata$is_normal),
  B1_tumor = which(!is.na(metadata$group) & metadata$group == "B1" & metadata$is_tumor),
  B1_normal = which(!is.na(metadata$group) & metadata$group == "B1" & metadata$is_normal)
)

# Filter groups with sufficient samples
groups <- groups[sapply(groups, length) > 3]
cat("\nGroups to process:", paste(names(groups), collapse = ", "), "\n")

# Unified gene filtering
cat("\nApplying unified gene filtering...\n")
all_samples <- unlist(groups, use.names = FALSE)
all_counts <- assay(se, "stranded_second")[, all_samples, drop = FALSE]

if (any(is.na(all_counts))) {
  cat("Warning: Found", sum(is.na(all_counts)), "NA values. Replacing with 0.\n")
  all_counts[is.na(all_counts)] <- 0
}

group_labels <- factor(rep(names(groups), sapply(groups, length)))
keep_genes <- filterByExpr(all_counts, group = group_labels)
cat("Common genes retained:", sum(keep_genes), "/", length(keep_genes), "\n")

# Process groups
run_one <- function(gn) {
  try(process_group(se, groups[[gn]], gn, keep_genes), silent = TRUE)
}

# Sort groups by size
order_g <- names(groups)[order(sapply(groups, length), decreasing = TRUE)]

cat("\nProcessing groups",
    if(.Platform$OS.type == "unix" && CONFIG$MC_CORES > 1L) {
      paste0(" in parallel (mc.cores=", CONFIG$MC_CORES, ")")
    } else {
      " sequentially"
    }, "...\n", sep = "")

if(.Platform$OS.type == "unix" && CONFIG$MC_CORES > 1L) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(CONFIG$PA_SEED)
  
  res_list <- mclapply(
    order_g, 
    run_one,
    mc.cores = CONFIG$MC_CORES,
    mc.set.seed = TRUE,
    mc.preschedule = FALSE
  )
  names(res_list) <- order_g
  results <- res_list[names(groups)]
} else {
  results <- setNames(lapply(names(groups), run_one), names(groups))
}

# Check for errors
failed <- sapply(results, function(x) inherits(x, "try-error"))
if(any(failed)) {
  warning("Failed groups: ", paste(names(results)[failed], collapse = ", "))
}

# Aggregate results
all_out_idx <- integer(0)
all_out_id  <- character(0)

for (group_name in names(results)) {
  if(!inherits(results[[group_name]], "try-error")) {
    res <- results[[group_name]]
    all_out_idx <- c(all_out_idx, res$out_idx)
    all_out_id  <- c(all_out_id,  res$out_id)
  }
}

# Summary
cat("\n=== Summary ===\n")
for (g in names(results)) {
  if(!inherits(results[[g]], "try-error")) {
    r <- results[[g]]
    cat(sprintf("%s: %d samples, K=%d, %d outliers (SD:%d, OD:%d)\n",
                g, length(r$samples_idx), r$K, r$n_outliers,
                r$diagnostics$hit_SD, r$diagnostics$hit_OD))
  }
}

all_out_idx <- sort(unique(all_out_idx))
all_out_id  <- unique(all_out_id)

cat("\nTotal unique outliers:", length(all_out_idx), "\n")
cat("Total samples processed:", length(all_samples), "\n")
cat("Retention rate:", sprintf("%.1f%%\n",
                               (length(all_samples) - length(all_out_idx)) / length(all_samples) * 100))

# Save results
cat("\nSaving results...\n")
saveRDS(results, paste0(paths$output, "stage1_pca_results.rds"))

# Create filtered SE
keep_samples_idx <- setdiff(seq_len(ncol(se)), all_out_idx)
se_filtered <- se[, keep_samples_idx]
saveRDS(se_filtered, paste0(paths$processed, "thyr_se_stage1_filtered.rds"))

# Save outlier list with diagnostics
outlier_df <- data.frame(
  sample_id = all_out_id,
  stage = "Stage1_PCA",
  stringsAsFactors = FALSE
)
write.csv(outlier_df, paste0(paths$output, "stage1_outliers.csv"), row.names = FALSE)

# Save diagnostic information
diagnostics_list <- lapply(results, function(r) {
  if(!inherits(r, "try-error") && !is.null(r$diagnostics)) {
    r$diagnostics
  }
})
saveRDS(diagnostics_list, paste0(paths$output, "stage1_diagnostics.rds"))

cat("\n=== Stage 1 Complete ===\n")
cat("Files created:\n")
cat("  - stage1_pca_results.rds\n")
cat("  - thyr_se_stage1_filtered.rds\n")
cat("  - stage1_outliers.csv\n")
cat("  - stage1_diagnostics.rds\n")
cat("\nConfiguration used:\n")
cat("  Method: ", ifelse(CONFIG$USE_ADAPTIVE, "Adaptive (IQR/MAD)", "Fixed quantile"), "\n")
cat("  IQR multiplier: ", CONFIG$IQR_MULTIPLIER, "\n")
cat("  MAD multiplier: ", CONFIG$MAD_MULTIPLIER, "\n")
cat("  Max outliers: ", CONFIG$MAX_OUTLIERS_PCT * 100, "%\n")
cat("\nNext: Run 08_purity_analysis.R\n")