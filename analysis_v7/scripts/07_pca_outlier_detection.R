# 07_pca_outlier_detection.R - Stage 1 PCA Outlier Detection (v7)
# Purpose: Remove technical outliers using CDM with PA-selected K
# Method: Unified gene set across all groups, parallel processing

source("analysis_v7/setup.R")

cat("\n=== Stage 1: PCA Outlier Detection (v7) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# Load packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(Rcpp)
  library(RcppArmadillo)
  library(parallel)
})

# --- single-thread lock & RNG kind (再現性の土台) ---
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

# Configuration
CONFIG <- list(
  PA_R = 200,
  PA_ALPHA = 0.05,
  PA_L = 8,
  PA_SEED = 123,
  OUTLIER_ALPHA = 0.975,
  USE_OD = TRUE,
  PRIOR_COUNT = 0.5,
  MC_CORES = min(3L, max(1L, parallel::detectCores() - 1L))
)

# Group-specific seed generation
seed_from_name <- function(name, base = 123L) {
  if (requireNamespace("digest", quietly = TRUE)) {
    # xxhash32で32bit hashを生成
    hx <- digest::digest(paste0(base, "::", name), algo = "xxhash32")
   
    # 8文字の16進数を2つの4文字に分割して処理（オーバーフロー回避）
    high <- strtoi(substr(hx, 1, 4), base = 16L)
    low <- strtoi(substr(hx, 5, 8), base = 16L)
   
    # 32bit整数に結合して31bitマスク
    seed_val <- bitwOr(bitwShiftL(high, 16L), low)
    seed_val <- bitwAnd(seed_val, 0x7FFFFFFF)
   
    # 正の整数を保証
    if(is.na(seed_val) || seed_val <= 0) {
      # フォールバック：名前の文字コードベース
      seed_val <- abs(sum(utf8ToInt(name)) * base) %% 2147483647L
    }
   
    as.integer(max(1L, seed_val))
  } else {
    # digestパッケージがない場合
    as.integer(max(1L, abs(sum(utf8ToInt(name)) * base) %% 2147483647L))
  }
}

# PA function with column centering
pa_select_K <- function(X, R = 200, alpha = 0.05, L = 8, seed = 123) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  X_t <- t(X)                            # samples x genes
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

# Outlier detection with centering
cdm_outlier_detection <- function(cdm_result, X, K, alpha = 0.975, use_OD = TRUE) {
  # X orientation → samples x genes
  if (ncol(X) == nrow(cdm_result$vectors)) Xuse <- X else Xuse <- t(X)
  n <- nrow(Xuse)
 
  # Column centering
  mu <- colMeans(Xuse)
  Xc <- sweep(Xuse, 2, mu, "-")

  K <- max(1L, min(K, n - 1L, ncol(cdm_result$vectors)))
  V <- qr.Q(qr(cdm_result$vectors[, 1:K, drop = FALSE]))
  S <- Xc %*% V                                           # Projection scores

  # SD calculation
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
  thr_SD <- quantile(SD, alpha, na.rm = TRUE)
  hit_SD <- SD > thr_SD

  # OD calculation
  hit_OD <- rep(FALSE, n)
  if (use_OD) {
    OD_raw <- sqrt(pmax(rowSums(Xc^2) - rowSums(S^2), 0))
    OD <- log1p(OD_raw)
    thr_OD <- quantile(OD, alpha, na.rm = TRUE)
    hit_OD <- OD > thr_OD
  }

  outliers <- hit_SD | hit_OD
  list(outliers = outliers, n_outliers = sum(outliers), K = ncol(V))
}

# Process one group
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

  # Outlier detection
  outlier_result <- cdm_outlier_detection(cdm_result, logcpm, K = pa_result$K,
                                          alpha = CONFIG$OUTLIER_ALPHA, use_OD = CONFIG$USE_OD)

  cat("  Outliers:", outlier_result$n_outliers, "/", length(group_samples), "\n")

  # Return both indices and IDs
  idx_out <- which(outlier_result$outliers)
  list(
    samples_idx = group_samples,
    samples_id  = colnames(se)[group_samples],
    out_idx     = group_samples[idx_out],
    out_id      = colnames(se)[group_samples[idx_out]],
    n_outliers  = outlier_result$n_outliers,
    K           = pa_result$K
  )
}

# Main analysis
cat("\nLoading data...\n")
se <- readRDS(paste0(paths$processed, "thyr_se_merged_err.rds"))
metadata <- as.data.frame(colData(se))

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

# Process groups (parallel or sequential)
run_one <- function(gn) {
  try(process_group(se, groups[[gn]], gn, keep_genes), silent = TRUE)
}

# Sort groups by size (largest first for better load balancing)
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
    mc.preschedule = FALSE  # Dynamic scheduling
  )
  names(res_list) <- order_g
  results <- res_list[names(groups)]  # Restore original order
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
    cat(sprintf("%s: %d samples, K=%d, %d outliers\n",
                g, length(results[[g]]$samples_idx), results[[g]]$K, results[[g]]$n_outliers))
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

# Create filtered SE (remove outliers by index)
keep_samples_idx <- setdiff(seq_len(ncol(se)), all_out_idx)
se_filtered <- se[, keep_samples_idx]
saveRDS(se_filtered, paste0(paths$processed, "thyr_se_stage1_filtered.rds"))

# Save outlier list (human-readable with IDs)
outlier_df <- data.frame(
  sample_id = all_out_id,
  stage = "Stage1_PCA",
  stringsAsFactors = FALSE
)
write.csv(outlier_df, paste0(paths$output, "stage1_outliers.csv"), row.names = FALSE)

cat("\n=== Stage 1 Complete ===\n")
cat("Files created:\n")
cat("  - stage1_pca_results.rds\n")
cat("  - thyr_se_stage1_filtered.rds\n")
cat("  - stage1_outliers.csv\n")
cat("\nNext: Run 08_purity_analysis.R\n")
