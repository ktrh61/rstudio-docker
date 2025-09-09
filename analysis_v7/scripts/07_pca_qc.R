# 07_pca_qc.R - Robust PCA analysis with outlier detection
# Purpose: CDM-based PCA with robust dimension selection and outlier detection

source("analysis_v7/setup.R")
source("analysis_v7/utils/robust_pca_utils.R")

cat("\n=== Robust PCA Quality Control (v7) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# Compile CDM
Rcpp::sourceCpp("analysis_v7/utils/CDM_fast3_arma_enhanced.cpp")

# Set RNG for reproducibility
RNGkind("L'Ecuyer-CMRG")
set.seed(12345)

# Control BLAS threads
old_omp <- Sys.getenv("OMP_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
if(requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
}

# Load data
se <- readRDS(paste0(paths$processed, "thyr_se_merged_err.rds"))
metadata <- as.data.frame(colData(se))

# Separate main/marker samples
valid_samples <- !is.na(metadata$group)
se_main <- se[, valid_samples]
se_marker <- se[, !valid_samples]
metadata_main <- metadata[valid_samples, ]

cat("\nSamples: Total", ncol(se), "| Main", ncol(se_main), "| Marker", ncol(se_marker), "\n")

saveRDS(se_marker, paste0(paths$processed, "thyr_se_marker.rds"))

# Extract counts
counts <- assay(se_main, "stranded_second")

# Filter genes
library(edgeR)
keep <- filterByExpr(counts, group = metadata_main$group)
counts_filtered <- counts[keep, ]
cat("Genes:", nrow(counts_filtered), "/", nrow(counts), "\n")

# Calculate logCPM
logcpm <- cpm(counts_filtered, log = TRUE, prior.count = 0.5)

# Define groups
pca_groups <- c("R0_normal", "R0_tumor", "R1_normal", "R1_tumor",
                "B0_normal", "B0_tumor", "B1_normal", "B1_tumor")

# Prepare parallel processing
library(parallel)
pca_results <- list()
outliers_all <- list()

# Function to process one group
process_group <- function(grp, seed_offset = 0) {
  parts <- strsplit(grp, "_")[[1]]
  group_name <- parts[1]
  sample_type <- parts[2]
  
  idx <- which(metadata_main$group == group_name & 
               grepl(sample_type, metadata_main$sample_type, ignore.case = TRUE))
  
  if(length(idx) < 3) {
    return(list(grp = grp, n_samples = length(idx), skip = TRUE))
  }
  
  # Run robust PCA
  X <- logcpm[, idx, drop = FALSE]
  robust_result <- cdm_robust_pca(t(X),  # Note: transpose for CDM
                                  alpha = 0.05,
                                  B = 100,
                                  seed = 12345 + seed_offset,
                                  verbose = FALSE)
  
  # Calculate outliers using spatial median and MAD
  K <- robust_result$K
  scores <- robust_result$cdm$scores[, 1:K, drop = FALSE]
  
  # MAD standardization
  scores_std <- apply(scores, 2, function(x) {
    med <- median(x)
    mad_val <- mad(x, constant = 1.4826)
    if(mad_val < 1e-8) mad_val <- IQR(x) / 1.349
    (x - med) / max(mad_val, 1e-8)
  })
  
  # Spatial median (simplified - use column medians)
  center <- apply(scores_std, 2, median)
  
  # Euclidean distances
  distances <- sqrt(rowSums((scores_std - rep(center, each = nrow(scores_std)))^2))
  
  # IQR outlier detection
  Q1 <- quantile(distances, 0.25)
  Q3 <- quantile(distances, 0.75)
  IQR <- Q3 - Q1
  threshold <- Q3 + 1.5 * IQR
  outliers <- which(distances > threshold)
  
  result <- list(
    grp = grp,
    n_samples = length(idx),
    K = K,
    median_stability = robust_result$final_median_cc,
    n_outliers = length(outliers),
    outlier_indices = idx[outliers],
    outlier_samples = metadata_main$sample_id[idx[outliers]],
    outlier_distances = distances[outliers],
    robust_result = robust_result
  )
  
  return(result)
}

# Process groups in parallel
# 07_pca_qc.R内の並列処理部分
cat("\nProcessing groups (3 cores)...\n")
results <- mclapply(seq_along(pca_groups), function(i) {
  process_group(pca_groups[i], seed_offset = i)
}, mc.cores = 3, mc.preschedule = TRUE, mc.set.seed = TRUE)

# Process results
for(res in results) {
  grp <- res$grp
  cat("\n", grp, ": n =", res$n_samples)
  
  if(!is.null(res$skip) && res$skip) {
    cat(" [SKIPPED - insufficient samples]")
    next
  }
  
  cat(", K =", res$K, ", stability =", sprintf("%.3f", res$median_stability))
  cat(", outliers =", res$n_outliers)
  
  if(res$n_outliers > 0) {
    cat("\n  Samples:", paste(res$outlier_samples, collapse = ", "))
    outliers_all[[grp]] <- data.frame(
      sample_id = res$outlier_samples,
      case_id = metadata_main$case_id[res$outlier_indices],
      group = grp,
      distance = res$outlier_distances,
      stringsAsFactors = FALSE
    )
  }
  
  pca_results[[grp]] <- res
}

# Save outlier information
if(length(outliers_all) > 0) {
  outliers_df <- do.call(rbind, outliers_all)
  cat("\n\n=== Total outliers:", nrow(outliers_df), "===\n")
  
  outlier_cases <- unique(outliers_df$case_id)
  cat("Cases with outliers:", paste(outlier_cases, collapse = ", "), "\n")
  
  data.table::fwrite(outliers_df, paste0(paths$processed, "pca_outliers_robust.tsv"), sep = "\t")
  
  metadata_main$is_outlier <- metadata_main$sample_id %in% outliers_df$sample_id
  metadata_main$outlier_case <- metadata_main$case_id %in% outlier_cases
} else {
  cat("\n\nNo outliers detected\n")
  metadata_main$is_outlier <- FALSE
  metadata_main$outlier_case <- FALSE
}

# Restore BLAS settings
Sys.setenv(OMP_NUM_THREADS = old_omp)

# Save results
colData(se_main) <- DataFrame(metadata_main)
saveRDS(se_main, paste0(paths$processed, "thyr_se_pca_robust.rds"))
saveRDS(pca_results, paste0(paths$processed, "pca_results_robust.rds"))

cat("\n=== Files created ===\n")
cat("  - thyr_se_pca_robust.rds\n")
cat("  - pca_results_robust.rds\n")
cat("  - pca_outliers_robust.tsv (if outliers found)\n")
cat("  - thyr_se_marker.rds\n")
