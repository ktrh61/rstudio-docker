# ==============================================================================
# REBC-THYR PCA Analysis Script v3 - New CDM with Pair Consistency
# 05_pca_analysis_v3.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(dplyr)
library(Rcpp)
library(RhpcBLASctl)

# Source the new CDM functions
source("./utils/with_openblas_threads.R")
sourceCpp("./utils/CDM_fast3_arma_enhanced.cpp")

cat("Starting PCA analysis v3 with new CDM and pair consistency...\n")

# ==============================================================================
# 1. Define CDM_fast Compatible Wrapper
# ==============================================================================

CDM_fast_compatible <- function(X, by_sample = FALSE, center = TRUE, scale. = TRUE, 
                                k = NULL, return_scores = TRUE, verbose = FALSE) {
  
  # サンプル名の保存
  if (by_sample) {
    sample_names <- rownames(X)
    X <- t(X)  # genes x samples に変換
  } else {
    sample_names <- colnames(X)
  }
  
  # 新CDMの実行
  result <- with_openblas_threads("auto-2", {
    CDM_fast3_arma(X, verbose = verbose)
  })
  
  # サンプル名とPC名の設定
  if (!is.null(result$scores) && !is.null(sample_names)) {
    rownames(result$scores) <- sample_names
    colnames(result$scores) <- paste0("PC", 1:ncol(result$scores))
  }
  
  # k制限の適用
  if (!is.null(k) && k < length(result$values)) {
    result$values <- result$values[1:k]
    result$vectors <- result$vectors[, 1:k, drop = FALSE]
    if (!is.null(result$scores)) {
      result$scores <- result$scores[, 1:k, drop = FALSE]
    }
  }
  
  # 分散説明率の計算
  variance_explained <- result$values^2 / sum(result$values^2)
  cumulative_variance <- cumsum(variance_explained)
  
  result$variance_explained <- variance_explained
  result$cumulative_variance <- cumulative_variance
  result$n_components <- length(result$values)
  
  return(result)
}

cat("CDM_fast compatible wrapper defined.\n")

# ==============================================================================
# 2. Load Data and Sample Lists
# ==============================================================================

cat("Loading data and sample lists...\n")

# Load sample lists from ERR integration (before any filtering)
load("./data/processed/sample_lists.rda")

# Load TPM data
load("./data/processed/thyr_tpm.rda")

cat("TPM data dimensions:", dim(tpm), "\n")
cat("Available groups:", names(sample_lists), "\n")

# ==============================================================================
# 3. Define Analysis Groups (Phase 1 - Individual Groups)
# ==============================================================================

cat("Setting up Phase 1 analysis groups...\n")

analysis_groups <- list(
  "R0_normal" = list(group = "R0", tissue = "normal", description = "RET Unexposed Normal"),
  "R0_tumor" = list(group = "R0", tissue = "tumor", description = "RET Unexposed Tumor"),
  "R1_normal" = list(group = "R1", tissue = "normal", description = "RET RadHigh Normal"),
  "R1_tumor" = list(group = "R1", tissue = "tumor", description = "RET RadHigh Tumor"),
  "B0_normal" = list(group = "B0", tissue = "normal", description = "BRAF Unexposed Normal"),
  "B0_tumor" = list(group = "B0", tissue = "tumor", description = "BRAF Unexposed Tumor"),
  "B1_normal" = list(group = "B1", tissue = "normal", description = "BRAF RadHigh Normal"),
  "B1_tumor" = list(group = "B1", tissue = "tumor", description = "BRAF RadHigh Tumor")
)

# ==============================================================================
# 4. Perform PCA Analysis with New CDM
# ==============================================================================

cat("Performing PCA analysis with new CDM for each group...\n")

pca_results <- list()

for (analysis_name in names(analysis_groups)) {
  cat(sprintf("\n--- Analyzing %s ---\n", analysis_name))
  
  group_info <- analysis_groups[[analysis_name]]
  group_name <- group_info$group
  tissue_type <- group_info$tissue
  description <- group_info$description
  
  # Get samples
  if (!group_name %in% names(sample_lists)) {
    cat(sprintf("Group %s not found, skipping...\n", group_name))
    next
  }
  
  if (tissue_type == "tumor") {
    sample_ids <- sample_lists[[group_name]]$tumor
  } else {
    sample_ids <- sample_lists[[group_name]]$normal
  }
  
  if (length(sample_ids) == 0) {
    cat(sprintf("No samples for %s, skipping...\n", analysis_name))
    next
  }
  
  # Check availability in TPM data
  available_samples <- sample_ids[sample_ids %in% colnames(tpm)]
  
  if (length(available_samples) < 3) {
    cat(sprintf("Insufficient samples (%d) for %s, skipping...\n", 
                length(available_samples), analysis_name))
    next
  }
  
  cat(sprintf("%s: %d samples\n", description, length(available_samples)))
  
  # Extract and filter TPM data
  group_tpm <- tpm[, available_samples, drop = FALSE]
  non_zero_genes <- rowSums(group_tpm > 0) > 0
  group_tpm_filtered <- group_tpm[non_zero_genes, ]
  
  cat(sprintf("Using %d genes for PCA (removed %d all-zero genes)\n", 
              nrow(group_tpm_filtered), sum(!non_zero_genes)))
  
  # Run PCA with new CDM
  tryCatch({
    pca_result <- CDM_fast_compatible(
      X = group_tpm_filtered,
      by_sample = FALSE,
      center = TRUE,
      scale. = TRUE,
      k = min(10, ncol(group_tpm_filtered) - 1),
      return_scores = TRUE,
      verbose = FALSE
    )
    
    # Store results
    pca_results[[analysis_name]] <- list(
      pca = pca_result,
      samples = available_samples,
      description = description,
      n_genes = nrow(group_tpm_filtered),
      n_samples = length(available_samples),
      analysis_name = analysis_name,
      group = group_name,
      tissue = tissue_type
    )
    
    cat(sprintf("PCA completed: %d components, %.1f%% variance explained by PC1-2\n",
                pca_result$n_components,
                sum(pca_result$variance_explained[1:min(2, length(pca_result$variance_explained))]) * 100))
    
  }, error = function(e) {
    cat(sprintf("Error in PCA for %s: %s\n", analysis_name, e$message))
    pca_results[[analysis_name]] <- NULL
  })
}

# ==============================================================================
# 5. Enhanced Outlier Detection Function
# ==============================================================================

enhanced_outlier_detection <- function(analysis_name, pca_data, verbose = TRUE) {
  
  if (verbose) cat(sprintf("\n--- Outlier detection for %s ---\n", analysis_name))
  
  scores <- pca_data$pca$scores
  if (is.null(scores) || nrow(scores) < 4) {
    if (verbose) cat(sprintf("Insufficient data for outlier detection\n"))
    return(NULL)
  }
  
  # Method 1: Mahalanobis distance
  n_pcs_for_outlier <- min(3, ncol(scores))
  pc_subset <- scores[, 1:n_pcs_for_outlier, drop = FALSE]
  
  center_point <- colMeans(pc_subset)
  cov_matrix <- cov(pc_subset)
  
  mahal_outliers <- integer(0)
  mahal_dist <- rep(NA, nrow(scores))
  threshold <- qchisq(0.95, df = n_pcs_for_outlier)
  
  if (det(cov_matrix) > .Machine$double.eps) {
    mahal_dist <- mahalanobis(pc_subset, center_point, cov_matrix)
    mahal_outliers <- which(mahal_dist > threshold)
  }
  
  # Method 2: IQR method (PC1 focused)
  pc1_values <- scores[, 1]
  q1 <- quantile(pc1_values, 0.25)
  q3 <- quantile(pc1_values, 0.75)
  iqr <- q3 - q1
  
  lower_fence <- q1 - 1.5 * iqr
  upper_fence <- q3 + 1.5 * iqr
  
  iqr_outliers <- which(pc1_values < lower_fence | pc1_values > upper_fence)
  
  # Combined outliers
  combined_outliers <- unique(c(mahal_outliers, iqr_outliers))
  combined_samples <- rownames(scores)[combined_outliers]
  
  if (verbose) {
    cat(sprintf("  Mahalanobis outliers: %d\n", length(mahal_outliers)))
    cat(sprintf("  IQR outliers: %d\n", length(iqr_outliers)))
    cat(sprintf("  Combined outliers: %d\n", length(combined_outliers)))
    if (length(combined_outliers) > 0) {
      cat(sprintf("  Outlier samples: %s\n", paste(combined_samples, collapse = ", ")))
    }
  }
  
  return(list(
    mahal_outliers = mahal_outliers,
    iqr_outliers = iqr_outliers,
    combined_outliers = combined_outliers,
    combined_samples = combined_samples,
    total_samples = nrow(scores),
    analysis_name = analysis_name
  ))
}

# ==============================================================================
# 6. Apply Outlier Detection to All Groups
# ==============================================================================

cat("\n==============================================\n")
cat("Outlier Detection Analysis\n")
cat("==============================================\n")

outlier_detection_results <- list()

for (analysis_name in names(pca_results)) {
  if (!is.null(pca_results[[analysis_name]])) {
    outlier_result <- enhanced_outlier_detection(
      analysis_name, 
      pca_results[[analysis_name]], 
      verbose = TRUE
    )
    
    if (!is.null(outlier_result)) {
      outlier_detection_results[[analysis_name]] <- outlier_result
    }
  }
}

# ==============================================================================
# 7. Create Filtered Sample Lists with Pair Consistency
# ==============================================================================

cat("\n==============================================\n")
cat("Creating filtered sample lists with pair consistency\n")
cat("==============================================\n")

# Initialize filtered sample lists
pca_filtered_sample_lists <- sample_lists

# First, collect all outliers by group and tissue
outliers_by_group <- list()

for (analysis_name in names(outlier_detection_results)) {
  if (is.null(outlier_detection_results[[analysis_name]])) next
  
  outlier_data <- outlier_detection_results[[analysis_name]]
  if (length(outlier_data$combined_samples) == 0) next
  
  # Parse group and tissue from analysis name
  pca_data <- pca_results[[analysis_name]]
  group_name <- pca_data$group
  tissue_type <- pca_data$tissue
  
  # Initialize group if not exists
  if (!group_name %in% names(outliers_by_group)) {
    outliers_by_group[[group_name]] <- list(tumor = character(0), normal = character(0))
  }
  
  # Add outliers to the appropriate tissue type
  outliers_by_group[[group_name]][[tissue_type]] <- outlier_data$combined_samples
}

# Apply pair consistency for each group
for (group_name in c("R0", "R1", "B0", "B1")) {
  if (!group_name %in% names(pca_filtered_sample_lists)) next
  
  # Get original samples and cases
  orig_tumor <- sample_lists[[group_name]]$tumor
  orig_normal <- sample_lists[[group_name]]$normal
  orig_cases <- sample_lists[[group_name]]$cases
  
  if (length(orig_cases) == 0 || length(orig_tumor) == 0 || length(orig_normal) == 0) {
    cat(sprintf("%s: No paired samples to process\n", group_name))
    next
  }
  
  # Collect outlier indices
  outlier_indices <- c()
  
  if (group_name %in% names(outliers_by_group)) {
    # Find indices of tumor outliers
    if (length(outliers_by_group[[group_name]]$tumor) > 0) {
      tumor_outlier_idx <- which(orig_tumor %in% outliers_by_group[[group_name]]$tumor)
      outlier_indices <- c(outlier_indices, tumor_outlier_idx)
      cat(sprintf("%s: Found %d tumor outlier(s)\n", group_name, length(tumor_outlier_idx)))
    }
    
    # Find indices of normal outliers
    if (length(outliers_by_group[[group_name]]$normal) > 0) {
      normal_outlier_idx <- which(orig_normal %in% outliers_by_group[[group_name]]$normal)
      outlier_indices <- c(outlier_indices, normal_outlier_idx)
      cat(sprintf("%s: Found %d normal outlier(s)\n", group_name, length(normal_outlier_idx)))
    }
  }
  
  # Remove duplicate indices
  outlier_indices <- unique(outlier_indices)
  
  if (length(outlier_indices) > 0) {
    # Get cases to exclude
    exclude_cases <- orig_cases[outlier_indices]
    cat(sprintf("%s: Excluding cases for pair consistency: %s\n", 
                group_name, paste(exclude_cases, collapse = ", ")))
    
    # Keep only non-outlier pairs
    keep_indices <- setdiff(1:length(orig_cases), outlier_indices)
    
    pca_filtered_sample_lists[[group_name]]$tumor <- orig_tumor[keep_indices]
    pca_filtered_sample_lists[[group_name]]$normal <- orig_normal[keep_indices]
    pca_filtered_sample_lists[[group_name]]$cases <- orig_cases[keep_indices]
    
    cat(sprintf("%s: %d pairs → %d pairs after filtering\n",
                group_name, length(orig_cases), length(keep_indices)))
  } else {
    cat(sprintf("%s: No outliers detected, keeping all %d pairs\n", 
                group_name, length(orig_cases)))
  }
}

# ==============================================================================
# 8. Summary Statistics
# ==============================================================================

cat("\n==============================================\n")
cat("PCA Phase 1 Summary (New CDM)\n")
cat("==============================================\n")

# Create summary table
pca_summary <- data.frame(
  Group = character(0),
  Tissue = character(0),
  Original_N = integer(0),
  Outliers = integer(0),
  Final_N = integer(0),
  PC1_Var = numeric(0),
  PC2_Var = numeric(0),
  stringsAsFactors = FALSE
)

for (analysis_name in names(analysis_groups)) {
  group_info <- analysis_groups[[analysis_name]]
  
  if (analysis_name %in% names(pca_results) && !is.null(pca_results[[analysis_name]])) {
    pca_data <- pca_results[[analysis_name]]
    
    n_outliers <- 0
    if (analysis_name %in% names(outlier_detection_results)) {
      n_outliers <- length(outlier_detection_results[[analysis_name]]$combined_outliers)
    }
    
    pc1_var <- pca_data$pca$variance_explained[1] * 100
    pc2_var <- ifelse(length(pca_data$pca$variance_explained) > 1, 
                      pca_data$pca$variance_explained[2] * 100, 0)
    
    pca_summary <- rbind(pca_summary, data.frame(
      Group = group_info$group,
      Tissue = group_info$tissue,
      Original_N = pca_data$n_samples,
      Outliers = n_outliers,
      Final_N = pca_data$n_samples - n_outliers,
      PC1_Var = round(pc1_var, 1),
      PC2_Var = round(pc2_var, 1),
      stringsAsFactors = FALSE
    ))
  }
}

print(pca_summary)

# Final pair counts after consistency filtering
cat("\nFinal paired sample counts after outlier removal:\n")
for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(pca_filtered_sample_lists)) {
    n_pairs <- length(pca_filtered_sample_lists[[group]]$cases)
    cat(sprintf("  %s: %d pairs\n", group, n_pairs))
  }
}

# ==============================================================================
# 9. Save Results
# ==============================================================================

cat("\nSaving PCA Phase 1 results with new CDM...\n")

# Create processed directory if it doesn't exist
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save comprehensive results
pca_phase1_results_v3 <- list(
  pca_results = pca_results,
  outlier_detection_results = outlier_detection_results,
  pca_filtered_sample_lists = pca_filtered_sample_lists,
  pca_summary = pca_summary,
  analysis_groups = analysis_groups,
  analysis_date = Sys.time(),
  analysis_version = "v3_new_cdm_with_pair_consistency",
  cdm_implementation = "CDM_fast3_arma_enhanced"
)

save(pca_phase1_results_v3, file = "./data/processed/pca_phase1_results_v3.rda")
save(pca_filtered_sample_lists, file = "./data/processed/pca_filtered_sample_lists_v3.rda")

cat("Results saved to ./data/processed/\n")

# ==============================================================================
# 10. Final Summary
# ==============================================================================

cat("\n==============================================\n")
cat("PCA Phase 1 Analysis v3 Complete\n")
cat("==============================================\n")

total_outliers <- sum(pca_summary$Outliers)
cat(sprintf("Total outliers detected: %d\n", total_outliers))

# Check minimum sample requirements
min_samples <- 5
viable_groups <- character(0)
for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(pca_filtered_sample_lists)) {
    n_pairs <- length(pca_filtered_sample_lists[[group]]$cases)
    if (n_pairs >= min_samples) {
      viable_groups <- c(viable_groups, group)
    }
  }
}

cat(sprintf("\nGroups viable for ContamDE analysis (≥%d pairs): %s\n",
            min_samples, paste(viable_groups, collapse = ", ")))

cat("\nNext steps:\n")
cat("1. ContamDE analysis with PCA-filtered samples (06_contamde_analysis_v3.R)\n")
cat("2. Continue with downstream analysis\n")
cat("3. Compare results with old CDM to assess impact\n")

cat("\nPCA Phase 1 with new CDM completed successfully!\n")
cat("==============================================\n")

