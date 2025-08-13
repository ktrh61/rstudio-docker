# ==============================================================================
# REBC-THYR PCA Analysis Script v2 - Phase 1 Outlier Detection Only
# 05_pca_analysis_v2.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(dplyr)

# Source the CDM_fast function
source("./utils/CDM_fast.R")

cat("Starting PCA analysis v2 for outlier detection (Phase 1 only)...\n")

# ==============================================================================
# 1. Load Data and Sample Lists
# ==============================================================================

cat("Loading data and sample lists...\n")

# Load sample lists from ERR integration (before ContamDE)
load("./data/processed/sample_lists.rda")

# Load TPM data
load("./data/processed/thyr_tpm.rda")

# Load SummarizedExperiment object for gene information
if (!exists("se_thyr")) {
  load("./data/raw/thyr_data.rda")
  se_thyr <- data
  rm(data)
}

gene_info <- rowData(se_thyr)

cat("TPM data dimensions:", dim(tpm), "\n")
cat("Available groups:", names(sample_lists), "\n")

# ==============================================================================
# 2. Prepare TPM Data for PCA (All Genes)
# ==============================================================================

cat("Preparing TPM data for PCA analysis (using all genes)...\n")

# Use all genes for PCA to preserve TPM relative value properties
cat("TPM data dimensions (all genes):", dim(tpm), "\n")
cat("Using all genes to preserve TPM compositional properties\n")

# ==============================================================================
# 3. Define Analysis Groups (Phase 1 Only - Individual Groups)
# ==============================================================================

cat("Setting up Phase 1 analysis groups (individual groups only)...\n")

# Define the 8 individual analysis groups for outlier detection
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

cat("Phase 1 analysis groups defined:", length(analysis_groups), "\n")

# ==============================================================================
# 4. Perform PCA Analysis for Each Individual Group
# ==============================================================================

cat("Performing Phase 1 PCA analysis for each individual group...\n")

# Initialize storage for results
pca_results <- list()

for (analysis_name in names(analysis_groups)) {
  cat(sprintf("\n--- Analyzing %s ---\n", analysis_name))
  
  group_info <- analysis_groups[[analysis_name]]
  group_name <- group_info$group
  tissue_type <- group_info$tissue
  description <- group_info$description
  
  # Check if group exists and has samples
  if (!group_name %in% names(sample_lists)) {
    cat(sprintf("Group %s not found in sample lists, skipping...\n", group_name))
    next
  }
  
  # Get samples for this group and tissue type
  if (tissue_type == "tumor") {
    sample_ids <- sample_lists[[group_name]]$tumor
  } else {
    sample_ids <- sample_lists[[group_name]]$normal
  }
  
  if (length(sample_ids) == 0) {
    cat(sprintf("No samples for %s, skipping...\n", analysis_name))
    next
  }
  
  cat(sprintf("%s: %d samples\n", description, length(sample_ids)))
  
  # Check sample availability in TPM data
  available_samples <- sample_ids[sample_ids %in% colnames(tpm)]
  if (length(available_samples) < length(sample_ids)) {
    missing_count <- length(sample_ids) - length(available_samples)
    cat(sprintf("Warning: %d samples missing from TPM data\n", missing_count))
  }
  
  if (length(available_samples) < 3) {
    cat(sprintf("Insufficient samples (%d) for PCA analysis, skipping...\n", length(available_samples)))
    next
  }
  
  # Extract TPM data for this group
  group_tpm <- tpm[, available_samples, drop = FALSE]
  
  # Minimal quality filtering: remove genes with all zeros
  non_zero_genes <- rowSums(group_tpm > 0) > 0
  group_tpm_filtered <- group_tpm[non_zero_genes, ]
  
  cat(sprintf("Using %d genes for PCA (removed %d all-zero genes)\n", 
              nrow(group_tpm_filtered), sum(!non_zero_genes)))
  
  # Run CDM_fast PCA
  cat("Running CDM_fast PCA...\n")
  
  tryCatch({
    # CDM_fast expects genes x samples format
    pca_result_temp <- CDM_fast(
      X = group_tpm_filtered,  # genes x samples
      by_sample = FALSE,       # X is already genes x samples
      center = TRUE,
      scale. = TRUE,
      k = min(10, ncol(group_tpm_filtered) - 1),  # Max 10 PCs or n-1
      return_scores = TRUE,
      verbose = FALSE
    )
    
    # Verify that scores match input samples
    expected_samples <- available_samples
    actual_score_names <- rownames(pca_result_temp$scores)
    
    if (!identical(sort(expected_samples), sort(actual_score_names))) {
      cat(sprintf("WARNING: Sample mismatch detected for %s\n", analysis_name))
      cat("Expected samples:", length(expected_samples), "\n")
      cat("Actual score samples:", length(actual_score_names), "\n")
    }
    
    # Store results with verification
    pca_results[[analysis_name]] <- list(
      pca = pca_result_temp,
      samples = available_samples,
      description = description,
      n_genes = nrow(group_tpm_filtered),
      n_samples = length(available_samples),
      analysis_name = analysis_name
    )
    
    cat(sprintf("PCA completed: %d components, %.1f%% variance explained by PC1-2\n",
                pca_result_temp$n_components,
                sum(pca_result_temp$variance_explained[1:2]) * 100))
    
  }, error = function(e) {
    cat(sprintf("Error in PCA for %s: %s\n", analysis_name, e$message))
    pca_results[[analysis_name]] <- NULL
  })
}

# ==============================================================================
# 5. Enhanced Outlier Detection (Mahalanobis + IQR)
# ==============================================================================

cat("\n==============================================\n")
cat("Enhanced Outlier Detection Analysis\n")
cat("==============================================\n")

# Enhanced outlier detection function
enhanced_outlier_detection <- function(analysis_name, pca_data, verbose = TRUE) {
  
  if (verbose) cat(sprintf("\n--- Enhanced outlier detection for %s ---\n", analysis_name))
  
  scores <- pca_data$pca$scores
  if (is.null(scores) || nrow(scores) < 4) {
    if (verbose) cat(sprintf("Insufficient data (%d samples) for outlier detection\n", 
                             ifelse(is.null(scores), 0, nrow(scores))))
    return(NULL)
  }
  
  if (verbose) cat(sprintf("%s: %d samples\n", pca_data$description, nrow(scores)))
  
  # Method 1: Mahalanobis distance (95% threshold)
  n_pcs_for_outlier <- min(3, ncol(scores))
  pc_subset <- scores[, 1:n_pcs_for_outlier, drop = FALSE]
  
  center_point <- colMeans(pc_subset)
  cov_matrix <- cov(pc_subset)
  
  # Check for singular covariance matrix
  if (det(cov_matrix) < .Machine$double.eps) {
    if (verbose) cat("Warning: Singular covariance matrix, using robust estimation\n")
    cov_matrix <- diag(diag(cov_matrix))
  }
  
  tryCatch({
    mahal_dist <- mahalanobis(pc_subset, center_point, cov_matrix)
    threshold <- qchisq(0.95, df = n_pcs_for_outlier)  # 95% threshold
    
    mahal_outliers <- which(mahal_dist > threshold)
    mahal_samples <- rownames(scores)[mahal_outliers]
    
    if (verbose) {
      cat(sprintf("Mahalanobis method (95%% threshold): %d outliers detected\n", length(mahal_outliers)))
      if (length(mahal_outliers) > 0) {
        for (i in mahal_outliers) {
          cat(sprintf("  %s: distance=%.2f\n", rownames(scores)[i], mahal_dist[i]))
        }
      }
    }
    
  }, error = function(e) {
    if (verbose) cat(sprintf("Error in Mahalanobis detection: %s\n", e$message))
    mahal_outliers <- integer(0)
    mahal_samples <- character(0)
    mahal_dist <- rep(NA, nrow(scores))
  })
  
  # Method 2: IQR method (PC1 focused)
  pc1_values <- scores[, 1]
  
  # Calculate IQR-based outliers
  q1 <- quantile(pc1_values, 0.25)
  q3 <- quantile(pc1_values, 0.75)
  iqr <- q3 - q1
  
  lower_fence <- q1 - 1.5 * iqr
  upper_fence <- q3 + 1.5 * iqr
  
  iqr_outliers <- which(pc1_values < lower_fence | pc1_values > upper_fence)
  iqr_samples <- rownames(scores)[iqr_outliers]
  
  if (verbose) {
    cat(sprintf("IQR method (PC1 focused): %d outliers detected\n", length(iqr_outliers)))
    cat(sprintf("  IQR range: [%.2f, %.2f], Outlier fences: [%.2f, %.2f]\n", 
                q1, q3, lower_fence, upper_fence))
    if (length(iqr_outliers) > 0) {
      for (i in iqr_outliers) {
        cat(sprintf("  %s: PC1=%.2f\n", rownames(scores)[i], pc1_values[i]))
      }
    }
  }
  
  # Method 3: Union approach (combined detection)
  combined_outliers <- unique(c(mahal_outliers, iqr_outliers))
  combined_samples <- rownames(scores)[combined_outliers]
  
  if (verbose) {
    cat(sprintf("Combined detection (Union approach): %d outliers\n", length(combined_outliers)))
    cat(sprintf("  Mahalanobis only: %d\n", length(setdiff(mahal_outliers, iqr_outliers))))
    cat(sprintf("  IQR only: %d\n", length(setdiff(iqr_outliers, mahal_outliers))))
    cat(sprintf("  Both methods: %d\n", length(intersect(mahal_outliers, iqr_outliers))))
    
    if (length(combined_outliers) > 0) {
      cat("Final outlier list:\n")
      for (i in combined_outliers) {
        sample_name <- rownames(scores)[i]
        in_mahal <- i %in% mahal_outliers
        in_iqr <- i %in% iqr_outliers
        method_str <- ifelse(in_mahal && in_iqr, "Both", 
                             ifelse(in_mahal, "Mahal", "IQR"))
        cat(sprintf("  %s (%s): PC1=%.2f\n", sample_name, method_str, pc1_values[i]))
      }
    }
  }
  
  # Return comprehensive results
  return(list(
    mahal_distances = mahal_dist,
    mahal_threshold = threshold,
    mahal_outliers = mahal_outliers,
    mahal_samples = mahal_samples,
    
    iqr_q1 = q1,
    iqr_q3 = q3,
    iqr_lower_fence = lower_fence,
    iqr_upper_fence = upper_fence,
    iqr_outliers = iqr_outliers,
    iqr_samples = iqr_samples,
    
    combined_outliers = combined_outliers,
    combined_samples = combined_samples,
    
    n_pcs_used = n_pcs_for_outlier,
    total_samples = nrow(scores),
    analysis_name = analysis_name
  ))
}

# Apply enhanced detection to all analysis groups
enhanced_outlier_detection_results <- list()

for (analysis_name in names(analysis_groups)) {
  cat(sprintf("\n--- Processing %s ---\n", analysis_name))
  
  if (analysis_name %in% names(pca_results) && 
      !is.null(pca_results[[analysis_name]])) {
    
    pca_data <- pca_results[[analysis_name]]
    enhanced_result <- enhanced_outlier_detection(analysis_name, pca_data, verbose = TRUE)
    
    if (!is.null(enhanced_result)) {
      enhanced_outlier_detection_results[[analysis_name]] <- enhanced_result
      cat(sprintf("%s processing completed successfully\n", analysis_name))
    } else {
      cat(sprintf("%s processing failed or insufficient data\n", analysis_name))
    }
  } else {
    cat(sprintf("%s data not available\n", analysis_name))
  }
}

# Create alias for backward compatibility
outlier_detection <- enhanced_outlier_detection_results

# ==============================================================================
# 6. Enhanced Summary and Quality Control
# ==============================================================================

cat("\n==============================================\n")
cat("Enhanced Outlier Detection Summary\n")
cat("==============================================\n")

# Create enhanced summary table
pca_summary <- data.frame(
  Analysis = character(0),
  Group = character(0),
  Tissue = character(0),
  N_Samples = integer(0),
  N_Genes = integer(0),
  PC1_Variance = numeric(0),
  PC2_Variance = numeric(0),
  Outliers = integer(0),
  Final_Samples = integer(0),
  stringsAsFactors = FALSE
)

for (analysis_name in names(analysis_groups)) {
  group_info <- analysis_groups[[analysis_name]]
  
  if (analysis_name %in% names(pca_results) && 
      !is.null(pca_results[[analysis_name]]) &&
      analysis_name %in% names(enhanced_outlier_detection_results)) {
    
    pca_data <- pca_results[[analysis_name]]
    outlier_data <- enhanced_outlier_detection_results[[analysis_name]]
    
    pc1_var <- pca_data$pca$variance_explained[1] * 100
    pc2_var <- pca_data$pca$variance_explained[2] * 100
    n_combined <- length(outlier_data$combined_outliers)
    final_samples <- pca_data$n_samples - n_combined
    
    pca_summary <- rbind(pca_summary, data.frame(
      Analysis = analysis_name,
      Group = group_info$group,
      Tissue = group_info$tissue,
      N_Samples = pca_data$n_samples,
      N_Genes = pca_data$n_genes,
      PC1_Variance = round(pc1_var, 1),
      PC2_Variance = round(pc2_var, 1),
      Outliers = n_combined,
      Final_Samples = final_samples,
      stringsAsFactors = FALSE
    ))
  } else {
    pca_summary <- rbind(pca_summary, data.frame(
      Analysis = analysis_name,
      Group = group_info$group,
      Tissue = group_info$tissue,
      N_Samples = 0,
      N_Genes = 0,
      PC1_Variance = NA,
      PC2_Variance = NA,
      Outliers = 0,
      Final_Samples = 0,
      stringsAsFactors = FALSE
    ))
  }
}

print(pca_summary)

# Check for excessive outlier detection
cat("\nOutlier detection assessment:\n")
excessive_detection <- pca_summary$Outliers / pca_summary$N_Samples > 0.25
if (any(excessive_detection, na.rm = TRUE)) {
  concerning_groups <- pca_summary$Analysis[excessive_detection & !is.na(excessive_detection)]
  cat("⚠️ Groups with high outlier rates (>25%):\n")
  for (group in concerning_groups) {
    row <- pca_summary[pca_summary$Analysis == group, ]
    rate <- round(row$Outliers / row$N_Samples * 100, 1)
    cat(sprintf("  %s: %d/%d (%s%%)\n", group, row$Outliers, row$N_Samples, rate))
  }
} else {
  cat("✅ Outlier detection rates are reasonable (<25% per group)\n")
}

# ==============================================================================
# 7. Create Filtered Sample Lists (Outliers Removed) with Pair Consistency
# ==============================================================================

cat("Creating sample lists with outliers removed (maintaining pair consistency)...\n")

# Create new sample lists excluding outliers
pca_filtered_sample_lists <- sample_lists

for (analysis_name in names(outlier_detection)) {
  if (is.null(outlier_detection[[analysis_name]])) next
  
  outlier_data <- outlier_detection[[analysis_name]]
  if (length(outlier_data$combined_samples) == 0) next
  
  # Parse analysis name to get group and tissue
  group_info <- analysis_groups[[analysis_name]]
  group_name <- group_info$group
  tissue_type <- group_info$tissue
  
  # Remove outliers from the appropriate sample list
  current_samples <- pca_filtered_sample_lists[[group_name]][[tissue_type]]
  outlier_samples <- outlier_data$combined_samples
  
  filtered_samples <- setdiff(current_samples, outlier_samples)
  pca_filtered_sample_lists[[group_name]][[tissue_type]] <- filtered_samples
  
  cat(sprintf("%s %s: removed %d outliers (%d -> %d samples)\n",
              group_name, tissue_type, 
              length(outlier_samples),
              length(current_samples),
              length(filtered_samples)))
}

# Special handling for pair consistency
cat("\nApplying pair consistency corrections...\n")

# For each group, ensure that cases with outliers in either tissue are completely removed
for (group_name in c("R0", "R1", "B0", "B1")) {
  if (!group_name %in% names(pca_filtered_sample_lists)) next
  
  group_samples <- pca_filtered_sample_lists[[group_name]]
  if (length(group_samples$tumor) == 0 || length(group_samples$normal) == 0) next
  
  # Get original sample lists for comparison
  orig_tumor <- sample_lists[[group_name]]$tumor
  orig_normal <- sample_lists[[group_name]]$normal
  orig_cases <- sample_lists[[group_name]]$cases
  
  # Find outlier samples for this group
  tumor_outliers <- c()
  normal_outliers <- c()
  
  tumor_analysis_name <- paste0(group_name, "_tumor")
  normal_analysis_name <- paste0(group_name, "_normal")
  
  if (tumor_analysis_name %in% names(outlier_detection) && 
      !is.null(outlier_detection[[tumor_analysis_name]])) {
    tumor_outliers <- outlier_detection[[tumor_analysis_name]]$combined_samples
  }
  
  if (normal_analysis_name %in% names(outlier_detection) && 
      !is.null(outlier_detection[[normal_analysis_name]])) {
    normal_outliers <- outlier_detection[[normal_analysis_name]]$combined_samples
  }
  
  if (length(tumor_outliers) > 0 || length(normal_outliers) > 0) {
    # Find indices of outlier samples
    normal_outlier_indices <- which(orig_normal %in% normal_outliers)
    tumor_outlier_indices <- which(orig_tumor %in% tumor_outliers)
    
    # Get corresponding case IDs
    normal_outlier_cases <- orig_cases[normal_outlier_indices]
    tumor_outlier_cases <- orig_cases[tumor_outlier_indices]
    
    # Cases to exclude completely (any outlier in either tissue)
    exclude_cases <- unique(c(normal_outlier_cases, tumor_outlier_cases))
    
    if (length(exclude_cases) > 0) {
      cat(sprintf("%s: Cases with outliers to exclude: %s\n", 
                  group_name, paste(exclude_cases, collapse = ", ")))
      
      # Find indices to keep (not in exclude list)
      exclude_indices <- which(orig_cases %in% exclude_cases)
      keep_indices <- setdiff(1:length(orig_cases), exclude_indices)
      
      # Create final sample lists with pair consistency
      final_tumor <- orig_tumor[keep_indices]
      final_normal <- orig_normal[keep_indices]
      final_cases <- orig_cases[keep_indices]
      
      # Update the filtered sample lists
      pca_filtered_sample_lists[[group_name]]$tumor <- final_tumor
      pca_filtered_sample_lists[[group_name]]$normal <- final_normal
      pca_filtered_sample_lists[[group_name]]$cases <- final_cases
      
      cat(sprintf("%s: Final paired samples: %d (removed %d pairs for consistency)\n",
                  group_name, length(final_cases), length(exclude_cases)))
    }
  }
}

# ==============================================================================
# 8. Save Results
# ==============================================================================

cat("Saving PCA Phase 1 analysis results...\n")

# Create processed directory if it doesn't exist
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save comprehensive PCA Phase 1 results
pca_phase1_results <- list(
  # Phase 1 results (individual group analysis)
  pca_results = pca_results,
  outlier_detection = outlier_detection,
  enhanced_outlier_detection = enhanced_outlier_detection_results,
  pca_summary = pca_summary,
  analysis_groups = analysis_groups,
  pca_filtered_sample_lists = pca_filtered_sample_lists,
  
  analysis_date = Sys.time(),
  analysis_version = "v2_phase1_only"
)

save(pca_phase1_results, file = "./data/processed/pca_phase1_results.rda")
save(pca_filtered_sample_lists, file = "./data/processed/pca_filtered_sample_lists.rda")

cat("Results saved to ./data/processed/\n")

# ==============================================================================
# 9. Final Summary
# ==============================================================================

cat("\n==============================================\n")
cat("PCA Phase 1 Analysis Final Summary\n")
cat("==============================================\n")

cat("Phase 1 - Individual Group Analysis:\n")
cat("Analysis completed for", nrow(pca_summary), "groups:\n")

# Summary by group
for (group in c("R0", "R1", "B0", "B1")) {
  group_data <- pca_summary[pca_summary$Group == group, ]
  if (nrow(group_data) > 0) {
    total_samples <- sum(group_data$N_Samples)
    total_outliers <- sum(group_data$Outliers)
    final_samples <- sum(group_data$Final_Samples)
    
    cat(sprintf("  %s: %d -> %d samples (-%d outliers)\n",
                group, total_samples, final_samples, total_outliers))
  }
}

# Check viability for next steps
min_samples <- 5
viable_groups <- pca_summary %>%
  group_by(Group) %>%
  summarise(
    min_final_samples = min(Final_Samples),
    .groups = "drop"
  ) %>%
  filter(min_final_samples >= min_samples) %>%
  pull(Group)

cat("\nGroups viable for ContamDE analysis (≥5 samples per tissue):\n")
if (length(viable_groups) > 0) {
  cat(paste(viable_groups, collapse = ", "), "\n")
} else {
  cat("⚠️  No groups meet minimum sample requirements\n")
}

cat("\nNext steps:\n")
cat("1. ContamDE analysis with PCA-filtered samples (06_contamde_analysis_v2.R)\n")
cat("2. Phase 2 PCA for group comparison visualization (07_pca_phase2_v2.R)\n")
cat("3. DEG analysis with high-quality samples (08_deg_analysis.R)\n")

total_cases <- sum(sapply(pca_filtered_sample_lists, function(x) length(x$tumor)))
cat(sprintf("\nTotal cases ready for ContamDE: %d\n", total_cases))

cat("PCA Phase 1 analysis (v2) completed successfully!\n")
cat("==============================================\n")