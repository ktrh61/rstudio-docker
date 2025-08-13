# ==============================================================================
# REBC-THYR PCA Phase 2 Analysis Script v2 - Group Comparison Visualization
# 07_pca_phase2_v2.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(dplyr)

# Source the CDM_fast function
source("./utils/CDM_fast.R")

cat("Starting PCA Phase 2 analysis v2 for group comparison visualization...\n")

# ==============================================================================
# 1. Load Data and Final Filtered Sample Lists
# ==============================================================================

cat("Loading data and final filtered sample lists...\n")

# Load final filtered sample lists from ContamDE v2
load("./data/processed/final_filtered_sample_lists.rda")

# Load comparison info
load("./data/processed/comparison_info_v2.rda")

# Load TPM data
load("./data/processed/thyr_tpm.rda")

cat("TPM data dimensions:", dim(tpm), "\n")
cat("Available groups:", names(final_filtered_sample_lists), "\n")

# Display final sample counts
cat("Final sample counts after dual filtering:\n")
for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(final_filtered_sample_lists)) {
    tumor_count <- length(final_filtered_sample_lists[[group]]$tumor)
    normal_count <- length(final_filtered_sample_lists[[group]]$normal)
    cat(sprintf("  %s: Tumor=%d, Normal=%d\n", group, tumor_count, normal_count))
  }
}

# ==============================================================================
# 2. Define Phase 2 Comparison Groups
# ==============================================================================

cat("Setting up Phase 2 comparison groups for group visualization...\n")

# Define the 4 comparison groups for Phase 2
comparison_groups <- list(
  "R0_R1_normal" = list(
    groups = c("R0", "R1"),
    tissue = "normal",
    description = "RET Normal: Unexposed vs RadHigh",
    title = "RET Normal Tissue (R0 vs R1)"
  ),
  "R0_R1_tumor" = list(
    groups = c("R0", "R1"),
    tissue = "tumor", 
    description = "RET Tumor: Unexposed vs RadHigh",
    title = "RET Tumor Tissue (R0 vs R1)"
  ),
  "B0_B1_normal" = list(
    groups = c("B0", "B1"),
    tissue = "normal",
    description = "BRAF Normal: Unexposed vs RadHigh", 
    title = "BRAF Normal Tissue (B0 vs B1)"
  ),
  "B0_B1_tumor" = list(
    groups = c("B0", "B1"),
    tissue = "tumor",
    description = "BRAF Tumor: Unexposed vs RadHigh",
    title = "BRAF Tumor Tissue (B0 vs B1)"
  )
)

cat("Phase 2 comparison groups defined:", length(comparison_groups), "\n")

# ==============================================================================
# 3. Phase 2: Integrated PCA Analysis for Group Comparison Visualization
# ==============================================================================

cat("Performing Phase 2 integrated PCA analysis...\n")

# Initialize storage for Phase 2 results
phase2_results <- list()

for (comp_name in names(comparison_groups)) {
  cat(sprintf("\n--- Phase 2 Analysis: %s ---\n", comp_name))
  
  comp_info <- comparison_groups[[comp_name]]
  groups <- comp_info$groups
  tissue <- comp_info$tissue
  description <- comp_info$description
  
  # Collect samples from both groups using final filtered lists
  all_samples <- c()
  all_labels <- c()
  
  for (group in groups) {
    if (group %in% names(final_filtered_sample_lists)) {
      group_samples <- final_filtered_sample_lists[[group]][[tissue]]
      
      if (length(group_samples) > 0) {
        all_samples <- c(all_samples, group_samples)
        all_labels <- c(all_labels, rep(group, length(group_samples)))
      }
    }
  }
  
  if (length(all_samples) < 6) {
    cat(sprintf("Insufficient samples (%d) for %s, skipping...\n", length(all_samples), comp_name))
    next
  }
  
  cat(sprintf("%s: %d samples total\n", description, length(all_samples)))
  
  # Display sample breakdown
  sample_counts <- table(all_labels)
  for (group in names(sample_counts)) {
    cat(sprintf("  %s: %d samples\n", group, sample_counts[group]))
  }
  
  # Extract TPM data for combined samples
  available_samples <- all_samples[all_samples %in% colnames(tpm)]
  if (length(available_samples) < length(all_samples)) {
    missing_count <- length(all_samples) - length(available_samples)
    cat(sprintf("Warning: %d samples missing from TPM data\n", missing_count))
  }
  
  if (length(available_samples) < 6) {
    cat(sprintf("Insufficient available samples (%d), skipping...\n", length(available_samples)))
    next
  }
  
  # Update labels to match available samples
  available_labels <- all_labels[all_samples %in% colnames(tpm)]
  
  # Extract TPM data
  combined_tpm <- tpm[, available_samples, drop = FALSE]
  
  # Minimal filtering (same as Phase 1)
  non_zero_genes <- rowSums(combined_tpm > 0) > 0
  combined_tpm_filtered <- combined_tpm[non_zero_genes, ]
  
  cat(sprintf("Using %d genes for PCA (removed %d all-zero genes)\n", 
              nrow(combined_tpm_filtered), sum(!non_zero_genes)))
  
  # Run CDM_fast PCA
  cat("Running integrated CDM_fast PCA...\n")
  
  tryCatch({
    pca_result <- CDM_fast(
      X = combined_tpm_filtered,
      by_sample = FALSE,
      center = TRUE,
      scale. = TRUE,
      k = min(10, ncol(combined_tpm_filtered) - 1),
      return_scores = TRUE,
      verbose = FALSE
    )
    
    # Calculate group separation distance
    unique_groups <- unique(available_labels)
    if (length(unique_groups) == 2) {
      group1_scores <- pca_result$scores[available_labels == unique_groups[1], 1:2, drop = FALSE]
      group2_scores <- pca_result$scores[available_labels == unique_groups[2], 1:2, drop = FALSE]
      
      group1_centroid <- colMeans(group1_scores)
      group2_centroid <- colMeans(group2_scores)
      
      # Euclidean distance between centroids
      separation_distance <- sqrt(sum((group1_centroid - group2_centroid)^2))
    } else {
      separation_distance <- NA
    }
    
    # Store results
    phase2_results[[comp_name]] <- list(
      pca = pca_result,
      samples = available_samples,
      labels = available_labels,
      groups = groups,
      tissue = tissue,
      description = description,
      title = comp_info$title,
      n_genes = nrow(combined_tpm_filtered),
      n_samples = length(available_samples),
      separation_distance = separation_distance,
      sample_counts = as.list(table(available_labels))
    )
    
    cat(sprintf("Integrated PCA completed: %d components, %.1f%% variance explained by PC1-2\n",
                pca_result$n_components,
                sum(pca_result$variance_explained[1:2]) * 100))
    
    if (!is.na(separation_distance)) {
      cat(sprintf("Group separation distance: %.3f\n", separation_distance))
    }
    
  }, error = function(e) {
    cat(sprintf("Error in integrated PCA for %s: %s\n", comp_name, e$message))
    phase2_results[[comp_name]] <- NULL
  })
}

# ==============================================================================
# 4. Phase 2: Generate Comparison Plots and Statistics
# ==============================================================================

cat("\nGenerating Phase 2 comparison plots and statistics...\n")

# Create plots directory if it doesn't exist
if (!dir.exists("./output/plots")) {
  dir.create("./output/plots", recursive = TRUE)
}

phase2_statistics <- list()

for (comp_name in names(phase2_results)) {
  if (is.null(phase2_results[[comp_name]])) next
  
  result <- phase2_results[[comp_name]]
  scores <- result$pca$scores
  labels <- result$labels
  groups <- result$groups
  title <- result$title
  
  if (is.null(scores)) {
    cat(sprintf("No scores available for %s\n", comp_name))
    next
  }
  
  cat(sprintf("Creating visualization and statistics for %s...\n", comp_name))
  
  # Create color mapping
  colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")  # Professional colors
  unique_groups <- unique(labels)
  color_map <- setNames(colors[1:length(unique_groups)], unique_groups)
  point_colors <- color_map[labels]
  
  # Enhanced plot with better visualization
  tryCatch({
    # Set up plot parameters
    pc1_var <- result$pca$variance_explained[1] * 100
    pc2_var <- result$pca$variance_explained[2] * 100
    
    # Create the plot
    plot(scores[, 1], scores[, 2],
         col = point_colors,
         pch = 19,
         cex = 1.2,
         xlab = sprintf("PC1 (%.1f%%)", pc1_var),
         ylab = sprintf("PC2 (%.1f%%)", pc2_var),
         main = title,
         cex.main = 1.2,
         cex.lab = 1.1)
    
    # Add legend
    legend("topright", 
           legend = unique_groups,
           col = color_map[unique_groups],
           pch = 19,
           cex = 1.0,
           bg = "white")
    
    # Add grid for better readability
    grid(lty = 3, col = "lightgray")
    
    # Add separation distance as text if available
    if (!is.na(result$separation_distance)) {
      mtext(sprintf("Separation Distance: %.3f", result$separation_distance), 
            side = 3, line = 0.5, cex = 0.8, col = "gray30")
    }
    
    cat(sprintf("Plot created successfully for %s\n", comp_name))
    
  }, error = function(e) {
    cat(sprintf("Error creating plot for %s: %s\n", comp_name, e$message))
  })
  
  # Generate detailed statistics
  cat(sprintf("\nDetailed statistics for %s:\n", comp_name))
  
  stats <- list(
    comparison = comp_name,
    description = result$description,
    total_samples = result$n_samples,
    n_genes = result$n_genes,
    pc1_variance = pc1_var,
    pc2_variance = pc2_var,
    cumulative_variance = pc1_var + pc2_var,
    separation_distance = result$separation_distance,
    sample_counts = result$sample_counts
  )
  
  # Group-specific statistics
  group_stats <- list()
  for (group in unique_groups) {
    group_indices <- which(labels == group)
    n_samples <- length(group_indices)
    pc1_mean <- mean(scores[group_indices, 1])
    pc2_mean <- mean(scores[group_indices, 2])
    pc1_sd <- sd(scores[group_indices, 1])
    pc2_sd <- sd(scores[group_indices, 2])
    
    group_stats[[group]] <- list(
      n_samples = n_samples,
      pc1_mean = pc1_mean,
      pc2_mean = pc2_mean,
      pc1_sd = pc1_sd,
      pc2_sd = pc2_sd
    )
    
    cat(sprintf("  %s: n=%d, PC1=%.2f±%.2f, PC2=%.2f±%.2f\n", 
                group, n_samples, pc1_mean, pc1_sd, pc2_mean, pc2_sd))
  }
  
  stats$group_statistics <- group_stats
  
  # Statistical significance assessment (informal)
  if (length(unique_groups) == 2 && !is.na(result$separation_distance)) {
    # Simple effect size assessment
    total_variance <- var(scores[, 1]) + var(scores[, 2])
    effect_size <- result$separation_distance / sqrt(total_variance)
    
    stats$effect_size <- effect_size
    
    if (effect_size > 2) {
      separation_quality <- "Excellent"
    } else if (effect_size > 1) {
      separation_quality <- "Good"
    } else if (effect_size > 0.5) {
      separation_quality <- "Moderate"
    } else {
      separation_quality <- "Poor"
    }
    
    stats$separation_quality <- separation_quality
    cat(sprintf("  Separation quality: %s (effect size: %.2f)\n", separation_quality, effect_size))
  }
  
  phase2_statistics[[comp_name]] <- stats
}

# ==============================================================================
# 5. Phase 2 Summary Table
# ==============================================================================

cat("\n==============================================\n")
cat("Phase 2 PCA Analysis Summary\n")
cat("==============================================\n")

# Create comprehensive summary table
phase2_summary <- data.frame(
  Comparison = character(0),
  Driver = character(0),
  Tissue = character(0),
  Group1_N = integer(0),
  Group2_N = integer(0),
  Total_N = integer(0),
  PC1_Var = numeric(0),
  PC2_Var = numeric(0),
  Cumulative_Var = numeric(0),
  Separation_Distance = numeric(0),
  Quality = character(0),
  stringsAsFactors = FALSE
)

for (comp_name in names(phase2_statistics)) {
  stats <- phase2_statistics[[comp_name]]
  
  # Parse comparison name for driver
  driver <- ifelse(grepl("R0_R1", comp_name), "RET", "BRAF")
  tissue <- stats$sample_counts
  group_names <- names(tissue)
  
  phase2_summary <- rbind(phase2_summary, data.frame(
    Comparison = comp_name,
    Driver = driver,
    Tissue = stringr::str_extract(comp_name, "(normal|tumor)"),
    Group1_N = as.integer(tissue[[group_names[1]]]),
    Group2_N = as.integer(tissue[[group_names[2]]]),
    Total_N = stats$total_samples,
    PC1_Var = round(stats$pc1_variance, 1),
    PC2_Var = round(stats$pc2_variance, 1),
    Cumulative_Var = round(stats$cumulative_variance, 1),
    Separation_Distance = round(stats$separation_distance, 3),
    Quality = ifelse(is.null(stats$separation_quality), "N/A", stats$separation_quality),
    stringsAsFactors = FALSE
  ))
}

print(phase2_summary)

# Identify best separations
cat("\nBest group separations (by distance):\n")
if (nrow(phase2_summary) > 0) {
  ordered_summary <- phase2_summary[order(phase2_summary$Separation_Distance, decreasing = TRUE), ]
  for (i in 1:min(3, nrow(ordered_summary))) {
    row <- ordered_summary[i, ]
    cat(sprintf("  %d. %s: %.3f (%s quality)\n", 
                i, row$Comparison, row$Separation_Distance, row$Quality))
  }
}

# ==============================================================================
# 6. Save Results
# ==============================================================================

cat("Saving Phase 2 PCA analysis results...\n")

# Create processed directory if it doesn't exist
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save comprehensive Phase 2 results
pca_phase2_results <- list(
  phase2_results = phase2_results,
  comparison_groups = comparison_groups,
  phase2_statistics = phase2_statistics,
  phase2_summary = phase2_summary,
  
  analysis_date = Sys.time(),
  analysis_version = "v2_phase2_post_dual_filtering"
)

save(pca_phase2_results, file = "./data/processed/pca_phase2_results.rda")

cat("Results saved to ./data/processed/\n")

# ==============================================================================
# 7. Final Summary and Recommendations
# ==============================================================================

cat("\n==============================================\n")
cat("Phase 2 PCA Analysis Final Summary\n")
cat("==============================================\n")

cat("Phase 2 - Group Comparison Analysis:\n")
cat(sprintf("Successfully analyzed %d group comparisons:\n", nrow(phase2_summary)))

for (i in 1:nrow(phase2_summary)) {
  row <- phase2_summary[i, ]
  cat(sprintf("  %s: %d vs %d samples, %.1f%% variance, %s separation\n",
              row$Comparison, row$Group1_N, row$Group2_N, 
              row$Cumulative_Var, row$Quality))
}

# Recommendations for DEG analysis
cat("\nRecommendations for DEG analysis:\n")

excellent_comparisons <- phase2_summary$Comparison[phase2_summary$Quality == "Excellent"]
good_comparisons <- phase2_summary$Comparison[phase2_summary$Quality == "Good"]

if (length(excellent_comparisons) > 0) {
  cat("🥇 Priority 1 (Excellent separation):\n")
  for (comp in excellent_comparisons) {
    desc <- phase2_statistics[[comp]]$description
    cat(sprintf("   - %s\n", desc))
  }
}

if (length(good_comparisons) > 0) {
  cat("🥈 Priority 2 (Good separation):\n")
  for (comp in good_comparisons) {
    desc <- phase2_statistics[[comp]]$description
    cat(sprintf("   - %s\n", desc))
  }
}

# Overall assessment
total_excellent <- length(excellent_comparisons)
total_good <- length(good_comparisons)
total_analyzed <- nrow(phase2_summary)

cat(sprintf("\nOverall quality assessment:\n"))
cat(sprintf("  Excellent separations: %d/%d\n", total_excellent, total_analyzed))
cat(sprintf("  Good+ separations: %d/%d\n", total_excellent + total_good, total_analyzed))

if (total_excellent > 0) {
  cat("✅ Ready for high-confidence DEG analysis\n")
} else if (total_good > 0) {
  cat("✅ Ready for standard DEG analysis\n")
} else {
  cat("⚠️ May need careful interpretation in DEG analysis\n")
}

cat("\nNext steps:\n")
cat("1. Review Phase 2 plots for visual confirmation\n")
cat("2. Proceed to DEG analysis with priority comparisons\n")
cat("3. Focus on RET comparisons if they show better separation\n")

cat("Phase 2 PCA analysis (v2) completed successfully!\n")
cat("==============================================\n")

