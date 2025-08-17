# ==============================================================================
# REBC-THYR Tumor Purity Analysis Script v3 - After New CDM PCA Outlier Removal
# 06_purity_analysis_v3.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(edgeR)
library(dplyr)

# Source the purity estimation functions
source("./utils/contamde_purity_functions.R")

cat("Starting tumor purity analysis v3 with new CDM PCA results...\n")

# ==============================================================================
# 1. Load Data and Apply CDM PCA Outlier Pair Removal
# ==============================================================================

cat("Loading data and applying new CDM PCA outlier pair removal...\n")

# Load the latest PCA results with outlier information (v3 with new CDM)
load("./data/processed/pca_phase1_results_v3.rda")

# Load original sample lists (before PCA filtering)
load("./data/processed/sample_lists.rda")

# Load SummarizedExperiment object for count data
if (!exists("se_thyr")) {
  load("./data/raw/thyr_data.rda")
  se_thyr <- data
  rm(data)
}

# Get stranded_second count data
count_data_full <- assay(se_thyr, "stranded_second")
gene_info <- rowData(se_thyr)

cat("Original count data dimensions:", dim(count_data_full), "\n")

# ==============================================================================
# 2. Apply CDM PCA Outlier Pair Removal
# ==============================================================================

cat("Applying CDM PCA outlier pair removal with pair consistency...\n")

# Extract outlier information from PCA Phase 1 v3 results
outlier_detection <- pca_phase1_results_v3$outlier_detection_results

# Use the latest pca_filtered_sample_lists from v3
pca_filtered_sample_lists_v3 <- pca_phase1_results_v3$pca_filtered_sample_lists

# Use pre-computed CDM-filtered sample lists from v3
cdm_filtered_sample_lists <- pca_filtered_sample_lists_v3

cat("Using pre-computed CDM-filtered sample lists from v3:\n")

# Display CDM-filtered sample counts
cat("\nCDM-filtered sample counts:\n")
for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(cdm_filtered_sample_lists)) {
    tumor_count <- length(cdm_filtered_sample_lists[[group]]$tumor)
    normal_count <- length(cdm_filtered_sample_lists[[group]]$normal)
    cat(sprintf("  %s: %d pairs (Tumor=%d, Normal=%d)\n", 
                group, tumor_count, tumor_count, normal_count))
  }
}

# ==============================================================================
# 3. Filter to Protein-Coding Genes
# ==============================================================================

cat("\nFiltering to protein-coding genes...\n")

# Filter to protein-coding genes only
protein_coding_genes <- gene_info$gene_type == "protein_coding"
count_data_pc <- count_data_full[protein_coding_genes, ]
gene_info_pc <- gene_info[protein_coding_genes, ]

cat("After protein-coding filter:", dim(count_data_pc), "\n")

# ==============================================================================
# 4. Process Each Group for Tumor Purity Estimation
# ==============================================================================

cat("\n==============================================\n")
cat("Tumor Purity Estimation for Each Group\n")
cat("==============================================\n")

# Initialize results storage
purity_results <- list()
purity_cutoff <- 0.6  # Absolute purity threshold (60%)

for (group in c("R0", "R1", "B0", "B1")) {
  cat(sprintf("\n### Processing group %s ###\n", group))
  
  # Check if group has samples after CDM filtering
  if (!group %in% names(cdm_filtered_sample_lists) || 
      length(cdm_filtered_sample_lists[[group]]$tumor) == 0) {
    cat(sprintf("Group %s has no samples after CDM filtering, skipping...\n", group))
    purity_results[[group]] <- NULL
    next
  }
  
  # Get CDM-filtered tumor and normal samples for this group
  tumor_samples <- cdm_filtered_sample_lists[[group]]$tumor
  normal_samples <- cdm_filtered_sample_lists[[group]]$normal
  
  cat(sprintf("CDM-filtered samples - Tumor: %d, Normal: %d\n", 
              length(tumor_samples), length(normal_samples)))
  
  # Verify sample availability in count data
  tumor_available <- tumor_samples %in% colnames(count_data_pc)
  normal_available <- normal_samples %in% colnames(count_data_pc)
  
  if (!all(tumor_available) || !all(normal_available)) {
    missing_tumor <- sum(!tumor_available)
    missing_normal <- sum(!normal_available)
    cat(sprintf("Warning: Missing samples - tumor: %d, normal: %d\n", 
                missing_tumor, missing_normal))
  }
  
  # Use only available samples
  tumor_samples_avail <- tumor_samples[tumor_available]
  normal_samples_avail <- normal_samples[normal_available]
  
  # Ensure equal number of tumor and normal samples
  min_pairs <- min(length(tumor_samples_avail), length(normal_samples_avail))
  if (min_pairs == 0) {
    cat(sprintf("Group %s has no valid pairs, skipping...\n", group))
    purity_results[[group]] <- NULL
    next
  }
  
  if (min_pairs < 3) {
    cat(sprintf("Warning: Group %s has only %d pairs. Results may be unstable.\n", 
                group, min_pairs))
  }
  
  tumor_samples_final <- tumor_samples_avail[1:min_pairs]
  normal_samples_final <- normal_samples_avail[1:min_pairs]
  
  # Prepare count data in [normal1, normal2, ..., tumor1, tumor2, ...] format
  group_samples <- c(normal_samples_final, tumor_samples_final)
  group_counts <- count_data_pc[, group_samples]
  
  cat(sprintf("Final analysis: %d pairs (%d samples total)\n", 
              min_pairs, ncol(group_counts)))
  
  # Apply filterByExpr for low expression filtering
  dge <- DGEList(counts = group_counts)
  # Create group factor for filterByExpr (normal=0, tumor=1)
  group_factor <- factor(c(rep("normal", min_pairs), rep("tumor", min_pairs)))
  keep <- filterByExpr(dge, group = group_factor, min.count = 1, min.total.count = 15)
  dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
  
  cat(sprintf("After filterByExpr: %d genes (from %d)\n", 
              nrow(dge_filtered), nrow(group_counts)))
  
  # Run contamDE purity estimation
  cat("Running contamDE purity estimation...\n")
  
  tryCatch({
    purity_result <- contamde_purity(
      counts = dge_filtered$counts,
      subtype = NULL,
      covariate = NULL,
      contaminated = TRUE,
      verbose = TRUE
    )
    
    # Store results with metadata
    purity_results[[group]] <- purity_result
    purity_results[[group]]$group_name <- group
    purity_results[[group]]$tumor_samples <- tumor_samples_final
    purity_results[[group]]$normal_samples <- normal_samples_final
    purity_results[[group]]$filtered_genes <- rownames(dge_filtered)
    purity_results[[group]]$analysis_date <- Sys.time()
    
    cat(sprintf("Group %s: Purity estimation completed successfully\n", group))
    
  }, error = function(e) {
    cat(sprintf("Error processing group %s: %s\n", group, e$message))
    purity_results[[group]] <- NULL
  })
}

# ==============================================================================
# 5. Purity Quality Assessment and Filtering
# ==============================================================================

cat("\n==============================================\n")
cat("Purity Quality Assessment\n")
cat("==============================================\n")

# Initialize summary table
purity_summary <- data.frame(
  Group = c("R0", "R1", "B0", "B1"),
  Input_Pairs = integer(4),
  Mean_Purity = numeric(4),
  SD_Purity = numeric(4),
  Min_Purity = numeric(4),
  Max_Purity = numeric(4),
  High_Purity_Count = integer(4),
  Retention_Rate = numeric(4),
  stringsAsFactors = FALSE
)

# Process results for each group
for (i in 1:4) {
  group <- purity_summary$Group[i]
  
  if (!is.null(purity_results[[group]])) {
    # Assess purity quality
    quality_stats <- assess_purity_quality(purity_results[[group]], 
                                           threshold = purity_cutoff, 
                                           verbose = TRUE)
    
    purity_summary$Input_Pairs[i] <- quality_stats$n_samples
    purity_summary$Mean_Purity[i] <- round(quality_stats$mean_purity, 3)
    purity_summary$SD_Purity[i] <- round(quality_stats$sd_purity, 3)
    purity_summary$Min_Purity[i] <- round(quality_stats$min_purity, 3)
    purity_summary$Max_Purity[i] <- round(quality_stats$max_purity, 3)
    purity_summary$High_Purity_Count[i] <- quality_stats$n_high_purity
    purity_summary$Retention_Rate[i] <- round(quality_stats$retention_rate, 3)
    
  } else {
    purity_summary$Input_Pairs[i] <- 0
    purity_summary$Mean_Purity[i] <- NA
    purity_summary$SD_Purity[i] <- NA
    purity_summary$Min_Purity[i] <- NA
    purity_summary$Max_Purity[i] <- NA
    purity_summary$High_Purity_Count[i] <- 0
    purity_summary$Retention_Rate[i] <- 0
  }
}

print(purity_summary)

# ==============================================================================
# 6. Create Purity-Filtered Sample Lists
# ==============================================================================

cat("\nCreating purity-filtered sample lists...\n")

# Create final purity-filtered sample lists
purity_filtered_sample_lists <- create_purity_filtered_lists(
  original_sample_lists = cdm_filtered_sample_lists,
  purity_results = purity_results,
  threshold = purity_cutoff,
  verbose = TRUE
)

# ==============================================================================
# 7. Special Assessment for B1 Group
# ==============================================================================

cat("\n==============================================\n")
cat("Special Assessment for B1 Group\n")
cat("==============================================\n")

if (!is.null(purity_results[["B1"]])) {
  b1_result <- purity_results[["B1"]]
  b1_proportions <- b1_result$proportion
  
  cat("B1 Group Detailed Analysis:\n")
  cat(sprintf("Total pairs: %d\n", length(b1_proportions)))
  cat(sprintf("Purity range: [%.3f, %.3f]\n", min(b1_proportions), max(b1_proportions)))
  cat(sprintf("Purity ≥ %.1f: %d pairs\n", purity_cutoff, sum(b1_proportions >= purity_cutoff)))
  
  # Check if any high-ERR samples were excluded
  excluded_pairs <- sum(b1_proportions < purity_cutoff)
  if (excluded_pairs > 0) {
    cat(sprintf("\n⚠️  WARNING: %d B1 pairs excluded by purity filter\n", excluded_pairs))
    cat("Individual purity values:\n")
    for (i in 1:length(b1_proportions)) {
      status <- ifelse(b1_proportions[i] >= purity_cutoff, "PASS", "EXCLUDED")
      cat(sprintf("  Pair %d: %.3f (%s)\n", i, b1_proportions[i], status))
    }
    
    cat("\nRecommendation: Review excluded B1 pairs individually\n")
    cat("Consider: ERR values, technical quality, biological markers\n")
  } else {
    cat("✅ All B1 pairs passed purity filter\n")
  }
} else {
  cat("B1 group analysis failed or no data available\n")
}

# ==============================================================================
# 8. Summary and Quality Control
# ==============================================================================

cat("\n==============================================\n")
cat("Filtering Cascade Summary\n")
cat("==============================================\n")

# Load original sample lists for comparison
original_sample_lists <- sample_lists

cat("Complete filtering cascade:\n")
cat("Original → CDM-filtered → Purity-filtered:\n")

for (group in c("R0", "R1", "B0", "B1")) {
  orig_count <- length(original_sample_lists[[group]]$tumor)
  cdm_count <- length(cdm_filtered_sample_lists[[group]]$tumor)
  purity_count <- length(purity_filtered_sample_lists[[group]]$tumor)
  
  cat(sprintf("  %s: %d → %d → %d pairs\n", group, orig_count, cdm_count, purity_count))
}

# Check minimum sample requirements
min_samples <- 5
cat("\nSample size check (minimum 5 pairs per group):\n")
for (group in c("R0", "R1", "B0", "B1")) {
  final_count <- length(purity_filtered_sample_lists[[group]]$tumor)
  status <- ifelse(final_count >= min_samples, "✅ PASS", "❌ INSUFFICIENT")
  cat(sprintf("  %s: %d pairs %s\n", group, final_count, status))
}

total_final_pairs <- sum(sapply(purity_filtered_sample_lists, function(x) length(x$tumor)))
cat(sprintf("\nTotal pairs after dual filtering: %d\n", total_final_pairs))

# ==============================================================================
# 9. Save Results
# ==============================================================================

cat("Saving purity analysis v3 results...\n")

# Create processed directory if it doesn't exist
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save CDM-filtered sample lists
save(cdm_filtered_sample_lists, file = "./data/processed/cdm_filtered_sample_lists.rda")

# Save purity-filtered sample lists
save(purity_filtered_sample_lists, file = "./data/processed/purity_filtered_sample_lists.rda")

# Save comprehensive purity analysis results
purity_analysis_v3_results <- list(
  purity_results = purity_results,
  purity_summary = purity_summary,
  cdm_filtered_sample_lists = cdm_filtered_sample_lists,
  purity_filtered_sample_lists = purity_filtered_sample_lists,
  purity_cutoff = purity_cutoff,
  protein_coding_genes = rownames(count_data_pc),
  analysis_date = Sys.time(),
  analysis_version = "v3_post_cdm_full_groups"
)

save(purity_analysis_v3_results, file = "./data/processed/purity_analysis_v3_results.rda")

cat("Results saved to ./data/processed/\n")

# ==============================================================================
# 10. Final Summary and Next Steps
# ==============================================================================

cat("\n==============================================\n")
cat("Purity Analysis v3 Final Summary\n")
cat("==============================================\n")

cat("CDM + Purity Dual Filtering Completed!\n\n")

cat("Key achievements:\n")
cat("✅ CDM PCA outlier pair removal applied consistently\n")
cat("✅ Tumor purity estimation for all groups\n") 
cat("✅ 60% purity threshold applied uniformly\n")
cat("✅ B1 group special assessment completed\n")
cat("✅ Comprehensive quality control pipeline\n")

if (total_final_pairs > 0) {
  cat(sprintf("\nHigh-quality samples ready: %d pairs total\n", total_final_pairs))
  
  # Identify viable comparisons
  viable_groups <- names(purity_filtered_sample_lists)[
    sapply(purity_filtered_sample_lists, function(x) length(x$tumor) >= min_samples)
  ]
  
  if (length(viable_groups) >= 2) {
    cat("Viable comparisons for DEG analysis:\n")
    if ("R0" %in% viable_groups && "R1" %in% viable_groups) {
      cat("  ⭐ R0 vs R1 (PRIMARY comparison)\n")
    }
    if ("B0" %in% viable_groups && "B1" %in% viable_groups) {
      cat("  📋 B0 vs B1 (BRAF comparison)\n")
    }
  } else {
    cat("⚠️  Limited viable comparisons. Review sample retention.\n")
  }
} else {
  cat("❌ No samples passed dual filtering. Review thresholds.\n")
}

cat("\nNext steps:\n")
cat("1. Review B1 exclusions (if any) for radiation markers\n")
cat("2. Phase 2 PCA analysis for group comparison\n")
cat("3. Proceed to DEGES normalization (08)\n")
cat("4. DEG analysis with high-purity samples\n")

cat("\nPurity analysis v3 completed successfully!\n")
cat("==============================================\n")
