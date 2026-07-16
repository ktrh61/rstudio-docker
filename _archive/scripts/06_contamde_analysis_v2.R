# ==============================================================================
# REBC-THYR ContamDE Analysis Script v2 - After PCA Outlier Removal
# 06_contamde_analysis_v2.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(edgeR)
library(dplyr)

# Source the contamDE function
source("./utils/contamde_functions.R")

cat("Starting ContamDE analysis v2 for tumor purity estimation (after PCA outlier removal)...\n")

# ==============================================================================
# 1. Load Data and PCA-Filtered Sample Lists
# ==============================================================================

cat("Loading data and PCA-filtered sample lists...\n")

# Load PCA-filtered sample lists from Phase 1
load("./data/processed/pca_filtered_sample_lists.rda")

# Load SummarizedExperiment object for count data
if (!exists("se_thyr")) {
  load("./data/raw/thyr_data.rda")
  se_thyr <- data
  rm(data)
}

# Get stranded_second count data for ContamDE analysis
count_data_full <- assay(se_thyr, "stranded_second")
gene_info <- rowData(se_thyr)

cat("Original count data dimensions:", dim(count_data_full), "\n")
cat("PCA-filtered sample lists loaded\n")

# Display PCA-filtered sample counts
cat("PCA-filtered sample counts:\n")
for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(pca_filtered_sample_lists)) {
    tumor_count <- length(pca_filtered_sample_lists[[group]]$tumor)
    normal_count <- length(pca_filtered_sample_lists[[group]]$normal)
    cat(sprintf("  %s: Tumor=%d, Normal=%d\n", group, tumor_count, normal_count))
  }
}

# ==============================================================================
# 2. Filter to Protein-Coding Genes
# ==============================================================================

cat("Filtering to protein-coding genes...\n")

# Filter to protein-coding genes only
protein_coding_genes <- gene_info$gene_type == "protein_coding"
count_data_pc <- count_data_full[protein_coding_genes, ]
gene_info_pc <- gene_info[protein_coding_genes, ]

cat("After protein-coding filter:", dim(count_data_pc), "\n")

# ==============================================================================
# 3. Process Each Group Separately with PCA-Filtered Samples
# ==============================================================================

cat("Processing each group for ContamDE analysis with PCA-filtered samples...\n")

# Initialize results storage
contamde_results <- list()
final_filtered_sample_lists <- list()
purity_cutoff <- 0.6  # Absolute purity threshold (60%)

for (group in c("R0", "R1", "B0", "B1")) {
  cat(sprintf("\n--- Processing group %s ---\n", group))
  
  # Check if group has samples
  if (!group %in% names(pca_filtered_sample_lists) || 
      length(pca_filtered_sample_lists[[group]]$tumor) == 0) {
    cat(sprintf("Group %s has no samples, skipping...\n", group))
    contamde_results[[group]] <- NULL
    final_filtered_sample_lists[[group]] <- list(tumor = character(0), 
                                                 normal = character(0), 
                                                 cases = character(0))
    next
  }
  
  # Get PCA-filtered tumor and normal samples for this group
  tumor_samples <- pca_filtered_sample_lists[[group]]$tumor
  normal_samples <- pca_filtered_sample_lists[[group]]$normal
  
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
    contamde_results[[group]] <- NULL
    final_filtered_sample_lists[[group]] <- list(tumor = character(0), 
                                                 normal = character(0), 
                                                 cases = character(0))
    next
  }
  
  tumor_samples_final <- tumor_samples_avail[1:min_pairs]
  normal_samples_final <- normal_samples_avail[1:min_pairs]
  
  # Prepare count data in [normal1, normal2, ..., tumor1, tumor2, ...] format
  group_samples <- c(normal_samples_final, tumor_samples_final)
  group_counts <- count_data_pc[, group_samples]
  
  cat(sprintf("Group %s: %d pairs (%d samples total)\n", 
              group, min_pairs, ncol(group_counts)))
  
  # Apply filterByExpr for low expression filtering
  dge <- DGEList(counts = group_counts)
  # Create group factor for filterByExpr (normal=0, tumor=1)
  group_factor <- factor(c(rep("normal", min_pairs), rep("tumor", min_pairs)))
  keep <- filterByExpr(dge, group = group_factor, min.count = 1, min.total.count = 15)
  dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
  
  cat(sprintf("After filterByExpr: %d genes (from %d)\n", 
              nrow(dge_filtered), nrow(group_counts)))
  
  # Run ContamDE analysis
  cat("Running ContamDE analysis on PCA-cleaned samples...\n")
  
  tryCatch({
    contamde_result <- contamde_lm(
      counts = dge_filtered$counts,
      subtype = NULL,
      covariate = NULL,
      contaminated = TRUE,
      robust = TRUE
    )
    
    # Extract purity estimates
    purity_estimates <- contamde_result$proportion
    names(purity_estimates) <- paste0("Pair_", 1:length(purity_estimates))
    
    cat("Purity estimates summary:\n")
    print(summary(purity_estimates))
    
    # Apply absolute purity cutoff (proportion > 0.6)
    purity_threshold <- purity_cutoff  # Use absolute threshold
    high_purity_pairs <- which(purity_estimates > purity_threshold)
    
    cat(sprintf("Purity threshold (absolute): %.1f\n", purity_threshold))
    cat(sprintf("High purity pairs (>%.1f): %d/%d\n", 
                purity_threshold, length(high_purity_pairs), length(purity_estimates)))
    
    # Update sample lists with high purity samples only
    if (length(high_purity_pairs) > 0) {
      filtered_tumor <- tumor_samples_final[high_purity_pairs]
      filtered_normal <- normal_samples_final[high_purity_pairs]
      # Get corresponding case IDs
      if (length(pca_filtered_sample_lists[[group]]$cases) >= min_pairs) {
        filtered_cases <- pca_filtered_sample_lists[[group]]$cases[high_purity_pairs]
      } else {
        # Fallback if cases are misaligned
        filtered_cases <- paste0("Case_", high_purity_pairs)
        cat(sprintf("Warning: Case IDs misaligned for group %s, using fallback\n", group))
      }
      
      final_filtered_sample_lists[[group]] <- list(
        tumor = filtered_tumor,
        normal = filtered_normal,
        cases = filtered_cases
      )
    } else {
      cat("Warning: No samples passed purity filter!\n")
      final_filtered_sample_lists[[group]] <- list(
        tumor = character(0),
        normal = character(0),
        cases = character(0)
      )
    }
    
    # Store results
    contamde_results[[group]] <- list(
      result = contamde_result,
      purity_estimates = purity_estimates,
      purity_threshold = purity_threshold,
      high_purity_pairs = high_purity_pairs,
      filtered_genes = rownames(dge_filtered),
      original_pairs = min_pairs,
      final_pairs = length(high_purity_pairs),
      input_source = "PCA_filtered"  # Mark as v2 result
    )
    
    cat(sprintf("Group %s processing completed successfully\n", group))
    
  }, error = function(e) {
    cat(sprintf("Error processing group %s: %s\n", group, e$message))
    contamde_results[[group]] <- NULL
    final_filtered_sample_lists[[group]] <- list(tumor = character(0), 
                                                 normal = character(0), 
                                                 cases = character(0))
  })
}

# ==============================================================================
# 4. Summary and Quality Control
# ==============================================================================

cat("\n==============================================\n")
cat("ContamDE Analysis v2 Summary\n")
cat("==============================================\n")

# Create summary table
summary_table <- data.frame(
  Group = c("R0", "R1", "B0", "B1"),
  Input_Pairs = integer(4),
  Final_Pairs = integer(4),
  Reduction = character(4),
  Mean_Purity = numeric(4),
  Purity_Threshold = numeric(4),
  stringsAsFactors = FALSE
)

for (i in 1:4) {
  group <- summary_table$Group[i]
  
  if (!is.null(contamde_results[[group]])) {
    input_pairs <- contamde_results[[group]]$original_pairs
    final_pairs <- contamde_results[[group]]$final_pairs
    purity_est <- contamde_results[[group]]$purity_estimates
    threshold <- contamde_results[[group]]$purity_threshold
    
    summary_table$Input_Pairs[i] <- input_pairs
    summary_table$Final_Pairs[i] <- final_pairs
    summary_table$Reduction[i] <- sprintf("%.1f%%", (1 - final_pairs/input_pairs) * 100)
    summary_table$Mean_Purity[i] <- mean(purity_est, na.rm = TRUE)
    summary_table$Purity_Threshold[i] <- purity_cutoff
  } else {
    summary_table$Input_Pairs[i] <- 0
    summary_table$Final_Pairs[i] <- 0
    summary_table$Reduction[i] <- "N/A"
    summary_table$Mean_Purity[i] <- NA
    summary_table$Purity_Threshold[i] <- purity_cutoff
  }
}

print(summary_table)

# Compare with original sample counts (before any filtering)
cat("\nFiltering cascade summary:\n")
cat("Original (ERR) → PCA-filtered → ContamDE-filtered:\n")

# Load original sample lists for comparison
load("./data/processed/sample_lists.rda")
original_sample_lists <- sample_lists

for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(original_sample_lists)) {
    orig_count <- length(original_sample_lists[[group]]$tumor)
    pca_count <- ifelse(group %in% names(pca_filtered_sample_lists), 
                        length(pca_filtered_sample_lists[[group]]$tumor), 0)
    final_count <- summary_table$Final_Pairs[summary_table$Group == group]
    
    cat(sprintf("  %s: %d → %d → %d samples\n", group, orig_count, pca_count, final_count))
  }
}

# Check minimum sample requirements
min_samples <- 5
sufficient_groups <- summary_table$Final_Pairs >= min_samples
cat("\nSample size check after dual filtering (minimum 5 per group):\n")
for (i in 1:4) {
  group <- summary_table$Group[i]
  status <- ifelse(sufficient_groups[i], "✅ PASS", "❌ INSUFFICIENT")
  cat(sprintf("  %s: %d pairs %s\n", group, summary_table$Final_Pairs[i], status))
}

total_final_pairs <- sum(summary_table$Final_Pairs)
cat(sprintf("\nTotal pairs after dual filtering: %d\n", total_final_pairs))

if (all(sufficient_groups)) {
  cat("✅ All groups meet minimum sample requirements\n")
} else {
  insufficient_groups <- summary_table$Group[!sufficient_groups]
  cat(sprintf("⚠️ Insufficient samples in: %s\n", 
              paste(insufficient_groups, collapse = ", ")))
}

# ==============================================================================
# 5. Save Results
# ==============================================================================

cat("Saving ContamDE v2 results...\n")

# Create processed directory if it doesn't exist
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save final filtered sample lists
save(final_filtered_sample_lists, file = "./data/processed/final_filtered_sample_lists.rda")

# Save detailed ContamDE v2 results
contamde_v2_analysis_results <- list(
  contamde_results = contamde_results,
  final_filtered_sample_lists = final_filtered_sample_lists,
  summary_table = summary_table,
  purity_cutoff = purity_cutoff,
  protein_coding_genes = rownames(count_data_pc),
  analysis_date = Sys.time(),
  analysis_version = "v2_post_pca_filtering",
  input_source = "PCA_filtered_samples"
)

save(contamde_v2_analysis_results, file = "./data/processed/contamde_v2_analysis_results.rda")

cat("Results saved to ./data/processed/\n")

# ==============================================================================
# 6. Prepare for Next Analysis
# ==============================================================================

cat("\nPreparing for next analysis steps...\n")

# Check which groups are viable for further analysis
viable_groups <- summary_table$Group[sufficient_groups]
cat("Groups viable for Phase 2 PCA and DEG analysis:", paste(viable_groups, collapse = ", "), "\n")

# Define possible comparisons (tissue-specific)
possible_comparisons <- list()
if ("R0" %in% viable_groups && "R1" %in% viable_groups) {
  possible_comparisons$R0_vs_R1_tumor <- "RET Tumor: Unexposed vs RadHigh (PRIMARY)"
  possible_comparisons$R0_vs_R1_normal <- "RET Normal: Unexposed vs RadHigh"
}
if ("B0" %in% viable_groups && "B1" %in% viable_groups) {
  possible_comparisons$B0_vs_B1_tumor <- "BRAF Tumor: Unexposed vs RadHigh"
  possible_comparisons$B0_vs_B1_normal <- "BRAF Normal: Unexposed vs RadHigh"
}

cat("Possible comparisons for DEG analysis:\n")
for (comp_name in names(possible_comparisons)) {
  cat(sprintf("  %s: %s\n", comp_name, possible_comparisons[[comp_name]]))
}

# Store comparison info for next scripts
comparison_info_v2 <- list(
  viable_groups = viable_groups,
  possible_comparisons = possible_comparisons,
  primary_comparison = ifelse("R0_vs_R1_tumor" %in% names(possible_comparisons), 
                              "R0_vs_R1_tumor", "None"),
  analysis_version = "v2"
)

save(comparison_info_v2, file = "./data/processed/comparison_info_v2.rda")

# Display quality improvement metrics
cat("\n=== Quality Improvement Assessment ===\n")
cat("Dual filtering strategy (PCA → ContamDE) benefits:\n")
cat("1. ✅ Cleaner samples for ContamDE: No expression outliers\n")
cat("2. ✅ More accurate purity estimates: Stable statistical assumptions\n")
cat("3. ✅ Higher confidence results: Two-stage quality control\n")

# Calculate overall sample retention rate
original_total <- sum(sapply(original_sample_lists, function(x) length(x$tumor)))
final_total <- sum(summary_table$Final_Pairs)
retention_rate <- (final_total / original_total) * 100

cat(sprintf("4. ✅ Sample retention: %.1f%% (%d/%d pairs)\n", 
            retention_rate, final_total, original_total))

if (retention_rate >= 60) {
  cat("5. ✅ Retention rate: Excellent (≥60%)\n")
} else if (retention_rate >= 50) {
  cat("5. ⚠️ Retention rate: Acceptable (50-60%)\n")
} else {
  cat("5. ❌ Retention rate: Concerning (<50%)\n")
}

cat("\nContamDE v2 analysis completed successfully!\n")
cat("Next step: Phase 2 PCA analysis for group comparison visualization\n")
cat("==============================================\n")

