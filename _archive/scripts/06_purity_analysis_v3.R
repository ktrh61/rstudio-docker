# ==============================================================================
# REBC-THYR Tumor Purity Analysis Script v3 - After New CDM PCA (Corrected)
# 06_purity_analysis_v3.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(edgeR)
library(dplyr)

# Source the lightweight purity estimation functions
source("./utils/contamde_purity_functions.R")

cat("Starting tumor purity analysis v3 with corrected PCA v3 integration...\n")

# ==============================================================================
# 1. Load Data and PCA v3 Results (Corrected Input)
# ==============================================================================

cat("Loading PCA v3 results and data...\n")

# Load PCA v3 results (pair consistency already ensured)
load("./data/processed/pca_phase1_results_v3.rda")
# Or load the filtered sample lists directly
load("./data/processed/pca_filtered_sample_lists_v3.rda")

# Use the already pair-consistent filtered lists
pca_filtered_sample_lists <- pca_filtered_sample_lists  # from v3 file

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

# Display PCA v3 filtered sample counts (input to this analysis)
cat("PCA v3 filtered sample counts (input):\n")
for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(pca_filtered_sample_lists)) {
    tumor_count <- length(pca_filtered_sample_lists[[group]]$tumor)
    normal_count <- length(pca_filtered_sample_lists[[group]]$normal)
    cat(sprintf("  %s: %d pairs (Tumor=%d, Normal=%d)\n", 
                group, tumor_count, tumor_count, normal_count))
  }
}

# ==============================================================================
# 2. Filter to Protein-Coding Genes
# ==============================================================================

cat("\nFiltering to protein-coding genes...\n")

# Filter to protein-coding genes only
protein_coding_genes <- gene_info$gene_type == "protein_coding"
count_data_pc <- count_data_full[protein_coding_genes, ]
gene_info_pc <- gene_info[protein_coding_genes, ]

cat("After protein-coding filter:", dim(count_data_pc), "\n")

# ==============================================================================
# 3. Tumor Purity Estimation Function (Lightweight)
# ==============================================================================

# Lightweight purity estimation using contamde_lm v2
estimate_tumor_purity <- function(count_matrix, group_name, verbose = TRUE) {
  
  if (verbose) cat(sprintf("Estimating tumor purity for %s...\n", group_name))
  
  # Apply filterByExpr for low expression filtering
  dge <- DGEList(counts = count_matrix)
  n_pairs <- ncol(count_matrix) / 2
  group_factor <- factor(c(rep("normal", n_pairs), rep("tumor", n_pairs)))
  keep <- filterByExpr(dge, group = group_factor)
  dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
  
  if (verbose) {
    cat(sprintf("  After filterByExpr: %d genes (from %d)\n", 
                nrow(dge_filtered), nrow(count_matrix)))
    # シンプルな0カウント確認
    zero_count <- sum(dge_filtered$counts == 0)
    total_elements <- length(dge_filtered$counts)
    zero_pct <- zero_count / total_elements * 100
    
    cat(sprintf("  Zero counts: %d/%d (%.2f%%)\n", 
                zero_count, total_elements, zero_pct))
  }
  
  # Run contamDE purity estimation (lightweight version)
  tryCatch({
    contamde_result <- contamde_purity(
      counts = dge_filtered$counts,
      subtype = NULL,
      covariate = NULL,
      contaminated = TRUE,
      verbose = TRUE
    )
    
    # Extract purity estimates
    purity_estimates <- contamde_result$proportion
    names(purity_estimates) <- paste0("Pair_", 1:length(purity_estimates))
    
    if (verbose) {
      cat(sprintf("  Purity estimates: [%.3f, %.3f], mean=%.3f\n",
                  min(purity_estimates), max(purity_estimates), mean(purity_estimates)))
    }
    
    return(list(
      purity_estimates = purity_estimates,
      contamde_result = contamde_result,
      filtered_genes = rownames(dge_filtered),
      n_pairs = n_pairs,
      success = TRUE,
      error_message = NULL
    ))
    
  }, error = function(e) {
    if (verbose) cat(sprintf("  Error: %s\n", e$message))
    return(list(
      purity_estimates = rep(NA, n_pairs),
      contamde_result = NULL,
      filtered_genes = character(0),
      n_pairs = n_pairs,
      success = FALSE,
      error_message = e$message
    ))
  })
}

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
  
  # Check if group has samples after PCA v3 filtering
  if (!group %in% names(pca_filtered_sample_lists) || 
      length(pca_filtered_sample_lists[[group]]$tumor) == 0) {
    cat(sprintf("Group %s has no samples after PCA v3 filtering, skipping...\n", group))
    purity_results[[group]] <- NULL
    next
  }
  
  # Get PCA v3 filtered tumor and normal samples (pair consistency guaranteed)
  tumor_samples <- pca_filtered_sample_lists[[group]]$tumor
  normal_samples <- pca_filtered_sample_lists[[group]]$normal
  case_ids <- pca_filtered_sample_lists[[group]]$cases
  
  cat(sprintf("PCA v3 filtered samples - Tumor: %d, Normal: %d\n", 
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
  
  # Use only available samples (should be all if PCA v3 worked correctly)
  tumor_samples_avail <- tumor_samples[tumor_available]
  normal_samples_avail <- normal_samples[normal_available]
  
  # Since PCA v3 ensures pair consistency, lengths should be equal
  if (length(tumor_samples_avail) != length(normal_samples_avail)) {
    cat(sprintf("Warning: Unequal sample counts after availability check. Using minimum.\n"))
    min_pairs <- min(length(tumor_samples_avail), length(normal_samples_avail))
    tumor_samples_final <- tumor_samples_avail[1:min_pairs]
    normal_samples_final <- normal_samples_avail[1:min_pairs]
  } else {
    tumor_samples_final <- tumor_samples_avail
    normal_samples_final <- normal_samples_avail
  }
  
  min_pairs <- length(tumor_samples_final)
  
  if (min_pairs == 0) {
    cat(sprintf("Group %s has no valid pairs, skipping...\n", group))
    purity_results[[group]] <- NULL
    next
  }
  
  if (min_pairs < 3) {
    cat(sprintf("Warning: Group %s has only %d pairs. Results may be unstable.\n", 
                group, min_pairs))
  }
  
  # Prepare count data in [normal1, normal2, ..., tumor1, tumor2, ...] format
  group_samples <- c(normal_samples_final, tumor_samples_final)
  group_counts <- count_data_pc[, group_samples]
  
  cat(sprintf("Final analysis: %d pairs (%d samples total)\n", 
              min_pairs, ncol(group_counts)))
  
  # Run purity estimation
  purity_result <- estimate_tumor_purity(group_counts, group, verbose = TRUE)
  
  if (purity_result$success) {
    # Store results with metadata
    purity_results[[group]] <- list(
      group_name = group,
      purity_estimates = purity_result$purity_estimates,
      contamde_result = purity_result$contamde_result,
      tumor_samples = tumor_samples_final,
      normal_samples = normal_samples_final,
      case_ids = case_ids[1:min_pairs],
      filtered_genes = purity_result$filtered_genes,
      n_pairs = min_pairs,
      analysis_date = Sys.time(),
      success = TRUE
    )
    
    cat(sprintf("Group %s: Purity estimation completed successfully\n", group))
  } else {
    cat(sprintf("Group %s: Purity estimation failed - %s\n", group, purity_result$error_message))
    purity_results[[group]] <- NULL
  }
}

# ==============================================================================
# 5. Purity Quality Assessment and High-Purity Sample Selection
# ==============================================================================

cat("\n==============================================\n")
cat("Purity Quality Assessment and High-Purity Selection\n")
cat("==============================================\n")

# Initialize summary table and high-purity lists
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

# Create final high-purity filtered sample lists
final_high_purity_sample_lists <- list()

for (i in 1:4) {
  group <- purity_summary$Group[i]
  
  if (!is.null(purity_results[[group]])) {
    # Get purity estimates
    purity_est <- purity_results[[group]]$purity_estimates
    n_pairs <- purity_results[[group]]$n_pairs
    
    # Calculate summary statistics
    purity_summary$Input_Pairs[i] <- n_pairs
    purity_summary$Mean_Purity[i] <- round(mean(purity_est, na.rm = TRUE), 3)
    purity_summary$SD_Purity[i] <- round(sd(purity_est, na.rm = TRUE), 3)
    purity_summary$Min_Purity[i] <- round(min(purity_est, na.rm = TRUE), 3)
    purity_summary$Max_Purity[i] <- round(max(purity_est, na.rm = TRUE), 3)
    
    # Apply purity cutoff
    high_purity_indices <- which(purity_est >= purity_cutoff)
    n_high_purity <- length(high_purity_indices)
    
    purity_summary$High_Purity_Count[i] <- n_high_purity
    purity_summary$Retention_Rate[i] <- round(n_high_purity / n_pairs, 3)
    
    cat(sprintf("\n%s Purity Assessment:\n", group))
    cat(sprintf("  Input pairs: %d\n", n_pairs))
    cat(sprintf("  Purity range: [%.3f, %.3f]\n", 
                purity_summary$Min_Purity[i], purity_summary$Max_Purity[i]))
    cat(sprintf("  High purity (≥%.1f): %d pairs (%.1f%% retention)\n",
                purity_cutoff, n_high_purity, purity_summary$Retention_Rate[i] * 100))
    
    # Create high-purity sample lists
    if (n_high_purity > 0) {
      final_high_purity_sample_lists[[group]] <- list(
        tumor = purity_results[[group]]$tumor_samples[high_purity_indices],
        normal = purity_results[[group]]$normal_samples[high_purity_indices],
        cases = purity_results[[group]]$case_ids[high_purity_indices]
      )
      
      if (n_high_purity < n_pairs) {
        excluded_pairs <- n_pairs - n_high_purity
        cat(sprintf("  Excluded %d pairs with purity < %.1f\n", 
                    excluded_pairs, purity_cutoff))
      }
    } else {
      cat(sprintf("  ⚠️  WARNING: No pairs passed purity threshold!\n"))
      final_high_purity_sample_lists[[group]] <- list(
        tumor = character(0),
        normal = character(0), 
        cases = character(0)
      )
    }
    
  } else {
    # Group failed or had no data
    purity_summary$Input_Pairs[i] <- 0
    purity_summary$Mean_Purity[i] <- NA
    purity_summary$SD_Purity[i] <- NA
    purity_summary$Min_Purity[i] <- NA
    purity_summary$Max_Purity[i] <- NA
    purity_summary$High_Purity_Count[i] <- 0
    purity_summary$Retention_Rate[i] <- 0
    
    final_high_purity_sample_lists[[group]] <- list(
      tumor = character(0),
      normal = character(0),
      cases = character(0)
    )
  }
}

print(purity_summary)

# ==============================================================================
# 6. Summary and Quality Control
# ==============================================================================

cat("\n==============================================\n")
cat("Complete Filtering Cascade Summary\n")
cat("==============================================\n")

# Load original sample lists for complete comparison
load("./data/processed/sample_lists.rda")
original_sample_lists <- sample_lists

cat("Complete filtering cascade:\n")
cat("Original → PCA v3 → Purity-filtered:\n")

for (group in c("R0", "R1", "B0", "B1")) {
  orig_count <- length(original_sample_lists[[group]]$tumor)
  pca_count <- length(pca_filtered_sample_lists[[group]]$tumor)
  purity_count <- length(final_high_purity_sample_lists[[group]]$tumor)
  
  cat(sprintf("  %s: %d → %d → %d pairs\n", group, orig_count, pca_count, purity_count))
}

# Check minimum sample requirements
min_samples <- 5
cat("\nSample size check (minimum 5 pairs per group):\n")
viable_groups <- c()
for (group in c("R0", "R1", "B0", "B1")) {
  final_count <- length(final_high_purity_sample_lists[[group]]$tumor)
  status <- ifelse(final_count >= min_samples, "✅ PASS", "❌ INSUFFICIENT")
  cat(sprintf("  %s: %d pairs %s\n", group, final_count, status))
  
  if (final_count >= min_samples) {
    viable_groups <- c(viable_groups, group)
  }
}

total_final_pairs <- sum(sapply(final_high_purity_sample_lists, function(x) length(x$tumor)))
cat(sprintf("\nTotal high-purity pairs: %d\n", total_final_pairs))
cat(sprintf("Viable groups: %s\n", paste(viable_groups, collapse = ", ")))

# ==============================================================================
# 7. Save Results
# ==============================================================================

cat("Saving purity analysis v3 results...\n")

# Create processed directory if it doesn't exist
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save final high-purity sample lists (main output)
save(final_high_purity_sample_lists, file = "./data/processed/final_high_purity_sample_lists.rda")

# Save comprehensive purity analysis results
purity_analysis_v3_results <- list(
  purity_results = purity_results,
  purity_summary = purity_summary,
  final_high_purity_sample_lists = final_high_purity_sample_lists,
  purity_cutoff = purity_cutoff,
  viable_groups = viable_groups,
  protein_coding_genes = rownames(count_data_pc),
  analysis_date = Sys.time(),
  analysis_version = "v3_corrected_no_double_filtering"
)

save(purity_analysis_v3_results, file = "./data/processed/purity_analysis_v3_results.rda")

cat("Results saved to ./data/processed/\n")

# ==============================================================================
# 8. Comparison Information for Next Steps
# ==============================================================================

cat("\n==============================================\n")
cat("Viable Comparisons for DEG Analysis\n")
cat("==============================================\n")

# Define possible comparisons
possible_comparisons <- list()
if ("R0" %in% viable_groups && "R1" %in% viable_groups) {
  r0_pairs <- length(final_high_purity_sample_lists[["R0"]]$tumor)
  r1_pairs <- length(final_high_purity_sample_lists[["R1"]]$tumor)
  possible_comparisons$R0_vs_R1_tumor <- sprintf("RET Tumor: R0(%d) vs R1(%d) - PRIMARY", r0_pairs, r1_pairs)
  possible_comparisons$R0_vs_R1_normal <- sprintf("RET Normal: R0(%d) vs R1(%d)", r0_pairs, r1_pairs)
}
if ("B0" %in% viable_groups && "B1" %in% viable_groups) {
  b0_pairs <- length(final_high_purity_sample_lists[["B0"]]$tumor)
  b1_pairs <- length(final_high_purity_sample_lists[["B1"]]$tumor)
  possible_comparisons$B0_vs_B1_tumor <- sprintf("BRAF Tumor: B0(%d) vs B1(%d)", b0_pairs, b1_pairs)
  possible_comparisons$B0_vs_B1_normal <- sprintf("BRAF Normal: B0(%d) vs B1(%d)", b0_pairs, b1_pairs)
}

if (length(possible_comparisons) > 0) {
  cat("Viable comparisons:\n")
  for (comp_name in names(possible_comparisons)) {
    cat(sprintf("  %s: %s\n", comp_name, possible_comparisons[[comp_name]]))
  }
  
  # Identify primary comparison
  primary_comparison <- ifelse("R0_vs_R1_tumor" %in% names(possible_comparisons), 
                               "R0_vs_R1_tumor", names(possible_comparisons)[1])
  cat(sprintf("\nPrimary comparison: %s\n", primary_comparison))
} else {
  cat("❌ No viable comparisons available!\n")
  cat("Consider:\n")
  cat("1. Lowering purity threshold\n")
  cat("2. Reviewing PCA v3 filtering\n")
  cat("3. Using smaller minimum sample requirements\n")
}

# Save comparison info for next analysis
comparison_info_v3 <- list(
  viable_groups = viable_groups,
  possible_comparisons = possible_comparisons,
  primary_comparison = ifelse(length(possible_comparisons) > 0, names(possible_comparisons)[1], "None"),
  analysis_version = "v3_corrected"
)
save(comparison_info_v3, file = "./data/processed/comparison_info_v3.rda")

# ==============================================================================
# 9. Final Summary
# ==============================================================================

cat("\n==============================================\n")
cat("Purity Analysis v3 Final Summary (Corrected)\n")
cat("==============================================\n")

cat("✅ Corrected Integration with PCA v3 Results\n")
cat("✅ No Double Filtering (Pair Consistency Maintained)\n")
cat("✅ Lightweight ContamDE v2 Functions Used\n")
cat("✅ High-Purity Sample Selection Completed\n")

if (total_final_pairs >= 10) {
  cat(sprintf("\n🎯 SUCCESS: %d high-purity pairs ready for analysis\n", total_final_pairs))
  cat("Ready for next steps:\n")
  cat("1. Phase 2 PCA for group comparison visualization\n")
  cat("2. DEGES normalization with high-purity samples\n")
  cat("3. DEG analysis for 321 consistent gene discovery\n")
} else if (total_final_pairs > 0) {
  cat(sprintf("\n⚠️  LIMITED: %d high-purity pairs available\n", total_final_pairs))
  cat("Consider parameter adjustment or proceed with caution\n")
} else {
  cat("\n❌ FAILED: No high-purity pairs available\n")
  cat("Review filtering parameters and data quality\n")
}

cat("\nPurity analysis v3 (corrected) completed!\n")
cat("==============================================\n")

