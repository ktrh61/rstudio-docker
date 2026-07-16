# ==============================================================================
# REBC-THYR Feature Selection Script - TPM Ratio Based v2 (Phase 2)
# 11_feature_selection_tpm_based_v2.R
# ==============================================================================
# Purpose: FoldChange control pattern testing and candidate count verification
# Strategy: Statistical performance absolute priority

# Required libraries
library(dplyr)

cat("Starting TPM ratio-based feature selection v2 (Phase 2: FoldChange testing)...\n")

# ==============================================================================
# 1. Load Phase 1 Results
# ==============================================================================

cat("Loading Phase 1 results...\n")

if (!file.exists("./data/processed/feature_selection_phase1_results.rda")) {
  stop("Phase 1 results not found. Please run 11_feature_selection_tpm_based_v1.R first.")
}
load("./data/processed/feature_selection_phase1_results.rda")

if (!file.exists("./data/processed/analysis_tpm_matrix.rda")) {
  stop("Analysis TPM matrix not found.")
}
load("./data/processed/analysis_tpm_matrix.rda")

cat("Phase 1 data loaded successfully\n")
cat(sprintf("  UP genes available: %d\n", length(phase1_results$up_genes_avail)))
cat(sprintf("  DOWN genes available: %d\n", length(phase1_results$down_genes_avail)))
cat(sprintf("  Total possible pairs: %s\n", 
            format(phase1_results$total_possible_pairs, big.mark=",")))

# Extract data from Phase 1
up_genes_avail <- phase1_results$up_genes_avail
down_genes_avail <- phase1_results$down_genes_avail
analysis_samples <- phase1_results$analysis_samples
outcome <- phase1_results$outcome
r0_indices <- phase1_results$r0_indices
r1_indices <- phase1_results$r1_indices

# ==============================================================================
# 2. Remove Zero Expression Genes (No Pseudocount Strategy)
# ==============================================================================

cat("Removing genes with zero expression (no pseudocount strategy)...\n")

# Filter genes with no zero values
up_genes_nonzero <- up_genes_avail[rowSums(analysis_tpm[up_genes_avail, ] == 0) == 0]
down_genes_nonzero <- down_genes_avail[rowSums(analysis_tpm[down_genes_avail, ] == 0) == 0]

cat(sprintf("Zero expression filtering:\n"))
cat(sprintf("  UP genes: %d → %d (removed %d)\n", 
            length(up_genes_avail), length(up_genes_nonzero), 
            length(up_genes_avail) - length(up_genes_nonzero)))
cat(sprintf("  DOWN genes: %d → %d (removed %d)\n", 
            length(down_genes_avail), length(down_genes_nonzero), 
            length(down_genes_avail) - length(down_genes_nonzero)))
cat(sprintf("  Pairs after zero filtering: %s\n", 
            format(length(up_genes_nonzero) * length(down_genes_nonzero), big.mark=",")))

# Update gene lists
up_genes_final <- up_genes_nonzero
down_genes_final <- down_genes_nonzero

# ==============================================================================
# 3. Define FoldChange Control Patterns
# ==============================================================================

cat("Defining FoldChange control patterns for testing...\n")

# Define various FoldChange control strategies
foldchange_patterns <- list(
  # All samples minimum threshold
  "all_above_1.2x" = list(
    method = "all_above",
    threshold = log2(1.2),
    description = "All samples >= 1.2-fold change"
  ),
  "all_above_1.3x" = list(
    method = "all_above", 
    threshold = log2(1.3),
    description = "All samples >= 1.3-fold change"
  ),
  "all_above_1.5x" = list(
    method = "all_above",
    threshold = log2(1.5), 
    description = "All samples >= 1.5-fold change"
  ),
  
  # Quantile-based thresholds
  "q75_1.5x" = list(
    method = "quantile",
    quantile = 0.75,
    threshold = log2(1.5),
    description = "75% of samples >= 1.5-fold change"
  ),
  "q75_2x" = list(
    method = "quantile",
    quantile = 0.75,
    threshold = log2(2),
    description = "75% of samples >= 2-fold change"
  ),
  "q75_2.5x" = list(
    method = "quantile",
    quantile = 0.75,
    threshold = log2(2.5),
    description = "75% of samples >= 2.5-fold change"
  ),
  
  # Median-based thresholds
  "median_1.5x" = list(
    method = "median",
    threshold = log2(1.5),
    description = "Median >= 1.5-fold change"
  ),
  "median_2x" = list(
    method = "median",
    threshold = log2(2),
    description = "Median >= 2-fold change"
  ),
  "median_2.5x" = list(
    method = "median",
    threshold = log2(2.5),
    description = "Median >= 2.5-fold change"
  )
)

cat(sprintf("Defined %d FoldChange control patterns\n", length(foldchange_patterns)))

# ==============================================================================
# 4. Test Function for Single Pair
# ==============================================================================

# Function to test a single pair against all criteria
test_single_pair <- function(up_gene, down_gene, patterns_to_test) {
  
  # Extract TPM values (no pseudocount)
  up_tpm <- as.numeric(analysis_tpm[up_gene, ])
  down_tpm <- as.numeric(analysis_tpm[down_gene, ])
  
  # Calculate ratios (no pseudocount - already filtered zero genes)
  ratios <- log2(up_tpm / down_tpm)
  
  # Split by groups
  r0_ratios <- ratios[r0_indices]
  r1_ratios <- ratios[r1_indices]
  
  # Basic criteria
  r0_median <- median(r0_ratios)
  r1_median <- median(r1_ratios)
  
  # Direction check
  direction_opposite <- r0_median != r1_median
  separation <- abs(r1_median - r0_median)
  
  # ERR=0 specificity check (90% minimum)
  if (r1_median > r0_median) {
    threshold_est <- max(r0_ratios) + 0.1
    err0_spec <- sum(r0_ratios <= threshold_est) / length(r0_ratios)
  } else {
    threshold_est <- min(r0_ratios) - 0.1
    err0_spec <- sum(r0_ratios >= threshold_est) / length(r0_ratios)
  }
  
  # Test each FoldChange pattern
  pattern_results <- list()
  
  for (pattern_name in names(patterns_to_test)) {
    pattern <- patterns_to_test[[pattern_name]]
    
    # Calculate absolute ratios for FoldChange test
    abs_ratios <- abs(ratios)
    
    # Apply pattern-specific test
    if (pattern$method == "all_above") {
      passes_foldchange <- all(abs_ratios >= pattern$threshold)
      
    } else if (pattern$method == "quantile") {
      q_value <- quantile(abs_ratios, pattern$quantile)
      passes_foldchange <- q_value >= pattern$threshold
      
    } else if (pattern$method == "median") {
      median_value <- median(abs_ratios)
      passes_foldchange <- median_value >= pattern$threshold
    }
    
    # Overall pass/fail for this pattern
    passes_all <- (err0_spec >= 0.90) && direction_opposite && passes_foldchange
    
    pattern_results[[pattern_name]] <- list(
      passes_foldchange = passes_foldchange,
      passes_all = passes_all
    )
  }
  
  return(list(
    up_gene = up_gene,
    down_gene = down_gene,
    err0_specificity = err0_spec,
    separation = separation,
    direction_opposite = direction_opposite,
    pattern_results = pattern_results
  ))
}

# ==============================================================================
# 5. Sample-Based Pattern Testing
# ==============================================================================

cat("Testing FoldChange patterns with sample pairs...\n")

# Test with a sample of pairs for efficiency
n_sample_pairs <- min(1000, length(up_genes_final) * length(down_genes_final))
set.seed(123)

# Generate sample pair indices
total_pairs <- length(up_genes_final) * length(down_genes_final)
sample_indices <- sample(total_pairs, n_sample_pairs)

cat(sprintf("Testing %s pairs out of %s total (%.1f%% sample)\n", 
            format(n_sample_pairs, big.mark=","),
            format(total_pairs, big.mark=","),
            n_sample_pairs/total_pairs*100))

# Test sample pairs
sample_results <- list()
pattern_pass_counts <- data.frame(
  pattern = names(foldchange_patterns),
  description = sapply(foldchange_patterns, function(x) x$description),
  foldchange_pass = 0,
  all_criteria_pass = 0,
  stringsAsFactors = FALSE
)

cat("Testing sample pairs... (this may take a moment)\n")
progress_interval <- max(1, floor(n_sample_pairs / 10))

for (i in seq_along(sample_indices)) {
  
  # Convert linear index to up/down indices
  pair_index <- sample_indices[i]
  up_idx <- ((pair_index - 1) %% length(up_genes_final)) + 1
  down_idx <- ceiling(pair_index / length(up_genes_final))
  
  up_gene <- up_genes_final[up_idx]
  down_gene <- down_genes_final[down_idx]
  
  # Test this pair
  tryCatch({
    result <- test_single_pair(up_gene, down_gene, foldchange_patterns)
    sample_results[[i]] <- result
    
    # Update counts
    for (j in 1:nrow(pattern_pass_counts)) {
      pattern_name <- pattern_pass_counts$pattern[j]
      if (pattern_name %in% names(result$pattern_results)) {
        if (result$pattern_results[[pattern_name]]$passes_foldchange) {
          pattern_pass_counts$foldchange_pass[j] <- pattern_pass_counts$foldchange_pass[j] + 1
        }
        if (result$pattern_results[[pattern_name]]$passes_all) {
          pattern_pass_counts$all_criteria_pass[j] <- pattern_pass_counts$all_criteria_pass[j] + 1
        }
      }
    }
  }, error = function(e) {
    cat(sprintf("Error with pair %d: %s\n", i, e$message))
  })
  
  # Progress reporting
  if (i %% progress_interval == 0) {
    cat(sprintf("  Processed %d/%d pairs (%.1f%%)\n", i, n_sample_pairs, i/n_sample_pairs*100))
  }
}

# Calculate percentages
pattern_pass_counts$foldchange_pass_pct <- pattern_pass_counts$foldchange_pass / n_sample_pairs * 100
pattern_pass_counts$all_criteria_pass_pct <- pattern_pass_counts$all_criteria_pass / n_sample_pairs * 100

# ==============================================================================
# 6. Results Analysis and Reporting
# ==============================================================================

cat("\n==============================================\n")
cat("FoldChange Pattern Testing Results\n")
cat("==============================================\n")

cat(sprintf("Sample testing completed: %d pairs tested\n\n", n_sample_pairs))

cat("Pattern Performance Summary:\n")
print(pattern_pass_counts[, c("pattern", "description", "foldchange_pass", 
                              "foldchange_pass_pct", "all_criteria_pass", 
                              "all_criteria_pass_pct")], row.names = FALSE)

# Identify viable patterns
viable_patterns <- pattern_pass_counts[pattern_pass_counts$all_criteria_pass > 0, ]
cat(sprintf("\nViable patterns (with candidates): %d/%d\n", 
            nrow(viable_patterns), nrow(pattern_pass_counts)))

if (nrow(viable_patterns) > 0) {
  cat("\nViable patterns ranked by candidate count:\n")
  viable_ranked <- viable_patterns[order(viable_patterns$all_criteria_pass, decreasing = TRUE), ]
  print(viable_ranked[, c("pattern", "all_criteria_pass", "all_criteria_pass_pct")], row.names = FALSE)
  
  # Estimate full dataset candidates
  cat("\nEstimated candidates for full dataset:\n")
  for (i in 1:nrow(viable_ranked)) {
    pattern_name <- viable_ranked$pattern[i]
    sample_rate <- viable_ranked$all_criteria_pass_pct[i] / 100
    estimated_full <- round(total_pairs * sample_rate)
    cat(sprintf("  %s: ~%s candidates (%.3f%%)\n", 
                pattern_name, format(estimated_full, big.mark=","), sample_rate*100))
  }
} else {
  cat("\n⚠️  No patterns produced viable candidates!\n")
  cat("Consider relaxing criteria:\n")
  cat("  - Lower ERR=0 specificity threshold (currently 90%)\n")
  cat("  - Reduce FoldChange requirements\n")
  cat("  - Alternative filtering strategies\n")
}

# ==============================================================================
# 7. Detailed Analysis of Best Patterns
# ==============================================================================

if (nrow(viable_patterns) > 0) {
  cat("\n==============================================\n")
  cat("Detailed Analysis of Top Patterns\n")
  cat("==============================================\n")
  
  # Analyze top 3 patterns
  top_patterns <- head(viable_ranked, 3)
  
  for (i in 1:nrow(top_patterns)) {
    pattern_name <- top_patterns$pattern[i]
    cat(sprintf("\n--- Pattern %d: %s ---\n", i, pattern_name))
    cat(sprintf("Description: %s\n", top_patterns$description[i]))
    cat(sprintf("Sample candidates: %d/%d (%.2f%%)\n", 
                top_patterns$all_criteria_pass[i], n_sample_pairs, 
                top_patterns$all_criteria_pass_pct[i]))
    
    # Find sample results for this pattern
    pattern_sample_results <- list()
    for (result in sample_results) {
      if (!is.null(result) && pattern_name %in% names(result$pattern_results)) {
        if (result$pattern_results[[pattern_name]]$passes_all) {
          pattern_sample_results <- append(pattern_sample_results, list(result), after = length(pattern_sample_results))
        }
      }
    }
    
    if (length(pattern_sample_results) > 0) {
      # Statistics for passing pairs
      separations <- sapply(pattern_sample_results, function(x) x$separation)
      err0_specs <- sapply(pattern_sample_results, function(x) x$err0_specificity)
      
      cat(sprintf("Passing pairs statistics:\n"))
      cat(sprintf("  Separation: mean=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
                  mean(separations), median(separations), min(separations), max(separations)))
      cat(sprintf("  ERR=0 spec: mean=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
                  mean(err0_specs), median(err0_specs), min(err0_specs), max(err0_specs)))
      
      # Show top 3 examples
      if (length(pattern_sample_results) >= 3) {
        top_by_separation <- pattern_sample_results[order(separations, decreasing = TRUE)[1:3]]
        cat(sprintf("Top 3 examples by separation:\n"))
        for (j in 1:3) {
          result <- top_by_separation[[j]]
          cat(sprintf("  %d. %s / %s: sep=%.3f, err0=%.1f%%\n", 
                      j, result$up_gene, result$down_gene, 
                      result$separation, result$err0_specificity*100))
        }
      }
    }
  }
}

# ==============================================================================
# 8. Save Phase 2 Results
# ==============================================================================

cat("\nSaving Phase 2 results...\n")

phase2_results <- list(
  # Input data
  up_genes_final = up_genes_final,
  down_genes_final = down_genes_final,
  total_pairs_after_zero_filter = length(up_genes_final) * length(down_genes_final),
  
  # Testing parameters
  sample_size = n_sample_pairs,
  sampling_rate = n_sample_pairs / total_pairs,
  
  # Pattern definitions
  foldchange_patterns = foldchange_patterns,
  
  # Results
  pattern_pass_counts = pattern_pass_counts,
  viable_patterns = viable_patterns,
  sample_results = sample_results,
  
  # Analysis data (from Phase 1)
  analysis_samples = analysis_samples,
  outcome = outcome,
  r0_indices = r0_indices,
  r1_indices = r1_indices,
  
  analysis_date = Sys.time(),
  analysis_version = "phase2_foldchange_testing"
)

save(phase2_results, file = "./data/processed/feature_selection_phase2_results.rda")

cat("Phase 2 results saved successfully\n")

# ==============================================================================
# 9. Recommendations for Phase 3
# ==============================================================================

cat("\n==============================================\n")
cat("Recommendations for Phase 3\n")
cat("==============================================\n")

if (nrow(viable_patterns) > 0) {
  best_pattern <- viable_ranked$pattern[1]
  cat(sprintf("✅ PROCEED TO PHASE 3\n"))
  cat(sprintf("Recommended pattern: %s\n", best_pattern))
  cat(sprintf("Expected candidates: ~%s\n", 
              format(round(total_pairs * viable_ranked$all_criteria_pass_pct[1] / 100), big.mark=",")))
  cat(sprintf("\nPhase 3 should implement:\n"))
  cat(sprintf("1. Full enumeration with %s pattern\n", best_pattern))
  cat(sprintf("2. Final pair selection (1-2 pairs)\n"))
  cat(sprintf("3. LOOCV and Bootstrap validation\n"))
  cat(sprintf("4. Performance optimization\n"))
} else {
  cat("❌ CANNOT PROCEED TO PHASE 3\n")
  cat("No viable patterns found. Consider:\n")
  cat("1. Relaxing ERR=0 specificity from 90% to 85%\n")
  cat("2. Reducing FoldChange thresholds\n")
  cat("3. Alternative filtering strategies\n")
  cat("4. Reviewing DEG quality and criteria\n")
}

cat("\nPhase 2 completed!\n")
cat("==============================================\n")

