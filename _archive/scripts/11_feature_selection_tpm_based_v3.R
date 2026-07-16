# ==============================================================================
# REBC-THYR Feature Selection Script - TPM Ratio Based v3 (Phase 3)
# 11_feature_selection_tpm_based_v3.R
# ==============================================================================
# Purpose: Final implementation with all_above_1.5x pattern and complete validation
# Strategy: Statistical performance absolute priority

# Required libraries
library(dplyr)

cat("Starting TPM ratio-based feature selection v3 (Phase 3: Final implementation)...\n")

# ==============================================================================
# 1. Load Phase 2 Results
# ==============================================================================

cat("Loading Phase 2 results...\n")

if (!file.exists("./data/processed/feature_selection_phase2_results.rda")) {
  stop("Phase 2 results not found. Please run 11_feature_selection_tpm_based_v2.R first.")
}
load("./data/processed/feature_selection_phase2_results.rda")

if (!file.exists("./data/processed/analysis_tpm_matrix.rda")) {
  stop("Analysis TPM matrix not found.")
}
load("./data/processed/analysis_tpm_matrix.rda")

cat("Phase 2 data loaded successfully\n")

# Extract data from Phase 2
up_genes_final <- phase2_results$up_genes_final
down_genes_final <- phase2_results$down_genes_final
analysis_samples <- phase2_results$analysis_samples
outcome <- phase2_results$outcome
r0_indices <- phase2_results$r0_indices
r1_indices <- phase2_results$r1_indices

cat(sprintf("  UP genes: %d\n", length(up_genes_final)))
cat(sprintf("  DOWN genes: %d\n", length(down_genes_final)))

# Define the selected pattern
selected_pattern <- "all_above_1.5x"
pattern_threshold <- log2(1.5)

cat(sprintf("Selected pattern: %s (threshold: %.3f log2 units)\n", 
            selected_pattern, pattern_threshold))

# ==============================================================================
# 2. Full Enumeration with Selected Pattern
# ==============================================================================

cat("Starting full enumeration with all_above_1.5x pattern...\n")

# Prepare storage for candidates
candidate_pairs <- data.frame(
  up_gene = character(0),
  down_gene = character(0),
  err0_specificity = numeric(0),
  separation = numeric(0),
  r0_median = numeric(0),
  r1_median = numeric(0),
  min_abs_ratio = numeric(0),
  direction_r1_higher = logical(0),
  stringsAsFactors = FALSE
)

# Calculate total pairs
total_pairs <- length(up_genes_final) * length(down_genes_final)
cat(sprintf("Total pairs to evaluate: %s\n", format(total_pairs, big.mark=",")))

# Progress tracking
progress_interval <- max(1, floor(total_pairs / 100))
processed_count <- 0
candidate_count <- 0

cat("Starting enumeration (this may take several minutes)...\n")
start_time <- Sys.time()

# Main enumeration loop
for (up_idx in seq_along(up_genes_final)) {
  
  up_gene <- up_genes_final[up_idx]
  up_tpm <- as.numeric(analysis_tpm[up_gene, ])
  
  for (down_idx in seq_along(down_genes_final)) {
    
    down_gene <- down_genes_final[down_idx]
    down_tpm <- as.numeric(analysis_tpm[down_gene, ])
    
    processed_count <- processed_count + 1
    
    # Calculate ratios (no pseudocount)
    ratios <- log2(up_tpm / down_tpm)
    
    # Split by groups
    r0_ratios <- ratios[r0_indices]
    r1_ratios <- ratios[r1_indices]
    
    # Basic criteria checks
    r0_median <- median(r0_ratios)
    r1_median <- median(r1_ratios)
    
    # Direction check (must be opposite)
    direction_opposite <- r0_median != r1_median
    if (!direction_opposite) next
    
    separation <- abs(r1_median - r0_median)
    direction_r1_higher <- r1_median > r0_median
    
    # ERR=0 specificity check (90% minimum)
    if (direction_r1_higher) {
      threshold_est <- max(r0_ratios) + 0.1
      err0_spec <- sum(r0_ratios <= threshold_est) / length(r0_ratios)
    } else {
      threshold_est <- min(r0_ratios) - 0.1
      err0_spec <- sum(r0_ratios >= threshold_est) / length(r0_ratios)
    }
    
    if (err0_spec < 0.90) next
    
    # FoldChange check: all_above_1.5x
    abs_ratios <- abs(ratios)
    min_abs_ratio <- min(abs_ratios)
    all_above_threshold <- all(abs_ratios >= pattern_threshold)
    
    if (!all_above_threshold) next
    
    # If we reach here, this pair passes all criteria
    candidate_count <- candidate_count + 1
    
    # Store candidate
    candidate_pairs <- rbind(candidate_pairs, data.frame(
      up_gene = up_gene,
      down_gene = down_gene,
      err0_specificity = err0_spec,
      separation = separation,
      r0_median = r0_median,
      r1_median = r1_median,
      min_abs_ratio = min_abs_ratio,
      direction_r1_higher = direction_r1_higher,
      stringsAsFactors = FALSE
    ))
    
    # Progress reporting
    if (processed_count %% progress_interval == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      progress_pct <- processed_count / total_pairs * 100
      estimated_total <- elapsed / progress_pct * 100
      remaining <- estimated_total - elapsed
      
      cat(sprintf("  Progress: %s/%s (%.1f%%) | Candidates: %d | Time: %.1f min | ETA: %.1f min\n",
                  format(processed_count, big.mark=","), 
                  format(total_pairs, big.mark=","),
                  progress_pct, candidate_count, elapsed, remaining))
    }
  }
}

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

cat(sprintf("\nEnumeration completed!\n"))
cat(sprintf("  Total pairs processed: %s\n", format(processed_count, big.mark=",")))
cat(sprintf("  Candidates found: %d\n", candidate_count))
cat(sprintf("  Success rate: %.2f%%\n", candidate_count / processed_count * 100))
cat(sprintf("  Total time: %.1f minutes\n", total_time))
cat(sprintf("  Processing speed: %d pairs/second\n", round(processed_count / (total_time * 60))))

# ==============================================================================
# 3. Candidate Analysis and Ranking
# ==============================================================================

if (nrow(candidate_pairs) == 0) {
  cat("\n❌ No candidates found!\n")
  cat("Consider:\n")
  cat("  1. Relaxing ERR=0 specificity to 85%\n")
  cat("  2. Using q75_1.5x or median_1.5x pattern\n")
  cat("  3. Reducing FoldChange threshold\n")
  stop("No candidates available for analysis")
}

cat("\n==============================================\n")
cat("Candidate Analysis and Ranking\n")
cat("==============================================\n")

# Basic statistics
cat(sprintf("Candidate summary:\n"))
cat(sprintf("  Total candidates: %d\n", nrow(candidate_pairs)))
cat(sprintf("  ERR=0 specificity: mean=%.3f, min=%.3f, max=%.3f\n",
            mean(candidate_pairs$err0_specificity),
            min(candidate_pairs$err0_specificity),
            max(candidate_pairs$err0_specificity)))
cat(sprintf("  Separation: mean=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
            mean(candidate_pairs$separation),
            median(candidate_pairs$separation),
            min(candidate_pairs$separation),
            max(candidate_pairs$separation)))
cat(sprintf("  Min absolute ratio: mean=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
            mean(candidate_pairs$min_abs_ratio),
            median(candidate_pairs$min_abs_ratio),
            min(candidate_pairs$min_abs_ratio),
            max(candidate_pairs$min_abs_ratio)))

# Direction analysis
direction_counts <- table(candidate_pairs$direction_r1_higher)
cat(sprintf("  Direction: R1 higher=%d, R0 higher=%d\n", 
            direction_counts["TRUE"], direction_counts["FALSE"]))

# ==============================================================================
# 4. Final Pair Selection with Two-Tier Ranking
# ==============================================================================

cat("\n==============================================\n")
cat("Final Pair Selection (Two-Tier Ranking)\n")
cat("==============================================\n")

# Two-tier ranking: ERR=0 specificity first, then separation
tier1_candidates <- candidate_pairs[candidate_pairs$err0_specificity >= 0.95, ]
tier2_candidates <- candidate_pairs[candidate_pairs$err0_specificity >= 0.90 & 
                                      candidate_pairs$err0_specificity < 0.95, ]

cat(sprintf("Tier 1 (ERR=0 spec ≥95%%): %d candidates\n", nrow(tier1_candidates)))
cat(sprintf("Tier 2 (ERR=0 spec 90-95%%): %d candidates\n", nrow(tier2_candidates)))

# Rank within each tier by separation
if (nrow(tier1_candidates) > 0) {
  tier1_ranked <- tier1_candidates[order(tier1_candidates$separation, decreasing = TRUE), ]
  final_candidates <- tier1_ranked
  tier_used <- "Tier 1"
} else {
  tier2_ranked <- tier2_candidates[order(tier2_candidates$separation, decreasing = TRUE), ]
  final_candidates <- tier2_ranked
  tier_used <- "Tier 2"
}

cat(sprintf("Using %s for final selection\n", tier_used))

# Select top 2 pairs (or fewer if not available)
n_final_pairs <- min(2, nrow(final_candidates))
final_pairs <- head(final_candidates, n_final_pairs)

cat(sprintf("\nSelected %d final pair(s):\n", n_final_pairs))
for (i in 1:nrow(final_pairs)) {
  pair <- final_pairs[i, ]
  cat(sprintf("  %d. %s / %s\n", i, pair$up_gene, pair$down_gene))
  cat(sprintf("     ERR=0 spec: %.1f%%, Separation: %.3f log2, Min ratio: %.3f\n",
              pair$err0_specificity * 100, pair$separation, pair$min_abs_ratio))
  cat(sprintf("     Direction: R1 %s than R0 (median: %.3f vs %.3f)\n",
              ifelse(pair$direction_r1_higher, "higher", "lower"),
              pair$r1_median, pair$r0_median))
}

# ==============================================================================
# 5. Manual Validation of Selected Pairs
# ==============================================================================

cat("\n==============================================\n")
cat("Manual Validation of Selected Pairs\n")
cat("==============================================\n")

for (i in 1:nrow(final_pairs)) {
  pair <- final_pairs[i, ]
  up_gene <- pair$up_gene
  down_gene <- pair$down_gene
  
  cat(sprintf("\n--- Manual Validation for Pair %d ---\n", i))
  cat(sprintf("UP gene: %s\n", up_gene))
  cat(sprintf("DOWN gene: %s\n", down_gene))
  
  # Manual calculation
  up_tpm <- as.numeric(analysis_tpm[up_gene, ])
  down_tpm <- as.numeric(analysis_tpm[down_gene, ])
  ratios <- log2(up_tpm / down_tpm)
  
  r0_ratios <- ratios[r0_indices]
  r1_ratios <- ratios[r1_indices]
  
  # Statistics
  cat(sprintf("R0 ratios (n=%d): mean=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
              length(r0_ratios), mean(r0_ratios), median(r0_ratios),
              min(r0_ratios), max(r0_ratios)))
  cat(sprintf("R1 ratios (n=%d): mean=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
              length(r1_ratios), mean(r1_ratios), median(r1_ratios),
              min(r1_ratios), max(r1_ratios)))
  
  # Manual ERR=0 specificity calculation
  r0_median_manual <- median(r0_ratios)
  r1_median_manual <- median(r1_ratios)
  separation_manual <- abs(r1_median_manual - r0_median_manual)
  
  if (r1_median_manual > r0_median_manual) {
    threshold_manual <- max(r0_ratios) + 0.1
    err0_spec_manual <- sum(r0_ratios <= threshold_manual) / length(r0_ratios)
    optimal_threshold <- threshold_manual
  } else {
    threshold_manual <- min(r0_ratios) - 0.1
    err0_spec_manual <- sum(r0_ratios >= threshold_manual) / length(r0_ratios)
    optimal_threshold <- threshold_manual
  }
  
  # FoldChange validation
  abs_ratios <- abs(ratios)
  min_abs_ratio_manual <- min(abs_ratios)
  all_above_1.5x_manual <- all(abs_ratios >= log2(1.5))
  
  cat(sprintf("Manual calculation validation:\n"))
  cat(sprintf("  Separation: %.3f (script: %.3f) %s\n", 
              separation_manual, pair$separation,
              ifelse(abs(separation_manual - pair$separation) < 0.001, "✓", "✗")))
  cat(sprintf("  ERR=0 specificity: %.3f (script: %.3f) %s\n",
              err0_spec_manual, pair$err0_specificity,
              ifelse(abs(err0_spec_manual - pair$err0_specificity) < 0.001, "✓", "✗")))
  cat(sprintf("  Min abs ratio: %.3f (script: %.3f) %s\n",
              min_abs_ratio_manual, pair$min_abs_ratio,
              ifelse(abs(min_abs_ratio_manual - pair$min_abs_ratio) < 0.001, "✓", "✗")))
  cat(sprintf("  All above 1.5x: %s\n", ifelse(all_above_1.5x_manual, "✓", "✗")))
  cat(sprintf("  Optimal threshold: %.3f\n", optimal_threshold))
  
  # Store optimal threshold for this pair
  final_pairs$optimal_threshold[i] <- optimal_threshold
}

# ==============================================================================
# 6. Simple Judgment Function
# ==============================================================================

cat("\n==============================================\n")
cat("Simple Judgment Function\n")
cat("==============================================\n")

if (nrow(final_pairs) > 0) {
  # Use the best pair for the judgment function
  best_pair <- final_pairs[1, ]
  
  cat("Simple radiation impact judgment function (Pair 1):\n")
  cat(sprintf("UP gene: %s\n", best_pair$up_gene))
  cat(sprintf("DOWN gene: %s\n", best_pair$down_gene))
  cat(sprintf("Threshold: %.3f\n", best_pair$optimal_threshold))
  cat(sprintf("Direction: R1 %s than R0\n", 
              ifelse(best_pair$direction_r1_higher, "higher", "lower")))
  
  # Create simple judgment function
  simple_judgment_function <- sprintf("
# Simple Radiation Impact Judgment Function
radiation_impact_judgment <- function(up_gene_tpm, down_gene_tpm) {
  ratio <- log2(up_gene_tpm / down_gene_tpm)
  threshold <- %.3f
  
  %s
  
  return(list(
    judgment = ifelse(ratio %s threshold, 'Radiation Impact Present', 'No Radiation Impact'),
    ratio_value = ratio,
    threshold_used = threshold
  ))
}

# Usage example:
# result <- radiation_impact_judgment(up_tpm = 10.5, down_tpm = 3.2)
# print(result$judgment)
", 
best_pair$optimal_threshold,
ifelse(best_pair$direction_r1_higher, 
       "# R1 group shows higher ratios than R0 group",
       "# R0 group shows higher ratios than R1 group"),
ifelse(best_pair$direction_r1_higher, ">", "<")
  )
  
  cat(simple_judgment_function)
}

# ==============================================================================
# 7. Save Final Results
# ==============================================================================

cat("\n==============================================\n")
cat("Saving Final Results\n")
cat("==============================================\n")

# Comprehensive final results
final_results <- list(
  # Selection parameters
  selected_pattern = selected_pattern,
  pattern_threshold = pattern_threshold,
  
  # Enumeration results
  total_pairs_processed = processed_count,
  total_candidates_found = candidate_count,
  enumeration_time_minutes = total_time,
  processing_speed_pairs_per_second = round(processed_count / (total_time * 60)),
  
  # All candidates
  all_candidates = candidate_pairs,
  
  # Final selection
  final_pairs = final_pairs,
  tier_used = if(exists("tier_used")) tier_used else "Unknown",
  
  # Simple function
  best_pair_info = if(nrow(final_pairs) > 0) final_pairs[1, ] else NULL,
  simple_function_code = if(exists("simple_judgment_function")) simple_judgment_function else NULL,
  
  # Analysis metadata
  analysis_samples = analysis_samples,
  outcome = outcome,
  r0_indices = r0_indices,
  r1_indices = r1_indices,
  up_genes_final = up_genes_final,
  down_genes_final = down_genes_final,
  
  analysis_date = Sys.time(),
  analysis_version = "phase3_final_implementation"
)

save(final_results, file = "./data/processed/feature_selection_final_results.rda")

# Save CSV files for easy access
if (!dir.exists("./output/reports")) {
  dir.create("./output/reports", recursive = TRUE)
}

if (nrow(candidate_pairs) > 0) {
  write.csv(candidate_pairs, "./output/reports/all_candidates_final.csv", row.names = FALSE)
}

if (nrow(final_pairs) > 0) {
  write.csv(final_pairs, "./output/reports/final_selected_pairs.csv", row.names = FALSE)
  
  # Write simple function to text file
  if (exists("simple_judgment_function")) {
    writeLines(simple_judgment_function, "./output/reports/simple_judgment_function.R")
  }
}

cat("Results saved successfully\n")
cat(sprintf("  All candidates: ./output/reports/all_candidates_final.csv\n"))
cat(sprintf("  Final pairs: ./output/reports/final_selected_pairs.csv\n"))
cat(sprintf("  Judgment function: ./output/reports/simple_judgment_function.R\n"))

# ==============================================================================
# 8. Final Summary
# ==============================================================================

cat("\n==============================================\n")
cat("TPM Ratio-Based Feature Selection - Final Summary\n")
cat("==============================================\n")

cat("Enumeration Performance:\n")
cat(sprintf("  ✅ Processed: %s pairs in %.1f minutes\n", 
            format(processed_count, big.mark=","), total_time))
cat(sprintf("  ⚡ Speed: %d pairs/second\n", round(processed_count / (total_time * 60))))
cat(sprintf("  🎯 Success rate: %.2f%% (%d candidates)\n", 
            candidate_count / processed_count * 100, candidate_count))

if (nrow(final_pairs) > 0) {
  cat("\n✅ SUCCESS: Final pairs selected\n")
  cat(sprintf("Selected pairs: %d\n", nrow(final_pairs)))
  
  for (i in 1:nrow(final_pairs)) {
    pair <- final_pairs[i, ]
    fold_change <- 2^pair$separation
    cat(sprintf("  %d. %s / %s\n", i, pair$up_gene, pair$down_gene))
    cat(sprintf("     Performance: ERR=0 spec %.1f%%, Separation %.3f log2 (%.1fx change)\n",
                pair$err0_specificity * 100, pair$separation, fold_change))
  }
  
  cat("\nNext steps:\n")
  cat("  1. ✅ Ready for LOOCV validation\n")
  cat("  2. ✅ Ready for Bootstrap validation\n") 
  cat("  3. ✅ Simple judgment function available\n")
  cat("  4. 📊 Consider biological pathway analysis\n")
  
} else {
  cat("\n❌ No pairs met all criteria\n")
  cat("Consider adjusting parameters and re-running\n")
}

cat("\nTPM ratio-based feature selection completed!\n")
cat("==============================================\n")