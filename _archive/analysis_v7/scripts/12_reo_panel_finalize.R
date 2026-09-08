#!/usr/bin/env Rscript
# =============================================================================
# 12_reo_panel_finalize.R
# REO Panel Construction - Finalization with Boundary Zone
# 
# Purpose: 
#   - Phase 5: Remove redundancy (gene overlap, high correlation)
#   - Phase 6: Evaluate separation and set boundary zone (data-driven)
#   - Phase 7: Define classification rules (negative/undetermined/positive)
#
# Input: 
#   - reo_candidate_pairs.rds (from 11_reo_pair_selection.R)
# Output:
#   - reo_final_panel.rds (final panel with boundary zone)
#   - reo_final_panel.csv
#   - reo_sample_classification.csv
#
# Version: v7.1
# Date: 2025-12-09
#
# Key Changes from v1:
#   - Ranking by median_diff (not separation score)
#   - Data-driven boundary zone: [max(R0)+1, min(R1)-1] = undetermined
#   - Three-class output: negative / undetermined / positive
# =============================================================================

source("analysis_v7/setup.R")

cat("\n=============================================================================\n")
cat("12_reo_panel_finalize.R - REO Panel Finalization (v7.1)\n")
cat("Date: 2025-12-09\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
PARAMS <- list(
  target_panel_size = 10,
  correlation_threshold = 0.75,    # Spearman correlation threshold
  dead_zone = log2(1.2)            # Must match 11_reo_pair_selection.R
)

cat("=== Parameters ===\n")
cat(sprintf("Target panel size: %d\n", PARAMS$target_panel_size))
cat(sprintf("Correlation threshold: %.2f (Spearman)\n", PARAMS$correlation_threshold))
cat(sprintf("Dead zone: |r| < %.4f -> non-reversal\n", PARAMS$dead_zone))
cat(sprintf("Selection: by median_diff (greedy, avoiding gene overlap)\n"))
cat("\n")

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
cat("=== Loading Data ===\n")

input <- readRDS(file.path(paths$processed, "reo_candidate_pairs.rds"))

candidates <- input$candidates
log2_tpm <- input$log2_tpm
r0_samples <- input$r0_samples
r1_samples <- input$r1_samples
n_r0 <- length(r0_samples)
n_r1 <- length(r1_samples)

cat(sprintf("Candidate pairs: %d\n", nrow(candidates)))
cat(sprintf("R0 samples: %d\n", n_r0))
cat(sprintf("R1 samples: %d\n", n_r1))
cat(sprintf("Input version: %s\n", input$version))
cat("\n")

# -----------------------------------------------------------------------------
# Phase 5a: Sort by Priority (median_diff, already sorted in 11_)
# -----------------------------------------------------------------------------
cat("=== Phase 5a: Priority Sorting ===\n")

# Ensure sorted by median_diff descending
candidates <- candidates[order(-candidates$median_diff), ]

cat(sprintf("Top candidate: %s x %s (median_diff=%.3f, rev=%d/%d)\n",
            candidates$up_name[1], candidates$down_name[1],
            candidates$median_diff[1], candidates$reversal_count[1], n_r1))
cat("\n")

# -----------------------------------------------------------------------------
# Phase 5b: Greedy Selection (Gene Non-overlap + Low Correlation)
# -----------------------------------------------------------------------------
cat("=== Phase 5b: Greedy Selection ===\n")
cat("Criteria: No gene overlap + Spearman correlation < threshold\n\n")

used_genes <- character(0)
selected_pairs <- data.frame()
selected_r_matrix <- NULL

for (k in 1:nrow(candidates)) {
  up_g <- candidates$up_gene[k]
  down_g <- candidates$down_gene[k]
  up_name <- candidates$up_name[k]
  down_name <- candidates$down_name[k]
  
  # Check gene overlap
  if ((up_g %in% used_genes) || (down_g %in% used_genes)) next
  
  # Compute r values for this pair (all samples)
  all_samples <- c(r0_samples, r1_samples)
  new_r <- log2_tpm[up_g, all_samples] - log2_tpm[down_g, all_samples]
  
  # Skip if contains -Inf or Inf
  if (any(is.infinite(new_r))) {
    cat(sprintf("  Skipped: %s x %s (contains Inf values)\n", up_name, down_name))
    next
  }
  
  # Check correlation with already selected pairs (Spearman)
  if (!is.null(selected_r_matrix)) {
    cors <- apply(selected_r_matrix, 2, function(x) {
      abs(cor(x, new_r, method = "spearman", use = "complete.obs"))
    })
    cors[is.na(cors)] <- 0
    if (any(cors >= PARAMS$correlation_threshold)) {
      max_cor_idx <- which.max(cors)
      cat(sprintf("  Skipped: %s x %s (cor=%.2f with P%d)\n", 
                  up_name, down_name, max(cors), max_cor_idx))
      next
    }
  }
  
  # Add to selected
  selected_pairs <- rbind(selected_pairs, candidates[k, ])
  used_genes <- c(used_genes, up_g, down_g)
  selected_r_matrix <- cbind(selected_r_matrix, new_r)
  
  cat(sprintf("  P%d: %s x %s (median_diff=%.3f, rev=%d/%d)\n",
              nrow(selected_pairs), up_name, down_name,
              candidates$median_diff[k], candidates$reversal_count[k], n_r1))
  
  # Stop when target reached
  if (nrow(selected_pairs) >= PARAMS$target_panel_size) break
}

cat(sprintf("\n=== Selected %d pairs ===\n\n", nrow(selected_pairs)))

# Check if enough pairs selected
if (nrow(selected_pairs) < PARAMS$target_panel_size) {
  cat(sprintf("WARNING: Only %d pairs selected (target: %d)\n", 
              nrow(selected_pairs), PARAMS$target_panel_size))
}

# Assign pair IDs
selected_pairs$pair_id <- paste0("P", 1:nrow(selected_pairs))
rownames(selected_pairs) <- selected_pairs$pair_id

# -----------------------------------------------------------------------------
# Display Final Panel
# -----------------------------------------------------------------------------
cat("=== Final Panel ===\n")
print(selected_pairs[, c("pair_id", "up_name", "down_name", "median_diff",
                         "reversal_count", "reversal_rate", "q10")])

# Correlation matrix (Spearman)
if (nrow(selected_pairs) > 1) {
  colnames(selected_r_matrix) <- selected_pairs$pair_id
  cor_matrix <- cor(selected_r_matrix, method = "spearman")
  
  cat("\n=== Pair Correlation Matrix (Spearman) ===\n")
  print(round(cor_matrix, 2))
  
  max_cor <- max(abs(cor_matrix[upper.tri(cor_matrix)]))
  cat(sprintf("\nMax off-diagonal correlation: %.2f\n", max_cor))
}

# -----------------------------------------------------------------------------
# Phase 6: Compute Reversal Counts for All Samples
# -----------------------------------------------------------------------------
cat("\n=== Phase 6: Sample-level Reversal Evaluation ===\n")

compute_sample_reversals <- function(selected_pairs, log2_tpm, 
                                     r0_samples, r1_samples, dead_zone) {
  
  n_pairs <- nrow(selected_pairs)
  all_samples <- c(r0_samples, r1_samples)
  n_all <- length(all_samples)
  
  # Matrix to store r values
  r_matrix <- matrix(NA, nrow = n_all, ncol = n_pairs)
  rownames(r_matrix) <- all_samples
  colnames(r_matrix) <- selected_pairs$pair_id
  
  # Matrix to store reversal status
  reversal_matrix <- matrix(0, nrow = n_all, ncol = n_pairs)
  rownames(reversal_matrix) <- all_samples
  colnames(reversal_matrix) <- selected_pairs$pair_id
  
  # R0 majority signs (pre-computed in selection)
  r0_majority_signs <- selected_pairs$r0_majority_sign
  names(r0_majority_signs) <- selected_pairs$pair_id
  
  for (k in 1:n_pairs) {
    up_g <- selected_pairs$up_gene[k]
    down_g <- selected_pairs$down_gene[k]
    pair_id <- selected_pairs$pair_id[k]
    
    # Compute r values
    r_all <- log2_tpm[up_g, all_samples] - log2_tpm[down_g, all_samples]
    r_matrix[, k] <- r_all
    
    # Reversal: outside dead zone AND opposite sign
    majority_sign <- r0_majority_signs[pair_id]
    outside_dz <- abs(r_all) >= dead_zone
    opposite <- sign(r_all) != majority_sign
    
    reversal_matrix[, k] <- as.numeric(outside_dz & opposite)
  }
  
  return(list(
    r_matrix = r_matrix,
    reversal_matrix = reversal_matrix,
    r0_majority_signs = r0_majority_signs
  ))
}

eval_result <- compute_sample_reversals(selected_pairs, log2_tpm, 
                                        r0_samples, r1_samples,
                                        PARAMS$dead_zone)
r_matrix <- eval_result$r_matrix
reversal_matrix <- eval_result$reversal_matrix

# Reversal counts
r0_reversals <- rowSums(reversal_matrix[r0_samples, , drop = FALSE])
r1_reversals <- rowSums(reversal_matrix[r1_samples, , drop = FALSE])

cat("\n--- R0 Sample Reversal Counts ---\n")
cat("(Expected: low, as these are sporadic cases)\n")
print(sort(r0_reversals, decreasing = TRUE))

cat("\n--- R1 Sample Reversal Counts ---\n")
cat("(Expected: variable, as R1 is a mixed group)\n")
print(sort(r1_reversals, decreasing = TRUE))

# Summary
cat("\n--- Summary ---\n")
cat(sprintf("R0: range=%d-%d, mean=%.1f, median=%.0f\n", 
            min(r0_reversals), max(r0_reversals), 
            mean(r0_reversals), median(r0_reversals)))
cat(sprintf("R1: range=%d-%d, mean=%.1f, median=%.0f\n", 
            min(r1_reversals), max(r1_reversals),
            mean(r1_reversals), median(r1_reversals)))

# -----------------------------------------------------------------------------
# Phase 7: Data-Driven Boundary Zone
# -----------------------------------------------------------------------------
cat("\n=== Phase 7: Boundary Zone Definition ===\n")

n_pairs <- nrow(selected_pairs)
max_r0 <- max(r0_reversals)
min_r1 <- min(r1_reversals)

cat(sprintf("R0 max reversal: %d\n", max_r0))
cat(sprintf("R1 min reversal: %d\n", min_r1))

# Define boundary zone
# negative: <= max_r0
# positive: >= min_r1
# undetermined: (max_r0, min_r1) if gap exists

if (max_r0 < min_r1) {
  # Clear separation - gap exists
  boundary <- list(
    negative_max = max_r0,
    positive_min = min_r1,
    undetermined_range = if (min_r1 - max_r0 > 1) (max_r0 + 1):(min_r1 - 1) else NULL,
    has_gap = TRUE,
    gap_size = min_r1 - max_r0 - 1
  )
  cat(sprintf("\n*** Clear separation detected ***\n"))
  cat(sprintf("Negative: reversal <= %d\n", boundary$negative_max))
  cat(sprintf("Positive: reversal >= %d\n", boundary$positive_min))
  if (!is.null(boundary$undetermined_range)) {
    cat(sprintf("Undetermined: reversal in {%s}\n", 
                paste(boundary$undetermined_range, collapse = ", ")))
  } else {
    cat("Undetermined: none (adjacent thresholds)\n")
  }
} else {
  # Overlap exists
  boundary <- list(
    negative_max = max_r0,
    positive_min = min_r1,
    undetermined_range = min_r1:max_r0,
    has_gap = FALSE,
    overlap_size = max_r0 - min_r1 + 1
  )
  cat(sprintf("\n*** WARNING: Overlap detected ***\n"))
  cat(sprintf("R0 max (%d) >= R1 min (%d)\n", max_r0, min_r1))
  cat(sprintf("Overlap range: {%s}\n", 
              paste(boundary$undetermined_range, collapse = ", ")))
  cat("Consider: values in overlap range = undetermined\n")
}

# Classification function
classify_sample <- function(reversal_count, boundary) {
  if (reversal_count <= boundary$negative_max && 
      (boundary$has_gap || reversal_count < boundary$positive_min)) {
    return("negative")
  } else if (reversal_count >= boundary$positive_min && 
             (boundary$has_gap || reversal_count > boundary$negative_max)) {
    return("positive")
  } else {
    return("undetermined")
  }
}

# Apply classification
r0_class <- sapply(r0_reversals, classify_sample, boundary = boundary)
r1_class <- sapply(r1_reversals, classify_sample, boundary = boundary)

cat("\n--- Classification Results ---\n")
cat("\nR0 (sporadic, expected: negative):\n")
print(table(r0_class))

cat("\nR1 (mixed, expected: mostly positive with some negative/undetermined):\n")
print(table(r1_class))

# -----------------------------------------------------------------------------
# Sensitivity Analysis: Threshold Sweep
# -----------------------------------------------------------------------------
cat("\n=== Threshold Sensitivity Analysis ===\n")
cat("(For reference - showing performance at different thresholds)\n\n")

cat("Threshold | R0 neg | R0 und | R0 pos | R1 neg | R1 und | R1 pos\n")
cat("--------- | ------ | ------ | ------ | ------ | ------ | ------\n")

threshold_results <- data.frame()

for (T in 0:n_pairs) {
  # Simple threshold: < T = negative, >= T = positive
  r0_pos_simple <- sum(r0_reversals >= T)
  r0_neg_simple <- sum(r0_reversals < T)
  r1_pos_simple <- sum(r1_reversals >= T)
  r1_neg_simple <- sum(r1_reversals < T)
  
  cat(sprintf("T = %2d    | %2d     | -      | %2d     | %2d     | -      | %2d\n",
              T, r0_neg_simple, r0_pos_simple, r1_neg_simple, r1_pos_simple))
  
  threshold_results <- rbind(threshold_results, data.frame(
    T = T,
    r0_negative = r0_neg_simple,
    r0_positive = r0_pos_simple,
    r1_negative = r1_neg_simple,
    r1_positive = r1_pos_simple,
    r0_specificity = r0_neg_simple / n_r0 * 100,
    r1_positive_rate = r1_pos_simple / n_r1 * 100
  ))
}

# Highlight data-driven boundary
cat(sprintf("\n--- Data-Driven Boundary ---\n"))
cat(sprintf("Negative threshold: <= %d (all R0 classified as negative)\n", max_r0))
cat(sprintf("Positive threshold: >= %d (based on R1 minimum)\n", min_r1))

# -----------------------------------------------------------------------------
# Notable Samples
# -----------------------------------------------------------------------------
cat("\n=== Notable Samples ===\n")

# R0 with highest reversal
cat("\nR0 samples with highest reversal:\n")
r0_sorted <- sort(r0_reversals, decreasing = TRUE)
for (i in 1:min(3, length(r0_sorted))) {
  s <- names(r0_sorted)[i]
  cat(sprintf("  %s: %d/%d reversals -> %s\n", 
              s, r0_reversals[s], n_pairs, r0_class[s]))
}

# R1 with lowest reversal
cat("\nR1 samples with lowest reversal:\n")
r1_sorted <- sort(r1_reversals, decreasing = FALSE)
for (i in 1:min(3, length(r1_sorted))) {
  s <- names(r1_sorted)[i]
  cat(sprintf("  %s: %d/%d reversals -> %s (likely sporadic in R1)\n", 
              s, r1_reversals[s], n_pairs, r1_class[s]))
}

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------
cat("\n=== Saving Results ===\n")

output <- list(
  # Version info
  version = "v7.1",
  date = "2025-12-09",
  
  # Parameters
  params = PARAMS,
  
  # Boundary zone (data-driven)
  boundary = boundary,
  
  # Panel definition
  selected_pairs = selected_pairs,
  correlation_matrix = if (exists("cor_matrix")) cor_matrix else NULL,
  
  # R0 majority signs
  r0_majority_signs = eval_result$r0_majority_signs,
  
  # Raw data matrices
  r_matrix = r_matrix,
  reversal_matrix = reversal_matrix,
  
  # Sample-level results
  r0_reversals = r0_reversals,
  r1_reversals = r1_reversals,
  r0_class = r0_class,
  r1_class = r1_class,
  
  # Threshold analysis
  threshold_results = threshold_results,
  
  # Expression data reference
  log2_tpm = log2_tpm,
  r0_samples = r0_samples,
  r1_samples = r1_samples,
  
  # Metadata
  timestamp = Sys.time()
)

saveRDS(output, file.path(paths$processed, "reo_final_panel.rds"))
cat(sprintf("Saved: %s\n", file.path(paths$processed, "reo_final_panel.rds")))

# Panel summary CSV
panel_csv <- selected_pairs[, c("pair_id", "up_name", "down_name", 
                                "up_gene", "down_gene",
                                "median_diff", "reversal_count", 
                                "reversal_rate", "q10")]
write.csv(panel_csv, file.path(paths$output, "reo_final_panel.csv"), row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_final_panel.csv")))

# Sample classification CSV
sample_results <- data.frame(
  sample_id = c(r0_samples, r1_samples),
  group = c(rep("R0", n_r0), rep("R1", n_r1)),
  reversal_count = c(r0_reversals, r1_reversals),
  classification = c(r0_class, r1_class)
)
write.csv(sample_results, file.path(paths$output, "reo_sample_classification.csv"), 
          row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_sample_classification.csv")))

# -----------------------------------------------------------------------------
# Final Summary
# -----------------------------------------------------------------------------
cat("\n=============================================================================\n")
cat("Phase 5-7 Complete\n")
cat("=============================================================================\n")

cat(sprintf("\nFinal panel: %d pairs\n", nrow(selected_pairs)))
cat(sprintf("Ranking: |median_R0 - median_R1| (top = %.3f)\n", 
            selected_pairs$median_diff[1]))

cat(sprintf("\nBoundary Zone:\n"))
cat(sprintf("  Negative: reversal <= %d\n", boundary$negative_max))
cat(sprintf("  Positive: reversal >= %d\n", boundary$positive_min))
if (boundary$has_gap && !is.null(boundary$undetermined_range)) {
  cat(sprintf("  Undetermined: {%s}\n", 
              paste(boundary$undetermined_range, collapse = ", ")))
}

cat(sprintf("\nClassification Summary:\n"))
cat(sprintf("  R0: %d negative, %d undetermined, %d positive\n",
            sum(r0_class == "negative"), 
            sum(r0_class == "undetermined"),
            sum(r0_class == "positive")))
cat(sprintf("  R1: %d negative, %d undetermined, %d positive\n",
            sum(r1_class == "negative"),
            sum(r1_class == "undetermined"),
            sum(r1_class == "positive")))

cat("\n=== Done ===\n")