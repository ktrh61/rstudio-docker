#!/usr/bin/env Rscript
# =============================================================================
# 12_reo_panel_finalize.R
# REO Panel Construction - Phase 5-7: Finalization
# 
# Purpose: 
#   - Phase 5: Remove redundancy (gene overlap, high correlation)
#   - Phase 6: Evaluate separation between R0 and R1
#   - Phase 7: Set voting threshold T
#
# Input: 
#   - reo_candidate_pairs.rds (from 11_reo_pair_selection.R)
# Output:
#   - reo_final_panel.rds (final panel with threshold, size depends on candidates)
#
# Reference: REO_Panel_Protocol_v2.md
# Date: 2025-12-08
# v1.1: Spearman correlation for pair independence check, >= threshold
# v1.2: Threshold T determined data-driven: max(R0 reversal) + 1
# v1.3: Sorting by separation score (from 11_reo_pair_selection.R)
# v1.4: Filter to separation > 0 only; max_panel_size as upper limit
# =============================================================================

source("analysis_v7/setup.R")

cat("\n=============================================================================\n")
cat("12_reo_panel_finalize.R - REO Panel Finalization\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
# Rationale:
# 1. Max panel size = 10: Upper limit, actual size depends on available candidates
# 2. Separation > 0: Only use pairs with R0/R1 distributions fully separated
# 3. Correlation threshold = 0.75: Spearman (rank-based, appropriate for REO)
# 4. Selection: By separation score (from 11_reo_pair_selection.R)
# 5. Threshold T: Data-driven, set to max(R0 reversal) + 1 to ensure 100% specificity

PARAMS <- list(
  max_panel_size = 10,
  correlation_threshold = 0.75,    # Spearman correlation (REO-appropriate)
  dead_zone_threshold = log2(1.5)  # For voting: |r| < this -> non-reversal
)

cat("=== Parameters ===\n")
cat(sprintf("Max panel size: %d\n", PARAMS$max_panel_size))
cat(sprintf("Separation requirement: > 0 (R0/R1 fully separated)\n"))
cat(sprintf("Correlation threshold: %.2f (Spearman)\n", PARAMS$correlation_threshold))
cat(sprintf("Selection: by separation score (greedy, avoiding gene overlap)\n"))
cat(sprintf("Dead zone for voting: |r| < %.3f -> non-reversal\n", 
            PARAMS$dead_zone_threshold))
cat("\n")

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
cat("=== Loading Data ===\n")

input <- readRDS(file.path(paths$processed, "reo_candidate_pairs.rds"))

candidates_all <- input$candidates
log2_tpm <- input$log2_tpm
r0_samples <- input$r0_samples
r1_samples <- input$r1_samples
n_r0 <- length(r0_samples)
n_r1 <- length(r1_samples)

cat(sprintf("Candidate pairs (total): %d\n", nrow(candidates_all)))

# Filter to separation > 0 only (R0/R1 distributions must be separated)
candidates <- candidates_all[candidates_all$separation > 0, ]
cat(sprintf("Candidate pairs (separation > 0): %d\n", nrow(candidates)))
cat(sprintf("R0 samples: %d\n", n_r0))
cat(sprintf("R1 samples: %d\n", n_r1))

if (nrow(candidates) == 0) {
  stop("No candidates with separation > 0. Consider relaxing selection criteria in 11_reo_pair_selection.R")
}
cat("\n")

# -----------------------------------------------------------------------------
# Phase 5a: Sort by Priority (Separation Score)
# -----------------------------------------------------------------------------
cat("=== Phase 5a: Priority Sorting ===\n")

# Already sorted by separation in 11_reo_pair_selection.R
# Re-sort to ensure consistency
candidates <- candidates[order(-candidates$separation), ]

cat(sprintf("Top candidate: %s x %s (separation=%.2f, rev=%d)\n",
            candidates$up_name[1], candidates$down_name[1],
            candidates$separation[1], candidates$reversal[1]))
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
  
  # Skip if contains -Inf or Inf (log2(0) issue)
  if (any(is.infinite(new_r))) {
    cat(sprintf("  Skipped: %s x %s (contains Inf values)\n", up_name, down_name))
    next
  }
  
  # Check correlation with already selected pairs (Spearman for REO consistency)
  if (!is.null(selected_r_matrix)) {
    cors <- apply(selected_r_matrix, 2, function(x) {
      abs(cor(x, new_r, method = "spearman", use = "complete.obs"))
    })
    # Handle NA correlations (treat as acceptable)
    cors[is.na(cors)] <- 0
    if (any(cors >= PARAMS$correlation_threshold)) {
      next
    }
  }
  
  # Add to selected
  selected_pairs <- rbind(selected_pairs, candidates[k, ])
  used_genes <- c(used_genes, up_g, down_g)
  selected_r_matrix <- cbind(selected_r_matrix, new_r)
  
  cat(sprintf("  P%d: %s x %s (sep=%.2f, rev=%d)\n",
              nrow(selected_pairs), up_name, down_name,
              candidates$separation[k], candidates$reversal[k]))
  
  # Stop when max size reached
  if (nrow(selected_pairs) >= PARAMS$max_panel_size) break
}

if (nrow(selected_pairs) < PARAMS$max_panel_size) {
  cat(sprintf("\n=== Selected %d pairs (all available with separation > 0) ===\n\n", 
              nrow(selected_pairs)))
} else {
  cat(sprintf("\n=== Selected %d pairs (reached max panel size) ===\n\n", 
              nrow(selected_pairs)))
}

# Assign pair IDs
selected_pairs$pair_id <- paste0("P", 1:nrow(selected_pairs))
rownames(selected_pairs) <- selected_pairs$pair_id

# -----------------------------------------------------------------------------
# Display Final Panel
# -----------------------------------------------------------------------------
cat("=== Final Panel ===\n")
print(selected_pairs[, c("pair_id", "up_name", "down_name", "separation",
                         "reversal", "r0_consistency", "dz_total")])

# Correlation matrix (Spearman)
if (nrow(selected_pairs) > 1) {
  colnames(selected_r_matrix) <- selected_pairs$pair_id
  cor_matrix <- cor(selected_r_matrix, method = "spearman")
  
  cat("\n=== Pair Correlation Matrix (Spearman) ===\n")
  print(round(cor_matrix, 2))
  
  max_cor <- max(abs(cor_matrix[upper.tri(cor_matrix)]))
  cat(sprintf("\nMax off-diagonal correlation: %.2f\n", max_cor))
  
  # Warn if any high correlations slipped through
  if (max_cor >= PARAMS$correlation_threshold) {
    cat("!!! WARNING: High correlation pairs detected !!!\n")
  }
}

# -----------------------------------------------------------------------------
# Phase 6: Evaluate Separation
# -----------------------------------------------------------------------------
cat("\n=== Phase 6: Separation Evaluation ===\n")

# Compute reversal count for each sample
# Note: |r| < dead_zone is treated as NON-REVERSAL (conservative)
compute_sample_reversals <- function(selected_pairs, log2_tpm, 
                                     r0_samples, r1_samples, dead_zone) {
  
  n_pairs <- nrow(selected_pairs)
  all_samples <- c(r0_samples, r1_samples)
  
  # Matrix to store reversal status (1 = reversed, 0 = not reversed)
  reversal_matrix <- matrix(0, nrow = length(all_samples), ncol = n_pairs)
  rownames(reversal_matrix) <- all_samples
  colnames(reversal_matrix) <- selected_pairs$pair_id
  
  # Matrix to store dead zone status (for reference)
  dz_matrix <- matrix(FALSE, nrow = length(all_samples), ncol = n_pairs)
  rownames(dz_matrix) <- all_samples
  colnames(dz_matrix) <- selected_pairs$pair_id
  
  for (k in 1:n_pairs) {
    up_g <- selected_pairs$up_gene[k]
    down_g <- selected_pairs$down_gene[k]
    
    # r values for all samples
    r_all <- log2_tpm[up_g, all_samples] - log2_tpm[down_g, all_samples]
    
    # R0 majority sign (baseline for this pair)
    r_r0 <- r_all[r0_samples]
    majority_sign <- ifelse(sum(r_r0 > 0) > sum(r_r0 < 0), 1, -1)
    
    # Dead zone check
    in_dead_zone <- abs(r_all) < dead_zone
    dz_matrix[, k] <- in_dead_zone
    
    # Reversal: sign differs from majority AND not in dead zone
    # If in dead zone -> treated as non-reversal (conservative)
    is_reversed <- (sign(r_all) != majority_sign) & !in_dead_zone
    reversal_matrix[, k] <- as.numeric(is_reversed)
  }
  
  return(list(reversal = reversal_matrix, dead_zone = dz_matrix))
}

eval_result <- compute_sample_reversals(selected_pairs, log2_tpm, 
                                        r0_samples, r1_samples,
                                        PARAMS$dead_zone_threshold)
reversal_matrix <- eval_result$reversal
dz_matrix <- eval_result$dead_zone

# R0 reversal counts
r0_reversals <- rowSums(reversal_matrix[r0_samples, , drop = FALSE])
cat("\n--- R0 Sample Reversal Counts ---\n")
cat("(Expected: low, as these are sporadic cases)\n")
print(sort(r0_reversals, decreasing = TRUE))

# R1 reversal counts  
r1_reversals <- rowSums(reversal_matrix[r1_samples, , drop = FALSE])
cat("\n--- R1 Sample Reversal Counts ---\n")
cat("(Expected: high for radiation-induced, variable for mixed group)\n")
print(sort(r1_reversals, decreasing = TRUE))

# Summary statistics
cat("\n--- Summary ---\n")
cat(sprintf("R0: range=%d-%d, mean=%.1f, median=%.0f\n", 
            min(r0_reversals), max(r0_reversals), 
            mean(r0_reversals), median(r0_reversals)))
cat(sprintf("R1: range=%d-%d, mean=%.1f, median=%.0f\n", 
            min(r1_reversals), max(r1_reversals),
            mean(r1_reversals), median(r1_reversals)))

# Dead zone summary
r0_dz_count <- rowSums(dz_matrix[r0_samples, , drop = FALSE])
r1_dz_count <- rowSums(dz_matrix[r1_samples, , drop = FALSE])
cat(sprintf("\nDead zone (|r| < %.3f) occurrences:\n", PARAMS$dead_zone_threshold))
cat(sprintf("  R0: mean=%.1f pairs/sample in DZ\n", mean(r0_dz_count)))
cat(sprintf("  R1: mean=%.1f pairs/sample in DZ\n", mean(r1_dz_count)))

# -----------------------------------------------------------------------------
# Phase 7: Threshold Setting
# -----------------------------------------------------------------------------
cat("\n=== Phase 7: Threshold Setting ===\n")

n_pairs <- nrow(selected_pairs)

# Data-driven threshold: max(R0 reversal) + 1 ensures 100% R0 specificity
max_r0_reversal <- max(r0_reversals)
threshold_T <- max_r0_reversal + 1

cat(sprintf("\n--- Data-Driven Threshold Selection ---\n"))
cat(sprintf("R0 max reversal: %d\n", max_r0_reversal))
cat(sprintf("Selected threshold T = %d (max R0 reversal + 1)\n", threshold_T))
cat("Rationale: Minimum threshold ensuring 100%% R0 specificity\n\n")

cat("--- Sensitivity Analysis ---\n")
cat("Threshold T | R0 Panel(+) | R1 Panel(-) | R0 Specificity | R1 Sensitivity\n")
cat("----------- | ----------- | ----------- | -------------- | --------------\n")

threshold_results <- data.frame()

for (T in 1:min(10, n_pairs)) {
  r0_positive <- sum(r0_reversals >= T)  # False positives
  r1_negative <- sum(r1_reversals < T)   # "Negatives" (may include sporadic in R1)
  r0_specificity <- (n_r0 - r0_positive) / n_r0 * 100
  r1_sensitivity <- (n_r1 - r1_negative) / n_r1 * 100  # Note: R1 is mixed
  
  marker <- ifelse(T == threshold_T, " <-- selected", "")
  
  cat(sprintf("T = %d       | %d/%d         | %d/%d         | %.1f%%          | %.1f%%%s\n", 
              T, r0_positive, n_r0, r1_negative, n_r1, 
              r0_specificity, r1_sensitivity, marker))
  
  threshold_results <- rbind(threshold_results, data.frame(
    T = T,
    r0_positive = r0_positive,
    r1_negative = r1_negative,
    r0_specificity = r0_specificity,
    r1_sensitivity = r1_sensitivity
  ))
}

cat(sprintf("\n--- Final Threshold: T = %d ---\n", threshold_T))
cat("Note: R1 'sensitivity' is approximate as R1 contains mixed cases\n")

# Final classification
cat("\n--- Classification with T = %d ---\n", threshold_T)

r0_panel_positive <- r0_reversals >= threshold_T
r1_panel_positive <- r1_reversals >= threshold_T

cat(sprintf("R0 (sporadic): %d Panel(+), %d Panel(-)\n",
            sum(r0_panel_positive), sum(!r0_panel_positive)))
cat(sprintf("R1 (high POC): %d Panel(+), %d Panel(-)\n",
            sum(r1_panel_positive), sum(!r1_panel_positive)))

# Identify notable samples
cat("\n--- Notable Samples ---\n")

# R0 with high reversal (potential false positives)
if (any(r0_reversals >= threshold_T - 1)) {
  cat("\nR0 samples near/above threshold:\n")
  high_r0 <- r0_samples[r0_reversals >= threshold_T - 1]
  for (s in high_r0) {
    status <- ifelse(r0_reversals[s] >= threshold_T, "Panel(+) - CAUTION", "Panel(-)")
    cat(sprintf("  %s: %d/%d reversals -> %s\n", s, r0_reversals[s], n_pairs, status))
  }
}

# R1 with low reversal (sporadic within R1)
if (any(r1_reversals < threshold_T + 1)) {
  cat("\nR1 samples near/below threshold:\n")
  low_r1 <- r1_samples[r1_reversals < threshold_T + 1]
  for (s in low_r1) {
    status <- ifelse(r1_reversals[s] >= threshold_T, "Panel(+)", "Panel(-) - likely sporadic")
    cat(sprintf("  %s: %d/%d reversals -> %s\n", s, r1_reversals[s], n_pairs, status))
  }
}

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------
cat("\n=== Saving Results ===\n")

output <- list(
  # Parameters
  params = PARAMS,
  threshold_T = threshold_T,
  
  # Panel definition
  selected_pairs = selected_pairs,
  correlation_matrix = cor_matrix,
  
  # Evaluation matrices
  reversal_matrix = reversal_matrix,
  dead_zone_matrix = dz_matrix,
  
  # Sample-level results
  r0_reversals = r0_reversals,
  r1_reversals = r1_reversals,
  r0_panel_positive = r0_panel_positive,
  r1_panel_positive = r1_panel_positive,
  
  # Threshold analysis
  threshold_results = threshold_results,
  
  # Data references
  log2_tpm = log2_tpm,
  r0_samples = r0_samples,
  r1_samples = r1_samples,
  
  # Metadata
  timestamp = Sys.time(),
  version = "v1.0"
)

saveRDS(output, file.path(paths$processed, "reo_final_panel.rds"))
cat(sprintf("Saved: %s\n", file.path(paths$processed, "reo_final_panel.rds")))

# Panel summary CSV
panel_csv <- selected_pairs[, c("pair_id", "up_name", "down_name", 
                                "up_gene", "down_gene",
                                "dz_total", "r0_consistency", 
                                "reversal", "wilson_lower")]
write.csv(panel_csv, file.path(paths$output, "reo_final_panel.csv"), row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_final_panel.csv")))

# Sample classification CSV
sample_results <- data.frame(
  sample_id = c(r0_samples, r1_samples),
  group = c(rep("R0", n_r0), rep("R1", n_r1)),
  reversal_count = c(r0_reversals, r1_reversals),
  panel_result = c(ifelse(r0_panel_positive, "Panel(+)", "Panel(-)"),
                   ifelse(r1_panel_positive, "Panel(+)", "Panel(-)"))
)
write.csv(sample_results, file.path(paths$output, "reo_sample_classification.csv"), 
          row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_sample_classification.csv")))

cat("\n=== Phase 5-7 Complete ===\n")
cat(sprintf("Final panel: %d pairs\n", nrow(selected_pairs)))
cat(sprintf("Voting threshold: T = %d\n", threshold_T))
cat(sprintf("R0 specificity: %.1f%%\n", 
            threshold_results$r0_specificity[threshold_results$T == threshold_T]))
cat(sprintf("R1 Panel(+) rate: %.1f%%\n", 
            threshold_results$r1_sensitivity[threshold_results$T == threshold_T]))