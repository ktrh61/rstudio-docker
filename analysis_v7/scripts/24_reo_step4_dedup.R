# 24_reo_step4_dedup.R - REO Step 4: Redundancy Removal (Deduplication)
# Purpose: Remove redundancy via maximum bipartite matching and correlation clustering
# Input: reo_step3b_data.rds
# Output: reo_step4_data.rds
# Version: v1.8 - CONFIG-based thresholds and controlled debug output
# Date: 2025-01-26 (Updated)
# Note: Lexicographic ordering for representative selection

source("analysis_v7/setup.R")

cat("\n=== REO STEP 4: Redundancy Removal (Deduplication) ===\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  # igraph not needed anymore - using greedy selection instead of bipartite matching
})

# --------------------------------------------------------------------------
# 1.0 Check for CONFIG parameters (defensive)
# --------------------------------------------------------------------------
# If CONFIG doesn't have correlation parameters, set defaults
# This ensures backward compatibility
DEFAULT_COR_THRESHOLD <- 0.9
DEFAULT_COR_METHOD <- "spearman"

# --------------------------------------------------------------------------
# 4.1 Load Step 3b data
# --------------------------------------------------------------------------

cat("--- 4.1 Loading Step 3b data ---\n")

step3b_data_path <- paste0(paths$processed, "reo_step3b_data.rds")
if (!file.exists(step3b_data_path)) {
  stop("ERROR: reo_step3b_data.rds not found. Please run 23_reo_step3b_proxy.R first.")
}

step3b_data <- readRDS(step3b_data_path)
cat("  Step 3b data loaded successfully\n")

# 薄い契約：最小限の存在チェック
stopifnot(all(c("passed_pairs", "reo_r0_passed", "reo_r1_passed", "config") %in% names(step3b_data)))

# CONFIGはRDSから取得
CONFIG <- step3b_data$config

# Set correlation parameters from CONFIG or defaults
cor_threshold <- ifelse(!is.null(CONFIG$cor_threshold), CONFIG$cor_threshold, DEFAULT_COR_THRESHOLD)
cor_method <- ifelse(!is.null(CONFIG$cor_method), CONFIG$cor_method, DEFAULT_COR_METHOD)

# Extract data
passed_pairs <- step3b_data$passed_pairs
reo_r0_passed <- step3b_data$reo_r0_passed
reo_r1_passed <- step3b_data$reo_r1_passed
r0_samples <- step3b_data$r0_samples
r1_samples <- step3b_data$r1_samples

cat(sprintf("  Input pairs from Step 3b: %d\n", nrow(passed_pairs)))
cat(sprintf("  Correlation settings (使用値=%s): threshold=%.2f\n", cor_method, cor_threshold))

# --------------------------------------------------------------------------
# 4.2 Prepare data for redundancy removal
# --------------------------------------------------------------------------

cat("\n--- 4.2 Preparing data for redundancy removal ---\n")

# Calculate a[2] (order statistic matching Step 3a definition) for each pair
# a[2] is the (n_R0-1)th smallest absolute value after dead zone filter
# This allows for at most 1 exception in R0 group
calculate_a2 <- function(reo_values, dead_zone, n_r0) {
  # Remove any NA values first (defensive programming)
  reo_values <- reo_values[!is.na(reo_values)]
  
  # Apply dead zone filter
  valid_values <- abs(reo_values[abs(reo_values) >= dead_zone])
  
  # Need at least (n_r0 - 1) valid samples to calculate a[2]
  if (length(valid_values) >= (n_r0 - 1)) {
    # Sort and take the (n_r0-1)th value (allowing 1 exception)
    # For n_r0=11, this takes the 10th smallest value
    # sort() removes NA by default, but we already removed them
    return(sort(valid_values, na.last = NA)[n_r0 - 1])
  } else {
    # Not enough valid samples - same as Step 3a rejection
    return(NA)
  }
}

# Calculate a[2] for all pairs
a2_values <- numeric(nrow(passed_pairs))
n_r0 <- length(r0_samples)  # Get R0 sample count
for (i in 1:nrow(passed_pairs)) {
  a2_values[i] <- calculate_a2(reo_r0_passed[i, ], CONFIG$dead_zone, n_r0)
}

# Add a[2] to passed_pairs dataframe
passed_pairs$a2_value <- a2_values

# Calculate missing rate for each pair
missing_rate <- function(r0_row, r1_row, dead_zone) {
  r0_invalid <- sum(abs(r0_row) < dead_zone)
  r1_invalid <- sum(abs(r1_row) < dead_zone)
  total_samples <- length(r0_row) + length(r1_row)
  return((r0_invalid + r1_invalid) / total_samples)
}

# Calculate missing rates
missing_rates <- numeric(nrow(passed_pairs))
for (i in 1:nrow(passed_pairs)) {
  missing_rates[i] <- missing_rate(reo_r0_passed[i, ], reo_r1_passed[i, ], CONFIG$dead_zone)
}

passed_pairs$missing_rate <- missing_rates

cat(sprintf("  Calculated metrics for %d pairs\n", nrow(passed_pairs)))
cat(sprintf("  a[2] range (strength, %dth order stat): %.2f - %.2f\n", 
            n_r0 - 1,
            min(a2_values, na.rm = TRUE), max(a2_values, na.rm = TRUE)))
cat(sprintf("  Missing rate range: %.1f%% - %.1f%%\n", 
            min(missing_rates) * 100, max(missing_rates) * 100))

# --------------------------------------------------------------------------
# 4.3 Step 4.1: Gene-unique selection (1 gene = 1 pair) 
# --------------------------------------------------------------------------

cat("\n--- 4.3 Gene-unique selection (1 gene = 1 pair) ---\n")

# Calculate lexicographic scores for prioritization
passed_pairs$lex_rank1 <- rank(-passed_pairs$a2_value, na.last = TRUE)  # Higher a[2] is better
passed_pairs$lex_rank2 <- rank(-passed_pairs$r1_reversal_rate, na.last = TRUE)  # Higher reversal is better
passed_pairs$lex_rank3 <- rank(passed_pairs$missing_rate, na.last = TRUE)  # Lower missing is better
passed_pairs$lex_rank4 <- rank(paste0(passed_pairs$gene_up, "_", passed_pairs$gene_down))  # Alphabetical

# Combine into single score (lexicographic)
n_pairs <- nrow(passed_pairs)
passed_pairs$lex_score <- (
  passed_pairs$lex_rank1 * (n_pairs^3) +
    passed_pairs$lex_rank2 * (n_pairs^2) +
    passed_pairs$lex_rank3 * n_pairs +
    passed_pairs$lex_rank4
)

# Simple greedy approach: Sort by priority and select non-overlapping pairs
# This is much simpler than bipartite matching and achieves the same goal

# Sort pairs by lexicographic score (best first)
sorted_indices <- order(passed_pairs$lex_score)
sorted_pairs <- passed_pairs[sorted_indices, ]
sorted_reo_r0 <- reo_r0_passed[sorted_indices, ]
sorted_reo_r1 <- reo_r1_passed[sorted_indices, ]

# Track which genes have been used
used_up_genes <- character(0)
used_down_genes <- character(0)
selected_indices <- integer(0)

# Greedy selection: iterate through sorted pairs
for (i in 1:nrow(sorted_pairs)) {
  current_up <- sorted_pairs$gene_up[i]
  current_down <- sorted_pairs$gene_down[i]
  
  # Check if either gene has been used
  if (!current_up %in% used_up_genes && !current_down %in% used_down_genes) {
    # Select this pair
    selected_indices <- c(selected_indices, i)
    used_up_genes <- c(used_up_genes, current_up)
    used_down_genes <- c(used_down_genes, current_down)
  }
}

# Extract selected pairs
if (length(selected_indices) > 0) {
  gene_unique_pairs <- sorted_pairs[selected_indices, ]
  gene_unique_reo_r0 <- sorted_reo_r0[selected_indices, ]
  gene_unique_reo_r1 <- sorted_reo_r1[selected_indices, ]
} else {
  gene_unique_pairs <- sorted_pairs[integer(0), ]
  gene_unique_reo_r0 <- sorted_reo_r0[integer(0), ]
  gene_unique_reo_r1 <- sorted_reo_r1[integer(0), ]
}

# Report results
cat(sprintf("  Input pairs: %d\n", nrow(passed_pairs)))
cat(sprintf("  Unique up genes: %d\n", length(unique(passed_pairs$gene_up))))
cat(sprintf("  Unique down genes: %d\n", length(unique(passed_pairs$gene_down))))
cat(sprintf("  Gene-unique pairs selected: %d\n", nrow(gene_unique_pairs)))
cat(sprintf("  Retention rate: %.1f%%\n", 
            nrow(gene_unique_pairs) / nrow(passed_pairs) * 100))

# Debug output
if (!is.null(CONFIG$VERBOSE) && CONFIG$VERBOSE && nrow(gene_unique_pairs) > 0) {
  cat("\n  DEBUG: Top 5 selected pairs:\n")
  for (i in 1:min(5, nrow(gene_unique_pairs))) {
    cat(sprintf("    %2d. %s vs %s (score=%.0f, a[2]=%.2f, rev=%.2f)\n",
                i,
                gene_unique_pairs$gene_up[i],
                gene_unique_pairs$gene_down[i],
                gene_unique_pairs$lex_score[i],
                gene_unique_pairs$a2_value[i],
                gene_unique_pairs$r1_reversal_rate[i]))
  }
  
  cat(sprintf("\n  Unique up genes in result: %d\n", 
              length(unique(gene_unique_pairs$gene_up))))
  cat(sprintf("  Unique down genes in result: %d\n", 
              length(unique(gene_unique_pairs$gene_down))))
}

# --------------------------------------------------------------------------
# 4.4 Step 4.2: Correlation-based clustering
# --------------------------------------------------------------------------

cat(sprintf("\n--- 4.4 Correlation-based clustering (使用値=%s, |ρ| > %.2f) ---\n", 
            cor_method, cor_threshold))

# Check if we have any pairs to cluster
if (nrow(gene_unique_pairs) == 0) {
  cat("  No pairs available for correlation clustering.\n")
  final_pairs <- gene_unique_pairs
  final_reo_r0 <- gene_unique_reo_r0
  final_reo_r1 <- gene_unique_reo_r1
} else if (nrow(gene_unique_pairs) == 1) {
  # Only one pair - no clustering needed
  cat("  Only one pair available. No clustering needed.\n")
  final_pairs <- gene_unique_pairs
  final_reo_r0 <- gene_unique_reo_r0
  final_reo_r1 <- gene_unique_reo_r1
} else {
  # Calculate pairwise correlations
  # Combine R0 and R1 REO values for correlation calculation
  combined_reo <- cbind(gene_unique_reo_r0, gene_unique_reo_r1)
  n_pairs_unique <- nrow(gene_unique_pairs)
  
  # Initialize correlation matrix
  cor_matrix <- matrix(0, nrow = n_pairs_unique, ncol = n_pairs_unique)
  diag(cor_matrix) <- 1
  
  # Calculate pairwise correlations
  for (i in 1:(n_pairs_unique - 1)) {
    for (j in (i + 1):n_pairs_unique) {
      # Correlation across all samples using CONFIG method
      cor_value <- cor(combined_reo[i, ], combined_reo[j, ], 
                       method = cor_method, use = "complete.obs")
      cor_matrix[i, j] <- cor_value
      cor_matrix[j, i] <- cor_value
    }
  }
  
  # Find pairs with high correlation
  high_cor_pairs <- which(abs(cor_matrix) > cor_threshold & 
                            upper.tri(cor_matrix), arr.ind = TRUE)
  
  cat(sprintf("  Pairs with |ρ| > %.2f: %d\n", cor_threshold, nrow(high_cor_pairs)))
  
  if (nrow(high_cor_pairs) == 0) {
    # No high correlations found - keep all pairs
    cat("  No highly correlated pairs found. Keeping all gene-unique pairs.\n")
    final_pairs <- gene_unique_pairs
    final_reo_r0 <- gene_unique_reo_r0
    final_reo_r1 <- gene_unique_reo_r1
  } else {
    # Simple clustering approach without igraph
    # Track which pairs belong to which cluster
    cluster_membership <- seq_len(n_pairs_unique)  # Start with each in own cluster
    
    # Merge clusters for highly correlated pairs
    for (i in 1:nrow(high_cor_pairs)) {
      idx1 <- high_cor_pairs[i, 1]
      idx2 <- high_cor_pairs[i, 2]
      
      # Merge clusters
      cluster1 <- cluster_membership[idx1]
      cluster2 <- cluster_membership[idx2]
      
      if (cluster1 != cluster2) {
        # Merge cluster2 into cluster1
        cluster_membership[cluster_membership == cluster2] <- cluster1
      }
    }
    
    # Renumber clusters to be sequential
    unique_clusters <- unique(cluster_membership)
    n_clusters <- length(unique_clusters)
    for (i in seq_along(unique_clusters)) {
      cluster_membership[cluster_membership == unique_clusters[i]] <- i
    }
    
    cat(sprintf("  Number of correlation clusters: %d\n", n_clusters))
    
    # Count cluster sizes
    cluster_sizes <- table(cluster_membership)
    cat(sprintf("  Singleton clusters: %d\n", sum(cluster_sizes == 1)))
    cat(sprintf("  Multi-member clusters: %d\n", sum(cluster_sizes > 1)))
    
    # Select representative from each cluster
    cluster_representatives <- integer(n_clusters)
    
    for (cluster_id in 1:n_clusters) {
      cluster_members <- which(cluster_membership == cluster_id)
      
      if (length(cluster_members) == 1) {
        # Singleton cluster
        cluster_representatives[cluster_id] <- cluster_members[1]
      } else {
        # Multi-member cluster - select based on lexicographic ordering
        cluster_subset <- gene_unique_pairs[cluster_members, ]
        
        # Find best member (lowest lex_score)
        best_idx <- cluster_members[which.min(cluster_subset$lex_score)]
        cluster_representatives[cluster_id] <- best_idx
      }
    }
    
    # Extract final deduplicated pairs
    final_pairs <- gene_unique_pairs[cluster_representatives, ]
    final_reo_r0 <- gene_unique_reo_r0[cluster_representatives, ]
    final_reo_r1 <- gene_unique_reo_r1[cluster_representatives, ]
    
    # Add cluster size information
    final_pairs$cluster_size <- as.numeric(cluster_sizes[cluster_membership[cluster_representatives]])
  }
  
  cat(sprintf("\n  Final deduplicated pairs: %d\n", nrow(final_pairs)))
  cat(sprintf("  Overall retention from Step 3b: %.1f%%\n", 
              nrow(final_pairs) / nrow(passed_pairs) * 100))
}

# --------------------------------------------------------------------------
# 4.5 Add selection metadata
# --------------------------------------------------------------------------

cat("\n--- 4.5 Adding selection metadata ---\n")

# Check if we have any pairs to process
if (nrow(final_pairs) > 0) {
  # Record selection reasons
  final_pairs$selection_step <- "correlation_clustering"
  final_pairs$dedup_method <- "greedy_selection_lexicographic"
  
  # Add cluster information if clustering was performed
  if (!"cluster_size" %in% names(final_pairs)) {
    final_pairs$cluster_size <- 1
  }
  final_pairs$from_multi_cluster <- final_pairs$cluster_size > 1
  
  cat(sprintf("  Pairs from multi-member clusters: %d\n", 
              sum(final_pairs$from_multi_cluster)))
  cat(sprintf("  Pairs from singleton clusters: %d\n", 
              sum(!final_pairs$from_multi_cluster)))
} else {
  cat("  WARNING: No pairs available for metadata annotation.\n")
  # Ensure final_pairs has the expected columns even when empty
  final_pairs$selection_step <- character(0)
  final_pairs$dedup_method <- character(0)
  final_pairs$cluster_size <- numeric(0)
  final_pairs$from_multi_cluster <- logical(0)
}

# --------------------------------------------------------------------------
# 4.6 Save Step 4 results
# --------------------------------------------------------------------------

cat("\n--- 4.6 Saving Step 4 results ---\n")

# Update CONFIG with actually used values (for traceability)
CONFIG$cor_threshold_used <- cor_threshold
CONFIG$cor_method_used <- cor_method

step4_data <- list(
  # Configuration (pass through with used values)
  config = CONFIG,
  
  # Final deduplicated results
  final_pairs = final_pairs,
  final_reo_r0 = final_reo_r0,
  final_reo_r1 = final_reo_r1,
  
  # Intermediate results for traceability
  gene_unique_pairs = gene_unique_pairs,
  n_gene_unique = nrow(gene_unique_pairs),
  
  # Summary statistics
  n_input = nrow(passed_pairs),
  n_after_gene_unique = nrow(gene_unique_pairs),
  n_final = nrow(final_pairs),
  retention_rate = nrow(final_pairs) / nrow(passed_pairs),
  
  # Deduplication metadata
  dedup_stats = list(
    gene_unique_retention = nrow(gene_unique_pairs) / nrow(passed_pairs),
    correlation_retention = nrow(final_pairs) / nrow(gene_unique_pairs),
    correlation_threshold_used = cor_threshold,
    correlation_method_used = cor_method,
    n_correlation_clusters = ifelse(exists("n_clusters"), n_clusters, NA)
  ),
  
  # Sample info
  r0_samples = r0_samples,
  r1_samples = r1_samples,
  
  # Metadata
  timestamp = Sys.time(),
  step = "Step 4 Complete"
)

saveRDS(step4_data, paste0(paths$processed, "reo_step4_data.rds"))
cat("  Step 4 data saved to: reo_step4_data.rds\n")

# Export final pairs to CSV (only if non-empty)
if (nrow(final_pairs) > 0) {
  write.csv(final_pairs[, c("pair_id", "gene_up", "gene_down", "a2_value", 
                            "r1_reversal_rate", "missing_rate", "cluster_size")],
            paste0(paths$output, "step4_deduplicated_pairs.csv"),
            row.names = FALSE)
  cat("  Deduplicated pairs exported to: step4_deduplicated_pairs.csv\n")
  
  # Export detailed tiebreak log
  tiebreak_log <- final_pairs[, c("pair_id", "gene_up", "gene_down",
                                  "lex_rank1", "lex_rank2", "lex_rank3", "lex_rank4",
                                  "lex_score", "selection_step", "cluster_size")]
  write.csv(tiebreak_log,
            paste0(paths$output, "step4_tiebreak_log.csv"),
            row.names = FALSE)
  cat("  Tiebreak log exported to: step4_tiebreak_log.csv\n")
} else {
  cat("  WARNING: No pairs to export to CSV files.\n")
}

# =============================================================================
# STEP 4 SUMMARY (固定フォーマット)
# =============================================================================

cat("\n=== STEP 4 SUMMARY ===\n")
cat(sprintf("Input pairs: %d\n", nrow(passed_pairs)))
cat(sprintf("After gene-unique selection: %d (%.1f%%)\n", 
            nrow(gene_unique_pairs),
            nrow(gene_unique_pairs) / nrow(passed_pairs) * 100))
cat(sprintf("After correlation clustering: %d (%.1f%%)\n", 
            nrow(final_pairs),
            nrow(final_pairs) / nrow(passed_pairs) * 100))
cat(sprintf("Final retention rate: %.1f%%\n", 
            nrow(final_pairs) / nrow(passed_pairs) * 100))
cat(sprintf("Correlation parameters: %s (threshold=%.2f)\n", cor_method, cor_threshold))
cat(sprintf("Timestamp: %s\n", format(step4_data$timestamp)))

# Final message
if (step4_data$n_final >= 10) {
  cat(sprintf("\n✔ Ready for Step 5: Panel Creation (%d pairs available)\n", 
              step4_data$n_final))
} else if (step4_data$n_final > 0) {
  cat(sprintf("\n⚠ Only %d pairs available. Consider relaxing thresholds if more pairs needed for panel.\n",
              step4_data$n_final))
} else {
  cat("\n⚠ No pairs survived deduplication. Review thresholds and data quality.\n")
}