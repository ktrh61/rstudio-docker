# 24_reo_step4_dedup.R - REO Step 4: Redundancy Removal (Deduplication)
# Purpose: Remove redundancy via maximum bipartite matching and correlation clustering
# Input: reo_step3b_data.rds
# Output: reo_step4_data.rds
# Version: v1.7 - Added detailed debugging for bipartite matching extraction
# Date: 2025-01-26
# Note: Lexicographic ordering for representative selection

source("analysis_v7/setup.R")

cat("\n=== REO STEP 4: Redundancy Removal (Deduplication) ===\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)  # For maximum matching
})

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

# Extract data
passed_pairs <- step3b_data$passed_pairs
reo_r0_passed <- step3b_data$reo_r0_passed
reo_r1_passed <- step3b_data$reo_r1_passed
r0_samples <- step3b_data$r0_samples
r1_samples <- step3b_data$r1_samples

cat(sprintf("  Input pairs from Step 3b: %d\n", nrow(passed_pairs)))

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
# 4.3 Step 4.1: Gene-unique maximum bipartite matching
# --------------------------------------------------------------------------

cat("\n--- 4.3 Gene-unique maximum bipartite matching (1 gene = 1 pair) ---\n")

# Create bipartite graph for maximum matching
# Use igraph's max_bipartite_match for optimal solution

# First, calculate lexicographic scores for prioritization
# We'll use weights to prefer better pairs in the matching
passed_pairs$lex_rank1 <- rank(-passed_pairs$a2_value, na.last = TRUE)  # Higher a[2] is better
passed_pairs$lex_rank2 <- rank(-passed_pairs$r1_reversal_rate, na.last = TRUE)  # Higher reversal is better
passed_pairs$lex_rank3 <- rank(passed_pairs$missing_rate, na.last = TRUE)  # Lower missing is better
passed_pairs$lex_rank4 <- rank(paste0(passed_pairs$gene_up, "_", passed_pairs$gene_down))  # Alphabetical

# Combine ranks with appropriate weighting to ensure lexicographic ordering
n_pairs <- nrow(passed_pairs)
passed_pairs$lex_score <- (
  passed_pairs$lex_rank1 * (n_pairs^3) +
    passed_pairs$lex_rank2 * (n_pairs^2) +
    passed_pairs$lex_rank3 * n_pairs +
    passed_pairs$lex_rank4
)

# Create bipartite graph
up_genes <- unique(passed_pairs$gene_up)
down_genes <- unique(passed_pairs$gene_down)
n_up <- length(up_genes)
n_down <- length(down_genes)

cat(sprintf("  Unique up genes in input: %d\n", n_up))
cat(sprintf("  Unique down genes in input: %d\n", n_down))
cat(sprintf("  Theoretical maximum matches: %d\n", min(n_up, n_down)))

# Create edge list with weights (lower lex_score = higher priority = higher weight)
edges <- data.frame(
  from = match(passed_pairs$gene_up, up_genes),
  to = match(passed_pairs$gene_down, down_genes) + n_up,
  weight = max(passed_pairs$lex_score) - passed_pairs$lex_score + 1,  # Invert for weight
  pair_idx = 1:nrow(passed_pairs)
)

cat(sprintf("  Total edges in bipartite graph: %d\n", nrow(edges)))

# Create bipartite graph
g <- graph_from_data_frame(edges, directed = FALSE, 
                           vertices = data.frame(
                             name = 1:(n_up + n_down),
                             type = c(rep(TRUE, n_up), rep(FALSE, n_down))
                           ))

# Find maximum weight bipartite matching
matching <- max_bipartite_match(g, weights = E(g)$weight)

cat(sprintf("  Matching algorithm completed. Matching size: %d\n", matching$matching_size))

# Debug the structure of matching result
cat("  DEBUG: Checking matching structure...\n")
cat(sprintf("    Length of matching$matching: %d\n", length(matching$matching)))
cat(sprintf("    Class of matching$matching: %s\n", class(matching$matching)))
if (length(matching$matching) > 0) {
  cat(sprintf("    First few elements: %s\n", 
              paste(head(matching$matching, 5), collapse=", ")))
  cat(sprintf("    Names (if any): %s\n", 
              paste(head(names(matching$matching), 5), collapse=", ")))
}

# Extract matched pairs
# matching$matching can be indexed by vertex ID or be a named vector
matched_pairs_list <- list()

# Try to extract matches - handle different possible structures
if (matching$matching_size > 0 && length(matching$matching) > 0) {
  # Check if it's a named vector
  if (!is.null(names(matching$matching))) {
    # Named vector case - names are vertex IDs as strings
    for (i in seq_along(matching$matching)) {
      v_name <- names(matching$matching)[i]
      v_idx <- as.integer(v_name)
      matched_to <- matching$matching[i]
      
      # Check if this is an up vertex that's matched
      if (!is.na(matched_to) && v_idx <= n_up) {
        matched_vertex_num <- as.integer(matched_to)
        
        if (matched_vertex_num > n_up && matched_vertex_num <= (n_up + n_down)) {
          down_idx <- matched_vertex_num - n_up
          up_gene <- up_genes[v_idx]
          down_gene <- down_genes[down_idx]
          
          pair_candidates <- which(
            passed_pairs$gene_up == up_gene & 
              passed_pairs$gene_down == down_gene
          )
          
          if (length(pair_candidates) > 0) {
            best_pair <- pair_candidates[which.min(passed_pairs$lex_score[pair_candidates])]
            matched_pairs_list[[length(matched_pairs_list) + 1]] <- best_pair
          }
        }
      }
    }
  } else {
    # Unnamed vector case - indexed directly
    for (v_idx in 1:min(n_up, length(matching$matching))) {
      matched_to <- matching$matching[v_idx]
      
      if (!is.na(matched_to)) {
        matched_vertex_num <- as.integer(matched_to)
        
        if (matched_vertex_num > n_up && matched_vertex_num <= (n_up + n_down)) {
          down_idx <- matched_vertex_num - n_up
          up_gene <- up_genes[v_idx]
          down_gene <- down_genes[down_idx]
          
          pair_candidates <- which(
            passed_pairs$gene_up == up_gene & 
              passed_pairs$gene_down == down_gene
          )
          
          if (length(pair_candidates) > 0) {
            best_pair <- pair_candidates[which.min(passed_pairs$lex_score[pair_candidates])]
            matched_pairs_list[[length(matched_pairs_list) + 1]] <- best_pair
          }
        }
      }
    }
  }
}

# Convert to vector
selected_pairs <- unlist(matched_pairs_list)

# Debug output
cat(sprintf("  Successfully extracted %d pairs from matching\n", length(selected_pairs)))

# Check if any pairs were matched
if (length(selected_pairs) == 0) {
  cat("  WARNING: No pairs could be matched in bipartite graph!\n")
  cat("  This may indicate a structural issue with the input pairs.\n")
  selected_pairs <- integer(0)  # Ensure it's an empty integer vector
}

# Extract selected pairs (handle empty case)
if (length(selected_pairs) > 0) {
  gene_unique_pairs <- passed_pairs[selected_pairs, ]
  gene_unique_reo_r0 <- reo_r0_passed[selected_pairs, ]
  gene_unique_reo_r1 <- reo_r1_passed[selected_pairs, ]
} else {
  # Create empty data frames with correct structure
  gene_unique_pairs <- passed_pairs[integer(0), ]
  gene_unique_reo_r0 <- reo_r0_passed[integer(0), ]
  gene_unique_reo_r1 <- reo_r1_passed[integer(0), ]
}

cat(sprintf("  Input pairs: %d\n", nrow(passed_pairs)))
cat(sprintf("  Maximum bipartite matching found: %d pairs\n", length(selected_pairs)))
cat(sprintf("  After bipartite matching: %d\n", nrow(gene_unique_pairs)))
cat(sprintf("  Retention rate: %.1f%%\n", 
            nrow(gene_unique_pairs) / nrow(passed_pairs) * 100))
cat(sprintf("  Unique up genes in result: %d\n", length(unique(gene_unique_pairs$gene_up))))
cat(sprintf("  Unique down genes in result: %d\n", length(unique(gene_unique_pairs$gene_down))))

# --------------------------------------------------------------------------
# 4.4 Step 4.2: Correlation-based clustering
# --------------------------------------------------------------------------

cat("\n--- 4.4 Correlation-based clustering (|ρ| > 0.9) ---\n")

# Check if we have any pairs to cluster
if (nrow(gene_unique_pairs) == 0) {
  cat("  No pairs available for correlation clustering.\n")
  final_pairs <- gene_unique_pairs
  final_reo_r0 <- gene_unique_reo_r0
  final_reo_r1 <- gene_unique_reo_r1
} else {
  # Calculate pairwise Spearman correlations
  # Combine R0 and R1 REO values for correlation calculation
  combined_reo <- cbind(gene_unique_reo_r0, gene_unique_reo_r1)
  n_pairs_unique <- nrow(gene_unique_pairs)
  
  # Initialize correlation matrix
  cor_matrix <- matrix(0, nrow = n_pairs_unique, ncol = n_pairs_unique)
  diag(cor_matrix) <- 1
  
  # Calculate pairwise correlations
  if (n_pairs_unique > 1) {
    for (i in 1:(n_pairs_unique - 1)) {
      for (j in (i + 1):n_pairs_unique) {
        # Spearman correlation across all samples
        cor_value <- cor(combined_reo[i, ], combined_reo[j, ], 
                         method = "spearman", use = "complete.obs")
        cor_matrix[i, j] <- cor_value
        cor_matrix[j, i] <- cor_value
      }
    }
    
    # Find pairs with high correlation (|ρ| > 0.9)
    high_cor_threshold <- 0.9
    high_cor_pairs <- which(abs(cor_matrix) > high_cor_threshold & 
                              upper.tri(cor_matrix), arr.ind = TRUE)
    
    cat(sprintf("  Pairs with |ρ| > %.1f: %d\n", high_cor_threshold, nrow(high_cor_pairs)))
    
    if (nrow(high_cor_pairs) > 0) {
      # Create graph with all nodes explicitly numbered 1:n_pairs_unique
      # This ensures proper index correspondence
      g <- make_empty_graph(n = n_pairs_unique, directed = FALSE)
      
      # Add edges for high correlation pairs
      edges_to_add <- as.vector(t(high_cor_pairs))
      g <- add_edges(g, edges_to_add)
      
      # Find connected components (correlation clusters)
      clusters <- components(g)
      n_clusters <- clusters$no
      
      cat(sprintf("  Number of correlation clusters: %d\n", n_clusters))
      cat(sprintf("  Singleton clusters: %d\n", sum(clusters$csize == 1)))
      cat(sprintf("  Multi-member clusters: %d\n", sum(clusters$csize > 1)))
      
      # Select representative from each cluster
      cluster_representatives <- integer(0)  # Initialize as empty vector
      
      for (cluster_id in 1:n_clusters) {
        cluster_members <- which(clusters$membership == cluster_id)
        
        if (length(cluster_members) == 1) {
          # Singleton cluster
          cluster_representatives <- c(cluster_representatives, cluster_members[1])
        } else {
          # Multi-member cluster - select based on lexicographic ordering
          cluster_subset <- gene_unique_pairs[cluster_members, ]
          
          # Recompute lexicographic ranks within cluster
          cluster_subset$cluster_lex_rank1 <- rank(-cluster_subset$a2_value, na.last = TRUE)
          cluster_subset$cluster_lex_rank2 <- rank(-cluster_subset$r1_reversal_rate, na.last = TRUE)
          cluster_subset$cluster_lex_rank3 <- rank(cluster_subset$missing_rate, na.last = TRUE)
          cluster_subset$cluster_lex_rank4 <- rank(paste0(cluster_subset$gene_up, "_", 
                                                          cluster_subset$gene_down))
          
          # Combine ranks
          n_cluster <- nrow(cluster_subset)
          cluster_subset$cluster_lex_score <- (
            cluster_subset$cluster_lex_rank1 * (n_cluster^3) +
              cluster_subset$cluster_lex_rank2 * (n_cluster^2) +
              cluster_subset$cluster_lex_rank3 * n_cluster +
              cluster_subset$cluster_lex_rank4
          )
          
          # Select best (lowest score)
          best_idx <- which.min(cluster_subset$cluster_lex_score)
          cluster_representatives <- c(cluster_representatives, cluster_members[best_idx])
        }
      }
      
      # Sort representatives to maintain consistent ordering
      cluster_representatives <- sort(cluster_representatives)
      
      # Extract final deduplicated pairs
      final_pairs <- gene_unique_pairs[cluster_representatives, ]
      final_reo_r0 <- gene_unique_reo_r0[cluster_representatives, ]
      final_reo_r1 <- gene_unique_reo_r1[cluster_representatives, ]
      
    } else {
      # No high correlations found - keep all pairs
      cat("  No highly correlated pairs found. Keeping all gene-unique pairs.\n")
      final_pairs <- gene_unique_pairs
      final_reo_r0 <- gene_unique_reo_r0
      final_reo_r1 <- gene_unique_reo_r1
    }
    
  } else {
    # Only one pair - no clustering needed
    cat("  Only one pair available. No clustering needed.\n")
    final_pairs <- gene_unique_pairs
    final_reo_r0 <- gene_unique_reo_r0
    final_reo_r1 <- gene_unique_reo_r1
  }
  
  cat(sprintf("\n  Final deduplicated pairs: %d\n", nrow(final_pairs)))
  cat(sprintf("  Overall retention from Step 3b: %.1f%%\n", 
              nrow(final_pairs) / nrow(passed_pairs) * 100))
}  # End of check for non-empty gene_unique_pairs

# --------------------------------------------------------------------------
# 4.5 Add selection metadata
# --------------------------------------------------------------------------

cat("\n--- 4.5 Adding selection metadata ---\n")

# Check if we have any pairs to process
if (nrow(final_pairs) > 0) {
  # Record selection reasons
  final_pairs$selection_step <- "correlation_clustering"
  final_pairs$dedup_method <- "bipartite_matching_lexicographic"
  
  # Add cluster information if clustering was performed
  if (exists("clusters") && exists("high_cor_pairs") && nrow(high_cor_pairs) > 0) {
    # cluster_representatives contains indices into gene_unique_pairs
    cluster_membership <- clusters$membership
    cluster_sizes <- numeric(length(cluster_representatives))
    for (i in seq_along(cluster_representatives)) {
      rep_idx <- cluster_representatives[i]
      cluster_id <- cluster_membership[rep_idx]
      cluster_sizes[i] <- clusters$csize[cluster_id]
    }
    final_pairs$cluster_size <- cluster_sizes
    final_pairs$from_multi_cluster <- cluster_sizes > 1
  } else {
    final_pairs$cluster_size <- 1
    final_pairs$from_multi_cluster <- FALSE
  }
  
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

step4_data <- list(
  # Configuration (pass through)
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
    correlation_threshold = ifelse(exists("high_cor_threshold"), high_cor_threshold, NA),
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
cat(sprintf("After bipartite matching: %d (%.1f%%)\n", 
            nrow(gene_unique_pairs),
            nrow(gene_unique_pairs) / nrow(passed_pairs) * 100))
cat(sprintf("After correlation clustering: %d (%.1f%%)\n", 
            nrow(final_pairs),
            nrow(final_pairs) / nrow(passed_pairs) * 100))
cat(sprintf("Final retention rate: %.1f%%\n", 
            nrow(final_pairs) / nrow(passed_pairs) * 100))
cat(sprintf("Timestamp: %s\n", format(step4_data$timestamp)))

# Final message
if (step4_data$n_final >= 10) {
  cat(sprintf("\n✓ Ready for Step 5: Panel Creation (%d pairs available)\n", 
              step4_data$n_final))
} else if (step4_data$n_final > 0) {
  cat(sprintf("\n⚠ Only %d pairs available. Consider relaxing thresholds if more pairs needed for panel.\n",
              step4_data$n_final))
} else {
  cat("\n⚠ No pairs survived deduplication. Review thresholds and data quality.\n")
}