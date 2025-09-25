# 24_reo_step4_dedup.R - REO Step 4: Redundancy Removal (Deduplication)
# Purpose: Remove redundancy via maximum bipartite matching and correlation clustering
# Input: reo_step3b_data.rds
# Output: reo_step4_data.rds
# Version: v2.0 - Complete rewrite with all red-level fixes
# Date: 2025-01-26 (Fully Revised)
# Note: Uses proper bipartite matching with simple indexing

source("analysis_v7/setup.R")

cat("\n=== REO STEP 4: Redundancy Removal (Deduplication) ===\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)  # Required for max_bipartite_match
})

# --------------------------------------------------------------------------
# 1.0 Configuration defaults and validation
# --------------------------------------------------------------------------
# Define defaults for correlation parameters
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

# Fix 赤-1: Robust parameter validation with NA/empty string protection
cor_threshold <- if (!is.null(CONFIG$cor_threshold) && 
                     !is.na(CONFIG$cor_threshold) && 
                     is.numeric(CONFIG$cor_threshold) &&
                     CONFIG$cor_threshold >= 0 && 
                     CONFIG$cor_threshold <= 1) {
  CONFIG$cor_threshold
} else {
  cat(sprintf("  WARNING: Invalid cor_threshold in CONFIG. Using default: %.2f\n", 
              DEFAULT_COR_THRESHOLD))
  DEFAULT_COR_THRESHOLD
}

cor_method <- if (!is.null(CONFIG$cor_method) && 
                  !is.na(CONFIG$cor_method) && 
                  nchar(CONFIG$cor_method) > 0 &&
                  CONFIG$cor_method %in% c("pearson", "spearman", "kendall")) {
  CONFIG$cor_method
} else {
  cat(sprintf("  WARNING: Invalid cor_method in CONFIG. Using default: %s\n", 
              DEFAULT_COR_METHOD))
  DEFAULT_COR_METHOD
}

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

# Fix 赤-3: Use a2_value from Step 3b instead of recalculating
if (!"a2_value" %in% names(passed_pairs)) {
  stop("ERROR: a2_value not found in Step 3b output. Data integrity issue.")
}

# Verify we have all required fields from Step 3b
required_fields <- c("pair_id", "gene_up", "gene_down", "a2_value", "r1_reversal_rate")
missing_fields <- setdiff(required_fields, names(passed_pairs))
if (length(missing_fields) > 0) {
  stop(sprintf("ERROR: Missing required fields from Step 3b: %s", 
               paste(missing_fields, collapse=", ")))
}

# Calculate missing rate if not already present
if (!"missing_rate" %in% names(passed_pairs)) {
  cat("  Calculating missing rates...\n")
  
  missing_rate <- function(r0_row, r1_row, dead_zone) {
    r0_invalid <- sum(abs(r0_row) < dead_zone)
    r1_invalid <- sum(abs(r1_row) < dead_zone)
    total_samples <- length(r0_row) + length(r1_row)
    return((r0_invalid + r1_invalid) / total_samples)
  }
  
  missing_rates <- numeric(nrow(passed_pairs))
  for (i in 1:nrow(passed_pairs)) {
    missing_rates[i] <- missing_rate(reo_r0_passed[i, ], reo_r1_passed[i, ], CONFIG$dead_zone)
  }
  
  passed_pairs$missing_rate <- missing_rates
} else {
  cat("  Missing rate already calculated in Step 3b\n")
}

cat(sprintf("  Using a[2] from Step 3b: range %.2f - %.2f\n", 
            min(passed_pairs$a2_value, na.rm = TRUE), 
            max(passed_pairs$a2_value, na.rm = TRUE)))
cat(sprintf("  Missing rate range: %.1f%% - %.1f%%\n", 
            min(passed_pairs$missing_rate) * 100, 
            max(passed_pairs$missing_rate) * 100))

# --------------------------------------------------------------------------
# 4.3 Step 4.1: Maximum bipartite matching (proper implementation)
# --------------------------------------------------------------------------

cat("\n--- 4.3 Maximum bipartite matching (1 gene = 1 pair) ---\n")

# Fix 赤-2: Proper bipartite matching implementation with simple indexing

# Calculate lexicographic scores for prioritization
passed_pairs$lex_rank1 <- rank(-passed_pairs$a2_value, na.last = TRUE)  
passed_pairs$lex_rank2 <- rank(-passed_pairs$r1_reversal_rate, na.last = TRUE)  
passed_pairs$lex_rank3 <- rank(passed_pairs$missing_rate, na.last = TRUE)  
passed_pairs$lex_rank4 <- rank(paste0(passed_pairs$gene_up, "_", passed_pairs$gene_down))

# Combine into single score (lower is better)
n_pairs <- nrow(passed_pairs)
passed_pairs$lex_score <- (
  passed_pairs$lex_rank1 * (n_pairs^3) +
    passed_pairs$lex_rank2 * (n_pairs^2) +
    passed_pairs$lex_rank3 * n_pairs +
    passed_pairs$lex_rank4
)

# Get unique genes
up_genes <- sort(unique(passed_pairs$gene_up))
down_genes <- sort(unique(passed_pairs$gene_down))
n_up <- length(up_genes)
n_down <- length(down_genes)

cat(sprintf("  Unique up genes: %d\n", n_up))
cat(sprintf("  Unique down genes: %d\n", n_down))

# Create adjacency matrix for bipartite graph
# Rows = up genes, Columns = down genes
# Values = best pair index for that combination (or 0 if no pair)
adj_matrix <- matrix(0, nrow = n_up, ncol = n_down)
weight_matrix <- matrix(Inf, nrow = n_up, ncol = n_down)

# Fill adjacency matrix with best pairs for each gene combination
for (i in 1:nrow(passed_pairs)) {
  up_idx <- which(up_genes == passed_pairs$gene_up[i])
  down_idx <- which(down_genes == passed_pairs$gene_down[i])
  
  # Store the best (lowest score) pair for this gene combination
  if (passed_pairs$lex_score[i] < weight_matrix[up_idx, down_idx]) {
    adj_matrix[up_idx, down_idx] <- i  # Store pair index
    weight_matrix[up_idx, down_idx] <- passed_pairs$lex_score[i]
  }
}

# Count actual edges
n_edges <- sum(adj_matrix > 0)
cat(sprintf("  Unique gene combinations: %d\n", n_edges))

# Create bipartite graph from adjacency matrix
# Use simple vertex numbering: 1:n_up for up genes, (n_up+1):(n_up+n_down) for down genes
g <- graph_from_incidence_matrix(adj_matrix > 0, directed = FALSE, weighted = NULL)

# Verify bipartite structure
if (!is_bipartite(g)) {
  stop("ERROR: Created graph is not bipartite. This should not happen.")
}

# Add weights (use negative of lex_score since igraph maximizes)
edge_weights <- c()
edges_mat <- as_edgelist(g)
for (i in 1:nrow(edges_mat)) {
  up_v <- edges_mat[i, 1]
  down_v <- edges_mat[i, 2] - n_up  # Adjust for indexing
  if (up_v <= n_up && down_v > 0 && down_v <= n_down) {
    edge_weights <- c(edge_weights, -weight_matrix[up_v, down_v])
  }
}
E(g)$weight <- edge_weights

# Debug output (if VERBOSE)
if (!is.null(CONFIG$VERBOSE) && CONFIG$VERBOSE) {
  cat(sprintf("  DEBUG: Graph has %d vertices and %d edges\n", vcount(g), ecount(g)))
  cat(sprintf("  DEBUG: Graph is bipartite: %s\n", is_bipartite(g)))
  
  components <- components(g)
  cat(sprintf("  DEBUG: Number of connected components: %d\n", components$no))
  
  # Test with simple case
  cat("  DEBUG: Testing max_bipartite_match with simple case...\n")
  test_g <- make_bipartite_graph(
    types = c(FALSE, FALSE, TRUE, TRUE),
    edges = c(1,3, 1,4, 2,3, 2,4),
    directed = FALSE
  )
  test_match <- max_bipartite_match(test_g)
  cat(sprintf("    Simple test matching size: %d (expected: 2)\n", test_match$matching_size))
}

# Perform maximum bipartite matching
cat("  Running maximum bipartite matching...\n")
matching <- tryCatch({
  max_bipartite_match(g, weights = E(g)$weight)
}, error = function(e) {
  cat(sprintf("  ERROR in weighted matching: %s\n", e$message))
  cat("  Attempting unweighted matching as fallback...\n")
  max_bipartite_match(g, weights = NULL)
})

cat(sprintf("  Matching found: %d pairs\n", matching$matching_size))

# Extract matched pairs
matched_pair_indices <- c()

if (matching$matching_size > 0) {
  # Process matching results
  for (up_v in 1:n_up) {
    matched_to <- matching$matching[up_v]
    
    if (!is.na(matched_to) && matched_to > n_up) {
      down_v <- matched_to - n_up
      
      # Get the pair index from adjacency matrix
      pair_idx <- adj_matrix[up_v, down_v]
      
      if (pair_idx > 0) {
        matched_pair_indices <- c(matched_pair_indices, pair_idx)
      }
    }
  }
}

# Extract selected pairs
if (length(matched_pair_indices) > 0) {
  gene_unique_pairs <- passed_pairs[matched_pair_indices, ]
  gene_unique_reo_r0 <- reo_r0_passed[matched_pair_indices, ]
  gene_unique_reo_r1 <- reo_r1_passed[matched_pair_indices, ]
} else {
  # Fallback to greedy if matching completely fails
  cat("  WARNING: Bipartite matching returned no results. Using greedy fallback...\n")
  
  # Sort by lex_score and greedily select
  sorted_indices <- order(passed_pairs$lex_score)
  sorted_pairs <- passed_pairs[sorted_indices, ]
  
  used_up <- character(0)
  used_down <- character(0)
  selected <- c()
  
  for (i in 1:nrow(sorted_pairs)) {
    if (!sorted_pairs$gene_up[i] %in% used_up && 
        !sorted_pairs$gene_down[i] %in% used_down) {
      selected <- c(selected, sorted_indices[i])
      used_up <- c(used_up, sorted_pairs$gene_up[i])
      used_down <- c(used_down, sorted_pairs$gene_down[i])
    }
  }
  
  gene_unique_pairs <- passed_pairs[selected, ]
  gene_unique_reo_r0 <- reo_r0_passed[selected, ]
  gene_unique_reo_r1 <- reo_r1_passed[selected, ]
}

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
  cat("  Only one pair available. No clustering needed.\n")
  final_pairs <- gene_unique_pairs
  final_reo_r0 <- gene_unique_reo_r0
  final_reo_r1 <- gene_unique_reo_r1
} else {
  # Calculate pairwise correlations
  combined_reo <- cbind(gene_unique_reo_r0, gene_unique_reo_r1)
  n_pairs_unique <- nrow(gene_unique_pairs)
  
  # Initialize correlation matrix
  cor_matrix <- matrix(0, nrow = n_pairs_unique, ncol = n_pairs_unique)
  diag(cor_matrix) <- 1
  
  # Calculate pairwise correlations
  for (i in 1:(n_pairs_unique - 1)) {
    for (j in (i + 1):n_pairs_unique) {
      cor_value <- cor(combined_reo[i, ], combined_reo[j, ], 
                       method = cor_method, use = "complete.obs")
      if (!is.na(cor_value)) {
        cor_matrix[i, j] <- cor_value
        cor_matrix[j, i] <- cor_value
      }
    }
  }
  
  # Find pairs with high correlation
  high_cor_pairs <- which(abs(cor_matrix) > cor_threshold & 
                            upper.tri(cor_matrix), arr.ind = TRUE)
  
  cat(sprintf("  Pairs with |ρ| > %.2f: %d\n", cor_threshold, nrow(high_cor_pairs)))
  
  if (nrow(high_cor_pairs) == 0) {
    # No high correlations found
    cat("  No highly correlated pairs found. Keeping all gene-unique pairs.\n")
    final_pairs <- gene_unique_pairs
    final_reo_r0 <- gene_unique_reo_r0
    final_reo_r1 <- gene_unique_reo_r1
  } else {
    # Create graph for correlation clustering
    cor_g <- make_empty_graph(n = n_pairs_unique, directed = FALSE)
    
    # Add edges for high correlation pairs
    edges_to_add <- as.vector(t(high_cor_pairs))
    cor_g <- add_edges(cor_g, edges_to_add)
    
    # Find connected components (correlation clusters)
    clusters <- components(cor_g)
    n_clusters <- clusters$no
    
    cat(sprintf("  Number of correlation clusters: %d\n", n_clusters))
    cat(sprintf("  Singleton clusters: %d\n", sum(clusters$csize == 1)))
    cat(sprintf("  Multi-member clusters: %d\n", sum(clusters$csize > 1)))
    
    # Select representative from each cluster
    cluster_representatives <- integer(n_clusters)
    
    for (cluster_id in 1:n_clusters) {
      cluster_members <- which(clusters$membership == cluster_id)
      
      if (length(cluster_members) == 1) {
        cluster_representatives[cluster_id] <- cluster_members[1]
      } else {
        # Select best member by lex_score
        cluster_subset <- gene_unique_pairs[cluster_members, ]
        best_idx <- cluster_members[which.min(cluster_subset$lex_score)]
        cluster_representatives[cluster_id] <- best_idx
      }
    }
    
    # Extract final deduplicated pairs
    final_pairs <- gene_unique_pairs[cluster_representatives, ]
    final_reo_r0 <- gene_unique_reo_r0[cluster_representatives, ]
    final_reo_r1 <- gene_unique_reo_r1[cluster_representatives, ]
    
    # Add cluster size information
    final_pairs$cluster_size <- clusters$csize[clusters$membership[cluster_representatives]]
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
  final_pairs$dedup_method <- "max_bipartite_matching"
  
  # Add cluster information if not already present
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
    correlation_retention = if(nrow(gene_unique_pairs) > 0) {
      nrow(final_pairs) / nrow(gene_unique_pairs)
    } else {
      NA
    },
    correlation_threshold_used = cor_threshold,
    correlation_method_used = cor_method,
    matching_type = ifelse(exists("matching") && matching$matching_size > 0, "bipartite", "greedy_fallback")
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
cat(sprintf("Correlation parameters: %s (threshold=%.2f)\n", cor_method, cor_threshold))
cat(sprintf("Matching method: %s\n", step4_data$dedup_stats$matching_type))
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