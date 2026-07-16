# ==============================================================================
# REBC-THYR Feature Selection Script v6 - High-Speed 3-Filter Strategy
# 11_feature_selection_v6.R
# ==============================================================================

library(data.table)
library(parallel)

cat("Starting feature selection v6 with 3-filter strategy...\n")
cat("Optimized for speed using data.table and parallel processing\n")

# ==============================================================================
# 1. Load Data
# ==============================================================================

cat("\nLoading data...\n")

# Load TPM v6
load("./data/processed/tpm_v6.rda")
tpm <- tpm_v6

# Load DEG results v6
load("./data/processed/final_deg_v6_results.rda")

# Load high-purity sample lists
load("./data/processed/final_high_purity_sample_lists_v6.rda")

# Extract R0_vs_R1_tumor DEGs
deg_result <- final_deg_v6_results$deg_analysis_results[["R0_vs_R1_tumor"]]
degs <- deg_result$deg_summary$results_df[deg_result$deg_summary$results_df$significant, ]

up_genes <- degs$gene_id[degs$direction == "UP"]
down_genes <- degs$gene_id[degs$direction == "DOWN"]

cat(sprintf("DEGs loaded: %d UP, %d DOWN\n", length(up_genes), length(down_genes)))

# Get R0 and R1 tumor samples
r0_samples <- final_high_purity_sample_lists_v6[["R0"]]$tumor
r1_samples <- final_high_purity_sample_lists_v6[["R1"]]$tumor

# Extract TPM for these samples
tpm_r0 <- tpm[, r0_samples]
tpm_r1 <- tpm[, r1_samples]

cat(sprintf("Samples: R0=%d, R1=%d\n", ncol(tpm_r0), ncol(tpm_r1)))

# ==============================================================================
# 2. Pre-filter Genes (Simplified)
# ==============================================================================

cat("\nChecking gene availability in TPM data...\n")

# Simply check if DEGs exist in TPM matrix
up_genes_available <- up_genes[up_genes %in% rownames(tpm)]
down_genes_available <- down_genes[down_genes %in% rownames(tpm)]

cat(sprintf("Available in TPM: %d/%d UP, %d/%d DOWN\n", 
            length(up_genes_available), length(up_genes),
            length(down_genes_available), length(down_genes)))

# Use all available DEGs without additional filtering
up_genes_filtered <- up_genes_available
down_genes_filtered <- down_genes_available

# ==============================================================================
# 3. Vectorized Filter Functions
# ==============================================================================

# Pre-calculate all ratios (vectorized)
calculate_all_ratios <- function(up_tpm, down_tpm) {
  # Add pseudocount
  pc <- 0.1
  log2((up_tpm + pc) / (down_tpm + pc))
}

# Filter 1: R0 consistency (90% same direction)
check_filter1 <- function(ratios_r0) {
  n_r0 <- ncol(ratios_r0)
  threshold <- ceiling(0.9 * n_r0)
  
  # Vectorized: count positive and negative per row
  n_positive <- rowSums(ratios_r0 > 0)
  n_negative <- rowSums(ratios_r0 < 0)
  
  (n_positive >= threshold) | (n_negative >= threshold)
}

# Filter 2: Direction reversal between R0 and R1 medians
check_filter2 <- function(ratios_r0, ratios_r1) {
  median_r0 <- apply(ratios_r0, 1, median)
  median_r1 <- apply(ratios_r1, 1, median)
  
  (median_r0 > 0 & median_r1 < 0) | (median_r0 < 0 & median_r1 > 0)
}

# Filter 3: Minimum absolute ratio
check_filter3 <- function(ratios_r0, ratios_r1, min_ratio = 1.5) {
  # Check if all samples meet minimum absolute ratio
  all_r0_pass <- rowSums(abs(ratios_r0) >= min_ratio) == ncol(ratios_r0)
  all_r1_pass <- rowSums(abs(ratios_r1) >= min_ratio) == ncol(ratios_r1)
  
  all_r0_pass & all_r1_pass
}

# ==============================================================================
# 4. High-Speed Pair Evaluation
# ==============================================================================

cat("\nEvaluating gene pairs with 3-filter strategy...\n")
cat(sprintf("Total pairs to evaluate: %d\n", 
            length(up_genes_filtered) * length(down_genes_filtered)))

# Extract TPM matrices
up_tpm_r0 <- tpm_r0[up_genes_filtered, , drop = FALSE]
up_tpm_r1 <- tpm_r1[up_genes_filtered, , drop = FALSE]
down_tpm_r0 <- tpm_r0[down_genes_filtered, , drop = FALSE]
down_tpm_r1 <- tpm_r1[down_genes_filtered, , drop = FALSE]

# Initialize results storage
results <- list()
batch_size <- 100  # Process in batches for memory efficiency

# Process in batches
n_up <- length(up_genes_filtered)
n_batches <- ceiling(n_up / batch_size)

for (batch in 1:n_batches) {
  # Define batch indices
  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- min(batch * batch_size, n_up)
  batch_up_genes <- up_genes_filtered[start_idx:end_idx]
  
  if (batch %% 10 == 1) {
    cat(sprintf("  Processing batch %d/%d...\n", batch, n_batches))
  }
  
  # Get batch TPM data
  batch_up_r0 <- up_tpm_r0[batch_up_genes, , drop = FALSE]
  batch_up_r1 <- up_tpm_r1[batch_up_genes, , drop = FALSE]
  
  # Calculate all ratios for this batch (vectorized)
  for (j in 1:length(down_genes_filtered)) {
    down_gene <- down_genes_filtered[j]
    
    # Expand down gene TPM to match batch size
    down_r0_expanded <- matrix(rep(down_tpm_r0[down_gene, ], nrow(batch_up_r0)), 
                               nrow = nrow(batch_up_r0), byrow = TRUE)
    down_r1_expanded <- matrix(rep(down_tpm_r1[down_gene, ], nrow(batch_up_r1)), 
                               nrow = nrow(batch_up_r1), byrow = TRUE)
    
    # Calculate ratios
    ratios_r0 <- calculate_all_ratios(batch_up_r0, down_r0_expanded)
    ratios_r1 <- calculate_all_ratios(batch_up_r1, down_r1_expanded)
    
    # Apply filters (vectorized)
    pass_filter1 <- check_filter1(ratios_r0)
    pass_filter2 <- check_filter2(ratios_r0, ratios_r1)
    pass_filter3 <- check_filter3(ratios_r0, ratios_r1)
    
    # Find pairs passing all filters
    pass_all <- pass_filter1 & pass_filter2 & pass_filter3
    
    if (any(pass_all)) {
      passing_indices <- which(pass_all)
      for (i in passing_indices) {
        up_gene <- batch_up_genes[i]
        
        # Calculate detailed statistics for passing pairs
        r0_ratios <- ratios_r0[i, ]
        r1_ratios <- ratios_r1[i, ]
        
        result <- list(
          up_gene = up_gene,
          down_gene = down_gene,
          r0_median = median(r0_ratios),
          r1_median = median(r1_ratios),
          r0_range = max(r0_ratios) - min(r0_ratios),
          r1_range = max(r1_ratios) - min(r1_ratios),
          separation = abs(median(r1_ratios) - median(r0_ratios))
        )
        
        results[[paste(up_gene, down_gene, sep = "_")]] <- result
      }
    }
  }
}

# Convert results to data frame
if (length(results) > 0) {
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  results_df <- results_df[order(results_df$separation, decreasing = TRUE), ]
  
  cat(sprintf("\nFound %d gene pairs passing all filters\n", nrow(results_df)))
  
  # Show top pairs
  cat("\nTop 10 pairs by separation:\n")
  print(head(results_df, 10))
  
  # Save results
  write.csv(results_df, "./output/reports/gene_pairs_3filter_v6.csv", row.names = FALSE)
  save(results_df, file = "./data/processed/gene_pairs_3filter_v6.rda")
  
} else {
  cat("\nNo pairs found passing all 3 filters\n")
  cat("Consider relaxing filter criteria:\n")
  cat("  - Reduce consistency threshold from 90% to 80%\n")
  cat("  - Reduce minimum ratio from 1.5 to 1.0\n")
  cat("  - Allow partial direction reversal\n")
}

# ==============================================================================
# 5. Summary
# ==============================================================================

cat("\n==============================================\n")
cat("Feature Selection v6 Summary\n")
cat("==============================================\n")

if (length(results) > 0) {
  cat(sprintf("Successfully identified %d gene pairs\n", nrow(results_df)))
  cat(sprintf("Best separation: %.2f\n", max(results_df$separation)))
  
  # Export top pair for validation
  top_pair <- results_df[1, ]
  cat(sprintf("\nRecommended pair for validation:\n"))
  cat(sprintf("  UP gene: %s\n", top_pair$up_gene))
  cat(sprintf("  DOWN gene: %s\n", top_pair$down_gene))
  cat(sprintf("  R0 median ratio: %.2f\n", top_pair$r0_median))
  cat(sprintf("  R1 median ratio: %.2f\n", top_pair$r1_median))
} else {
  cat("No pairs met all criteria - parameter adjustment needed\n")
}

cat("\nFeature selection v6 completed!\n")