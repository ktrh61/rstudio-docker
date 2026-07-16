# ==============================================================================
# REBC-THYR Feature Selection Script v6 - Option C with Zero Filtering
# 11_feature_selection_v6_optionC_cleaned.R
# ==============================================================================

library(data.table)

cat("Starting feature selection v6 Option C with zero filtering...\n")
cat("Strategy: R0 consistency + R0/R1 separation, excluding sparse genes\n\n")

# ==============================================================================
# 1. Load Data
# ==============================================================================

cat("Loading data...\n")

# Load TPM v6
load("./data/processed/tpm_v6.rda")
tpm <- tpm_v6

# Load DEG results v6
load("./data/processed/final_deg_v6_results.rda")

# Load high-purity sample lists
load("./data/processed/final_high_purity_sample_lists_v6.rda")

# Get R0 and R1 tumor samples
r0_samples <- final_high_purity_sample_lists_v6[["R0"]]$tumor
r1_samples <- final_high_purity_sample_lists_v6[["R1"]]$tumor

# Extract TPM for these samples
tpm_r0 <- tpm[, r0_samples]
tpm_r1 <- tpm[, r1_samples]

cat(sprintf("Samples: R0=%d, R1=%d\n", ncol(tpm_r0), ncol(tpm_r1)))

# ==============================================================================
# 2. Extract and Clean DEG Lists
# ==============================================================================

cat("\nExtracting DEG lists...\n")

# Extract R0_vs_R1_tumor DEGs
deg_result <- final_deg_v6_results$deg_analysis_results[["R0_vs_R1_tumor"]]
degs <- deg_result$deg_summary$results_df[deg_result$deg_summary$results_df$significant, ]

# Get CORRECT up and down gene lists
up_genes_all <- degs$gene_id[degs$direction == "UP"]
down_genes_all <- degs$gene_id[degs$direction == "DOWN"]

cat(sprintf("Initial DEGs: %d UP, %d DOWN\n", length(up_genes_all), length(down_genes_all)))

# Ensure genes exist in TPM data
up_genes_all <- up_genes_all[up_genes_all %in% rownames(tpm)]
down_genes_all <- down_genes_all[down_genes_all %in% rownames(tpm)]

cat(sprintf("Available in TPM: %d UP, %d DOWN\n", length(up_genes_all), length(down_genes_all)))

# ==============================================================================
# 3. Filter Genes with Excessive Zero Values
# ==============================================================================

cat("\nFiltering genes with excessive zero values...\n")

# Function to check zero proportion
check_zeros <- function(tpm_matrix, gene_list, max_zero_prop = 0.3) {
  if (length(gene_list) == 0) return(character(0))
  
  gene_list_valid <- gene_list[gene_list %in% rownames(tpm_matrix)]
  if (length(gene_list_valid) == 0) return(character(0))
  
  zero_props <- apply(tpm_matrix[gene_list_valid, , drop = FALSE], 1, 
                      function(x) sum(x == 0) / length(x))
  gene_list_valid[zero_props <= max_zero_prop]
}

# Filter UP genes (max 30% zeros in both R0 and R1)
up_genes_r0_ok <- check_zeros(tpm_r0, up_genes_all, 0.3)
up_genes_r1_ok <- check_zeros(tpm_r1, up_genes_all, 0.3)
up_genes_filtered <- intersect(up_genes_r0_ok, up_genes_r1_ok)

# Filter DOWN genes (max 30% zeros in both R0 and R1)
down_genes_r0_ok <- check_zeros(tpm_r0, down_genes_all, 0.3)
down_genes_r1_ok <- check_zeros(tpm_r1, down_genes_all, 0.3)
down_genes_filtered <- intersect(down_genes_r0_ok, down_genes_r1_ok)

cat(sprintf("After zero filtering (≤30%% zeros): %d UP, %d DOWN\n", 
            length(up_genes_filtered), length(down_genes_filtered)))

# ==============================================================================
# 4. Additional Expression Filter
# ==============================================================================

cat("\nApplying minimum expression filter...\n")

# Require median TPM ≥ 5 in at least one group
min_median_tpm <- 5

# UP genes
if (length(up_genes_filtered) > 0) {
  median_r0_up <- apply(tpm_r0[up_genes_filtered, , drop = FALSE], 1, median)
  median_r1_up <- apply(tpm_r1[up_genes_filtered, , drop = FALSE], 1, median)
  up_genes_expressed <- up_genes_filtered[(median_r0_up >= min_median_tpm) | 
                                            (median_r1_up >= min_median_tpm)]
} else {
  up_genes_expressed <- character(0)
}

# DOWN genes
if (length(down_genes_filtered) > 0) {
  median_r0_down <- apply(tpm_r0[down_genes_filtered, , drop = FALSE], 1, median)
  median_r1_down <- apply(tpm_r1[down_genes_filtered, , drop = FALSE], 1, median)
  down_genes_expressed <- down_genes_filtered[(median_r0_down >= min_median_tpm) | 
                                                (median_r1_down >= min_median_tpm)]
} else {
  down_genes_expressed <- character(0)
}

cat(sprintf("After expression filter (median TPM ≥ %d): %d UP, %d DOWN\n", 
            min_median_tpm, length(up_genes_expressed), length(down_genes_expressed)))
cat(sprintf("Total pairs to evaluate: %d\n", 
            length(up_genes_expressed) * length(down_genes_expressed)))

# ==============================================================================
# 5. Option C Analysis: R0 Consistency + Separation
# ==============================================================================

if (length(up_genes_expressed) == 0 || length(down_genes_expressed) == 0) {
  cat("\nNo genes remaining after filtering. Consider relaxing criteria.\n")
  stop("Insufficient genes for analysis")
}

cat("\nEvaluating gene pairs with Option C criteria...\n")

# Parameters
min_separation <- 2.0
consistency_threshold <- 0.9  # 90% of R0 samples

# Extract matrices
up_tpm_r0_mat <- as.matrix(tpm_r0[up_genes_expressed, ])
up_tpm_r1_mat <- as.matrix(tpm_r1[up_genes_expressed, ])
down_tpm_r0_mat <- as.matrix(tpm_r0[down_genes_expressed, ])
down_tpm_r1_mat <- as.matrix(tpm_r1[down_genes_expressed, ])

pc <- 0.1  # Pseudocount
n_up <- nrow(up_tpm_r0_mat)
n_down <- nrow(down_tpm_r0_mat)

# Process in chunks
results_list <- list()
chunk_size <- 50
n_chunks <- ceiling(n_up / chunk_size)

cat(sprintf("Processing %d chunks...\n", n_chunks))

for (chunk in 1:n_chunks) {
  if (chunk %% 10 == 1 || chunk == n_chunks) {
    cat(sprintf("  Chunk %d/%d (%.1f%%)...\n", chunk, n_chunks, chunk/n_chunks*100))
  }
  
  # Define chunk range
  start_idx <- (chunk - 1) * chunk_size + 1
  end_idx <- min(chunk * chunk_size, n_up)
  chunk_indices <- start_idx:end_idx
  
  # Get chunk data
  up_r0_chunk <- up_tpm_r0_mat[chunk_indices, , drop = FALSE]
  up_r1_chunk <- up_tpm_r1_mat[chunk_indices, , drop = FALSE]
  
  for (i in 1:length(chunk_indices)) {
    up_idx <- chunk_indices[i]
    
    # Calculate ratios for all down genes
    ratios_r0 <- sweep(down_tpm_r0_mat + pc, 2, up_r0_chunk[i, ] + pc, "/")
    ratios_r0 <- log2(1/ratios_r0)  # Flip to get UP/DOWN
    
    ratios_r1 <- sweep(down_tpm_r1_mat + pc, 2, up_r1_chunk[i, ] + pc, "/")
    ratios_r1 <- log2(1/ratios_r1)  # Flip to get UP/DOWN
    
    # Check Filter 1: R0 consistency
    n_positive_r0 <- rowSums(ratios_r0 > 0)
    n_negative_r0 <- rowSums(ratios_r0 < 0)
    n_r0_consistent <- pmax(n_positive_r0, n_negative_r0)
    pass_consistency <- n_r0_consistent >= (consistency_threshold * ncol(ratios_r0))
    
    # Check Filter 2: Separation
    median_r0 <- apply(ratios_r0, 1, median)
    median_r1 <- apply(ratios_r1, 1, median)
    separation <- abs(median_r1 - median_r0)
    pass_separation <- separation >= min_separation
    
    # Find passing pairs
    pass_both <- pass_consistency & pass_separation
    
    if (any(pass_both)) {
      passing_idx <- which(pass_both)
      
      for (j in passing_idx) {
        pair_key <- paste(up_genes_expressed[up_idx], 
                          down_genes_expressed[j], sep = "_")
        
        results_list[[pair_key]] <- data.frame(
          up_gene = up_genes_expressed[up_idx],
          down_gene = down_genes_expressed[j],
          r0_median = median_r0[j],
          r1_median = median_r1[j],
          separation = separation[j],
          r0_consistency = n_r0_consistent[j] / ncol(ratios_r0),
          r0_direction = ifelse(n_positive_r0[j] > n_negative_r0[j], "positive", "negative"),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

# ==============================================================================
# 6. Compile and Save Results
# ==============================================================================

cat("\nCompiling results...\n")

if (length(results_list) > 0) {
  results_df <- do.call(rbind, results_list)
  results_df <- results_df[order(results_df$separation, decreasing = TRUE), ]
  rownames(results_df) <- NULL
  
  cat(sprintf("Found %d gene pairs passing Option C criteria\n", nrow(results_df)))
  
  # Show top results
  cat("\nTop 20 pairs by separation:\n")
  print(head(results_df, 20))
  
  # Filter for high-confidence pairs
  high_conf <- results_df[results_df$separation >= 3.0, ]
  cat(sprintf("\nHigh-confidence pairs (separation ≥ 3.0): %d\n", nrow(high_conf)))
  
  if (nrow(high_conf) > 0 && nrow(high_conf) <= 50) {
    cat("\nAll high-confidence pairs:\n")
    print(high_conf[, c("up_gene", "down_gene", "r0_median", "r1_median", "separation")])
  }
  
  # Save results
  write.csv(results_df, "./output/reports/gene_pairs_optionC_cleaned_v6.csv", row.names = FALSE)
  save(results_df, file = "./data/processed/gene_pairs_optionC_cleaned_v6.rda")
  
  cat("\nResults saved to:\n")
  cat("  ./output/reports/gene_pairs_optionC_cleaned_v6.csv\n")
  cat("  ./data/processed/gene_pairs_optionC_cleaned_v6.rda\n")
  
} else {
  cat("\nNo pairs found meeting criteria.\n")
  cat("Consider:\n")
  cat("  - Reducing separation threshold (currently 2.0)\n")
  cat("  - Allowing more zeros (currently 30%)\n")
  cat("  - Lowering minimum expression (currently median TPM ≥ 5)\n")
}

cat("\nFeature selection v6 Option C completed!\n")

