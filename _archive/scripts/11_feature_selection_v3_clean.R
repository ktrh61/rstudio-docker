# ==============================================================================
# REBC-THYR Feature Selection Script v3 - Clean Version
# 11_feature_selection_v3_clean.R
# ==============================================================================

# Required libraries
library(dplyr)
library(pROC)

cat("Starting feature selection v3: Clean implementation...\n")

# ==============================================================================
# 1. Load Data and DEG Results
# ==============================================================================

cat("Loading DEG v3 results and TPM data...\n")

# Load DEG analysis v3 results
load("./data/processed/final_deg_v3_results.rda")

# Load TPM data
load("./data/processed/thyr_tpm.rda")

# Load high-purity sample lists for R0/R1 groups
load("./data/processed/final_high_purity_sample_lists.rda")
high_purity_sample_lists <- final_high_purity_sample_lists

cat("Data loaded successfully\n")

# ==============================================================================
# 2. Extract DEG Results and Prepare Data
# ==============================================================================

cat("Extracting R0_vs_R1_tumor DEG results...\n")

# Check if R0_vs_R1_tumor comparison is available
if (!"R0_vs_R1_tumor" %in% names(final_deg_v3_results$deg_analysis_results)) {
  stop("R0_vs_R1_tumor comparison not found in DEG results")
}

# Get tumor DEG results
tumor_deg_result <- final_deg_v3_results$deg_analysis_results[["R0_vs_R1_tumor"]]
tumor_degs <- tumor_deg_result$deg_summary$results_df

# Extract significant UP and DOWN genes
up_genes <- tumor_degs$gene_id[tumor_degs$significance_category == "Upregulated"]
down_genes <- tumor_degs$gene_id[tumor_degs$significance_category == "Downregulated"]

cat(sprintf("DEG genes: %d UP, %d DOWN\n", length(up_genes), length(down_genes)))

# ==============================================================================
# 3. Prepare Sample Data for Analysis
# ==============================================================================

cat("Preparing R0/R1 sample data...\n")

# Get R0 and R1 tumor samples
r0_tumor_samples <- high_purity_sample_lists$R0$tumor
r1_tumor_samples <- high_purity_sample_lists$R1$tumor

# Combine samples and create labels
all_tumor_samples <- c(r0_tumor_samples, r1_tumor_samples)
sample_labels <- c(rep("R0", length(r0_tumor_samples)), rep("R1", length(r1_tumor_samples)))

# Check sample availability in TPM data
available_samples <- all_tumor_samples[all_tumor_samples %in% colnames(tpm)]
available_labels <- sample_labels[all_tumor_samples %in% colnames(tpm)]

# Extract TPM data for analysis
analysis_tpm <- tpm[, available_samples]

# Create binary outcome (R0=0, R1=1)
outcome <- ifelse(available_labels == "R0", 0, 1)
n_r0 <- sum(outcome == 0)
n_r1 <- sum(outcome == 1)

cat(sprintf("Analysis samples: %d total (R0=%d, R1=%d)\n", length(available_samples), n_r0, n_r1))

# Filter genes to those available in TPM data
up_genes_avail <- up_genes[up_genes %in% rownames(analysis_tpm)]
down_genes_avail <- down_genes[down_genes %in% rownames(analysis_tpm)]

cat(sprintf("Available genes: %d UP, %d DOWN\n", length(up_genes_avail), length(down_genes_avail)))

total_pairs <- length(up_genes_avail) * length(down_genes_avail)
cat(sprintf("Total pairs to evaluate: %d\n", total_pairs))

# ==============================================================================
# 4. Set Filter Criteria
# ==============================================================================

# Filter settings (based on diagnostic results)
r0_specificity_threshold <- 0.95  # 95% (high specificity required - near 100%)
separation_threshold <- 1.0       # 1.0 log2 units
r0_consistency_threshold <- 2.0   # 2.0 log2 units

cat(sprintf("Filter criteria:\n"))
cat(sprintf("  R0 specificity >= %.1f%%\n", r0_specificity_threshold * 100))
cat(sprintf("  Separation >= %.1f log2 units\n", separation_threshold))
cat(sprintf("  R0 range <= %.1f log2 units\n", r0_consistency_threshold))

# ==============================================================================
# 5. Feature Selection with Simple Loop
# ==============================================================================

cat("Starting feature selection...\n")

# Initialize storage for candidates
candidate_pairs <- data.frame(
  up_gene = character(0),
  down_gene = character(0),
  r0_mean = numeric(0),
  r1_mean = numeric(0),
  separation = numeric(0),
  r0_range = numeric(0),
  r0_specificity = numeric(0),
  stringsAsFactors = FALSE
)

# Progress tracking
progress_interval <- max(1, floor(length(up_genes_avail) / 20))
candidate_count <- 0
pseudocount <- 0.1

start_time <- Sys.time()

for (i in seq_along(up_genes_avail)) {
  
  up_gene <- up_genes_avail[i]
  up_tpm_values <- as.numeric(analysis_tpm[up_gene, ])
  
  for (j in seq_along(down_genes_avail)) {
    
    down_gene <- down_genes_avail[j]
    down_tpm_values <- as.numeric(analysis_tpm[down_gene, ])
    
    # Calculate log2 ratio
    ratio_values <- log2((up_tpm_values + pseudocount) / (down_tpm_values + pseudocount))
    
    # Split by outcome
    r0_ratios <- ratio_values[outcome == 0]
    r1_ratios <- ratio_values[outcome == 1]
    
    # Calculate statistics
    r0_mean <- mean(r0_ratios)
    r1_mean <- mean(r1_ratios)
    separation <- abs(r1_mean - r0_mean)
    r0_range <- max(r0_ratios) - min(r0_ratios)
    
    # Calculate R0 specificity
    threshold_est <- (r0_mean + r1_mean) / 2
    if (r0_mean < r1_mean) {
      r0_correct <- sum(r0_ratios <= threshold_est)
    } else {
      r0_correct <- sum(r0_ratios >= threshold_est)
    }
    specificity <- r0_correct / length(r0_ratios)
    
    # Apply filters
    filter1_pass <- specificity >= r0_specificity_threshold
    filter2_pass <- separation > 0  # Direction clear (always true)
    filter3_pass <- separation >= separation_threshold
    filter4_pass <- r0_range <= r0_consistency_threshold
    
    all_filters_pass <- filter1_pass & filter2_pass & filter3_pass & filter4_pass
    
    # Store candidates that pass all filters
    if (all_filters_pass) {
      new_candidate <- data.frame(
        up_gene = up_gene,
        down_gene = down_gene,
        r0_mean = r0_mean,
        r1_mean = r1_mean,
        separation = separation,
        r0_range = r0_range,
        r0_specificity = specificity,
        stringsAsFactors = FALSE
      )
      
      candidate_pairs <- rbind(candidate_pairs, new_candidate)
      candidate_count <- candidate_count + 1
    }
  }
  
  # Progress tracking
  if (i %% progress_interval == 0) {
    elapsed <- as.numeric(Sys.time() - start_time)
    progress_pct <- i / length(up_genes_avail) * 100
    pairs_processed <- i * length(down_genes_avail)
    cat(sprintf("  Progress: %.1f%% (%d UP genes, %d candidates, %.1fs)\n", 
                progress_pct, i, candidate_count, elapsed))
  }
}

end_time <- Sys.time()
total_time <- as.numeric(end_time - start_time)

cat(sprintf("Feature selection completed in %.1f seconds\n", total_time))
cat(sprintf("Candidates found: %d/%d (%.4f%%)\n", 
            nrow(candidate_pairs), total_pairs, nrow(candidate_pairs)/total_pairs*100))

# ==============================================================================
# 6. ROC Analysis and Final Selection
# ==============================================================================

if (nrow(candidate_pairs) > 0) {
  cat("Performing ROC analysis...\n")
  
  # Function for ROC analysis
  calculate_roc_metrics <- function(i) {
    up_gene <- candidate_pairs$up_gene[i]
    down_gene <- candidate_pairs$down_gene[i]
    
    # Extract TPM values
    up_tpm <- as.numeric(analysis_tpm[up_gene, ])
    down_tpm <- as.numeric(analysis_tpm[down_gene, ])
    
    # Calculate ratio
    ratio <- log2((up_tpm + pseudocount) / (down_tpm + pseudocount))
    
    # ROC analysis
    tryCatch({
      roc_result <- roc(outcome, ratio, quiet = TRUE)
      auc_value <- as.numeric(auc(roc_result))
      
      # Get optimal threshold
      coords_result <- coords(roc_result, "best", best.method = "youden")
      
      list(
        auc = auc_value,
        threshold = coords_result$threshold,
        sensitivity = coords_result$sensitivity,
        specificity = coords_result$specificity
      )
    }, error = function(e) {
      list(auc = NA, threshold = NA, sensitivity = NA, specificity = NA)
    })
  }
  
  # Apply ROC analysis to all candidates
  cat("Computing ROC metrics for all candidates...\n")
  roc_results <- lapply(1:nrow(candidate_pairs), calculate_roc_metrics)
  
  # Extract metrics
  candidate_pairs$auc <- sapply(roc_results, function(x) x$auc)
  candidate_pairs$threshold <- sapply(roc_results, function(x) x$threshold)
  candidate_pairs$sensitivity <- sapply(roc_results, function(x) x$sensitivity)
  candidate_pairs$specificity <- sapply(roc_results, function(x) x$specificity)
  
  # Remove failed ROC analyses
  valid_candidates <- candidate_pairs[!is.na(candidate_pairs$auc), ]
  
  cat(sprintf("Valid ROC results: %d/%d candidates\n", nrow(valid_candidates), nrow(candidate_pairs)))
  
  if (nrow(valid_candidates) > 0) {
    # Sort by AUC (descending)
    valid_candidates <- valid_candidates[order(valid_candidates$auc, decreasing = TRUE), ]
    
    # Select top candidates (limit to 2 pairs to avoid overfitting)
    final_pairs <- head(valid_candidates, 2)
    
    cat("\n=== FINAL SELECTED PAIRS ===\n")
    for (i in 1:nrow(final_pairs)) {
      pair <- final_pairs[i, ]
      cat(sprintf("%d. %s / %s\n", i, pair$up_gene, pair$down_gene))
      cat(sprintf("   AUC=%.3f, R0_Spec=%.1f%%, Separation=%.2f, R0_Range=%.2f\n",
                  pair$auc, pair$r0_specificity * 100, pair$separation, pair$r0_range))
    }
    
  } else {
    cat("No valid ROC results obtained\n")
    final_pairs <- data.frame()
  }
  
} else {
  cat("No candidates found meeting all filter criteria\n")
  final_pairs <- data.frame()
}

# ==============================================================================
# 7. Save Results
# ==============================================================================

cat("Saving results...\n")

# Create output directories
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}
if (!dir.exists("./output/reports")) {
  dir.create("./output/reports", recursive = TRUE)
}

# Save comprehensive results
feature_selection_v3_results <- list(
  analysis_parameters = list(
    r0_specificity_threshold = r0_specificity_threshold,
    separation_threshold = separation_threshold,
    r0_consistency_threshold = r0_consistency_threshold,
    total_pairs_evaluated = total_pairs,
    processing_time_seconds = total_time,
    analysis_version = "v3_clean_simple_loop"
  ),
  input_data = list(
    up_genes_total = length(up_genes),
    down_genes_total = length(down_genes),
    up_genes_available = length(up_genes_avail),
    down_genes_available = length(down_genes_avail),
    r0_samples = n_r0,
    r1_samples = n_r1
  ),
  all_candidates = candidate_pairs,
  final_selected_pairs = if(exists("final_pairs")) final_pairs else data.frame(),
  performance_summary = list(
    candidates_found = nrow(candidate_pairs),
    candidate_rate = nrow(candidate_pairs) / total_pairs,
    final_pairs_count = if(exists("final_pairs")) nrow(final_pairs) else 0,
    processing_efficiency = paste0(round(total_pairs / total_time), " pairs/second")
  ),
  analysis_date = Sys.time()
)

save(feature_selection_v3_results, file = "./data/processed/feature_selection_v3_results.rda")

# Save CSV files
if (nrow(candidate_pairs) > 0) {
  write.csv(candidate_pairs, "./output/reports/all_candidates_v3.csv", row.names = FALSE)
}

if (exists("final_pairs") && nrow(final_pairs) > 0) {
  write.csv(final_pairs, "./output/reports/final_selected_pairs_v3.csv", row.names = FALSE)
}

# ==============================================================================
# 8. Final Summary
# ==============================================================================

cat("\n==============================================\n")
cat("Feature Selection v3 Clean - Final Summary\n")
cat("==============================================\n")

cat("Simple Loop Implementation Completed!\n\n")

cat("Performance metrics:\n")
cat(sprintf("  Total pairs evaluated: %d\n", total_pairs))
cat(sprintf("  Processing time: %.1f seconds\n", total_time))
cat(sprintf("  Processing speed: %d pairs/second\n", round(total_pairs / total_time)))
cat(sprintf("  Candidates found: %d (%.4f%%)\n", nrow(candidate_pairs), nrow(candidate_pairs)/total_pairs*100))

if (exists("final_pairs") && nrow(final_pairs) > 0) {
  cat("\n✅ SUCCESS: Final pairs selected\n")
  cat(sprintf("  Final pairs: %d\n", nrow(final_pairs)))
} else {
  cat("\n❌ No suitable pairs found\n")
  cat("Consider relaxing filter criteria if needed\n")
}

cat("\nNext steps:\n")
cat("1. Review selected pairs for biological relevance\n")
cat("2. Proceed to validation experiment (12_validation_experiment.R)\n")
cat("3. Test performance across ERR gradient\n")

cat("\nFeature selection v3 clean completed!\n")
cat("==============================================\n")