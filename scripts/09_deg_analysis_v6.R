# ==============================================================================
# REBC-THYR DEG Analysis Script v6 - Brunner-Munzel + Storey Method
# 09_deg_analysis_v6.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(edgeR)
library(qvalue)
library(dplyr)
library(ggplot2)
library(brunnermunzel)

cat("Starting DEG analysis v6 with Brunner-Munzel + Storey method...\n")
cat("Focus: R0_vs_R1_tumor comparison (most promising for DEGs)\n")

# ==============================================================================
# 1. Load DEGES Normalization Results v6
# ==============================================================================

cat("\nLoading DEGES normalization v6 results...\n")

# Load DEGES normalization results from 08_deges_normalization_v6.R
load("./data/processed/deges_normalization_v6_results.rda")

cat("DEGES normalization v6 results loaded\n")
cat("Available comparisons:", names(deges_normalization_v6_results$deges_results), "\n")

# Display summary
summary_table <- deges_normalization_v6_results$summary_table
print(summary_table)

# ==============================================================================
# 2. Helper Functions for DEG Analysis v6
# ==============================================================================

# Enhanced filterByExpr for statistical testing
apply_statistical_filter_v6 <- function(original_counts, normalized_counts, sample_groups) {
  cat("Applying enhanced filterByExpr for statistical testing...\n")
  
  # Create DGEList with original counts for filtering decision
  dgelist <- DGEList(counts = original_counts, group = factor(sample_groups))
  
  # Apply filterByExpr with enhanced parameters
  keep <- filterByExpr(dgelist, group = factor(sample_groups))
  
  cat(sprintf("Enhanced filterByExpr: %d/%d genes passed (%.1f%%)\n",
              sum(keep), length(keep), sum(keep)/length(keep)*100))
  
  # Apply filter to normalized data
  filtered_normalized <- normalized_counts[keep, ]
  
  return(list(
    filtered_normalized = filtered_normalized,
    filter_passed = keep,
    n_filtered = sum(keep),
    gene_names = rownames(filtered_normalized)
  ))
}

# Enhanced Brunner-Munzel test with progress tracking
perform_brunner_munzel_test_v6 <- function(normalized_data, sample_groups, group1_name, group2_name) {
  cat("Performing enhanced Brunner-Munzel tests...\n")
  
  # Identify group indices
  group1_indices <- which(sample_groups == group1_name)
  group2_indices <- which(sample_groups == group2_name)
  
  cat(sprintf("Testing %s (n=%d) vs %s (n=%d)\n", 
              group1_name, length(group1_indices),
              group2_name, length(group2_indices)))
  
  # Initialize results vectors
  n_genes <- nrow(normalized_data)
  pvalues <- rep(NA, n_genes)
  statistics <- rep(NA, n_genes)
  fold_changes <- rep(NA, n_genes)
  group1_means <- rep(NA, n_genes)
  group2_means <- rep(NA, n_genes)
  group1_medians <- rep(NA, n_genes)
  group2_medians <- rep(NA, n_genes)
  
  # Progress tracking
  progress_interval <- max(1, floor(n_genes / 20))
  
  # Perform test for each gene
  for (i in 1:n_genes) {
    if (i %% progress_interval == 0) {
      cat(sprintf("  Progress: %d/%d genes (%.1f%%)\n", i, n_genes, i/n_genes*100))
    }
    
    # Extract expression values for each group
    group1_values <- as.numeric(normalized_data[i, group1_indices])
    group2_values <- as.numeric(normalized_data[i, group2_indices])
    
    # Calculate means, medians and log2 fold change
    mean1 <- mean(group1_values, na.rm = TRUE)
    mean2 <- mean(group2_values, na.rm = TRUE)
    median1 <- median(group1_values, na.rm = TRUE)
    median2 <- median(group2_values, na.rm = TRUE)
    
    group1_means[i] <- mean1
    group2_means[i] <- mean2
    group1_medians[i] <- median1
    group2_medians[i] <- median2
    
    # Log2 fold change (group2 vs group1)
    # Add small pseudocount to avoid log of zero
    pseudocount <- 1
    fold_changes[i] <- log2(mean2 + pseudocount) - log2(mean1 + pseudocount)
    
    # Brunner-Munzel test
    tryCatch({
      # Check for sufficient variation
      combined_values <- c(group1_values, group2_values)
      if (length(unique(combined_values)) > 1 && 
          !all(is.na(combined_values)) && 
          var(combined_values, na.rm = TRUE) > 0) {
        
        bm_result <- brunnermunzel::brunnermunzel.test(group1_values, group2_values)
        pvalues[i] <- bm_result$p.value
        statistics[i] <- bm_result$statistic
      } else {
        # No variation - set p-value to 1
        pvalues[i] <- 1.0
        statistics[i] <- 0
      }
    }, error = function(e) {
      # Handle any errors in BM test
      pvalues[i] <- 1.0
      statistics[i] <- 0
    })
  }
  
  cat("Enhanced Brunner-Munzel testing completed\n")
  
  # Return results
  return(list(
    pvalues = pvalues,
    statistics = statistics,
    fold_changes = fold_changes,
    group1_means = group1_means,
    group2_means = group2_means,
    group1_medians = group1_medians,
    group2_medians = group2_medians,
    gene_names = rownames(normalized_data),
    group1_name = group1_name,
    group2_name = group2_name,
    n_group1 = length(group1_indices),
    n_group2 = length(group2_indices)
  ))
}

# Enhanced Storey method with bootstrap validation
apply_storey_correction_v6 <- function(pvalues, alpha = 0.05) {
  cat("Applying enhanced Storey method for multiple testing correction...\n")
  
  # Remove NA p-values for qvalue calculation
  valid_pvals <- !is.na(pvalues) & is.finite(pvalues)
  
  if (sum(valid_pvals) == 0) {
    cat("No valid p-values for correction\n")
    return(list(
      qvalues = rep(NA, length(pvalues)),
      pi0 = NA,
      significant = rep(FALSE, length(pvalues)),
      n_significant = 0
    ))
  }
  
  cat(sprintf("Valid p-values for correction: %d/%d\n", sum(valid_pvals), length(pvalues)))
  
  # Apply Storey method with error handling
  tryCatch({
    qvalue_result <- qvalue(pvalues[valid_pvals], pi0.method = "bootstrap")
    
    # Create full qvalue vector
    qvalues <- rep(NA, length(pvalues))
    qvalues[valid_pvals] <- qvalue_result$qvalues
    
    # Identify significant genes
    significant <- qvalues < alpha & !is.na(qvalues)
    n_significant <- sum(significant)
    
    cat(sprintf("Storey method results:\n"))
    cat(sprintf("  pi0 (proportion of non-DE genes): %.3f\n", qvalue_result$pi0))
    cat(sprintf("  Significant genes (q < %.2f): %d/%d (%.2f%%)\n",
                alpha, n_significant, sum(valid_pvals), 
                n_significant/sum(valid_pvals)*100))
    
    # Quality check: pi0 should be reasonable
    if (qvalue_result$pi0 < 0.3) {
      cat("Note: Low pi0 (<0.3), high proportion of DE genes detected\n")
    } else if (qvalue_result$pi0 > 0.95) {
      cat("Note: High pi0 (>0.95), very few DE genes detected\n")
    } else {
      cat("pi0 value appears reasonable\n")
    }
    
    return(list(
      qvalues = qvalues,
      pi0 = qvalue_result$pi0,
      significant = significant,
      n_significant = n_significant,
      alpha = alpha,
      qvalue_object = qvalue_result
    ))
    
  }, error = function(e) {
    cat(sprintf("Error in Storey method: %s\n", e$message))
    cat("Falling back to Benjamini-Hochberg method...\n")
    
    # Fallback to BH method
    qvalues <- rep(NA, length(pvalues))
    qvalues[valid_pvals] <- p.adjust(pvalues[valid_pvals], method = "BH")
    significant <- qvalues < alpha & !is.na(qvalues)
    n_significant <- sum(significant)
    
    cat(sprintf("Benjamini-Hochberg fallback results:\n"))
    cat(sprintf("  Significant genes (BH q < %.2f): %d/%d (%.2f%%)\n",
                alpha, n_significant, sum(valid_pvals), 
                n_significant/sum(valid_pvals)*100))
    
    return(list(
      qvalues = qvalues,
      pi0 = NA,
      significant = significant,
      n_significant = n_significant,
      alpha = alpha,
      method_used = "BH_fallback"
    ))
  })
}

# Enhanced results summary with gene annotation
# Enhanced results summary with gene annotation
create_deg_summary_v6 <- function(bm_result, storey_result, comparison_name) {
  
  # Create comprehensive results data frame
  results_df <- data.frame(
    gene_id = bm_result$gene_names,
    pvalue = bm_result$pvalues,
    qvalue = storey_result$qvalues,
    log2FC = bm_result$fold_changes,
    statistic = bm_result$statistics,
    group1_mean = bm_result$group1_means,
    group2_mean = bm_result$group2_means,
    group1_median = bm_result$group1_medians,
    group2_median = bm_result$group2_medians,
    significant = storey_result$significant,
    stringsAsFactors = FALSE
  )
  
  # Add fold change magnitude and direction
  results_df$abs_log2FC <- abs(results_df$log2FC)
  results_df$direction <- ifelse(results_df$log2FC > 0, "UP", "DOWN")
  
  # Add significance categories
  results_df$significance_category <- "Non-significant"
  results_df$significance_category[results_df$significant & results_df$log2FC > 0] <- "Upregulated"
  results_df$significance_category[results_df$significant & results_df$log2FC < 0] <- "Downregulated"
  
  # Load gene info if available
  if (!exists("gene_info")) {
    if (!exists("se_thyr")) {
      load("./data/raw/thyr_data.rda")
      se_thyr <- data
      rm(data)
    }
    gene_info <- rowData(se_thyr)
  }
  
  # Add gene annotation if available - FIXED S4 handling
  if (any(results_df$gene_id %in% rownames(gene_info))) {
    matching_genes <- results_df$gene_id[results_df$gene_id %in% rownames(gene_info)]
    
    # Extract annotation data properly from S4 object
    # Convert to data.frame if it's an S4 object
    if (is(gene_info, "DataFrame")) {
      # It's a DataFrame from S4Vectors
      gene_annotation <- as.data.frame(gene_info[matching_genes, c("gene_name", "gene_type")])
    } else {
      # It's already a regular data.frame
      gene_annotation <- gene_info[matching_genes, c("gene_name", "gene_type")]
    }
    
    # Initialize columns
    results_df$gene_name <- NA
    results_df$gene_type <- NA
    
    # Match and assign annotations
    for (i in 1:nrow(results_df)) {
      gene_id <- results_df$gene_id[i]
      if (gene_id %in% rownames(gene_annotation)) {
        results_df$gene_name[i] <- as.character(gene_annotation[gene_id, "gene_name"])
        results_df$gene_type[i] <- as.character(gene_annotation[gene_id, "gene_type"])
      }
    }
  } else {
    # No matching genes - add empty columns
    results_df$gene_name <- NA
    results_df$gene_type <- NA
  }
  
  # Sort by q-value
  results_df <- results_df[order(results_df$qvalue, na.last = TRUE), ]
  
  # Summary statistics
  summary_stats <- list(
    comparison = comparison_name,
    total_genes_tested = nrow(results_df),
    valid_tests = sum(!is.na(results_df$pvalue)),
    significant_genes = storey_result$n_significant,
    pi0 = storey_result$pi0,
    alpha = storey_result$alpha,
    upregulated = sum(results_df$significance_category == "Upregulated", na.rm = TRUE),
    downregulated = sum(results_df$significance_category == "Downregulated", na.rm = TRUE),
    group1_name = bm_result$group1_name,
    group2_name = bm_result$group2_name,
    n_group1 = bm_result$n_group1,
    n_group2 = bm_result$n_group2,
    method_used = ifelse(is.null(storey_result$method_used), "Storey", storey_result$method_used)
  )
  
  return(list(
    results_df = results_df,
    summary_stats = summary_stats
  ))
}
# ==============================================================================
# 3. Process Each Comparison for DEG Analysis
# ==============================================================================

cat("\n==============================================\n")
cat("DEG Analysis for Each Comparison\n")
cat("==============================================\n")

# Initialize results storage
deg_analysis_v6_results <- list()

# Process each comparison from DEGES v6 results
for (comp_name in names(deges_normalization_v6_results$deges_results)) {
  cat(sprintf("\n### Processing %s ###\n", comp_name))
  
  # Get DEGES results for this comparison
  deges_comp_result <- deges_normalization_v6_results$deges_results[[comp_name]]
  
  # Extract key information
  normalized_counts <- deges_comp_result$final_normalized_counts
  samples <- deges_comp_result$samples
  groups <- deges_comp_result$groups
  comp_info <- deges_comp_result$comparison_info
  
  cat(sprintf("DEGES normalized data: %d genes, %d samples\n",
              nrow(normalized_counts), ncol(normalized_counts)))
  
  # Get original count data for the same genes and samples
  if (!exists("se_thyr")) {
    load("./data/raw/thyr_data.rda")
    se_thyr <- data
    rm(data)
  }
  
  count_data_full <- assay(se_thyr, "stranded_second")
  common_genes <- intersect(rownames(count_data_full), rownames(normalized_counts))
  original_counts_subset <- count_data_full[common_genes, samples]
  normalized_counts_subset <- normalized_counts[common_genes, ]
  
  cat(sprintf("Matched original counts: %d genes\n", length(common_genes)))
  
  # Apply enhanced statistical filtering
  filter_result <- apply_statistical_filter_v6(
    original_counts_subset, 
    normalized_counts_subset, 
    groups
  )
  
  # Perform enhanced Brunner-Munzel test
  bm_result <- perform_brunner_munzel_test_v6(
    filter_result$filtered_normalized,
    groups,
    comp_info$group1,
    comp_info$group2
  )
  
  # Apply enhanced Storey correction
  storey_result <- apply_storey_correction_v6(bm_result$pvalues)
  
  # Create comprehensive summary
  deg_summary <- create_deg_summary_v6(bm_result, storey_result, comp_name)
  
  # Store results
  deg_analysis_v6_results[[comp_name]] <- list(
    comparison_info = comp_info,
    filter_result = filter_result,
    bm_result = bm_result,
    storey_result = storey_result,
    deg_summary = deg_summary,
    samples_used = samples,
    groups_used = groups,
    analysis_date = Sys.time(),
    analysis_version = "v6_enhanced_bm_storey"
  )
  
  # Print summary
  cat(sprintf("\n%s Results Summary:\n", comp_name))
  cat(sprintf("  Genes tested: %d\n", deg_summary$summary_stats$total_genes_tested))
  cat(sprintf("  Significant DEGs: %d (%.2f%%)\n", 
              deg_summary$summary_stats$significant_genes,
              deg_summary$summary_stats$significant_genes / deg_summary$summary_stats$total_genes_tested * 100))
  cat(sprintf("  Upregulated: %d\n", deg_summary$summary_stats$upregulated))
  cat(sprintf("  Downregulated: %d\n", deg_summary$summary_stats$downregulated))
  if (!is.na(deg_summary$summary_stats$pi0)) {
    cat(sprintf("  Pi0 estimate: %.3f\n", deg_summary$summary_stats$pi0))
  }
  cat(sprintf("  Method used: %s\n", deg_summary$summary_stats$method_used))
  
  cat(sprintf("%s: DEG analysis completed\n", comp_name))
}

# ==============================================================================
# 4. Create Overall Summary
# ==============================================================================

cat("\n==============================================\n")
cat("Overall DEG Analysis v6 Summary\n")
cat("==============================================\n")

# Create summary table
deg_summary_table_v6 <- data.frame(
  Comparison = character(0),
  Group1 = character(0),
  Group2 = character(0),
  Tissue = character(0),
  N1 = integer(0),
  N2 = integer(0),
  Genes_Tested = integer(0),
  DEGs_Total = integer(0),
  DEGs_Up = integer(0),
  DEGs_Down = integer(0),
  DEG_Rate = numeric(0),
  Pi0 = numeric(0),
  Method = character(0),
  stringsAsFactors = FALSE
)

for (comp_name in names(deg_analysis_v6_results)) {
  result <- deg_analysis_v6_results[[comp_name]]
  summary_stats <- result$deg_summary$summary_stats
  
  deg_summary_table_v6 <- rbind(deg_summary_table_v6, data.frame(
    Comparison = comp_name,
    Group1 = summary_stats$group1_name,
    Group2 = summary_stats$group2_name,
    Tissue = result$comparison_info$tissue,
    N1 = summary_stats$n_group1,
    N2 = summary_stats$n_group2,
    Genes_Tested = summary_stats$total_genes_tested,
    DEGs_Total = summary_stats$significant_genes,
    DEGs_Up = summary_stats$upregulated,
    DEGs_Down = summary_stats$downregulated,
    DEG_Rate = round(summary_stats$significant_genes / summary_stats$total_genes_tested * 100, 2),
    Pi0 = round(summary_stats$pi0, 3),
    Method = summary_stats$method_used,
    stringsAsFactors = FALSE
  ))
}

print(deg_summary_table_v6)

# ==============================================================================
# 5. Focus on Primary Comparison (R0_vs_R1_tumor)
# ==============================================================================

cat("\n=== Primary Comparison Analysis (R0_vs_R1_tumor) ===\n")

if ("R0_vs_R1_tumor" %in% names(deg_analysis_v6_results)) {
  primary_result <- deg_analysis_v6_results[["R0_vs_R1_tumor"]]
  primary_degs <- primary_result$deg_summary$results_df
  
  # Extract significant DEGs
  sig_degs <- primary_degs[primary_degs$significant, ]
  
  cat(sprintf("R0_vs_R1_tumor detailed results:\n"))
  cat(sprintf("  Total significant DEGs: %d\n", nrow(sig_degs)))
  
  if (nrow(sig_degs) > 0) {
    # Top upregulated genes
    top_up <- head(sig_degs[sig_degs$direction == "UP", ], 10)
    cat(sprintf("\n  Top %d upregulated genes:\n", min(10, nrow(top_up))))
    if (nrow(top_up) > 0) {
      for (i in 1:min(5, nrow(top_up))) {
        gene <- top_up[i, ]
        cat(sprintf("    %s: log2FC=%.2f, q=%.3e\n", 
                    gene$gene_id, gene$log2FC, gene$qvalue))
      }
    }
    
    # Top downregulated genes
    top_down <- head(sig_degs[sig_degs$direction == "DOWN", ], 10)
    cat(sprintf("\n  Top %d downregulated genes:\n", min(10, nrow(top_down))))
    if (nrow(top_down) > 0) {
      for (i in 1:min(5, nrow(top_down))) {
        gene <- top_down[i, ]
        cat(sprintf("    %s: log2FC=%.2f, q=%.3e\n", 
                    gene$gene_id, gene$log2FC, gene$qvalue))
      }
    }
    
    # Expression level distribution
    cat("\n  Expression level distribution of DEGs:\n")
    mean_expr <- rowMeans(cbind(sig_degs$group1_mean, sig_degs$group2_mean))
    expr_quantiles <- quantile(mean_expr, probs = c(0, 0.25, 0.5, 0.75, 1))
    cat(sprintf("    Min: %.1f, Q1: %.1f, Median: %.1f, Q3: %.1f, Max: %.1f\n",
                expr_quantiles[1], expr_quantiles[2], expr_quantiles[3], 
                expr_quantiles[4], expr_quantiles[5]))
  }
} else {
  cat("R0_vs_R1_tumor comparison not available\n")
}

# ==============================================================================
# 6. Save Results
# ==============================================================================

cat("\nSaving DEG analysis v6 results...\n")

# Create processed directory if it doesn't exist
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save comprehensive DEG v6 results
final_deg_v6_results <- list(
  deg_analysis_results = deg_analysis_v6_results,
  deg_summary_table = deg_summary_table_v6,
  analysis_parameters = list(
    statistical_test = "Brunner-Munzel",
    multiple_correction = "Storey",
    significance_threshold = 0.05,
    filter_method = "enhanced_filterByExpr",
    analysis_version = "v6"
  ),
  total_comparisons = nrow(deg_summary_table_v6),
  analysis_date = Sys.time(),
  analysis_version = "v6_enhanced_bm_storey"
)

save(final_deg_v6_results, file = "./data/processed/final_deg_v6_results.rda")

# Save individual comparison results for easy access
if (!dir.exists("./output/reports")) {
  dir.create("./output/reports", recursive = TRUE)
}

for (comp_name in names(deg_analysis_v6_results)) {
  results_df <- deg_analysis_v6_results[[comp_name]]$deg_summary$results_df
  write.csv(results_df, 
            file = paste0("./output/reports/deg_results_v6_", comp_name, ".csv"),
            row.names = FALSE)
}

cat("Results saved to ./data/processed/ and ./output/reports/\n")

# ==============================================================================
# 7. Final Summary and Next Steps
# ==============================================================================

cat("\n==============================================\n")
cat("DEG Analysis v6 Final Summary\n")
cat("==============================================\n")

cat("Enhanced Brunner-Munzel + Storey Method Analysis Completed!\n\n")

cat("Key achievements:\n")
cat("- DEGES v6 normalized data integration\n")
cat("- Robust Brunner-Munzel test implementation\n")
cat("- Storey method with bootstrap validation\n")
cat("- Comprehensive gene annotation\n")

total_degs <- sum(deg_summary_table_v6$DEGs_Total)
total_genes_tested <- sum(deg_summary_table_v6$Genes_Tested)

cat(sprintf("\nOverall Results:\n"))
cat(sprintf("  Total comparisons: %d\n", nrow(deg_summary_table_v6)))
cat(sprintf("  Total genes tested: %d\n", total_genes_tested))
cat(sprintf("  Total DEGs identified: %d\n", total_degs))

if (nrow(deg_summary_table_v6) > 0) {
  cat(sprintf("\nComparison-wise DEG counts:\n"))
  for (i in 1:nrow(deg_summary_table_v6)) {
    row <- deg_summary_table_v6[i, ]
    cat(sprintf("  %s: %d DEGs (%.1f%%) [pi0=%.3f]\n", 
                row$Comparison, row$DEGs_Total, row$DEG_Rate, row$Pi0))
  }
}

cat("\nNext steps:\n")
cat("1. Review DEG results, especially R0_vs_R1_tumor\n")
cat("2. Consider consistency analysis if needed (10_deg_overlap_analysis_v6.R)\n")
cat("3. Proceed to feature selection based on tumor DEGs\n")

cat("\nDEG analysis v6 completed successfully!\n")
cat("==============================================\n")