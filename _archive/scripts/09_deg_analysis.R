# ==============================================================================
# REBC-THYR DEG Analysis Script
# 09_deg_analysis.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(edgeR)
library(qvalue)
library(dplyr)
library(ggplot2)
library(brunnermunzel)

cat("Starting DEG analysis with Brunner-Munzel test + Storey method...\n")

# ==============================================================================
# 1. Load Data and DEGES Normalization Results
# ==============================================================================

cat("Loading DEGES normalization results...\n")

# Load DEGES normalization results from 08
load("./data/processed/deges_normalization_results.rda")

# Load original count data for filterByExpr
if (!exists("se_thyr")) {
  load("./data/raw/thyr_data.rda")
  se_thyr <- data
  rm(data)
}

count_data_full <- assay(se_thyr, "stranded_second")
gene_info <- rowData(se_thyr)

cat("DEGES normalization results loaded\n")
cat("Available comparisons:", names(deges_normalization_results$deges_results), "\n")

# ==============================================================================
# 2. Helper Functions
# ==============================================================================

# Function to apply filterByExpr for statistical testing
apply_statistical_filter <- function(original_counts, normalized_counts, sample_groups) {
  cat("Applying filterByExpr for statistical testing...\n")
  
  # Create DGEList with original counts for filtering decision
  dgelist <- DGEList(counts = original_counts, group = factor(sample_groups))
  
  # Apply filterByExpr
  keep <- filterByExpr(dgelist, min.count = 1, min.total.count = 15)
  
  cat(sprintf("FilterByExpr: %d/%d genes passed (%.1f%%)\n",
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

# Function to perform Brunner-Munzel test
perform_brunner_munzel_test <- function(normalized_data, sample_groups, group1_name, group2_name) {
  cat("Performing Brunner-Munzel tests...\n")
  
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
    
    # Calculate means and log2 fold change
    mean1 <- mean(group1_values, na.rm = TRUE)
    mean2 <- mean(group2_values, na.rm = TRUE)
    
    group1_means[i] <- mean1
    group2_means[i] <- mean2
    
    # Log2 fold change (group2 vs group1)
    # Add small pseudocount to avoid log(0)
    pseudocount <- 1e-6
    fold_changes[i] <- log2((mean2 + pseudocount) / (mean1 + pseudocount))
    
    # Brunner-Munzel test
    tryCatch({
      # Check for sufficient variation
      if (length(unique(c(group1_values, group2_values))) > 1) {
        bm_result <- brunnermunzel::brunnermunzel.test(group1_values, group2_values)
        pvalues[i] <- bm_result$p.value
        statistics[i] <- bm_result$statistic
      } else {
        # No variation - set p-value to 1
        pvalues[i] <- 1.0
        statistics[i] <- 0
      }
    }, error = function(e) {
      cat(sprintf("Error in Brunner-Munzel test for gene %d: %s\n", i, e$message))
      stop("Brunner-Munzel test failed. Please review data or test parameters.")
    })
  }
  
  cat("Brunner-Munzel testing completed\n")
  
  # Return results
  return(list(
    pvalues = pvalues,
    statistics = statistics,
    fold_changes = fold_changes,
    group1_means = group1_means,
    group2_means = group2_means,
    gene_names = rownames(normalized_data),
    group1_name = group1_name,
    group2_name = group2_name,
    n_group1 = length(group1_indices),
    n_group2 = length(group2_indices)
  ))
}

# Function to apply Storey method for multiple testing correction
apply_storey_correction <- function(pvalues, alpha = 0.05) {
  cat("Applying Storey method for multiple testing correction...\n")
  
  # Remove NA p-values for qvalue calculation
  valid_pvals <- !is.na(pvalues)
  
  if (sum(valid_pvals) == 0) {
    cat("No valid p-values for correction\n")
    return(list(
      qvalues = rep(NA, length(pvalues)),
      pi0 = NA,
      significant = rep(FALSE, length(pvalues)),
      n_significant = 0
    ))
  }
  
  # Apply Storey method
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
    
    return(list(
      qvalues = qvalues,
      pi0 = qvalue_result$pi0,
      significant = significant,
      n_significant = n_significant,
      alpha = alpha
    ))
    
  }, error = function(e) {
    cat(sprintf("Error in Storey method: %s\n", e$message))
    stop("Storey method failed. Please review p-value distribution or method parameters.")
  })
}

# Function to create results summary
create_deg_summary <- function(bm_result, storey_result, comparison_name) {
  
  # Create comprehensive results data frame
  results_df <- data.frame(
    gene_id = bm_result$gene_names,
    pvalue = bm_result$pvalues,
    qvalue = storey_result$qvalues,
    log2FC = bm_result$fold_changes,
    statistic = bm_result$statistics,
    group1_mean = bm_result$group1_means,
    group2_mean = bm_result$group2_means,
    significant = storey_result$significant,
    stringsAsFactors = FALSE
  )
  
  # Add gene annotation if available
  if (exists("gene_info")) {
    gene_annotation <- gene_info[results_df$gene_id, c("gene_name", "gene_type")]
    results_df <- cbind(results_df, gene_annotation)
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
    upregulated = sum(results_df$significant & results_df$log2FC > 0, na.rm = TRUE),
    downregulated = sum(results_df$significant & results_df$log2FC < 0, na.rm = TRUE),
    group1_name = bm_result$group1_name,
    group2_name = bm_result$group2_name,
    n_group1 = bm_result$n_group1,
    n_group2 = bm_result$n_group2
  )
  
  return(list(
    results_df = results_df,
    summary_stats = summary_stats
  ))
}

# ==============================================================================
# 3. Process Each Comparison
# ==============================================================================

cat("\n==============================================\n")
cat("DEG Analysis for Each Comparison\n")
cat("==============================================\n")

# Initialize results storage
deg_analysis_results <- list()

# Process each comparison from DEGES results
for (comp_name in names(deges_normalization_results$deges_results)) {
  cat(sprintf("\n### Processing %s ###\n", comp_name))
  
  # Get DEGES results for this comparison
  deges_comp_result <- deges_normalization_results$deges_results[[comp_name]]
  
  # Extract key information
  normalized_counts <- deges_comp_result$final_normalized_counts
  samples <- deges_comp_result$samples
  groups <- deges_comp_result$groups
  comp_info <- deges_comp_result$comparison_info
  
  cat(sprintf("DEGES normalized data: %d genes, %d samples\n",
              nrow(normalized_counts), ncol(normalized_counts)))
  
  # Get original count data for the same genes and samples
  common_genes <- intersect(rownames(count_data_full), rownames(normalized_counts))
  original_counts_subset <- count_data_full[common_genes, samples]
  normalized_counts_subset <- normalized_counts[common_genes, ]
  
  cat(sprintf("Matched original counts: %d genes\n", length(common_genes)))
  
  # Apply statistical filtering (filterByExpr)
  filter_result <- apply_statistical_filter(
    original_counts_subset, 
    normalized_counts_subset, 
    groups
  )
  
  # Perform Brunner-Munzel test
  bm_result <- perform_brunner_munzel_test(
    filter_result$filtered_normalized,
    groups,
    comp_info$group1,
    comp_info$group2
  )
  
  # Apply Storey correction
  storey_result <- apply_storey_correction(bm_result$pvalues)
  
  # Create comprehensive summary
  deg_summary <- create_deg_summary(bm_result, storey_result, comp_name)
  
  # Store results
  deg_analysis_results[[comp_name]] <- list(
    comparison_info = comp_info,
    filter_result = filter_result,
    bm_result = bm_result,
    storey_result = storey_result,
    deg_summary = deg_summary,
    samples_used = samples,
    groups_used = groups,
    analysis_date = Sys.time()
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
  
  cat(sprintf("%s: DEG analysis completed\n", comp_name))
}

# ==============================================================================
# 4. Create Overall Summary
# ==============================================================================

cat("\n==============================================\n")
cat("Overall DEG Analysis Summary\n")
cat("==============================================\n")

# Create summary table
deg_summary_table <- data.frame(
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
  stringsAsFactors = FALSE
)

for (comp_name in names(deg_analysis_results)) {
  result <- deg_analysis_results[[comp_name]]
  summary_stats <- result$deg_summary$summary_stats
  
  deg_summary_table <- rbind(deg_summary_table, data.frame(
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
    stringsAsFactors = FALSE
  ))
}

print(deg_summary_table)

# ==============================================================================
# 5. Quality Assessment and Validation
# ==============================================================================

cat("\nQuality Assessment:\n")

# Check for expected patterns based on research hypothesis
primary_comparison <- "R0_vs_R1_tumor"
if (primary_comparison %in% names(deg_analysis_results)) {
  primary_result <- deg_analysis_results[[primary_comparison]]
  primary_summary <- primary_result$deg_summary$summary_stats
  
  cat(sprintf("\nPrimary comparison (%s) assessment:\n", primary_comparison))
  cat(sprintf("  Expected: ~3000 DEGs based on research hypothesis\n"))
  cat(sprintf("  Observed: %d DEGs\n", primary_summary$significant_genes))
  
  if (primary_summary$significant_genes >= 1000) {
    cat("  ✅ DEG count within expected range\n")
  } else if (primary_summary$significant_genes >= 500) {
    cat("  ⚠️ DEG count lower than expected but reasonable\n")
  } else {
    cat("  ❌ DEG count much lower than expected - review parameters\n")
  }
}

# Cross-comparison validation
ret_tumor_degs <- ifelse("R0_vs_R1_tumor" %in% names(deg_analysis_results),
                         deg_analysis_results[["R0_vs_R1_tumor"]]$deg_summary$summary_stats$significant_genes, 0)
ret_normal_degs <- ifelse("R0_vs_R1_normal" %in% names(deg_analysis_results),
                          deg_analysis_results[["R0_vs_R1_normal"]]$deg_summary$summary_stats$significant_genes, 0)
braf_normal_degs <- ifelse("B0_vs_B1_normal" %in% names(deg_analysis_results),
                           deg_analysis_results[["B0_vs_B1_normal"]]$deg_summary$summary_stats$significant_genes, 0)

cat(sprintf("\nCross-comparison validation:\n"))
cat(sprintf("  RET Tumor DEGs: %d\n", ret_tumor_degs))
cat(sprintf("  RET Normal DEGs: %d\n", ret_normal_degs))
cat(sprintf("  BRAF Normal DEGs: %d\n", braf_normal_degs))

if (ret_tumor_degs > ret_normal_degs && ret_tumor_degs > braf_normal_degs) {
  cat("  ✅ Expected pattern: RET Tumor shows most DEGs\n")
} else {
  cat("  ⚠️ Unexpected pattern in DEG distribution\n")
}

# ==============================================================================
# 6. Create Visualizations
# ==============================================================================

cat("Creating visualization plots...\n")

# Create plots directory if it doesn't exist
if (!dir.exists("./output/plots")) {
  dir.create("./output/plots", recursive = TRUE)
}

# Volcano plots for each comparison
for (comp_name in names(deg_analysis_results)) {
  result <- deg_analysis_results[[comp_name]]
  results_df <- result$deg_summary$results_df
  
  # Skip if no valid results
  if (nrow(results_df) == 0 || all(is.na(results_df$pvalue))) {
    cat(sprintf("Skipping volcano plot for %s (no valid results)\n", comp_name))
    next
  }
  
  # Create volcano plot data
  plot_data <- results_df[!is.na(results_df$pvalue) & !is.na(results_df$log2FC), ]
  
  if (nrow(plot_data) < 10) {
    cat(sprintf("Skipping volcano plot for %s (insufficient data points)\n", comp_name))
    next
  }
  
  # Add significance categories
  plot_data$category <- "Non-significant"
  plot_data$category[plot_data$significant & plot_data$log2FC > 0] <- "Upregulated"
  plot_data$category[plot_data$significant & plot_data$log2FC < 0] <- "Downregulated"
  
  # Create volcano plot
  volcano_plot <- ggplot(plot_data, aes(x = log2FC, y = -log10(pvalue), color = category)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("Upregulated" = "red", 
                                  "Downregulated" = "blue", 
                                  "Non-significant" = "gray")) +
    labs(title = paste("Volcano Plot:", comp_name),
         subtitle = paste0("Significant DEGs: ", sum(plot_data$significant)),
         x = "Log2 Fold Change",
         y = "-Log10 P-value",
         color = "Category") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5)
  
  # Save plot
  ggsave(filename = paste0("./output/plots/volcano_", comp_name, ".png"),
         plot = volcano_plot, width = 8, height = 6, dpi = 300)
  
  cat(sprintf("Volcano plot saved for %s\n", comp_name))
}

# ==============================================================================
# 7. Save Results
# ==============================================================================

cat("Saving DEG analysis results...\n")

# Create processed directory if it doesn't exist
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save comprehensive DEG results
final_deg_results <- list(
  deg_analysis_results = deg_analysis_results,
  deg_summary_table = deg_summary_table,
  analysis_parameters = list(
    statistical_test = "Brunner-Munzel",
    multiple_correction = "Storey",
    significance_threshold = 0.05,
    filter_method = "filterByExpr"
  ),
  total_comparisons = nrow(deg_summary_table),
  analysis_date = Sys.time(),
  analysis_version = "v1_bm_storey_integration"
)

save(final_deg_results, file = "./data/processed/final_deg_results.rda")

# Save individual comparison results for easy access
for (comp_name in names(deg_analysis_results)) {
  results_df <- deg_analysis_results[[comp_name]]$deg_summary$results_df
  write.csv(results_df, 
            file = paste0("./output/reports/deg_results_", comp_name, ".csv"),
            row.names = FALSE)
}

cat("Results saved to ./data/processed/ and ./output/reports/\n")

# ==============================================================================
# 8. Final Summary and Next Steps
# ==============================================================================

cat("\n==============================================\n")
cat("DEG Analysis Final Summary\n")
cat("==============================================\n")

cat("Brunner-Munzel + Storey Method Analysis Completed Successfully!\n\n")

cat("Key achievements:\n")
cat("✅ DEGES normalized data integration\n")
cat("✅ Statistical filtering (filterByExpr) applied\n")
cat("✅ Brunner-Munzel test for robust comparison\n")
cat("✅ Storey method for mixed distribution correction\n")
cat("✅ Comprehensive quality assessment\n")

total_degs <- sum(deg_summary_table$DEGs_Total)
total_genes_tested <- sum(deg_summary_table$Genes_Tested)

cat(sprintf("\nOverall Results:\n"))
cat(sprintf("  Total comparisons: %d\n", nrow(deg_summary_table)))
cat(sprintf("  Total genes tested: %d\n", total_genes_tested))
cat(sprintf("  Total DEGs identified: %d\n", total_degs))
cat(sprintf("  Overall DEG rate: %.2f%%\n", total_degs / total_genes_tested * 100))

if (nrow(deg_summary_table) > 0) {
  cat(sprintf("\nComparison-wise DEG counts:\n"))
  for (i in 1:nrow(deg_summary_table)) {
    row <- deg_summary_table[i, ]
    cat(sprintf("  %s: %d DEGs (%.1f%%)\n", 
                row$Comparison, row$DEGs_Total, row$DEG_Rate))
  }
}

cat("\nOutput files created:\n")
cat("  ./data/processed/final_deg_results.rda - Complete analysis results\n")
cat("  ./output/reports/deg_results_*.csv - Individual comparison results\n")
cat("  ./output/plots/volcano_*.png - Volcano plots\n")

cat("\nNext steps:\n")
cat("1. Review DEG results and volcano plots\n")
cat("2. Feature selection for 12 → 1-2 gene pairs\n")
cat("3. Pathway enrichment analysis (optional)\n")
cat("4. Predictive model development\n")

cat("\nDEG analysis (09) completed successfully!\n")
cat("Two-stage analysis (08: DEGES + 09: BM+Storey) achieved!\n")
cat("==============================================\n")

