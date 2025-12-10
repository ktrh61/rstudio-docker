# 09_deg_analysis.R - DEG Analysis with Brunner-Munzel and Storey Method
# Purpose: Perform differential expression analysis on DEGES-normalized data
# Method: Brunner-Munzel test with Storey (qvalue, lambda=0.5) correction
# Input: analysis_dgelist_*.rds (from 07_deges_normalization.R)
# Output: thyr_deg_results.rds with DEG lists and consistency evaluation
# Version: v7.6 - Fixed prior.count to 0 (consistent with 07)
# Date: 2025-12-09

source("analysis_v7/setup.R")

cat("\n=== DEG Analysis with Brunner-Munzel + Storey Method (v7.6) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Method: Brunner-Munzel test with Storey correction (lambda=0.5)\n")
cat("Effect size: PI (primary) with Cliff's delta (derived)\n")
cat("Focus: R0 vs R1, B0 vs B1 comparisons (tumor/normal)\n")

# Load packages
suppressPackageStartupMessages({
  library(edgeR)
  library(qvalue)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(brunnermunzel)
})

# ============================================================================
# Configuration
# ============================================================================

CONFIG <- list(
  # Statistical testing
  ALPHA = 0.05,               # Significance threshold for q-values
  
  # Output control
  VERBOSE = TRUE              # Verbose output
)

cat("\nConfiguration:\n")
cat("  Significance threshold (q-value):", CONFIG$ALPHA, "\n")

# ============================================================================
# Helper functions for DEG analysis
# ============================================================================

# Prepare normalized counts from DGEList
prepare_normalized_counts <- function(dgelist) {
  # Calculate normalized counts (CPM)
  # norm.factors already in dgelist$samples
  normalized_counts <- cpm(dgelist, normalized.lib.sizes = TRUE, 
                           prior.count = 0, log = FALSE)
  
  return(normalized_counts)
}

# Perform Brunner-Munzel test
perform_brunner_munzel_test <- function(normalized_data, sample_groups, 
                                        group1_name, group2_name) {
  cat("  Performing Brunner-Munzel tests...\n")
  
  # Identify group indices
  group1_indices <- which(sample_groups == group1_name)
  group2_indices <- which(sample_groups == group2_name)
  
  cat(sprintf("    Testing %s (n=%d) vs %s (n=%d)\n", 
              group1_name, length(group1_indices),
              group2_name, length(group2_indices)))
  
  # Initialize results
  n_genes <- nrow(normalized_data)
  pvalues <- rep(NA_real_, n_genes)
  statistics <- rep(NA_real_, n_genes)
  fold_changes <- rep(NA_real_, n_genes)
  PI <- rep(NA_real_, n_genes)  # Probability of Superiority: P(Y>X) + 0.5*P(Y=X)
  group1_means <- rep(NA_real_, n_genes)
  group2_means <- rep(NA_real_, n_genes)
  
  # Progress tracking
  progress_interval <- max(1, floor(n_genes / 10))
  
  # Perform test for each gene
  for (i in seq_len(n_genes)) {
    if (i %% progress_interval == 0 && CONFIG$VERBOSE) {
      cat(sprintf("      Progress: %d/%d genes (%.0f%%)\n", 
                  i, n_genes, i/n_genes*100))
    }
    
    # Extract expression values
    group1_values <- as.numeric(normalized_data[i, group1_indices])
    group2_values <- as.numeric(normalized_data[i, group2_indices])
    
    # Calculate means
    mean1 <- mean(group1_values, na.rm = TRUE)
    mean2 <- mean(group2_values, na.rm = TRUE)
    group1_means[i] <- mean1
    group2_means[i] <- mean2
    
    # Log2 fold change (group2 vs group1)
    pseudocount <- 1
    fold_changes[i] <- log2(mean2 + pseudocount) - log2(mean1 + pseudocount)
    
    # Brunner-Munzel test
    # Check for sufficient variation
    if (length(unique(c(group1_values, group2_values))) > 1 && 
        var(c(group1_values, group2_values), na.rm = TRUE) > 0) {
      
      bm_out <- brunnermunzel::brunnermunzel.test(group1_values, group2_values)
      pvalues[i] <- bm_out$p.value
      statistics[i] <- bm_out$statistic
      PI[i] <- bm_out$estimate  # P(Y>X) + 0.5*P(Y=X), primary output
    } else {
      # No variation: set neutral values
      pvalues[i] <- 1.0
      statistics[i] <- 0
      PI[i] <- 0.5
    }
  }
  
  cat("    Brunner-Munzel testing completed\n")
  
  return(list(
    pvalues = pvalues,
    statistics = statistics,
    fold_changes = fold_changes,
    PI = PI,  # Primary output from Brunner-Munzel
    group1_means = group1_means,
    group2_means = group2_means,
    gene_names = rownames(normalized_data),
    group1_name = group1_name,
    group2_name = group2_name,
    n_group1 = length(group1_indices),
    n_group2 = length(group2_indices)
  ))
}

# Apply Storey correction
apply_storey_correction <- function(pvalues, alpha = 0.05) {
  cat("  Applying Storey method for multiple testing correction...\n")
  
  # Remove NA p-values
  valid_pvals <- !is.na(pvalues) & is.finite(pvalues)
  
  if (sum(valid_pvals) == 0) {
    cat("    No valid p-values for correction\n")
    return(list(
      qvalues = rep(NA_real_, length(pvalues)),
      pi0 = NA_real_,
      significant = rep(FALSE, length(pvalues)),
      n_significant = 0L
    ))
  }
  
  cat(sprintf("    Valid p-values: %d/%d\n", sum(valid_pvals), length(pvalues)))
  
  # Apply Storey method (lambda = 0.5 fixed for reproducibility)
  qvalue_result <- qvalue(pvalues[valid_pvals], lambda = 0.5)
  
  # Create full qvalue vector
  qvalues <- rep(NA_real_, length(pvalues))
  qvalues[valid_pvals] <- qvalue_result$qvalues
  
  # Identify significant genes
  significant <- qvalues < alpha & !is.na(qvalues)
  n_significant <- sum(significant)
  
  cat(sprintf("    pi0 estimate: %.3f (lambda=0.5)\n", qvalue_result$pi0))
  cat(sprintf("    Significant genes (q < %.2f): %d (%.1f%%)\n",
              alpha, n_significant, n_significant/sum(valid_pvals)*100))
  
  return(list(
    qvalues = qvalues,
    pi0 = qvalue_result$pi0,
    significant = significant,
    n_significant = n_significant,
    alpha = alpha
  ))
}

# Create DEG results summary
create_deg_summary <- function(bm_result, storey_result, comparison_name, gene_info_df) {
  
  # Create results data frame
  # PI is primary output, cliffs_delta is derived (single conversion)
  results_df <- data.frame(
    gene_id = bm_result$gene_names,
    pvalue = bm_result$pvalues,
    qvalue = storey_result$qvalues,
    log2FC = bm_result$fold_changes,
    PI = bm_result$PI,                      # Primary output from Brunner-Munzel
    cliffs_delta = 2 * bm_result$PI - 1,    # Derived: delta = 2*PI - 1
    statistic = bm_result$statistics,
    group1_mean = bm_result$group1_means,
    group2_mean = bm_result$group2_means,
    significant = storey_result$significant,
    stringsAsFactors = FALSE
  )
  
  # Add fold change magnitude and direction
  results_df$abs_log2FC <- abs(results_df$log2FC)
  results_df$abs_cliffs_delta <- abs(results_df$cliffs_delta)
  results_df$direction <- ifelse(results_df$log2FC > 0, "UP", "DOWN")
  
  # Add gene annotation
  if (!is.null(gene_info_df) && nrow(gene_info_df) > 0) {
    # Match gene annotations
    matching_idx <- match(results_df$gene_id, rownames(gene_info_df))
    results_df$gene_name <- gene_info_df$gene_name[matching_idx]
    results_df$gene_type <- gene_info_df$gene_type[matching_idx]
  } else {
    results_df$gene_name <- NA_character_
    results_df$gene_type <- NA_character_
  }
  
  # Sort by q-value
  results_df <- results_df[order(results_df$qvalue, na.last = TRUE), ]
  
  # Helper function for descriptive statistics
  calc_effect_size_stats <- function(delta_values) {
    delta_values <- delta_values[!is.na(delta_values)]
    if (length(delta_values) == 0) {
      return(list(
        n = 0L,
        min = NA_real_, q1 = NA_real_, median = NA_real_, 
        q3 = NA_real_, max = NA_real_,
        mean = NA_real_, sd = NA_real_, iqr = NA_real_
      ))
    }
    quants <- quantile(delta_values, probs = c(0, 0.25, 0.5, 0.75, 1))
    list(
      n = length(delta_values),
      min = quants[1], q1 = quants[2], median = quants[3], 
      q3 = quants[4], max = quants[5],
      mean = mean(delta_values), 
      sd = sd(delta_values),
      iqr = quants[4] - quants[2]
    )
  }
  
  # Effect size statistics for ALL genes (PI-based, primary)
  # Using cliffs_delta for compatibility with existing visualizations
  all_delta_stats <- calc_effect_size_stats(results_df$cliffs_delta)
  # Also compute absolute delta median for reference
  all_abs_delta_median <- median(results_df$abs_cliffs_delta, na.rm = TRUE)
  
  # PI statistics for ALL genes
  all_pi_stats <- calc_effect_size_stats(results_df$PI)
  
  # Effect size statistics for DEGs
  sig_genes <- results_df[results_df$significant, ]
  deg_delta_stats <- calc_effect_size_stats(sig_genes$cliffs_delta)
  deg_abs_delta_median <- if (nrow(sig_genes) > 0) {
    median(sig_genes$abs_cliffs_delta, na.rm = TRUE)
  } else {
    NA_real_
  }
  deg_pi_stats <- calc_effect_size_stats(sig_genes$PI)
  
  # Summary statistics
  summary_stats <- list(
    comparison = comparison_name,
    total_genes_tested = nrow(results_df),
    valid_tests = sum(!is.na(results_df$pvalue)),
    significant_genes = storey_result$n_significant,
    pi0 = storey_result$pi0,
    alpha = storey_result$alpha,
    upregulated = sum(results_df$significant & results_df$direction == "UP"),
    downregulated = sum(results_df$significant & results_df$direction == "DOWN"),
    group1_name = bm_result$group1_name,
    group2_name = bm_result$group2_name,
    n_group1 = bm_result$n_group1,
    n_group2 = bm_result$n_group2,
    # PI descriptive statistics (primary)
    all_pi_stats = all_pi_stats,
    deg_pi_stats = deg_pi_stats,
    # Cliff's delta descriptive statistics (derived, for compatibility)
    all_delta_stats = all_delta_stats,
    all_abs_delta_median = all_abs_delta_median,
    deg_delta_stats = deg_delta_stats,
    deg_abs_delta_median = deg_abs_delta_median
  )
  
  return(list(
    results_df = results_df,
    summary_stats = summary_stats
  ))
}

# ============================================================================
# Process each comparison
# ============================================================================

cat("\n--- Processing DEG comparisons ---\n")

thyr_deg_results <- list()

# Define comparisons
comparisons <- list(
  R0_vs_R1 = c("R0", "R1"),
  B0_vs_B1 = c("B0", "B1")
)

for (comp_name in names(comparisons)) {
  groups <- comparisons[[comp_name]]
  
  cat(sprintf("\n=== Comparison: %s ===\n", comp_name))
  
  # Process tumor and normal separately
  for (tissue_type in c("tumor", "normal")) {
    cat(sprintf("\n--- %s %s ---\n", comp_name, tissue_type))
    
    # Load DGEList for this comparison
    comp_tissue <- paste(comp_name, tissue_type, sep = "_")
    dgelist_file <- paste0(paths$processed, "analysis_dgelist_", comp_tissue, ".rds")
    
    if (!file.exists(dgelist_file)) {
      cat(sprintf("  DGEList file not found: %s\n", basename(dgelist_file)))
      cat("  Skipping this comparison\n")
      next
    }
    
    # Load DGEList
    dgelist <- readRDS(dgelist_file)
    cat(sprintf("  Loaded DGEList: %d genes, %d samples\n", 
                nrow(dgelist), ncol(dgelist)))
    
    # Extract information from DGEList
    count_matrix <- dgelist$counts
    sample_groups <- as.character(dgelist$samples$group)
    norm_factors <- dgelist$samples$norm.factors
    gene_info_subset <- dgelist$genes
    
    # Report sample distribution
    group_table <- table(sample_groups)
    cat(sprintf("  Samples: %s=%d, %s=%d\n", 
                groups[1], group_table[groups[1]],
                groups[2], group_table[groups[2]]))
    
    # Report normalization factors
    cat(sprintf("  Norm factors range: [%.3f, %.3f], median=%.3f\n",
                min(norm_factors), max(norm_factors), median(norm_factors)))
    
    # Apply filterByExpr (already filtered, but reapply to be sure)
    keep <- filterByExpr(dgelist, group = dgelist$samples$group)
    dgelist_filtered <- dgelist[keep, , keep.lib.sizes = FALSE]
    
    cat(sprintf("  Genes after filterByExpr: %d\n", nrow(dgelist_filtered)))
    
    # Prepare normalized counts
    normalized_counts <- prepare_normalized_counts(dgelist_filtered)
    
    # Perform Brunner-Munzel test
    bm_result <- perform_brunner_munzel_test(
      normalized_counts,
      as.character(dgelist_filtered$samples$group),
      groups[1],
      groups[2]
    )
    
    # Apply Storey correction
    storey_result <- apply_storey_correction(bm_result$pvalues, CONFIG$ALPHA)
    
    # Create summary
    deg_summary <- create_deg_summary(bm_result, storey_result, comp_tissue, 
                                      dgelist_filtered$genes)
    
    # Store results
    thyr_deg_results[[comp_tissue]] <- list(
      comparison = comp_name,
      tissue = tissue_type,
      groups = groups,
      samples = list(
        group1 = rownames(dgelist_filtered$samples)[dgelist_filtered$samples$group == groups[1]],
        group2 = rownames(dgelist_filtered$samples)[dgelist_filtered$samples$group == groups[2]]
      ),
      bm_result = bm_result,
      storey_result = storey_result,
      deg_summary = deg_summary,
      norm_factors_used = norm_factors,
      analysis_date = Sys.time()
    )
    
    # Print summary
    cat(sprintf("\nResults for %s:\n", comp_tissue))
    cat(sprintf("  Genes tested: %d\n", deg_summary$summary_stats$total_genes_tested))
    cat(sprintf("  Significant DEGs: %d (%.1f%%)\n", 
                deg_summary$summary_stats$significant_genes,
                deg_summary$summary_stats$significant_genes / 
                  deg_summary$summary_stats$total_genes_tested * 100))
    cat(sprintf("  Upregulated: %d\n", deg_summary$summary_stats$upregulated))
    cat(sprintf("  Downregulated: %d\n", deg_summary$summary_stats$downregulated))
    if (!is.na(deg_summary$summary_stats$pi0)) {
      cat(sprintf("  Pi0 estimate: %.3f\n", deg_summary$summary_stats$pi0))
    }
    
    # PI statistics (primary)
    pi_all <- deg_summary$summary_stats$all_pi_stats
    cat(sprintf("  PI (ALL genes, n=%d):\n", pi_all$n))
    cat(sprintf("    5-num: %.3f, %.3f, %.3f, %.3f, %.3f (min, Q1, med, Q3, max)\n", 
                pi_all$min, pi_all$q1, pi_all$median, pi_all$q3, pi_all$max))
    
    # Cliff's delta statistics (derived, for reference)
    es_all <- deg_summary$summary_stats$all_delta_stats
    cat(sprintf("  Cliff's delta (derived, ALL genes):\n"))
    cat(sprintf("    5-num: %.3f, %.3f, %.3f, %.3f, %.3f\n", 
                es_all$min, es_all$q1, es_all$median, es_all$q3, es_all$max))
    cat(sprintf("    Median |d|=%.3f (absolute)\n", 
                deg_summary$summary_stats$all_abs_delta_median))
    
    # Effect size for DEGs (if any)
    if (deg_summary$summary_stats$significant_genes > 0) {
      pi_deg <- deg_summary$summary_stats$deg_pi_stats
      es_deg <- deg_summary$summary_stats$deg_delta_stats
      cat(sprintf("  DEGs: PI median=%.3f, delta median=%.3f, |d| median=%.3f\n", 
                  pi_deg$median, es_deg$median, 
                  deg_summary$summary_stats$deg_abs_delta_median))
    }
  }
}

# ============================================================================
# Consistency analysis between tumor and normal
# ============================================================================

cat("\n--- Evaluating tumor-normal consistency ---\n")

thyr_consistency_results <- list()

for (comp_name in names(comparisons)) {
  tumor_name <- paste(comp_name, "tumor", sep = "_")
  normal_name <- paste(comp_name, "normal", sep = "_")
  
  if (!tumor_name %in% names(thyr_deg_results) || 
      !normal_name %in% names(thyr_deg_results)) {
    cat(sprintf("%s: Missing tumor or normal results\n", comp_name))
    next
  }
  
  cat(sprintf("\n%s consistency analysis:\n", comp_name))
  
  # Get DEG results
  tumor_degs <- thyr_deg_results[[tumor_name]]$deg_summary$results_df
  normal_degs <- thyr_deg_results[[normal_name]]$deg_summary$results_df
  
  # Find common genes
  common_genes <- intersect(tumor_degs$gene_id, normal_degs$gene_id)
  cat(sprintf("  Common genes tested: %d\n", length(common_genes)))
  
  # Extract significant genes
  tumor_sig <- tumor_degs[tumor_degs$significant, ]
  normal_sig <- normal_degs[normal_degs$significant, ]
  
  # Find overlapping significant genes
  overlap_sig <- intersect(tumor_sig$gene_id, normal_sig$gene_id)
  cat(sprintf("  Overlapping significant genes: %d\n", length(overlap_sig)))
  
  if (length(overlap_sig) > 0) {
    # Check direction consistency (simple sign check)
    consistent_genes <- character()
    inconsistent_genes <- character()
    
    for (gene in overlap_sig) {
      tumor_fc <- tumor_sig$log2FC[tumor_sig$gene_id == gene]
      normal_fc <- normal_sig$log2FC[normal_sig$gene_id == gene]
      
      # Check if same direction
      if (sign(tumor_fc) == sign(normal_fc)) {
        consistent_genes <- c(consistent_genes, gene)
      } else {
        inconsistent_genes <- c(inconsistent_genes, gene)
      }
    }
    
    cat(sprintf("  Consistent direction: %d\n", length(consistent_genes)))
    cat(sprintf("  Inconsistent direction: %d\n", length(inconsistent_genes)))
    
    # Store consistency results
    thyr_consistency_results[[comp_name]] <- list(
      common_genes = common_genes,
      tumor_sig_count = nrow(tumor_sig),
      normal_sig_count = nrow(normal_sig),
      overlap_sig = overlap_sig,
      consistent_genes = consistent_genes,
      inconsistent_genes = inconsistent_genes,
      consistency_rate = if(length(overlap_sig) > 0) length(consistent_genes) / length(overlap_sig) else NA
    )
    
    # Show top consistent genes
    if (length(consistent_genes) > 0) {
      cat("\n  Top consistent genes:\n")
      for (i in seq_len(min(5, length(consistent_genes)))) {
        gene <- consistent_genes[i]
        tumor_row <- tumor_sig[tumor_sig$gene_id == gene, ]
        normal_row <- normal_sig[normal_sig$gene_id == gene, ]
        gene_name <- ifelse(is.na(tumor_row$gene_name[1]), gene, tumor_row$gene_name[1])
        cat(sprintf("    %s: tumor FC=%.2f (q=%.3e), normal FC=%.2f (q=%.3e)\n",
                    gene_name, 
                    tumor_row$log2FC[1], tumor_row$qvalue[1],
                    normal_row$log2FC[1], normal_row$qvalue[1]))
      }
    }
  } else {
    thyr_consistency_results[[comp_name]] <- list(
      common_genes = common_genes,
      tumor_sig_count = nrow(tumor_sig),
      normal_sig_count = nrow(normal_sig),
      overlap_sig = character(),
      consistent_genes = character(),
      inconsistent_genes = character(),
      consistency_rate = NA_real_
    )
  }
}

# ============================================================================
# Create overall summary
# ============================================================================

cat("\n--- Creating summary ---\n")

summary_data <- data.frame()

for (comp_tissue in names(thyr_deg_results)) {
  result <- thyr_deg_results[[comp_tissue]]
  summary_stats <- result$deg_summary$summary_stats
  pi_all <- summary_stats$all_pi_stats
  es_all <- summary_stats$all_delta_stats
  pi_deg <- summary_stats$deg_pi_stats
  es_deg <- summary_stats$deg_delta_stats
  
  summary_row <- data.frame(
    comparison = result$comparison,
    tissue = result$tissue,
    n_group1 = summary_stats$n_group1,
    n_group2 = summary_stats$n_group2,
    genes_tested = summary_stats$total_genes_tested,
    degs_total = summary_stats$significant_genes,
    degs_up = summary_stats$upregulated,
    degs_down = summary_stats$downregulated,
    deg_rate = round(summary_stats$significant_genes / 
                       summary_stats$total_genes_tested * 100, 2),
    pi0 = round(summary_stats$pi0, 3),
    # PI statistics (primary)
    all_pi_median = round(pi_all$median, 3),
    deg_pi_median = round(pi_deg$median, 3),
    # Cliff's delta statistics (derived)
    all_delta_min = round(es_all$min, 3),
    all_delta_q1 = round(es_all$q1, 3),
    all_delta_median = round(es_all$median, 3),
    all_delta_q3 = round(es_all$q3, 3),
    all_delta_max = round(es_all$max, 3),
    all_delta_mean = round(es_all$mean, 3),
    all_delta_sd = round(es_all$sd, 3),
    all_abs_delta_median = round(summary_stats$all_abs_delta_median, 3),
    deg_delta_median = round(es_deg$median, 3),
    deg_abs_delta_median = round(summary_stats$deg_abs_delta_median, 3),
    stringsAsFactors = FALSE
  )
  
  summary_data <- rbind(summary_data, summary_row)
}

print(summary_data)

# ============================================================================
# Save results
# ============================================================================

cat("\n--- Saving results ---\n")

# Save main DEG results
deg_output <- list(
  date = Sys.Date(),
  config = CONFIG,
  deg_results = thyr_deg_results,
  consistency_results = thyr_consistency_results,
  summary = summary_data,
  version = "v7.5"
)

saveRDS(deg_output, paste0(paths$processed, "thyr_deg_results.rds"))
cat("  DEG results saved: thyr_deg_results.rds\n")

# Export summary as CSV
write.csv(summary_data, 
          paste0(paths$output, "deg_analysis_summary.csv"),
          row.names = FALSE)
cat("  Summary CSV saved: deg_analysis_summary.csv\n")

# Export individual DEG lists
for (comp_tissue in names(thyr_deg_results)) {
  results_df <- thyr_deg_results[[comp_tissue]]$deg_summary$results_df
  
  # Export all genes
  write.csv(results_df,
            paste0(paths$output, sprintf("deg_results_%s_all.csv", comp_tissue)),
            row.names = FALSE)
  
  # Export significant genes only
  sig_df <- results_df[results_df$significant, ]
  if (nrow(sig_df) > 0) {
    write.csv(sig_df,
              paste0(paths$output, sprintf("deg_results_%s_significant.csv", comp_tissue)),
              row.names = FALSE)
  }
}
cat("  Individual DEG lists exported to output/\n")

# ============================================================================
# Effect size distribution visualization
# ============================================================================

cat("\n--- Generating effect size distribution plots ---\n")

# Create individual histograms for each comparison
effect_size_plots <- list()

for (comp_tissue in names(thyr_deg_results)) {
  result <- thyr_deg_results[[comp_tissue]]
  results_df <- result$deg_summary$results_df
  es_stats <- result$deg_summary$summary_stats$all_delta_stats
  n_degs <- result$deg_summary$summary_stats$significant_genes
  
  # Create histogram (signed delta to preserve direction)
  # binwidth aligned with tie-parity discrete steps (2/(n1*n2)):
  # R0/R1: 8/143 ~= 0.056 (4x of 2/143), B0/B1: 14/260 ~= 0.054 (7x of 2/260)
  bw <- if (grepl("^R0", comp_tissue)) 8/143 else 14/260
  
  p <- ggplot(results_df, aes(x = cliffs_delta)) +
    geom_histogram(binwidth = bw, boundary = 0, fill = "steelblue", color = "white", alpha = 0.7) +
    geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 0.5) +
    geom_vline(xintercept = median(results_df$cliffs_delta, na.rm = TRUE), 
               color = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = sprintf("%s (n=%d genes, %d DEGs)", 
                      comp_tissue, nrow(results_df), n_degs),
      subtitle = sprintf("Median d = %.3f (red), 0 line (black)", 
                         median(results_df$cliffs_delta, na.rm = TRUE)),
      x = "Cliff's delta (d): negative = down, positive = up",
      y = "Count"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.title = element_text(size = 10)
    ) +
    coord_cartesian(xlim = c(-1, 1))
  
  effect_size_plots[[comp_tissue]] <- p
  
  # Display in RStudio
  print(p)
  
  # Save individual plot
  ggsave(
    filename = paste0(paths$output, sprintf("effect_size_hist_%s.pdf", comp_tissue)),
    plot = p, width = 6, height = 4
  )
  cat(sprintf("  Saved: effect_size_hist_%s.pdf\n", comp_tissue))
}

# Create combined figure with unified Y-axis
if (length(effect_size_plots) >= 2) {
  cat("  Creating combined figure with unified Y-axis...\n")
  
  # Fixed Y-axis: range up to 1500, tick marks up to 1500
  y_max <- 1500
  y_breaks <- seq(0, 1500, 500)
  
  # Recreate plots with unified Y-axis for combined figure
  effect_size_plots_unified <- list()
  
  for (comp_tissue in names(thyr_deg_results)) {
    result <- thyr_deg_results[[comp_tissue]]
    results_df <- result$deg_summary$results_df
    n_degs <- result$deg_summary$summary_stats$significant_genes
    
    # binwidth aligned with tie-parity discrete steps (2/(n1*n2)):
    # R0/R1: 8/143 ~= 0.056, B0/B1: 14/260 ~= 0.054
    bw <- if (grepl("^R0", comp_tissue)) 8/143 else 14/260
    
    p <- ggplot(results_df, aes(x = cliffs_delta)) +
      geom_histogram(binwidth = bw, boundary = 0, fill = "steelblue", color = "white", alpha = 0.7) +
      geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 0.5) +
      geom_vline(xintercept = median(results_df$cliffs_delta, na.rm = TRUE), 
                 color = "red", linetype = "dashed", linewidth = 1) +
      labs(
        title = sprintf("%s (n=%d genes, %d DEGs)", 
                        comp_tissue, nrow(results_df), n_degs),
        subtitle = sprintf("Median d = %.3f (red), 0 line (black)", 
                           median(results_df$cliffs_delta, na.rm = TRUE)),
        x = "Cliff's delta (d): negative = down, positive = up",
        y = "Count"
      ) +
      scale_y_continuous(breaks = y_breaks) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        axis.title = element_text(size = 9)
      ) +
      coord_cartesian(xlim = c(-1, 1), ylim = c(0, y_max))
    
    effect_size_plots_unified[[comp_tissue]] <- p
  }
  
  # Arrange combined plot
  if (length(effect_size_plots_unified) >= 4) {
    combined_plot <- gridExtra::grid.arrange(
      effect_size_plots_unified[[1]], effect_size_plots_unified[[2]],
      effect_size_plots_unified[[3]], effect_size_plots_unified[[4]],
      ncol = 2, nrow = 2,
      top = grid::textGrob("Effect Size Distribution (Cliff's delta) by Comparison", 
                           gp = grid::gpar(fontsize = 14, fontface = "bold"))
    )
    ggsave(
      filename = paste0(paths$output, "effect_size_distributions_combined.pdf"),
      plot = combined_plot, width = 12, height = 10
    )
  } else {
    combined_plot <- gridExtra::grid.arrange(
      grobs = effect_size_plots_unified,
      ncol = 2,
      top = grid::textGrob("Effect Size Distribution (Cliff's delta) by Comparison", 
                           gp = grid::gpar(fontsize = 14, fontface = "bold"))
    )
    ggsave(
      filename = paste0(paths$output, "effect_size_distributions_combined.pdf"),
      plot = combined_plot, width = 12, height = 6
    )
  }
  cat(sprintf("  Saved: effect_size_distributions_combined.pdf (Y-axis: 0-%d, ticks to 1000)\n", y_max))
}

# ============================================================================
# Density plot overlay (4 comparisons in 1 panel)
# ============================================================================

cat("\n--- Generating density overlay plot ---\n")

# Combine all delta values into one data.frame
all_deltas <- do.call(rbind, lapply(names(thyr_deg_results), function(comp_tissue) {
  result <- thyr_deg_results[[comp_tissue]]
  deltas <- result$deg_summary$results_df$cliffs_delta
  n_degs <- result$deg_summary$summary_stats$significant_genes
  # Create label with DEG count
  label <- sprintf("%s (%d DEGs)", comp_tissue, n_degs)
  data.frame(
    comparison = comp_tissue,
    label = label,
    delta = deltas,
    stringsAsFactors = FALSE
  )
}))

# Order for legend: tumor first, then normal; within each, R0 then B0
ordered_labels <- c()
for (comp in c("R0_vs_R1_tumor", "B0_vs_B1_tumor", "R0_vs_R1_normal", "B0_vs_B1_normal")) {
  if (comp %in% names(thyr_deg_results)) {
    n <- thyr_deg_results[[comp]]$deg_summary$summary_stats$significant_genes
    ordered_labels <- c(ordered_labels, sprintf("%s (%d DEGs)", comp, n))
  }
}
all_deltas$label <- factor(all_deltas$label, levels = ordered_labels)

# Custom colors: tumor = warm (red/orange), normal = cool (blue/purple)
custom_colors <- c(
  ordered_labels[1],  # R0_vs_R1_tumor - red
  ordered_labels[2],  # B0_vs_B1_tumor - orange
  ordered_labels[3],  # R0_vs_R1_normal - blue
  ordered_labels[4]   # B0_vs_B1_normal - purple
)
names(custom_colors) <- ordered_labels
color_values <- c("#E41A1C", "#FF7F00", "#377EB8", "#984EA3")
names(color_values) <- ordered_labels

# Create density overlay plot
density_plot <- ggplot(all_deltas, aes(x = delta, color = label, fill = label)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 0.5) +
  labs(
    title = "Effect Size Distribution Comparison",
    subtitle = "Kernel density estimation of Cliff's delta across all comparisons",
    x = "Cliff's delta (d): negative = down, positive = up",
    y = "Density",
    color = "Comparison",
    fill = "Comparison"
  ) +
  scale_color_manual(values = color_values) +
  scale_fill_manual(values = color_values) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  ) +
  coord_cartesian(xlim = c(-1, 1))

# Display
print(density_plot)

# Save
ggsave(
  filename = paste0(paths$output, "effect_size_density_overlay.pdf"),
  plot = density_plot, width = 10, height = 6
)
cat("  Saved: effect_size_density_overlay.pdf\n")

# ============================================================================
# Volcano plots (Cliff's delta and log2FC versions)
# ============================================================================

cat("\n--- Generating volcano plots ---\n")

volcano_plots_delta <- list()
volcano_plots_log2fc <- list()

for (comp_tissue in names(thyr_deg_results)) {
  result <- thyr_deg_results[[comp_tissue]]
  results_df <- result$deg_summary$results_df
  n_degs <- result$deg_summary$summary_stats$significant_genes
  
  # Add -log10(qvalue) column, handling NA and 0 values
  # Cap at 10 for visualization
  y_cap <- 10
  results_df$neg_log10_q <- pmin(-log10(pmax(results_df$qvalue, 1e-300)), y_cap)
  
  # Color by significance and direction
  results_df$deg_status <- ifelse(!results_df$significant, "non-DEG",
                                  ifelse(results_df$direction == "UP", "Up", "Down"))
  results_df$deg_status <- factor(results_df$deg_status, levels = c("Up", "Down", "non-DEG"))
  
  # Cliff's delta version
  p_delta <- ggplot(results_df, aes(x = cliffs_delta, y = neg_log10_q, color = deg_status)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray60", linewidth = 0.3) +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "non-DEG" = "gray70")) +
    labs(
      title = sprintf("%s", comp_tissue),
      subtitle = sprintf("%d DEGs (Up: %d, Down: %d)", 
                         n_degs,
                         sum(results_df$deg_status == "Up"),
                         sum(results_df$deg_status == "Down")),
      x = "Cliff's delta",
      y = "-log10(q-value)",
      color = "Status"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "right"
    ) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, y_cap))
  
  volcano_plots_delta[[comp_tissue]] <- p_delta
  
  # log2FC version
  p_log2fc <- ggplot(results_df, aes(x = log2FC, y = neg_log10_q, color = deg_status)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray60", linewidth = 0.3) +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "non-DEG" = "gray70")) +
    labs(
      title = sprintf("%s", comp_tissue),
      subtitle = sprintf("%d DEGs (Up: %d, Down: %d)", 
                         n_degs,
                         sum(results_df$deg_status == "Up"),
                         sum(results_df$deg_status == "Down")),
      x = "log2 Fold Change",
      y = "-log10(q-value)",
      color = "Status"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "right"
    ) +
    coord_cartesian(ylim = c(0, y_cap))
  
  volcano_plots_log2fc[[comp_tissue]] <- p_log2fc
}

# Create combined volcano plots (2x2 for each version)
if (length(volcano_plots_delta) >= 4) {
  # Cliff's delta version
  combined_volcano_delta <- gridExtra::grid.arrange(
    volcano_plots_delta[[1]], volcano_plots_delta[[2]],
    volcano_plots_delta[[3]], volcano_plots_delta[[4]],
    ncol = 2, nrow = 2,
    top = grid::textGrob("Volcano Plot (Cliff's delta vs q-value)", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  ggsave(
    filename = paste0(paths$output, "volcano_cliffs_delta.pdf"),
    plot = combined_volcano_delta, width = 12, height = 10
  )
  cat("  Saved: volcano_cliffs_delta.pdf\n")
  
  # log2FC version
  combined_volcano_log2fc <- gridExtra::grid.arrange(
    volcano_plots_log2fc[[1]], volcano_plots_log2fc[[2]],
    volcano_plots_log2fc[[3]], volcano_plots_log2fc[[4]],
    ncol = 2, nrow = 2,
    top = grid::textGrob("Volcano Plot (log2 Fold Change vs q-value)", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  ggsave(
    filename = paste0(paths$output, "volcano_log2fc.pdf"),
    plot = combined_volcano_log2fc, width = 12, height = 10
  )
  cat("  Saved: volcano_log2fc.pdf\n")
}

# ============================================================================
# Final report
# ============================================================================

cat("\n=== DEG Analysis Complete (v7.5) ===\n")
cat("Configuration:\n")
cat("  Test: Brunner-Munzel\n")
cat("  Correction: Storey method (qvalue, lambda=0.5)\n")
cat("  Effect size: PI (primary), Cliff's delta (derived)\n")
cat("  Significance: q <", CONFIG$ALPHA, "\n")
cat("\nProcessed comparisons:\n")

for (comp_tissue in names(thyr_deg_results)) {
  result <- thyr_deg_results[[comp_tissue]]
  pi_all <- result$deg_summary$summary_stats$all_pi_stats
  es_all <- result$deg_summary$summary_stats$all_delta_stats
  abs_med <- result$deg_summary$summary_stats$all_abs_delta_median
  n_genes <- result$deg_summary$summary_stats$total_genes_tested
  n_degs <- result$deg_summary$summary_stats$significant_genes
  cat(sprintf("  %s: %d DEGs from %d genes\n",
              comp_tissue,
              n_degs,
              n_genes))
  cat(sprintf("    PI: median=%.3f, Cliff's delta: median=%.3f, |d| median=%.3f\n",
              pi_all$median, es_all$median, abs_med))
}

if (length(thyr_consistency_results) > 0) {
  cat("\nConsistency analysis:\n")
  for (comp_name in names(thyr_consistency_results)) {
    cons <- thyr_consistency_results[[comp_name]]
    if (length(cons$consistent_genes) > 0) {
      cat(sprintf("  %s: %d consistent genes (%.0f%% of overlap)\n",
                  comp_name,
                  length(cons$consistent_genes),
                  cons$consistency_rate * 100))
    } else {
      cat(sprintf("  %s: No overlapping DEGs\n", comp_name))
    }
  }
}

cat("\nOutputs:\n")
cat("  Main: thyr_deg_results.rds\n")
cat("  Summary: deg_analysis_summary.csv\n")
cat("  DEG lists: deg_results_*_all.csv, deg_results_*_significant.csv\n")
cat("  Effect size plots: effect_size_hist_*.pdf, effect_size_distributions_combined.pdf\n")

# Highlight key findings
cat("\n=== Key Findings ===\n")

# Find comparison with most DEGs
max_degs_idx <- which.max(summary_data$degs_total)
if (length(max_degs_idx) > 0) {
  max_row <- summary_data[max_degs_idx, ]
  cat(sprintf("Most DEGs: %s_%s with %d genes (%.1f%%)\n",
              max_row$comparison, max_row$tissue,
              max_row$degs_total, max_row$deg_rate))
}

# Check if R0_vs_R1_tumor has results
if ("R0_vs_R1_tumor" %in% names(thyr_deg_results)) {
  r0r1_tumor <- thyr_deg_results[["R0_vs_R1_tumor"]]
  cat(sprintf("\nR0_vs_R1_tumor (primary target): %d DEGs\n",
              r0r1_tumor$deg_summary$summary_stats$significant_genes))
  
  # Show top genes if available
  sig_genes <- r0r1_tumor$deg_summary$results_df[r0r1_tumor$deg_summary$results_df$significant, ]
  if (nrow(sig_genes) >= 5) {
    cat("  Top 5 genes by q-value:\n")
    for (i in 1:5) {
      gene_name <- ifelse(is.na(sig_genes$gene_name[i]), 
                          sig_genes$gene_id[i], 
                          sig_genes$gene_name[i])
      cat(sprintf("    %s: FC=%.2f, PI=%.3f, q=%.3e\n",
                  gene_name, sig_genes$log2FC[i], sig_genes$PI[i], sig_genes$qvalue[i]))
    }
  }
}

cat("\nNext steps:\n")
if (sum(summary_data$degs_total) > 20) {
  cat("  1. Proceed to enrichment analysis (sufficient DEGs)\n")
  cat("  2. Feature selection for biomarker development\n")
} else if (sum(summary_data$degs_total) > 0) {
  cat("  1. Limited DEGs - consider GSEA or direct biomarker selection\n")
  cat("  2. Review individual genes for biological relevance\n")
} else {
  cat("  1. No significant DEGs detected\n")
  cat("  2. Consider adjusting parameters or alternative approaches\n")
}

# Clean up (preserve key objects)
rm(list = setdiff(ls(), c("paths", "thyr_deg_results", 
                          "thyr_consistency_results")))
gc()

cat("\nAnalysis completed successfully!\n")