# 09_deg_analysis.R - DEG Analysis with Brunner-Munzel and Storey Method
# Purpose: Perform differential expression analysis on DEGES-normalized CPM data
# Method: Brunner-Munzel test on CPM values with Storey (qvalue) correction
# Input: analysis_dgelist_*.rds, analysis_cpm_*.rds (from 08_deges_normalization.R)
# Output: thyr_deg_results.rds with DEG lists and enhanced consistency evaluation
# Version: v7.3 - Enhanced consistency analysis (tumor-normal + cross-driver)
# Date: 2025-01-21

source("analysis_v7/setup.R")

cat("\n=== DEG Analysis with Brunner-Munzel + Storey Method (v7.3) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Method: Brunner-Munzel test on normalized CPM (no log)\n")
cat("Correction: Storey method (q < 0.05)\n")
cat("Enhancement: Cross-driver consistency evaluation\n")

# Load packages
suppressPackageStartupMessages({
  library(edgeR)
  library(qvalue)
  library(dplyr)
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
cat("  Input: Normalized CPM (no log transformation)\n")
cat("  No pseudocount added\n")

# ============================================================================
# Helper functions for DEG analysis
# ============================================================================

# Perform Brunner-Munzel test on CPM data
perform_brunner_munzel_test <- function(cpm_data, sample_groups, 
                                        group1_name, group2_name) {
  cat("  Performing Brunner-Munzel tests...\n")
  
  # Identify group indices
  group1_indices <- which(sample_groups == group1_name)
  group2_indices <- which(sample_groups == group2_name)
  
  cat(sprintf("    Testing %s (n=%d) vs %s (n=%d)\n", 
              group1_name, length(group1_indices),
              group2_name, length(group2_indices)))
  
  # Initialize results
  n_genes <- nrow(cpm_data)
  pvalues <- rep(NA_real_, n_genes)
  statistics <- rep(NA_real_, n_genes)
  fold_changes <- rep(NA_real_, n_genes)
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
    
    # Extract CPM values (no log transformation)
    group1_values <- as.numeric(cpm_data[i, group1_indices])
    group2_values <- as.numeric(cpm_data[i, group2_indices])
    
    # Calculate means
    mean1 <- mean(group1_values, na.rm = TRUE)
    mean2 <- mean(group2_values, na.rm = TRUE)
    group1_means[i] <- mean1
    group2_means[i] <- mean2
    
    # Log2 fold change (group2 vs group1) - no pseudocount
    fold_changes[i] <- log2(mean2) - log2(mean1)
    
    # Brunner-Munzel test
    tryCatch({
      # Check for sufficient variation
      if (length(unique(c(group1_values, group2_values))) > 1 && 
          var(c(group1_values, group2_values), na.rm = TRUE) > 0) {
        
        bm_result <- brunnermunzel::brunnermunzel.test(group1_values, group2_values)
        pvalues[i] <- bm_result$p.value
        statistics[i] <- bm_result$statistic
      } else {
        pvalues[i] <- 1.0
        statistics[i] <- 0
      }
    }, error = function(e) {
      pvalues[i] <- 1.0
      statistics[i] <- 0
    })
  }
  
  cat("    Brunner-Munzel testing completed\n")
  
  return(list(
    pvalues = pvalues,
    statistics = statistics,
    fold_changes = fold_changes,
    group1_means = group1_means,
    group2_means = group2_means,
    gene_names = rownames(cpm_data),
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
  
  # Apply Storey method
  tryCatch({
    qvalue_result <- qvalue(pvalues[valid_pvals], pi0.method = "bootstrap")
    
    # Create full qvalue vector
    qvalues <- rep(NA_real_, length(pvalues))
    qvalues[valid_pvals] <- qvalue_result$qvalues
    
    # Identify significant genes (q < 0.05, no FC cutoff)
    significant <- qvalues < alpha & !is.na(qvalues)
    n_significant <- sum(significant)
    
    cat(sprintf("    pi0 estimate: %.3f\n", qvalue_result$pi0))
    cat(sprintf("    Significant genes (q < %.2f): %d (%.1f%%)\n",
                alpha, n_significant, n_significant/sum(valid_pvals)*100))
    
    return(list(
      qvalues = qvalues,
      pi0 = qvalue_result$pi0,
      significant = significant,
      n_significant = n_significant,
      alpha = alpha
    ))
    
  }, error = function(e) {
    cat(sprintf("    Storey method failed: %s\n", e$message))
    cat("    Falling back to Benjamini-Hochberg\n")
    
    # Fallback to BH
    qvalues <- rep(NA_real_, length(pvalues))
    qvalues[valid_pvals] <- p.adjust(pvalues[valid_pvals], method = "BH")
    significant <- qvalues < alpha & !is.na(qvalues)
    n_significant <- sum(significant)
    
    cat(sprintf("    Significant genes (BH q < %.2f): %d\n", alpha, n_significant))
    
    return(list(
      qvalues = qvalues,
      pi0 = NA_real_,
      significant = significant,
      n_significant = n_significant,
      alpha = alpha,
      method = "BH"
    ))
  })
}

# Create DEG results summary
create_deg_summary <- function(bm_result, storey_result, comparison_name, gene_info_df) {
  
  # Create results data frame
  results_df <- data.frame(
    gene_id = bm_result$gene_names,
    pvalue = bm_result$pvalues,
    qvalue = storey_result$qvalues,
    log2FC = bm_result$fold_changes,
    statistic = bm_result$statistics,
    group1_mean_cpm = bm_result$group1_means,
    group2_mean_cpm = bm_result$group2_means,
    significant = storey_result$significant,
    stringsAsFactors = FALSE
  )
  
  # Add fold change magnitude and direction
  results_df$abs_log2FC <- abs(results_df$log2FC)
  results_df$direction <- ifelse(results_df$log2FC > 0, "UP", "DOWN")
  
  # Add gene annotation if available
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
    n_group2 = bm_result$n_group2
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
    
    # Load normalized CPM for this comparison
    comp_tissue <- paste(comp_name, tissue_type, sep = "_")
    cpm_file <- paste0(paths$processed, "analysis_cpm_", comp_tissue, ".rds")
    dgelist_file <- paste0(paths$processed, "analysis_dgelist_", comp_tissue, ".rds")
    
    if (!file.exists(cpm_file) || !file.exists(dgelist_file)) {
      cat(sprintf("  Files not found for %s:\n", comp_tissue))
      cat(sprintf("    CPM file exists: %s\n", file.exists(cpm_file)))
      cat(sprintf("    DGEList file exists: %s\n", file.exists(dgelist_file)))
      cat("  Possible causes:\n")
      cat("    - Insufficient samples after quality filtering\n")
      cat("    - All genes removed by Cook's distance filtering\n")
      cat("    - filterByExpr removed all genes\n")
      cat("  Check 08_deges_normalization.R output for details\n")
      cat("  Skipping this comparison\n")
      next
    }
    
    # Load normalized CPM
    normalized_cpm <- readRDS(cpm_file)
    cat(sprintf("  Loaded CPM: %d genes, %d samples\n", 
                nrow(normalized_cpm), ncol(normalized_cpm)))
    
    # Load DGEList for sample information and gene info
    dgelist <- readRDS(dgelist_file)
    sample_groups <- as.character(dgelist$samples$group)
    gene_info_subset <- dgelist$genes
    
    # Report sample distribution
    group_table <- table(sample_groups)
    cat(sprintf("  Samples: %s=%d, %s=%d\n", 
                groups[1], group_table[groups[1]],
                groups[2], group_table[groups[2]]))
    
    # Perform Brunner-Munzel test on CPM
    bm_result <- perform_brunner_munzel_test(
      normalized_cpm,
      sample_groups,
      groups[1],
      groups[2]
    )
    
    # Apply Storey correction
    storey_result <- apply_storey_correction(bm_result$pvalues, CONFIG$ALPHA)
    
    # Create summary
    deg_summary <- create_deg_summary(bm_result, storey_result, comp_tissue, 
                                      gene_info_subset)
    
    # Store results
    thyr_deg_results[[comp_tissue]] <- list(
      comparison = comp_name,
      tissue = tissue_type,
      groups = groups,
      samples = list(
        group1 = rownames(dgelist$samples)[dgelist$samples$group == groups[1]],
        group2 = rownames(dgelist$samples)[dgelist$samples$group == groups[2]]
      ),
      bm_result = bm_result,
      storey_result = storey_result,
      deg_summary = deg_summary,
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
  }
}

# ============================================================================
# Enhanced consistency analysis
# ============================================================================

cat("\n=== ENHANCED CONSISTENCY ANALYSIS ===\n")

# Initialize comprehensive consistency results
thyr_consistency_results <- list()

# ----------------------------------------------------------------------------
# 1. Tumor-Normal consistency (within same comparison)
# ----------------------------------------------------------------------------

cat("\n--- 1. Tumor-Normal Consistency ---\n")

for (comp_name in names(comparisons)) {
  tumor_name <- paste(comp_name, "tumor", sep = "_")
  normal_name <- paste(comp_name, "normal", sep = "_")
  
  if (!tumor_name %in% names(thyr_deg_results) || 
      !normal_name %in% names(thyr_deg_results)) {
    cat(sprintf("%s: Missing tumor or normal results\n", comp_name))
    next
  }
  
  cat(sprintf("\n%s tumor-normal consistency:\n", comp_name))
  
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
    # Check direction consistency
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
    
    # Store tumor-normal consistency
    thyr_consistency_results[[paste0(comp_name, "_tumor_normal")]] <- list(
      type = "tumor_normal",
      comparison = comp_name,
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
      cat("\n  Top consistent genes (up to 5):\n")
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
    thyr_consistency_results[[paste0(comp_name, "_tumor_normal")]] <- list(
      type = "tumor_normal",
      comparison = comp_name,
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

# ----------------------------------------------------------------------------
# 2. Cross-driver consistency (RET vs BRAF)
# ----------------------------------------------------------------------------

cat("\n--- 2. Cross-Driver Consistency (RET vs BRAF) ---\n")

# Tumor comparison
if ("R0_vs_R1_tumor" %in% names(thyr_deg_results) && 
    "B0_vs_B1_tumor" %in% names(thyr_deg_results)) {
  
  cat("\nTumor cross-driver consistency:\n")
  
  ret_tumor <- thyr_deg_results[["R0_vs_R1_tumor"]]$deg_summary$results_df
  braf_tumor <- thyr_deg_results[["B0_vs_B1_tumor"]]$deg_summary$results_df
  
  # Common genes
  common_genes <- intersect(ret_tumor$gene_id, braf_tumor$gene_id)
  cat(sprintf("  Common genes tested: %d\n", length(common_genes)))
  
  # Significant genes
  ret_sig <- ret_tumor[ret_tumor$significant, ]
  braf_sig <- braf_tumor[braf_tumor$significant, ]
  
  cat(sprintf("  RET significant: %d\n", nrow(ret_sig)))
  cat(sprintf("  BRAF significant: %d\n", nrow(braf_sig)))
  
  # Overlapping significant
  overlap_sig <- intersect(ret_sig$gene_id, braf_sig$gene_id)
  cat(sprintf("  Overlapping significant: %d\n", length(overlap_sig)))
  
  if (length(overlap_sig) > 0) {
    # Direction consistency
    consistent_genes <- character()
    inconsistent_genes <- character()
    
    for (gene in overlap_sig) {
      ret_fc <- ret_sig$log2FC[ret_sig$gene_id == gene]
      braf_fc <- braf_sig$log2FC[braf_sig$gene_id == gene]
      
      if (sign(ret_fc) == sign(braf_fc)) {
        consistent_genes <- c(consistent_genes, gene)
      } else {
        inconsistent_genes <- c(inconsistent_genes, gene)
      }
    }
    
    cat(sprintf("  Consistent direction: %d\n", length(consistent_genes)))
    cat(sprintf("  Inconsistent direction: %d\n", length(inconsistent_genes)))
    
    # Store cross-driver tumor consistency
    thyr_consistency_results[["cross_driver_tumor"]] <- list(
      type = "cross_driver",
      tissue = "tumor",
      ret_sig_count = nrow(ret_sig),
      braf_sig_count = nrow(braf_sig),
      overlap_sig = overlap_sig,
      consistent_genes = consistent_genes,
      inconsistent_genes = inconsistent_genes,
      consistency_rate = length(consistent_genes) / length(overlap_sig)
    )
    
    # Show top consistent cross-driver genes
    if (length(consistent_genes) > 0) {
      cat("\n  Top cross-driver consistent genes (radiation response candidates):\n")
      for (i in seq_len(min(5, length(consistent_genes)))) {
        gene <- consistent_genes[i]
        ret_row <- ret_sig[ret_sig$gene_id == gene, ]
        braf_row <- braf_sig[braf_sig$gene_id == gene, ]
        gene_name <- ifelse(is.na(ret_row$gene_name[1]), gene, ret_row$gene_name[1])
        cat(sprintf("    %s: RET FC=%.2f (q=%.3e), BRAF FC=%.2f (q=%.3e)\n",
                    gene_name,
                    ret_row$log2FC[1], ret_row$qvalue[1],
                    braf_row$log2FC[1], braf_row$qvalue[1]))
      }
    }
  } else {
    thyr_consistency_results[["cross_driver_tumor"]] <- list(
      type = "cross_driver",
      tissue = "tumor",
      ret_sig_count = nrow(ret_sig),
      braf_sig_count = nrow(braf_sig),
      overlap_sig = character(),
      consistent_genes = character(),
      inconsistent_genes = character(),
      consistency_rate = NA_real_
    )
  }
}

# Normal comparison
if ("R0_vs_R1_normal" %in% names(thyr_deg_results) && 
    "B0_vs_B1_normal" %in% names(thyr_deg_results)) {
  
  cat("\nNormal cross-driver consistency:\n")
  
  ret_normal <- thyr_deg_results[["R0_vs_R1_normal"]]$deg_summary$results_df
  braf_normal <- thyr_deg_results[["B0_vs_B1_normal"]]$deg_summary$results_df
  
  # Significant genes
  ret_sig <- ret_normal[ret_normal$significant, ]
  braf_sig <- braf_normal[braf_normal$significant, ]
  
  cat(sprintf("  RET significant: %d\n", nrow(ret_sig)))
  cat(sprintf("  BRAF significant: %d\n", nrow(braf_sig)))
  
  # Overlapping significant
  overlap_sig <- intersect(ret_sig$gene_id, braf_sig$gene_id)
  cat(sprintf("  Overlapping significant: %d\n", length(overlap_sig)))
  
  if (length(overlap_sig) > 0) {
    # Direction consistency
    consistent_genes <- character()
    inconsistent_genes <- character()
    
    for (gene in overlap_sig) {
      ret_fc <- ret_sig$log2FC[ret_sig$gene_id == gene]
      braf_fc <- braf_sig$log2FC[braf_sig$gene_id == gene]
      
      if (sign(ret_fc) == sign(braf_fc)) {
        consistent_genes <- c(consistent_genes, gene)
      } else {
        inconsistent_genes <- c(inconsistent_genes, gene)
      }
    }
    
    cat(sprintf("  Consistent direction: %d\n", length(consistent_genes)))
    cat(sprintf("  Inconsistent direction: %d\n", length(inconsistent_genes)))
    
    # Store cross-driver normal consistency
    thyr_consistency_results[["cross_driver_normal"]] <- list(
      type = "cross_driver",
      tissue = "normal",
      ret_sig_count = nrow(ret_sig),
      braf_sig_count = nrow(braf_sig),
      overlap_sig = overlap_sig,
      consistent_genes = consistent_genes,
      inconsistent_genes = inconsistent_genes,
      consistency_rate = if(length(overlap_sig) > 0) length(consistent_genes) / length(overlap_sig) else NA
    )
  } else {
    thyr_consistency_results[["cross_driver_normal"]] <- list(
      type = "cross_driver",
      tissue = "normal",
      ret_sig_count = nrow(ret_sig),
      braf_sig_count = nrow(braf_sig),
      overlap_sig = character(),
      consistent_genes = character(),
      inconsistent_genes = character(),
      consistency_rate = NA_real_
    )
  }
}

# ----------------------------------------------------------------------------
# 3. Core radiation response genes (4-way consistency)
# ----------------------------------------------------------------------------

cat("\n--- 3. Core Radiation Response Genes (4-way analysis) ---\n")

# Collect all significant gene lists
all_sig_lists <- list()
for (comp_tissue in names(thyr_deg_results)) {
  sig_genes <- thyr_deg_results[[comp_tissue]]$deg_summary$results_df
  sig_genes <- sig_genes[sig_genes$significant, ]
  all_sig_lists[[comp_tissue]] <- sig_genes$gene_id
}

if (length(all_sig_lists) >= 2) {
  # Find genes appearing in multiple comparisons
  all_sig_genes <- unlist(all_sig_lists)
  gene_freq <- table(all_sig_genes)
  
  # Genes appearing in 2+ comparisons
  multi_comp_genes <- names(gene_freq[gene_freq >= 2])
  cat(sprintf("  Genes significant in 2+ comparisons: %d\n", length(multi_comp_genes)))
  
  # Genes appearing in 3+ comparisons
  three_comp_genes <- names(gene_freq[gene_freq >= 3])
  cat(sprintf("  Genes significant in 3+ comparisons: %d\n", length(three_comp_genes)))
  
  # Genes appearing in all 4 comparisons
  four_comp_genes <- names(gene_freq[gene_freq == 4])
  cat(sprintf("  Genes significant in all 4 comparisons: %d\n", length(four_comp_genes)))
  
  # Store multi-comparison results
  thyr_consistency_results[["multi_comparison"]] <- list(
    type = "multi_comparison",
    two_plus = multi_comp_genes,
    three_plus = three_comp_genes,
    four_way = four_comp_genes,
    gene_frequencies = gene_freq[gene_freq >= 2]
  )
  
  # Show top multi-comparison genes
  if (length(three_comp_genes) > 0) {
    cat("\n  Top multi-comparison genes (3+ comparisons):\n")
    for (i in seq_len(min(5, length(three_comp_genes)))) {
      gene <- three_comp_genes[i]
      cat(sprintf("    %s appears in %d comparisons\n", gene, gene_freq[gene]))
      
      # Show details from each comparison
      for (comp_tissue in names(all_sig_lists)) {
        if (gene %in% all_sig_lists[[comp_tissue]]) {
          deg_row <- thyr_deg_results[[comp_tissue]]$deg_summary$results_df
          deg_row <- deg_row[deg_row$gene_id == gene & deg_row$significant, ]
          if (nrow(deg_row) > 0) {
            cat(sprintf("      %s: FC=%.2f, q=%.3e\n",
                        comp_tissue, deg_row$log2FC[1], deg_row$qvalue[1]))
          }
        }
      }
    }
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
  method = "Brunner-Munzel on CPM (no log)",
  deg_results = thyr_deg_results,
  consistency_results = thyr_consistency_results,
  summary = summary_data,
  version = "v7.3_enhanced"
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
    cat(sprintf("  %s: exported %d significant genes\n", comp_tissue, nrow(sig_df)))
  }
}

# Export consistency analysis results
if (length(thyr_consistency_results) > 0) {
  # Create consistency summary
  consistency_summary <- data.frame()
  
  for (cons_name in names(thyr_consistency_results)) {
    cons <- thyr_consistency_results[[cons_name]]
    if (cons$type == "tumor_normal") {
      consistency_summary <- rbind(consistency_summary, data.frame(
        analysis = cons_name,
        type = cons$type,
        comparison = cons$comparison,
        overlap_count = length(cons$overlap_sig),
        consistent_count = length(cons$consistent_genes),
        consistency_rate = cons$consistency_rate,
        stringsAsFactors = FALSE
      ))
    } else if (cons$type == "cross_driver") {
      consistency_summary <- rbind(consistency_summary, data.frame(
        analysis = cons_name,
        type = cons$type,
        comparison = cons$tissue,
        overlap_count = length(cons$overlap_sig),
        consistent_count = length(cons$consistent_genes),
        consistency_rate = cons$consistency_rate,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  write.csv(consistency_summary,
            paste0(paths$output, "deg_consistency_summary.csv"),
            row.names = FALSE)
  cat("  Consistency summary saved: deg_consistency_summary.csv\n")
}

# ============================================================================
# Final report
# ============================================================================

cat("\n=== DEG Analysis Complete (v7.3 Enhanced) ===\n")
cat("Configuration:\n")
cat("  Test: Brunner-Munzel on normalized CPM\n")
cat("  Correction: Storey method (qvalue)\n")
cat("  Significance: q <", CONFIG$ALPHA, "\n")
cat("  No log transformation of CPM\n")
cat("  No pseudocount added\n")
cat("\nProcessed comparisons:\n")

for (comp_tissue in names(thyr_deg_results)) {
  result <- thyr_deg_results[[comp_tissue]]
  cat(sprintf("  %s: %d DEGs from %d genes (%.1f%%)\n",
              comp_tissue,
              result$deg_summary$summary_stats$significant_genes,
              result$deg_summary$summary_stats$total_genes_tested,
              result$deg_summary$summary_stats$significant_genes / 
                result$deg_summary$summary_stats$total_genes_tested * 100))
}

cat("\n=== Enhanced Consistency Analysis Results ===\n")

# Tumor-normal consistency
if ("R0_vs_R1_tumor_normal" %in% names(thyr_consistency_results)) {
  cons <- thyr_consistency_results[["R0_vs_R1_tumor_normal"]]
  if (length(cons$consistent_genes) > 0) {
    cat(sprintf("R0_vs_R1 tumor-normal: %d consistent genes (%.0f%% of overlap)\n",
                length(cons$consistent_genes),
                cons$consistency_rate * 100))
  }
}

if ("B0_vs_B1_tumor_normal" %in% names(thyr_consistency_results)) {
  cons <- thyr_consistency_results[["B0_vs_B1_tumor_normal"]]
  if (length(cons$consistent_genes) > 0) {
    cat(sprintf("B0_vs_B1 tumor-normal: %d consistent genes (%.0f%% of overlap)\n",
                length(cons$consistent_genes),
                cons$consistency_rate * 100))
  }
}

# Cross-driver consistency
if ("cross_driver_tumor" %in% names(thyr_consistency_results)) {
  cons <- thyr_consistency_results[["cross_driver_tumor"]]
  cat(sprintf("\nCross-driver tumor: %d consistent genes\n",
              length(cons$consistent_genes)))
  if (length(cons$consistent_genes) > 0) {
    cat("  → Potential universal radiation response markers\n")
  }
}

if ("cross_driver_normal" %in% names(thyr_consistency_results)) {
  cons <- thyr_consistency_results[["cross_driver_normal"]]
  cat(sprintf("Cross-driver normal: %d consistent genes\n",
              length(cons$consistent_genes)))
}

# Multi-comparison summary
if ("multi_comparison" %in% names(thyr_consistency_results)) {
  multi <- thyr_consistency_results[["multi_comparison"]]
  cat(sprintf("\nMulti-comparison genes:\n"))
  cat(sprintf("  2+ comparisons: %d genes\n", length(multi$two_plus)))
  cat(sprintf("  3+ comparisons: %d genes\n", length(multi$three_plus)))
  cat(sprintf("  All 4 comparisons: %d genes\n", length(multi$four_way)))
  
  if (length(multi$four_way) > 0) {
    cat("  → Highest confidence radiation markers\n")
  }
}

cat("\nOutputs:\n")
cat("  Main: thyr_deg_results.rds\n")
cat("  Summary: deg_analysis_summary.csv\n")
cat("  Consistency: deg_consistency_summary.csv\n")
cat("  DEG lists: deg_results_*_all.csv, deg_results_*_significant.csv\n")

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

# Check if R0_vs_R1_tumor has results (primary target)
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
      cat(sprintf("    %s: FC=%.2f, q=%.3e\n",
                  gene_name, sig_genes$log2FC[i], sig_genes$qvalue[i]))
    }
  }
}

cat("\nNext steps:\n")
if (sum(summary_data$degs_total) > 20) {
  cat("  1. Proceed to enrichment analysis (sufficient DEGs)\n")
  cat("  2. Feature selection focusing on cross-driver consistent genes\n")
  cat("  3. Validate multi-comparison genes as biomarker candidates\n")
} else if (sum(summary_data$degs_total) > 0) {
  cat("  1. Limited DEGs - consider GSEA or direct biomarker selection\n")
  cat("  2. Focus on consistency-validated genes\n")
} else {
  cat("  1. No significant DEGs detected\n")
  cat("  2. Consider adjusting parameters or alternative approaches\n")
}

# Clean up (preserve key objects)
rm(list = setdiff(ls(), c("paths", "thyr_deg_results", 
                          "thyr_consistency_results")))
gc()

cat("\nAnalysis completed successfully!\n")