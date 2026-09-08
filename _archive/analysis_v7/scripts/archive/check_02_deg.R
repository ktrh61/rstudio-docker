# check_02_deg_revised.R
# Purpose: Check DEG analysis results for manuscript (REVISED)
# Section 2: Normalization and Differential Expression Analysis
# Date: 2025-01-26
# Fixed: Proper extraction of DEGES information

# Setup
source("analysis_v7/setup.R")

cat("=====================================\n")
cat("DEG ANALYSIS RESULTS CHECK (REVISED)\n")
cat("=====================================\n\n")

# ============================================================================
# Gene Filtering Information
# ============================================================================

cat("--- Checking initial gene counts ---\n")

# Check for filtered SE
se_filtered_file <- paste0(paths$processed, "thyr_se_strand2_nonzero.rds")
if (file.exists(se_filtered_file)) {
  se_filtered <- readRDS(se_filtered_file)
  cat("□ Genes after strand selection & zero removal:", nrow(se_filtered), "\n")
  rm(se_filtered)  # Free memory
}

# ============================================================================
# DEGES Normalization Information
# ============================================================================

cat("\n--- DEGES-MUREN Normalization ---\n")

deges_file <- paste0(paths$processed, "analysis_deges_results.rds")
if (file.exists(deges_file)) {
  deges_results <- readRDS(deges_file)
  cat("✓ DEGES results loaded\n")
  
  # Check structure
  cat("\nDEGES object structure:\n")
  cat("  Available elements:", names(deges_results), "\n")
  
  # Extract information from different possible locations
  if ("summary" %in% names(deges_results)) {
    cat("\n--- DEGES Summary ---\n")
    if (is.list(deges_results$summary)) {
      # Print summary information
      for (key in names(deges_results$summary)) {
        cat(sprintf("  %s: %s\n", key, deges_results$summary[[key]]))
      }
    }
  }
  
  if ("results" %in% names(deges_results)) {
    cat("\n--- DEGES Results ---\n")
    
    # Check for iterations
    if ("iterations" %in% names(deges_results$results)) {
      cat("□ Number of iterations:", deges_results$results$iterations, "\n")
    }
    
    # Check for filtered genes
    if ("n_genes" %in% names(deges_results$results)) {
      cat("□ Genes after filterByExpr:", deges_results$results$n_genes, "\n")
    }
    
    # Check for outliers
    if ("outliers_removed" %in% names(deges_results$results)) {
      cat("□ Outliers removed (Cook's):", deges_results$results$outliers_removed, "\n")
    }
  }
  
  # Try to find normalization factors
  if ("norm.factors" %in% names(deges_results)) {
    nf <- deges_results$norm.factors
    cat("\n--- Normalization factors ---\n")
    cat("□ Number of samples:", length(nf), "\n")
    cat("□ Range:", sprintf("%.3f - %.3f\n", min(nf), max(nf)))
  }
  
  # Check config for parameters
  if ("config" %in% names(deges_results)) {
    cat("\n--- DEGES Configuration ---\n")
    if ("iteration" %in% names(deges_results$config)) {
      cat("□ Max iterations:", deges_results$config$iteration, "\n")
    }
    if ("FDR" %in% names(deges_results$config)) {
      cat("□ FDR threshold:", deges_results$config$FDR, "\n")
    }
  }
}

# Alternative: Check norm factors file directly
norm_factors_file <- paste0(paths$processed, "analysis_norm_factors.rds")
if (file.exists(norm_factors_file)) {
  norm_factors <- readRDS(norm_factors_file)
  cat("\n--- Normalization Factors (direct) ---\n")
  
  if (is.list(norm_factors)) {
    for (comparison in names(norm_factors)) {
      nf <- norm_factors[[comparison]]
      cat(sprintf("  %s: %d samples, range %.3f-%.3f\n", 
                  comparison, length(nf), min(nf), max(nf)))
    }
  }
}

# ============================================================================
# DEG Results
# ============================================================================

cat("\n--- Loading DEG results ---\n")

deg_file <- paste0(paths$processed, "thyr_deg_results.rds")
if (!file.exists(deg_file)) {
  stop("DEG results file not found")
}

deg_results <- readRDS(deg_file)
cat("✓ DEG results loaded\n")

# Check structure
cat("\nAvailable comparisons:\n")
print(names(deg_results$deg_results))

# Process each comparison
comparisons <- names(deg_results$deg_results)

for (comp in comparisons) {
  cat("\n=====================================\n")
  cat(toupper(gsub("_", " ", comp)), "\n")
  cat("=====================================\n")
  
  comp_data <- deg_results$deg_results[[comp]]
  
  if ("summary_stats" %in% names(comp_data$deg_summary)) {
    stats <- comp_data$deg_summary$summary_stats
    
    cat("\n--- Summary statistics ---\n")
    cat("□ Total genes tested:", stats$total_genes_tested, "\n")
    cat("□ Valid tests:", stats$valid_tests, "\n")
    cat("□ Significant DEGs (q<0.05):", stats$significant_genes, "\n")
    
    if ("upregulated" %in% names(stats)) {
      cat("□ Up-regulated:", stats$upregulated, "\n")
      cat("□ Down-regulated:", stats$downregulated, "\n")
    }
    
    if ("pi0" %in% names(stats)) {
      cat("□ Pi0 estimate:", round(stats$pi0, 3), "\n")
    }
    
    if ("alpha" %in% names(stats)) {
      cat("□ Alpha used:", stats$alpha, "\n")
    }
    
    if ("n_group1" %in% names(stats)) {
      cat("□ Group 1 samples:", stats$n_group1, "\n")
      cat("□ Group 2 samples:", stats$n_group2, "\n")
    }
  }
  
  # Fold change statistics for tumor comparisons
  if (grepl("tumor", comp) && "results_df" %in% names(comp_data$deg_summary)) {
    results <- comp_data$deg_summary$results_df
    sig_results <- results[results$significant, ]
    
    if (nrow(sig_results) > 0) {
      cat("\n--- Fold change statistics (significant DEGs) ---\n")
      cat("□ Max log2FC:", round(max(sig_results$log2FC), 2), "\n")
      cat("□ Min log2FC:", round(min(sig_results$log2FC), 2), "\n")
      cat("□ Median |log2FC|:", round(median(abs(sig_results$log2FC)), 2), "\n")
    }
  }
}

# ============================================================================
# Cross-driver Analysis
# ============================================================================

cat("\n=====================================\n")
cat("CROSS-DRIVER ANALYSIS\n")
cat("=====================================\n")

if ("consistency_results" %in% names(deg_results)) {
  consistency <- deg_results$consistency_results
  
  if ("cross_driver_tumor" %in% names(consistency)) {
    cross_driver <- consistency$cross_driver_tumor
    
    cat("\n--- Cross-driver genes (tumor) ---\n")
    cat("□ Total consistent genes:", length(cross_driver$consistent_genes), "\n")
    
    # Direction analysis
    if ("R0_vs_R1_tumor" %in% names(deg_results$deg_results) && 
        "B0_vs_B1_tumor" %in% names(deg_results$deg_results)) {
      
      r0r1_df <- deg_results$deg_results$R0_vs_R1_tumor$deg_summary$results_df
      b0b1_df <- deg_results$deg_results$B0_vs_B1_tumor$deg_summary$results_df
      
      # Get cross-driver genes
      cross_genes <- cross_driver$consistent_genes
      
      # Check directions
      r0r1_cross <- r0r1_df[r0r1_df$gene_id %in% cross_genes, ]
      b0b1_cross <- b0b1_df[b0b1_df$gene_id %in% cross_genes, ]
      
      # Match genes properly
      merged <- merge(
        r0r1_cross[, c("gene_id", "log2FC")],
        b0b1_cross[, c("gene_id", "log2FC")],
        by = "gene_id",
        suffixes = c("_r", "_b")
      )
      
      both_up <- sum(merged$log2FC_r > 0 & merged$log2FC_b > 0)
      both_down <- sum(merged$log2FC_r < 0 & merged$log2FC_b < 0)
      mixed <- nrow(merged) - both_up - both_down
      
      cat("□ Both up-regulated:", both_up, "\n")
      cat("□ Both down-regulated:", both_down, "\n")
      cat("□ Mixed direction:", mixed, "\n")
    }
  }
}

# ============================================================================
# Tumor-Normal Consistency
# ============================================================================

cat("\n=====================================\n")
cat("TUMOR-NORMAL CONSISTENCY\n")
cat("=====================================\n")

if ("consistency_results" %in% names(deg_results)) {
  
  # R0 vs R1 consistency
  if ("R0_vs_R1_consistency" %in% names(deg_results$consistency_results)) {
    r0r1_consist <- deg_results$consistency_results$R0_vs_R1_consistency
    cat("\n--- R0 vs R1 consistency ---\n")
    cat("□ Tumor-specific DEGs:", r0r1_consist$tumor_only_count, "\n")
    cat("□ Normal-specific DEGs:", r0r1_consist$normal_only_count, "\n")
    cat("□ Consistent DEGs (same direction):", r0r1_consist$consistent_count, "\n")
    cat("□ Opposite direction:", r0r1_consist$opposite_count, "\n")
    
    # Calculate actual overlap
    total_tumor <- r0r1_consist$tumor_only_count + r0r1_consist$consistent_count + r0r1_consist$opposite_count
    overlap_pct <- round((r0r1_consist$consistent_count + r0r1_consist$opposite_count) / total_tumor * 100, 1)
    cat(sprintf("□ Tumor-normal overlap: %.1f%%\n", overlap_pct))
  }
  
  # B0 vs B1 consistency
  if ("B0_vs_B1_consistency" %in% names(deg_results$consistency_results)) {
    b0b1_consist <- deg_results$consistency_results$B0_vs_B1_consistency
    cat("\n--- B0 vs B1 consistency ---\n")
    cat("□ Tumor-specific DEGs:", b0b1_consist$tumor_only_count, "\n")
    cat("□ Normal-specific DEGs:", b0b1_consist$normal_only_count, "\n")
    cat("□ Consistent DEGs (same direction):", b0b1_consist$consistent_count, "\n")
    cat("□ Opposite direction:", b0b1_consist$opposite_count, "\n")
    
    # Calculate actual overlap
    total_tumor <- b0b1_consist$tumor_only_count + b0b1_consist$consistent_count + b0b1_consist$opposite_count
    if (total_tumor > 0) {
      overlap_pct <- round((b0b1_consist$consistent_count + b0b1_consist$opposite_count) / total_tumor * 100, 1)
      cat(sprintf("□ Tumor-normal overlap: %.1f%%\n", overlap_pct))
    }
  }
}

# ============================================================================
# SUMMARY FOR MANUSCRIPT
# ============================================================================

cat("\n=====================================\n")
cat("SUMMARY FOR MANUSCRIPT\n")
cat("=====================================\n")

cat("\nSection 2 key numbers:\n")
cat("----------------------------------\n")

# Gene filtering
if (exists("se_filtered")) {
  cat("Initial genes (after zero removal):", nrow(se_filtered), "\n")
}

# DEGES information
if (exists("deges_results")) {
  if ("results" %in% names(deges_results)) {
    if ("n_genes" %in% names(deges_results$results)) {
      cat("Genes after filterByExpr:", deges_results$results$n_genes, "\n")
    }
    if ("iterations" %in% names(deges_results$results)) {
      cat("DEGES iterations:", deges_results$results$iterations, "\n")
    }
  }
}

# DEG results
if (exists("deg_results")) {
  # R0 vs R1
  if ("R0_vs_R1_tumor" %in% names(deg_results$deg_results)) {
    stats <- deg_results$deg_results$R0_vs_R1_tumor$deg_summary$summary_stats
    cat("\nR0 vs R1 (tumor):\n")
    cat("  Total DEGs:", stats$significant_genes, "\n")
    cat("  Up in R1:", stats$upregulated, "\n")
    cat("  Down in R1:", stats$downregulated, "\n")
    cat("  q-value cutoff:", stats$alpha, "\n")
  }
  
  # B0 vs B1
  if ("B0_vs_B1_tumor" %in% names(deg_results$deg_results)) {
    stats <- deg_results$deg_results$B0_vs_B1_tumor$deg_summary$summary_stats
    cat("\nB0 vs B1 (tumor):\n")
    cat("  Total DEGs:", stats$significant_genes, "\n")
    cat("  Up in B1:", stats$upregulated, "\n")
    cat("  Down in B1:", stats$downregulated, "\n")
  }
  
  # Cross-driver
  if ("cross_driver_tumor" %in% names(deg_results$consistency_results)) {
    cat("\nCross-driver genes:", 
        length(deg_results$consistency_results$cross_driver_tumor$consistent_genes), "\n")
  }
}

cat("\n=== END OF DEG CHECK (REVISED) ===\n")