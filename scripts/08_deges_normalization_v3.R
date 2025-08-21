# ==============================================================================
# REBC-THYR DEGES Normalization Script v3 - Pure Normalization Approach
# 08_deges_normalization_v3.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(edgeR)
library(DESeq2)
library(qvalue)
library(dplyr)
library(MUREN)

cat("Starting DEGES normalization v3 with pure normalization approach...\n")

# ==============================================================================
# 1. Load Data and High-Purity Sample Lists
# ==============================================================================

cat("Loading high-purity sample lists and data...\n")

# Load high-purity sample lists from 06_purity_analysis_v3.R
load("./data/processed/final_high_purity_sample_lists.rda")

# Load SummarizedExperiment object
if (!exists("se_thyr")) {
  load("./data/raw/thyr_data.rda")
  se_thyr <- data
  rm(data)
}

# Get count data
count_data_full <- assay(se_thyr, "stranded_second")
gene_info <- rowData(se_thyr)

cat("Original count data dimensions:", dim(count_data_full), "\n")

# Use the correct variable name from 06
high_purity_sample_lists <- final_high_purity_sample_lists

cat("High-purity sample lists loaded:\n")
for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(high_purity_sample_lists)) {
    n_pairs <- length(high_purity_sample_lists[[group]]$tumor)
    cat(sprintf("  %s: %d pairs\n", group, n_pairs))
  }
}

# ==============================================================================
# 2. Helper Functions
# ==============================================================================

# DESeq2-based Cook's distance outlier detection
deseq2_cook_outlier_detection <- function(count_matrix, sample_groups, cutoff_percentile = 0.99) {
  cat("Performing DESeq2-based Cook's distance outlier detection...\n")
  
  # Create sample metadata
  sample_info <- data.frame(
    sample_id = colnames(count_matrix),
    group = factor(sample_groups),
    stringsAsFactors = FALSE
  )
  rownames(sample_info) <- sample_info$sample_id
  
  # Check minimum group sizes
  group_sizes <- table(sample_groups)
  min_group_size <- min(group_sizes)
  if (min_group_size < 3) {
    cat(sprintf("Warning: Minimum group size is %d. Cook's distance may be less reliable.\n", min_group_size))
  }
  
  cat("Group sizes:\n")
  print(group_sizes)
  
  # Create DESeq2 object
  tryCatch({
    dds <- DESeqDataSetFromMatrix(
      countData = count_matrix,
      colData = sample_info,
      design = ~ group
    )
    
    # Run DESeq2 for Cook's distance calculation
    cat("Running DESeq2 for Cook's distance calculation...\n")
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    
    # Get Cook's distances
    cooks_distances <- assays(dds)[["cooks"]]
    
    # Calculate cutoff
    m <- ncol(count_matrix)  # number of samples
    p <- length(resultsNames(dds))  # number of parameters
    cooks_cutoff <- qf(cutoff_percentile, df1 = p, df2 = m - p)
    
    cat(sprintf("Cook's distance cutoff (%.1f%% percentile): %.3f\n", 
                cutoff_percentile * 100, cooks_cutoff))
    
    # Find maximum Cook's distance for each gene
    max_cooks <- apply(cooks_distances, 1, function(x) {
      if (all(is.na(x))) {
        return(0)  # Return 0 for genes with all NA Cook's distances
      } else {
        return(max(x, na.rm = TRUE))
      }
    })
    
    # Identify outlier genes
    outlier_genes <- names(max_cooks)[max_cooks > cooks_cutoff & max_cooks != -Inf & !is.na(max_cooks)]
    
    cat(sprintf("Genes flagged as Cook's outliers: %d/%d (%.2f%%)\n",
                length(outlier_genes), nrow(count_matrix),
                length(outlier_genes) / nrow(count_matrix) * 100))
    
    # Return results
    return(list(
      outlier_genes = outlier_genes,
      max_cooks = max_cooks,
      cooks_cutoff = cooks_cutoff,
      cooks_distances = cooks_distances,
      outlier_count = length(outlier_genes),
      dds = dds  # Return DESeq2 object for reference
    ))
    
  }, error = function(e) {
    cat(sprintf("Error in Cook's distance calculation: %s\n", e$message))
    cat("Proceeding without Cook's distance filtering...\n")
    return(list(
      outlier_genes = character(0),
      max_cooks = rep(0, nrow(count_matrix)),
      cooks_cutoff = Inf,
      cooks_distances = NULL,
      outlier_count = 0,
      dds = NULL
    ))
  })
}

# DEGES iteration function (pure normalization)
perform_deges_iteration <- function(count_matrix, sample_groups, iteration = 0, max_iterations = 3) {
  cat(sprintf("\n--- DEGES Iteration %d ---\n", iteration))
  
  # Create DGEList
  dgelist <- DGEList(counts = count_matrix, group = factor(sample_groups))
  
  # Apply filterByExpr for low expression filtering
  keep <- filterByExpr(dgelist, group = factor(sample_groups))
  dgelist_filtered <- dgelist[keep, , keep.lib.sizes = FALSE]
  
  cat(sprintf("After filterByExpr: %d genes (from %d)\n", 
              nrow(dgelist_filtered), nrow(count_matrix)))
  
  # Normalization approach depends on iteration
  if (iteration == 0) {
    # First iteration: use calcNormFactors (TMM)
    dgelist_normalized <- calcNormFactors(dgelist_filtered)
  } else {
    # Subsequent iterations: use calcNormFactors
    dgelist_normalized <- calcNormFactors(dgelist_filtered)
  }
  
  # MUREN normalization
  cat("Applying MUREN normalization...\n")
  tryCatch({
    muren_coeff <- MUREN::muren_norm(
      data.frame(dgelist_normalized$counts),
      workers = max(1, parallel::detectCores() - 1),
      res_return = "scaling_coeff"
    )
    
    # Apply MUREN coefficients
    dgelist_muren <- dgelist_normalized
    dgelist_muren$samples$norm.factors <- muren_coeff
    
  }, error = function(e) {
    cat(sprintf("MUREN normalization failed: %s\n", e$message))
    cat("Using TMM normalization only...\n")
    dgelist_muren <- dgelist_normalized
    muren_coeff <- dgelist_normalized$samples$norm.factors
  })
  
  # Create design matrix
  design <- model.matrix(~ factor(sample_groups))
  
  # GLM fitting and testing for DEGES
  cat("Performing GLM analysis for DEGES...\n")
  dgelist_estimated <- estimateDisp(dgelist_muren, design)
  dgelist_fitted <- glmFit(dgelist_estimated, design)
  dgelist_tested <- glmLRT(dgelist_fitted, coef = 2)
  
  # Get p-values and apply Storey method
  pvalues <- dgelist_tested$table$PValue
  
  cat("Applying Storey method for multiple testing correction...\n")
  tryCatch({
    qvalue_result <- qvalue(pvalues, pi0.method = "bootstrap")
    qvalues <- qvalue_result$qvalues
    pi0 <- qvalue_result$pi0
    cat(sprintf("Estimated pi0 (proportion of non-DE genes): %.3f\n", pi0))
  }, error = function(e) {
    cat(sprintf("Error in qvalue calculation: %s\n", e$message))
    cat("Using Benjamini-Hochberg method as fallback...\n")
    qvalues <- p.adjust(pvalues, method = "BH")
    pi0 <- NA
  })
  
  # Determine potential DEG threshold
  q_threshold <- 0.10
  potential_deg_candidates <- qvalues < q_threshold
  potential_deg_count <- sum(potential_deg_candidates, na.rm = TRUE)
  
  cat(sprintf("Genes with q-value < %.2f: %d/%d (%.2f%%)\n",
              q_threshold, potential_deg_count, length(qvalues),
              potential_deg_count / length(qvalues) * 100))
  
  # Apply DEGES criteria for potential DEG identification
  if (potential_deg_count / length(qvalues) > 0.05) {
    # Use q-value threshold
    non_deg_indices <- qvalues >= q_threshold
    cat("Using q-value >= 0.10 criterion for non-DEG selection\n")
  } else {
    # Use 5% percentile
    percentile_threshold <- sort(qvalues)[round(length(qvalues) * 0.05) + 1]
    non_deg_indices <- qvalues >= percentile_threshold
    cat(sprintf("Using 5%% percentile criterion (q >= %.4f) for non-DEG selection\n", percentile_threshold))
  }
  
  # Extract non-DEG counts for next iteration
  non_deg_counts <- dgelist_normalized$counts[non_deg_indices, ]
  
  cat(sprintf("Non-DEG genes selected for next iteration: %d/%d\n",
              nrow(non_deg_counts), length(qvalues)))
  
  # Check convergence (no potential DEGs remaining)
  converged <- (potential_deg_count == 0) || (nrow(non_deg_counts) == length(qvalues))
  
  # Return results
  return(list(
    dgelist_normalized = dgelist_normalized,
    dgelist_muren = dgelist_muren,
    muren_coeff = muren_coeff,
    pvalues = pvalues,
    qvalues = qvalues,
    pi0 = pi0,
    non_deg_counts = non_deg_counts,
    potential_deg_count = potential_deg_count,
    converged = converged,
    gene_names = rownames(dgelist_normalized),
    iteration = iteration
  ))
}

# Complete DEGES workflow
run_deges_workflow <- function(count_matrix, sample_groups, max_iterations = 3) {
  cat("Starting DEGES workflow...\n")
  
  iteration_results <- list()
  current_counts <- count_matrix
  iteration <- 0
  converged <- FALSE
  
  while (iteration <= max_iterations && !converged) {
    # Run DEGES iteration
    result <- perform_deges_iteration(current_counts, sample_groups, iteration, max_iterations)
    iteration_results[[paste0("iteration_", iteration)]] <- result
    
    # Check convergence
    converged <- result$converged
    
    if (converged) {
      cat(sprintf("DEGES converged at iteration %d\n", iteration))
      break
    } else if (iteration == max_iterations) {
      cat(sprintf("DEGES reached maximum iterations (%d)\n", max_iterations))
      break
    } else {
      # Prepare for next iteration
      current_counts <- result$non_deg_counts
      iteration <- iteration + 1
    }
  }
  
  # Get final normalization coefficients
  final_iteration <- iteration_results[[length(iteration_results)]]
  final_muren_coeff <- final_iteration$muren_coeff
  
  # Apply final normalization to original count matrix (excluding outliers)
  cat("Applying final MUREN coefficients to all genes...\n")
  final_normalized_counts <- t(t(count_matrix) / final_muren_coeff)
  
  # Summary
  cat("\n=== DEGES Workflow Summary ===\n")
  cat(sprintf("Total iterations completed: %d\n", iteration + 1))
  cat(sprintf("Convergence status: %s\n", ifelse(converged, "Converged", "Max iterations reached")))
  cat(sprintf("Final normalization applied to %d genes\n", nrow(final_normalized_counts)))
  
  return(list(
    iteration_results = iteration_results,
    final_muren_coeff = final_muren_coeff,
    final_normalized_counts = final_normalized_counts,
    total_iterations = iteration + 1,
    converged = converged
  ))
}

# ==============================================================================
# 3. Define Analysis Comparisons
# ==============================================================================

cat("Setting up analysis comparisons...\n")

# Define comparisons based on available high-purity groups
comparisons <- list()

# Check which groups have sufficient samples
viable_groups <- c()
min_samples <- 5

for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(high_purity_sample_lists)) {
    n_pairs <- length(high_purity_sample_lists[[group]]$tumor)
    if (n_pairs >= min_samples) {
      viable_groups <- c(viable_groups, group)
      cat(sprintf("  %s: %d pairs (viable)\n", group, n_pairs))
    } else {
      cat(sprintf("  %s: %d pairs (insufficient)\n", group, n_pairs))
    }
  }
}

# Define comparisons based on viable groups
if ("R0" %in% viable_groups && "R1" %in% viable_groups) {
  comparisons[["R0_vs_R1_tumor"]] <- list(
    group1 = "R0", group2 = "R1", tissue = "tumor",
    description = "RET Tumor: Unexposed vs RadHigh (PRIMARY)",
    priority = 1
  )
  comparisons[["R0_vs_R1_normal"]] <- list(
    group1 = "R0", group2 = "R1", tissue = "normal", 
    description = "RET Normal: Unexposed vs RadHigh",
    priority = 2
  )
}

if ("B0" %in% viable_groups && "B1" %in% viable_groups) {
  comparisons[["B0_vs_B1_tumor"]] <- list(
    group1 = "B0", group2 = "B1", tissue = "tumor",
    description = "BRAF Tumor: Unexposed vs RadHigh",
    priority = 3
  )
  comparisons[["B0_vs_B1_normal"]] <- list(
    group1 = "B0", group2 = "B1", tissue = "normal",
    description = "BRAF Normal: Unexposed vs RadHigh",
    priority = 4
  )
}

cat("Defined comparisons:\n")
for (comp_name in names(comparisons)) {
  comp <- comparisons[[comp_name]]
  cat(sprintf("  %s: %s (Priority %d)\n", comp_name, comp$description, comp$priority))
}

if (length(comparisons) == 0) {
  stop("No viable comparisons found. Check sample counts and filtering parameters.")
}

# ==============================================================================
# 4. Process Each Comparison
# ==============================================================================

cat("\n==============================================\n")
cat("DEGES Normalization for Each Comparison\n")
cat("==============================================\n")

# Initialize results storage
deges_results <- list()

for (comp_name in names(comparisons)) {
  cat(sprintf("\n### Processing %s ###\n", comp_name))
  
  comp_info <- comparisons[[comp_name]]
  group1 <- comp_info$group1
  group2 <- comp_info$group2
  tissue <- comp_info$tissue
  
  # Get samples for this comparison
  samples1 <- high_purity_sample_lists[[group1]][[tissue]]
  samples2 <- high_purity_sample_lists[[group2]][[tissue]]
  
  cat(sprintf("%s %s: %d samples\n", group1, tissue, length(samples1)))
  cat(sprintf("%s %s: %d samples\n", group2, tissue, length(samples2)))
  
  if (length(samples1) < 3 || length(samples2) < 3) {
    cat(sprintf("Insufficient samples for comparison (minimum 3 per group), skipping...\n"))
    next
  }
  
  # Combine samples and create group labels
  all_samples <- c(samples1, samples2)
  sample_groups <- c(rep(group1, length(samples1)), rep(group2, length(samples2)))
  
  # Check sample availability in count data
  available_samples <- all_samples[all_samples %in% colnames(count_data_full)]
  if (length(available_samples) < length(all_samples)) {
    missing_count <- length(all_samples) - length(available_samples)
    cat(sprintf("Warning: %d samples missing from count data\n", missing_count))
  }
  
  if (length(available_samples) < 6) {
    cat(sprintf("Insufficient available samples (%d), skipping...\n", length(available_samples)))
    next
  }
  
  # Update group labels to match available samples
  available_groups <- sample_groups[all_samples %in% colnames(count_data_full)]
  
  # Extract count data for this comparison
  comparison_counts <- count_data_full[, available_samples]
  
  # Filter to protein-coding genes
  protein_coding <- gene_info$gene_type == "protein_coding"
  comparison_counts_pc <- comparison_counts[protein_coding, ]
  
  cat(sprintf("Filtered to protein-coding genes: %d -> %d\n", 
              nrow(comparison_counts), nrow(comparison_counts_pc)))
  
  # Step 1: Cook's distance outlier detection using DESeq2
  cook_result <- deseq2_cook_outlier_detection(comparison_counts_pc, available_groups)
  
  # Remove Cook's outlier genes
  if (cook_result$outlier_count > 0) {
    clean_genes <- setdiff(rownames(comparison_counts_pc), cook_result$outlier_genes)
    comparison_counts_clean <- comparison_counts_pc[clean_genes, ]
    cat(sprintf("After Cook's outlier removal: %d genes (removed %d outliers)\n",
                nrow(comparison_counts_clean), cook_result$outlier_count))
  } else {
    comparison_counts_clean <- comparison_counts_pc
    cat("No Cook's outliers detected, proceeding with all genes\n")
  }
  
  # Step 2: Run DEGES workflow
  deges_result <- run_deges_workflow(comparison_counts_clean, available_groups)
  
  # Store comprehensive results
  deges_results[[comp_name]] <- list(
    comparison_info = comp_info,
    samples = available_samples,
    groups = available_groups,
    cook_result = cook_result,
    deges_result = deges_result,
    final_normalized_counts = deges_result$final_normalized_counts,
    original_gene_count = nrow(comparison_counts),
    protein_coding_gene_count = nrow(comparison_counts_pc),
    final_gene_count = nrow(deges_result$final_normalized_counts),
    analysis_date = Sys.time()
  )
  
  cat(sprintf("%s: DEGES normalization completed successfully\n", comp_name))
}

# ==============================================================================
# 5. Summary and Quality Control
# ==============================================================================

cat("\n==============================================\n")
cat("DEGES Normalization v3 Summary\n")
cat("==============================================\n")

# Create summary table
summary_table <- data.frame(
  Comparison = character(0),
  Group1 = character(0),
  Group2 = character(0),
  Tissue = character(0),
  N1 = integer(0),
  N2 = integer(0),
  Original_Genes = integer(0),
  Cook_Outliers = integer(0),
  Final_Genes = integer(0),
  DEGES_Iterations = integer(0),
  Converged = logical(0),
  stringsAsFactors = FALSE
)

for (comp_name in names(deges_results)) {
  result <- deges_results[[comp_name]]
  comp_info <- result$comparison_info
  
  group_counts <- table(result$groups)
  
  summary_table <- rbind(summary_table, data.frame(
    Comparison = comp_name,
    Group1 = comp_info$group1,
    Group2 = comp_info$group2,
    Tissue = comp_info$tissue,
    N1 = as.integer(group_counts[comp_info$group1]),
    N2 = as.integer(group_counts[comp_info$group2]),
    Original_Genes = result$original_gene_count,
    Cook_Outliers = result$cook_result$outlier_count,
    Final_Genes = result$final_gene_count,
    DEGES_Iterations = result$deges_result$total_iterations,
    Converged = result$deges_result$converged,
    stringsAsFactors = FALSE
  ))
}

print(summary_table)

# Quality assessment
cat("\nQuality Assessment:\n")
total_comparisons <- nrow(summary_table)
converged_comparisons <- sum(summary_table$Converged)
cat(sprintf("Successfully processed: %d comparisons\n", total_comparisons))
cat(sprintf("DEGES converged: %d/%d comparisons\n", converged_comparisons, total_comparisons))

if (total_comparisons > 0) {
  avg_iterations <- mean(summary_table$DEGES_Iterations)
  avg_cook_outliers <- mean(summary_table$Cook_Outliers)
  avg_gene_retention <- mean(summary_table$Final_Genes / summary_table$Original_Genes)
  
  cat(sprintf("Average DEGES iterations: %.1f\n", avg_iterations))
  cat(sprintf("Average Cook's outliers per comparison: %.1f\n", avg_cook_outliers))
  cat(sprintf("Average gene retention rate: %.1f%%\n", avg_gene_retention * 100))
} else {
  cat("No comparisons were successfully processed\n")
}

# ==============================================================================
# 6. Save Results
# ==============================================================================

cat("Saving DEGES normalization v3 results...\n")

# Create processed directory if it doesn't exist
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save comprehensive results
deges_normalization_v3_results <- list(
  deges_results = deges_results,
  summary_table = summary_table,
  comparisons = comparisons,
  viable_groups = viable_groups,
  high_purity_sample_lists = high_purity_sample_lists,
  total_comparisons_processed = total_comparisons,
  analysis_date = Sys.time(),
  analysis_version = "v3_high_purity_pure_normalization"
)

save(deges_normalization_v3_results, file = "./data/processed/deges_normalization_v3_results.rda")

cat("Results saved to ./data/processed/deges_normalization_v3_results.rda\n")

# ==============================================================================
# 7. Final Summary and Next Steps
# ==============================================================================

cat("\n==============================================\n")
cat("DEGES Normalization v3 Final Summary\n")
cat("==============================================\n")

cat("Pure Normalization Approach Completed Successfully!\n\n")

cat("Key achievements:\n")
cat("✅ High-purity sample integration from 06_purity_analysis_v3\n")
cat("✅ DESeq2-based Cook's distance outlier detection\n")
cat("✅ DEGES iterative normalization with convergence detection\n")
cat("✅ MUREN scaling factor integration\n")
cat("✅ No contamDE dependency - pure normalization focus\n")

if (total_comparisons > 0) {
  primary_comparison <- "R0_vs_R1_tumor"
  if (primary_comparison %in% names(deges_results)) {
    primary_result <- deges_results[[primary_comparison]]
    cat(sprintf("\nPrimary comparison (%s) results:\n", primary_comparison))
    cat(sprintf("  Final genes: %d\n", primary_result$final_gene_count))
    cat(sprintf("  DEGES iterations: %d\n", primary_result$deges_result$total_iterations))
    cat(sprintf("  Converged: %s\n", ifelse(primary_result$deges_result$converged, "Yes", "No")))
  }
}

cat("\nNormalized data structure:\n")
cat("  Each comparison contains:\n")
cat("    $final_normalized_counts: Ready for statistical testing\n")
cat("    $cook_result: Cook's distance outlier information\n")
cat("    $deges_result: Complete DEGES iteration history\n")

cat("\nNext steps:\n")
cat("1. Review normalization quality (iteration convergence, outlier counts)\n")
cat("2. Proceed to statistical testing (09_deg_analysis_v3.R)\n")
cat("3. Apply Brunner-Munzel test + Storey method to normalized data\n")
cat("4. Target 321 consistent gene discovery\n")

cat("\nDEGES normalization v3 completed successfully!\n")
cat("==============================================\n")

