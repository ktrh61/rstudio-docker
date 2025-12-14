# 07_deges_normalization.R - DEGES Iterative Normalization
# Purpose: Apply DEGES-MUREN normalization to high-purity paired samples
# Method: filterByExpr -> Cook's distance -> MUREN (LTS) + Brunner-Munzel iteration
# Input: thyr_case_master_stage2_filtered, thyr_se_strand2_nonzero  
# Output: Normalized CPM values and DGEList objects with DEGES-MUREN factors
# Version: v7.9 - Separated MA plot generation from DEGES processing
#                 MA plots can now be re-generated independently
#                 Added memory cleanup between processing and visualization
# Date: 2025-12-14

source("analysis_v7/setup.R")

cat("\n=== DEGES Normalization (v7.9) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Method: filterByExpr -> Cook's -> DEGES-MUREN (Brunner-Munzel) -> CPM output\n")
cat("Groups: R0/R1/B0/B1 high-purity pairs only\n")

# Load packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(DESeq2)
  library(brunnermunzel)
  library(qvalue)
  library(dplyr)
})

# Load utility functions
source("utils/utils_improved.R")
source("utils/norm_improved.R")

# ============================================================================
# Configuration
# ============================================================================

CONFIG <- list(
  # DEGES parameters
  MAX_ITERATIONS = 3,         # Maximum DEGES iterations
  CONVERGENCE_THRESHOLD = 0.05, # 5% DEG threshold for floor processing
  
  # Cook's distance
  COOKS_QUANTILE = 0.99,     # Cook's distance quantile threshold
  
  # MUREN parameters
  MUREN_METHOD = "lts",       # LTS robust regression
  MUREN_WORKERS = 3,          # Fixed 3 cores for MUREN
  
  # Thread control
  VERBOSE = TRUE              # Verbose output
)

cat("\nConfiguration:\n")
cat("  Max iterations:", CONFIG$MAX_ITERATIONS, "\n")
cat("  Floor threshold:", sprintf("%.0f%%", CONFIG$CONVERGENCE_THRESHOLD * 100), "\n")
cat("  Cook's quantile:", sprintf("%.0f%%", CONFIG$COOKS_QUANTILE * 100), "\n")
cat("  MUREN method:", CONFIG$MUREN_METHOD, "\n")
cat("  MUREN workers:", CONFIG$MUREN_WORKERS, "\n")
cat("  DEG screening: Brunner-Munzel test + Storey q-value (lambda=0.5)\n")

# Thread control
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  cat("  BLAS/OMP threads: 1\n")
}

deges_processing_log <- list(
  date = Sys.Date(),
  config = CONFIG,
  comparisons = list()
)

# ============================================================================
# Load data
# ============================================================================

cat("\n--- Loading data ---\n")

# Load SummarizedExperiment
if (exists("thyr_se_strand2_nonzero")) {
  cat("Using existing thyr_se_strand2_nonzero\n")
  se <- thyr_se_strand2_nonzero
} else {
  se_path <- paste0(paths$processed, "thyr_se_strand2_nonzero.rds")
  if (!file.exists(se_path)) {
    stop("thyr_se_strand2_nonzero.rds not found. Please run 03_prepare_counts.R first.")
  }
  cat("Loading SE from file...\n")
  se <- readRDS(se_path)
}
cat("SE dimensions:", format(nrow(se), big.mark=","), "genes x ", 
    format(ncol(se), big.mark=","), "samples\n")

# Load case_master_stage2_filtered
if (exists("thyr_case_master_stage2_filtered")) {
  cat("Using existing thyr_case_master_stage2_filtered\n")
  case_master <- thyr_case_master_stage2_filtered
} else {
  case_master_path <- paste0(paths$processed, "thyr_case_master_stage2_filtered.rds")
  if (!file.exists(case_master_path)) {
    stop("thyr_case_master_stage2_filtered.rds not found. Please run 06_purity_analysis.R first.")
  }
  cat("Loading case_master from file...\n")
  case_master <- readRDS(case_master_path)
}
cat("Case master loaded:", nrow(case_master), "cases\n")

# Get sample metadata
sample_metadata <- as.data.frame(colData(se))

# ============================================================================
# Filter to high-purity clean cases
# ============================================================================

cat("\n--- Filtering to high-purity clean cases ---\n")

# Apply all quality filters
high_purity_cases <- case_master %>%
  filter(
    has_outlier_tumor == 0 &    # No tumor outliers
      has_outlier_normal == 0 &   # No normal outliers  
      low_purity == 0              # High tumor purity
  )

cat("High-purity clean cases:", nrow(high_purity_cases), "/", nrow(case_master), 
    sprintf("(%.1f%%)\n", nrow(high_purity_cases)/nrow(case_master)*100))

# Group distribution
high_purity_summary <- high_purity_cases %>%
  group_by(group) %>%
  summarise(n = n(), .groups = "drop")
print(high_purity_summary)

# ============================================================================
# Helper function: Extract paired samples (from 05-07 pattern)
# ============================================================================

extract_paired_samples <- function(group_name, case_df, sample_meta) {
  # Get cases for this group
  group_cases <- case_df$case_submitter_id[case_df$group == group_name]
  
  if (length(group_cases) == 0) {
    return(list(tumor = character(0), normal = character(0), cases = character(0)))
  }
  
  # Initialize vectors to maintain pairing
  tumor_ids <- character(length(group_cases))
  normal_ids <- character(length(group_cases))
  
  for (i in seq_along(group_cases)) {
    case_id <- group_cases[i]
    
    # Find tumor sample (Primary Tumor, _merged)
    tumor_sample <- sample_meta$sample_submitter_id[
      sample_meta$case_submitter_id == case_id &
        sample_meta$sample_type == "Primary Tumor" &
        grepl("_merged", sample_meta$sample_submitter_id)
    ]
    
    # Find normal sample (Solid Tissue Normal, _merged)
    normal_sample <- sample_meta$sample_submitter_id[
      sample_meta$case_submitter_id == case_id &
        sample_meta$sample_type == "Solid Tissue Normal" &
        grepl("_merged", sample_meta$sample_submitter_id)
    ]
    
    # Take first if multiple (shouldn't happen with _merged)
    tumor_ids[i] <- tumor_sample[1]
    normal_ids[i] <- normal_sample[1]
  }
  
  # Remove any NA (shouldn't happen with high-purity cases)
  valid <- !is.na(tumor_ids) & !is.na(normal_ids)
  
  return(list(
    tumor = tumor_ids[valid],
    normal = normal_ids[valid],
    cases = group_cases[valid]
  ))
}

# ============================================================================
# Cook's distance outlier detection
# ============================================================================

detect_cook_outliers <- function(count_matrix, sample_groups, quantile_cutoff = 0.99) {
  cat("  Cook's distance outlier detection...\n")
  
  # Create DESeq2 dataset
  coldata <- data.frame(
    group = factor(sample_groups),
    row.names = colnames(count_matrix)
  )
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = ~ group
  )
  
  # Estimate size factors and dispersions
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  
  # Get Cook's distances
  cooks <- assays(dds)[["cooks"]]
  
  # F-distribution based cutoff (DESeq2 default)
  n_samples <- ncol(count_matrix)
  n_params <- ncol(model.matrix(~ sample_groups))
  f_cutoff <- qf(quantile_cutoff, df1 = n_params, df2 = n_samples - n_params)
  
  # Determine outliers
  max_cooks <- apply(cooks, 1, function(x) {
    if (all(is.na(x))) {
      return(NA)
    } else {
      return(max(x, na.rm = TRUE))
    }
  })
  
  outliers <- !is.na(max_cooks) & max_cooks > f_cutoff
  
  cat(sprintf("    F-distribution cutoff: %.3f\n", f_cutoff))
  cat(sprintf("    Outlier genes: %d (%.2f%%)\n", 
              sum(outliers), sum(outliers) / length(outliers) * 100))
  
  outlier_gene_ids <- rownames(count_matrix)[outliers]
  
  return(list(
    outliers = outliers,
    threshold = f_cutoff,
    outlier_count = sum(outliers),
    outlier_gene_ids = outlier_gene_ids
  ))
}

# ============================================================================
# Brunner-Munzel test for DEG screening
# ============================================================================

perform_bm_test <- function(cpm_matrix, sample_groups, group_levels) {
  # Perform Brunner-Munzel test for each gene
  # Input: normalized CPM matrix (genes x samples)
  # Returns: vector of p-values
  
  group1_idx <- which(sample_groups == group_levels[1])
  group2_idx <- which(sample_groups == group_levels[2])
  
  n_genes <- nrow(cpm_matrix)
  pvalues <- numeric(n_genes)
  
  for (i in seq_len(n_genes)) {
    x <- cpm_matrix[i, group1_idx]
    y <- cpm_matrix[i, group2_idx]
    
    if (length(unique(c(x, y))) > 1 && var(c(x, y), na.rm = TRUE) > 0) {
      result <- brunnermunzel::brunnermunzel.test(x, y)
      pvalues[i] <- result$p.value
    } else {
      pvalues[i] <- 1.0
    }
  }
  
  return(pvalues)
}

# ============================================================================
# DEGES iteration function (Brunner-Munzel-based)
# ============================================================================

perform_deges_iteration <- function(count_matrix, sample_groups, iteration = 0, 
                                    config = CONFIG) {
  cat(sprintf("    Iteration %d:\n", iteration))
  
  group_levels <- levels(factor(sample_groups))
  
  # Create DGEList (no additional filtering here - already filtered)
  dgelist <- DGEList(counts = count_matrix, group = factor(sample_groups))
  
  # MUREN normalization
  cat("      Applying MUREN normalization...\n")
  
  muren_coeff <- muren_norm(
    dgelist$counts,
    refs = "saturated",
    pairwise_method = config$MUREN_METHOD,
    single_param = TRUE,
    res_return = "scaling_coeff",
    workers = config$MUREN_WORKERS
  )
  
  # Validate coefficients
  if (any(!is.finite(muren_coeff)) || any(muren_coeff <= 0)) {
    stop("MUREN scaling coefficients contain non-finite or non-positive values.")
  }
  
  # Apply MUREN factors
  dgelist$samples$norm.factors <- muren_coeff / mean(muren_coeff)
  
  # Calculate normalized CPM for Brunner-Munzel test
  normalized_cpm <- cpm(dgelist, normalized.lib.sizes = TRUE, 
                        prior.count = 0, log = FALSE)
  
  # Perform Brunner-Munzel test
  cat("      Performing Brunner-Munzel tests...\n")
  pvalues <- perform_bm_test(normalized_cpm, sample_groups, group_levels)
  
  # Apply Storey q-value (lambda = 0.5 fixed)
  qval_result <- qvalue(pvalues, lambda = 0.5)
  qvalues <- qval_result$qvalues
  pi0_estimate <- qval_result$pi0
  
  cat(sprintf("      Estimated pi0: %.3f\n", pi0_estimate))
  
  # Count potential DEGs
  deg_count <- sum(qvalues < 0.10, na.rm = TRUE)
  deg_proportion <- deg_count / length(qvalues)
  
  cat(sprintf("      DEGs: %d (%.2f%%)\n", deg_count, deg_proportion * 100))
  
  # Determine potential DEGs and non-DEGs
  n_genes <- length(qvalues)
  exclusion_method <- "none"  # Track which method was used
  
  if (deg_proportion > config$CONVERGENCE_THRESHOLD) {
    # More than 5% DEGs: use q < 0.10 as potential DEGs
    potential_deg_indices <- which(qvalues < 0.10)
    non_deg_indices <- which(qvalues >= 0.10)
    exclusion_method <- "qvalue"
    
  } else {
    # 5% or fewer DEGs: apply floor processing
    sorted_qvalues <- sort(qvalues)
    n_5pct <- floor(n_genes * 0.05)
    
    if (n_5pct > 0) {
      q_threshold <- sorted_qvalues[n_5pct]
      potential_deg_indices <- which(qvalues < q_threshold)
      non_deg_indices <- which(qvalues >= q_threshold)
      
      if (length(potential_deg_indices) > 0) {
        exclusion_method <- "floor"
        cat(sprintf("      Floor processing: selected %d genes (%.2f%%)\n", 
                    length(potential_deg_indices), 
                    length(potential_deg_indices) / n_genes * 100))
      }
    } else {
      potential_deg_indices <- integer(0)
      non_deg_indices <- seq_len(n_genes)
    }
  }
  
  # Check termination
  terminate <- (length(potential_deg_indices) == 0 || n_genes == 0)
  
  norm_factors_with_names <- dgelist$samples$norm.factors
  names(norm_factors_with_names) <- rownames(dgelist$samples)
  
  return(list(
    norm_factors = norm_factors_with_names,
    potential_deg_indices = potential_deg_indices,
    non_deg_indices = non_deg_indices,
    terminate = terminate,
    deg_count = deg_count,
    deg_proportion = deg_proportion,
    n_excluded = length(potential_deg_indices),
    exclusion_method = exclusion_method,
    pi0_estimate = pi0_estimate,
    n_genes_input = n_genes
  ))
}

# ============================================================================
# Extract samples for each group
# ============================================================================

cat("\n--- Extracting paired samples for each group ---\n")

# Create sample lists using case-based extraction
sample_lists <- list()

for (group_name in c("R0", "R1", "B0", "B1")) {
  paired <- extract_paired_samples(group_name, high_purity_cases, sample_metadata)
  
  if (length(paired$cases) > 0) {
    sample_lists[[group_name]] <- paired
    cat(sprintf("  %s: %d pairs\n", group_name, length(paired$cases)))
  } else {
    cat(sprintf("  %s: No pairs available\n", group_name))
  }
}

saveRDS(sample_lists, paste0(paths$processed, "analysis_sample_lists.rds"))

# ============================================================================
# Process each comparison
# ============================================================================

cat("\n--- Processing comparisons ---\n")

thyr_deges_results <- list()
thyr_norm_factors <- list()

# Define comparisons
comparisons <- list(
  R0_vs_R1 = c("R0", "R1"),
  B0_vs_B1 = c("B0", "B1")
)

for (comp_name in names(comparisons)) {
  groups <- comparisons[[comp_name]]
  
  cat(sprintf("\n=== Comparison: %s ===\n", comp_name))
  
  # Check if both groups have samples
  if (!all(groups %in% names(sample_lists))) {
    cat("Skipping: Missing groups\n")
    next
  }
  
  if (length(sample_lists[[groups[1]]]$cases) == 0 || 
      length(sample_lists[[groups[2]]]$cases) == 0) {
    cat("Skipping: Empty groups\n")
    next
  }
  
  # Process tumor and normal separately
  for (tissue_type in c("tumor", "normal")) {
    cat(sprintf("\n--- %s %s ---\n", comp_name, tissue_type))
    
    # Get samples
    if (tissue_type == "tumor") {
      samples1 <- sample_lists[[groups[1]]]$tumor
      samples2 <- sample_lists[[groups[2]]]$tumor
    } else {
      samples1 <- sample_lists[[groups[1]]]$normal
      samples2 <- sample_lists[[groups[2]]]$normal
    }
    
    # Combine samples
    all_samples <- c(samples1, samples2)
    sample_groups <- c(rep(groups[1], length(samples1)),
                       rep(groups[2], length(samples2)))
    
    cat(sprintf("  Samples: %s=%d, %s=%d\n", 
                groups[1], length(samples1),
                groups[2], length(samples2)))
    
    # Extract counts
    count_matrix <- assay(se)[, all_samples, drop = FALSE]
    
    # Filter to protein_coding genes only
    gene_info <- as.data.frame(rowData(se))
    is_protein_coding <- gene_info$gene_type == "protein_coding"
    count_matrix <- count_matrix[is_protein_coding, , drop = FALSE]
    gene_info_pc <- gene_info[is_protein_coding, , drop = FALSE]
    
    cat(sprintf("  Protein coding genes: %d\n", nrow(count_matrix)))
    
    # ========================================================================
    # WORKFLOW: filterByExpr first
    # ========================================================================
    
    # Step 1: Apply filterByExpr
    cat("  Applying filterByExpr...\n")
    dgelist_initial <- DGEList(counts = count_matrix, group = factor(sample_groups))
    keep <- filterByExpr(dgelist_initial, group = dgelist_initial$samples$group)
    count_matrix_filtered <- count_matrix[keep, , drop = FALSE]
    gene_info_filtered <- gene_info_pc[keep, , drop = FALSE]
    cat(sprintf("    After filterByExpr: %d genes\n", nrow(count_matrix_filtered)))
    
    # Store the filtered gene set (this will be our consistent reference)
    filtered_gene_set <- which(keep)
    
    # Step 2: Apply Cook's distance outlier removal (for normalization only)
    cook_result <- detect_cook_outliers(count_matrix_filtered, sample_groups, 
                                        CONFIG$COOKS_QUANTILE)
    
    if (cook_result$outlier_count > 0) {
      count_matrix_for_norm <- count_matrix_filtered[!cook_result$outliers, , drop = FALSE]
      cat(sprintf("    For normalization: %d genes (after Cook's removal)\n", 
                  nrow(count_matrix_for_norm)))
    } else {
      count_matrix_for_norm <- count_matrix_filtered
      cat("    No Cook's outliers detected\n")
    }
    
    # Step 3: DEGES iterations (on Cook's-filtered subset)
    cat("  Starting DEGES iterations...\n")
    current_counts <- count_matrix_for_norm
    iteration_results <- list()
    
    for (iter in 0:(CONFIG$MAX_ITERATIONS - 1)) {
      result <- perform_deges_iteration(
        current_counts, 
        sample_groups, 
        iteration = iter,
        config = CONFIG
      )
      
      iteration_results[[paste0("iteration_", iter)]] <- result
      
      # Check termination
      if (result$terminate) {
        cat(sprintf("    DEGES terminated at iteration %d\n", iter))
        break
      }
      
      if (iter == CONFIG$MAX_ITERATIONS - 1) {
        cat(sprintf("    DEGES reached maximum iterations (%d)\n", CONFIG$MAX_ITERATIONS))
        break
      }
      
      # Prepare for next iteration: use only non-DEG genes
      if (length(result$non_deg_indices) > 0) {
        current_counts <- count_matrix_for_norm[result$non_deg_indices, , drop = FALSE]
        cat(sprintf("    Genes for next iteration: %d\n", nrow(current_counts)))
      } else {
        # No non-DEG genes left (all genes are DEGs)
        cat("    No non-DEG genes for next iteration, stopping\n")
        break
      }
    }
    
    # ========================================================================
    # Step 4: Apply norm.factors to filterByExpr-filtered data (before Cook's)
    # ========================================================================
    
    final_result <- iteration_results[[length(iteration_results)]]
    
    # Create DGEList with filterByExpr-filtered genes (not Cook's filtered)
    # NOTE: norm.factors must be set AFTER DGEList creation (edgeR ignores it in samples arg)
    dgelist_final <- DGEList(
      counts = count_matrix_filtered,  # After filterByExpr, before Cook's
      group = factor(sample_groups),
      genes = gene_info_filtered
    )
    
    # Apply DEGES-MUREN normalization factors
    dgelist_final$samples$norm.factors <- final_result$norm_factors[all_samples]
    
    # Verify norm.factors were applied correctly
    nf_range <- range(dgelist_final$samples$norm.factors)
    cat(sprintf("  Norm factors applied: range [%.4f, %.4f]\n", nf_range[1], nf_range[2]))
    
    # Calculate normalized CPM values
    cat("  Calculating normalized CPM values...\n")
    normalized_cpm <- cpm(dgelist_final, normalized.lib.sizes = TRUE, 
                          prior.count = 0, log = FALSE)
    
    # Store results
    comp_tissue <- paste(comp_name, tissue_type, sep = "_")
    
    # Save DGEList
    dgelist_filename <- paste0(paths$processed, "analysis_dgelist_", comp_tissue, ".rds")
    saveRDS(dgelist_final, dgelist_filename)
    cat(sprintf("  DGEList saved: %s\n", basename(dgelist_filename)))
    
    # Save normalized CPM
    cpm_filename <- paste0(paths$processed, "analysis_cpm_", comp_tissue, ".rds")
    saveRDS(normalized_cpm, cpm_filename)
    cat(sprintf("  Normalized CPM saved: %s\n", basename(cpm_filename)))
    
    # Store summary information
    thyr_deges_results[[comp_tissue]] <- list(
      comparison = comp_name,
      tissue = tissue_type,
      groups = groups,
      n_samples = list(
        group1 = length(samples1),
        group2 = length(samples2)
      ),
      n_genes_initial = sum(is_protein_coding),
      n_genes_filtered = length(filtered_gene_set),
      cook_outliers = cook_result$outlier_count,
      iterations = length(iteration_results),
      iteration_degs = sapply(iteration_results, function(x) x$deg_count),
      iteration_excluded = unname(sapply(iteration_results, function(x) x$n_excluded)),
      iteration_methods = unname(sapply(iteration_results, function(x) x$exclusion_method)),
      final_deg_count = final_result$deg_count,
      normalized_cpm = normalized_cpm
    )
    
    # Save normalization factors
    thyr_norm_factors[[comp_tissue]] <- final_result$norm_factors
    
    # Log processing details (expanded)
    deges_processing_log$comparisons[[comp_tissue]] <- list(
      # Basic info
      comparison = comp_name,
      tissue = tissue_type,
      date = Sys.Date(),
      
      # Sample info
      n_samples = list(
        group1 = length(samples1), 
        group2 = length(samples2)
      ),
      sample_ids = list(
        group1 = samples1,
        group2 = samples2
      ),
      
      # Gene filtering
      n_protein_coding = sum(is_protein_coding),
      n_after_filterByExpr = length(filtered_gene_set),
      
      # Cook's distance details
      cook_outliers = list(
        n_outliers = cook_result$outlier_count,
        threshold = cook_result$threshold,
        outlier_gene_ids = cook_result$outlier_gene_ids
      ),
      
      # DEGES iterations (detailed)
      iterations_summary = lapply(names(iteration_results), function(iter_name) {
        x <- iteration_results[[iter_name]]
        list(
          iteration = iter_name,
          n_genes_input = x$n_genes_input,
          deg_count = x$deg_count,
          n_excluded = x$n_excluded,
          exclusion_method = x$exclusion_method,
          deg_proportion = x$deg_proportion,
          non_deg_count = length(x$non_deg_indices),
          pi0_estimate = x$pi0_estimate,
          terminated = x$terminate
        )
      }),
      
      # Final results
      final_iteration = length(iteration_results),
      final_deg_count = final_result$deg_count,
      final_pi0 = final_result$pi0_estimate,
      n_genes_output = nrow(count_matrix_filtered)
    )
  }
}

# ============================================================================
# Create summary
# ============================================================================

cat("\n--- Creating summary ---\n")

summary_data <- data.frame()

for (comp_name in names(thyr_deges_results)) {
  result <- thyr_deges_results[[comp_name]]
  iter_excluded <- result$iteration_excluded
  iter_methods <- result$iteration_methods
  
  summary_row <- data.frame(
    comparison = result$comparison,
    tissue = result$tissue,
    n_group1 = result$n_samples$group1,
    n_group2 = result$n_samples$group2,
    genes_initial = result$n_genes_initial,
    genes_filtered = result$n_genes_filtered,
    cook_outliers = result$cook_outliers,
    iterations = result$iterations,
    iter0_excluded = iter_excluded[1],
    iter0_method = iter_methods[1],
    iter1_excluded = ifelse(length(iter_excluded) > 1, iter_excluded[2], NA),
    iter1_method = ifelse(length(iter_methods) > 1, iter_methods[2], NA),
    final_degs = result$final_deg_count,
    stringsAsFactors = FALSE
  )
  
  summary_data <- rbind(summary_data, summary_row)
}

print(summary_data)

# ============================================================================
# Save results
# ============================================================================

cat("\n--- Saving results ---\n")

# Save main DEGES results
deges_output <- list(
  date = Sys.Date(),
  config = CONFIG,
  case_master = high_purity_cases,
  sample_lists = sample_lists,
  results = thyr_deges_results,
  summary = summary_data,
  version = "v7.8"
)

saveRDS(deges_output, paste0(paths$processed, "analysis_deges_results.rds"))
cat("  DEGES results saved: analysis_deges_results.rds\n")

# Save normalization factors
saveRDS(thyr_norm_factors, paste0(paths$processed, "analysis_norm_factors.rds"))
cat("  Normalization factors saved: analysis_norm_factors.rds\n")

# Export summary as CSV
write.csv(summary_data, 
          paste0(paths$output, "analysis_deges_summary.csv"),
          row.names = FALSE)
cat("  Summary CSV saved: analysis_deges_summary.csv\n")

# Save processing log
log_file <- paste0(paths$logs, "deges_processing_info_", 
                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
saveRDS(deges_processing_log, log_file)
cat("  Processing log saved to logs/\n")

# ============================================================================
# DEGES Processing Complete - Clean up before visualization
# ============================================================================

cat("\n=== DEGES Processing Complete ===\n")
cat("All results saved. Cleaning up intermediate variables...\n")

# Keep only essential objects for next steps
rm(list = setdiff(ls(), c("paths", "thyr_case_master_stage2_filtered", 
                          "thyr_se_strand2_nonzero", "CONFIG")))
gc()

cat("Memory cleaned. MA Plot generation can be run independently from here.\n")

# ============================================================================
# MA Plots for normalization quality assessment
# (This section can be re-run independently after DEGES processing)
# ============================================================================

cat("\n--- Generating MA Plots ---\n")

# Ensure paths is available (for independent execution)
if (!exists("paths")) {
  source("analysis_v7/setup.R")
}

# Load saved results
cat("Loading saved DEGES results...\n")
deges_output <- readRDS(paste0(paths$processed, "analysis_deges_results.rds"))
comparison_names <- names(deges_output$results)
cat(sprintf("  Found %d comparisons: %s\n", 
            length(comparison_names), 
            paste(comparison_names, collapse = ", ")))

library(ggplot2)

ma_plot_data <- list()

for (comp_tissue in comparison_names) {
  cat(sprintf("  Processing %s...\n", comp_tissue))
  
  # Load DGEList
  dgelist_file <- paste0(paths$processed, "analysis_dgelist_", comp_tissue, ".rds")
  if (!file.exists(dgelist_file)) {
    cat(sprintf("    Skipping: DGEList not found\n"))
    next
  }
  
  dgelist <- readRDS(dgelist_file)
  
  # Calculate normalized log2 CPM
  cpm_log2 <- cpm(dgelist, normalized.lib.sizes = TRUE, prior.count = 1, log = TRUE)
  
  # Get group information
  groups <- levels(dgelist$samples$group)
  group1_samples <- rownames(dgelist$samples)[dgelist$samples$group == groups[1]]
  group2_samples <- rownames(dgelist$samples)[dgelist$samples$group == groups[2]]
  
  # Calculate M and A values
  mean_g1 <- rowMeans(cpm_log2[, group1_samples, drop = FALSE])
  mean_g2 <- rowMeans(cpm_log2[, group2_samples, drop = FALSE])
  
  M <- mean_g2 - mean_g1  # log2 FC (group2 vs group1)
  A <- (mean_g2 + mean_g1) / 2
  
  # Store statistics
  ma_plot_data[[comp_tissue]] <- list(
    M_median = median(M),
    M_mean = mean(M),
    pct_positive = sum(M > 0) / length(M) * 100,
    n_genes = length(M),
    nf_g1_median = median(dgelist$samples$norm.factors[dgelist$samples$group == groups[1]]),
    nf_g2_median = median(dgelist$samples$norm.factors[dgelist$samples$group == groups[2]])
  )
  
  # Create MA plot
  ma_df <- data.frame(A = A, M = M)
  
  p <- ggplot(ma_df, aes(x = A, y = M)) +
    geom_point(alpha = 0.2, size = 0.5, color = "gray30") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = median(M), color = "blue", linewidth = 0.8) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, 
                color = "darkgreen", linewidth = 1) +
    annotate("text", x = max(A) - 2, y = 2.5, 
             label = sprintf("M median = %.3f", median(M)), 
             color = "blue", size = 4) +
    labs(
      title = sprintf("MA Plot: %s (DEGES-MUREN normalized)", comp_tissue),
      subtitle = sprintf("n=%d genes | Red=0, Blue=median, Green=loess | %s vs %s",
                         nrow(ma_df), groups[2], groups[1]),
      x = "A (Average log2 expression)",
      y = sprintf("M (log2 FC: %s vs %s)", groups[2], groups[1])
    ) +
    theme_bw() +
    coord_cartesian(ylim = c(-3, 3))
  
  # Display in RStudio
  print(p)
  
  # Save PDF
  pdf_file <- paste0(paths$output, "maplot_", comp_tissue, ".pdf")
  ggsave(pdf_file, p, width = 8, height = 6)
  cat(sprintf("    Saved: %s\n", basename(pdf_file)))
}

# Print MA plot summary
cat("\n=== MA Plot Summary ===\n")
cat(sprintf("%-20s %10s %10s %12s %10s %10s\n", 
            "Comparison", "M_median", "M_mean", "Pct_M>0", "NF_G1_med", "NF_G2_med"))
cat(paste(rep("-", 75), collapse = ""), "\n")

for (comp_tissue in names(ma_plot_data)) {
  stats <- ma_plot_data[[comp_tissue]]
  cat(sprintf("%-20s %10.4f %10.4f %11.1f%% %10.3f %10.3f\n",
              comp_tissue, stats$M_median, stats$M_mean, stats$pct_positive,
              stats$nf_g1_median, stats$nf_g2_median))
}

# ============================================================================
# Final report
# ============================================================================

cat("\n=== DEGES Normalization Complete (v7.9) ===\n")
cat("Configuration:\n")
cat("  Workflow: filterByExpr -> Cook's -> DEGES iterations\n")
cat("  DEG screening: Brunner-Munzel test + Storey q-value (lambda=0.5)\n")
cat("  Gene set: Consistent (filterByExpr-filtered) throughout\n")
cat("  Output: Normalized CPM values (prior.count = 0)\n")
cat("\nProcessed comparisons:\n")

for (comp_name in comparison_names) {
  result <- deges_output$results[[comp_name]]
  cat(sprintf("  %s: %d genes, %d iterations\n",
              comp_name, 
              result$n_genes_filtered,
              result$iterations))
}

cat("\nOutputs:\n")
cat("  Main: analysis_deges_results.rds\n")
cat("  Factors: analysis_norm_factors.rds\n")
cat("  DGELists: analysis_dgelist_*.rds (4 files)\n")
cat("  CPM values: analysis_cpm_*.rds (4 files)\n")
cat("  MA plots: maplot_*.pdf (4 files)\n")
cat("  Summary: analysis_deges_summary.csv\n")
cat("\nNext: Run 09_deg_analysis.R for differential expression\n")

# Clean up (keep SE and case_master for potential subsequent scripts)
rm(list = setdiff(ls(), c("paths", "thyr_case_master_stage2_filtered", 
                          "thyr_se_strand2_nonzero")))
gc()