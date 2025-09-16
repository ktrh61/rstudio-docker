# 06_purity_analysis.R - Stage 2 Tumor Purity Estimation
# Purpose: Estimate tumor purity using ContamDE for clean paired samples
# Method: MUREN normalization with contamDE purity estimation  
# Input: thyr_case_master_stage1_filtered, thyr_se_strand2_nonzero
# Output: thyr_case_master_stage2_filtered with tumor_purity and low_purity flags (clean cases only)
# Version: v7.1 - Serial outer loop, parallel MUREN (3 threads)
# Date: 2025-01-20

source("analysis_v7/setup.R")

cat("\n=== Stage 2: Tumor Purity Analysis (v7.1) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Analysis: R0/R1/B0/B1 groups with clean paired samples\n")
cat("Method: ContamDE with MUREN normalization\n")
cat("Processing: Serial group processing, parallel MUREN (3 threads)\n")

# Load packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
  library(qvalue)
  library(dplyr)
})

# Load utility functions
source("utils/utils_improved.R")
source("utils/norm_improved.R")
source("utils/contamde_purity_functions.R")

# ============================================================================
# Configuration
# ============================================================================

CONFIG <- list(
  PURITY_THRESHOLD = 0.6,      # 60% minimum purity
  PAIRWISE_METHOD = "lts",     # MUREN method: lts, median, trim10
  WORKERS = 3,                 # Fixed 3 cores for MUREN internal parallelization
  PRIOR_COUNT = 0.5,           # Pseudocount for log ratio stability
  VERBOSE = TRUE               # Verbose output for contamDE
)

cat("\nConfiguration:\n")
cat("  Purity threshold:", sprintf("%.0f%%", CONFIG$PURITY_THRESHOLD * 100), "\n")
cat("  MUREN method:", CONFIG$PAIRWISE_METHOD, "\n")
cat("  MUREN workers:", CONFIG$WORKERS, "(internal parallelization)\n")
cat("  Prior count:", CONFIG$PRIOR_COUNT, "\n")
cat("  Processing: Serial (group-by-group)\n")

# Thread control for BLAS
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  cat("  BLAS/OMP threads: 1 (to avoid nested parallelization)\n")
}

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
cat("SE dimensions:", format(nrow(se), big.mark=","), "genes × ", 
    format(ncol(se), big.mark=","), "samples\n")

# Load case_master_stage1_filtered
if (exists("thyr_case_master_stage1_filtered")) {
  cat("Using existing thyr_case_master_stage1_filtered\n")
  case_master <- thyr_case_master_stage1_filtered
} else {
  case_master_path <- paste0(paths$processed, "thyr_case_master_stage1_filtered.rds")
  if (!file.exists(case_master_path)) {
    stop("thyr_case_master_stage1_filtered.rds not found. Please run 05_pca_outlier_detection.R first.")
  }
  cat("Loading case_master from file...\n")
  case_master <- readRDS(case_master_path)
}
cat("Case master loaded:", nrow(case_master), "cases\n")

# Get sample metadata
sample_metadata <- as.data.frame(colData(se))

# Get gene information
gene_info <- as.data.frame(rowData(se))

# ============================================================================
# Filter to clean cases only (no outliers)
# ============================================================================

cat("\n--- Filtering to clean cases (no outliers) ---\n")

# Keep only cases without outliers - these proceed to Stage 2
clean_cases <- case_master %>%
  filter(has_outlier_tumor == 0 & has_outlier_normal == 0)

cat("Clean cases:", nrow(clean_cases), "/", nrow(case_master), 
    sprintf("(%.1f%%)\n", nrow(clean_cases)/nrow(case_master)*100))

# Group distribution of clean cases
clean_summary <- clean_cases %>%
  group_by(group) %>%
  summarise(n = n(), .groups = "drop")

cat("\nClean cases by group:\n")
print(clean_summary)

# ============================================================================
# Helper functions
# ============================================================================

# Extract paired samples for a group (from clean cases only)
extract_paired_samples <- function(group_name, clean_cases_df, sample_meta) {
  # Get cases for this group (already filtered for clean)
  group_cases <- clean_cases_df$case_submitter_id[clean_cases_df$group == group_name]
  
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
    if (length(tumor_sample) > 0) tumor_ids[i] <- tumor_sample[1]
    if (length(normal_sample) > 0) normal_ids[i] <- normal_sample[1]
  }
  
  # Remove any NA (shouldn't happen with clean paired cases)
  valid <- !is.na(tumor_ids) & !is.na(normal_ids)
  
  return(list(
    tumor = tumor_ids[valid],
    normal = normal_ids[valid],
    cases = group_cases[valid]
  ))
}

# Gene filtering for purity estimation
filter_genes_for_purity <- function(se, tumor_ids, normal_ids, gene_info) {
  # 1. Protein coding genes only
  is_protein_coding <- gene_info$gene_type == "protein_coding"
  
  # 2. Extract counts for the paired samples
  all_sample_ids <- c(normal_ids, tumor_ids)  # Normal first, then tumor (contamDE convention)
  counts_group <- assay(se)[, all_sample_ids, drop = FALSE]
  
  # 3. Create group factor for filterByExpr (normal vs tumor)
  n_pairs <- length(tumor_ids)
  group <- factor(c(rep("Normal", n_pairs), rep("Tumor", n_pairs)))
  
  # 4. Apply filterByExpr
  keep_expr <- edgeR::filterByExpr(counts_group, group = group)
  
  # 5. Combine conditions
  keep_genes <- is_protein_coding & keep_expr
  
  return(keep_genes)
}

# Prepare counts matrix for ContamDE (normal first, then tumor)
prepare_contamde_matrix <- function(se, normal_ids, tumor_ids, keep_genes) {
  # Extract counts with gene filtering
  counts <- cbind(
    assay(se)[keep_genes, normal_ids, drop = FALSE],
    assay(se)[keep_genes, tumor_ids, drop = FALSE]
  )
  
  # Set informative column names
  n_pairs <- length(normal_ids)
  colnames(counts) <- c(
    paste0("Normal_", seq_len(n_pairs)),
    paste0("Tumor_", seq_len(n_pairs))
  )
  
  return(counts)
}

# ============================================================================
# Process each group (SERIAL)
# ============================================================================

cat("\n--- Processing groups (serial) ---\n")
cat("Note: Each group processed sequentially, MUREN uses internal parallelization\n")

# Initialize results storage
purity_results <- list()
processing_times <- list()

# Process each group sequentially
for (group_name in c("R0", "R1", "B0", "B1")) {
  cat(sprintf("\n[%s] Processing %s...\n", Sys.time(), group_name))
  start_time <- Sys.time()
  
  # Extract paired samples (clean cases only)
  paired <- extract_paired_samples(group_name, clean_cases, sample_metadata)
  
  if (length(paired$cases) == 0) {
    cat("  No clean paired samples available\n")
    processing_times[[group_name]] <- 0
    next
  }
  
  cat(sprintf("  Clean pairs: %d\n", length(paired$cases)))
  
  # Gene filtering
  keep_genes <- filter_genes_for_purity(se, paired$tumor, paired$normal, gene_info)
  cat(sprintf("  Genes: %d total, %d protein coding, %d after filterByExpr\n",
              nrow(se), sum(gene_info$gene_type == "protein_coding"), sum(keep_genes)))
  
  # Prepare counts matrix
  counts <- prepare_contamde_matrix(se, paired$normal, paired$tumor, keep_genes)
  cat(sprintf("  Count matrix: %d genes × %d samples\n", nrow(counts), ncol(counts)))
  
  # Run ContamDE purity estimation (with internal MUREN parallelization)
  cat("  Running ContamDE with MUREN (3 parallel workers)...\n")
  
  tryCatch({
    purity_result <- contamde_purity(
      counts = counts,
      subtype = NULL,
      covariate = NULL,
      contaminated = TRUE,
      pairwise_method = CONFIG$PAIRWISE_METHOD,
      workers = CONFIG$WORKERS,  # 3 workers for internal MUREN parallelization
      prior.count = CONFIG$PRIOR_COUNT,
      verbose = FALSE  # Suppress internal messages
    )
    
    # Store results with case IDs
    purity_results[[group_name]] <- list(
      cases = paired$cases,
      purity = purity_result$proportion,
      n_pairs = length(paired$cases),
      n_genes = sum(keep_genes),
      informative_genes = length(purity_result$informative_genes$up) + 
        length(purity_result$informative_genes$down)
    )
    
    # Summary statistics
    purity_vals <- purity_result$proportion
    cat(sprintf("  Purity: mean=%.3f, median=%.3f, sd=%.3f\n",
                mean(purity_vals), median(purity_vals), sd(purity_vals)))
    cat(sprintf("  Range: [%.3f, %.3f]\n", min(purity_vals), max(purity_vals)))
    cat(sprintf("  Above %.0f%%: %d/%d (%.1f%%)\n",
                CONFIG$PURITY_THRESHOLD * 100,
                sum(purity_vals >= CONFIG$PURITY_THRESHOLD),
                length(purity_vals),
                sum(purity_vals >= CONFIG$PURITY_THRESHOLD) / length(purity_vals) * 100))
    
  }, error = function(e) {
    cat("  Error in purity estimation:", e$message, "\n")
    purity_results[[group_name]] <- NULL
  })
  
  # Record processing time
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  processing_times[[group_name]] <- elapsed
  cat(sprintf("  Processing time: %.1f seconds\n", elapsed))
}

# Report total processing time
total_time <- sum(unlist(processing_times))
cat(sprintf("\nTotal processing time: %.1f seconds\n", total_time))

# ============================================================================
# Create case_master_stage2_filtered (CLEAN CASES ONLY)
# ============================================================================

cat("\n--- Creating case_master_stage2_filtered (clean cases only) ---\n")

# Start with clean cases only - Stage 2 contains ONLY non-outlier cases
thyr_case_master_stage2_filtered <- clean_cases

# Initialize new columns for purity information
thyr_case_master_stage2_filtered$tumor_purity <- NA_real_
thyr_case_master_stage2_filtered$low_purity <- 0  # Default 0 for clean cases

# Reorder columns to place purity info after has_outlier columns
col_order <- names(thyr_case_master_stage2_filtered)
has_outlier_idx <- which(col_order == "has_outlier_normal")
if (length(has_outlier_idx) > 0) {
  new_order <- c(
    col_order[1:has_outlier_idx],
    "tumor_purity", "low_purity",
    col_order[(has_outlier_idx + 1):(length(col_order) - 2)]  # Exclude the two purity columns at the end
  )
  thyr_case_master_stage2_filtered <- thyr_case_master_stage2_filtered[, new_order, with = FALSE]
}

# Fill in purity values for clean cases
for (group_name in names(purity_results)) {
  if (!is.null(purity_results[[group_name]])) {
    result <- purity_results[[group_name]]
    
    for (i in seq_along(result$cases)) {
      case_id <- result$cases[i]
      purity_value <- result$purity[i]
      
      # Find the case in case_master
      idx <- which(thyr_case_master_stage2_filtered$case_submitter_id == case_id)
      
      if (length(idx) == 1) {
        thyr_case_master_stage2_filtered$tumor_purity[idx] <- purity_value
        
        # Set low_purity flag
        if (!is.na(purity_value) && purity_value < CONFIG$PURITY_THRESHOLD) {
          thyr_case_master_stage2_filtered$low_purity[idx] <- 1
        }
      }
    }
  }
}

cat(sprintf("Stage 2 filtered dataset: %d clean cases (outlier cases excluded)\n",
            nrow(thyr_case_master_stage2_filtered)))

# ============================================================================
# Summary statistics (clean cases only)
# ============================================================================

cat("\n--- Summary by group (clean cases only) ---\n")

summary_stats <- thyr_case_master_stage2_filtered %>%
  group_by(group) %>%
  summarise(
    n_clean_cases = n(),
    n_purity_measured = sum(!is.na(tumor_purity)),
    mean_purity = mean(tumor_purity, na.rm = TRUE),
    median_purity = median(tumor_purity, na.rm = TRUE),
    sd_purity = sd(tumor_purity, na.rm = TRUE),
    n_low_purity = sum(low_purity == 1, na.rm = TRUE),
    n_high_purity = sum(low_purity == 0 & !is.na(tumor_purity), na.rm = TRUE),
    .groups = "drop"
  )

print(as.data.frame(summary_stats))

# Overall statistics
total_clean <- nrow(thyr_case_master_stage2_filtered)
high_purity <- sum(thyr_case_master_stage2_filtered$low_purity == 0 & 
                     !is.na(thyr_case_master_stage2_filtered$tumor_purity))
measured <- sum(!is.na(thyr_case_master_stage2_filtered$tumor_purity))

cat(sprintf("\nOverall (clean cases only):\n"))
cat(sprintf("  Total clean cases: %d\n", total_clean))
cat(sprintf("  Purity measured: %d (%.1f%%)\n", measured, measured/total_clean*100))
cat(sprintf("  High purity: %d (%.1f%% of measured)\n", 
            high_purity, high_purity/measured*100))

# ============================================================================
# Save results
# ============================================================================

cat("\n--- Saving results ---\n")

# Save updated case_master (CLEAN CASES ONLY)
output_file <- paste0(paths$processed, "thyr_case_master_stage2_filtered.rds")
saveRDS(thyr_case_master_stage2_filtered, output_file)
cat("  Case master saved:", basename(output_file), "\n")
cat(sprintf("    Contains: %d clean cases with purity information\n", 
            nrow(thyr_case_master_stage2_filtered)))
cat("    Note: Outlier cases from Stage 1 are excluded\n")

# Save purity details for reference
purity_details <- list(
  date = Sys.Date(),
  config = CONFIG,
  results = purity_results,
  summary = summary_stats,
  processing_times = processing_times,
  total_time = total_time,
  processing_mode = "Serial groups, parallel MUREN"
)
saveRDS(purity_details, paste0(paths$output, "stage2_purity_details.rds"))
cat("  Purity details saved: stage2_purity_details.rds\n")

# Export summary as CSV
write.csv(summary_stats, 
          paste0(paths$output, "stage2_purity_summary.csv"),
          row.names = FALSE)
cat("  Summary CSV saved: stage2_purity_summary.csv\n")

# Export case-level purity data for review (clean cases only)
case_purity_export <- thyr_case_master_stage2_filtered %>%
  select(case_submitter_id, group, tumor_purity, low_purity) %>%
  arrange(group, case_submitter_id)

write.csv(case_purity_export,
          paste0(paths$output, "stage2_case_purity.csv"),
          row.names = FALSE)
cat("  Case purity CSV saved: stage2_case_purity.csv\n")

# ============================================================================
# Final report
# ============================================================================

cat("\n=== Stage 2 Complete ===\n")
cat("Configuration:\n")
cat("  Processing: Serial (group-by-group)\n")
cat("  MUREN: 3 parallel workers (internal)\n")
cat("  Purity threshold:", sprintf("%.0f%%", CONFIG$PURITY_THRESHOLD * 100), "\n")
cat("  Gene filtering: protein_coding + filterByExpr\n")
cat("\nResults:\n")
cat("  Input cases (Stage 1):", nrow(case_master), "\n")
cat("  Clean cases (no outliers):", total_clean, sprintf("(%.1f%%)\n", total_clean/nrow(case_master)*100))
cat("  Purity measured:", measured, "\n")
cat("  High purity (≥60%):", high_purity, sprintf("(%.1f%% of measured)\n", high_purity/measured*100))
cat("\nProcessing performance:\n")
cat("  Total time:", sprintf("%.1f seconds\n", total_time))
cat("  Average per group:", sprintf("%.1f seconds\n", total_time/4))
cat("\nOutputs:\n")
cat("  Main: thyr_case_master_stage2_filtered.rds (clean cases only)\n")
cat("  Details: stage2_purity_details.rds\n")
cat("  Summaries: stage2_purity_summary.csv, stage2_case_purity.csv\n")
cat("\nNext steps:\n")
cat("  - Use high purity cases (low_purity == 0) for downstream analysis\n")
cat("  - Stage 2 dataset contains ONLY clean cases from Stage 1\n")

# Restore SE to original name
thyr_se_strand2_nonzero <- se

# Clean up
rm(list = setdiff(ls(), c("paths", "thyr_case_master_stage2_filtered", "thyr_se_strand2_nonzero")))
gc()