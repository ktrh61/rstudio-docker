# 06_purity_analysis.R - Stage 2 Tumor Purity Estimation
# Purpose: Estimate relative tumor purity scores using ContamDE for clean paired samples
# Method: MUREN normalization with contamDE purity estimation
# Note: ContamDE outputs normalized proportions reflecting relative contamination intensities,
#       not absolute tumor purity percentages. The highest-purity sample is scaled to ~1.0.
# Input: thyr_case_master_stage1_filtered, thyr_se_strand2_nonzero
# Output: thyr_case_master_stage2_filtered with tumor_purity (relative score) and low_purity flags
# Version: v7.4 - Fixed dplyr::select namespace conflict
#                 New group design (R0/R_Low/R_Mid/R1/B0/B1)
#                 B_Low/B_Mid excluded (not used in validation)
#                 Output includes ALL cases (not just clean) for layer-based validation
# Date: 2025-12-07

source("analysis_v7/setup.R")

cat("\n=== Stage 2: Tumor Purity Analysis (v7.4) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Analysis: R0/R_Low/R_Mid/R1/B0/B1 groups with clean paired samples\n")
cat("Note: B_Low/B_Mid excluded (not used in validation)\n")
cat("Method: ContamDE with MUREN normalization\n")
cat("Output: Relative purity scores (normalized, not absolute %)\n")

# Load packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
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
  PURITY_THRESHOLD = 0.6,      # Relative score threshold (ContamDE outputs normalized scores, not absolute %)
  PAIRWISE_METHOD = "lts",     # MUREN method: lts, median, trim10
  WORKERS = 3,                 # Fixed 3 cores for MUREN internal parallelization
  PRIOR_COUNT = 1.0,           # Pseudocount for log ratio stability
  VERBOSE = TRUE               # Verbose output for contamDE
)

cat("\nConfiguration:\n")
cat("  Purity score threshold:", CONFIG$PURITY_THRESHOLD, "(relative, not %)\n")
cat("  MUREN method:", CONFIG$PAIRWISE_METHOD, "\n")
cat("  Prior count:", CONFIG$PRIOR_COUNT, "\n")

# Thread control
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  cat("  BLAS/OMP threads: 1\n")
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
cat("SE dimensions:", format(nrow(se), big.mark=","), "genes x",
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
# Filter to clean cases only (FIXED: removed low_purity reference)
# ============================================================================

cat("\n--- Filtering to clean cases ---\n")

# Keep only cases without outliers
# NOTE: low_purity will be determined later in this script
clean_cases <- case_master %>%
  filter(has_outlier_tumor == 0 & has_outlier_normal == 0)

cat("Clean cases:", nrow(clean_cases), "/", nrow(case_master), 
    sprintf("(%.1f%%)\n", nrow(clean_cases)/nrow(case_master)*100))

# Group distribution
clean_summary <- clean_cases %>%
  group_by(group) %>%
  summarise(n = n(), .groups = "drop")
print(clean_summary)

# ============================================================================
# Helper functions
# ============================================================================

# Extract paired samples for a group
extract_paired_samples <- function(group_name, clean_cases_df, sample_meta) {
  # Get paired cases for this group (explicit has_pair check for safety)
  group_cases <- clean_cases_df$case_submitter_id[
    clean_cases_df$group == group_name & clean_cases_df$has_pair
  ]
  
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
  
  # Remove any NA (shouldn't happen with clean cases)
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
# Process each group
# ============================================================================

cat("\n--- Processing groups ---\n")

# Initialize results storage
purity_results <- list()

for (group_name in c("R0", "R_Low", "R_Mid", "R1", "B0", "B1")) {
  cat(sprintf("\nProcessing %s...\n", group_name))
  
  # Extract paired samples
  paired <- extract_paired_samples(group_name, clean_cases, sample_metadata)
  
  if (length(paired$cases) == 0) {
    cat("  No clean paired samples available\n")
    next
  }
  
  cat(sprintf("  Pairs: %d\n", length(paired$cases)))
  
  # Gene filtering
  keep_genes <- filter_genes_for_purity(se, paired$tumor, paired$normal, gene_info)
  cat(sprintf("  Genes: %d total, %d protein coding, %d after filterByExpr\n",
              nrow(se), sum(gene_info$gene_type == "protein_coding"), sum(keep_genes)))
  
  # Prepare counts matrix
  counts <- prepare_contamde_matrix(se, paired$normal, paired$tumor, keep_genes)
  cat(sprintf("  Count matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))
  
  # Run ContamDE purity estimation
  cat("  Running ContamDE...\n")
  
  tryCatch({
    purity_result <- contamde_purity(
      counts = counts,
      subtype = NULL,
      covariate = NULL,
      contaminated = TRUE,
      pairwise_method = CONFIG$PAIRWISE_METHOD,
      workers = CONFIG$WORKERS,
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
    cat(sprintf("  Relative score: mean=%.3f, median=%.3f, sd=%.3f\n",
                mean(purity_vals), median(purity_vals), sd(purity_vals)))
    cat(sprintf("  Range: [%.3f, %.3f]\n", min(purity_vals), max(purity_vals)))
    cat(sprintf("  Above threshold (%.2f): %d/%d (%.1f%%)\n",
                CONFIG$PURITY_THRESHOLD,
                sum(purity_vals >= CONFIG$PURITY_THRESHOLD),
                length(purity_vals),
                sum(purity_vals >= CONFIG$PURITY_THRESHOLD) / length(purity_vals) * 100))
    
  }, error = function(e) {
    cat("  Error in purity estimation:", e$message, "\n")
    purity_results[[group_name]] <- NULL
  })
}

# ============================================================================
# Update case_master with purity information
# ============================================================================

cat("\n--- Updating case_master with purity results ---\n")

# Start with ALL cases (not just clean cases)
# This allows downstream filtering by different criteria
thyr_case_master_stage2_filtered <- case_master

# Initialize new columns (NA by default)
thyr_case_master_stage2_filtered$tumor_purity <- NA_real_
thyr_case_master_stage2_filtered$low_purity <- NA_integer_

# Reorder columns to place purity info after has_outlier columns
col_order <- names(thyr_case_master_stage2_filtered)
has_outlier_idx <- which(col_order == "has_outlier_normal")
if (length(has_outlier_idx) > 0) {
  new_order <- c(
    col_order[1:has_outlier_idx],
    "tumor_purity", "low_purity",
    col_order[(has_outlier_idx + 1):(length(col_order) - 2)]  # Exclude the two purity columns at the end
  )
  # Use dplyr select for data.frame compatibility
  thyr_case_master_stage2_filtered <- thyr_case_master_stage2_filtered %>%
    dplyr::select(all_of(new_order))
}

# Fill in purity values
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
        
        # Set low_purity flag (0 or 1, never NA since these are clean cases)
        if (!is.na(purity_value)) {
          if (purity_value < CONFIG$PURITY_THRESHOLD) {
            thyr_case_master_stage2_filtered$low_purity[idx] <- 1
          } else {
            thyr_case_master_stage2_filtered$low_purity[idx] <- 0
          }
        }
      }
    }
  }
}

# ============================================================================
# Summary statistics
# ============================================================================

cat("\n--- Summary by group ---\n")

summary_stats <- thyr_case_master_stage2_filtered %>%
  group_by(group) %>%
  summarise(
    n_cases = n(),
    n_paired = sum(has_pair),
    n_qc_clean = sum(has_outlier_tumor == 0 & has_outlier_normal == 0, na.rm = TRUE),
    n_qc_outlier = sum(has_outlier_tumor == 1 | has_outlier_normal == 1, na.rm = TRUE),
    n_no_qc = sum(is.na(has_outlier_tumor)),
    n_purity_measured = sum(!is.na(tumor_purity)),
    mean_purity = mean(tumor_purity, na.rm = TRUE),
    median_purity = median(tumor_purity, na.rm = TRUE),
    n_low_purity = sum(low_purity == 1, na.rm = TRUE),
    n_high_purity = sum(low_purity == 0, na.rm = TRUE),
    .groups = "drop"
  )

print(as.data.frame(summary_stats))

# Overall statistics
total_cases <- nrow(thyr_case_master_stage2_filtered)
usable_cases <- sum(
  thyr_case_master_stage2_filtered$has_outlier_tumor == 0 & 
    thyr_case_master_stage2_filtered$has_outlier_normal == 0 & 
    thyr_case_master_stage2_filtered$low_purity == 0,
  na.rm = TRUE
)

cat(sprintf("\nOverall: %d/%d cases usable (%.1f%%)\n",
            usable_cases, total_cases, usable_cases/total_cases*100))

# ============================================================================
# Save results
# ============================================================================

cat("\n--- Saving results ---\n")

# Save updated case_master
output_file <- paste0(paths$processed, "thyr_case_master_stage2_filtered.rds")
saveRDS(thyr_case_master_stage2_filtered, output_file)
cat("  Case master saved:", basename(output_file), "\n")
cat(sprintf("    Contains: %d cases with purity information\n", 
            nrow(thyr_case_master_stage2_filtered)))

# Save purity details for reference
purity_details <- list(
  date = Sys.Date(),
  config = CONFIG,
  results = purity_results,
  summary = summary_stats
)
saveRDS(purity_details, paste0(paths$output, "stage2_purity_details.rds"))
cat("  Purity details saved: stage2_purity_details.rds\n")

# Export summary as CSV
write.csv(summary_stats, 
          paste0(paths$output, "stage2_purity_summary.csv"),
          row.names = FALSE)
cat("  Summary CSV saved: stage2_purity_summary.csv\n")

# Export case-level purity data for review
case_purity_export <- thyr_case_master_stage2_filtered %>%
  dplyr::select(case_submitter_id, group, has_pair, has_outlier_tumor, has_outlier_normal, 
                tumor_purity, low_purity) %>%
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
cat("  Purity score threshold:", CONFIG$PURITY_THRESHOLD, "(relative, not %)\n")
cat("  Gene filtering: protein_coding + filterByExpr\n")
cat("  MUREN method:", CONFIG$PAIRWISE_METHOD, "\n")
cat("\nResults:\n")
cat("  Total cases:", total_cases, "\n")

# Layer statistics for validation
n_qc_clean <- sum(thyr_case_master_stage2_filtered$has_outlier_tumor == 0 & 
                    thyr_case_master_stage2_filtered$has_outlier_normal == 0, na.rm = TRUE)
n_high_purity <- sum(thyr_case_master_stage2_filtered$low_purity == 0, na.rm = TRUE)
n_no_qc <- sum(is.na(thyr_case_master_stage2_filtered$has_outlier_tumor))

cat("  Layer 1 (QC clean + high score):", usable_cases, "\n")
cat("  Layer 2 (QC clean, any score):", n_qc_clean, "\n")
cat("  Layer 3 (No QC - unpaired/B_Low/B_Mid):", n_no_qc, "\n")
cat("\nOutputs:\n")
cat("  Main: thyr_case_master_stage2_filtered.rds\n")
cat("  Details: stage2_purity_details.rds\n")
cat("  Summaries: stage2_purity_summary.csv, stage2_case_purity.csv\n")
cat("\nFiltering guide:\n")
cat("  Layer 1: has_outlier_tumor==0 & has_outlier_normal==0 & low_purity==0\n")
cat("  Layer 2: has_outlier_tumor==0 & has_outlier_normal==0\n")
cat("  Layer 3: is.na(has_outlier_tumor)\n")

# Restore SE to original name
thyr_se_strand2_nonzero <- se

# Clean up
rm(list = setdiff(ls(), c("paths", "thyr_case_master_stage2_filtered", "thyr_se_strand2_nonzero")))
gc()