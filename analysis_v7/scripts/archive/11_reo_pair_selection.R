# 11_reo_pair_selection.R - REO-based Biomarker Pair Selection
# Purpose: Select stable REO pairs for radiation exposure detection
# Method: REO (Relative Expression Orderings) analysis on R0 vs R1 DEGs
# Input: thyr_deg_results.rds, normalized count data
# Output: Top REO pairs for clinical validation
# Version: v1.2 - Phase 1: Data Preparation, Phase 2: REO Pair Generation
# Date: 2025-01-23

source("analysis_v7/setup.R")

cat("\n=== REO-based Biomarker Pair Selection (v1.2) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Phase 1: Data Preparation\n")
cat("Phase 2: REO Pair Generation (if Phase 1 data exists)\n")

# Load required packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
  library(qvalue)
})

# ============================================================================
# Configuration
# ============================================================================

CONFIG <- list(
  # Input/Output
  VERBOSE = TRUE,
  
  # Filtering
  MIN_TPM = 0  # Exclude complete zeros only (no pseudocount)
)

cat("\nConfiguration:\n")
cat("  TPM calculation: Using GENCODE v36 gene lengths\n")
cat("  Filter: All samples must have TPM > 0\n")

# ============================================================================
# Phase 1: Data Preparation
# ============================================================================

cat("\n--- Phase 1: Data Preparation ---\n")

# --------------------------------------------------------------------------
# 1.1 Load DEG results
# --------------------------------------------------------------------------

cat("\n1.1 Loading DEG results...\n")

deg_results_path <- paste0(paths$processed, "thyr_deg_results.rds")
if (!file.exists(deg_results_path)) {
  stop("DEG results not found. Please run 09_deg_analysis.R first.")
}

deg_output <- readRDS(deg_results_path)

# Extract R0_vs_R1_tumor results
# Structure: deg_output$deg_results$R0_vs_R1_tumor$deg_summary$results_df
r0r1_tumor_result <- deg_output$deg_results$R0_vs_R1_tumor

if (is.null(r0r1_tumor_result)) {
  stop("R0_vs_R1_tumor results not found in DEG output.")
}

# Get the results dataframe
results_df <- r0r1_tumor_result$deg_summary$results_df

# Get significant DEGs (already marked in the significant column)
sig_degs <- results_df %>%
  dplyr::filter(significant == TRUE)

# Separate up and down regulated genes based on log2FC
up_genes <- sig_degs %>%
  dplyr::filter(log2FC > 0) %>%
  dplyr::pull(gene_id)

down_genes <- sig_degs %>%
  dplyr::filter(log2FC < 0) %>%
  dplyr::pull(gene_id)

cat(sprintf("  Total DEGs: %d (q < 0.05)\n", nrow(sig_degs)))
cat(sprintf("  Upregulated: %d\n", length(up_genes)))
cat(sprintf("  Downregulated: %d\n", length(down_genes)))

# --------------------------------------------------------------------------
# 1.2 Load count data and clinical info
# --------------------------------------------------------------------------

cat("\n1.2 Loading count data...\n")

# Load the SummarizedExperiment
se_path <- paste0(paths$processed, "thyr_se_strand2_nonzero.rds")
if (!file.exists(se_path)) {
  stop("SummarizedExperiment not found. Please run 03_prepare_counts.R first.")
}

se <- readRDS(se_path)

# Load case master for sample grouping
case_master_path <- paste0(paths$processed, "thyr_case_master_stage2_filtered.rds")
if (!file.exists(case_master_path)) {
  stop("Case master not found. Please run 06_purity_analysis.R first.")
}

case_master <- readRDS(case_master_path)

# Load sample lists (high purity tumor samples)
sample_lists_path <- paste0(paths$processed, "analysis_sample_lists.rds")
if (!file.exists(sample_lists_path)) {
  stop("Sample lists not found. Please run 07_pca_group_visualization.R first.")
}

sample_lists <- readRDS(sample_lists_path)

# Get R0 and R1 tumor samples
r0_tumor_samples <- sample_lists$R0$tumor
r1_tumor_samples <- sample_lists$R1$tumor

cat(sprintf("  R0 tumor samples: %d\n", length(r0_tumor_samples)))
cat(sprintf("  R1 tumor samples: %d\n", length(r1_tumor_samples)))

# --------------------------------------------------------------------------
# 1.3 Extract raw counts for DEGs
# --------------------------------------------------------------------------

cat("\n1.3 Extracting raw counts for DEGs...\n")

# Combine all samples
all_samples <- c(r0_tumor_samples, r1_tumor_samples)

# Get all significant genes
all_sig_genes <- unique(c(up_genes, down_genes))

# Check gene availability in SE
available_genes <- intersect(all_sig_genes, rownames(se))
missing_genes <- setdiff(all_sig_genes, rownames(se))

if (length(missing_genes) > 0) {
  cat(sprintf("  Warning: %d genes not found in count matrix\n", length(missing_genes)))
}

# Extract counts
raw_counts <- assay(se, "counts")[available_genes, all_samples]

cat(sprintf("  Count matrix: %d genes × %d samples\n", nrow(raw_counts), ncol(raw_counts)))

# Update gene lists to available genes only
up_genes <- intersect(up_genes, available_genes)
down_genes <- intersect(down_genes, available_genes)

cat(sprintf("  Available upregulated: %d\n", length(up_genes)))
cat(sprintf("  Available downregulated: %d\n", length(down_genes)))

# --------------------------------------------------------------------------
# 1.4 Calculate TPM
# --------------------------------------------------------------------------

cat("\n1.4 Calculating TPM values...\n")

# Get gene lengths from rowData (added by 03_prepare_counts.R)
if (!"gene_length" %in% colnames(rowData(se))) {
  stop("gene_length not found in rowData. Please ensure 03_prepare_counts.R has been run with gene length addition.")
}

# Extract gene lengths for available genes
gene_lengths <- rowData(se)[available_genes, "gene_length"]

# Check for any NA values
na_lengths <- is.na(gene_lengths)
if (any(na_lengths)) {
  cat(sprintf("  Warning: %d genes have NA lengths, removing them...\n", sum(na_lengths)))
  
  # Remove genes with NA lengths
  available_genes <- available_genes[!na_lengths]
  raw_counts <- raw_counts[available_genes, ]  # Use gene names directly
  gene_lengths <- gene_lengths[!na_lengths]
  
  # Update gene lists
  up_genes <- intersect(up_genes, available_genes)
  down_genes <- intersect(down_genes, available_genes)
  
  cat(sprintf("  Genes after NA removal: %d\n", length(available_genes)))
  cat(sprintf("  Updated upregulated: %d\n", length(up_genes)))
  cat(sprintf("  Updated downregulated: %d\n", length(down_genes)))
}

cat(sprintf("  Using GENCODE v36 gene lengths for %d genes\n", length(gene_lengths)))

# Calculate TPM
calculate_tpm <- function(counts, lengths) {
  # RPK (Reads Per Kilobase)
  rpk <- counts / (lengths / 1000)
  
  # TPM (Transcripts Per Million)
  scaling_factor <- colSums(rpk) / 1e6
  tpm <- sweep(rpk, 2, scaling_factor, "/")
  
  return(tpm)
}

tpm_matrix <- calculate_tpm(raw_counts, gene_lengths)

# Verify TPM sums
tpm_sums <- colSums(tpm_matrix)
cat(sprintf("  TPM sum range: %.2f - %.2f (expected: ~1M)\n", 
            min(tpm_sums), max(tpm_sums)))

# --------------------------------------------------------------------------
# 1.5 Filter genes with zero expression
# --------------------------------------------------------------------------

cat("\n1.5 Filtering zero expression genes...\n")

# Identify genes with TPM > 0 in all samples
non_zero_genes <- rownames(tpm_matrix)[rowSums(tpm_matrix > CONFIG$MIN_TPM) == ncol(tpm_matrix)]

cat(sprintf("  Genes before filtering: %d\n", nrow(tpm_matrix)))
cat(sprintf("  Genes after filtering (TPM > 0 in all): %d\n", length(non_zero_genes)))

# Update gene lists - intersect with available_genes to maintain consistency
up_genes_filtered <- intersect(up_genes, non_zero_genes)
down_genes_filtered <- intersect(down_genes, non_zero_genes)

cat(sprintf("  Filtered upregulated: %d (removed %d)\n", 
            length(up_genes_filtered), 
            length(up_genes) - length(up_genes_filtered)))
cat(sprintf("  Filtered downregulated: %d (removed %d)\n", 
            length(down_genes_filtered),
            length(down_genes) - length(down_genes_filtered)))

# --------------------------------------------------------------------------
# 1.6 Prepare final data structures
# --------------------------------------------------------------------------

cat("\n1.6 Preparing final data structures...\n")

# Create filtered TPM matrices
tpm_r0 <- tpm_matrix[non_zero_genes, r0_tumor_samples]
tpm_r1 <- tpm_matrix[non_zero_genes, r1_tumor_samples]

# Create combined data structure for Phase 2
reo_data <- list(
  # TPM matrices
  tpm_r0 = tpm_r0,
  tpm_r1 = tpm_r1,
  
  # Gene lists
  up_genes = up_genes_filtered,
  down_genes = down_genes_filtered,
  
  # Sample info
  r0_samples = r0_tumor_samples,
  r1_samples = r1_tumor_samples,
  
  # Metadata
  total_pairs_possible = length(up_genes_filtered) * length(down_genes_filtered),
  
  # Original DEG info for reference
  deg_info = sig_degs %>%
    dplyr::filter(gene_id %in% non_zero_genes) %>%
    .[, c("gene_id", "log2FC", "pvalue", "qvalue")]
)

cat(sprintf("\nData preparation complete:\n"))
cat(sprintf("  R0 samples: %d\n", length(r0_tumor_samples)))
cat(sprintf("  R1 samples: %d\n", length(r1_tumor_samples)))
cat(sprintf("  Up genes (filtered): %d\n", length(up_genes_filtered)))
cat(sprintf("  Down genes (filtered): %d\n", length(down_genes_filtered)))
cat(sprintf("  Potential REO pairs: %s\n", format(reo_data$total_pairs_possible, big.mark=",")))

# --------------------------------------------------------------------------
# 1.7 Save intermediate results
# --------------------------------------------------------------------------

cat("\n1.7 Saving intermediate results...\n")

output_path <- paste0(paths$processed, "reo_phase1_data.rds")
saveRDS(reo_data, output_path)
cat(sprintf("  Saved to: %s\n", basename(output_path)))

# --------------------------------------------------------------------------
# Summary statistics for verification
# --------------------------------------------------------------------------

cat("\n--- Phase 1 Summary ---\n")

# TPM distribution summary
tpm_all <- tpm_matrix[non_zero_genes, ]
cat("\nTPM distribution (all samples):\n")
cat(sprintf("  Median: %.2f\n", median(as.vector(tpm_all))))
cat(sprintf("  Mean: %.2f\n", mean(as.vector(tpm_all))))
cat(sprintf("  Q1-Q3: %.2f - %.2f\n", 
            quantile(as.vector(tpm_all), 0.25),
            quantile(as.vector(tpm_all), 0.75)))

# Gene-wise TPM ranges
cat("\nGene-wise median TPM ranges:\n")
gene_medians <- apply(tpm_all, 1, median)
cat(sprintf("  Min median TPM: %.4f\n", min(gene_medians)))
cat(sprintf("  Max median TPM: %.2f\n", max(gene_medians)))
cat("  Note: TPM calculated using GENCODE v36 exon union lengths\n")

cat("\n=== Phase 1 Complete ===\n")

# ============================================================================
# Phase 2: REO Pair Generation and Matrix Construction
# ============================================================================

# Uncomment below to run Phase 2
# source("analysis_v7/setup.R")
# reo_data <- readRDS(paste0(paths$processed, "reo_phase1_data.rds"))

if (exists("reo_data")) {
  
  cat("\n=== Phase 2: REO Pair Generation ===\n")
  
  # --------------------------------------------------------------------------
  # 2.1 Generate all up/down pairs
  # --------------------------------------------------------------------------
  
  cat("\n2.1 Generating up/down gene pairs...\n")
  
  # Create pair combinations
  pair_grid <- expand.grid(
    gene_up = reo_data$up_genes,
    gene_down = reo_data$down_genes,
    stringsAsFactors = FALSE
  )
  
  n_pairs <- nrow(pair_grid)
  cat(sprintf("  Total pairs generated: %s\n", format(n_pairs, big.mark = ",")))
  
  # Add DEG fold change information for each gene
  # Create lookup tables for FC values
  deg_fc_lookup <- setNames(reo_data$deg_info$log2FC, reo_data$deg_info$gene_id)
  
  # Add FC information with clear naming
  pair_grid$gene_up_deg_fc <- deg_fc_lookup[pair_grid$gene_up]
  pair_grid$gene_down_deg_fc <- deg_fc_lookup[pair_grid$gene_down]
  
  cat("  Added DEG fold change information (R1 vs R0)\n")
  
  # --------------------------------------------------------------------------
  # 2.2 Calculate REO ratios for all pairs
  # --------------------------------------------------------------------------
  
  cat("\n2.2 Calculating REO ratios (log2 scale)...\n")
  
  # Extract TPM matrices
  tpm_r0 <- reo_data$tpm_r0
  tpm_r1 <- reo_data$tpm_r1
  
  # Calculate REO ratios (vectorized for efficiency)
  cat("  Calculating R0 group ratios...\n")
  
  # Vectorized calculation for R0
  up_tpm_r0 <- tpm_r0[pair_grid$gene_up, ]
  down_tpm_r0 <- tpm_r0[pair_grid$gene_down, ]
  reo_matrix_r0 <- log2(up_tpm_r0 / down_tpm_r0)
  
  # Ensure correct dimensions
  if (!is.matrix(reo_matrix_r0)) {
    reo_matrix_r0 <- matrix(reo_matrix_r0, nrow = n_pairs, ncol = length(reo_data$r0_samples))
  }
  colnames(reo_matrix_r0) <- reo_data$r0_samples
  
  cat("  Calculating R1 group ratios...\n")
  
  # Vectorized calculation for R1
  up_tpm_r1 <- tpm_r1[pair_grid$gene_up, ]
  down_tpm_r1 <- tpm_r1[pair_grid$gene_down, ]
  reo_matrix_r1 <- log2(up_tpm_r1 / down_tpm_r1)
  
  # Ensure correct dimensions
  if (!is.matrix(reo_matrix_r1)) {
    reo_matrix_r1 <- matrix(reo_matrix_r1, nrow = n_pairs, ncol = length(reo_data$r1_samples))
  }
  colnames(reo_matrix_r1) <- reo_data$r1_samples
  
  cat(sprintf("  REO matrices created: %d pairs × 24 samples\n", n_pairs))
  
  # --------------------------------------------------------------------------
  # 2.3 Calculate group-wise statistics
  # --------------------------------------------------------------------------
  
  cat("\n2.3 Calculating group-wise statistics...\n")
  
  # Calculate medians for each pair
  pair_grid$median_reo_r0 <- apply(reo_matrix_r0, 1, median)
  pair_grid$median_reo_r1 <- apply(reo_matrix_r1, 1, median)
  
  # Calculate MAD (Median Absolute Deviation) for robustness assessment
  pair_grid$mad_reo_r0 <- apply(reo_matrix_r0, 1, mad)
  pair_grid$mad_reo_r1 <- apply(reo_matrix_r1, 1, mad)
  
  # Calculate difference in medians (potential reversal indicator)
  pair_grid$median_diff <- pair_grid$median_reo_r1 - pair_grid$median_reo_r0
  
  cat("  Group statistics calculated\n")
  
  # --------------------------------------------------------------------------
  # 2.4 Summary statistics
  # --------------------------------------------------------------------------
  
  cat("\n2.4 Summary of REO pairs:\n")
  
  # Distribution of median REO values
  cat("\nR0 group median REO distribution:\n")
  r0_med_summary <- summary(pair_grid$median_reo_r0)
  print(r0_med_summary)
  
  cat("\nR1 group median REO distribution:\n")
  r1_med_summary <- summary(pair_grid$median_reo_r1)
  print(r1_med_summary)
  
  cat("\nMedian difference (R1 - R0) distribution:\n")
  diff_summary <- summary(pair_grid$median_diff)
  print(diff_summary)
  
  # Count potential reversals (sign change in median)
  potential_reversals <- sum(
    (pair_grid$median_reo_r0 > 0 & pair_grid$median_reo_r1 < 0) |
      (pair_grid$median_reo_r0 < 0 & pair_grid$median_reo_r1 > 0)
  )
  
  cat(sprintf("\nPotential reversals (sign change in median): %d (%.2f%%)\n",
              potential_reversals,
              potential_reversals / n_pairs * 100))
  
  # --------------------------------------------------------------------------
  # 2.5 Save Phase 2 results
  # --------------------------------------------------------------------------
  
  cat("\n2.5 Saving Phase 2 results...\n")
  
  reo_phase2_data <- list(
    # Pair information with statistics
    pair_info = pair_grid,
    
    # REO matrices (can be large)
    reo_matrix_r0 = reo_matrix_r0,
    reo_matrix_r1 = reo_matrix_r1,
    
    # Sample information
    r0_samples = reo_data$r0_samples,
    r1_samples = reo_data$r1_samples,
    
    # Metadata
    n_pairs = n_pairs,
    n_up_genes = length(reo_data$up_genes),
    n_down_genes = length(reo_data$down_genes),
    date_created = Sys.Date()
  )
  
  output_path <- paste0(paths$processed, "reo_phase2_data.rds")
  saveRDS(reo_phase2_data, output_path)
  cat(sprintf("  Saved to: %s\n", basename(output_path)))
  
  # Save size information
  file_size_mb <- file.info(output_path)$size / 1024^2
  cat(sprintf("  File size: %.1f MB\n", file_size_mb))
  
  # --------------------------------------------------------------------------
  # Phase 2 Summary
  # --------------------------------------------------------------------------
  
  cat("\n=== Phase 2 Complete ===\n")
  cat("Summary:\n")
  cat(sprintf("  Total REO pairs: %s\n", format(n_pairs, big.mark = ",")))
  cat(sprintf("  REO matrices: %d pairs × %d samples (R0: 11, R1: 13)\n", 
              n_pairs, length(reo_data$r0_samples) + length(reo_data$r1_samples)))
  cat(sprintf("  Memory usage: ~%.1f MB per matrix\n", 
              object.size(reo_matrix_r0) / 1024^2))
  cat("\nNext: Run Phase 3 for filtering and statistical testing\n")
  
  # Clean up large objects
  rm(reo_matrix_r0, reo_matrix_r1)
  gc()
  
} else {
  cat("\nTo run Phase 2, first load Phase 1 data:\n")
  cat("  reo_data <- readRDS(paste0(paths$processed, 'reo_phase1_data.rds'))\n")
  cat("  Then re-run this script\n")
}