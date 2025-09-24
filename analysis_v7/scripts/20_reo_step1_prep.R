# 20_reo_step1_prep.R - REO Step 1: Input and Preprocessing
# Purpose: Load DEG results, calculate TPM, prepare zero-free expression matrices
# Input: thyr_deg_results.rds, thyr_case_master_stage2_filtered.rds
# Output: reo_step1_data.rds
# Version: v1.0
# Date: 2025-01-26

source("analysis_v7/setup.R")
source("analysis_v7/config/reo_config.R")  # Step1のみCONFIGを読む

cat("\n=== REO STEP 1: Input and Preprocessing ===\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
})

# --------------------------------------------------------------------------
# 1.1 Load DEG results from 09_deg_analysis.R
# --------------------------------------------------------------------------

cat("--- 1.1 Loading DEG results ---\n")

deg_results_path <- paste0(paths$processed, "thyr_deg_results.rds")
if (!file.exists(deg_results_path)) {
  stop("ERROR: thyr_deg_results.rds not found. Please run 09_deg_analysis.R first.")
}

deg_output <- readRDS(deg_results_path)
cat("  DEG results loaded successfully\n")

# Extract R0_vs_R1_tumor results
r0r1_tumor <- deg_output$deg_results$R0_vs_R1_tumor
if (is.null(r0r1_tumor)) {
  stop("ERROR: R0_vs_R1_tumor results not found in DEG output.")
}

# Get DEG dataframe
deg_df <- r0r1_tumor$deg_summary$results_df

# Get significant DEGs
sig_degs <- deg_df[deg_df$significant == TRUE, ]

cat(sprintf("  R0_vs_R1_tumor DEGs: %d genes (q < 0.05)\n", nrow(sig_degs)))
cat(sprintf("    Upregulated: %d\n", sum(sig_degs$log2FC > 0)))
cat(sprintf("    Downregulated: %d\n", sum(sig_degs$log2FC < 0)))

# --------------------------------------------------------------------------
# 1.2 Load case master and identify RET samples
# --------------------------------------------------------------------------

cat("\n--- 1.2 Identifying RET fusion samples ---\n")

# Load case master
case_master_path <- paste0(paths$processed, "thyr_case_master_stage2_filtered.rds")
if (!file.exists(case_master_path)) {
  stop("ERROR: Case master not found. Please run 06_purity_analysis.R first.")
}

case_master <- readRDS(case_master_path)

# Filter for RET fusion cases with high purity
ret_cases <- case_master %>%
  filter(driver == "RET" & 
           has_outlier_tumor == 0 & 
           has_outlier_normal == 0 & 
           low_purity == 0)

cat(sprintf("  RET fusion cases (high purity): %d\n", nrow(ret_cases)))

# Separate R0 and R1 cases
r0_cases <- ret_cases %>% filter(group == "R0")
r1_cases <- ret_cases %>% filter(group == "R1")

cat(sprintf("  R0 cases (POC=0%%): %d\n", nrow(r0_cases)))
cat(sprintf("  R1 cases (POC≥66.6%%): %d\n", nrow(r1_cases)))

# Check sample counts
if (nrow(r0_cases) < 5 || nrow(r1_cases) < 5) {
  stop("ERROR: Insufficient samples. Need at least 5 samples per group.")
}

# --------------------------------------------------------------------------
# 1.3 Get tumor sample IDs (using paired samples only)
# --------------------------------------------------------------------------

cat("\n--- 1.3 Getting tumor sample IDs ---\n")

# Load sample metadata
se_path <- paste0(paths$processed, "thyr_se_strand2_nonzero.rds")
if (!file.exists(se_path)) {
  stop("ERROR: SummarizedExperiment not found. Please run 03_prepare_counts.R first.")
}

se <- readRDS(se_path)
sample_metadata <- as.data.frame(colData(se))

# Load paired sample lists from previous analysis
sample_lists_path <- paste0(paths$processed, "analysis_sample_lists.rds")
if (file.exists(sample_lists_path)) {
  # Use pre-computed paired sample lists
  sample_lists <- readRDS(sample_lists_path)
  
  # Get RET group paired tumor samples (correct structure)
  r0_samples <- sample_lists$R0$tumor
  r1_samples <- sample_lists$R1$tumor
  
  cat("  Using pre-computed paired sample lists\n")
} else {
  # Fallback: compute samples from case master and RET filter
  cat("  Computing RET tumor samples from case master\n")
  
  # Get RET cases only
  ret_r0_cases <- r0_cases$case_submitter_id
  ret_r1_cases <- r1_cases$case_submitter_id
  
  # Find tumor samples for RET cases (with _merged suffix)
  r0_samples <- sample_metadata %>%
    filter(case_submitter_id %in% ret_r0_cases &
             sample_type == "Primary Tumor" &
             grepl("_merged", sample_submitter_id)) %>%
    pull(sample_submitter_id)
  
  r1_samples <- sample_metadata %>%
    filter(case_submitter_id %in% ret_r1_cases &
             sample_type == "Primary Tumor" &
             grepl("_merged", sample_submitter_id)) %>%
    pull(sample_submitter_id)
}

cat(sprintf("  R0 tumor samples: %d\n", length(r0_samples)))
cat(sprintf("  R1 tumor samples: %d\n", length(r1_samples)))

# Verify we have the expected counts
expected_r0 <- 11
expected_r1 <- 13

if (length(r0_samples) != expected_r0) {
  warning(sprintf("  WARNING: Expected %d R0 samples, found %d", expected_r0, length(r0_samples)))
}
if (length(r1_samples) != expected_r1) {
  warning(sprintf("  WARNING: Expected %d R1 samples, found %d", expected_r1, length(r1_samples)))
}

# --------------------------------------------------------------------------
# 1.4 Calculate TPM from count data
# --------------------------------------------------------------------------

cat("\n--- 1.4 Calculating TPM ---\n")

# Extract raw counts
raw_counts <- assay(se, "counts")

# Check for gene_length in rowData
if (!"gene_length" %in% colnames(rowData(se))) {
  stop("ERROR: gene_length not found in rowData. Please check 03_prepare_counts.R")
}

# Get gene lengths
gene_lengths <- rowData(se)$gene_length
names(gene_lengths) <- rownames(se)

# Validate gene lengths
cat("  Validating gene lengths...\n")
invalid_lengths <- which(is.na(gene_lengths) | gene_lengths <= 0)
if (length(invalid_lengths) > 0) {
  cat(sprintf("  WARNING: Found %d genes with invalid lengths (NA or ≤0)\n", 
              length(invalid_lengths)))
  
  # Remove genes with invalid lengths
  valid_genes <- setdiff(rownames(raw_counts), rownames(raw_counts)[invalid_lengths])
  raw_counts <- raw_counts[valid_genes, ]
  gene_lengths <- gene_lengths[valid_genes]
  
  cat(sprintf("  Proceeding with %d genes with valid lengths\n", length(valid_genes)))
}

# Function to calculate TPM
calculate_tpm <- function(counts, lengths) {
  # Final safety check
  if (any(lengths <= 0, na.rm = TRUE)) {
    stop("ERROR: Invalid gene lengths detected in TPM calculation")
  }
  
  # Ensure matching dimensions
  if (length(lengths) != nrow(counts)) {
    stop("ERROR: Gene length vector doesn't match count matrix rows")
  }
  
  # Calculate RPK (Reads Per Kilobase)
  rpk <- counts / (lengths / 1000)
  
  # Check for non-finite values
  if (any(!is.finite(rpk))) {
    n_inf <- sum(is.infinite(rpk))
    n_nan <- sum(is.nan(rpk))
    stop(sprintf("ERROR: Non-finite values in RPK calculation (Inf: %d, NaN: %d)", 
                 n_inf, n_nan))
  }
  
  # Calculate scaling factors
  scaling_factor <- colSums(rpk) / 1e6
  
  # Calculate TPM
  tpm <- sweep(rpk, 2, scaling_factor, "/")
  
  return(tpm)
}

# Calculate TPM for all samples
tpm_matrix <- calculate_tpm(raw_counts, gene_lengths)

# Verify TPM sums (should be ~1 million)
tpm_sums <- colSums(tpm_matrix)
cat(sprintf("  TPM sum range: %.2f - %.2f\n", min(tpm_sums), max(tpm_sums)))

if (any(abs(tpm_sums - 1e6) > 1)) {
  warning("  WARNING: Some samples have TPM sums deviating from 1M")
}

# --------------------------------------------------------------------------
# 1.5 Filter for zero-free genes
# --------------------------------------------------------------------------

cat("\n--- 1.5 Filtering for zero-free genes ---\n")

# Combine R0 and R1 samples
analysis_samples <- c(r0_samples, r1_samples)

# Extract TPM for analysis samples
tpm_analysis <- tpm_matrix[, analysis_samples]

# Find genes with TPM > 0 in all samples
zero_free_genes <- rownames(tpm_analysis)[rowSums(tpm_analysis > 0) == ncol(tpm_analysis)]

cat(sprintf("  Total genes: %d\n", nrow(tpm_analysis)))
cat(sprintf("  Zero-free genes: %d\n", length(zero_free_genes)))
cat(sprintf("  Genes with zeros removed: %d\n", nrow(tpm_analysis) - length(zero_free_genes)))

# --------------------------------------------------------------------------
# 1.6 Apply log2 transformation (no pseudocount)
# --------------------------------------------------------------------------

cat("\n--- 1.6 Log2 transformation ---\n")

# Filter TPM matrix to zero-free genes
tpm_zero_free <- tpm_analysis[zero_free_genes, ]

# Apply log2 transformation
log2_tpm <- log2(tpm_zero_free)

# Check for any NA, NaN, or Inf values
if (any(!is.finite(log2_tpm))) {
  stop("ERROR: Non-finite values detected after log2 transformation")
}

cat("  Log2 transformation completed successfully\n")
cat(sprintf("  Matrix dimensions: %d genes × %d samples\n", 
            nrow(log2_tpm), ncol(log2_tpm)))

# --------------------------------------------------------------------------
# 1.7 Split into R0 and R1 matrices
# --------------------------------------------------------------------------

cat("\n--- 1.7 Creating group-specific matrices ---\n")

# Create R0 and R1 matrices
log2_tpm_r0 <- log2_tpm[, r0_samples]
log2_tpm_r1 <- log2_tpm[, r1_samples]

cat(sprintf("  R0 matrix: %d genes × %d samples\n", 
            nrow(log2_tpm_r0), ncol(log2_tpm_r0)))
cat(sprintf("  R1 matrix: %d genes × %d samples\n", 
            nrow(log2_tpm_r1), ncol(log2_tpm_r1)))

# --------------------------------------------------------------------------
# 1.8 Summary statistics
# --------------------------------------------------------------------------

cat("\n--- 1.8 Summary statistics ---\n")

# R0 group statistics
r0_means <- rowMeans(log2_tpm_r0)
cat("\nR0 group (log2 TPM):\n")
cat(sprintf("  Mean expression range: %.2f to %.2f\n", 
            min(r0_means), max(r0_means)))
cat(sprintf("  Median of means: %.2f\n", median(r0_means)))

# R1 group statistics
r1_means <- rowMeans(log2_tpm_r1)
cat("\nR1 group (log2 TPM):\n")
cat(sprintf("  Mean expression range: %.2f to %.2f\n", 
            min(r1_means), max(r1_means)))
cat(sprintf("  Median of means: %.2f\n", median(r1_means)))

# --------------------------------------------------------------------------
# 1.9 Save Step 1 results
# --------------------------------------------------------------------------

cat("\n--- 1.9 Saving Step 1 results ---\n")

# Check how many DEGs are in zero-free genes
degs_in_zero_free <- intersect(sig_degs$gene_id, zero_free_genes)
cat(sprintf("  DEGs in zero-free genes: %d / %d (%.1f%%)\n", 
            length(degs_in_zero_free), nrow(sig_degs),
            length(degs_in_zero_free) / nrow(sig_degs) * 100))

# Separate up/down DEGs that are zero-free
up_degs_zero_free <- intersect(
  sig_degs$gene_id[sig_degs$log2FC > 0],
  zero_free_genes
)
down_degs_zero_free <- intersect(
  sig_degs$gene_id[sig_degs$log2FC < 0], 
  zero_free_genes
)

cat(sprintf("    Up-regulated in zero-free: %d\n", length(up_degs_zero_free)))
cat(sprintf("    Down-regulated in zero-free: %d\n", length(down_degs_zero_free)))

step1_data <- list(
  # Configuration (スナップショット)
  config = CONFIG,
  
  # Sample information
  r0_samples = r0_samples,
  r1_samples = r1_samples,
  n_r0 = length(r0_samples),
  n_r1 = length(r1_samples),
  
  # Gene information
  zero_free_genes = zero_free_genes,
  degs_in_zero_free = degs_in_zero_free,
  up_degs_zero_free = up_degs_zero_free,
  down_degs_zero_free = down_degs_zero_free,
  
  # Expression matrices (log2 TPM)
  log2_tpm_r0 = log2_tpm_r0,
  log2_tpm_r1 = log2_tpm_r1,
  
  # DEG information
  deg_df = deg_df,
  sig_degs = sig_degs,
  
  # Metadata
  timestamp = Sys.time(),
  step = "Step 1 Complete"
)

saveRDS(step1_data, paste0(paths$processed, "reo_step1_data.rds"))
cat("  Step 1 data saved to: reo_step1_data.rds\n")

# =============================================================================
# STEP 1 SUMMARY (固定フォーマット)
# =============================================================================

cat("\n=== STEP 1 SUMMARY ===\n")
cat(sprintf("R0 samples: %d\n", length(r0_samples)))
cat(sprintf("R1 samples: %d\n", length(r1_samples)))
cat(sprintf("Zero-free genes: %d\n", length(zero_free_genes)))
cat(sprintf("DEGs in zero-free: %d (%.1f%%)\n", 
            length(degs_in_zero_free),
            length(degs_in_zero_free) / nrow(sig_degs) * 100))
cat(sprintf("Timestamp: %s\n", format(step1_data$timestamp)))