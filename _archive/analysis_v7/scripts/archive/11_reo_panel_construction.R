# 11_reo_panel_construction.R - REO Panel Construction following Strategy Document
# Purpose: Construct REO panel for radiation exposure detection
# Method: 6-step process as defined in REO_Panel_Strategy.md
# Version: v1.0
# Date: 2025-01-25

source("analysis_v7/setup.R")

cat("\n=== REO Panel Construction (v1.0) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Following REO_Panel_Strategy.md specifications\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
})

# ============================================================================
# Configuration (from REO_Panel_Strategy.md)
# ============================================================================

CONFIG <- list(
  # Basic parameters
  tau_strength = log2(1.5),        # Strength criterion
  dead_zone = log2(1.2),           # Dead zone threshold
  trim_low = 0.07,                 # Lower trim rate (ADJUSTED from 0.10)
  trim_high = 0.10,                # Upper trim rate
  M = 12,                          # Maximum pairs per gene (ADJUSTED from 10)
  N = 10,                          # Panel size (default)
  T_ratio = 0.4,                   # Decision threshold ratio
  CI_method = "Wilson",            # CI calculation method (no continuity correction)
  
  # Stopping condition parameters
  min_candidates_abs = 100,        # Minimum candidate count (absolute)
  min_candidates_ratio = 0.02,     # Minimum candidate ratio (relative to initial)
  
  # Processing parameters
  VERBOSE = TRUE                   # Verbose output
)

# Log configuration
cat("=== Configuration Parameters ===\n")
cat("  NOTE: trim_low adjusted to 0.07, M adjusted to 12 due to candidate shortage\n")
for(param in names(CONFIG)) {
  if(is.numeric(CONFIG[[param]])) {
    cat(sprintf("  %s: %.4f\n", param, CONFIG[[param]]))
  } else {
    cat(sprintf("  %s: %s\n", param, CONFIG[[param]]))
  }
}
cat("\n")

# ============================================================================
# STEP 1: Input and Preprocessing
# ============================================================================

cat("=== STEP 1: Input and Preprocessing ===\n\n")

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

# Display actual column names for transparency
cat("  DEG result columns: ", paste(colnames(deg_df)[1:min(8, ncol(deg_df))], collapse=", "))
if (ncol(deg_df) > 8) cat(", ...")
cat("\n")

# Get significant DEGs
sig_degs <- deg_df[deg_df$significant == TRUE, ]

# Additional check: ensure we have data
if (nrow(sig_degs) == 0) {
  warning("WARNING: No significant DEGs found (significant == TRUE)")
  cat("  Check: unique values in 'significant' column:", paste(unique(deg_df$significant), collapse=", "), "\n")
}

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
cat(sprintf("  R1 cases (POCâ‰¥66.6%%): %d\n", nrow(r1_cases)))

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
  cat(sprintf("  WARNING: Found %d genes with invalid lengths (NA or â‰¤0)\n", 
              length(invalid_lengths)))
  
  # Show examples
  if (length(invalid_lengths) <= 5) {
    for (i in invalid_lengths) {
      cat(sprintf("    %s: length = %s\n", 
                  names(gene_lengths)[i], 
                  ifelse(is.na(gene_lengths[i]), "NA", as.character(gene_lengths[i]))))
    }
  } else {
    cat(sprintf("    First 5: %s\n", 
                paste(names(gene_lengths)[invalid_lengths[1:5]], collapse=", ")))
  }
  
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
cat(sprintf("  Matrix dimensions: %d genes Ã— %d samples\n", 
            nrow(log2_tpm), ncol(log2_tpm)))

# --------------------------------------------------------------------------
# 1.7 Split into R0 and R1 matrices
# --------------------------------------------------------------------------

cat("\n--- 1.7 Creating group-specific matrices ---\n")

# Create R0 and R1 matrices
log2_tpm_r0 <- log2_tpm[, r0_samples]
log2_tpm_r1 <- log2_tpm[, r1_samples]

cat(sprintf("  R0 matrix: %d genes Ã— %d samples\n", 
            nrow(log2_tpm_r0), ncol(log2_tpm_r0)))
cat(sprintf("  R1 matrix: %d genes Ã— %d samples\n", 
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
# Note: zero_free_genes are gene IDs (rownames), sig_degs$gene_id are also gene IDs
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
  # Configuration
  config = CONFIG,
  
  # Sample information
  r0_samples = r0_samples,
  r1_samples = r1_samples,
  n_r0 = length(r0_samples),
  n_r1 = length(r1_samples),
  
  # Gene information
  zero_free_genes = zero_free_genes,
  degs_in_zero_free = degs_in_zero_free,  # All DEGs that are zero-free
  up_degs_zero_free = up_degs_zero_free,  # Up-regulated DEGs that are zero-free
  down_degs_zero_free = down_degs_zero_free,  # Down-regulated DEGs that are zero-free
  
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

cat("\n=== STEP 1 COMPLETE ===\n")
cat("Ready for Step 2: Candidate Pair Generation\n")

# Final check
cat("\n--- Final Status ---\n")
cat(sprintf("  âœ“ R0 samples: %d (expected: 11)\n", length(r0_samples)))
cat(sprintf("  âœ“ R1 samples: %d (expected: 13)\n", length(r1_samples)))
cat(sprintf("  âœ“ Zero-free genes: %d\n", length(zero_free_genes)))
cat(sprintf("  âœ“ DEGs total: %d\n", nrow(sig_degs)))
cat(sprintf("  âœ“ DEGs in zero-free: %d (%.1f%%)\n", 
            length(degs_in_zero_free),
            length(degs_in_zero_free) / nrow(sig_degs) * 100))
cat("  âœ“ Data saved successfully\n")
cat("\nNote: DEG filtering will be applied in Step 2 for candidate pair generation\n")

# ============================================================================
# STEP 2: Candidate Pair Generation (RET-limited)
# ============================================================================

cat("\n\n=== STEP 2: Candidate Pair Generation ===\n\n")

# --------------------------------------------------------------------------
# 2.1 Load Step 1 data
# --------------------------------------------------------------------------

cat("--- 2.1 Loading Step 1 data ---\n")

step1_data <- readRDS(paste0(paths$processed, "reo_step1_data.rds"))
cat("  Step 1 data loaded successfully\n")

# Extract necessary data
up_degs <- step1_data$up_degs_zero_free
down_degs <- step1_data$down_degs_zero_free
log2_tpm_r0 <- step1_data$log2_tpm_r0
log2_tpm_r1 <- step1_data$log2_tpm_r1

cat(sprintf("  Up-regulated DEGs (zero-free): %d\n", length(up_degs)))
cat(sprintf("  Down-regulated DEGs (zero-free): %d\n", length(down_degs)))
cat(sprintf("  Initial pairs possible: %s\n", 
            format(length(up_degs) * length(down_degs), big.mark = ",")))

# Store initial pair count for ratio calculation
initial_pairs_raw <- length(up_degs) * length(down_degs)

# --------------------------------------------------------------------------
# 2.2 Calculate gene median TPM for filtering
# --------------------------------------------------------------------------

cat("\n--- 2.2 Calculating median TPM across R0+R1 ---\n")

# Combine R0 and R1 log2 TPM matrices
all_samples_log2tpm <- cbind(log2_tpm_r0, log2_tpm_r1)
cat(sprintf("  Combined matrix: %d genes Ã— %d samples\n", 
            nrow(all_samples_log2tpm), ncol(all_samples_log2tpm)))

# Calculate median log2 TPM for each gene across all samples
gene_median_log2tpm <- apply(all_samples_log2tpm, 1, median)

# Filter to DEGs only
up_medians <- gene_median_log2tpm[up_degs]
down_medians <- gene_median_log2tpm[down_degs]

cat("  Median log2 TPM calculated for all genes\n")

# --------------------------------------------------------------------------
# 2.3 Apply expression level filter (trim upper/lower)
# --------------------------------------------------------------------------

cat("\n--- 2.3 Applying expression level filter ---\n")

# Use CONFIG values for trimming
trim_low <- CONFIG$trim_low    # Now 0.08 (adjusted)
trim_high <- CONFIG$trim_high  # Still 0.10

cat(sprintf("  Trim settings: lower %.0f%%, upper %.0f%%\n", 
            trim_low * 100, trim_high * 100))

# For up-regulated genes
up_quantiles <- quantile(up_medians, probs = c(trim_low, 1 - trim_high))
up_keep <- up_medians >= up_quantiles[1] & up_medians <= up_quantiles[2]
up_filtered <- names(up_medians)[up_keep]

cat(sprintf("  Up-regulated genes:\n"))
cat(sprintf("    Before trim: %d\n", length(up_medians)))
cat(sprintf("    Lower %.0f%% cutoff: %.2f\n", trim_low * 100, up_quantiles[1]))
cat(sprintf("    Upper %.0f%% cutoff: %.2f\n", trim_high * 100, up_quantiles[2]))
cat(sprintf("    After trim: %d (removed %d)\n", 
            length(up_filtered), length(up_medians) - length(up_filtered)))

# For down-regulated genes
down_quantiles <- quantile(down_medians, probs = c(trim_low, 1 - trim_high))
down_keep <- down_medians >= down_quantiles[1] & down_medians <= down_quantiles[2]
down_filtered <- names(down_medians)[down_keep]

cat(sprintf("  Down-regulated genes:\n"))
cat(sprintf("    Before trim: %d\n", length(down_medians)))
cat(sprintf("    Lower %.0f%% cutoff: %.2f\n", trim_low * 100, down_quantiles[1]))
cat(sprintf("    Upper %.0f%% cutoff: %.2f\n", trim_high * 100, down_quantiles[2]))
cat(sprintf("    After trim: %d (removed %d)\n", 
            length(down_filtered), length(down_medians) - length(down_filtered)))

# Store for logging
initial_pairs_after_trim <- length(up_filtered) * length(down_filtered)
cat(sprintf("\n  Pairs after trim: %s\n", 
            format(initial_pairs_after_trim, big.mark = ",")))

# --------------------------------------------------------------------------
# 2.4 Apply one-sided multiplicity limit (M=10)
# --------------------------------------------------------------------------

cat("\n--- 2.4 Applying multiplicity limit (M=", CONFIG$M, ") ---\n", sep="")

# Get DEG ranks for selection priority
deg_df <- step1_data$deg_df
sig_degs <- step1_data$sig_degs

# Create priority ranking: |logFC| desc, then qvalue asc, then gene_id asc
sig_degs$abs_logFC <- abs(sig_degs$log2FC)
sig_degs <- sig_degs[order(-sig_degs$abs_logFC, sig_degs$qvalue, sig_degs$gene_id), ]

# Create rank lookup
gene_ranks <- setNames(1:nrow(sig_degs), sig_degs$gene_id)

# Function to select top M partners based on rank
select_top_partners <- function(genes, partner_genes, M = 10) {
  # For each gene, select top M partners
  selected_pairs <- list()
  genes_hitting_limit <- 0
  
  for (gene in genes) {
    # Get all possible partners
    partners <- partner_genes
    
    # Rank partners
    partner_ranks <- gene_ranks[partners]
    
    # Count how many hit the M limit
    if (length(partner_ranks) >= M) {
      genes_hitting_limit <- genes_hitting_limit + 1
    }
    
    # Select top M (lowest rank = highest priority)
    top_partners <- names(sort(partner_ranks)[1:min(M, length(partner_ranks))])
    
    # Store pairs with proper naming
    for (partner in top_partners) {
      pair_name <- paste(gene, partner, sep = "_")
      selected_pairs[[pair_name]] <- c(gene, partner)
    }
  }
  
  return(list(pairs = selected_pairs, hitting_limit = genes_hitting_limit))
}

# Apply M limit for up genes (each up gene pairs with at most M down genes)
up_result <- select_top_partners(up_filtered, down_filtered, M = CONFIG$M)
up_limited_pairs <- up_result$pairs  # Extract pairs from result
up_hitting_limit <- up_result$hitting_limit

# Apply M limit for down genes (each down gene pairs with at most M up genes)  
down_result <- select_top_partners(down_filtered, up_filtered, M = CONFIG$M)
down_limited_pairs <- down_result$pairs  # Extract pairs from result
down_hitting_limit <- down_result$hitting_limit

# Diagnostic information
cat(sprintf("  Up genes hitting M limit: %d / %d (%.1f%%)\n",
            up_hitting_limit, length(up_filtered),
            up_hitting_limit / length(up_filtered) * 100))
cat(sprintf("  Down genes hitting M limit: %d / %d (%.1f%%)\n",
            down_hitting_limit, length(down_filtered),
            down_hitting_limit / length(down_filtered) * 100))

# Combine and deduplicate
all_pairs_list <- c(up_limited_pairs, down_limited_pairs)

# Deduplicate while preserving structure
# unique() drops names, so we need a different approach
pair_names <- names(all_pairs_list)
unique_pair_names <- unique(pair_names)
unique_pairs <- all_pairs_list[unique_pair_names]

cat(sprintf("\n  Pairs before M limit: %s\n", 
            format(initial_pairs_after_trim, big.mark = ",")))
cat(sprintf("  Pairs after M=%d limit: %d\n", CONFIG$M, length(unique_pairs)))
cat(sprintf("  Reduction: %.1f%%\n", 
            (1 - length(unique_pairs) / initial_pairs_after_trim) * 100))

# Additional diagnostic
if (up_hitting_limit + down_hitting_limit == 0) {
  cat("\n  WARNING: No genes hitting M limit. Increasing M will have no effect!\n")
} else {
  potential_gain <- (up_hitting_limit + down_hitting_limit) * 2  # Rough estimate
  cat(sprintf("  Potential gain from M increase: ~%d pairs\n", potential_gain))
}

final_candidates <- length(unique_pairs)

# --------------------------------------------------------------------------
# 2.5 Check stopping conditions
# --------------------------------------------------------------------------

cat("\n--- 2.5 Checking stopping conditions ---\n")

# Calculate minimum candidates threshold
min_candidates_abs <- CONFIG$min_candidates_abs
min_candidates_ratio <- CONFIG$min_candidates_ratio
min_candidates_dynamic <- ceiling(min_candidates_ratio * initial_pairs_raw)
min_candidates_required <- max(min_candidates_abs, min_candidates_dynamic)

cat(sprintf("  Minimum candidates (absolute): %d\n", min_candidates_abs))
cat(sprintf("  Minimum candidates (%.1f%% of initial): %d\n", 
            min_candidates_ratio * 100, min_candidates_dynamic))
cat(sprintf("  Required threshold: %d\n", min_candidates_required))
cat(sprintf("  Actual candidates: %d\n", final_candidates))

if (final_candidates < min_candidates_required) {
  cat("\n  WARNING: Insufficient candidate pairs!\n")
  
  # Check if we already adjusted parameters
  if (CONFIG$trim_low != 0.10 || CONFIG$M != 10) {
    cat("  NOTE: Already adjusted parameters:\n")
    if (CONFIG$trim_low != 0.10) {
      cat(sprintf("    - trim_low: 0.10 â†’ %.2f\n", CONFIG$trim_low))
    }
    if (CONFIG$M != 10) {
      cat(sprintf("    - M: 10 â†’ %d\n", CONFIG$M))
    }
  }
  
  cat("\n  Consider further adjustments:\n")
  cat("    - Further reduce trim thresholds\n")
  cat("    - Further increase M limit\n")
  cat("    - Check DEG quality\n")
  stop(sprintf("STOPPING: Only %d pairs, need at least %d", 
               final_candidates, min_candidates_required))
} else {
  cat(sprintf("\n  âœ“ Sufficient candidates: %d â‰¥ %d\n", 
              final_candidates, min_candidates_required))
}

# --------------------------------------------------------------------------
# 2.6 Create pair matrix
# --------------------------------------------------------------------------

cat("\n--- 2.6 Creating pair matrix ---\n")

# Check structure
if (length(unique_pairs) == 0) {
  stop("ERROR: No unique pairs generated")
}

# Verify names exist
if (is.null(names(unique_pairs))) {
  stop("ERROR: Pair names lost during processing")
}

# Convert list to data frame - more explicit approach
pairs_df <- data.frame(
  pair_id = names(unique_pairs),
  gene_up = sapply(unique_pairs, function(x) x[1]),
  gene_down = sapply(unique_pairs, function(x) x[2]),
  stringsAsFactors = FALSE
)

# Verify structure
if (nrow(pairs_df) != length(unique_pairs)) {
  stop("ERROR: Data frame creation failed")
}

# Add median TPM info
pairs_df$up_median_log2tpm <- gene_median_log2tpm[pairs_df$gene_up]
pairs_df$down_median_log2tpm <- gene_median_log2tpm[pairs_df$gene_down]

# Add DEG info (logFC, qvalue)
up_deg_info <- sig_degs[match(pairs_df$gene_up, sig_degs$gene_id), ]
down_deg_info <- sig_degs[match(pairs_df$gene_down, sig_degs$gene_id), ]

pairs_df$up_logFC <- up_deg_info$log2FC
pairs_df$up_qvalue <- up_deg_info$qvalue
pairs_df$down_logFC <- down_deg_info$log2FC
pairs_df$down_qvalue <- down_deg_info$qvalue

cat(sprintf("  Pair matrix created: %d pairs\n", nrow(pairs_df)))
cat("  Columns: pair_id, gene_up, gene_down, expression levels, DEG stats\n")

# --------------------------------------------------------------------------
# 2.7 Save Step 2 results
# --------------------------------------------------------------------------

cat("\n--- 2.7 Saving Step 2 results ---\n")

step2_data <- list(
  # Configuration
  config = CONFIG,
  
  # Pair information
  pairs_df = pairs_df,
  n_pairs = nrow(pairs_df),
  
  # Gene lists after filtering
  up_genes_filtered = up_filtered,
  down_genes_filtered = down_filtered,
  
  # Expression data (passed through)
  log2_tpm_r0 = log2_tpm_r0,
  log2_tpm_r1 = log2_tpm_r1,
  r0_samples = step1_data$r0_samples,
  r1_samples = step1_data$r1_samples,
  
  # Logging info
  initial_pairs_raw = initial_pairs_raw,
  initial_pairs_after_trim = initial_pairs_after_trim,
  final_candidates = final_candidates,
  
  # Metadata
  timestamp = Sys.time(),
  step = "Step 2 Complete"
)

saveRDS(step2_data, paste0(paths$processed, "reo_step2_data.rds"))
cat("  Step 2 data saved to: reo_step2_data.rds\n")

cat("\n=== STEP 2 COMPLETE ===\n")
cat("Ready for Step 3: Adoption Criteria\n")

# Final summary
cat("\n--- Step 2 Summary ---\n")
cat(sprintf("  Initial pairs (up Ã— down): %s\n", 
            format(initial_pairs_raw, big.mark = ",")))
cat(sprintf("  After expression trim: %s\n", 
            format(initial_pairs_after_trim, big.mark = ",")))
cat(sprintf("  After M=%d limit: %d\n", CONFIG$M, final_candidates))
cat(sprintf("  Retention rate: %.1f%%\n", 
            final_candidates / initial_pairs_raw * 100))
cat("  âœ“ All criteria satisfied\n")

# ============================================================================
# STEP 3a: Adoption Criteria - Basic Judgment
# ============================================================================

cat("\n\n=== STEP 3a: Adoption Criteria (Basic Judgment) ===\n\n")

# --------------------------------------------------------------------------
# 3a.1 Load Step 2 data and calculate REO values
# --------------------------------------------------------------------------

cat("--- 3a.1 Loading data and calculating REO values ---\n")

step2_data <- readRDS(paste0(paths$processed, "reo_step2_data.rds"))
pairs_df <- step2_data$pairs_df
log2_tpm_r0 <- step2_data$log2_tpm_r0
log2_tpm_r1 <- step2_data$log2_tpm_r1

cat(sprintf("  Loaded %d candidate pairs\n", nrow(pairs_df)))
cat(sprintf("  R0 samples: %d, R1 samples: %d\n", 
            ncol(log2_tpm_r0), ncol(log2_tpm_r1)))

# Calculate REO values (r = log2(up/down) = log2(up) - log2(down))
# Since we already have log2 TPM, just subtract
cat("\n  Calculating REO values (r = log2(up/down))...\n")

# For R0 samples
reo_r0 <- matrix(NA, nrow = nrow(pairs_df), ncol = ncol(log2_tpm_r0))
colnames(reo_r0) <- colnames(log2_tpm_r0)
for (i in 1:nrow(pairs_df)) {
  up_values <- log2_tpm_r0[pairs_df$gene_up[i], ]
  down_values <- log2_tpm_r0[pairs_df$gene_down[i], ]
  reo_r0[i, ] <- up_values - down_values
}

# For R1 samples
reo_r1 <- matrix(NA, nrow = nrow(pairs_df), ncol = ncol(log2_tpm_r1))
colnames(reo_r1) <- colnames(log2_tpm_r1)
for (i in 1:nrow(pairs_df)) {
  up_values <- log2_tpm_r1[pairs_df$gene_up[i], ]
  down_values <- log2_tpm_r1[pairs_df$gene_down[i], ]
  reo_r1[i, ] <- up_values - down_values
}

cat(sprintf("  REO matrices created: %d pairs Ã— %d R0 + %d R1 samples\n",
            nrow(pairs_df), ncol(reo_r0), ncol(reo_r1)))

# --------------------------------------------------------------------------
# 3a.2 Apply five judgment criteria
# --------------------------------------------------------------------------

cat("\n--- 3a.2 Applying judgment criteria ---\n")

# Initialize results
judgment_results <- data.frame(
  pair_id = pairs_df$pair_id,
  # Criterion 1: R0 directionality
  r0_consistent_direction = NA,
  r0_direction = NA,
  r0_exceptions = NA,
  r0_valid_samples = NA,
  criterion1_pass = FALSE,
  # Criterion 2: Strength
  r0_a2_value = NA,
  r0_strong_samples = NA,
  criterion2_pass = FALSE,
  # Criterion 3: R1 reversal
  r1_reversal_rate = NA,
  r1_reversal_count = NA,
  r1_valid_samples = NA,
  r1_ci_lower = NA,
  criterion3_pass = FALSE,
  # Criterion 4: 100% reversal exclusion
  criterion4_pass = TRUE,  # Default true, set false if 100% reversal
  # Criterion 5: Dead zone (implicit in other criteria)
  # Overall
  all_criteria_pass = FALSE,
  stringsAsFactors = FALSE
)

# Process each pair
cat("  Processing criteria for each pair...\n")
pb <- txtProgressBar(min = 0, max = nrow(pairs_df), style = 3)

for (i in 1:nrow(pairs_df)) {
  # Get REO values for this pair
  r0_values <- reo_r0[i, ]
  r1_values <- reo_r1[i, ]
  
  # ===== Criterion 1: R0 Directionality =====
  # Only consider |r| >= dead_zone
  r0_valid <- abs(r0_values) >= CONFIG$dead_zone
  r0_valid_values <- r0_values[r0_valid]
  n_r0_valid <- length(r0_valid_values)
  
  if (n_r0_valid >= (length(r0_values) - 1)) {
    # Enough valid samples (at least n-1)
    r0_signs <- sign(r0_valid_values)
    r0_positive <- sum(r0_signs > 0)
    r0_negative <- sum(r0_signs < 0)
    
    # Determine majority direction
    if (r0_positive >= r0_negative) {
      r0_direction <- 1
      r0_exceptions <- r0_negative
    } else {
      r0_direction <- -1
      r0_exceptions <- r0_positive
    }
    
    # Check if consistent (at most 1 exception)
    criterion1 <- r0_exceptions <= 1
    
    judgment_results$r0_direction[i] <- r0_direction
    judgment_results$r0_exceptions[i] <- r0_exceptions
    judgment_results$r0_valid_samples[i] <- n_r0_valid
    judgment_results$criterion1_pass[i] <- criterion1
  } else {
    # Not enough valid samples
    judgment_results$r0_valid_samples[i] <- n_r0_valid
    judgment_results$criterion1_pass[i] <- FALSE
    # Skip to next pair
    setTxtProgressBar(pb, i)
    next
  }
  
  # ===== Criterion 2: Strength (a[2] >= tau_strength) =====
  # Remove |r| < dead_zone
  r0_strong <- abs(r0_values) >= CONFIG$dead_zone
  r0_strong_values <- abs(r0_values[r0_strong])  # Use absolute values for strength
  n_r0_strong <- length(r0_strong_values)
  
  if (n_r0_strong >= (length(r0_values) - 1)) {
    # Sort and get a[2] (second smallest)
    r0_sorted <- sort(r0_strong_values)
    if (length(r0_sorted) >= 2) {
      a2_value <- r0_sorted[2]  # Second smallest
      criterion2 <- a2_value >= CONFIG$tau_strength
      
      judgment_results$r0_a2_value[i] <- a2_value
      judgment_results$r0_strong_samples[i] <- n_r0_strong
      judgment_results$criterion2_pass[i] <- criterion2
    } else {
      judgment_results$r0_strong_samples[i] <- n_r0_strong
      judgment_results$criterion2_pass[i] <- FALSE
    }
  } else {
    judgment_results$r0_strong_samples[i] <- n_r0_strong
    judgment_results$criterion2_pass[i] <- FALSE
  }
  
  # ===== Criterion 3: R1 Reversal =====
  # Remove |r| < dead_zone for R1 too
  r1_valid <- abs(r1_values) >= CONFIG$dead_zone
  r1_valid_values <- r1_values[r1_valid]
  n_r1_valid <- length(r1_valid_values)
  
  if (n_r1_valid > 0 && !is.na(r0_direction)) {
    # Count reversals (opposite sign to R0 expectation)
    if (r0_direction == 1) {
      # R0 expects positive, count negatives in R1
      r1_reversals <- sum(r1_valid_values < 0)
    } else {
      # R0 expects negative, count positives in R1
      r1_reversals <- sum(r1_valid_values > 0)
    }
    
    r1_reversal_rate <- r1_reversals / n_r1_valid
    
    # Wilson confidence interval (no continuity correction)
    if (n_r1_valid > 0) {
      # Wilson score interval
      z <- qnorm(0.975)  # 95% CI
      p_hat <- r1_reversal_rate
      n <- n_r1_valid
      
      denominator <- 1 + z^2/n
      center <- (p_hat + z^2/(2*n)) / denominator
      margin <- z * sqrt(p_hat * (1-p_hat) / n + z^2/(4*n^2)) / denominator
      
      ci_lower <- center - margin
      ci_upper <- center + margin
      
      # Ensure bounds are [0,1]
      ci_lower <- max(0, ci_lower)
      ci_upper <- min(1, ci_upper)
      
      criterion3 <- ci_lower > 0.5
      
      judgment_results$r1_reversal_rate[i] <- r1_reversal_rate
      judgment_results$r1_reversal_count[i] <- r1_reversals
      judgment_results$r1_valid_samples[i] <- n_r1_valid
      judgment_results$r1_ci_lower[i] <- ci_lower
      judgment_results$criterion3_pass[i] <- criterion3
    }
  } else {
    judgment_results$r1_valid_samples[i] <- n_r1_valid
    judgment_results$criterion3_pass[i] <- FALSE
  }
  
  # ===== Criterion 4: 100% Reversal Exclusion =====
  if (!is.na(judgment_results$r1_reversal_rate[i])) {
    if (judgment_results$r1_reversal_rate[i] == 1.0) {
      judgment_results$criterion4_pass[i] <- FALSE
    }
  }
  
  # ===== Overall Pass/Fail =====
  judgment_results$all_criteria_pass[i] <- 
    judgment_results$criterion1_pass[i] &
    judgment_results$criterion2_pass[i] &
    judgment_results$criterion3_pass[i] &
    judgment_results$criterion4_pass[i]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# --------------------------------------------------------------------------
# 3a.3 Summary of results
# --------------------------------------------------------------------------

cat("\n--- 3a.3 Summary of judgment results ---\n")

# Count passes for each criterion
cat("\nCriterion pass rates:\n")
cat(sprintf("  1. R0 directionality: %d / %d (%.1f%%)\n",
            sum(judgment_results$criterion1_pass),
            nrow(judgment_results),
            sum(judgment_results$criterion1_pass) / nrow(judgment_results) * 100))

cat(sprintf("  2. Strength (a[2]â‰¥%.2f): %d / %d (%.1f%%)\n",
            CONFIG$tau_strength,
            sum(judgment_results$criterion2_pass),
            nrow(judgment_results),
            sum(judgment_results$criterion2_pass) / nrow(judgment_results) * 100))

cat(sprintf("  3. R1 reversal (CI>0.5): %d / %d (%.1f%%)\n",
            sum(judgment_results$criterion3_pass),
            nrow(judgment_results),
            sum(judgment_results$criterion3_pass) / nrow(judgment_results) * 100))

cat(sprintf("  4. Not 100%% reversal: %d / %d (%.1f%%)\n",
            sum(judgment_results$criterion4_pass),
            nrow(judgment_results),
            sum(judgment_results$criterion4_pass) / nrow(judgment_results) * 100))

cat(sprintf("\nOverall pass (all criteria): %d / %d (%.1f%%)\n",
            sum(judgment_results$all_criteria_pass),
            nrow(judgment_results),
            sum(judgment_results$all_criteria_pass) / nrow(judgment_results) * 100))

# Filter to passed pairs
passed_pairs <- judgment_results[judgment_results$all_criteria_pass, ]
cat(sprintf("\nâœ“ Pairs passing Step 3a: %d\n", nrow(passed_pairs)))

# --------------------------------------------------------------------------
# 3a.4 Save Step 3a results
# --------------------------------------------------------------------------

cat("\n--- 3a.4 Saving Step 3a results ---\n")

# Combine with original pair info
step3a_results <- cbind(
  pairs_df,
  judgment_results[, -1]  # Remove duplicate pair_id
)

step3a_data <- list(
  # Configuration
  config = CONFIG,
  
  # Full results
  all_results = step3a_results,
  passed_pairs = step3a_results[step3a_results$all_criteria_pass, ],
  
  # REO matrices for passed pairs
  reo_r0_passed = reo_r0[judgment_results$all_criteria_pass, ],
  reo_r1_passed = reo_r1[judgment_results$all_criteria_pass, ],
  
  # Summary statistics
  n_input = nrow(pairs_df),
  n_passed = sum(judgment_results$all_criteria_pass),
  pass_rate = sum(judgment_results$all_criteria_pass) / nrow(pairs_df),
  
  # Sample info
  r0_samples = step2_data$r0_samples,
  r1_samples = step2_data$r1_samples,
  
  # Metadata
  timestamp = Sys.time(),
  step = "Step 3a Complete"
)

saveRDS(step3a_data, paste0(paths$processed, "reo_step3a_data.rds"))
cat("  Step 3a data saved to: reo_step3a_data.rds\n")

# Export summary CSV
write.csv(step3a_results[step3a_results$all_criteria_pass, ],
          paste0(paths$output, "step3a_passed_pairs.csv"),
          row.names = FALSE)
cat("  Passed pairs exported to: step3a_passed_pairs.csv\n")

cat("\n=== STEP 3a COMPLETE ===\n")

if (nrow(passed_pairs) > 0) {
  cat("Ready for Step 3b: Proxy Check (Fusion subtype)\n")
} else {
  cat("WARNING: No pairs passed Step 3a criteria!\n")
  cat("Consider adjusting parameters:\n")
  cat("  - tau_strength (currently", CONFIG$tau_strength, ")\n")
  cat("  - dead_zone (currently", CONFIG$dead_zone, ")\n")
}