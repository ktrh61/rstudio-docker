# 21_reo_step2_pairs.R - REO Step 2: Candidate Pair Generation
# Purpose: Generate candidate gene pairs with expression filtering and multiplicity limit
# Input: reo_step1_data.rds
# Output: reo_step2_data.rds
# Version: v1.0
# Date: 2025-01-26

source("analysis_v7/setup.R")

cat("\n=== REO STEP 2: Candidate Pair Generation ===\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
})

# --------------------------------------------------------------------------
# 2.1 Load Step 1 data
# --------------------------------------------------------------------------

cat("--- 2.1 Loading Step 1 data ---\n")

step1_data_path <- paste0(paths$processed, "reo_step1_data.rds")
if (!file.exists(step1_data_path)) {
  stop("ERROR: reo_step1_data.rds not found. Please run 20_reo_step1_prep.R first.")
}

step1_data <- readRDS(step1_data_path)
cat("  Step 1 data loaded successfully\n")

# 薄い契約：最小限の存在チェック
stopifnot(all(c("log2_tpm_r0", "log2_tpm_r1", "sig_degs", "config") %in% names(step1_data)))

# CONFIGはRDSから取得（reo_config.Rは読まない）
CONFIG <- step1_data$config

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
cat(sprintf("  Combined matrix: %d genes × %d samples\n", 
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
trim_low <- CONFIG$trim_low    # Now 0.07 (adjusted)
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
# 2.4 Apply one-sided multiplicity limit (M=12)
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
up_limited_pairs <- up_result$pairs
up_hitting_limit <- up_result$hitting_limit

# Apply M limit for down genes (each down gene pairs with at most M up genes)  
down_result <- select_top_partners(down_filtered, up_filtered, M = CONFIG$M)
down_limited_pairs <- down_result$pairs
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
      cat(sprintf("    - trim_low: 0.10 → %.2f\n", CONFIG$trim_low))
    }
    if (CONFIG$M != 10) {
      cat(sprintf("    - M: 10 → %d\n", CONFIG$M))
    }
  }
  
  cat("\n  Consider further adjustments:\n")
  cat("    - Further reduce trim thresholds\n")
  cat("    - Further increase M limit\n")
  cat("    - Check DEG quality\n")
  stop(sprintf("STOPPING: Only %d pairs, need at least %d", 
               final_candidates, min_candidates_required))
} else {
  cat(sprintf("\n  ✓ Sufficient candidates: %d ≥ %d\n", 
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

# Convert list to data frame
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
  # Configuration (pass through)
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

# =============================================================================
# STEP 2 SUMMARY (固定フォーマット)
# =============================================================================

cat("\n=== STEP 2 SUMMARY ===\n")
cat(sprintf("Initial pairs: %s\n", format(initial_pairs_raw, big.mark = ",")))
cat(sprintf("After expression trim: %s\n", format(initial_pairs_after_trim, big.mark = ",")))
cat(sprintf("After M=%d limit: %d\n", CONFIG$M, final_candidates))
cat(sprintf("Retention rate: %.1f%%\n", final_candidates / initial_pairs_raw * 100))
cat(sprintf("Timestamp: %s\n", format(step2_data$timestamp)))