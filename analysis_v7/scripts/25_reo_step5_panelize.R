# 25_reo_step5_panelize.R - REO Step 5: Panel Creation and QC
# Purpose: Create REO panel with quality control anchors
# Input: reo_step4_data.rds, thyr_deg_results.rds, thyr_tpm.rds
# Output: reo_step5_data.rds
# Version: v1.0
# Date: 2025-01-26
# Note: Implements N_eff-based quality control and QC anchors

source("analysis_v7/setup.R")

cat("\n=== REO STEP 5: Panel Creation and Quality Control ===\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(matrixStats)
})

# --------------------------------------------------------------------------
# 5.1 Configuration for Panel
# --------------------------------------------------------------------------

cat("--- 5.1 Panel Configuration ---\n")

# Panel-specific parameters (extend CONFIG)
PANEL_CONFIG <- list(
  # Panel size parameters
  panel_size_target = 10,        # Target panel size
  panel_size_min = 8,            # Minimum acceptable size
  panel_size_max = 12,           # Maximum panel size
  
  # Quality thresholds
  max_missing_rate = 0.20,       # Maximum missing rate for panel pairs
  n_eff_threshold = 0.7,         # Minimum effective rate for sample judgment
  judgment_threshold_ratio = 0.4, # Ratio for T calculation
  
  # QC anchor parameters
  anchor_delta_max = log2(1.2),  # Maximum group difference for anchors (same value as dead_zone but different meaning)
  anchor_count_target = 300,     # Number of QC anchor genes
  anchor_percentile_low = 0.2,   # Lower bound for expression (20th percentile)
  anchor_percentile_high = 0.8   # Upper bound for expression (80th percentile)
)

# Display configuration
for (param in names(PANEL_CONFIG)) {
  cat(sprintf("  %s: %s\n", param, PANEL_CONFIG[[param]]))
}

# --------------------------------------------------------------------------
# 5.2 Load Step 4 data
# --------------------------------------------------------------------------

cat("\n--- 5.2 Loading Step 4 data ---\n")

step4_data_path <- paste0(paths$processed, "reo_step4_data.rds")
if (!file.exists(step4_data_path)) {
  stop("ERROR: reo_step4_data.rds not found. Please run 24_reo_step4_dedup.R first.")
}

step4_data <- readRDS(step4_data_path)
cat("  Step 4 data loaded successfully\n")

# Extract data
final_pairs <- step4_data$final_pairs
final_reo_r0 <- step4_data$final_reo_r0
final_reo_r1 <- step4_data$final_reo_r1
CONFIG <- step4_data$config  # Get original CONFIG
r0_samples <- step4_data$r0_samples
r1_samples <- step4_data$r1_samples

cat(sprintf("  Input pairs from Step 4: %d\n", nrow(final_pairs)))

# --------------------------------------------------------------------------
# 5.3 Select panel pairs with quality filtering
# --------------------------------------------------------------------------

cat("\n--- 5.3 Selecting panel pairs ---\n")

# Filter by missing rate and keep track of original indices
quality_idx <- which(final_pairs$missing_rate <= PANEL_CONFIG$max_missing_rate)
quality_pairs <- final_pairs[quality_idx, ]
cat(sprintf("  After missing rate filter (≤%.1f%%): %d pairs\n", 
            PANEL_CONFIG$max_missing_rate * 100, nrow(quality_pairs)))

# Check minimum panel size
if (nrow(quality_pairs) < PANEL_CONFIG$panel_size_min) {
  stop(sprintf("ERROR: Only %d pairs available after quality filter. Minimum required: %d\n",
               nrow(quality_pairs), PANEL_CONFIG$panel_size_min))
}

# Determine actual panel size
N <- min(nrow(quality_pairs), PANEL_CONFIG$panel_size_target)
if (N > PANEL_CONFIG$panel_size_max) {
  N <- PANEL_CONFIG$panel_size_max
}

cat(sprintf("  Panel size N: %d\n", N))

# Select top N pairs using CORRECT indices from original data
# quality_idx contains the row numbers from final_pairs that passed the filter
panel_idx <- quality_idx[1:N]  # Get the first N indices from quality_pairs
panel_pairs <- final_pairs[panel_idx, ]
panel_reo_r0 <- final_reo_r0[panel_idx, ]
panel_reo_r1 <- final_reo_r1[panel_idx, ]

# Calculate judgment threshold
T <- ceiling(N * PANEL_CONFIG$judgment_threshold_ratio)
cat(sprintf("  Judgment threshold T: %d (%.0f%% of %d)\n", 
            T, PANEL_CONFIG$judgment_threshold_ratio * 100, N))

# Display panel summary
cat("\n  Panel composition:\n")
for (i in 1:N) {
  cat(sprintf("    Pair %2d: %s vs %s (a[2]=%.2f, rev_rate=%.2f, miss_rate=%.1f%%)\n",
              i,
              panel_pairs$gene_up[i],
              panel_pairs$gene_down[i],
              panel_pairs$a2_value[i],
              panel_pairs$r1_reversal_rate[i],
              panel_pairs$missing_rate[i] * 100))
}

# --------------------------------------------------------------------------
# 5.4 Calculate REO scores for all samples
# --------------------------------------------------------------------------

cat("\n--- 5.4 Calculating REO scores ---\n")

# Function to calculate REO score for one sample
calculate_reo_score <- function(reo_values, r0_directions_vec, dead_zone, N) {
  # Count reversals
  # A reversal is when the sign differs from R0 direction
  # Only count if |r| >= dead_zone
  
  valid_mask <- abs(reo_values) >= dead_zone
  n_eff <- sum(valid_mask, na.rm = TRUE)
  
  if (n_eff == 0) {
    return(list(
      reo_score = NA,
      n_reversals = NA,
      n_eff = 0,
      judgment_possible = FALSE
    ))
  }
  
  # Get valid values and their corresponding R0 directions
  valid_values <- reo_values[valid_mask]
  valid_directions <- sign(valid_values)
  valid_r0_directions <- r0_directions_vec[valid_mask]
  
  # Count reversals (opposite to R0 direction)
  n_reversals <- sum(valid_directions != valid_r0_directions, na.rm = TRUE)
  
  # REO score = reversals / N (denominator fixed)
  reo_score <- n_reversals / N
  
  # Check if judgment is possible
  judgment_possible <- n_eff >= ceiling(PANEL_CONFIG$n_eff_threshold * N)
  
  return(list(
    reo_score = reo_score,
    n_reversals = n_reversals,
    n_eff = n_eff,
    judgment_possible = judgment_possible
  ))
}

# Get R0 directions for each pair (from Step 3a results)
r0_directions <- panel_pairs$r0_direction

# Calculate scores for R0 samples
r0_scores <- data.frame(
  sample_id = r0_samples,
  group = "R0",
  reo_score = NA_real_,
  n_reversals = NA_integer_,
  n_eff = NA_integer_,
  judgment_possible = FALSE,
  stringsAsFactors = FALSE
)

for (i in seq_along(r0_samples)) {
  sample_reo <- panel_reo_r0[, i]
  result <- calculate_reo_score(sample_reo, r0_directions, CONFIG$dead_zone, N)
  r0_scores[i, c("reo_score", "n_reversals", "n_eff", "judgment_possible")] <- result
}

# Calculate scores for R1 samples
r1_scores <- data.frame(
  sample_id = r1_samples,
  group = "R1",
  reo_score = NA_real_,
  n_reversals = NA_integer_,
  n_eff = NA_integer_,
  judgment_possible = FALSE,
  stringsAsFactors = FALSE
)

for (i in seq_along(r1_samples)) {
  sample_reo <- panel_reo_r1[, i]
  result <- calculate_reo_score(sample_reo, r0_directions, CONFIG$dead_zone, N)
  r1_scores[i, c("reo_score", "n_reversals", "n_eff", "judgment_possible")] <- result
}

# Combine scores
all_scores <- rbind(r0_scores, r1_scores)

# Summary statistics
cat("\n  REO Score Summary:\n")
cat("  R0 group:\n")
r0_valid <- r0_scores[r0_scores$judgment_possible, ]
cat(sprintf("    Samples with valid judgment: %d/%d (%.1f%%)\n",
            nrow(r0_valid), nrow(r0_scores),
            nrow(r0_valid) / nrow(r0_scores) * 100))
if (nrow(r0_valid) > 0) {
  cat(sprintf("    Mean REO score: %.3f\n", mean(r0_valid$reo_score)))
  cat(sprintf("    Median REO score: %.3f\n", median(r0_valid$reo_score)))
  cat(sprintf("    Range: %.3f - %.3f\n", 
              min(r0_valid$reo_score), max(r0_valid$reo_score)))
}

cat("\n  R1 group:\n")
r1_valid <- r1_scores[r1_scores$judgment_possible, ]
cat(sprintf("    Samples with valid judgment: %d/%d (%.1f%%)\n",
            nrow(r1_valid), nrow(r1_scores),
            nrow(r1_valid) / nrow(r1_scores) * 100))
if (nrow(r1_valid) > 0) {
  cat(sprintf("    Mean REO score: %.3f\n", mean(r1_valid$reo_score)))
  cat(sprintf("    Median REO score: %.3f\n", median(r1_valid$reo_score)))
  cat(sprintf("    Range: %.3f - %.3f\n", 
              min(r1_valid$reo_score), max(r1_valid$reo_score)))
}

# Check if we have enough valid samples
min_valid_samples <- 5  # Minimum for meaningful statistics
if (nrow(r0_valid) < min_valid_samples || nrow(r1_valid) < min_valid_samples) {
  warning("WARNING: Very few samples with valid judgment. Consider relaxing n_eff_threshold.")
}

# --------------------------------------------------------------------------
# 5.5 Select QC anchor genes
# --------------------------------------------------------------------------

cat("\n--- 5.5 Selecting QC anchor genes ---\n")

# Load TPM data for anchor selection
tpm_path <- paste0(paths$processed, "thyr_tpm.rds")
if (!file.exists(tpm_path)) {
  cat("  TPM data not found. Attempting to create from counts...\n")
  
  # Load counts and create TPM (simplified version)
  se_path <- paste0(paths$processed, "thyr_se_strand2_nonzero.rds")
  if (file.exists(se_path)) {
    se <- readRDS(se_path)
    counts <- assay(se, "counts")  # Fixed: use "counts" instead of "stranded_second"
    
    # Simple TPM calculation (without effective length)
    rpk <- counts / 1000  # Assume 1kb length for all
    scaling_factor <- colSums(rpk) / 1e6
    tpm <- sweep(rpk, 2, scaling_factor, "/")
    
    saveRDS(tpm, tpm_path)
    cat("  TPM data created and saved.\n")
  } else {
    cat("  WARNING: Cannot find count data. Skipping QC anchor selection.\n")
    tpm <- NULL
  }
} else {
  tpm <- readRDS(tpm_path)
  cat("  TPM data loaded successfully.\n")
}

if (!is.null(tpm)) {
  # Get samples used in analysis
  all_analysis_samples <- c(r0_samples, r1_samples)
  tpm_subset <- tpm[, all_analysis_samples, drop = FALSE]
  
  # Calculate log2 TPM
  log2_tpm <- log2(tpm_subset + 1)
  
  # Load DEG results to identify non-DEGs
  deg_path <- paste0(paths$processed, "thyr_deg_results.rds")
  if (file.exists(deg_path)) {
    deg_results <- readRDS(deg_path)
    
    # Get DEGs from R0 vs R1 comparison
    if ("R0_vs_R1_tumor" %in% names(deg_results)) {
      degs <- deg_results$R0_vs_R1_tumor$gene_id[deg_results$R0_vs_R1_tumor$qvalue < 0.05]
    } else {
      degs <- character(0)
    }
    
    # Identify non-DEGs
    all_genes <- rownames(log2_tpm)
    non_degs <- setdiff(all_genes, degs)
    cat(sprintf("  Non-DEG genes available: %d\n", length(non_degs)))
    
    # Calculate median expression across all samples
    gene_medians <- rowMedians(as.matrix(log2_tpm))
    names(gene_medians) <- rownames(log2_tpm)
    
    # Get percentile thresholds
    all_medians <- gene_medians[non_degs]
    percentile_low <- quantile(all_medians, PANEL_CONFIG$anchor_percentile_low)
    percentile_high <- quantile(all_medians, PANEL_CONFIG$anchor_percentile_high)
    
    cat(sprintf("  Expression thresholds at %dth and %dth percentiles: %.2f, %.2f (log2 TPM)\n",
                PANEL_CONFIG$anchor_percentile_low * 100,
                PANEL_CONFIG$anchor_percentile_high * 100,
                percentile_low, percentile_high))
    
    # Filter candidates
    candidate_anchors <- non_degs[
      gene_medians[non_degs] >= percentile_low &
        gene_medians[non_degs] <= percentile_high
    ]
    cat(sprintf("  Candidates within expression thresholds: %d\n", length(candidate_anchors)))
    
    # Calculate group differences
    r0_medians <- rowMedians(as.matrix(log2_tpm[, r0_samples, drop = FALSE]))
    r1_medians <- rowMedians(as.matrix(log2_tpm[, r1_samples, drop = FALSE]))
    names(r0_medians) <- rownames(log2_tpm)
    names(r1_medians) <- rownames(log2_tpm)
    
    group_diffs <- abs(r0_medians - r1_medians)
    
    # Filter by group difference
    stable_anchors <- candidate_anchors[
      group_diffs[candidate_anchors] <= PANEL_CONFIG$anchor_delta_max
    ]
    cat(sprintf("  After group difference filter (≤%.3f): %d\n",
                PANEL_CONFIG$anchor_delta_max, length(stable_anchors)))
    
    if (length(stable_anchors) >= PANEL_CONFIG$anchor_count_target) {
      # Calculate MAD for each candidate
      mad_values <- apply(log2_tpm[stable_anchors, , drop = FALSE], 1, mad)
      
      # Select top K by MAD
      qc_anchors <- names(sort(mad_values))[1:min(PANEL_CONFIG$anchor_count_target, 
                                                  length(stable_anchors))]
      
      cat(sprintf("  QC anchors selected: %d genes\n", length(qc_anchors)))
      
      # Create anchor dataframe
      anchor_df <- data.frame(
        gene = qc_anchors,
        median_expression = gene_medians[qc_anchors],
        group_diff = group_diffs[qc_anchors],
        mad = mad_values[qc_anchors],
        stringsAsFactors = FALSE
      ) %>%
        arrange(mad)
      
    } else {
      cat(sprintf("  WARNING: Only %d stable anchors available (target: %d)\n",
                  length(stable_anchors), PANEL_CONFIG$anchor_count_target))
      qc_anchors <- stable_anchors
      
      # Create anchor dataframe with available anchors
      mad_values <- apply(log2_tpm[stable_anchors, , drop = FALSE], 1, mad)
      anchor_df <- data.frame(
        gene = stable_anchors,
        median_expression = gene_medians[stable_anchors],
        group_diff = group_diffs[stable_anchors],
        mad = mad_values,
        stringsAsFactors = FALSE
      ) %>%
        arrange(mad)
    }
  } else {
    cat("  WARNING: Cannot identify non-DEGs for anchors.\n")
    qc_anchors <- character(0)
    anchor_df <- data.frame()
  }
} else {
  qc_anchors <- character(0)
  anchor_df <- data.frame()
}

# --------------------------------------------------------------------------
# 5.6 Create panel classification function
# --------------------------------------------------------------------------

cat("\n--- 5.6 Creating classification function ---\n")

# Classification function for new samples
classify_sample <- function(sample_reo, panel_info, config) {
  N <- panel_info$N
  T <- panel_info$T
  r0_directions <- panel_info$r0_directions
  dead_zone <- config$dead_zone
  n_eff_threshold <- config$n_eff_threshold
  
  # Check input dimensions
  if (length(sample_reo) != N) {
    stop(sprintf("Input must have %d REO values", N))
  }
  
  # Calculate N_eff and reversals
  valid_mask <- abs(sample_reo) >= dead_zone
  n_eff <- sum(valid_mask, na.rm = TRUE)
  
  # Check if judgment is possible
  if (n_eff < ceiling(n_eff_threshold * N)) {
    return(list(
      classification = "Undetermined",
      reo_score = NA,
      n_reversals = NA,
      n_eff = n_eff,
      reason = "Insufficient valid pairs"
    ))
  }
  
  # Count reversals
  valid_values <- sample_reo[valid_mask]
  valid_directions <- sign(valid_values)
  expected_directions <- r0_directions[valid_mask]
  n_reversals <- sum(valid_directions != expected_directions, na.rm = TRUE)
  
  # Calculate REO score
  reo_score <- n_reversals / N
  
  # Classify based on threshold
  if (n_reversals >= T) {
    classification <- "R1-like"
  } else {
    classification <- "R0-like"
  }
  
  return(list(
    classification = classification,
    reo_score = reo_score,
    n_reversals = n_reversals,
    n_eff = n_eff,
    reason = "Valid classification"
  ))
}

cat("  Classification function created.\n")
cat(sprintf("  Decision rule: ≥%d reversals → R1-like, <%d → R0-like\n", T, T))

# --------------------------------------------------------------------------
# 5.7 Save Step 5 results
# --------------------------------------------------------------------------

cat("\n--- 5.7 Saving Step 5 results ---\n")

step5_data <- list(
  # Configuration
  config = CONFIG,
  panel_config = PANEL_CONFIG,
  
  # Panel information
  panel_pairs = panel_pairs,
  panel_reo_r0 = panel_reo_r0,
  panel_reo_r1 = panel_reo_r1,
  N = N,
  T = T,
  
  # REO scores
  sample_scores = all_scores,
  
  # QC anchors
  qc_anchors = qc_anchors,
  anchor_df = anchor_df,
  
  # Classification function and info
  panel_info = list(
    N = N,
    T = T,
    r0_directions = r0_directions,
    gene_up = panel_pairs$gene_up,
    gene_down = panel_pairs$gene_down
  ),
  
  # Summary statistics
  summary_stats = list(
    n_input_pairs = nrow(final_pairs),
    n_quality_pairs = nrow(quality_pairs),
    n_panel_pairs = N,
    threshold_T = T,
    n_r0_valid = nrow(r0_valid),
    n_r1_valid = nrow(r1_valid),
    r0_mean_score = ifelse(nrow(r0_valid) > 0, mean(r0_valid$reo_score), NA),
    r1_mean_score = ifelse(nrow(r1_valid) > 0, mean(r1_valid$reo_score), NA),
    n_qc_anchors = length(qc_anchors)
  ),
  
  # Sample info
  r0_samples = r0_samples,
  r1_samples = r1_samples,
  
  # Metadata
  timestamp = Sys.time(),
  step = "Step 5 Complete"
)

saveRDS(step5_data, paste0(paths$processed, "reo_step5_data.rds"))
cat("  Step 5 data saved to: reo_step5_data.rds\n")

# Export panel definition
panel_def <- data.frame(
  pair_idx = 1:N,
  gene_up = panel_pairs$gene_up,
  gene_down = panel_pairs$gene_down,
  r0_direction = r0_directions,
  a2_value = panel_pairs$a2_value,
  r1_reversal_rate = panel_pairs$r1_reversal_rate,
  missing_rate = panel_pairs$missing_rate
)
write.csv(panel_def, 
          paste0(paths$output, "step5_reo_panel_definition.csv"),
          row.names = FALSE)
cat("  Panel definition exported to: step5_reo_panel_definition.csv\n")

# Export sample scores
write.csv(all_scores,
          paste0(paths$output, "step5_sample_scores.csv"),
          row.names = FALSE)
cat("  Sample scores exported to: step5_sample_scores.csv\n")

# Export QC anchors if available
if (length(qc_anchors) > 0) {
  write.csv(anchor_df,
            paste0(paths$output, "step5_qc_anchors.csv"),
            row.names = FALSE)
  cat("  QC anchors exported to: step5_qc_anchors.csv\n")
}

# =============================================================================
# STEP 5 SUMMARY (固定フォーマット)
# =============================================================================

cat("\n=== STEP 5 SUMMARY ===\n")
cat(sprintf("Input pairs: %d\n", nrow(final_pairs)))
cat(sprintf("After quality filter: %d\n", nrow(quality_pairs)))
cat(sprintf("Panel size N: %d\n", N))
cat(sprintf("Judgment threshold T: %d (≥%.0f%% reversals)\n", 
            T, PANEL_CONFIG$judgment_threshold_ratio * 100))
cat(sprintf("R0 valid judgments: %d/%d\n", nrow(r0_valid), length(r0_samples)))
cat(sprintf("R1 valid judgments: %d/%d\n", nrow(r1_valid), length(r1_samples)))
cat(sprintf("QC anchors selected: %d\n", length(qc_anchors)))
cat(sprintf("Timestamp: %s\n", format(step5_data$timestamp)))

# Classification performance preview
if (nrow(r0_valid) > 0 && nrow(r1_valid) > 0) {
  r0_classified_r1 <- sum(r0_valid$n_reversals >= T)
  r1_classified_r1 <- sum(r1_valid$n_reversals >= T)
  
  cat("\nClassification preview:\n")
  cat(sprintf("  R0 samples classified as R0-like: %d/%d (%.1f%%)\n",
              nrow(r0_valid) - r0_classified_r1, nrow(r0_valid),
              (nrow(r0_valid) - r0_classified_r1) / nrow(r0_valid) * 100))
  cat(sprintf("  R1 samples classified as R1-like: %d/%d (%.1f%%)\n",
              r1_classified_r1, nrow(r1_valid),
              r1_classified_r1 / nrow(r1_valid) * 100))
}

# Final message
if (N >= PANEL_CONFIG$panel_size_min) {
  cat(sprintf("\n✓ Panel created successfully. Ready for Step 6: R_test validation\n"))
} else {
  cat(sprintf("\n⚠ Panel size %d is below minimum. Review data quality.\n", N))
}