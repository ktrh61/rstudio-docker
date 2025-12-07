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
# Initialize from CONFIG if available, otherwise use defaults
default_panel_config <- list(
  # Panel size parameters (matching reo_config.R naming)
  size_target = 10,              # Target panel size
  size_min = 8,                  # Minimum acceptable size
  size_max = 12,                 # Maximum panel size
  t_ratio = 0.4,                 # Threshold ratio for T calculation
  max_missing_rate = 0.20,       # Maximum missing rate for panel pairs
  
  # Additional parameters not in CONFIG$panel_defaults
  n_eff_threshold = 0.7,         # Minimum effective rate for sample judgment
  
  # QC anchor parameters (from CONFIG$qc_anchor if available)
  anchor_delta_max = log2(1.2),  # Maximum group difference for anchors
  anchor_count_target = 300,     # Number of QC anchor genes
  anchor_percentile_low = 0.2,   # Lower bound for expression (20th percentile)
  anchor_percentile_high = 0.8   # Upper bound for expression (80th percentile)
)

if (exists("panel_defaults", where = CONFIG)) {
  cat("  Initializing from CONFIG$panel_defaults\n")
  PANEL_CONFIG <- CONFIG$panel_defaults
  
  # Add n_eff_threshold from main CONFIG if available
  if (exists("n_eff_threshold", where = CONFIG)) {
    PANEL_CONFIG$n_eff_threshold <- CONFIG$n_eff_threshold
  } else {
    PANEL_CONFIG$n_eff_threshold <- default_panel_config$n_eff_threshold
  }
  
  # Add QC anchor parameters from CONFIG$qc_anchor if available
  if (exists("qc_anchor", where = CONFIG)) {
    PANEL_CONFIG$anchor_delta_max <- CONFIG$qc_anchor$delta_max
    PANEL_CONFIG$anchor_count_target <- CONFIG$qc_anchor$n
    PANEL_CONFIG$anchor_percentile_low <- CONFIG$qc_anchor$expr_pct_low
    PANEL_CONFIG$anchor_percentile_high <- CONFIG$qc_anchor$expr_pct_high
  } else {
    # Use defaults for QC anchor parameters
    PANEL_CONFIG$anchor_delta_max <- default_panel_config$anchor_delta_max
    PANEL_CONFIG$anchor_count_target <- default_panel_config$anchor_count_target
    PANEL_CONFIG$anchor_percentile_low <- default_panel_config$anchor_percentile_low
    PANEL_CONFIG$anchor_percentile_high <- default_panel_config$anchor_percentile_high
  }
  
  # Ensure all required parameters exist
  for (param in names(default_panel_config)) {
    if (!exists(param, where = PANEL_CONFIG)) {
      PANEL_CONFIG[[param]] <- default_panel_config[[param]]
      cat(sprintf("    Added missing parameter: %s = %s\n", param, default_panel_config[[param]]))
    }
  }
} else {
  cat("  Using default panel configuration\n")
  PANEL_CONFIG <- default_panel_config
}

# Validate critical parameters (using correct names)
stopifnot("size_min must exist" = !is.null(PANEL_CONFIG$size_min))
stopifnot("size_max must exist" = !is.null(PANEL_CONFIG$size_max))
stopifnot("size_target must exist" = !is.null(PANEL_CONFIG$size_target))
stopifnot("max_missing_rate must exist" = !is.null(PANEL_CONFIG$max_missing_rate))
stopifnot("t_ratio must exist" = !is.null(PANEL_CONFIG$t_ratio))

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

# Validate required components
required_components <- c("final_pairs", "final_reo_r0", "final_reo_r1", 
                         "config", "r0_samples", "r1_samples")
missing_components <- setdiff(required_components, names(step4_data))
if (length(missing_components) > 0) {
  stop(sprintf("ERROR: Missing required components in step4_data: %s",
               paste(missing_components, collapse = ", ")))
}

# Extract data with validation
final_pairs <- step4_data$final_pairs
if (!is.data.frame(final_pairs)) {
  stop("ERROR: final_pairs is not a data frame")
}

final_reo_r0 <- step4_data$final_reo_r0
final_reo_r1 <- step4_data$final_reo_r1
CONFIG <- step4_data$config  # Get original CONFIG
r0_samples <- step4_data$r0_samples
r1_samples <- step4_data$r1_samples

# Validate final_pairs structure
required_cols <- c("gene_up", "gene_down", "missing_rate", "a2_value", 
                   "r1_reversal_rate", "r0_direction")
missing_cols <- setdiff(required_cols, names(final_pairs))
if (length(missing_cols) > 0) {
  cat("  WARNING: Missing expected columns in final_pairs:", 
      paste(missing_cols, collapse = ", "), "\n")
}

cat(sprintf("  Input pairs from Step 4: %d\n", nrow(final_pairs)))
cat(sprintf("  Data dimensions - REO_R0: %d x %d, REO_R1: %d x %d\n",
            nrow(final_reo_r0), ncol(final_reo_r0),
            nrow(final_reo_r1), ncol(final_reo_r1)))

# --------------------------------------------------------------------------
# 5.3 Select panel pairs with quality filtering
# --------------------------------------------------------------------------

cat("\n--- 5.3 Selecting panel pairs ---\n")

# Ensure final_pairs is a data frame
if (!is.data.frame(final_pairs)) {
  cat("  WARNING: Converting final_pairs to data frame\n")
  final_pairs <- as.data.frame(final_pairs)
}

# Filter by missing rate and keep track of original indices
quality_idx <- which(final_pairs$missing_rate <= PANEL_CONFIG$max_missing_rate)
cat(sprintf("  Pairs passing missing rate filter (≤%.1f%%): %d\n", 
            PANEL_CONFIG$max_missing_rate * 100, length(quality_idx)))

# Check if we have any pairs
if (length(quality_idx) == 0) {
  stop("ERROR: No pairs passed the missing rate filter")
}

# Create quality_pairs subset
quality_pairs <- final_pairs[quality_idx, , drop = FALSE]

# Verify quality_pairs is valid
if (!is.data.frame(quality_pairs) || nrow(quality_pairs) == 0) {
  stop(sprintf("ERROR: Failed to create quality_pairs subset. quality_idx length: %d", 
               length(quality_idx)))
}

cat(sprintf("  Quality pairs created: %d rows\n", nrow(quality_pairs)))

# Check minimum panel size
if (nrow(quality_pairs) < PANEL_CONFIG$size_min) {
  stop(sprintf("ERROR: Only %d pairs available after quality filter. Minimum required: %d\n",
               nrow(quality_pairs), PANEL_CONFIG$size_min))
}

# Determine actual panel size
N <- min(nrow(quality_pairs), PANEL_CONFIG$size_target)
if (N > PANEL_CONFIG$size_max) {
  N <- PANEL_CONFIG$size_max
}

cat(sprintf("  Panel size N: %d\n", N))

# Select top N pairs 
# Note: quality_idx contains the row indices from final_pairs
# We need the first N indices that passed the quality filter
panel_idx <- quality_idx[1:N]  # Indices in original final_pairs
panel_pairs <- final_pairs[panel_idx, , drop = FALSE]

# Validate matrix subsetting
if (!is.matrix(final_reo_r0) && !is.data.frame(final_reo_r0)) {
  stop("ERROR: final_reo_r0 is neither matrix nor data.frame")
}
if (!is.matrix(final_reo_r1) && !is.data.frame(final_reo_r1)) {
  stop("ERROR: final_reo_r1 is neither matrix nor data.frame")
}

# Create panel REO matrices
if (is.matrix(final_reo_r0)) {
  panel_reo_r0 <- final_reo_r0[panel_idx, , drop = FALSE]
} else {
  panel_reo_r0 <- as.matrix(final_reo_r0[panel_idx, , drop = FALSE])
}

if (is.matrix(final_reo_r1)) {
  panel_reo_r1 <- final_reo_r1[panel_idx, , drop = FALSE]
} else {
  panel_reo_r1 <- as.matrix(final_reo_r1[panel_idx, , drop = FALSE])
}

# Validate dimensions
cat(sprintf("  Panel dimensions - Pairs: %d, R0 samples: %d, R1 samples: %d\n",
            nrow(panel_pairs), ncol(panel_reo_r0), ncol(panel_reo_r1)))

# Calculate judgment threshold
T <- ceiling(N * PANEL_CONFIG$t_ratio)
cat(sprintf("  Judgment threshold T: %d (%.0f%% of %d)\n", 
            T, PANEL_CONFIG$t_ratio * 100, N))

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
calculate_reo_score <- function(reo_values, r0_directions_vec, dead_zone, N, T) {
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
      judgment_possible = FALSE,
      boundary_flag = NA,
      decision_margin = NA,
      boundary_type = "undetermined"
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
  
  # Boundary zone analysis (only if judgment possible)
  if (judgment_possible) {
    # Decision margin: distance from threshold T
    decision_margin <- abs(n_reversals - T)
    
    # Boundary flag: TRUE if in boundary zone (T-1 or T reversals)
    boundary_flag <- (n_reversals == T - 1) || (n_reversals == T)
    
    # Boundary type classification
    if (n_reversals < T - 1) {
      boundary_type <- "clear_R0"
    } else if (n_reversals == T - 1 || n_reversals == T) {
      boundary_type <- "boundary"
    } else {
      boundary_type <- "clear_R1"
    }
  } else {
    decision_margin <- NA
    boundary_flag <- NA
    boundary_type <- "undetermined"
  }
  
  return(list(
    reo_score = reo_score,
    n_reversals = n_reversals,
    n_eff = n_eff,
    judgment_possible = judgment_possible,
    boundary_flag = boundary_flag,
    decision_margin = decision_margin,
    boundary_type = boundary_type
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
  boundary_flag = NA,
  decision_margin = NA_integer_,
  boundary_type = character(length(r0_samples)),
  stringsAsFactors = FALSE
)

for (i in seq_along(r0_samples)) {
  sample_reo <- panel_reo_r0[, i]
  result <- calculate_reo_score(sample_reo, r0_directions, CONFIG$dead_zone, N, T)
  r0_scores[i, c("reo_score", "n_reversals", "n_eff", "judgment_possible", 
                 "boundary_flag", "decision_margin", "boundary_type")] <- result
}

# Calculate scores for R1 samples
r1_scores <- data.frame(
  sample_id = r1_samples,
  group = "R1",
  reo_score = NA_real_,
  n_reversals = NA_integer_,
  n_eff = NA_integer_,
  judgment_possible = FALSE,
  boundary_flag = NA,
  decision_margin = NA_integer_,
  boundary_type = character(length(r1_samples)),
  stringsAsFactors = FALSE
)

for (i in seq_along(r1_samples)) {
  sample_reo <- panel_reo_r1[, i]
  result <- calculate_reo_score(sample_reo, r0_directions, CONFIG$dead_zone, N, T)
  r1_scores[i, c("reo_score", "n_reversals", "n_eff", "judgment_possible",
                 "boundary_flag", "decision_margin", "boundary_type")] <- result
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
  cat(sprintf("    Mean decision margin: %.1f reversals\n", mean(r0_valid$decision_margin)))
  
  # Boundary type distribution
  boundary_dist <- table(r0_valid$boundary_type)
  cat("    Distribution:\n")
  for (bt in names(boundary_dist)) {
    cat(sprintf("      %s: %d (%.1f%%)\n", bt, boundary_dist[bt], 
                boundary_dist[bt] / nrow(r0_valid) * 100))
  }
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
  cat(sprintf("    Mean decision margin: %.1f reversals\n", mean(r1_valid$decision_margin)))
  
  # Boundary type distribution
  boundary_dist <- table(r1_valid$boundary_type)
  cat("    Distribution:\n")
  for (bt in names(boundary_dist)) {
    cat(sprintf("      %s: %d (%.1f%%)\n", bt, boundary_dist[bt], 
                boundary_dist[bt] / nrow(r1_valid) * 100))
  }
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

# Check if log2 TPM data exists in Step 4 data
if (!all(c("log2_tpm_r0", "log2_tpm_r1") %in% names(step4_data))) {
  cat("  WARNING: TPM data not found in Step 4 data.\n")
  cat("  QC anchor selection requires TPM data from upstream steps.\n")
  cat("  Skipping QC anchor selection.\n")
  qc_anchors <- character(0)
  anchor_df <- data.frame()
} else if (is.null(step4_data$log2_tpm_r0) || is.null(step4_data$log2_tpm_r1)) {
  cat("  WARNING: TPM data is NULL in Step 4 data.\n")
  cat("  This indicates TPM was not available in Step 3b.\n")
  cat("  Skipping QC anchor selection.\n")
  qc_anchors <- character(0)
  anchor_df <- data.frame()
} else {
  # Get log2 TPM matrices from Step 4 (originally from Step 1->2->3->4)
  log2_tpm_r0 <- step4_data$log2_tpm_r0
  log2_tpm_r1 <- step4_data$log2_tpm_r1
  
  # Combine R0 and R1 for full analysis
  log2_tpm <- cbind(log2_tpm_r0, log2_tpm_r1)
  cat(sprintf("  TPM data loaded from Step 4: %d genes × %d samples\n", 
              nrow(log2_tpm), ncol(log2_tpm)))
  cat("  Note: This is zero-free log2 TPM from Step 1 processing\n")
  
  # TPM data from Step 4 already has correct samples (zero-free and no need for subsetting)
  # Verify data integrity
  if (any(!is.finite(as.matrix(log2_tpm)))) {
    cat("  WARNING: Non-finite values found in log2 TPM.\n")
    cat("  This should not happen with zero-free data from Step 1. Investigating...\n")
    n_inf <- sum(!is.finite(as.matrix(log2_tpm)))
    cat(sprintf("  Non-finite values: %d\n", n_inf))
  }
  
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
    
    # Identify non-DEGs (only from genes remaining after filtering)
    all_genes <- rownames(log2_tpm)
    non_degs <- setdiff(all_genes, degs)
    cat(sprintf("  Non-DEG genes available: %d\n", length(non_degs)))
    
    if (length(non_degs) == 0) {
      cat("  WARNING: No non-DEG genes found after filtering.\n")
      qc_anchors <- character(0)
      anchor_df <- data.frame()
    } else {
      # Calculate median expression across all samples
      gene_medians <- rowMedians(as.matrix(log2_tpm))
      names(gene_medians) <- rownames(log2_tpm)
      
      # Get percentile thresholds (only for non-DEGs)
      all_medians <- gene_medians[non_degs]
      
      # Check for valid values
      if (any(!is.finite(all_medians))) {
        cat("  WARNING: Non-finite median values found. Removing.\n")
        all_medians <- all_medians[is.finite(all_medians)]
      }
      
      if (length(all_medians) == 0) {
        cat("  ERROR: No valid median values for non-DEGs.\n")
        qc_anchors <- character(0)
        anchor_df <- data.frame()
      } else {
        percentile_low <- quantile(all_medians, PANEL_CONFIG$anchor_percentile_low, na.rm = TRUE)
        percentile_high <- quantile(all_medians, PANEL_CONFIG$anchor_percentile_high, na.rm = TRUE)
        
        cat(sprintf("  Expression thresholds at %dth and %dth percentiles: %.2f, %.2f (log2 TPM)\n",
                    PANEL_CONFIG$anchor_percentile_low * 100,
                    PANEL_CONFIG$anchor_percentile_high * 100,
                    percentile_low, percentile_high))
        
        # Filter candidates (ensure finite thresholds)
        if (is.finite(percentile_low) && is.finite(percentile_high)) {
          candidate_anchors <- non_degs[
            gene_medians[non_degs] >= percentile_low &
              gene_medians[non_degs] <= percentile_high &
              is.finite(gene_medians[non_degs])
          ]
        } else {
          cat("  WARNING: Non-finite percentile thresholds. Cannot filter candidates.\n")
          candidate_anchors <- character(0)
        }
        
        cat(sprintf("  Candidates within expression thresholds: %d\n", length(candidate_anchors)))
        
        if (length(candidate_anchors) == 0) {
          cat("  WARNING: No candidate anchors found.\n")
          qc_anchors <- character(0)
          anchor_df <- data.frame()
        } else {
          # Calculate group differences using R0 and R1 separately
          # (log2_tpm_r0 and log2_tpm_r1 from Step 4 already have correct samples)
          r0_medians <- rowMedians(as.matrix(log2_tpm_r0))
          r1_medians <- rowMedians(as.matrix(log2_tpm_r1))
          names(r0_medians) <- rownames(log2_tpm_r0)
          names(r1_medians) <- rownames(log2_tpm_r1)
          
          group_diff <- abs(r0_medians - r1_medians)
          
          # Filter by group difference
          stable_anchors <- candidate_anchors[
            group_diff[candidate_anchors] <= PANEL_CONFIG$anchor_delta_max &
              is.finite(group_diff[candidate_anchors])
          ]
          
          cat(sprintf("  Stable anchors (|∆| ≤ %.3f): %d\n", 
                      PANEL_CONFIG$anchor_delta_max, length(stable_anchors)))
          
          if (length(stable_anchors) > 0) {
            # Calculate MAD for stable anchors
            anchor_mad <- apply(log2_tpm[stable_anchors, , drop = FALSE], 1, mad)
            
            # Remove any genes with NA MAD
            valid_mad <- !is.na(anchor_mad) & is.finite(anchor_mad)
            stable_anchors <- stable_anchors[valid_mad]
            anchor_mad <- anchor_mad[valid_mad]
            
            if (length(stable_anchors) == 0) {
              cat("  WARNING: No anchors with valid MAD values.\n")
              qc_anchors <- character(0)
              anchor_df <- data.frame()
            } else {
              # Select top K by lowest MAD
              n_anchors <- min(length(stable_anchors), PANEL_CONFIG$anchor_count_target)
              qc_anchors <- names(sort(anchor_mad))[1:n_anchors]
              
              cat(sprintf("  QC anchors selected: %d\n", length(qc_anchors)))
              cat(sprintf("  MAD range: %.3f - %.3f\n", 
                          min(anchor_mad[qc_anchors]), max(anchor_mad[qc_anchors])))
              
              # Create anchor summary
              anchor_df <- data.frame(
                gene_id = qc_anchors,
                median_expression = gene_medians[qc_anchors],
                group_diff = group_diff[qc_anchors],
                mad = anchor_mad[qc_anchors],
                stringsAsFactors = FALSE
              )
            }
          } else {
            cat("  WARNING: No stable anchors found. Consider adjusting thresholds.\n")
            qc_anchors <- character(0)
            anchor_df <- data.frame()
          }
        }
      }
    }
  } else {
    cat("  WARNING: Cannot identify non-DEGs for anchors (DEG results not found).\n")
    qc_anchors <- character(0)
    anchor_df <- data.frame()
  }
}

# --------------------------------------------------------------------------
# 5.6 Create panel classification function
# --------------------------------------------------------------------------

cat("\n--- 5.6 Creating classification function ---\n")

# Classification function for new samples
classify_sample <- function(sample_reo, panel_info, panel_config) {
  N <- panel_info$N
  T <- panel_info$T
  r0_directions <- panel_info$r0_directions
  dead_zone <- panel_info$dead_zone  # From panel_info (self-contained)
  n_eff_threshold <- panel_config$n_eff_threshold  # From PANEL_CONFIG
  
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
      boundary_flag = NA,
      decision_margin = NA,
      boundary_type = "undetermined",
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
  
  # Boundary analysis
  decision_margin <- abs(n_reversals - T)
  boundary_flag <- (n_reversals == T - 1) || (n_reversals == T)
  
  # Classify based on threshold
  if (n_reversals >= T) {
    classification <- "R1-like"
  } else {
    classification <- "R0-like"
  }
  
  # Determine boundary type
  if (n_reversals < T - 1) {
    boundary_type <- "clear_R0"
  } else if (boundary_flag) {
    boundary_type <- "boundary"
  } else {
    boundary_type <- "clear_R1"
  }
  
  return(list(
    classification = classification,
    reo_score = reo_score,
    n_reversals = n_reversals,
    n_eff = n_eff,
    boundary_flag = boundary_flag,
    decision_margin = decision_margin,
    boundary_type = boundary_type,
    reason = ifelse(boundary_flag, 
                    sprintf("In boundary zone (%d-%d reversals)", T-1, T),
                    "Clear classification")
  ))
}

cat("  Classification function created.\n")
cat(sprintf("  Decision rule: ≥%d reversals → R1-like, <%d → R0-like\n", T, T))
cat(sprintf("  Boundary zone: %d-%d reversals (flagged for sensitivity analysis)\n", T-1, T))
cat("  Usage: classify_sample(sample_reo, panel_info, panel_config)\n")

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
    dead_zone = CONFIG$dead_zone,  # Include dead_zone for self-contained classification
    gene_up = panel_pairs$gene_up,
    gene_down = panel_pairs$gene_down,
    panel_config = PANEL_CONFIG,  # Include panel configuration
    boundary_zone = list(
      lower = T - 1,
      upper = T,
      label = sprintf("%d-%d reversals", T-1, T)
    )
  ),
  
  # Summary statistics
  summary_stats = list(
    n_input_pairs = nrow(final_pairs),
    n_quality_pairs = nrow(quality_pairs),
    n_panel_pairs = N,
    threshold_T = T,
    threshold_ratio = PANEL_CONFIG$t_ratio,
    n_r0_valid = nrow(r0_valid),
    n_r1_valid = nrow(r1_valid),
    r0_mean_score = ifelse(nrow(r0_valid) > 0, mean(r0_valid$reo_score), NA),
    r1_mean_score = ifelse(nrow(r1_valid) > 0, mean(r1_valid$reo_score), NA),
    n_qc_anchors = length(qc_anchors),
    # Boundary zone statistics
    boundary_zone = c(T - 1, T),
    boundary_label = sprintf("%d-%d reversals", T - 1, T),
    n_r0_boundary = ifelse(nrow(r0_valid) > 0, 
                           sum(r0_valid$boundary_flag, na.rm = TRUE), 
                           NA),
    n_r1_boundary = ifelse(nrow(r1_valid) > 0,
                           sum(r1_valid$boundary_flag, na.rm = TRUE),
                           NA),
    mean_r0_margin = ifelse(nrow(r0_valid) > 0,
                            mean(r0_valid$decision_margin),
                            NA),
    mean_r1_margin = ifelse(nrow(r1_valid) > 0,
                            mean(r1_valid$decision_margin),
                            NA)
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
  missing_rate = panel_pairs$missing_rate,
  # Add boundary zone info
  threshold_T = T,
  boundary_lower = T - 1,
  boundary_upper = T
)
write.csv(panel_def, 
          paste0(paths$output, "step5_reo_panel_definition.csv"),
          row.names = FALSE)
cat("  Panel definition exported to: step5_reo_panel_definition.csv\n")
cat(sprintf("    (includes boundary zone: %d-%d reversals)\n", T-1, T))

# Export sample scores
write.csv(all_scores,
          paste0(paths$output, "step5_sample_scores.csv"),
          row.names = FALSE)
cat("  Sample scores exported to: step5_sample_scores.csv\n")
cat("    (includes boundary_flag, decision_margin, and boundary_type columns)\n")

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
            T, PANEL_CONFIG$t_ratio * 100))
cat(sprintf("R0 valid judgments: %d/%d\n", nrow(r0_valid), length(r0_samples)))
cat(sprintf("R1 valid judgments: %d/%d\n", nrow(r1_valid), length(r1_samples)))
cat(sprintf("QC anchors selected: %d\n", length(qc_anchors)))
cat(sprintf("Boundary zone: %d-%d reversals (%.0f%%-%.0f%% of panel)\n",
            T - 1, T, (T - 1) / N * 100, T / N * 100))
cat(sprintf("Boundary cases: R0=%d, R1=%d\n", 
            ifelse(is.na(step5_data$summary_stats$n_r0_boundary), 0, step5_data$summary_stats$n_r0_boundary),
            ifelse(is.na(step5_data$summary_stats$n_r1_boundary), 0, step5_data$summary_stats$n_r1_boundary)))
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
  
  # Boundary zone analysis using boundary_flag
  cat("\nBoundary zone analysis:\n")
  
  # Define boundary zone
  boundary_zone <- list(
    lower = T - 1,
    upper = T,
    label = sprintf("Boundary (%d-%d reversals)", T-1, T)
  )
  
  cat(sprintf("  %s: %.0f%%-%.0f%% of panel\n",
              boundary_zone$label,
              boundary_zone$lower / N * 100, 
              boundary_zone$upper / N * 100))
  
  # Count samples in each zone using boundary_type
  r0_boundary_count <- sum(r0_valid$boundary_flag, na.rm = TRUE)
  r1_boundary_count <- sum(r1_valid$boundary_flag, na.rm = TRUE)
  
  cat(sprintf("  R0 samples in boundary: %d/%d (%.1f%%)\n",
              r0_boundary_count, nrow(r0_valid),
              r0_boundary_count / nrow(r0_valid) * 100))
  cat(sprintf("  R1 samples in boundary: %d/%d (%.1f%%)\n",
              r1_boundary_count, nrow(r1_valid),
              r1_boundary_count / nrow(r1_valid) * 100))
  
  # Mean decision margins
  cat("\nDecision margins (distance from threshold):\n")
  cat(sprintf("  R0 mean margin: %.1f reversals\n", mean(r0_valid$decision_margin)))
  cat(sprintf("  R1 mean margin: %.1f reversals\n", mean(r1_valid$decision_margin)))
  
  # Detailed distribution using boundary_type
  cat("\nDetailed classification zones:\n")
  r0_types <- table(r0_valid$boundary_type)
  r1_types <- table(r1_valid$boundary_type)
  
  # Ensure all types are represented
  all_types <- c("clear_R0", "boundary", "clear_R1")
  for (type in all_types) {
    r0_count <- ifelse(type %in% names(r0_types), r0_types[type], 0)
    r1_count <- ifelse(type %in% names(r1_types), r1_types[type], 0)
    cat(sprintf("  %s: R0=%d, R1=%d\n", type, r0_count, r1_count))
  }
  
  # Sensitivity preview (excluding boundary cases)
  r0_clear <- r0_valid[!r0_valid$boundary_flag, ]
  r1_clear <- r1_valid[!r1_valid$boundary_flag, ]
  
  if (nrow(r0_clear) > 0 && nrow(r1_clear) > 0) {
    cat("\nPerformance excluding boundary cases:\n")
    r0_clear_correct <- sum(r0_clear$n_reversals < T)
    r1_clear_correct <- sum(r1_clear$n_reversals >= T)
    cat(sprintf("  R0 accuracy (clear cases): %d/%d (%.1f%%)\n",
                r0_clear_correct, nrow(r0_clear),
                r0_clear_correct / nrow(r0_clear) * 100))
    cat(sprintf("  R1 accuracy (clear cases): %d/%d (%.1f%%)\n",
                r1_clear_correct, nrow(r1_clear),
                r1_clear_correct / nrow(r1_clear) * 100))
  }
}

# Final message
if (N >= PANEL_CONFIG$size_min) {
  cat(sprintf("\n✓ Panel created successfully. Ready for Step 6: R_test validation\n"))
} else {
  cat(sprintf("\n⚠ Panel size %d is below minimum. Review data quality.\n", N))
}