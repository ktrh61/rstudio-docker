# 22_reo_step3a_basic.R - REO Step 3a: Adoption Criteria (Basic Judgment)
# Purpose: Apply five judgment criteria to candidate pairs
# Input: reo_step2_data.rds
# Output: reo_step3a_data.rds
# Version: v1.0
# Date: 2025-01-26
# Note: Wilson CI calculation corrected

source("analysis_v7/setup.R")

cat("\n=== REO STEP 3a: Adoption Criteria (Basic Judgment) ===\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
})

# --------------------------------------------------------------------------
# 3a.1 Load Step 2 data and calculate REO values
# --------------------------------------------------------------------------

cat("--- 3a.1 Loading data and calculating REO values ---\n")

step2_data_path <- paste0(paths$processed, "reo_step2_data.rds")
if (!file.exists(step2_data_path)) {
  stop("ERROR: reo_step2_data.rds not found. Please run 21_reo_step2_pairs.R first.")
}

step2_data <- readRDS(step2_data_path)
cat("  Step 2 data loaded successfully\n")

# 薄い契約：最小限の存在チェック
stopifnot(all(c("pairs_df", "log2_tpm_r0", "log2_tpm_r1", "config") %in% names(step2_data)))

# CONFIGはRDSから取得（reo_config.Rは読まない）
CONFIG <- step2_data$config

# Extract data
pairs_df <- step2_data$pairs_df
log2_tpm_r0 <- step2_data$log2_tpm_r0
log2_tpm_r1 <- step2_data$log2_tpm_r1

cat(sprintf("  Loaded %d candidate pairs\n", nrow(pairs_df)))
cat(sprintf("  R0 samples: %d, R1 samples: %d\n", 
            ncol(log2_tpm_r0), ncol(log2_tpm_r1)))

# Calculate REO values (r = log2(up/down) = log2(up) - log2(down))
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

cat(sprintf("  REO matrices created: %d pairs × %d R0 + %d R1 samples\n",
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
  
  # ===== Criterion 3: R1 Reversal (Wilson CI - CORRECTED) =====
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
      # CORRECTED: sqrt calculated first, then divided by denominator
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

cat(sprintf("  2. Strength (a[2]≥%.2f): %d / %d (%.1f%%)\n",
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
cat(sprintf("\n✓ Pairs passing Step 3a: %d\n", nrow(passed_pairs)))

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
  # Configuration (pass through)
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

# =============================================================================
# STEP 3a SUMMARY (固定フォーマット)
# =============================================================================

cat("\n=== STEP 3a SUMMARY ===\n")
cat(sprintf("R0 directionality: %d / %d (%.1f%%)\n",
            sum(judgment_results$criterion1_pass),
            nrow(judgment_results),
            sum(judgment_results$criterion1_pass) / nrow(judgment_results) * 100))
cat(sprintf("Strength (a[2]≥τ): %d / %d (%.1f%%)\n",
            sum(judgment_results$criterion2_pass),
            nrow(judgment_results),
            sum(judgment_results$criterion2_pass) / nrow(judgment_results) * 100))
cat(sprintf("R1 reversal (CI>0.5): %d / %d (%.1f%%)\n",
            sum(judgment_results$criterion3_pass),
            nrow(judgment_results),
            sum(judgment_results$criterion3_pass) / nrow(judgment_results) * 100))
cat(sprintf("Not 100%% reversal: %d / %d (%.1f%%)\n",
            sum(judgment_results$criterion4_pass),
            nrow(judgment_results),
            sum(judgment_results$criterion4_pass) / nrow(judgment_results) * 100))
cat(sprintf("Overall pass: %d / %d (%.1f%%)\n",
            sum(judgment_results$all_criteria_pass),
            nrow(judgment_results),
            sum(judgment_results$all_criteria_pass) / nrow(judgment_results) * 100))