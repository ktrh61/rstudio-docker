# 23_reo_step3b_proxy.R - REO Step 3b: Proxy Check (Fusion Subtype)
# Purpose: Evaluate pairs for fusion subtype-specific effects (CCDC6 vs non-CCDC6)
# Input: reo_step3a_data.rds, thyr_case_master_stage2_filtered.rds
# Output: reo_step3b_data.rds
# Version: v1.2 - Complete transparent mode with CI boundary fixes
# Date: 2025-01-26
# Note: Modified criteria - R0 allows 20% exceptions, R1 uses Wilson CI with [0,1] clipping

source("analysis_v7/setup.R")

cat("\n=== REO STEP 3b: Proxy Check (Fusion Subtype) ===\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(SummarizedExperiment)
})

# --------------------------------------------------------------------------
# 3b.1 Load Step 3a data
# --------------------------------------------------------------------------

cat("--- 3b.1 Loading Step 3a data ---\n")

step3a_data_path <- paste0(paths$processed, "reo_step3a_data.rds")
if (!file.exists(step3a_data_path)) {
  stop("ERROR: reo_step3a_data.rds not found. Please run 22_reo_step3a_basic.R first.")
}

step3a_data <- readRDS(step3a_data_path)
cat("  Step 3a data loaded successfully\n")

# 薄い契約：最小限の存在チェック
stopifnot(all(c("passed_pairs", "reo_r0_passed", "reo_r1_passed", "config") %in% names(step3a_data)))

# CONFIGはRDSから取得
CONFIG <- step3a_data$config

# Extract data
passed_pairs <- step3a_data$passed_pairs
reo_r0_passed <- step3a_data$reo_r0_passed
reo_r1_passed <- step3a_data$reo_r1_passed
r0_samples <- step3a_data$r0_samples
r1_samples <- step3a_data$r1_samples

cat(sprintf("  Input pairs from Step 3a: %d\n", nrow(passed_pairs)))

# --------------------------------------------------------------------------
# 3b.2 Load case master for fusion subtype information
# --------------------------------------------------------------------------

cat("\n--- 3b.2 Loading case master for fusion subtypes ---\n")

case_master_path <- paste0(paths$processed, "thyr_case_master_stage2_filtered.rds")
if (!file.exists(case_master_path)) {
  stop("ERROR: thyr_case_master_stage2_filtered.rds not found. Please run 06_purity_analysis.R first.")
}

case_master <- readRDS(case_master_path)
cat("  Case master loaded successfully\n")

# Filter for RET cases in our analysis (with purity filters)
ret_cases <- case_master %>%
  filter(driver == "RET" & 
           group %in% c("R0", "R1") &
           has_outlier_tumor == 0 & 
           has_outlier_normal == 0 & 
           low_purity == 0)

# Extract fusion subtype from Designated_Driver column
# Expected format: "CCDC6.FusionRET" or similar
ret_cases$fusion_partner <- gsub("\\.Fusion.*", "", ret_cases$Designated_Driver)
ret_cases$is_ccdc6 <- grepl("CCDC6", ret_cases$fusion_partner, ignore.case = TRUE)

# Count by subtype
cat("\n  Fusion partner distribution:\n")
fusion_table <- table(ret_cases$fusion_partner, useNA = "ifany")
print(fusion_table)

cat(sprintf("\n  CCDC6-RET: %d cases\n", sum(ret_cases$is_ccdc6)))
cat(sprintf("  Non-CCDC6-RET: %d cases\n", sum(!ret_cases$is_ccdc6)))

# --------------------------------------------------------------------------
# 3b.3 Create sample-to-subtype mapping
# --------------------------------------------------------------------------

cat("\n--- 3b.3 Mapping samples to fusion subtypes ---\n")

# Get case IDs by subtype and group
r0_cases_ccdc6 <- ret_cases %>% 
  filter(group == "R0" & is_ccdc6) %>% 
  pull(case_submitter_id)
r0_cases_non_ccdc6 <- ret_cases %>% 
  filter(group == "R0" & !is_ccdc6) %>% 
  pull(case_submitter_id)

r1_cases_ccdc6 <- ret_cases %>% 
  filter(group == "R1" & is_ccdc6) %>% 
  pull(case_submitter_id)
r1_cases_non_ccdc6 <- ret_cases %>% 
  filter(group == "R1" & !is_ccdc6) %>% 
  pull(case_submitter_id)

# Load sample metadata to map case IDs to sample IDs
se_path <- paste0(paths$processed, "thyr_se_strand2_nonzero.rds")
se <- readRDS(se_path)
sample_metadata <- as.data.frame(colData(se))

# Map samples to subtypes
r0_samples_ccdc6 <- sample_metadata %>%
  filter(case_submitter_id %in% r0_cases_ccdc6 &
           sample_type == "Primary Tumor" &
           sample_submitter_id %in% r0_samples) %>%
  pull(sample_submitter_id)

r0_samples_non_ccdc6 <- sample_metadata %>%
  filter(case_submitter_id %in% r0_cases_non_ccdc6 &
           sample_type == "Primary Tumor" &
           sample_submitter_id %in% r0_samples) %>%
  pull(sample_submitter_id)

r1_samples_ccdc6 <- sample_metadata %>%
  filter(case_submitter_id %in% r1_cases_ccdc6 &
           sample_type == "Primary Tumor" &
           sample_submitter_id %in% r1_samples) %>%
  pull(sample_submitter_id)

r1_samples_non_ccdc6 <- sample_metadata %>%
  filter(case_submitter_id %in% r1_cases_non_ccdc6 &
           sample_type == "Primary Tumor" &
           sample_submitter_id %in% r1_samples) %>%
  pull(sample_submitter_id)

# Report sample counts by subtype
cat("\nSample counts by fusion subtype:\n")
cat(sprintf("  R0 CCDC6: %d samples\n", length(r0_samples_ccdc6)))
cat(sprintf("  R0 non-CCDC6: %d samples\n", length(r0_samples_non_ccdc6)))
cat(sprintf("  R1 CCDC6: %d samples\n", length(r1_samples_ccdc6)))
cat(sprintf("  R1 non-CCDC6: %d samples\n", length(r1_samples_non_ccdc6)))

# Verify total samples match
total_mapped <- length(r0_samples_ccdc6) + length(r0_samples_non_ccdc6) + 
  length(r1_samples_ccdc6) + length(r1_samples_non_ccdc6)
total_expected <- length(r0_samples) + length(r1_samples)

if (total_mapped != total_expected) {
  warning(sprintf("WARNING: Sample mapping mismatch. Expected %d, got %d", 
                  total_expected, total_mapped))
}

# --------------------------------------------------------------------------
# 3b.4 Apply Proxy Check criteria (Modified version)
# --------------------------------------------------------------------------

cat("\n--- 3b.4 Applying Proxy Check criteria (Modified) ---\n")
cat("  R0: exceptions ≤ min(max(1, ceil(0.2*n)), 3)\n")
cat("  R1: Wilson CI lower ≥ 0.6\n")
cat("  Relaxed logic: One subtype OK + other not contradicting (CI upper > 0.5)\n")
cat("  Extreme case filter: n≥5 and one≥0.9 and other≤0.1\n\n")

# Initialize results
proxy_results <- data.frame(
  pair_id = passed_pairs$pair_id,
  # R0 subtype evaluations
  r0_ccdc6_n = NA,
  r0_ccdc6_exceptions = NA,
  r0_ccdc6_exception_rate = NA,
  r0_ccdc6_pass = NA,
  r0_non_ccdc6_n = NA,
  r0_non_ccdc6_exceptions = NA,
  r0_non_ccdc6_exception_rate = NA,
  r0_non_ccdc6_pass = NA,
  # R1 subtype evaluations  
  r1_ccdc6_n = NA,
  r1_ccdc6_reversal_count = NA,
  r1_ccdc6_reversal_rate = NA,
  r1_ccdc6_wilson_lower = NA,
  r1_ccdc6_wilson_upper = NA,  # Add upper bound column
  r1_ccdc6_pass = NA,
  r1_non_ccdc6_n = NA,
  r1_non_ccdc6_reversal_count = NA,
  r1_non_ccdc6_reversal_rate = NA,
  r1_non_ccdc6_wilson_lower = NA,
  r1_non_ccdc6_wilson_upper = NA,  # Add upper bound column
  r1_non_ccdc6_pass = NA,
  # Overall proxy status
  extreme_case_flag = FALSE,
  proxy_status = NA_character_,
  proxy_pass = FALSE,
  pass_mode = NA_character_,  # Add pass mode tracking
  stringsAsFactors = FALSE
)

# Process each pair
cat("  Processing proxy check for each pair...\n")
pb <- txtProgressBar(min = 0, max = nrow(passed_pairs), style = 3)

for (i in 1:nrow(passed_pairs)) {
  # Get REO values for this pair
  r0_values <- reo_r0_passed[i, ]
  r1_values <- reo_r1_passed[i, ]
  
  # Get R0 direction from Step 3a results
  r0_direction <- passed_pairs$r0_direction[i]
  
  # === R0 CCDC6 evaluation ===
  if (length(r0_samples_ccdc6) >= 3) {
    r0_ccdc6_values <- r0_values[r0_samples_ccdc6]
    # Apply dead zone
    r0_ccdc6_valid <- abs(r0_ccdc6_values) >= CONFIG$dead_zone
    r0_ccdc6_valid_values <- r0_ccdc6_values[r0_ccdc6_valid]
    n_r0_ccdc6_valid <- length(r0_ccdc6_valid_values)
    
    if (n_r0_ccdc6_valid >= 3) {
      # Count exceptions (opposite to expected direction)
      if (r0_direction == 1) {
        r0_ccdc6_exceptions <- sum(r0_ccdc6_valid_values < 0)
      } else {
        r0_ccdc6_exceptions <- sum(r0_ccdc6_valid_values > 0)
      }
      
      # Modified threshold: min(max(1, ceil(0.2*n)), 3)
      max_allowed_exceptions <- min(max(1, ceiling(0.2 * n_r0_ccdc6_valid)), 3)
      r0_ccdc6_pass <- r0_ccdc6_exceptions <= max_allowed_exceptions
      
      proxy_results$r0_ccdc6_n[i] <- n_r0_ccdc6_valid
      proxy_results$r0_ccdc6_exceptions[i] <- r0_ccdc6_exceptions
      proxy_results$r0_ccdc6_exception_rate[i] <- r0_ccdc6_exceptions / n_r0_ccdc6_valid
      proxy_results$r0_ccdc6_pass[i] <- r0_ccdc6_pass
    }
  }
  
  # === R0 non-CCDC6 evaluation ===
  if (length(r0_samples_non_ccdc6) >= 3) {
    r0_non_ccdc6_values <- r0_values[r0_samples_non_ccdc6]
    # Apply dead zone
    r0_non_ccdc6_valid <- abs(r0_non_ccdc6_values) >= CONFIG$dead_zone
    r0_non_ccdc6_valid_values <- r0_non_ccdc6_values[r0_non_ccdc6_valid]
    n_r0_non_ccdc6_valid <- length(r0_non_ccdc6_valid_values)
    
    if (n_r0_non_ccdc6_valid >= 3) {
      # Count exceptions
      if (r0_direction == 1) {
        r0_non_ccdc6_exceptions <- sum(r0_non_ccdc6_valid_values < 0)
      } else {
        r0_non_ccdc6_exceptions <- sum(r0_non_ccdc6_valid_values > 0)
      }
      
      # Modified threshold
      max_allowed_exceptions <- min(max(1, ceiling(0.2 * n_r0_non_ccdc6_valid)), 3)
      r0_non_ccdc6_pass <- r0_non_ccdc6_exceptions <= max_allowed_exceptions
      
      proxy_results$r0_non_ccdc6_n[i] <- n_r0_non_ccdc6_valid
      proxy_results$r0_non_ccdc6_exceptions[i] <- r0_non_ccdc6_exceptions
      proxy_results$r0_non_ccdc6_exception_rate[i] <- r0_non_ccdc6_exceptions / n_r0_non_ccdc6_valid
      proxy_results$r0_non_ccdc6_pass[i] <- r0_non_ccdc6_pass
    }
  }
  
  # === R1 CCDC6 evaluation ===
  if (length(r1_samples_ccdc6) >= 3) {
    r1_ccdc6_values <- r1_values[r1_samples_ccdc6]
    # Apply dead zone
    r1_ccdc6_valid <- abs(r1_ccdc6_values) >= CONFIG$dead_zone
    r1_ccdc6_valid_values <- r1_ccdc6_values[r1_ccdc6_valid]
    n_r1_ccdc6_valid <- length(r1_ccdc6_valid_values)
    
    if (n_r1_ccdc6_valid >= 3) {
      # Count reversals (opposite to R0 direction)
      if (r0_direction == 1) {
        r1_ccdc6_reversals <- sum(r1_ccdc6_valid_values < 0)
      } else {
        r1_ccdc6_reversals <- sum(r1_ccdc6_valid_values > 0)
      }
      
      # Wilson confidence interval
      r1_ccdc6_reversal_rate <- r1_ccdc6_reversals / n_r1_ccdc6_valid
      z <- qnorm(0.975)  # 95% CI (consistent with Step 3a)
      n <- n_r1_ccdc6_valid
      p_hat <- r1_ccdc6_reversal_rate
      
      # Wilson CI calculation (corrected formula)
      denominator <- 1 + z^2/n
      center <- (p_hat + z^2/(2*n)) / denominator
      margin <- z * sqrt(p_hat * (1-p_hat) / n + z^2/(4*n^2)) / denominator
      
      ci_lower <- center - margin
      ci_upper <- center + margin  # Calculate upper bound
      
      # Clip to [0, 1] range (consistent with Step 3a)
      ci_lower <- max(0, ci_lower)
      ci_upper <- min(1, ci_upper)
      
      proxy_results$r1_ccdc6_n[i] <- n_r1_ccdc6_valid
      proxy_results$r1_ccdc6_reversal_count[i] <- r1_ccdc6_reversals
      proxy_results$r1_ccdc6_reversal_rate[i] <- r1_ccdc6_reversal_rate
      proxy_results$r1_ccdc6_wilson_lower[i] <- ci_lower
      proxy_results$r1_ccdc6_wilson_upper[i] <- ci_upper  # Store upper bound
      proxy_results$r1_ccdc6_pass[i] <- ci_lower >= 0.6
    }
  }
  
  # === R1 non-CCDC6 evaluation ===
  if (length(r1_samples_non_ccdc6) >= 3) {
    r1_non_ccdc6_values <- r1_values[r1_samples_non_ccdc6]
    # Apply dead zone
    r1_non_ccdc6_valid <- abs(r1_non_ccdc6_values) >= CONFIG$dead_zone
    r1_non_ccdc6_valid_values <- r1_non_ccdc6_values[r1_non_ccdc6_valid]
    n_r1_non_ccdc6_valid <- length(r1_non_ccdc6_valid_values)
    
    if (n_r1_non_ccdc6_valid >= 3) {
      # Count reversals
      if (r0_direction == 1) {
        r1_non_ccdc6_reversals <- sum(r1_non_ccdc6_valid_values < 0)
      } else {
        r1_non_ccdc6_reversals <- sum(r1_non_ccdc6_valid_values > 0)
      }
      
      # Wilson confidence interval
      r1_non_ccdc6_reversal_rate <- r1_non_ccdc6_reversals / n_r1_non_ccdc6_valid
      z <- qnorm(0.975)  # 95% CI (consistent with Step 3a)
      n <- n_r1_non_ccdc6_valid
      p_hat <- r1_non_ccdc6_reversal_rate
      
      # Wilson CI calculation (corrected formula)
      denominator <- 1 + z^2/n
      center <- (p_hat + z^2/(2*n)) / denominator
      margin <- z * sqrt(p_hat * (1-p_hat) / n + z^2/(4*n^2)) / denominator
      
      ci_lower <- center - margin
      ci_upper <- center + margin  # Calculate upper bound
      
      # Clip to [0, 1] range (consistent with Step 3a)
      ci_lower <- max(0, ci_lower)
      ci_upper <- min(1, ci_upper)
      
      proxy_results$r1_non_ccdc6_n[i] <- n_r1_non_ccdc6_valid
      proxy_results$r1_non_ccdc6_reversal_count[i] <- r1_non_ccdc6_reversals
      proxy_results$r1_non_ccdc6_reversal_rate[i] <- r1_non_ccdc6_reversal_rate
      proxy_results$r1_non_ccdc6_wilson_lower[i] <- ci_lower
      proxy_results$r1_non_ccdc6_wilson_upper[i] <- ci_upper  # Store upper bound
      proxy_results$r1_non_ccdc6_pass[i] <- ci_lower >= 0.6
    }
  }
  
  # === Check extreme case (safety valve) ===
  # Only apply if both subtypes have n≥5
  if (!is.na(proxy_results$r1_ccdc6_n[i]) && !is.na(proxy_results$r1_non_ccdc6_n[i]) &&
      proxy_results$r1_ccdc6_n[i] >= 5 && proxy_results$r1_non_ccdc6_n[i] >= 5) {
    
    rate_ccdc6 <- proxy_results$r1_ccdc6_reversal_rate[i]
    rate_non_ccdc6 <- proxy_results$r1_non_ccdc6_reversal_rate[i]
    
    if ((rate_ccdc6 >= 0.9 && rate_non_ccdc6 <= 0.1) ||
        (rate_ccdc6 <= 0.1 && rate_non_ccdc6 >= 0.9)) {
      proxy_results$extreme_case_flag[i] <- TRUE
    }
  }
  
  # === Determine overall proxy status (COMPLETE TRANSPARENT MODE LOGIC) ===
  ccdc6_r0_evaluable <- !is.na(proxy_results$r0_ccdc6_pass[i])
  ccdc6_r1_evaluable <- !is.na(proxy_results$r1_ccdc6_pass[i])
  non_ccdc6_r0_evaluable <- !is.na(proxy_results$r0_non_ccdc6_pass[i])
  non_ccdc6_r1_evaluable <- !is.na(proxy_results$r1_non_ccdc6_pass[i])
  
  ccdc6_evaluable <- ccdc6_r0_evaluable || ccdc6_r1_evaluable
  non_ccdc6_evaluable <- non_ccdc6_r0_evaluable || non_ccdc6_r1_evaluable
  
  # Check extreme case first
  if (proxy_results$extreme_case_flag[i]) {
    proxy_results$proxy_status[i] <- "extreme_case_fail"
    proxy_results$proxy_pass[i] <- FALSE
    proxy_results$pass_mode[i] <- NA_character_  # No pass mode for failures
    
  } else if (!ccdc6_evaluable && !non_ccdc6_evaluable) {
    # Both subtypes not evaluable (n<3 for both)
    proxy_results$proxy_status[i] <- "skipped_both"
    proxy_results$proxy_pass[i] <- TRUE  # Continue to adoption
    proxy_results$pass_mode[i] <- "skipped_both"
    
  } else if (ccdc6_evaluable && !non_ccdc6_evaluable) {
    # Only CCDC6 evaluable
    ccdc6_r0_ok <- !ccdc6_r0_evaluable || proxy_results$r0_ccdc6_pass[i]
    ccdc6_r1_ok <- !ccdc6_r1_evaluable || proxy_results$r1_ccdc6_pass[i]
    
    if (ccdc6_r0_ok && ccdc6_r1_ok) {
      proxy_results$proxy_status[i] <- "one_pass_one_na"
      proxy_results$proxy_pass[i] <- TRUE
      proxy_results$pass_mode[i] <- "one_pass_one_na"
    } else {
      proxy_results$proxy_status[i] <- "fail_on_one"
      proxy_results$proxy_pass[i] <- FALSE
      proxy_results$pass_mode[i] <- NA_character_
    }
    
  } else if (!ccdc6_evaluable && non_ccdc6_evaluable) {
    # Only non-CCDC6 evaluable
    non_ccdc6_r0_ok <- !non_ccdc6_r0_evaluable || proxy_results$r0_non_ccdc6_pass[i]
    non_ccdc6_r1_ok <- !non_ccdc6_r1_evaluable || proxy_results$r1_non_ccdc6_pass[i]
    
    if (non_ccdc6_r0_ok && non_ccdc6_r1_ok) {
      proxy_results$proxy_status[i] <- "one_pass_one_na"
      proxy_results$proxy_pass[i] <- TRUE
      proxy_results$pass_mode[i] <- "one_pass_one_na"
    } else {
      proxy_results$proxy_status[i] <- "fail_on_one"
      proxy_results$proxy_pass[i] <- FALSE
      proxy_results$pass_mode[i] <- NA_character_
    }
    
  } else {
    # Both subtypes evaluable - use relaxed criteria
    ccdc6_r0_ok <- !ccdc6_r0_evaluable || proxy_results$r0_ccdc6_pass[i]
    ccdc6_r1_ok <- !ccdc6_r1_evaluable || proxy_results$r1_ccdc6_pass[i]
    non_ccdc6_r0_ok <- !non_ccdc6_r0_evaluable || proxy_results$r0_non_ccdc6_pass[i]
    non_ccdc6_r1_ok <- !non_ccdc6_r1_evaluable || proxy_results$r1_non_ccdc6_pass[i]
    
    # Combined OK status for each subtype
    ccdc6_both_ok <- ccdc6_r0_ok && ccdc6_r1_ok
    non_ccdc6_both_ok <- non_ccdc6_r0_ok && non_ccdc6_r1_ok
    
    # Check for contradictions (R1 CI upper ≤ 0.5 means statistically confident non-reversal)
    ccdc6_contra <- FALSE
    if (ccdc6_r1_evaluable && !is.na(proxy_results$r1_ccdc6_wilson_upper[i])) {
      ccdc6_contra <- proxy_results$r1_ccdc6_wilson_upper[i] <= 0.5
    }
    
    non_ccdc6_contra <- FALSE
    if (non_ccdc6_r1_evaluable && !is.na(proxy_results$r1_non_ccdc6_wilson_upper[i])) {
      non_ccdc6_contra <- proxy_results$r1_non_ccdc6_wilson_upper[i] <= 0.5
    }
    
    # Apply complete pass mode logic
    if (ccdc6_both_ok && non_ccdc6_both_ok) {
      # Both subtypes meet all criteria
      proxy_results$proxy_status[i] <- "pass"
      proxy_results$proxy_pass[i] <- TRUE
      proxy_results$pass_mode[i] <- "both_ok"
      
    } else if ((ccdc6_both_ok && !non_ccdc6_contra) || (non_ccdc6_both_ok && !ccdc6_contra)) {
      # One subtype meets criteria, other not contradicting
      proxy_results$proxy_status[i] <- "pass"
      proxy_results$proxy_pass[i] <- TRUE
      proxy_results$pass_mode[i] <- "one_ok_not_contra"
      
    } else {
      # Neither condition met - fail
      proxy_results$proxy_status[i] <- "fail"
      proxy_results$proxy_pass[i] <- FALSE
      proxy_results$pass_mode[i] <- NA_character_
    }
  }
  
  setTxtProgressBar(pb, i)
}
close(pb)

# --------------------------------------------------------------------------
# 3b.5 Summary of results
# --------------------------------------------------------------------------

cat("\n--- 3b.5 Summary of proxy check results ---\n")

# Status distribution
cat("\nProxy status distribution:\n")
status_table <- table(proxy_results$proxy_status)
for (status in names(status_table)) {
  cat(sprintf("  %s: %d pairs\n", status, status_table[status]))
}

# Pass mode distribution (for passed pairs only)
cat("\nPass mode distribution (passed pairs only):\n")
pass_mode_table <- table(proxy_results$pass_mode[proxy_results$proxy_pass])
for (mode in names(pass_mode_table)) {
  cat(sprintf("  %s: %d pairs\n", mode, pass_mode_table[mode]))
}

# Overall pass rate
cat(sprintf("\nOverall proxy pass rate: %d / %d (%.1f%%)\n",
            sum(proxy_results$proxy_pass),
            nrow(proxy_results),
            sum(proxy_results$proxy_pass) / nrow(proxy_results) * 100))

# Extreme cases
if (sum(proxy_results$extreme_case_flag) > 0) {
  cat(sprintf("  Extreme cases detected: %d pairs\n", 
              sum(proxy_results$extreme_case_flag)))
}

# Filter to passed pairs
proxy_passed_pairs <- passed_pairs[proxy_results$proxy_pass, ]
cat(sprintf("\n✓ Pairs passing Step 3b: %d\n", nrow(proxy_passed_pairs)))

# --------------------------------------------------------------------------
# 3b.6 Save Step 3b results
# --------------------------------------------------------------------------

cat("\n--- 3b.6 Saving Step 3b results ---\n")

# Combine with Step 3a results
step3b_results <- cbind(
  passed_pairs,
  proxy_results[, -1]  # Remove duplicate pair_id
)

step3b_data <- list(
  # Configuration (pass through)
  config = CONFIG,
  
  # Full results
  all_results = step3b_results,
  passed_pairs = step3b_results[step3b_results$proxy_pass, ],
  
  # REO matrices for passed pairs (subset from 3a)
  reo_r0_passed = reo_r0_passed[proxy_results$proxy_pass, ],
  reo_r1_passed = reo_r1_passed[proxy_results$proxy_pass, ],
  
  # Sample subtype mapping (for reference)
  sample_mapping = list(
    r0_ccdc6 = r0_samples_ccdc6,
    r0_non_ccdc6 = r0_samples_non_ccdc6,
    r1_ccdc6 = r1_samples_ccdc6,
    r1_non_ccdc6 = r1_samples_non_ccdc6
  ),
  
  # Summary statistics
  n_input = nrow(passed_pairs),
  n_passed = sum(proxy_results$proxy_pass),
  pass_rate = sum(proxy_results$proxy_pass) / nrow(passed_pairs),
  
  # Pass mode summary
  pass_mode_summary = table(proxy_results$pass_mode[proxy_results$proxy_pass]),
  
  # Sample info
  r0_samples = r0_samples,
  r1_samples = r1_samples,
  
  # Metadata
  timestamp = Sys.time(),
  step = "Step 3b Complete"
)

saveRDS(step3b_data, paste0(paths$processed, "reo_step3b_data.rds"))
cat("  Step 3b data saved to: reo_step3b_data.rds\n")

# Export detailed results to CSV for review
write.csv(step3b_results,
          paste0(paths$output, "step3b_proxy_check_results.csv"),
          row.names = FALSE)
cat("  Full proxy check results exported to: step3b_proxy_check_results.csv\n")

# Export passed pairs only
write.csv(step3b_results[step3b_results$proxy_pass, ],
          paste0(paths$output, "step3b_passed_pairs.csv"),
          row.names = FALSE)
cat("  Passed pairs exported to: step3b_passed_pairs.csv\n")

# =============================================================================
# STEP 3b SUMMARY (固定フォーマット)
# =============================================================================

cat("\n=== STEP 3b SUMMARY ===\n")
cat(sprintf("Input pairs: %d\n", nrow(passed_pairs)))
cat(sprintf("Proxy pass: %d (%.1f%%)\n",
            sum(proxy_results$proxy_pass),
            sum(proxy_results$proxy_pass) / nrow(proxy_results) * 100))
cat(sprintf("Skipped (both n<3): %d\n", 
            sum(proxy_results$proxy_status == "skipped_both", na.rm = TRUE)))
cat(sprintf("Extreme case fail: %d\n",
            sum(proxy_results$extreme_case_flag)))

# Display pass mode breakdown
cat("\nPass modes:\n")
if (length(pass_mode_table) > 0) {
  for (mode in names(pass_mode_table)) {
    cat(sprintf("  %s: %d\n", mode, pass_mode_table[mode]))
  }
}

cat(sprintf("Timestamp: %s\n", format(step3b_data$timestamp)))

# Final message
if (step3b_data$n_passed > 0) {
  cat(sprintf("\n✓ Ready for Step 4: Redundancy Removal (%d pairs available)\n", 
              step3b_data$n_passed))
} else {
  cat("\n⚠ No pairs passed Step 3b. Review parameters or data quality.\n")
}