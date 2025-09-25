# 23_reo_step3b_proxy.R - REO Step 3b: Proxy Check (Fusion Subtype) + POC Trend Diagnosis
# Purpose: Evaluate pairs for fusion subtype-specific effects and add POC trend diagnostic tags
# Input: reo_step3a_data.rds, thyr_case_master_stage2_filtered.rds
# Output: reo_step3b_data.rds
# Version: v2.1 - Modified CR1 to use Bayesian criterion with Jeffreys prior
# Date: 2025-01-27
# Note: POC trend tags are diagnostic only - not used for pair selection/judgment
# Change: CR1 now uses Pr(p > 0.5 | data) >= 0.95 with Jeffreys prior Beta(0.5, 0.5)
#         instead of Wilson CI lower >= 0.6 to better accommodate small samples

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

# Get sample metadata from SummarizedExperiment
se_path <- paste0(paths$processed, "thyr_se_strand2_nonzero.rds")
if (!file.exists(se_path)) {
  stop("ERROR: thyr_se_strand2_nonzero.rds not found. Please run 03_prepare_counts.R first.")
}
se <- readRDS(se_path)
sample_metadata <- as.data.frame(colData(se))

# Map samples to fusion subtypes
r0_sample_info <- sample_metadata %>%
  filter(sample_submitter_id %in% r0_samples) %>%
  left_join(case_master %>% select(case_submitter_id, Designated_Driver),
            by = "case_submitter_id")

r1_sample_info <- sample_metadata %>%
  filter(sample_submitter_id %in% r1_samples) %>%
  left_join(case_master %>% select(case_submitter_id, Designated_Driver),
            by = "case_submitter_id")

# Convert Designated_Driver to fusion partner
# Format: "CCDC6-RET", "NCOA4-RET", "RET-OTHER"
get_fusion_partner <- function(designated_driver) {
  ifelse(grepl("CCDC6", designated_driver), "CCDC6",
         ifelse(grepl("NCOA4", designated_driver), "NCOA4", "non-CCDC6"))
}

r0_sample_info$fusion_partner <- get_fusion_partner(r0_sample_info$Designated_Driver)
r1_sample_info$fusion_partner <- get_fusion_partner(r1_sample_info$Designated_Driver)

# Classify samples by subtype
r0_samples_ccdc6 <- r0_sample_info %>%
  filter(fusion_partner == "CCDC6") %>%
  pull(sample_submitter_id)

r0_samples_non_ccdc6 <- r0_sample_info %>%
  filter(fusion_partner != "CCDC6") %>%
  pull(sample_submitter_id)

r1_samples_ccdc6 <- r1_sample_info %>%
  filter(fusion_partner == "CCDC6") %>%
  pull(sample_submitter_id)

r1_samples_non_ccdc6 <- r1_sample_info %>%
  filter(fusion_partner != "CCDC6") %>%
  pull(sample_submitter_id)

cat("\nSubtype breakdown:\n")
cat(sprintf("  R0: CCDC6=%d, non-CCDC6=%d\n", 
            length(r0_samples_ccdc6), length(r0_samples_non_ccdc6)))
cat(sprintf("  R1: CCDC6=%d, non-CCDC6=%d\n", 
            length(r1_samples_ccdc6), length(r1_samples_non_ccdc6)))

# --------------------------------------------------------------------------
# 3b.3 Helper functions for proxy evaluation
# --------------------------------------------------------------------------

cat("\n--- 3b.3 Setting up proxy evaluation functions ---\n")

# Function to calculate allowed exceptions for R0
calc_max_exceptions <- function(n) {
  # Allow min(max(1, ceil(0.2*n)), 3)
  return(min(max(1, ceiling(0.2 * n)), 3))
}

# Function to calculate Bayesian posterior probability
# Using Jeffreys prior Beta(0.5, 0.5)
calc_bayes_prob_greater_half <- function(successes, n) {
  # Posterior is Beta(k + 0.5, n - k + 0.5)
  # Return Pr(p > 0.5 | data)
  return(pbeta(0.5, successes + 0.5, n - successes + 0.5, lower.tail = FALSE))
}

# Function to calculate Bayesian one-sided lower credible limit
# Using Jeffreys prior Beta(0.5, 0.5)
calc_bayes_lcl <- function(successes, n, level = 0.95) {
  # Posterior is Beta(k + 0.5, n - k + 0.5)
  # Return one-sided lower credible limit (level% from below)
  alpha <- 1 - level  # For one-sided 95% LCL, use 5% quantile
  return(qbeta(alpha, successes + 0.5, n - successes + 0.5))
}

# --------------------------------------------------------------------------
# 3b.4 Apply Proxy Check criteria (Modified with Bayesian approach)
# --------------------------------------------------------------------------

cat("\n--- 3b.4 Applying Proxy Check criteria (Modified with Bayesian) ---\n")
cat("  R0: exceptions ≤ min(max(1, ceil(0.2*n)), 3)\n")
cat("  R1: Pr(p > 0.5 | data) ≥ 0.95 AND n≥6 AND point estimate≥0.7\n")
cat("  Relaxed logic: One subtype OK + other not contradicting\n")
cat("  Non-contradiction: Pr(p > 0.5) ≥ 0.20 AND Wilson CI upper ≥ 0.55 (n≥5 required)\n")
cat("  Veto: one-sided LCL_ok ≥ 0.7 AND UCL_other ≤ 0.5 → fail (both n≥5)\n\n")

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
  # R1 subtype evaluations with Bayesian  
  r1_ccdc6_n = NA,
  r1_ccdc6_reversal_count = NA,
  r1_ccdc6_reversal_rate = NA,
  r1_ccdc6_bayes_prob = NA,  # Pr(p > 0.5 | data)
  r1_ccdc6_wilson_lower = NA,  # Keep for reference
  r1_ccdc6_wilson_upper = NA,
  r1_ccdc6_pass = NA,
  r1_non_ccdc6_n = NA,
  r1_non_ccdc6_reversal_count = NA,
  r1_non_ccdc6_reversal_rate = NA,
  r1_non_ccdc6_bayes_prob = NA,  # Pr(p > 0.5 | data)
  r1_non_ccdc6_wilson_lower = NA,  # Keep for reference
  r1_non_ccdc6_wilson_upper = NA,
  r1_non_ccdc6_pass = NA,
  # Overall proxy status
  extreme_case_flag = FALSE,
  proxy_status = NA_character_,
  proxy_pass = FALSE,
  pass_mode = NA_character_,
  pass_subtype = NA_character_,  # Track which subtype passed in one_ok_not_contra
  # POC trend diagnosis (diagnostic tags only - not used for selection)
  rho_poc_flip = NA_real_,
  p_poc_flip = NA_real_,
  n_valid_poc = NA_integer_,
  stable = NA,
  trend_tag = NA_character_,
  trend_method = NA_character_,
  partner_hetero_p = NA_real_,
  partner_k = NA_integer_,
  per_level_n = NA_character_,
  rank_penalty_trend = 0L,
  rank_penalty_partner = 0L,
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
      # Count exceptions (wrong sign)
      if (r0_direction == 1) {
        r0_ccdc6_exceptions <- sum(r0_ccdc6_valid_values < 0)
      } else {
        r0_ccdc6_exceptions <- sum(r0_ccdc6_valid_values > 0)
      }
      
      max_allowed <- calc_max_exceptions(n_r0_ccdc6_valid)
      
      proxy_results$r0_ccdc6_n[i] <- n_r0_ccdc6_valid
      proxy_results$r0_ccdc6_exceptions[i] <- r0_ccdc6_exceptions
      proxy_results$r0_ccdc6_exception_rate[i] <- r0_ccdc6_exceptions / n_r0_ccdc6_valid
      proxy_results$r0_ccdc6_pass[i] <- r0_ccdc6_exceptions <= max_allowed
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
      
      max_allowed <- calc_max_exceptions(n_r0_non_ccdc6_valid)
      
      proxy_results$r0_non_ccdc6_n[i] <- n_r0_non_ccdc6_valid
      proxy_results$r0_non_ccdc6_exceptions[i] <- r0_non_ccdc6_exceptions
      proxy_results$r0_non_ccdc6_exception_rate[i] <- r0_non_ccdc6_exceptions / n_r0_non_ccdc6_valid
      proxy_results$r0_non_ccdc6_pass[i] <- r0_non_ccdc6_exceptions <= max_allowed
    }
  }
  
  # === R1 CCDC6 evaluation (MODIFIED: Bayesian approach with stricter criteria) ===
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
      
      r1_ccdc6_reversal_rate <- r1_ccdc6_reversals / n_r1_ccdc6_valid
      
      # Bayesian approach with Jeffreys prior Beta(0.5, 0.5)
      bayes_prob <- calc_bayes_prob_greater_half(r1_ccdc6_reversals, n_r1_ccdc6_valid)
      
      # Calculate one-sided Bayesian LCL for veto check
      bayes_lcl <- calc_bayes_lcl(r1_ccdc6_reversals, n_r1_ccdc6_valid, 0.95)
      
      # Also calculate Wilson CI for reference and non-contradiction check
      z <- qnorm(0.975)  # 95% CI
      p_hat <- r1_ccdc6_reversal_rate
      n <- n_r1_ccdc6_valid
      
      denominator <- 1 + z^2/n
      center <- (p_hat + z^2/(2*n)) / denominator
      margin <- z * sqrt(p_hat * (1-p_hat) / n + z^2/(4*n^2)) / denominator
      
      ci_lower <- max(0, center - margin)
      ci_upper <- min(1, center + margin)
      
      proxy_results$r1_ccdc6_n[i] <- n_r1_ccdc6_valid
      proxy_results$r1_ccdc6_reversal_count[i] <- r1_ccdc6_reversals
      proxy_results$r1_ccdc6_reversal_rate[i] <- r1_ccdc6_reversal_rate
      proxy_results$r1_ccdc6_bayes_prob[i] <- bayes_prob
      proxy_results$r1_ccdc6_wilson_lower[i] <- ci_lower
      proxy_results$r1_ccdc6_wilson_upper[i] <- ci_upper
      
      # Stricter criteria: Pr(p > 0.5) ≥ 0.95 AND n≥6 AND point estimate≥0.7
      proxy_results$r1_ccdc6_pass[i] <- (bayes_prob >= 0.95 && 
                                           n_r1_ccdc6_valid >= 6 && 
                                           r1_ccdc6_reversal_rate >= 0.7)
    }
  }
  
  # === R1 non-CCDC6 evaluation (MODIFIED: Bayesian approach with stricter criteria) ===
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
      
      r1_non_ccdc6_reversal_rate <- r1_non_ccdc6_reversals / n_r1_non_ccdc6_valid
      
      # Bayesian approach with Jeffreys prior Beta(0.5, 0.5)
      bayes_prob <- calc_bayes_prob_greater_half(r1_non_ccdc6_reversals, n_r1_non_ccdc6_valid)
      
      # Calculate one-sided Bayesian LCL for veto check
      bayes_lcl <- calc_bayes_lcl(r1_non_ccdc6_reversals, n_r1_non_ccdc6_valid, 0.95)
      
      # Also calculate Wilson CI for reference and non-contradiction check
      z <- qnorm(0.975)
      p_hat <- r1_non_ccdc6_reversal_rate
      n <- n_r1_non_ccdc6_valid
      
      denominator <- 1 + z^2/n
      center <- (p_hat + z^2/(2*n)) / denominator
      margin <- z * sqrt(p_hat * (1-p_hat) / n + z^2/(4*n^2)) / denominator
      
      ci_lower <- max(0, center - margin)
      ci_upper <- min(1, center + margin)
      
      proxy_results$r1_non_ccdc6_n[i] <- n_r1_non_ccdc6_valid
      proxy_results$r1_non_ccdc6_reversal_count[i] <- r1_non_ccdc6_reversals
      proxy_results$r1_non_ccdc6_reversal_rate[i] <- r1_non_ccdc6_reversal_rate
      proxy_results$r1_non_ccdc6_bayes_prob[i] <- bayes_prob
      proxy_results$r1_non_ccdc6_wilson_lower[i] <- ci_lower
      proxy_results$r1_non_ccdc6_wilson_upper[i] <- ci_upper
      
      # Stricter criteria: Pr(p > 0.5) ≥ 0.95 AND n≥6 AND point estimate≥0.7
      proxy_results$r1_non_ccdc6_pass[i] <- (bayes_prob >= 0.95 && 
                                               n_r1_non_ccdc6_valid >= 6 && 
                                               r1_non_ccdc6_reversal_rate >= 0.7)
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
  
  # === Determine overall proxy status (STRICTER CRITERIA) ===
  ccdc6_r0_evaluable <- !is.na(proxy_results$r0_ccdc6_pass[i])
  ccdc6_r1_evaluable <- !is.na(proxy_results$r1_ccdc6_pass[i])
  non_ccdc6_r0_evaluable <- !is.na(proxy_results$r0_non_ccdc6_pass[i])
  non_ccdc6_r1_evaluable <- !is.na(proxy_results$r1_non_ccdc6_pass[i])
  
  ccdc6_evaluable <- ccdc6_r0_evaluable || ccdc6_r1_evaluable
  non_ccdc6_evaluable <- non_ccdc6_r0_evaluable || non_ccdc6_r1_evaluable
  
  # Check extreme case first (legacy, rarely triggers)
  if (proxy_results$extreme_case_flag[i]) {
    proxy_results$proxy_status[i] <- "extreme_case_fail"
    proxy_results$proxy_pass[i] <- FALSE
    proxy_results$pass_mode[i] <- NA_character_
    
    # === NEW: Veto check (LCL_ok ≥ 0.7 AND UCL_other ≤ 0.5) ===
  } else if (!is.na(proxy_results$r1_ccdc6_n[i]) && !is.na(proxy_results$r1_non_ccdc6_n[i]) &&
             proxy_results$r1_ccdc6_n[i] >= 5 && proxy_results$r1_non_ccdc6_n[i] >= 5) {
    
    # Calculate one-sided LCLs for veto check
    ccdc6_lcl <- calc_bayes_lcl(proxy_results$r1_ccdc6_reversal_count[i], 
                                proxy_results$r1_ccdc6_n[i], 0.95)
    non_ccdc6_lcl <- calc_bayes_lcl(proxy_results$r1_non_ccdc6_reversal_count[i], 
                                    proxy_results$r1_non_ccdc6_n[i], 0.95)
    
    ccdc6_ucl <- proxy_results$r1_ccdc6_wilson_upper[i]
    non_ccdc6_ucl <- proxy_results$r1_non_ccdc6_wilson_upper[i]
    
    # Veto if strong contradiction detected
    if ((ccdc6_lcl >= 0.7 && non_ccdc6_ucl <= 0.5) || 
        (non_ccdc6_lcl >= 0.7 && ccdc6_ucl <= 0.5)) {
      proxy_results$proxy_status[i] <- "veto_contradiction"
      proxy_results$proxy_pass[i] <- FALSE
      proxy_results$pass_mode[i] <- NA_character_
    } else {
      # Continue with regular logic (see below)
      veto_applied <- FALSE
    }
  } else {
    veto_applied <- FALSE
  }
  
  # Continue only if not vetoed
  if (is.na(proxy_results$proxy_status[i]) || proxy_results$proxy_status[i] == "") {
    
    if (!ccdc6_evaluable && !non_ccdc6_evaluable) {
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
      # Both subtypes evaluable - use stricter criteria
      ccdc6_r0_ok <- !ccdc6_r0_evaluable || proxy_results$r0_ccdc6_pass[i]
      ccdc6_r1_ok <- !ccdc6_r1_evaluable || proxy_results$r1_ccdc6_pass[i]
      non_ccdc6_r0_ok <- !non_ccdc6_r0_evaluable || proxy_results$r0_non_ccdc6_pass[i]
      non_ccdc6_r1_ok <- !non_ccdc6_r1_evaluable || proxy_results$r1_non_ccdc6_pass[i]
      
      # Combined OK status for each subtype
      ccdc6_both_ok <- ccdc6_r0_ok && ccdc6_r1_ok
      non_ccdc6_both_ok <- non_ccdc6_r0_ok && non_ccdc6_r1_ok
      
      # === STRICTER non-contradiction check ===
      # Non-contradiction requires: Pr(p > 0.5) ≥ 0.20 AND Wilson upper ≥ 0.55
      # AND n≥5 (otherwise cannot evaluate non-contradiction)
      ccdc6_contra <- FALSE
      ccdc6_can_eval_contra <- FALSE
      if (ccdc6_r1_evaluable && !is.na(proxy_results$r1_ccdc6_n[i]) && 
          proxy_results$r1_ccdc6_n[i] >= 5) {
        ccdc6_can_eval_contra <- TRUE
        # Contradiction if: Pr < 0.20 OR Wilson upper < 0.55
        if (proxy_results$r1_ccdc6_bayes_prob[i] < 0.20 || 
            proxy_results$r1_ccdc6_wilson_upper[i] < 0.55) {
          ccdc6_contra <- TRUE
        }
      }
      
      non_ccdc6_contra <- FALSE
      non_ccdc6_can_eval_contra <- FALSE
      if (non_ccdc6_r1_evaluable && !is.na(proxy_results$r1_non_ccdc6_n[i]) && 
          proxy_results$r1_non_ccdc6_n[i] >= 5) {
        non_ccdc6_can_eval_contra <- TRUE
        # Contradiction if: Pr < 0.20 OR Wilson upper < 0.55
        if (proxy_results$r1_non_ccdc6_bayes_prob[i] < 0.20 || 
            proxy_results$r1_non_ccdc6_wilson_upper[i] < 0.55) {
          non_ccdc6_contra <- TRUE
        }
      }
      
      # Apply complete pass mode logic with stricter criteria
      if (ccdc6_both_ok && non_ccdc6_both_ok) {
        # Both subtypes meet all criteria
        proxy_results$proxy_status[i] <- "pass"
        proxy_results$proxy_pass[i] <- TRUE
        proxy_results$pass_mode[i] <- "both_ok"
        
      } else if (ccdc6_both_ok && non_ccdc6_can_eval_contra && !non_ccdc6_contra) {
        # CCDC6 OK, non-CCDC6 not contradicting (with n≥5)
        proxy_results$proxy_status[i] <- "pass"
        proxy_results$proxy_pass[i] <- TRUE
        proxy_results$pass_mode[i] <- "one_ok_not_contra"
        proxy_results$pass_subtype[i] <- "CCDC6"
        
      } else if (non_ccdc6_both_ok && ccdc6_can_eval_contra && !ccdc6_contra) {
        # non-CCDC6 OK, CCDC6 not contradicting (with n≥5)
        proxy_results$proxy_status[i] <- "pass"
        proxy_results$proxy_pass[i] <- TRUE
        proxy_results$pass_mode[i] <- "one_ok_not_contra"
        proxy_results$pass_subtype[i] <- "non-CCDC6"
        
      } else {
        # Neither condition met - fail
        proxy_results$proxy_status[i] <- "fail"
        proxy_results$proxy_pass[i] <- FALSE
        proxy_results$pass_mode[i] <- NA_character_
      }
    }
  }
  
  setTxtProgressBar(pb, i)
}
close(pb)

# --------------------------------------------------------------------------
# 3b.4.5 POC Trend Diagnosis (R1 samples only - diagnostic tags)
# --------------------------------------------------------------------------

cat("\n--- 3b.4.5 POC Trend Diagnosis (diagnostic tags only) ---\n")
cat("  Note: POC trend analysis is for interpretation assistance only.\n")
cat("  These tags DO NOT affect pair selection or judgment logic.\n\n")

# Extract POC values for R1 samples
# POC is at case level, need to map from samples to cases
r1_case_ids <- sample_metadata %>%
  filter(sample_submitter_id %in% r1_samples) %>%
  pull(case_submitter_id)

# Check for POC column name (could be 'poc' or 'POC')
if ("POC" %in% colnames(case_master)) {
  poc_col <- "POC"
} else if ("poc" %in% colnames(case_master)) {
  poc_col <- "poc"
} else {
  cat("  WARNING: POC column not found in case_master\n")
  poc_col <- NULL
}

if (!is.null(poc_col)) {
  r1_poc_values <- case_master %>%
    filter(case_submitter_id %in% r1_case_ids) %>%
    select(case_submitter_id, !!sym(poc_col)) %>%
    rename(poc = !!sym(poc_col))
  
  # Map POC values to R1 samples
  poc_r1 <- numeric(length(r1_samples))
  names(poc_r1) <- r1_samples
  
  for (samp in r1_samples) {
    case_id <- sample_metadata %>%
      filter(sample_submitter_id == samp) %>%
      pull(case_submitter_id)
    
    if (length(case_id) == 1 && case_id %in% r1_poc_values$case_submitter_id) {
      poc_r1[samp] <- r1_poc_values %>%
        filter(case_submitter_id == case_id) %>%
        pull(poc)
    } else {
      poc_r1[samp] <- NA
    }
  }
  
  cat(sprintf("  POC values available for %d/%d R1 samples\n",
              sum(!is.na(poc_r1)), length(r1_samples)))
  
  # Helper function for leave-one-out correlation stability
  check_loo_stability <- function(poc_valid, rev_valid) {
    if (length(poc_valid) < 10) return(NA)
    
    base_cor <- suppressWarnings(cor(poc_valid, rev_valid, method = "spearman", use = "complete.obs"))
    if (is.na(base_cor)) return(NA)
    
    base_sign <- sign(base_cor)
    sign_consistency <- 0
    
    for (i in seq_along(poc_valid)) {
      loo_cor <- suppressWarnings(cor(poc_valid[-i], rev_valid[-i], method = "spearman", use = "complete.obs"))
      if (!is.na(loo_cor) && sign(loo_cor) == base_sign) {
        sign_consistency <- sign_consistency + 1
      }
    }
    
    # Stable if sign consistent in at least (n_valid - 1) iterations
    return(sign_consistency >= (length(poc_valid) - 1))
  }
  
  # Process POC trend for each pair
  cat("  Computing POC-reversal correlations for each pair...\n")
  pb2 <- txtProgressBar(min = 0, max = nrow(passed_pairs), style = 3)
  
  for (i in 1:nrow(passed_pairs)) {
    # Only process if pair passed proxy check
    if (!proxy_results$proxy_pass[i]) {
      setTxtProgressBar(pb2, i)
      next
    }
    
    # Get R1 REO values for this pair
    r1_values_for_poc <- reo_r1_passed[i, ]
    
    # Get R0 direction from original results
    r0_direction_poc <- passed_pairs$r0_direction[i]
    
    # Determine reversal flags (1 = reversed, 0 = not reversed)
    # Only for samples outside dead zone
    dead_zone_mask <- abs(r1_values_for_poc) >= CONFIG$dead_zone
    reversal_flags <- numeric(length(r1_values_for_poc))
    
    if (r0_direction_poc == 1) {
      # R0 expects positive, reversal is negative
      reversal_flags[dead_zone_mask] <- as.numeric(r1_values_for_poc[dead_zone_mask] < 0)
    } else {
      # R0 expects negative, reversal is positive
      reversal_flags[dead_zone_mask] <- as.numeric(r1_values_for_poc[dead_zone_mask] > 0)
    }
    
    # Check if we have enough valid samples
    valid_for_poc <- !is.na(poc_r1) & dead_zone_mask
    n_valid_poc <- sum(valid_for_poc)
    
    if (n_valid_poc >= 10) {
      # Compute Spearman correlation
      poc_valid <- poc_r1[valid_for_poc]
      rev_valid <- reversal_flags[valid_for_poc]
      
      # Check for zero variance
      if (sd(rev_valid) == 0 || sd(poc_valid) == 0) {
        # Cannot compute correlation - all values are the same
        proxy_results$n_valid_poc[i] <- n_valid_poc
        proxy_results$trend_tag[i] <- "NA"
        proxy_results$trend_method[i] <- "zero_variance"
        proxy_results$rank_penalty_trend[i] <- 0
      } else {
        # Spearman correlation with exact test if possible
        cor_result <- suppressWarnings(
          cor.test(poc_valid, rev_valid, 
                   method = "spearman", 
                   alternative = "greater",  # One-sided test for positive correlation
                   exact = (n_valid_poc <= 30))  # Use exact test for small samples
        )
        
        proxy_results$rho_poc_flip[i] <- cor_result$estimate
        proxy_results$p_poc_flip[i] <- cor_result$p.value
        proxy_results$n_valid_poc[i] <- n_valid_poc
        proxy_results$trend_method[i] <- ifelse(n_valid_poc <= 30, "exact", "asymptotic")
        
        # Check leave-one-out stability
        proxy_results$stable[i] <- check_loo_stability(poc_valid, rev_valid)
        
        # Assign trend tags based on correlation and significance
        rho <- cor_result$estimate
        p_val <- cor_result$p.value
        
        if (abs(rho) >= 0.5 && p_val < 0.05) {
          if (rho > 0) {
            proxy_results$trend_tag[i] <- "POC_positive"
            proxy_results$rank_penalty_trend[i] <- 1  # Penalty for positive trend
          } else {
            proxy_results$trend_tag[i] <- "POC_negative"
            proxy_results$rank_penalty_trend[i] <- -1  # Bonus for negative trend
          }
        } else if (abs(rho) < 0.3) {
          proxy_results$trend_tag[i] <- "POC_neutral"
          proxy_results$rank_penalty_trend[i] <- 0
        } else {
          proxy_results$trend_tag[i] <- "POC_weak"
          proxy_results$rank_penalty_trend[i] <- 0
        }
      }
      
      # Check for partner heterogeneity (supplementary)
      if (length(unique(r1_sample_info$fusion_partner)) > 1) {
        # Kruskal-Wallis test for differences between partners
        partner_groups <- r1_sample_info %>%
          filter(sample_submitter_id %in% names(which(valid_for_poc))) %>%
          pull(fusion_partner)
        
        if (length(unique(partner_groups)) > 1 && length(partner_groups) == sum(valid_for_poc)) {
          kw_test <- suppressWarnings(kruskal.test(reversal_flags[valid_for_poc] ~ partner_groups))
          proxy_results$partner_hetero_p[i] <- kw_test$p.value
          proxy_results$partner_k[i] <- length(unique(partner_groups))
          
          # Count per partner group
          per_level_counts <- table(partner_groups)
          proxy_results$per_level_n[i] <- paste(per_level_counts, collapse = ",")
          
          if (kw_test$p.value < 0.05) {
            proxy_results$rank_penalty_partner[i] <- 1  # Penalty for heterogeneity
          }
        }
      }
    } else if (n_valid_poc > 0) {
      # Not enough samples for correlation
      proxy_results$n_valid_poc[i] <- n_valid_poc
      proxy_results$trend_tag[i] <- "insufficient_n"
      proxy_results$trend_method[i] <- "none"
    }
    
    setTxtProgressBar(pb2, i)
  }
  close(pb2)
  
  # Summary of POC trends
  trend_table <- table(proxy_results$trend_tag[proxy_results$proxy_pass])
  cat("\n  POC trend distribution (passed pairs only):\n")
  for (tag in names(trend_table)) {
    cat(sprintf("    %s: %d\n", tag, trend_table[tag]))
  }
  
  # Count stable correlations
  n_stable <- sum(proxy_results$stable == TRUE & proxy_results$proxy_pass, na.rm = TRUE)
  n_evaluated <- sum(!is.na(proxy_results$stable) & proxy_results$proxy_pass)
  cat(sprintf("  Stable correlations (LOO): %d/%d (%.1f%%)\n",
              n_stable, n_evaluated, 
              ifelse(n_evaluated > 0, n_stable/n_evaluated * 100, 0)))
} else {
  cat("  POC analysis skipped (POC column not found)\n")
  trend_table <- NULL
}

# --------------------------------------------------------------------------
# 3b.5 Summary statistics
# --------------------------------------------------------------------------

cat("\n--- 3b.5 Summary of proxy check results ---\n")

# Pass mode summary
pass_mode_table <- table(proxy_results$pass_mode[proxy_results$proxy_pass])
cat("\nPass modes for successful pairs:\n")
for (mode in names(pass_mode_table)) {
  cat(sprintf("  %s: %d\n", mode, pass_mode_table[mode]))
}

# Detailed analysis of one_ok_not_contra cases
if ("one_ok_not_contra" %in% names(pass_mode_table) && pass_mode_table["one_ok_not_contra"] > 0) {
  cat("\n--- Detailed analysis of one_ok_not_contra cases ---\n")
  
  # Which subtype was OK?
  one_ok_cases <- proxy_results[proxy_results$pass_mode == "one_ok_not_contra" & !is.na(proxy_results$pass_mode), ]
  subtype_counts <- table(one_ok_cases$pass_subtype)
  
  cat("Which subtype passed:\n")
  if ("CCDC6" %in% names(subtype_counts)) {
    cat(sprintf("  CCDC6 OK: %d cases\n", subtype_counts["CCDC6"]))
  } else {
    cat("  CCDC6 OK: 0 cases\n")
  }
  if ("non-CCDC6" %in% names(subtype_counts)) {
    cat(sprintf("  non-CCDC6 OK: %d cases\n", subtype_counts["non-CCDC6"]))
  } else {
    cat("  non-CCDC6 OK: 0 cases\n")
  }
  
  # CI upper bounds of the non-OK side
  cat("\nCI upper bounds of non-OK side:\n")
  
  # For CCDC6 OK cases, get non-CCDC6 CI upper
  ccdc6_ok_indices <- which(proxy_results$pass_mode == "one_ok_not_contra" & 
                              proxy_results$pass_subtype == "CCDC6" & 
                              !is.na(proxy_results$pass_subtype))
  if (length(ccdc6_ok_indices) > 0) {
    non_ccdc6_uppers <- proxy_results$r1_non_ccdc6_wilson_upper[ccdc6_ok_indices]
    non_ccdc6_uppers <- non_ccdc6_uppers[!is.na(non_ccdc6_uppers)]
    if (length(non_ccdc6_uppers) > 0) {
      cat(sprintf("  When CCDC6 OK, non-CCDC6 CI upper: median=%.3f, range=[%.3f, %.3f]\n",
                  median(non_ccdc6_uppers),
                  min(non_ccdc6_uppers),
                  max(non_ccdc6_uppers)))
    }
  }
  
  # For non-CCDC6 OK cases, get CCDC6 CI upper
  non_ccdc6_ok_indices <- which(proxy_results$pass_mode == "one_ok_not_contra" & 
                                  proxy_results$pass_subtype == "non-CCDC6" & 
                                  !is.na(proxy_results$pass_subtype))
  if (length(non_ccdc6_ok_indices) > 0) {
    ccdc6_uppers <- proxy_results$r1_ccdc6_wilson_upper[non_ccdc6_ok_indices]
    ccdc6_uppers <- ccdc6_uppers[!is.na(ccdc6_uppers)]
    if (length(ccdc6_uppers) > 0) {
      cat(sprintf("  When non-CCDC6 OK, CCDC6 CI upper: median=%.3f, range=[%.3f, %.3f]\n",
                  median(ccdc6_uppers),
                  min(ccdc6_uppers),
                  max(ccdc6_uppers)))
    }
  }
}

# Overall pass rate
cat(sprintf("\nOverall proxy check pass rate: %d/%d (%.1f%%)\n",
            sum(proxy_results$proxy_pass),
            nrow(proxy_results),
            sum(proxy_results$proxy_pass) / nrow(proxy_results) * 100))

# Extreme cases
if (sum(proxy_results$extreme_case_flag) > 0) {
  cat(sprintf("  Extreme cases detected: %d pairs\n", 
              sum(proxy_results$extreme_case_flag)))
}

# Display Bayesian vs Wilson comparison for debugging
cat("\n--- Bayesian vs Wilson CI comparison (first 10 passed pairs) ---\n")
comparison_data <- proxy_results %>%
  filter(proxy_pass) %>%
  select(pair_id, 
         r1_ccdc6_n, r1_ccdc6_reversal_rate, 
         r1_ccdc6_bayes_prob, r1_ccdc6_wilson_lower,
         r1_non_ccdc6_n, r1_non_ccdc6_reversal_rate,
         r1_non_ccdc6_bayes_prob, r1_non_ccdc6_wilson_lower) %>%
  head(10)

if (nrow(comparison_data) > 0) {
  cat("CCDC6 subtype:\n")
  for (i in 1:min(5, nrow(comparison_data))) {
    if (!is.na(comparison_data$r1_ccdc6_n[i])) {
      cat(sprintf("  Pair %s: n=%d, rate=%.2f, Pr(p>0.5)=%.3f, Wilson_lower=%.3f\n",
                  comparison_data$pair_id[i],
                  comparison_data$r1_ccdc6_n[i],
                  comparison_data$r1_ccdc6_reversal_rate[i],
                  comparison_data$r1_ccdc6_bayes_prob[i],
                  comparison_data$r1_ccdc6_wilson_lower[i]))
    }
  }
  
  cat("\nnon-CCDC6 subtype:\n")
  for (i in 1:min(5, nrow(comparison_data))) {
    if (!is.na(comparison_data$r1_non_ccdc6_n[i])) {
      cat(sprintf("  Pair %s: n=%d, rate=%.2f, Pr(p>0.5)=%.3f, Wilson_lower=%.3f\n",
                  comparison_data$pair_id[i],
                  comparison_data$r1_non_ccdc6_n[i],
                  comparison_data$r1_non_ccdc6_reversal_rate[i],
                  comparison_data$r1_non_ccdc6_bayes_prob[i],
                  comparison_data$r1_non_ccdc6_wilson_lower[i]))
    }
  }
}

# Filter to passed pairs
proxy_passed_pairs <- passed_pairs[proxy_results$proxy_pass, ]
cat(sprintf("\n✔ Pairs passing Step 3b: %d\n", nrow(proxy_passed_pairs)))

# --------------------------------------------------------------------------
# 3b.6 Save Step 3b results
# --------------------------------------------------------------------------

cat("\n--- 3b.6 Saving Step 3b results ---\n")

# Combine with Step 3a results
step3b_results <- cbind(
  passed_pairs,
  proxy_results[, -1]  # Remove duplicate pair_id
)

# Add columns to clearly identify which subtype was OK in one_ok_not_contra cases
# and provide the decisive metrics
step3b_results$ok_subtype_detail <- NA_character_
step3b_results$ok_subtype_metrics <- NA_character_

for (i in 1:nrow(step3b_results)) {
  if (step3b_results$pass_mode[i] == "one_ok_not_contra" && !is.na(step3b_results$pass_mode[i])) {
    if (step3b_results$pass_subtype[i] == "CCDC6" && !is.na(step3b_results$pass_subtype[i])) {
      # CCDC6 was OK
      step3b_results$ok_subtype_detail[i] <- "CCDC6"
      step3b_results$ok_subtype_metrics[i] <- sprintf(
        "CCDC6: n=%d, rate=%.3f, Pr=%.3f, CI=[%.3f,%.3f] | non-CCDC6: n=%d, CI_upper=%.3f",
        step3b_results$r1_ccdc6_n[i],
        step3b_results$r1_ccdc6_reversal_rate[i],
        step3b_results$r1_ccdc6_bayes_prob[i],
        step3b_results$r1_ccdc6_wilson_lower[i],
        step3b_results$r1_ccdc6_wilson_upper[i],
        ifelse(is.na(step3b_results$r1_non_ccdc6_n[i]), 0, step3b_results$r1_non_ccdc6_n[i]),
        ifelse(is.na(step3b_results$r1_non_ccdc6_wilson_upper[i]), NA, step3b_results$r1_non_ccdc6_wilson_upper[i])
      )
    } else if (step3b_results$pass_subtype[i] == "non-CCDC6" && !is.na(step3b_results$pass_subtype[i])) {
      # non-CCDC6 was OK
      step3b_results$ok_subtype_detail[i] <- "non-CCDC6"
      step3b_results$ok_subtype_metrics[i] <- sprintf(
        "non-CCDC6: n=%d, rate=%.3f, Pr=%.3f, CI=[%.3f,%.3f] | CCDC6: n=%d, CI_upper=%.3f",
        step3b_results$r1_non_ccdc6_n[i],
        step3b_results$r1_non_ccdc6_reversal_rate[i],
        step3b_results$r1_non_ccdc6_bayes_prob[i],
        step3b_results$r1_non_ccdc6_wilson_lower[i],
        step3b_results$r1_non_ccdc6_wilson_upper[i],
        ifelse(is.na(step3b_results$r1_ccdc6_n[i]), 0, step3b_results$r1_ccdc6_n[i]),
        ifelse(is.na(step3b_results$r1_ccdc6_wilson_upper[i]), NA, step3b_results$r1_ccdc6_wilson_upper[i])
      )
    }
  } else if (step3b_results$pass_mode[i] == "both_ok" && !is.na(step3b_results$pass_mode[i])) {
    step3b_results$ok_subtype_detail[i] <- "both"
    step3b_results$ok_subtype_metrics[i] <- sprintf(
      "CCDC6: n=%d, rate=%.3f, Pr=%.3f | non-CCDC6: n=%d, rate=%.3f, Pr=%.3f",
      ifelse(is.na(step3b_results$r1_ccdc6_n[i]), 0, step3b_results$r1_ccdc6_n[i]),
      ifelse(is.na(step3b_results$r1_ccdc6_reversal_rate[i]), 0, step3b_results$r1_ccdc6_reversal_rate[i]),
      ifelse(is.na(step3b_results$r1_ccdc6_bayes_prob[i]), 0, step3b_results$r1_ccdc6_bayes_prob[i]),
      ifelse(is.na(step3b_results$r1_non_ccdc6_n[i]), 0, step3b_results$r1_non_ccdc6_n[i]),
      ifelse(is.na(step3b_results$r1_non_ccdc6_reversal_rate[i]), 0, step3b_results$r1_non_ccdc6_reversal_rate[i]),
      ifelse(is.na(step3b_results$r1_non_ccdc6_bayes_prob[i]), 0, step3b_results$r1_non_ccdc6_bayes_prob[i])
    )
  }
}

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
  pass_rate = sum(proxy_results$proxy_pass) / nrow(proxy_results),
  
  # Pass mode summary
  pass_mode_summary = table(proxy_results$pass_mode[proxy_results$proxy_pass]),
  
  # POC trend diagnosis summary (diagnostic only)
  poc_trend_summary = list(
    trend_distribution = table(proxy_results$trend_tag[proxy_results$proxy_pass]),
    n_stable = sum(proxy_results$stable == TRUE & proxy_results$proxy_pass, na.rm = TRUE),
    n_evaluated = sum(!is.na(proxy_results$stable) & proxy_results$proxy_pass),
    n_hetero = sum(proxy_results$partner_hetero_p < 0.05 & proxy_results$proxy_pass, na.rm = TRUE),
    diagnostic_note = "POC trend tags are for interpretation assistance only and do not affect pair selection"
  ),
  
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

# Export decision rationale summary for passed pairs
# Ensure pass_subtype exists (for compatibility)
if (!"pass_subtype" %in% names(step3b_results)) {
  step3b_results$pass_subtype <- NA_character_
}

decision_summary <- step3b_results[step3b_results$proxy_pass, 
                                   c("pair_id", "gene_up", "gene_down", 
                                     "proxy_status", "pass_mode",
                                     "ok_subtype_detail", "ok_subtype_metrics",
                                     "r1_ccdc6_n", "r1_ccdc6_reversal_rate", 
                                     "r1_ccdc6_bayes_prob", "r1_ccdc6_wilson_lower",
                                     "r1_non_ccdc6_n", "r1_non_ccdc6_reversal_rate",
                                     "r1_non_ccdc6_bayes_prob", "r1_non_ccdc6_wilson_lower")]
write.csv(decision_summary,
          paste0(paths$output, "step3b_decision_rationale.csv"),
          row.names = FALSE)
cat("  Decision rationale exported to: step3b_decision_rationale.csv\n")

# Quick summary of one_ok_not_contra decisions
if (sum(step3b_results$pass_mode == "one_ok_not_contra", na.rm = TRUE) > 0) {
  cat("\n  One-OK decision breakdown in CSV:\n")
  one_ok_ccdc6 <- sum(step3b_results$ok_subtype_detail == "CCDC6", na.rm = TRUE)
  one_ok_non <- sum(step3b_results$ok_subtype_detail == "non-CCDC6", na.rm = TRUE)
  cat(sprintf("    CCDC6 was decisive: %d pairs\n", one_ok_ccdc6))
  cat(sprintf("    non-CCDC6 was decisive: %d pairs\n", one_ok_non))
}

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
cat(sprintf("Veto contradiction: %d\n",
            sum(proxy_results$proxy_status == "veto_contradiction", na.rm = TRUE)))
cat(sprintf("Extreme case fail: %d\n",
            sum(proxy_results$extreme_case_flag)))

# Display pass mode breakdown
cat("\nPass modes:\n")
if (length(pass_mode_table) > 0) {
  for (mode in names(pass_mode_table)) {
    cat(sprintf("  %s: %d\n", mode, pass_mode_table[mode]))
  }
}

# One-line summary of one_ok_not_contra
if ("one_ok_not_contra" %in% names(pass_mode_table) && pass_mode_table["one_ok_not_contra"] > 0) {
  one_ok_summary <- proxy_results[proxy_results$pass_mode == "one_ok_not_contra" & !is.na(proxy_results$pass_mode), ]
  ccdc6_ok_count <- sum(one_ok_summary$pass_subtype == "CCDC6", na.rm = TRUE)
  non_ccdc6_ok_count <- sum(one_ok_summary$pass_subtype == "non-CCDC6", na.rm = TRUE)
  cat(sprintf("  → one_ok breakdown: CCDC6 OK=%d, non-CCDC6 OK=%d\n", 
              ccdc6_ok_count, non_ccdc6_ok_count))
}

# Display Bayesian criterion summary
cat("\nBayesian criterion (Jeffreys prior):\n")
cat(sprintf("  R1 OK criterion: Pr(p > 0.5) ≥ 0.95 AND n≥6 AND rate≥0.7\n"))
cat(sprintf("  Non-contradiction: Pr(p > 0.5) ≥ 0.20 AND Wilson upper ≥ 0.55 (n≥5 required)\n"))
cat(sprintf("  Veto: one-sided LCL≥0.7 AND UCL_other≤0.5 (both n≥5)\n"))

# POC trend diagnosis summary (diagnostic only)
cat("\nPOC trend diagnosis (R1, diagnostic only):\n")
if (exists("trend_table") && length(trend_table) > 0) {
  for (tag in names(trend_table)) {
    cat(sprintf("  %s: %d\n", tag, trend_table[tag]))
  }
}
cat(sprintf("  Method: Spearman correlation (one-sided, exact when n≤30)\n"))
cat(sprintf("  Min n_valid: 10, effect size threshold: |ρ|≥0.5\n"))
cat(sprintf("  Note: These tags DO NOT affect selection/judgment\n"))

cat(sprintf("\nTimestamp: %s\n", format(step3b_data$timestamp)))

# Final message
if (step3b_data$n_passed > 0) {
  cat(sprintf("\n✔ Ready for Step 4: Redundancy Removal (%d pairs available)\n", 
              step3b_data$n_passed))
} else {
  cat("\n⚠ No pairs passed Step 3b. Review parameters or data quality.\n")
}