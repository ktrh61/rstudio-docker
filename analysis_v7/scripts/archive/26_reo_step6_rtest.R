# 26_reo_step6_rtest.R - REO Step 6: R_test Validation and Robustness
# Purpose: Validate REO panel on R_test samples with comprehensive diagnostics
# Input: reo_step5_data.rds, thyr_case_master_full.rds, thyr_se_strand2_nonzero.rds
# Output: reo_step6_validation.rds, step6_validation_summary.csv
# Version: v2.0
# Date: 2025-01-27
# Note: This panel assesses the presence of RET fusion REO signature. It does NOT directly estimate radiation exposure.
#
# Classification Categories:
# - Positive: n_reversals >= T+1 (clear positive signal)
# - Negative: n_reversals <= T-2 (clear negative signal)  
# - Boundary: n_reversals == T-1 or T (ambiguous zone)
# - Undetermined: n_eff < ceil(0.7×N) (insufficient data)

source("analysis_v7/setup.R")

cat("\n=== REO STEP 6: R_test Validation and Robustness ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("\n*** This panel assesses RET fusion REO signature presence. ***\n")
cat("*** It does NOT directly estimate radiation exposure. ***\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
  library(matrixStats)
})

# Load utility functions for contamDE
source("utils/utils_improved.R")
source("utils/norm_improved.R")
source("utils/contamde_purity_functions.R")

# --------------------------------------------------------------------------
# 6.1 Configuration
# --------------------------------------------------------------------------

cat("--- 6.1 Configuration ---\n")

VALIDATION_CONFIG <- list(
  # Purity thresholds
  purity_threshold = 0.6,        # High purity threshold (60%)
  
  # POC stratification (2-layer analysis)
  poc_cutoff = 33.3,             # POC stratification point
  
  # N_eff threshold
  n_eff_threshold = 0.7,          # Proportion of N required for valid judgment
  
  # Minimum sample sizes
  min_samples_per_stratum = 5,   # Minimum for CI calculation
  
  # Jackknife settings
  jackknife_enabled = TRUE,      # Perform N-1 jackknife
  
  # BRAF validation (internal QC only)
  braf_validation = TRUE,        # Test on BRAF samples
  
  # ContamDE settings
  pairwise_method = "lts",       # MUREN method
  workers = 3,                   # Parallel workers
  prior_count = 1.0,             # Pseudocount
  
  # CI method
  ci_method = "wilson",          # Wilson CI (no continuity correction)
  ci_level = 0.95                # 95% confidence interval
)

# Display configuration
for (param in names(VALIDATION_CONFIG)) {
  cat(sprintf("  %s: %s\n", param, VALIDATION_CONFIG[[param]]))
}

# --------------------------------------------------------------------------
# Helper function: REO classification (B-scheme)
# --------------------------------------------------------------------------

classify_reo <- function(n_reversals, n_eff, N, T, n_eff_threshold = 0.7) {
  # Step 1: Check if judgment is possible
  if (n_eff < ceiling(n_eff_threshold * N)) {
    return(list(
      RET_REO_call = "Undetermined",
      judged = FALSE
    ))
  }
  
  # Step 2: Apply B-scheme classification
  if (n_reversals >= T + 1) {
    call <- "Positive"
  } else if (n_reversals == T || n_reversals == T - 1) {
    call <- "Boundary"
  } else {  # n_reversals <= T - 2
    call <- "Negative"
  }
  
  return(list(
    RET_REO_call = call,
    judged = TRUE
  ))
}

# --------------------------------------------------------------------------
# 6.2 Load Step 5 results and data
# --------------------------------------------------------------------------

cat("\n--- 6.2 Loading data ---\n")

# Load Step 5 results
step5_data_path <- paste0(paths$processed, "reo_step5_data.rds")
if (!file.exists(step5_data_path)) {
  stop("ERROR: reo_step5_data.rds not found. Please run 25_reo_step5_panelize.R first.")
}

step5_data <- readRDS(step5_data_path)
cat("  Step 5 data loaded successfully\n")

# Extract panel information
panel_info <- step5_data$panel_info
panel_pairs <- step5_data$panel_pairs
N <- step5_data$N
T <- step5_data$T
CONFIG <- step5_data$config  # Original CONFIG
qc_anchors <- step5_data$qc_anchors

cat(sprintf("  Panel: N=%d pairs, threshold T=%d\n", N, T))
cat(sprintf("  Boundary rule: n_reversals in {%d, %d}\n", T-1, T))
cat(sprintf("  N_eff threshold: n_eff >= ceil(%.1f×%d) = %d\n", 
            VALIDATION_CONFIG$n_eff_threshold, N, 
            ceiling(VALIDATION_CONFIG$n_eff_threshold * N)))

# Load case master
case_master_path <- paste0(paths$processed, "thyr_case_master_full.rds")
if (!file.exists(case_master_path)) {
  stop("ERROR: thyr_case_master_full.rds not found. Please run 04_clinical_integration.R first.")
}

case_master <- readRDS(case_master_path)
cat("  Case master loaded successfully\n")

# Load SummarizedExperiment
se_path <- paste0(paths$processed, "thyr_se_strand2_nonzero.rds")
if (!file.exists(se_path)) {
  stop("ERROR: thyr_se_strand2_nonzero.rds not found. Please run 03_prepare_counts.R first.")
}

se <- readRDS(se_path)
cat("  SummarizedExperiment loaded successfully\n")

# --------------------------------------------------------------------------
# 6.3 Extract R_test cases and samples
# --------------------------------------------------------------------------

cat("\n--- 6.3 Extracting R_test cases ---\n")

# Get R_test cases
r_test_cases <- case_master[case_master$group == "R_test", ]
cat(sprintf("  R_test cases found: %d\n", nrow(r_test_cases)))

# Function to extract tumor samples
extract_r_test_samples <- function(case_ids, metadata, sample_type = "Primary Tumor") {
  samples <- c()
  cases <- c()
  
  for (case_id in case_ids) {
    case_samples <- metadata$sample_submitter_id[
      metadata$case_submitter_id == case_id & 
        metadata$sample_type == sample_type
    ]
    
    if (length(case_samples) > 0) {
      # Take first tumor sample if multiple exist
      samples <- c(samples, case_samples[1])
      cases <- c(cases, case_id)
    }
  }
  
  return(list(samples = samples, cases = cases))
}

# Extract R_test tumor samples
sample_metadata <- colData(se)
r_test_all <- extract_r_test_samples(
  r_test_cases$case_submitter_id,
  sample_metadata,
  "Primary Tumor"
)

cat(sprintf("  R_test tumor samples extracted: %d\n", length(r_test_all$samples)))

# Identify paired and unpaired samples
r_test_paired <- data.frame()
r_test_unpaired <- data.frame()

for (i in seq_along(r_test_all$cases)) {
  case_id <- r_test_all$cases[i]
  
  # Check for normal sample
  normal_samples <- sample_metadata$sample_submitter_id[
    sample_metadata$case_submitter_id == case_id & 
      sample_metadata$sample_type == "Solid Tissue Normal"
  ]
  
  if (length(normal_samples) > 0) {
    r_test_paired <- rbind(r_test_paired, data.frame(
      case_id = case_id,
      tumor_sample = r_test_all$samples[i],
      normal_sample = normal_samples[1],
      stringsAsFactors = FALSE
    ))
  } else {
    r_test_unpaired <- rbind(r_test_unpaired, data.frame(
      case_id = case_id,
      tumor_sample = r_test_all$samples[i],
      stringsAsFactors = FALSE
    ))
  }
}

cat(sprintf("\n  Sample breakdown:\n"))
cat(sprintf("    Paired (T/N): %d cases\n", nrow(r_test_paired)))
cat(sprintf("    Unpaired (T only): %d cases\n", nrow(r_test_unpaired)))

# --------------------------------------------------------------------------
# 6.4 Estimate tumor purity (ONE TIME using contamDE)
# --------------------------------------------------------------------------

cat("\n--- 6.4 Estimating tumor purity (contamDE) ---\n")
cat("  Note: Purity estimation performed ONCE for all R_test samples\n")

if (nrow(r_test_paired) > 0) {
  # Prepare paired samples for contamDE
  paired_samples <- list(
    cases = r_test_paired$case_id,
    tumors = r_test_paired$tumor_sample,
    normals = r_test_paired$normal_sample
  )
  
  cat(sprintf("  Using %d paired samples for purity estimation\n", 
              nrow(r_test_paired)))
  
  # Extract counts
  counts_for_purity <- cbind(
    assay(se, "counts")[, paired_samples$normals, drop = FALSE],
    assay(se, "counts")[, paired_samples$tumors, drop = FALSE]
  )
  
  # Filter genes
  keep_genes <- rowSums(counts_for_purity > 0) >= ncol(counts_for_purity) * 0.5
  
  n_pairs <- length(paired_samples$cases)
  if (sum(keep_genes) >= 1000 && n_pairs >= 3) {
    cat(sprintf("  Genes retained for contamDE: %d\n", sum(keep_genes)))
    
    # Prepare contamDE matrix
    counts_contamde <- counts_for_purity[keep_genes, , drop = FALSE]
    
    # Set informative column names
    colnames(counts_contamde) <- c(
      paste0("Normal_", seq_len(n_pairs)),
      paste0("Tumor_", seq_len(n_pairs))
    )
    
    # Run contamDE
    cat("  Running contamDE purity estimation...\n")
    
    purity_result <- tryCatch({
      contamde_purity(
        counts = counts_contamde,
        subtype = NULL,
        covariate = NULL,
        contaminated = TRUE,
        pairwise_method = VALIDATION_CONFIG$pairwise_method,
        workers = VALIDATION_CONFIG$workers,
        prior.count = VALIDATION_CONFIG$prior_count,
        verbose = TRUE
      )
    }, error = function(e) {
      cat("  WARNING: ContamDE failed:", conditionMessage(e), "\n")
      NULL
    })
    
    if (!is.null(purity_result)) {
      # Assess purity quality
      quality_stats <- assess_purity_quality(
        purity_result, 
        threshold = VALIDATION_CONFIG$purity_threshold,
        verbose = TRUE
      )
      
      # Create purity mapping
      purity_values <- purity_result$proportion
      names(purity_values) <- paired_samples$cases
      
    } else {
      purity_values <- NULL
      quality_stats <- NULL
    }
  } else {
    cat("  Insufficient data for contamDE\n")
    purity_values <- NULL
    quality_stats <- NULL
  }
} else {
  cat("  No paired samples in R_test\n")
  purity_values <- NULL
  quality_stats <- NULL
}

# --------------------------------------------------------------------------
# 6.5 Prepare sample information
# --------------------------------------------------------------------------

cat("\n--- 6.5 Preparing sample information ---\n")

# Create sample classification dataframe
r_test_sample_info <- data.frame(
  case_id = r_test_all$cases,
  sample_id = r_test_all$samples,
  stringsAsFactors = FALSE
)

# Add POC values (continuous, for stratification)
r_test_sample_info$poc <- r_test_cases$POC[
  match(r_test_sample_info$case_id, r_test_cases$case_submitter_id)
]

# Add purity values and classification
r_test_sample_info$purity <- NA_real_
r_test_sample_info$purity_class <- "Unknown"

if (!is.null(purity_values)) {
  # Match purity values
  idx <- match(r_test_sample_info$case_id, names(purity_values))
  r_test_sample_info$purity[!is.na(idx)] <- purity_values[idx[!is.na(idx)]]
  
  # Classify purity
  r_test_sample_info$purity_class[!is.na(r_test_sample_info$purity) & 
                                    r_test_sample_info$purity >= VALIDATION_CONFIG$purity_threshold] <- "High"
  r_test_sample_info$purity_class[!is.na(r_test_sample_info$purity) & 
                                    r_test_sample_info$purity < VALIDATION_CONFIG$purity_threshold] <- "Low"
}

# POC stratification (2-layer for main analysis)
r_test_sample_info$poc_stratum <- ifelse(
  r_test_sample_info$poc <= VALIDATION_CONFIG$poc_cutoff,
  "Low_POC",
  "High_POC"
)

# Summary
cat("\nPurity classification:\n")
print(table(r_test_sample_info$purity_class))

cat("\nPOC stratification (2-layer):\n")
print(table(r_test_sample_info$poc_stratum))

cat("\nCross-tabulation (Purity × POC):\n")
print(table(r_test_sample_info$purity_class, r_test_sample_info$poc_stratum))

# --------------------------------------------------------------------------
# 6.6 Calculate REO scores for all R_test samples
# --------------------------------------------------------------------------

cat("\n--- 6.6 Calculating REO scores for R_test samples ---\n")
cat("  Note: REO-score = n_reversals / N (denominator is fixed at N)\n")

# Get counts and calculate TPM
r_test_counts <- assay(se, "counts")[, r_test_sample_info$sample_id, drop = FALSE]

# Simple TPM calculation (as in Step 5)
rpk <- r_test_counts / 1000  # Assume 1kb length for all
scaling_factor <- colSums(rpk) / 1e6
tpm_r_test <- sweep(rpk, 2, scaling_factor, "/")
log2_tpm_r_test <- log2(tpm_r_test + 1)

# Calculate REO values for each sample
calculate_sample_reo <- function(sample_idx) {
  sample_tpm <- log2_tpm_r_test[, sample_idx]
  
  # Calculate REO for each pair
  reo_values <- numeric(N)
  for (i in 1:N) {
    gene_up <- panel_pairs$gene_up[i]
    gene_down <- panel_pairs$gene_down[i]
    
    if (gene_up %in% rownames(log2_tpm_r_test) && 
        gene_down %in% rownames(log2_tpm_r_test)) {
      reo_values[i] <- sample_tpm[gene_up] - sample_tpm[gene_down]
    } else {
      reo_values[i] <- NA_real_
    }
  }
  
  return(reo_values)
}

# Calculate for all samples
reo_matrix_r_test <- matrix(NA_real_, nrow = N, ncol = nrow(r_test_sample_info))
colnames(reo_matrix_r_test) <- r_test_sample_info$sample_id

for (i in seq_len(nrow(r_test_sample_info))) {
  reo_matrix_r_test[, i] <- calculate_sample_reo(i)
}

# Apply classification using B-scheme
r_test_classifications <- data.frame(
  sample_id = r_test_sample_info$sample_id,
  case_id = r_test_sample_info$case_id,
  RET_REO_call = character(nrow(r_test_sample_info)),
  reo_score = numeric(nrow(r_test_sample_info)),
  n_reversals = integer(nrow(r_test_sample_info)),
  n_eff = integer(nrow(r_test_sample_info)),
  judged = logical(nrow(r_test_sample_info)),
  stringsAsFactors = FALSE
)

# Classify each sample
for (i in seq_len(nrow(r_test_sample_info))) {
  sample_reo <- reo_matrix_r_test[, i]
  
  # Calculate n_eff and n_reversals
  valid_mask <- abs(sample_reo) >= CONFIG$dead_zone
  n_eff <- sum(valid_mask, na.rm = TRUE)
  
  if (n_eff >= ceiling(VALIDATION_CONFIG$n_eff_threshold * N)) {
    valid_values <- sample_reo[valid_mask]
    valid_directions <- sign(valid_values)
    valid_r0_directions <- panel_info$r0_directions[valid_mask]
    
    n_reversals <- sum(valid_directions != valid_r0_directions, na.rm = TRUE)
  } else {
    n_reversals <- NA_integer_
  }
  
  # Calculate REO-score (always n_reversals/N)
  reo_score <- ifelse(!is.na(n_reversals), n_reversals / N, NA_real_)
  
  # Apply classification
  classification_result <- classify_reo(
    n_reversals, n_eff, N, T, 
    n_eff_threshold = VALIDATION_CONFIG$n_eff_threshold
  )
  
  r_test_classifications[i, c("RET_REO_call", "reo_score", "n_reversals", 
                              "n_eff", "judged")] <- 
    list(classification_result$RET_REO_call, reo_score, n_reversals, 
         n_eff, classification_result$judged)
}

# Merge with sample info
r_test_results <- cbind(r_test_sample_info, r_test_classifications[, -c(1:2)])

cat(sprintf("  Classifications complete: %d samples\n", nrow(r_test_results)))
cat("\n  Overall classification summary:\n")
cat(sprintf("    Judged: %d, Undetermined: %d\n", 
            sum(r_test_results$judged), 
            sum(!r_test_results$judged)))

# Show breakdown for judged samples
judged_results <- r_test_results[r_test_results$judged, ]
if (nrow(judged_results) > 0) {
  cat("\n  Among judged samples:\n")
  classification_table <- table(judged_results$RET_REO_call)
  for (call in c("Positive", "Boundary", "Negative")) {
    count <- ifelse(call %in% names(classification_table), 
                    classification_table[call], 0)
    cat(sprintf("    %s: %d\n", call, count))
  }
}

# --------------------------------------------------------------------------
# Part A: Main Analysis - Stratified by POC (High purity only)
# --------------------------------------------------------------------------

cat("\n=== PART A: Main R_test Analysis (High Purity Only) ===\n")

# Filter to high purity
high_purity_results <- r_test_results[r_test_results$purity_class == "High", ]

cat(sprintf("\nHigh purity samples: %d/%d (%.1f%%)\n",
            nrow(high_purity_results), nrow(r_test_results),
            nrow(high_purity_results)/nrow(r_test_results)*100))

if (nrow(high_purity_results) > 0) {
  # Stratify by POC (2-layer)
  strata <- list(
    Low_POC = high_purity_results[high_purity_results$poc_stratum == "Low_POC", ],
    High_POC = high_purity_results[high_purity_results$poc_stratum == "High_POC", ]
  )
  
  # Calculate performance by stratum
  stratum_results <- list()
  
  for (stratum_name in names(strata)) {
    stratum_data <- strata[[stratum_name]]
    stratum_judged <- stratum_data[stratum_data$judged, ]
    
    cat(sprintf("\n%s stratum (≤33.3%% vs >33.3%%):\n", stratum_name))
    cat(sprintf("  Total: %d\n", nrow(stratum_data)))
    cat(sprintf("  Judged: %d\n", nrow(stratum_judged)))
    cat(sprintf("  Excluded (Undetermined): %d\n", 
                sum(!stratum_data$judged)))
    
    if (nrow(stratum_judged) > 0) {
      # Count each category
      n_positive <- sum(stratum_judged$RET_REO_call == "Positive")
      n_negative <- sum(stratum_judged$RET_REO_call == "Negative")
      n_boundary <- sum(stratum_judged$RET_REO_call == "Boundary")
      
      # Display breakdown
      cat(sprintf("  Positive/Boundary/Negative: %d/%d/%d\n", 
                  n_positive, n_boundary, n_negative))
      
      # Calculate positive rate (excluding Boundary from denominator)
      n_denominator <- n_positive + n_negative
      
      if (n_denominator >= VALIDATION_CONFIG$min_samples_per_stratum) {
        pos_rate <- n_positive / n_denominator
        
        # Wilson CI (no continuity correction)
        wilson_ci <- prop.test(n_positive, n_denominator, 
                               conf.level = VALIDATION_CONFIG$ci_level,
                               correct = FALSE)$conf.int
        
        cat(sprintf("  Positive rate: %d/%d = %.1f%% [%.1f%%, %.1f%%]\n",
                    n_positive, n_denominator, pos_rate*100,
                    wilson_ci[1]*100, wilson_ci[2]*100))
        cat(sprintf("  (Note: Denominator = Positive + Negative, excludes Boundary)\n"))
        
        stratum_results[[stratum_name]] <- list(
          n_total = nrow(stratum_data),
          n_judged = nrow(stratum_judged),
          n_positive = n_positive,
          n_negative = n_negative,
          n_boundary = n_boundary,
          n_undetermined = sum(!stratum_data$judged),
          n_denominator = n_denominator,
          positive_rate = pos_rate,
          ci_lower = wilson_ci[1],
          ci_upper = wilson_ci[2],
          sufficient_samples = TRUE
        )
      } else {
        cat(sprintf("  Insufficient samples for CI calculation (n<%d)\n", 
                    VALIDATION_CONFIG$min_samples_per_stratum))
        cat("  Descriptive statistics only - not used for conclusions\n")
        
        stratum_results[[stratum_name]] <- list(
          n_total = nrow(stratum_data),
          n_judged = nrow(stratum_judged),
          n_positive = n_positive,
          n_negative = n_negative,
          n_boundary = n_boundary,
          n_undetermined = sum(!stratum_data$judged),
          n_denominator = n_denominator,
          sufficient_samples = FALSE
        )
      }
    } else {
      cat("  No judged samples in this stratum\n")
      stratum_results[[stratum_name]] <- list(
        n_total = nrow(stratum_data),
        n_judged = 0,
        sufficient_samples = FALSE
      )
    }
  }
} else {
  cat("  No high purity samples available for main analysis\n")
  stratum_results <- NULL
}

# --------------------------------------------------------------------------
# Part B: Diagnostic Checks
# --------------------------------------------------------------------------

cat("\n=== PART B: Diagnostic Checks ===\n")

# B1. Classification Symmetry (FPR vs FNR on Judged basis)
cat("\n--- B1. Classification Symmetry (Training Data) ---\n")

# Use training data from Step 5 for symmetry check
training_scores <- step5_data$sample_scores

# Apply B-scheme classification to training data
training_scores$RET_REO_call <- NA_character_
for (i in 1:nrow(training_scores)) {
  if (training_scores$judgment_possible[i]) {
    class_result <- classify_reo(
      training_scores$n_reversals[i],
      training_scores$n_eff[i],
      N, T,
      n_eff_threshold = VALIDATION_CONFIG$n_eff_threshold
    )
    training_scores$RET_REO_call[i] <- class_result$RET_REO_call
  } else {
    training_scores$RET_REO_call[i] <- "Undetermined"
  }
}

# Filter to judged samples only (excluding Boundary and Undetermined)
r0_judged <- training_scores[training_scores$group == "R0" & 
                               training_scores$RET_REO_call %in% c("Positive", "Negative"), ]
r1_judged <- training_scores[training_scores$group == "R1" & 
                               training_scores$RET_REO_call %in% c("Positive", "Negative"), ]

# Calculate error rates (Judged basis)
if (nrow(r0_judged) > 0) {
  r0_fpr <- sum(r0_judged$RET_REO_call == "Positive") / nrow(r0_judged)
  r0_fpr_ci <- prop.test(sum(r0_judged$RET_REO_call == "Positive"), 
                         nrow(r0_judged),
                         conf.level = 0.95, correct = FALSE)$conf.int
  cat(sprintf("  R0 FPR (Positive rate): %.1f%% [%.1f%%, %.1f%%]\n", 
              r0_fpr*100, r0_fpr_ci[1]*100, r0_fpr_ci[2]*100))
  cat(sprintf("    (Based on %d judged R0 samples)\n", nrow(r0_judged)))
} else {
  r0_fpr <- NA
  r0_fpr_ci <- c(NA, NA)
  cat("  R0 FPR: Cannot calculate (no judged samples)\n")
}

if (nrow(r1_judged) > 0) {
  r1_fnr <- sum(r1_judged$RET_REO_call == "Negative") / nrow(r1_judged)
  r1_fnr_ci <- prop.test(sum(r1_judged$RET_REO_call == "Negative"), 
                         nrow(r1_judged),
                         conf.level = 0.95, correct = FALSE)$conf.int
  cat(sprintf("  R1 FNR (Negative rate): %.1f%% [%.1f%%, %.1f%%]\n",
              r1_fnr*100, r1_fnr_ci[1]*100, r1_fnr_ci[2]*100))
  cat(sprintf("    (Based on %d judged R1 samples)\n", nrow(r1_judged)))
} else {
  r1_fnr <- NA
  r1_fnr_ci <- c(NA, NA)
  cat("  R1 FNR: Cannot calculate (no judged samples)\n")
}

# Check CI overlap
if (!is.na(r0_fpr) && !is.na(r1_fnr)) {
  ci_overlap <- (r0_fpr_ci[1] <= r1_fnr_ci[2]) && (r1_fnr_ci[1] <= r0_fpr_ci[2])
  cat(sprintf("  CI overlap: %s\n", ifelse(ci_overlap, "Yes (balanced)", "No (imbalanced)")))
} else {
  ci_overlap <- NA
  cat("  CI overlap: Cannot assess\n")
}

# B2. N_eff Influence
cat("\n--- B2. N_eff Influence ---\n")

if (sum(r_test_results$judged) > 5) {
  # Check correlation between N_eff and REO score
  valid_results <- r_test_results[r_test_results$judged, ]
  
  if (nrow(valid_results) > 2) {
    cor_n_eff <- cor(valid_results$n_eff, valid_results$reo_score, 
                     method = "spearman", use = "complete.obs")
    cat(sprintf("  Spearman correlation (N_eff vs REO-score): %.3f\n", cor_n_eff))
    
    if (abs(cor_n_eff) > 0.4) {
      cat("  ⚠ Warning: Moderate to strong correlation detected\n")
    } else {
      cat("  ✔ Weak correlation - N_eff influence minimal\n")
    }
    
    # Check boundary samples
    n_boundary_samples <- sum(r_test_results$RET_REO_call == "Boundary", na.rm = TRUE)
    n_low_neff <- sum(r_test_results$n_eff < ceiling(VALIDATION_CONFIG$n_eff_threshold * N), na.rm = TRUE)
    
    cat(sprintf("  Boundary samples (n_reversals in {%d, %d}): %d\n", 
                T-1, T, n_boundary_samples))
    cat(sprintf("  Low N_eff samples (Undetermined): %d\n", n_low_neff))
  }
} else {
  cat("  Insufficient judged samples for N_eff analysis\n")
}

# B3. N-1 Jackknife (Stability for Boundary samples)
cat("\n--- B3. Jackknife Stability (Boundary Samples) ---\n")

if (VALIDATION_CONFIG$jackknife_enabled && nrow(r_test_results) > 0) {
  # Identify boundary samples
  boundary_samples <- r_test_results[r_test_results$RET_REO_call == "Boundary", ]
  
  if (nrow(boundary_samples) > 0) {
    cat(sprintf("  Boundary samples found: %d\n", nrow(boundary_samples)))
    
    # Perform jackknife for up to 3 boundary samples
    jackknife_results <- list()
    n_to_test <- min(3, nrow(boundary_samples))
    
    for (i in seq_len(n_to_test)) {
      sample_idx <- which(r_test_results$sample_id == boundary_samples$sample_id[i])
      sample_reo <- reo_matrix_r_test[, sample_idx]
      original_call <- boundary_samples$RET_REO_call[i]
      
      # N-1 classifications
      jackknife_calls <- character(N)
      
      for (j in 1:N) {
        # Remove pair j
        sample_reo_j <- sample_reo[-j]
        r0_directions_j <- panel_info$r0_directions[-j]
        N_j <- N - 1
        T_j <- ceiling(N_j * 0.4)
        
        # Recalculate with N-1 pairs
        valid_mask_j <- abs(sample_reo_j) >= CONFIG$dead_zone
        n_eff_j <- sum(valid_mask_j, na.rm = TRUE)
        
        if (n_eff_j >= ceiling(VALIDATION_CONFIG$n_eff_threshold * N_j)) {
          n_reversals_j <- sum(sign(sample_reo_j[valid_mask_j]) != 
                                 r0_directions_j[valid_mask_j], na.rm = TRUE)
        } else {
          n_reversals_j <- NA_integer_
        }
        
        # Apply B-scheme classification
        class_result_j <- classify_reo(n_reversals_j, n_eff_j, N_j, T_j,
                                       n_eff_threshold = VALIDATION_CONFIG$n_eff_threshold)
        jackknife_calls[j] <- class_result_j$RET_REO_call
      }
      
      # Calculate consistency
      consistency <- sum(jackknife_calls == original_call) / N
      
      cat(sprintf("    Sample %d: %.0f%% consistent (stays in Boundary)\n", 
                  i, consistency * 100))
      
      jackknife_results[[i]] <- list(
        sample_id = boundary_samples$sample_id[i],
        original = original_call,
        consistency = consistency,
        jackknife_calls = jackknife_calls
      )
    }
  } else {
    cat("  No boundary samples found\n")
    jackknife_results <- NULL
  }
} else {
  cat("  Jackknife analysis skipped\n")
  jackknife_results <- NULL
}

# --------------------------------------------------------------------------
# Part C: External Specificity (BRAF samples - Internal QC only)
# --------------------------------------------------------------------------

cat("\n=== PART C: External Specificity (BRAF - Internal QC) ===\n")
cat("  Note: This is for internal validation only, not for publication\n")

if (VALIDATION_CONFIG$braf_validation) {
  # Get BRAF test samples
  braf_cases <- case_master[case_master$group == "B_test", ]
  
  if (nrow(braf_cases) > 0) {
    cat(sprintf("  BRAF test cases: %d\n", nrow(braf_cases)))
    
    # Extract BRAF samples
    braf_samples <- extract_r_test_samples(
      braf_cases$case_submitter_id,
      sample_metadata,
      "Primary Tumor"
    )
    
    cat(sprintf("  BRAF tumor samples found: %d\n", length(braf_samples$samples)))
    
    if (length(braf_samples$samples) > 0) {
      # Calculate REO scores for BRAF samples
      braf_counts <- assay(se, "counts")[, braf_samples$samples, drop = FALSE]
      
      # TPM calculation
      rpk_braf <- braf_counts / 1000
      scaling_factor_braf <- colSums(rpk_braf) / 1e6
      tpm_braf <- sweep(rpk_braf, 2, scaling_factor_braf, "/")
      log2_tpm_braf <- log2(tpm_braf + 1)
      
      # Calculate REO and classify
      braf_results <- data.frame(
        case_id = braf_samples$cases,
        sample_id = braf_samples$samples,
        RET_REO_call = character(length(braf_samples$samples)),
        reo_score = numeric(length(braf_samples$samples)),
        n_eff = integer(length(braf_samples$samples)),
        n_reversals = integer(length(braf_samples$samples)),
        judged = logical(length(braf_samples$samples)),
        stringsAsFactors = FALSE
      )
      
      for (i in seq_along(braf_samples$samples)) {
        # Calculate REO for each pair
        sample_tpm <- log2_tpm_braf[, i]
        reo_values <- numeric(N)
        
        for (j in 1:N) {
          gene_up <- panel_pairs$gene_up[j]
          gene_down <- panel_pairs$gene_down[j]
          
          if (gene_up %in% rownames(log2_tpm_braf) && 
              gene_down %in% rownames(log2_tpm_braf)) {
            reo_values[j] <- sample_tpm[gene_up] - sample_tpm[gene_down]
          } else {
            reo_values[j] <- NA_real_
          }
        }
        
        # Calculate metrics
        valid_mask <- abs(reo_values) >= CONFIG$dead_zone
        n_eff <- sum(valid_mask, na.rm = TRUE)
        
        if (n_eff >= ceiling(VALIDATION_CONFIG$n_eff_threshold * N)) {
          n_reversals <- sum(sign(reo_values[valid_mask]) != 
                               panel_info$r0_directions[valid_mask], na.rm = TRUE)
        } else {
          n_reversals <- NA_integer_
        }
        
        reo_score <- ifelse(!is.na(n_reversals), n_reversals / N, NA_real_)
        
        # Apply B-scheme classification
        class_result <- classify_reo(n_reversals, n_eff, N, T,
                                     n_eff_threshold = VALIDATION_CONFIG$n_eff_threshold)
        
        braf_results[i, c("RET_REO_call", "reo_score", "n_eff", 
                          "n_reversals", "judged")] <- 
          list(class_result$RET_REO_call, reo_score, n_eff, n_reversals,
               class_result$judged)
      }
      
      # Calculate RET-REO positive rate in BRAF
      braf_judged <- braf_results[braf_results$judged & 
                                    braf_results$RET_REO_call %in% c("Positive", "Negative"), ]
      
      cat(sprintf("\n  BRAF samples classification:\n"))
      cat(sprintf("    Total: %d\n", nrow(braf_results)))
      cat(sprintf("    Judged (Pos/Neg only): %d\n", nrow(braf_judged)))
      
      if (nrow(braf_judged) > 0) {
        n_positive_braf <- sum(braf_judged$RET_REO_call == "Positive")
        pos_rate_braf <- n_positive_braf / nrow(braf_judged)
        
        cat(sprintf("  RET-REO Positive rate in BRAF: %d/%d = %.1f%%\n",
                    n_positive_braf, nrow(braf_judged), pos_rate_braf*100))
        cat("  (Exploratory finding for internal QC)\n")
        
        # Wilson CI if sufficient samples
        if (nrow(braf_judged) >= 5) {
          wilson_ci_braf <- prop.test(n_positive_braf, nrow(braf_judged),
                                      conf.level = 0.95, correct = FALSE)$conf.int
          cat(sprintf("  95%% CI: [%.1f%%, %.1f%%]\n",
                      wilson_ci_braf[1]*100, wilson_ci_braf[2]*100))
        }
        
        # Show full breakdown
        cat("\n  Full BRAF classification:\n")
        print(table(braf_results$RET_REO_call))
      }
    } else {
      braf_results <- NULL
      cat("  No BRAF samples available for testing\n")
    }
  } else {
    braf_results <- NULL
    cat("  No BRAF test cases in dataset\n")
  }
} else {
  braf_results <- NULL
  cat("  BRAF validation skipped\n")
}

# --------------------------------------------------------------------------
# Part D: Sensitivity Analysis (Purity stratification)
# --------------------------------------------------------------------------

cat("\n=== PART D: Sensitivity Analysis ===\n")

# D1. Performance by purity class
cat("\n--- D1. REO-score by purity class ---\n")

purity_classes <- unique(r_test_results$purity_class)
purity_performance <- list()

for (pc in purity_classes) {
  pc_data <- r_test_results[r_test_results$purity_class == pc & 
                              r_test_results$judged, ]
  
  if (nrow(pc_data) > 0) {
    cat(sprintf("\n%s purity:\n", pc))
    cat(sprintf("  n = %d\n", nrow(pc_data)))
    cat(sprintf("  REO-score: mean=%.3f, median=%.3f, sd=%.3f\n",
                mean(pc_data$reo_score, na.rm=TRUE),
                median(pc_data$reo_score, na.rm=TRUE),
                sd(pc_data$reo_score, na.rm=TRUE)))
    
    # Classification breakdown
    class_table <- table(pc_data$RET_REO_call)
    cat("  Classification: ")
    cat(paste(names(class_table), "=", class_table, collapse=", "))
    cat("\n")
    
    purity_performance[[pc]] <- list(
      n = nrow(pc_data),
      mean_score = mean(pc_data$reo_score, na.rm=TRUE),
      median_score = median(pc_data$reo_score, na.rm=TRUE),
      sd_score = sd(pc_data$reo_score, na.rm=TRUE),
      classifications = class_table
    )
  }
}

# Check for score stability across purity levels
if (length(purity_performance) > 1) {
  score_diff <- abs(purity_performance[["High"]]$median_score - 
                      purity_performance[["Low"]]$median_score)
  if (score_diff < 0.1) {
    cat("  ✔ REO-scores appear stable across purity levels\n")
  } else {
    cat("  ⚠ REO-scores differ between purity levels\n")
  }
}

# D2. Threshold sensitivity (T±1)
cat("\n--- D2. Threshold sensitivity analysis ---\n")
cat("  Testing alternative thresholds: T-1, T, T+1\n")

for (delta in c(-1, 0, 1)) {
  T_alt <- T + delta
  
  # Recount with alternative threshold
  n_pos_alt <- sum(r_test_results$n_reversals >= T_alt + 1 & 
                     r_test_results$judged, na.rm=TRUE)
  n_neg_alt <- sum(r_test_results$n_reversals <= T_alt - 2 & 
                     r_test_results$judged, na.rm=TRUE)
  n_total_alt <- sum(r_test_results$judged)
  
  if (n_total_alt > 0) {
    pos_rate_alt <- n_pos_alt / n_total_alt
    cat(sprintf("  T=%d: Positive rate = %.1f%% (%d/%d)\n",
                T_alt, pos_rate_alt*100, n_pos_alt, n_total_alt))
  }
}

# --------------------------------------------------------------------------
# 6.7 Save Step 6 results
# --------------------------------------------------------------------------

cat("\n--- 6.7 Saving Step 6 results ---\n")

step6_data <- list(
  # Configuration
  config = CONFIG,
  validation_config = VALIDATION_CONFIG,
  
  # Panel info from Step 5
  panel_info = panel_info,
  panel_pairs = panel_pairs,
  N = N,
  T = T,
  
  # R_test results
  r_test_results = r_test_results,
  reo_matrix_r_test = reo_matrix_r_test,
  
  # Purity analysis
  purity_results = if (!is.null(purity_values)) {
    list(
      purity_values = purity_values,
      quality_stats = quality_stats
    )
  } else NULL,
  
  # Main analysis (stratified)
  main_analysis = stratum_results,
  
  # Diagnostic checks
  diagnostics = list(
    symmetry = list(
      r0_fpr = r0_fpr,
      r0_fpr_ci = r0_fpr_ci,
      r1_fnr = r1_fnr,
      r1_fnr_ci = r1_fnr_ci,
      ci_overlap = ci_overlap
    ),
    n_eff_influence = if (exists("cor_n_eff")) {
      list(correlation = cor_n_eff)
    } else NULL,
    jackknife = jackknife_results
  ),
  
  # External validation
  braf_validation = braf_results,
  
  # Sensitivity analysis
  purity_performance = purity_performance,
  
  # Summary statistics
  summary_stats = list(
    n_r_test_total = nrow(r_test_cases),
    n_r_test_analyzed = nrow(r_test_results),
    n_high_purity = sum(r_test_results$purity_class == "High"),
    n_judged = sum(r_test_results$judged),
    n_positive = sum(r_test_results$RET_REO_call == "Positive", na.rm=TRUE),
    n_negative = sum(r_test_results$RET_REO_call == "Negative", na.rm=TRUE),
    n_boundary = sum(r_test_results$RET_REO_call == "Boundary", na.rm=TRUE),
    n_undetermined = sum(r_test_results$RET_REO_call == "Undetermined", na.rm=TRUE),
    timestamp = Sys.time()
  ),
  
  # Metadata
  step = "Step 6 Complete",
  classification_scheme = "B-scheme (n_reversals >= T+1 → Positive)"
)

saveRDS(step6_data, paste0(paths$processed, "reo_step6_validation.rds"))
cat("  Step 6 data saved to: reo_step6_validation.rds\n")

# Export summary CSV
validation_summary <- data.frame(
  metric = character(),
  value = character(),
  stringsAsFactors = FALSE
)

# Add main results
if (!is.null(stratum_results)) {
  for (stratum in names(stratum_results)) {
    sr <- stratum_results[[stratum]]
    if (!is.null(sr$sufficient_samples) && sr$sufficient_samples) {
      validation_summary <- rbind(validation_summary,
                                  data.frame(
                                    metric = paste0(stratum, "_positive_rate"),
                                    value = sprintf("%.1f%% [%.1f-%.1f]", 
                                                    sr$positive_rate * 100,
                                                    sr$ci_lower * 100,
                                                    sr$ci_upper * 100)
                                  ),
                                  data.frame(
                                    metric = paste0(stratum, "_breakdown"),
                                    value = sprintf("Pos=%d, Neg=%d, Bnd=%d, Und=%d",
                                                    sr$n_positive, sr$n_negative, 
                                                    sr$n_boundary, sr$n_undetermined)
                                  )
      )
    } else if (!is.null(sr$n_judged)) {
      validation_summary <- rbind(validation_summary,
                                  data.frame(
                                    metric = paste0(stratum, "_status"),
                                    value = sprintf("n=%d (insufficient for CI)", sr$n_judged)
                                  )
      )
    }
  }
}

# Add diagnostic results
if (!is.na(r0_fpr)) {
  validation_summary <- rbind(validation_summary,
                              data.frame(
                                metric = "R0_FPR",
                                value = sprintf("%.1f%% [%.1f-%.1f]", 
                                                r0_fpr * 100, r0_fpr_ci[1] * 100, r0_fpr_ci[2] * 100)
                              )
  )
}

if (!is.na(r1_fnr)) {
  validation_summary <- rbind(validation_summary,
                              data.frame(
                                metric = "R1_FNR", 
                                value = sprintf("%.1f%% [%.1f-%.1f]",
                                                r1_fnr * 100, r1_fnr_ci[1] * 100, r1_fnr_ci[2] * 100)
                              )
  )
}

if (!is.na(ci_overlap)) {
  validation_summary <- rbind(validation_summary,
                              data.frame(
                                metric = "CI_overlap",
                                value = ifelse(ci_overlap, "Yes", "No")
                              )
  )
}

# Add overall classification summary
validation_summary <- rbind(validation_summary,
                            data.frame(
                              metric = "Total_samples",
                              value = as.character(nrow(r_test_results))
                            ),
                            data.frame(
                              metric = "Classification_summary",
                              value = sprintf("Pos=%d, Neg=%d, Bnd=%d, Und=%d",
                                              step6_data$summary_stats$n_positive,
                                              step6_data$summary_stats$n_negative,
                                              step6_data$summary_stats$n_boundary,
                                              step6_data$summary_stats$n_undetermined)
                            )
)

write.csv(validation_summary, 
          paste0(paths$output, "step6_validation_summary.csv"),
          row.names = FALSE)
cat("  Validation summary exported to: step6_validation_summary.csv\n")

# =============================================================================
# STEP 6 SUMMARY
# =============================================================================

cat("\n=== STEP 6 SUMMARY ===\n")
cat(sprintf("R_test cases analyzed: %d\n", nrow(r_test_cases)))
cat(sprintf("Samples with purity estimation: %d\n", 
            sum(!is.na(r_test_results$purity))))
cat(sprintf("High purity samples: %d\n", 
            sum(r_test_results$purity_class == "High")))
cat(sprintf("\nClassification breakdown (all samples):\n"))
cat(sprintf("  Positive: %d\n", step6_data$summary_stats$n_positive))
cat(sprintf("  Negative: %d\n", step6_data$summary_stats$n_negative))
cat(sprintf("  Boundary: %d\n", step6_data$summary_stats$n_boundary))
cat(sprintf("  Undetermined: %d\n", step6_data$summary_stats$n_undetermined))
cat(sprintf("  Judged (Pos+Neg): %d\n", 
            step6_data$summary_stats$n_positive + step6_data$summary_stats$n_negative))

# Performance summary
if (!is.null(stratum_results) && length(stratum_results) > 0) {
  cat("\nMain Analysis Results (High Purity, 2-layer POC):\n")
  for (stratum in names(stratum_results)) {
    sr <- stratum_results[[stratum]]
    if (!is.null(sr$sufficient_samples) && sr$sufficient_samples) {
      cat(sprintf("  %s: Positive rate = %.1f%% [%.1f%%, %.1f%%]\n",
                  stratum, sr$positive_rate*100, 
                  sr$ci_lower*100, sr$ci_upper*100))
    } else {
      cat(sprintf("  %s: n=%d (insufficient for conclusions)\n", 
                  stratum, ifelse(!is.null(sr$n_judged), sr$n_judged, 0)))
    }
  }
}

# Diagnostic summary
cat("\nDiagnostic Checks:\n")
if (!is.na(ci_overlap)) {
  cat(sprintf("  Classification symmetry: %s\n", 
              ifelse(ci_overlap, "Balanced", "Imbalanced")))
} else {
  cat("  Classification symmetry: Cannot assess\n")
}

cat(sprintf("  External validation (BRAF): %s\n",
            ifelse(!is.null(braf_results), "Completed (internal QC)", "Not performed")))
cat(sprintf("  Jackknife stability: %s\n",
            ifelse(!is.null(jackknife_results), "Assessed", "Not assessed")))

cat(sprintf("\nTimestamp: %s\n", format(step6_data$summary_stats$timestamp)))
cat("\n✔ Step 6 validation complete. REO panel performance validated.\n")
cat("✔ Classification scheme: B-scheme (Positive/Negative/Boundary/Undetermined)\n")