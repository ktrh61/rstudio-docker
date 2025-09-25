# 26_reo_step6_rtest.R - REO Step 6: R_test Validation and Robustness
# Purpose: Validate REO panel on R_test samples with comprehensive diagnostics
# Input: reo_step5_data.rds, thyr_case_master_full.rds, thyr_se_strand2_nonzero.rds
# Output: reo_step6_validation.rds
# Version: v1.0
# Date: 2025-01-26
# Note: High-purity focus with multi-layer diagnostics

source("analysis_v7/setup.R")

cat("\n=== REO STEP 6: R_test Validation and Robustness ===\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

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
  
  # POC stratification
  poc_cutoff = 33.3,             # POC stratification point
  
  # Minimum sample sizes
  min_samples_per_stratum = 5,   # Minimum for CI calculation
  
  # Jackknife settings
  jackknife_enabled = TRUE,      # Perform N-1 jackknife
  
  # BRAF validation
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

cat(sprintf("  Panel: %d pairs, threshold T=%d\n", N, T))

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

# Get sample metadata
sample_metadata <- as.data.frame(colData(se))

# --------------------------------------------------------------------------
# 6.3 Extract R_test samples (following safe pattern from 06_purity_analysis.R)
# --------------------------------------------------------------------------

cat("\n--- 6.3 Extracting R_test samples ---\n")

# Get R_test cases
r_test_cases <- case_master[case_master$group == "R_test", ]
cat(sprintf("  R_test cases: %d\n", nrow(r_test_cases)))

# Separate paired and unpaired
r_test_paired <- r_test_cases[r_test_cases$has_pair == TRUE, ]
r_test_unpaired <- r_test_cases[r_test_cases$has_pair == FALSE, ]

cat(sprintf("  Paired: %d, Unpaired: %d\n", 
            nrow(r_test_paired), nrow(r_test_unpaired)))

# Safe sample extraction function (from 06_purity_analysis.R pattern)
extract_r_test_samples <- function(case_ids, sample_meta, sample_type = "Primary Tumor") {
  sample_ids <- character(length(case_ids))
  
  for (i in seq_along(case_ids)) {
    case_id <- case_ids[i]
    
    # Find sample with explicit conditions
    sample_id <- sample_meta$sample_submitter_id[
      sample_meta$case_submitter_id == case_id &
        sample_meta$sample_type == sample_type &
        grepl("_merged", sample_meta$sample_submitter_id)
    ]
    
    # Take first if multiple (shouldn't happen with _merged)
    if (length(sample_id) > 0) {
      sample_ids[i] <- sample_id[1]
    } else {
      sample_ids[i] <- NA_character_
    }
  }
  
  # Return only valid samples
  valid_samples <- sample_ids[!is.na(sample_ids)]
  valid_cases <- case_ids[!is.na(sample_ids)]
  
  return(list(
    samples = valid_samples,
    cases = valid_cases
  ))
}

# Extract all R_test tumor samples
r_test_all <- extract_r_test_samples(
  r_test_cases$case_submitter_id, 
  sample_metadata, 
  "Primary Tumor"
)

cat(sprintf("  R_test tumor samples found: %d\n", length(r_test_all$samples)))

# --------------------------------------------------------------------------
# 6.4 Tumor purity estimation for paired samples
# --------------------------------------------------------------------------

cat("\n--- 6.4 Estimating tumor purity (paired samples only) ---\n")

if (nrow(r_test_paired) > 0) {
  # Extract paired samples (following 06_purity_analysis.R exactly)
  extract_paired_for_contamde <- function(case_ids, sample_meta) {
    tumor_ids <- character(length(case_ids))
    normal_ids <- character(length(case_ids))
    
    for (i in seq_along(case_ids)) {
      case_id <- case_ids[i]
      
      # Find tumor sample
      tumor_sample <- sample_meta$sample_submitter_id[
        sample_meta$case_submitter_id == case_id &
          sample_meta$sample_type == "Primary Tumor" &
          grepl("_merged", sample_meta$sample_submitter_id)
      ]
      
      # Find normal sample
      normal_sample <- sample_meta$sample_submitter_id[
        sample_meta$case_submitter_id == case_id &
          sample_meta$sample_type == "Solid Tissue Normal" &
          grepl("_merged", sample_meta$sample_submitter_id)
      ]
      
      # Store (may be NA)
      tumor_ids[i] <- if(length(tumor_sample) > 0) tumor_sample[1] else NA
      normal_ids[i] <- if(length(normal_sample) > 0) normal_sample[1] else NA
    }
    
    # Keep only valid pairs
    valid <- !is.na(tumor_ids) & !is.na(normal_ids)
    
    return(list(
      tumor = tumor_ids[valid],
      normal = normal_ids[valid],
      cases = case_ids[valid]
    ))
  }
  
  # Get paired samples
  paired_samples <- extract_paired_for_contamde(
    r_test_paired$case_submitter_id,
    sample_metadata
  )
  
  cat(sprintf("  Valid pairs for contamDE: %d\n", length(paired_samples$cases)))
  
  if (length(paired_samples$cases) > 0) {
    # Gene filtering (from 06_purity_analysis.R)
    gene_info <- as.data.frame(rowData(se))
    
    # Protein coding only
    is_protein_coding <- gene_info$gene_type == "protein_coding"
    
    # Extract counts (Normal first, then Tumor - contamDE convention)
    all_sample_ids <- c(paired_samples$normal, paired_samples$tumor)
    counts_for_purity <- assay(se, "counts")[, all_sample_ids, drop = FALSE]
    
    # Apply filterByExpr
    n_pairs <- length(paired_samples$tumor)
    group <- factor(c(rep("Normal", n_pairs), rep("Tumor", n_pairs)))
    keep_expr <- edgeR::filterByExpr(counts_for_purity, group = group)
    
    # Combine filters
    keep_genes <- is_protein_coding & keep_expr
    
    cat(sprintf("  Genes for purity: %d protein coding, %d after filter\n",
                sum(is_protein_coding), sum(keep_genes)))
    
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
    cat("  No valid pairs for contamDE\n")
    purity_values <- NULL
    quality_stats <- NULL
  }
} else {
  cat("  No paired samples in R_test\n")
  purity_values <- NULL
  quality_stats <- NULL
}

# --------------------------------------------------------------------------
# 6.5 Classify samples by purity
# --------------------------------------------------------------------------

cat("\n--- 6.5 Classifying samples by purity ---\n")

# Create sample classification
r_test_sample_info <- data.frame(
  case_id = r_test_all$cases,
  sample_id = r_test_all$samples,
  stringsAsFactors = FALSE
)

# Add POC values
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

# POC classification
r_test_sample_info$poc_class <- ifelse(
  r_test_sample_info$poc <= VALIDATION_CONFIG$poc_cutoff,
  "Low_POC",
  "High_POC"
)

# Summary
cat("\nPurity classification:\n")
print(table(r_test_sample_info$purity_class))

cat("\nPOC classification:\n")
print(table(r_test_sample_info$poc_class))

cat("\nCross-tabulation (Purity × POC):\n")
print(table(r_test_sample_info$purity_class, r_test_sample_info$poc_class))

# --------------------------------------------------------------------------
# 6.6 Calculate REO scores for all R_test samples
# --------------------------------------------------------------------------

cat("\n--- 6.6 Calculating REO scores for R_test samples ---\n")

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

# Apply classification function from Step 5
r_test_classifications <- data.frame(
  sample_id = r_test_sample_info$sample_id,
  case_id = r_test_sample_info$case_id,
  classification = character(nrow(r_test_sample_info)),
  reo_score = numeric(nrow(r_test_sample_info)),
  n_reversals = integer(nrow(r_test_sample_info)),
  n_eff = integer(nrow(r_test_sample_info)),
  judgment_possible = logical(nrow(r_test_sample_info)),
  stringsAsFactors = FALSE
)

# Classify each sample
for (i in seq_len(nrow(r_test_sample_info))) {
  sample_reo <- reo_matrix_r_test[, i]
  
  # Manual classification (matching Step 5 logic)
  valid_mask <- abs(sample_reo) >= CONFIG$dead_zone
  n_eff <- sum(valid_mask, na.rm = TRUE)
  
  if (n_eff >= ceiling(0.7 * N)) {  # N_eff threshold from Step 5
    valid_values <- sample_reo[valid_mask]
    valid_directions <- sign(valid_values)
    valid_r0_directions <- panel_info$r0_directions[valid_mask]
    
    n_reversals <- sum(valid_directions != valid_r0_directions, na.rm = TRUE)
    reo_score <- n_reversals / N
    
    classification <- ifelse(n_reversals >= T, "R1-like", "R0-like")
    judgment_possible <- TRUE
  } else {
    classification <- "Undetermined"
    reo_score <- NA_real_
    n_reversals <- NA_integer_
    judgment_possible <- FALSE
  }
  
  r_test_classifications[i, c("classification", "reo_score", "n_reversals", 
                              "n_eff", "judgment_possible")] <- 
    list(classification, reo_score, n_reversals, n_eff, judgment_possible)
}

# Merge with sample info
r_test_results <- cbind(r_test_sample_info, r_test_classifications[, -c(1:2)])

cat(sprintf("  Classifications complete: %d samples\n", nrow(r_test_results)))
cat("  Overall classification:\n")
print(table(r_test_results$classification))

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
  # Stratify by POC
  strata <- list(
    Low_POC = high_purity_results[high_purity_results$poc_class == "Low_POC", ],
    High_POC = high_purity_results[high_purity_results$poc_class == "High_POC", ]
  )
  
  # Calculate performance by stratum
  stratum_results <- list()
  
  for (stratum_name in names(strata)) {
    stratum_data <- strata[[stratum_name]]
    stratum_judged <- stratum_data[stratum_data$judgment_possible, ]
    
    cat(sprintf("\n%s stratum:\n", stratum_name))
    cat(sprintf("  Total: %d, Judged: %d\n", 
                nrow(stratum_data), nrow(stratum_judged)))
    
    if (nrow(stratum_judged) >= VALIDATION_CONFIG$min_samples_per_stratum) {
      # Calculate positivity rate
      n_positive <- sum(stratum_judged$classification == "R1-like")
      n_total <- nrow(stratum_judged)
      pos_rate <- n_positive / n_total
      
      # Wilson CI (no continuity correction)
      wilson_ci <- prop.test(n_positive, n_total, 
                             conf.level = VALIDATION_CONFIG$ci_level,
                             correct = FALSE)$conf.int
      
      cat(sprintf("  R1-like rate: %d/%d = %.1f%% [%.1f%%, %.1f%%]\n",
                  n_positive, n_total, pos_rate*100,
                  wilson_ci[1]*100, wilson_ci[2]*100))
      
      stratum_results[[stratum_name]] <- list(
        n = n_total,
        n_positive = n_positive,
        rate = pos_rate,
        ci_lower = wilson_ci[1],
        ci_upper = wilson_ci[2],
        insufficient = FALSE  # Add this field
      )
    } else {
      cat("  Insufficient samples for CI calculation (descriptive only)\n")
      if (nrow(stratum_judged) > 0) {
        cat(sprintf("  Classifications: %s\n", 
                    paste(table(stratum_judged$classification), collapse=", ")))
      }
      stratum_results[[stratum_name]] <- list(
        n = nrow(stratum_judged),
        insufficient = TRUE
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

# B1. Classification Symmetry (FPR vs FNR)
cat("\n--- B1. Classification Symmetry ---\n")

# Use training data from Step 5 for symmetry check
training_scores <- step5_data$sample_scores
r0_training <- training_scores[training_scores$group == "R0" & 
                                 training_scores$judgment_possible, ]
r1_training <- training_scores[training_scores$group == "R1" & 
                                 training_scores$judgment_possible, ]

# Calculate error rates
r0_fpr <- sum(r0_training$n_reversals >= T) / nrow(r0_training)  # False Positive Rate
r1_fnr <- sum(r1_training$n_reversals < T) / nrow(r1_training)   # False Negative Rate

# Wilson CIs
r0_fpr_ci <- prop.test(sum(r0_training$n_reversals >= T), nrow(r0_training),
                       conf.level = 0.95, correct = FALSE)$conf.int
r1_fnr_ci <- prop.test(sum(r1_training$n_reversals < T), nrow(r1_training),
                       conf.level = 0.95, correct = FALSE)$conf.int

cat(sprintf("  R0 FPR: %.1f%% [%.1f%%, %.1f%%]\n", 
            r0_fpr*100, r0_fpr_ci[1]*100, r0_fpr_ci[2]*100))
cat(sprintf("  R1 FNR: %.1f%% [%.1f%%, %.1f%%]\n",
            r1_fnr*100, r1_fnr_ci[1]*100, r1_fnr_ci[2]*100))

# Check CI overlap
ci_overlap <- (r0_fpr_ci[1] <= r1_fnr_ci[2]) && (r1_fnr_ci[1] <= r0_fpr_ci[2])
cat(sprintf("  CI overlap: %s\n", ifelse(ci_overlap, "Yes (balanced)", "No (imbalanced)")))

# B2. N_eff Influence
cat("\n--- B2. N_eff Influence ---\n")

if (nrow(r_test_results) > 5) {
  # Check correlation between N_eff and REO score
  valid_results <- r_test_results[r_test_results$judgment_possible, ]
  
  if (nrow(valid_results) > 2) {
    cor_n_eff <- cor(valid_results$n_eff, valid_results$reo_score, 
                     method = "spearman", use = "complete.obs")
    cat(sprintf("  Spearman correlation (N_eff vs REO-score): %.3f\n", cor_n_eff))
    
    if (abs(cor_n_eff) > 0.4) {
      cat("  ⚠ Warning: Moderate to strong correlation detected\n")
    } else {
      cat("  ✓ Weak correlation - N_eff influence minimal\n")
    }
    
    # Check boundary samples
    n_eff_low <- sum(r_test_results$n_eff < 7, na.rm = TRUE)
    n_boundary <- sum(abs(r_test_results$reo_score - 0.4) < 0.1, na.rm = TRUE)
    cat(sprintf("  Low N_eff (<7): %d samples\n", n_eff_low))
    cat(sprintf("  Boundary samples (score 0.3-0.5): %d samples\n", n_boundary))
  }
} else {
  cat("  Insufficient samples for N_eff analysis\n")
}

# B3. N-1 Jackknife (Stability)
cat("\n--- B3. Jackknife Stability (Boundary Samples) ---\n")

if (VALIDATION_CONFIG$jackknife_enabled && nrow(r_test_results) > 0) {
  # Identify boundary samples (REO score near threshold)
  boundary_samples <- r_test_results[
    r_test_results$judgment_possible &
      abs(r_test_results$n_reversals - T) <= 1,  # Within 1 reversal of threshold
  ]
  
  if (nrow(boundary_samples) > 0) {
    cat(sprintf("  Boundary samples found: %d\n", nrow(boundary_samples)))
    
    # Perform jackknife for each boundary sample
    jackknife_results <- list()
    
    for (i in seq_len(min(3, nrow(boundary_samples)))) {  # Test up to 3 samples
      sample_idx <- which(r_test_results$sample_id == boundary_samples$sample_id[i])
      sample_reo <- reo_matrix_r_test[, sample_idx]
      original_class <- boundary_samples$classification[i]
      
      # N-1 classifications
      jackknife_classes <- character(N)
      
      for (j in 1:N) {
        # Remove pair j
        sample_reo_j <- sample_reo[-j]
        r0_directions_j <- panel_info$r0_directions[-j]
        N_j <- N - 1
        T_j <- ceiling(N_j * 0.4)
        
        # Reclassify
        valid_mask_j <- abs(sample_reo_j) >= CONFIG$dead_zone
        n_eff_j <- sum(valid_mask_j, na.rm = TRUE)
        
        if (n_eff_j >= ceiling(0.7 * N_j)) {
          n_reversals_j <- sum(sign(sample_reo_j[valid_mask_j]) != 
                                 r0_directions_j[valid_mask_j], na.rm = TRUE)
          jackknife_classes[j] <- ifelse(n_reversals_j >= T_j, "R1-like", "R0-like")
        } else {
          jackknife_classes[j] <- "Undetermined"
        }
      }
      
      # Calculate stability
      consistency <- sum(jackknife_classes == original_class) / N
      
      cat(sprintf("    Sample %d: %.0f%% consistent\n", i, consistency * 100))
      
      jackknife_results[[i]] <- list(
        sample_id = boundary_samples$sample_id[i],
        original = original_class,
        consistency = consistency
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
# Part C: External Specificity (BRAF samples)
# --------------------------------------------------------------------------

cat("\n=== PART C: External Specificity (BRAF) ===\n")

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
        classification = character(length(braf_samples$samples)),
        reo_score = numeric(length(braf_samples$samples)),
        n_eff = integer(length(braf_samples$samples)),
        judgment_possible = logical(length(braf_samples$samples)),
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
        
        # Classify
        valid_mask <- abs(reo_values) >= CONFIG$dead_zone
        n_eff <- sum(valid_mask, na.rm = TRUE)
        
        if (n_eff >= ceiling(0.7 * N)) {
          n_reversals <- sum(sign(reo_values[valid_mask]) != 
                               panel_info$r0_directions[valid_mask], na.rm = TRUE)
          reo_score <- n_reversals / N
          classification <- ifelse(n_reversals >= T, "R1-like", "R0-like")
          judgment_possible <- TRUE
        } else {
          classification <- "Undetermined"
          reo_score <- NA_real_
          judgment_possible <- FALSE
        }
        
        braf_results[i, c("classification", "reo_score", "n_eff", "judgment_possible")] <- 
          list(classification, reo_score, n_eff, judgment_possible)
      }
      
      # Calculate specificity
      braf_judged <- braf_results[braf_results$judgment_possible, ]
      
      cat(sprintf("\n  BRAF samples judged: %d/%d\n", 
                  nrow(braf_judged), nrow(braf_results)))
      
      if (nrow(braf_judged) > 0) {
        n_positive_braf <- sum(braf_judged$classification == "R1-like")
        pos_rate_braf <- n_positive_braf / nrow(braf_judged)
        
        # Wilson CI
        if (nrow(braf_judged) >= 5) {
          wilson_ci_braf <- prop.test(n_positive_braf, nrow(braf_judged),
                                      conf.level = 0.95, correct = FALSE)$conf.int
          cat(sprintf("  False positive rate: %d/%d = %.1f%% [%.1f%%, %.1f%%]\n",
                      n_positive_braf, nrow(braf_judged), pos_rate_braf*100,
                      wilson_ci_braf[1]*100, wilson_ci_braf[2]*100))
        } else {
          cat(sprintf("  False positive rate: %d/%d = %.1f%% (CI not calculated, n<%d)\n",
                      n_positive_braf, nrow(braf_judged), pos_rate_braf*100,
                      VALIDATION_CONFIG$min_samples_per_stratum))
        }
        
        # Classification table
        cat("\n  BRAF classifications:\n")
        print(table(braf_judged$classification))
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
# Part D: QC Anchor Drift
# --------------------------------------------------------------------------

cat("\n=== PART D: QC Anchor Drift ===\n")

if (length(qc_anchors) > 0) {
  # Calculate anchor expression in R_test samples
  anchor_expression <- log2_tpm_r_test[qc_anchors, , drop = FALSE]
  
  # Group by purity class
  purity_groups <- unique(r_test_results$purity_class)
  
  cat("  Anchor drift by purity class:\n")
  
  for (purity_class in purity_groups) {
    samples_in_class <- r_test_results$sample_id[
      r_test_results$purity_class == purity_class
    ]
    
    if (length(samples_in_class) > 1) {
      class_expression <- anchor_expression[, samples_in_class, drop = FALSE]
      
      # Calculate median absolute deviation across samples
      anchor_mads <- apply(class_expression, 1, mad)
      
      cat(sprintf("    %s purity: median MAD = %.3f, 95th percentile = %.3f\n",
                  purity_class, median(anchor_mads), quantile(anchor_mads, 0.95)))
      
      # Check against threshold
      if (median(anchor_mads) > CONFIG$dead_zone / 2) {
        cat("      ⚠ Warning: Potential drift detected\n")
      }
    }
  }
} else {
  cat("  No QC anchors available\n")
}

# --------------------------------------------------------------------------
# Sensitivity Analysis: Performance by purity class
# --------------------------------------------------------------------------

cat("\n=== Sensitivity Analysis: Performance by Purity Class ===\n")

purity_classes <- c("High", "Low", "Unknown")
purity_performance <- list()

for (pc in purity_classes) {
  pc_results <- r_test_results[r_test_results$purity_class == pc, ]
  pc_judged <- pc_results[pc_results$judgment_possible, ]
  
  cat(sprintf("\n%s purity:\n", pc))
  cat(sprintf("  Total: %d, Judged: %d (%.1f%%)\n",
              nrow(pc_results), nrow(pc_judged),
              nrow(pc_judged)/max(1, nrow(pc_results))*100))
  
  if (nrow(pc_judged) > 0) {
    # REO score statistics
    cat(sprintf("  REO score: mean=%.3f, median=%.3f\n",
                mean(pc_judged$reo_score), median(pc_judged$reo_score)))
    
    # Classification distribution
    cat("  Classifications:\n")
    print(table(pc_judged$classification))
    
    # Store results
    purity_performance[[pc]] <- list(
      n_total = nrow(pc_results),
      n_judged = nrow(pc_judged),
      judgment_rate = nrow(pc_judged) / max(1, nrow(pc_results)),
      mean_score = mean(pc_judged$reo_score),
      median_score = median(pc_judged$reo_score),
      classifications = table(pc_judged$classification)
    )
  } else {
    purity_performance[[pc]] <- list(
      n_total = nrow(pc_results),
      n_judged = 0,
      judgment_rate = 0
    )
  }
}

# Check monotonicity
if (length(purity_performance) >= 2 && 
    purity_performance[["High"]]$n_judged > 0 &&
    purity_performance[["Low"]]$n_judged > 0) {
  
  score_diff <- purity_performance[["High"]]$median_score - 
    purity_performance[["Low"]]$median_score
  
  cat(sprintf("\n  Score monotonicity (High - Low): %.3f\n", score_diff))
  if (abs(score_diff) < 0.1) {
    cat("  ✓ Scores stable across purity levels\n")
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
    n_r_test_paired = nrow(r_test_paired),
    n_r_test_unpaired = nrow(r_test_unpaired),
    n_samples_analyzed = nrow(r_test_results),
    n_high_purity = sum(r_test_results$purity_class == "High"),
    n_judged = sum(r_test_results$judgment_possible),
    timestamp = Sys.time()
  ),
  
  # Metadata
  step = "Step 6 Complete"
)

saveRDS(step6_data, paste0(paths$processed, "reo_step6_validation.rds"))
cat("  Step 6 data saved to: reo_step6_validation.rds\n")

# Export key results
validation_summary <- data.frame(
  metric = character(),
  value = character(),
  stringsAsFactors = FALSE
)

# Add main results
if (!is.null(stratum_results)) {
  for (stratum in names(stratum_results)) {
    # Check if insufficient field exists and is FALSE
    if (is.null(stratum_results[[stratum]]$insufficient) || 
        !stratum_results[[stratum]]$insufficient) {
      validation_summary <- rbind(validation_summary,
                                  data.frame(
                                    metric = paste0(stratum, "_rate"),
                                    value = sprintf("%.1f%% [%.1f-%.1f]", 
                                                    stratum_results[[stratum]]$rate * 100,
                                                    stratum_results[[stratum]]$ci_lower * 100,
                                                    stratum_results[[stratum]]$ci_upper * 100)
                                  )
      )
    }
  }
}

# Add diagnostic results
validation_summary <- rbind(validation_summary,
                            data.frame(
                              metric = c("R0_FPR", "R1_FNR", "CI_overlap"),
                              value = c(
                                sprintf("%.1f%%", r0_fpr * 100),
                                sprintf("%.1f%%", r1_fnr * 100),
                                ifelse(ci_overlap, "Yes", "No")
                              )
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
cat(sprintf("Samples with valid judgment: %d\n", 
            sum(r_test_results$judgment_possible)))
cat(sprintf("Timestamp: %s\n", format(step6_data$summary_stats$timestamp)))

# Performance summary
if (!is.null(stratum_results) && length(stratum_results) > 0) {
  cat("\nMain Analysis Results (High Purity):\n")
  for (stratum in names(stratum_results)) {
    sr <- stratum_results[[stratum]]
    # Check if insufficient field exists
    is_insufficient <- !is.null(sr$insufficient) && sr$insufficient
    if (!is_insufficient && !is.null(sr$rate)) {
      cat(sprintf("  %s: R1-like rate = %.1f%% [%.1f%%, %.1f%%]\n",
                  stratum, sr$rate*100, sr$ci_lower*100, sr$ci_upper*100))
    } else {
      cat(sprintf("  %s: n=%d (insufficient for CI)\n", stratum, sr$n))
    }
  }
}

# Diagnostic summary
cat("\nDiagnostic Checks:\n")
cat(sprintf("  Classification symmetry: %s\n", 
            ifelse(ci_overlap, "Balanced", "Imbalanced")))
cat(sprintf("  External specificity (BRAF): %s\n",
            ifelse(!is.null(braf_results), "Tested", "Not tested")))
cat(sprintf("  Jackknife stability: %s\n",
            ifelse(!is.null(jackknife_results), "Assessed", "Not assessed")))

# Final assessment
cat("\n✓ Step 6 validation complete. REO panel performance validated.\n")