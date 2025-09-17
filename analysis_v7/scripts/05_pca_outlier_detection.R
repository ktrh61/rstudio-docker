# 05_pca_outlier_detection.R - Stage 1 PCA Outlier Detection (v7.9)
# Purpose: Detect technical outliers using CDM and update case_master
# Method: Cross-Data Matrix with Permutation Analysis and adaptive thresholds
# Version: v7.9 - Case-master centric approach with outlier flags
#          Outputs thyr_case_master_stage1_filtered with has_outlier_tumor/normal flags
# Date: 2025-01-20

source("analysis_v7/setup.R")

cat("\n=== Stage 1: PCA Outlier Detection (v7.9) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Analysis: 8 groups (R0/R1/B0/B1 Ã— tumor/normal)\n")
cat("Sample types: Primary Tumor and Solid Tissue Normal only\n")
cat("Sample filtering: _merged samples only\n")
cat("Output: case_master with outlier flags (0/1)\n")

# Load packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(Rcpp)
  library(RcppArmadillo)
  library(parallel)
  library(dplyr)
})

# ============================================================================
# Configuration
# ============================================================================

CONFIG <- list(
  # Permutation Analysis parameters
  PA_R = 200,           # Permutation replicates
  PA_ALPHA = 0.05,      # Significance level for PA
  PA_L = 8,             # Maximum components to evaluate
  PA_SEED = 123,        # Base seed for reproducibility
  
  # Adaptive outlier detection parameters
  MIN_OUTLIERS_PCT = 0.00,    # Minimum outlier percentage (0% = no forced outliers)
  MAX_OUTLIERS_PCT = 0.20,    # Maximum outlier percentage (20% cap)
  IQR_MULTIPLIER = 3.0,       # IQR multiplier for outlier detection
  MAD_MULTIPLIER = 5.0,       # MAD multiplier for outlier detection
  USE_ADAPTIVE = TRUE,        # Use adaptive thresholds
  USE_OD = TRUE,              # Use Orthogonal Distance
  
  # Processing parameters
  PRIOR_COUNT = 0.5,          # Pseudocount for logCPM
  MC_CORES = min(3L, max(1L, parallel::detectCores() - 1L)),  # Conservative setting
  VERBOSE_CDM = FALSE         # Disable verbose in CDM for parallel processing
)

cat("\nConfiguration:\n")
cat("  Parallel cores:", CONFIG$MC_CORES, "\n")
cat("  Adaptive thresholds:", CONFIG$USE_ADAPTIVE, "\n")
cat("  Max outliers per group:", sprintf("%.0f%%", CONFIG$MAX_OUTLIERS_PCT * 100), "\n")

# ============================================================================
# Thread control utilities
# ============================================================================

# Source BLAS thread control utility
if (file.exists("utils/with_openblas_threads.R")) {
  source("utils/with_openblas_threads.R")
  cat("  BLAS control: with_openblas_threads loaded\n")
} else {
  # Fallback: manual thread control
  cat("  BLAS control: Manual (with_openblas_threads.R not found)\n")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  } else {
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      BLIS_NUM_THREADS = "1"
    )
  }
}

# Set RNG for reproducibility
RNGkind("L'Ecuyer-CMRG")

# ============================================================================
# Compile CDM implementation
# ============================================================================

cat("\nCompiling CDM implementation...\n")
sourceCpp("utils/CDM_fast3_arma_enhanced.cpp")
cat("  CDM compiled successfully\n")

# ============================================================================
# Load data
# ============================================================================

cat("\n--- Loading data ---\n")

# Load processed SummarizedExperiment
if (exists("thyr_se_strand2_nonzero")) {
  cat("Using existing thyr_se_strand2_nonzero\n")
  se <- thyr_se_strand2_nonzero
} else {
  se_path <- paste0(paths$processed, "thyr_se_strand2_nonzero.rds")
  if (!file.exists(se_path)) {
    stop("thyr_se_strand2_nonzero.rds not found. Please run 03_prepare_counts.R first.")
  }
  cat("Loading SE from file...\n")
  se <- readRDS(se_path)
}
cat("SE dimensions:", format(nrow(se), big.mark=","), "genes Ã—", 
    format(ncol(se), big.mark=","), "samples\n")

# Load case master for group definitions
case_master_path <- paste0(paths$processed, "thyr_case_master.rds")
if (!file.exists(case_master_path)) {
  stop("thyr_case_master.rds not found. Please run 04_clinical_integration.R first.")
}

thyr_case_master <- readRDS(case_master_path)
cat("Case master loaded:", nrow(thyr_case_master), "cases\n")

# Filter to main analysis groups only
thyr_case_master_stage1_filtered <- thyr_case_master %>%
  filter(group %in% c("R0", "R1", "B0", "B1"))

cat("Main analysis groups:", nrow(thyr_case_master_stage1_filtered), "cases\n")

# Get sample metadata
sample_metadata <- as.data.frame(colData(se))

# ============================================================================
# Define analysis groups
# ============================================================================

cat("\n--- Defining analysis groups ---\n")
cat("Using only Primary Tumor and Solid Tissue Normal samples\n")
cat("Filtering to _merged samples only for proper pairing\n")

# Function to extract samples for a group
extract_group_samples <- function(group_name, sample_type_filter, case_master_df, sample_meta) {
  # Get cases for this group
  group_cases <- case_master_df$case_submitter_id[case_master_df$group == group_name]
  
  if (length(group_cases) == 0) {
    return(character(0))
  }
  
  # Find samples for these cases
  sample_ids <- sample_meta$sample_submitter_id[
    sample_meta$case_submitter_id %in% group_cases &
      sample_meta$sample_type == sample_type_filter &
      grepl("_merged", sample_meta$sample_submitter_id)
  ]
  
  return(sample_ids)
}

# Create groups as lists of sample_submitter_ids
groups_ids <- list(
  R0_tumor = extract_group_samples("R0", "Primary Tumor", thyr_case_master_stage1_filtered, sample_metadata),
  R0_normal = extract_group_samples("R0", "Solid Tissue Normal", thyr_case_master_stage1_filtered, sample_metadata),
  R1_tumor = extract_group_samples("R1", "Primary Tumor", thyr_case_master_stage1_filtered, sample_metadata),
  R1_normal = extract_group_samples("R1", "Solid Tissue Normal", thyr_case_master_stage1_filtered, sample_metadata),
  B0_tumor = extract_group_samples("B0", "Primary Tumor", thyr_case_master_stage1_filtered, sample_metadata),
  B0_normal = extract_group_samples("B0", "Solid Tissue Normal", thyr_case_master_stage1_filtered, sample_metadata),
  B1_tumor = extract_group_samples("B1", "Primary Tumor", thyr_case_master_stage1_filtered, sample_metadata),
  B1_normal = extract_group_samples("B1", "Solid Tissue Normal", thyr_case_master_stage1_filtered, sample_metadata)
)

# Remove groups with insufficient samples (less than 3)
groups_ids <- groups_ids[sapply(groups_ids, length) >= 3]

# Report group sizes
cat("\nGroups to process:\n")
for (g in names(groups_ids)) {
  cat(sprintf("  %s: %d samples\n", g, length(groups_ids[[g]])))
}
cat("Total groups:", length(groups_ids), "\n")

# ============================================================================
# Unified gene filtering
# ============================================================================

cat("\n--- Applying unified gene filtering ---\n")

# Get all sample IDs from analysis groups
all_analysis_sample_ids <- unlist(groups_ids, use.names = FALSE)

# Extract counts for these samples
all_analysis_counts <- assay(se)[, all_analysis_sample_ids, drop = FALSE]

# Handle any NA values
if (any(is.na(all_analysis_counts))) {
  cat("  Warning: Found", sum(is.na(all_analysis_counts)), "NA values. Replacing with 0.\n")
  all_analysis_counts[is.na(all_analysis_counts)] <- 0
}

# Create group labels for filterByExpr
group_labels <- factor(rep(names(groups_ids), sapply(groups_ids, length)))

# Apply filterByExpr with group consideration
keep_genes <- filterByExpr(all_analysis_counts, group = group_labels)

cat("Gene filtering results:\n")
cat("  Starting genes:", nrow(se), "\n")
cat("  After filterByExpr:", sum(keep_genes), "\n")
cat("  Retention rate:", sprintf("%.1f%%\n", sum(keep_genes) / nrow(se) * 100))

# ============================================================================
# Helper functions
# ============================================================================

# Deterministic seed generation from group name
seed_from_name <- function(name, base = 123L) {
  as.integer(max(1L, abs(sum(utf8ToInt(name)) * base) %% 2147483647L))
}

# Permutation Analysis function
pa_select_K <- function(X, R = 200, alpha = 0.05, L = 8, seed = 123) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  
  X_t <- t(X)  # samples Ã— genes
  n <- nrow(X_t)
  p <- ncol(X_t)
  L <- min(L, n - 1L, p)
  
  if (L < 1L) {
    return(list(K = 1L, real = numeric(0), threshold = numeric(0)))
  }
  
  # Column centering
  mu <- colMeans(X_t)
  center <- function(M) sweep(M, 2, mu, "-")
  
  # Function to get top L eigenvalues
  topL <- function(M) {
    sv <- svd(center(M), nu = 0, nv = 0)$d^2
    sv[seq_len(min(L, length(sv)))]
  }
  
  # Real eigenvalues
  real_sv <- topL(X_t)
  actual_L <- length(real_sv)
  
  # Permutation null distribution
  null_sv <- replicate(R, {
    Xp <- apply(X_t, 2, base::sample)
    topL(Xp)
  })
  
  if (is.null(dim(null_sv))) {
    null_sv <- matrix(null_sv, nrow = 1L)
  }
  
  # Calculate threshold
  thr <- apply(null_sv, 1, quantile, probs = 1 - alpha, na.rm = TRUE)
  K <- max(1L, sum(real_sv > thr * 1.00001))
  
  list(K = K, real = real_sv, threshold = thr)
}

# Adaptive outlier detection function
cdm_outlier_detection <- function(cdm_result, X, K, config = CONFIG) {
  # Ensure X is samples Ã— genes
  if (ncol(X) == nrow(cdm_result$vectors)) {
    Xuse <- X
  } else {
    Xuse <- t(X)
  }
  n <- nrow(Xuse)
  
  # Define outlier constraints
  min_outliers <- ceiling(n * config$MIN_OUTLIERS_PCT)
  max_outliers <- floor(n * config$MAX_OUTLIERS_PCT)
  
  # Center data
  mu <- colMeans(Xuse)
  Xc <- sweep(Xuse, 2, mu, "-")
  
  # Extract eigenvectors and compute scores
  K <- max(1L, min(K, n - 1L, ncol(cdm_result$vectors)))
  V <- qr.Q(qr(cdm_result$vectors[, 1:K, drop = FALSE]))
  S <- Xc %*% V
  
  # Calculate Score Distance (SD) with robust normalization
  Z <- apply(S, 2, function(x) {
    m <- median(x)
    s <- mad(x, constant = 1.4826)
    if (!is.finite(s) || s <= 1e-8) {
      s <- IQR(x) / 1.349
    }
    if (!is.finite(s) || s <= 1e-8) {
      return(rep(NA_real_, length(x)))
    }
    (x - m) / s
  })
  
  # Filter out PCs with zero variance
  keep_pc <- which(colSums(!is.na(Z)) == n)
  if (length(keep_pc) == 0L) {
    stop("All PCs dropped due to zero variance")
  }
  
  Z <- Z[, keep_pc, drop = FALSE]
  SD <- sqrt(rowSums(Z^2))
  
  # Calculate adaptive threshold for SD
  if (config$USE_ADAPTIVE) {
    q1 <- quantile(SD, 0.25, na.rm = TRUE)
    q3 <- quantile(SD, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    
    thr_SD_iqr <- if (is.finite(iqr) && iqr > 1e-8) {
      q3 + config$IQR_MULTIPLIER * iqr
    } else {
      Inf
    }
    
    med_SD <- median(SD, na.rm = TRUE)
    mad_SD <- mad(SD, constant = 1.4826, na.rm = TRUE)
    
    thr_SD_mad <- if (is.finite(mad_SD) && mad_SD > 1e-8) {
      med_SD + config$MAD_MULTIPLIER * mad_SD
    } else {
      Inf
    }
    
    thr_SD <- max(thr_SD_iqr, thr_SD_mad, na.rm = TRUE)
    
    if (!is.finite(thr_SD)) {
      thr_SD <- quantile(SD, 0.975, na.rm = TRUE)
    }
  } else {
    thr_SD <- quantile(SD, 0.975, na.rm = TRUE)
  }
  
  # Apply max outlier constraint to SD threshold
  if (max_outliers > 0 && sum(SD > thr_SD) > max_outliers) {
    sorted_SD <- sort(SD, decreasing = TRUE)
    thr_SD <- sorted_SD[max_outliers] + 1e-10
  }
  
  hit_SD <- SD > thr_SD
  
  # Calculate Orthogonal Distance (OD) if enabled
  hit_OD <- rep(FALSE, n)
  thr_OD <- NA_real_
  OD <- NULL
  
  if (config$USE_OD) {
    OD_raw <- sqrt(pmax(rowSums(Xc^2) - rowSums(S^2), 0))
    OD <- log1p(OD_raw)
    
    if (config$USE_ADAPTIVE) {
      q1_OD <- quantile(OD, 0.25, na.rm = TRUE)
      q3_OD <- quantile(OD, 0.75, na.rm = TRUE)
      iqr_OD <- q3_OD - q1_OD
      
      thr_OD_iqr <- if (is.finite(iqr_OD) && iqr_OD > 1e-8) {
        q3_OD + config$IQR_MULTIPLIER * iqr_OD
      } else {
        Inf
      }
      
      med_OD <- median(OD, na.rm = TRUE)
      mad_OD <- mad(OD, constant = 1.4826, na.rm = TRUE)
      
      thr_OD_mad <- if (is.finite(mad_OD) && mad_OD > 1e-8) {
        med_OD + config$MAD_MULTIPLIER * mad_OD
      } else {
        Inf
      }
      
      thr_OD <- max(thr_OD_iqr, thr_OD_mad, na.rm = TRUE)
      
      if (!is.finite(thr_OD)) {
        thr_OD <- quantile(OD, 0.975, na.rm = TRUE)
      }
    } else {
      thr_OD <- quantile(OD, 0.975, na.rm = TRUE)
    }
    
    # Apply max outlier constraint to OD threshold
    if (max_outliers > 0 && sum(OD > thr_OD) > max_outliers) {
      sorted_OD <- sort(OD, decreasing = TRUE)
      thr_OD <- sorted_OD[max_outliers] + 1e-10
    }
    
    hit_OD <- OD > thr_OD
  }
  
  # Combine outliers
  outliers <- hit_SD | hit_OD
  
  # Final enforcement of maximum outlier constraint
  if (max_outliers > 0L && sum(outliers) > max_outliers) {
    # Rank by combined deviation
    mad_SD_global <- mad(SD, constant = 1.4826, na.rm = TRUE)
    zSD <- if (mad_SD_global > 1e-8) {
      (SD - thr_SD) / mad_SD_global
    } else {
      SD - thr_SD
    }
    
    if (config$USE_OD && !is.null(OD)) {
      mad_OD_global <- mad(OD, constant = 1.4826, na.rm = TRUE)
      zOD <- if (mad_OD_global > 1e-8) {
        (OD - thr_OD) / mad_OD_global
      } else {
        OD - thr_OD
      }
    } else {
      zOD <- rep(-Inf, length(zSD))
    }
    
    fused <- pmax(zSD, zOD, na.rm = TRUE)
    keep_idx <- order(fused, decreasing = TRUE)[seq_len(max_outliers)]
    
    new_outliers <- rep(FALSE, length(outliers))
    new_outliers[keep_idx] <- TRUE
    outliers <- new_outliers
  }
  
  # Return results
  list(
    outliers = outliers,
    n_outliers = sum(outliers),
    K = ncol(V),
    SD = SD,
    OD = OD,
    thr_SD = thr_SD,
    thr_OD = thr_OD,
    hit_SD = sum(hit_SD),
    hit_OD = if(config$USE_OD) sum(hit_OD) else 0
  )
}

# Process one group function
process_group <- function(group_name, group_sample_ids, se, keep_genes, config) {
  # Ensure single thread in worker
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
  
  start_time <- Sys.time()
  
  # Extract counts using sample_ids directly
  counts <- assay(se)[keep_genes, group_sample_ids, drop = FALSE]
  if (any(is.na(counts))) {
    counts[is.na(counts)] <- 0
  }
  
  # logCPM transformation
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge)
  logcpm <- cpm(dge, log = TRUE, prior.count = config$PRIOR_COUNT)
  
  # Permutation Analysis for K selection
  pa_seed <- seed_from_name(group_name, config$PA_SEED)
  pa_result <- pa_select_K(
    logcpm,
    R = config$PA_R,
    alpha = config$PA_ALPHA,
    L = config$PA_L,
    seed = pa_seed
  )
  
  # CDM computation
  cdm_result <- CDM_fast3_arma(logcpm, verbose = config$VERBOSE_CDM)
  
  # Outlier detection
  outlier_result <- cdm_outlier_detection(
    cdm_result, 
    logcpm, 
    K = pa_result$K, 
    config = config
  )
  
  # Get outlier sample_submitter_ids directly
  outlier_ids <- if(any(outlier_result$outliers)) {
    group_sample_ids[outlier_result$outliers]
  } else {
    character(0)
  }
  
  # Processing time
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Return results
  list(
    group = group_name,
    n_samples = length(group_sample_ids),
    n_outliers = outlier_result$n_outliers,
    outlier_ids = outlier_ids,
    K = pa_result$K,
    elapsed_time = elapsed,
    diagnostics = list(
      hit_SD = outlier_result$hit_SD,
      hit_OD = outlier_result$hit_OD,
      thr_SD = outlier_result$thr_SD,
      thr_OD = outlier_result$thr_OD
    )
  )
}

# ============================================================================
# Main processing
# ============================================================================

cat("\n--- Processing groups ---\n")

overall_start <- Sys.time()

# Process groups in parallel
if (CONFIG$MC_CORES > 1L && .Platform$OS.type == "unix") {
  cat("Starting parallel processing with", CONFIG$MC_CORES, "cores...\n")
  
  # Set seed for reproducibility
  RNGkind("L'Ecuyer-CMRG")
  set.seed(CONFIG$PA_SEED)
  
  # Create wrapper with error handling
  process_wrapper <- function(gn) {
    tryCatch({
      if (exists("with_openblas_threads") && is.function(with_openblas_threads)) {
        with_openblas_threads(
          threads = 1,
          expr = process_group(gn, groups_ids[[gn]], se, keep_genes, CONFIG)
        )
      } else {
        process_group(gn, groups_ids[[gn]], se, keep_genes, CONFIG)
      }
    }, error = function(e) {
      list(
        group = gn,
        n_samples = length(groups_ids[[gn]]),
        n_outliers = 0,
        outlier_ids = character(0),
        K = NA,
        elapsed_time = NA,
        error = e$message
      )
    })
  }
  
  # Run parallel processing
  results <- mclapply(
    names(groups_ids),
    process_wrapper,
    mc.cores = CONFIG$MC_CORES,
    mc.preschedule = FALSE
  )
  names(results) <- names(groups_ids)
  
} else {
  # Sequential processing
  cat("Starting sequential processing...\n")
  results <- list()
  
  for (gn in names(groups_ids)) {
    cat(sprintf("Processing %s...\n", gn))
    
    tryCatch({
      if (exists("with_openblas_threads") && is.function(with_openblas_threads)) {
        results[[gn]] <- with_openblas_threads(
          threads = 1,
          expr = process_group(gn, groups_ids[[gn]], se, keep_genes, CONFIG)
        )
      } else {
        results[[gn]] <- process_group(gn, groups_ids[[gn]], se, keep_genes, CONFIG)
      }
    }, error = function(e) {
      cat("  Error:", e$message, "\n")
      results[[gn]] <- list(
        group = gn,
        n_samples = length(groups_ids[[gn]]),
        n_outliers = 0,
        outlier_ids = character(0),
        K = NA,
        elapsed_time = NA,
        error = e$message
      )
    })
  }
}

overall_end <- Sys.time()
total_time <- difftime(overall_end, overall_start, units = "secs")

cat(sprintf("\nAll groups processed in %.1f seconds\n", as.numeric(total_time)))

# ============================================================================
# Create outlier mapping by sample_id
# ============================================================================

cat("\n--- Creating outlier mapping ---\n")

# Collect all outlier sample IDs
outlier_samples <- unlist(lapply(results, function(r) r$outlier_ids))
outlier_samples <- unique(outlier_samples)

cat("Total unique outlier samples:", length(outlier_samples), "\n")

# ============================================================================
# Update case_master with outlier flags
# ============================================================================

cat("\n--- Updating case_master with outlier flags ---\n")

# Initialize outlier flag columns right after has_pair
thyr_case_master_stage1_filtered$has_outlier_tumor <- 0
thyr_case_master_stage1_filtered$has_outlier_normal <- 0

# Reorder columns to place outlier flags after has_pair
col_order <- names(thyr_case_master_stage1_filtered)
has_pair_idx <- which(col_order == "has_pair")
if (length(has_pair_idx) > 0) {
  # Insert outlier columns after has_pair
  new_order <- c(
    col_order[1:has_pair_idx],
    "has_outlier_tumor", "has_outlier_normal",
    col_order[(has_pair_idx + 1):(length(col_order) - 2)]  # Exclude the two outlier columns at the end
  )
  thyr_case_master_stage1_filtered <- thyr_case_master_stage1_filtered[, ..new_order]
}

# For each case, check if its samples are outliers
for (i in seq_len(nrow(thyr_case_master_stage1_filtered))) {
  case_id <- thyr_case_master_stage1_filtered$case_submitter_id[i]
  
  # Find tumor samples for this case
  tumor_samples <- sample_metadata$sample_submitter_id[
    sample_metadata$case_submitter_id == case_id &
      sample_metadata$sample_type == "Primary Tumor" &
      grepl("_merged", sample_metadata$sample_submitter_id)
  ]
  
  # Find normal samples for this case
  normal_samples <- sample_metadata$sample_submitter_id[
    sample_metadata$case_submitter_id == case_id &
      sample_metadata$sample_type == "Solid Tissue Normal" &
      grepl("_merged", sample_metadata$sample_submitter_id)
  ]
  
  # Check if any tumor sample is an outlier
  if (length(tumor_samples) > 0 && any(tumor_samples %in% outlier_samples)) {
    thyr_case_master_stage1_filtered$has_outlier_tumor[i] <- 1
  }
  
  # Check if any normal sample is an outlier
  if (length(normal_samples) > 0 && any(normal_samples %in% outlier_samples)) {
    thyr_case_master_stage1_filtered$has_outlier_normal[i] <- 1
  }
}

# Summary of outliers by group
cat("\n--- Outlier summary by group ---\n")
outlier_summary <- thyr_case_master_stage1_filtered %>%
  group_by(group) %>%
  summarise(
    n_cases = n(),
    n_tumor_outliers = sum(has_outlier_tumor),
    n_normal_outliers = sum(has_outlier_normal),
    n_both_outliers = sum(has_outlier_tumor == 1 & has_outlier_normal == 1),
    n_clean = sum(has_outlier_tumor == 0 & has_outlier_normal == 0),
    .groups = "drop"
  )

print(outlier_summary)

# Overall statistics
total_cases <- nrow(thyr_case_master_stage1_filtered)
clean_cases <- sum(thyr_case_master_stage1_filtered$has_outlier_tumor == 0 & 
                     thyr_case_master_stage1_filtered$has_outlier_normal == 0)

cat(sprintf("\nOverall: %d/%d cases clean (%.1f%%)\n",
            clean_cases, total_cases, clean_cases/total_cases*100))

# ============================================================================
# Processing summary
# ============================================================================

cat("\n=== Processing Summary ===\n")
cat(sprintf("%-12s %6s %4s %8s %6s %6s %8s\n", 
            "Group", "N_Samp", "K", "Outliers", "SD_Out", "OD_Out", "Time(s)"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (g in names(results)) {
  r <- results[[g]]
  if (!is.null(r$error)) {
    cat(sprintf("%-12s %6d %4s %8s %6s %6s %8s ERROR\n",
                g, r$n_samples, "NA", "NA", "NA", "NA", "NA"))
  } else {
    cat(sprintf("%-12s %6d %4d %8d %6d %6d %8.1f\n",
                g, r$n_samples, r$K, r$n_outliers,
                r$diagnostics$hit_SD, r$diagnostics$hit_OD,
                r$elapsed_time))
  }
}

cat(paste(rep("-", 60), collapse = ""), "\n")

# ============================================================================
# Save results
# ============================================================================

cat("\n--- Saving results ---\n")

# Save the updated case_master
output_file <- paste0(paths$processed, "thyr_case_master_stage1_filtered.rds")
saveRDS(thyr_case_master_stage1_filtered, output_file)
cat("  Updated case_master saved:", basename(output_file), "\n")
cat(sprintf("    Contains: %d cases (R0/R1/B0/B1 only)\n", nrow(thyr_case_master_stage1_filtered)))
cat("    New columns: has_outlier_tumor, has_outlier_normal (0/1 flags)\n")

# Save processing details for diagnostic purposes
processing_info <- list(
  date = Sys.Date(),
  config = CONFIG,
  results = results,
  outlier_samples = outlier_samples,
  processing_time = total_time,
  summary = outlier_summary
)
saveRDS(processing_info, paste0(paths$output, "stage1_pca_info.rds"))
cat("  Processing info saved: stage1_pca_info.rds\n")

# Export summary as CSV for easy review
write.csv(outlier_summary, 
          paste0(paths$output, "stage1_outlier_summary.csv"),
          row.names = FALSE)
cat("  Summary CSV saved: stage1_outlier_summary.csv\n")

# ============================================================================
# Final report
# ============================================================================

cat("\n=== Stage 1 Complete ===\n")
cat("Configuration:\n")
cat("  Method: Adaptive thresholds (IQR/MAD)\n")
cat("  IQR multiplier:", CONFIG$IQR_MULTIPLIER, "\n")
cat("  MAD multiplier:", CONFIG$MAD_MULTIPLIER, "\n")
cat("  Max outliers:", sprintf("%.0f%%", CONFIG$MAX_OUTLIERS_PCT * 100), "\n")
cat("\nProcessing:\n")
cat("  Groups analyzed:", length(groups_ids), "(R0/R1/B0/B1 Ã— tumor/normal)\n")
cat("  Total time:", sprintf("%.1f seconds\n", as.numeric(total_time)))
cat("\nResults:\n")
cat("  Total outlier samples:", length(outlier_samples), "\n")
cat("  Clean cases:", clean_cases, "/", total_cases, 
    sprintf("(%.1f%%)\n", clean_cases/total_cases*100))
cat("\nOutputs:\n")
cat("  Main: thyr_case_master_stage1_filtered.rds\n")
cat("  Info: stage1_pca_info.rds, stage1_outlier_summary.csv\n")
cat("\nNext: Use thyr_case_master_stage1_filtered for downstream analysis\n")
cat("      Filter cases with: has_outlier_tumor == 0 & has_outlier_normal == 0\n")

# Restore SE to original variable name
thyr_se_strand2_nonzero <- se

# Clean up
rm(list = setdiff(ls(), c("paths", "thyr_case_master", "thyr_case_master_stage1_filtered", "thyr_se_strand2_nonzero")))
gc()