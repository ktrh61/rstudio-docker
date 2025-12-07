# 05_pca_outlier_detection.R - PCA Outlier Detection (v7.12)
# Purpose: Detect technical outliers using CDM-PCA and update case_master
#          - Protein coding genes only (consistent with downstream analysis)
#          - Per-group filterByExpr (group-specific expression filtering)
#          - MUREN normalization (consistent with 08_deges_normalization.R)
# Method: Cross-Data Matrix with Permutation Analysis and adaptive thresholds
# Version: v7.12 - Protein coding, per-group filterByExpr, MUREN normalization
#                  SD vs OD diagnostic plots
# Date: 2025-12-07

source("analysis_v7/setup.R")

cat("\n=== PCA Outlier Detection (v7.12) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Analysis: 12 groups (R0/R_Low/R_High/R1/B0/B1 × tumor/normal)\n")
cat("Note: B_Low/B_High excluded (not used in validation)\n")
cat("Gene filtering: Protein coding only, per-group filterByExpr\n")
cat("Normalization: MUREN (consistent with downstream)\n")
cat("Sample filtering: _merged samples only, paired cases only\n")
cat("Output: case_master with outlier flags (0/1/NA)\n")

# Load packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(Rcpp)
  library(RcppArmadillo)
  library(parallel)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})

# ============================================================================
# Configuration
# ============================================================================

CONFIG <- list(
  # Permutation Analysis parameters
  PA_R = 200,           # Permutation replicates
  PA_ALPHA = 0.05,      # Significance level for PA
  PA_L = 8,             # Maximum components to evaluate
  PA_SEED = 1986,       # Base seed for reproducibility
  
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
  VERBOSE_CDM = FALSE,        # Disable verbose in CDM for parallel processing
  
  # MUREN parameters
  MUREN_METHOD = "lts",       # LTS robust regression
  MUREN_WORKERS = 3           # Workers for MUREN
)

cat("\nConfiguration:\n")
cat("  Parallel cores:", CONFIG$MC_CORES, "\n")
cat("  Adaptive thresholds:", CONFIG$USE_ADAPTIVE, "\n")
cat("  Max outliers per group:", sprintf("%.0f%%", CONFIG$MAX_OUTLIERS_PCT * 100), "\n")
cat("  MUREN method:", CONFIG$MUREN_METHOD, "\n")
cat("  Seed:", CONFIG$PA_SEED, "\n")

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

# Load utility functions (including MUREN)
source("utils/utils_improved.R")
source("utils/norm_improved.R")
cat("  Utilities loaded: utils_improved.R, norm_improved.R\n")

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
cat("SE dimensions:", format(nrow(se), big.mark=","), "genes ×", 
    format(ncol(se), big.mark=","), "samples\n")

# Load case master for group definitions
case_master_path <- paste0(paths$processed, "thyr_case_master.rds")
if (!file.exists(case_master_path)) {
  stop("thyr_case_master.rds not found. Please run 04_clinical_integration.R first.")
}

thyr_case_master <- readRDS(case_master_path)
cat("Case master loaded:", nrow(thyr_case_master), "cases\n")

# Filter to all analysis groups (paired samples only for PCA QC)
# Note: All groups included, but PCA QC only runs on paired samples
thyr_case_master_stage1_filtered <- thyr_case_master %>%
  filter(group %in% c("R0", "R_Low", "R_High", "R1", "B0", "B_Low", "B_High", "B1"))

cat("Analysis groups:", nrow(thyr_case_master_stage1_filtered), "cases\n")
cat("  Paired:", sum(thyr_case_master_stage1_filtered$has_pair), "cases (PCA QC target)\n")
cat("  Unpaired:", sum(!thyr_case_master_stage1_filtered$has_pair), "cases (no PCA QC)\n")

# Get sample metadata
sample_metadata <- as.data.frame(colData(se))

# ============================================================================
# Protein coding gene identification
# ============================================================================

cat("\n--- Identifying protein coding genes ---\n")

gene_info <- as.data.frame(rowData(se))
is_protein_coding <- gene_info$gene_type == "protein_coding"

cat(sprintf("Protein coding genes: %d / %d (%.1f%%)\n", 
            sum(is_protein_coding), nrow(se), 
            sum(is_protein_coding) / nrow(se) * 100))

# ============================================================================
# Define analysis groups
# ============================================================================

cat("\n--- Defining analysis groups ---\n")
cat("Using only Primary Tumor and Solid Tissue Normal samples\n")
cat("Filtering to _merged samples only for proper pairing\n")
cat("PCA QC: paired samples only\n")

# Function to extract samples for a group (paired cases only for PCA QC)
extract_group_samples <- function(group_name, sample_type_filter, case_master_df, sample_meta) {
  # Get paired cases for this group (PCA QC requires paired samples)
  group_cases <- case_master_df$case_submitter_id[
    case_master_df$group == group_name & case_master_df$has_pair
  ]
  
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

# Create groups as lists of sample_submitter_ids (paired only)
# Note: B_Low/B_High excluded from PCA QC (not used in validation)
groups_ids <- list(
  # RET groups (all for validation)
  R0_tumor = extract_group_samples("R0", "Primary Tumor", thyr_case_master_stage1_filtered, sample_metadata),
  R0_normal = extract_group_samples("R0", "Solid Tissue Normal", thyr_case_master_stage1_filtered, sample_metadata),
  R_Low_tumor = extract_group_samples("R_Low", "Primary Tumor", thyr_case_master_stage1_filtered, sample_metadata),
  R_Low_normal = extract_group_samples("R_Low", "Solid Tissue Normal", thyr_case_master_stage1_filtered, sample_metadata),
  R_High_tumor = extract_group_samples("R_High", "Primary Tumor", thyr_case_master_stage1_filtered, sample_metadata),
  R_High_normal = extract_group_samples("R_High", "Solid Tissue Normal", thyr_case_master_stage1_filtered, sample_metadata),
  R1_tumor = extract_group_samples("R1", "Primary Tumor", thyr_case_master_stage1_filtered, sample_metadata),
  R1_normal = extract_group_samples("R1", "Solid Tissue Normal", thyr_case_master_stage1_filtered, sample_metadata),
  # BRAF groups (B0/B1 only for driver generalization validation)
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
# Helper functions
# ============================================================================

# Deterministic seed generation from group name
seed_from_name <- function(name, base = 1986L) {
  as.integer(max(1L, abs(sum(utf8ToInt(name)) * base) %% 2147483647L))
}

# Permutation Analysis function
pa_select_K <- function(X, R = 200, alpha = 0.05, L = 8, seed = 1986) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  
  X_t <- t(X)  # samples × genes
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
  # Ensure X is samples × genes
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
process_group <- function(group_name, group_sample_ids, se, is_protein_coding, config) {
  # Ensure single thread in worker
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
  
  start_time <- Sys.time()
  
  # Extract counts for protein coding genes only
  counts_pc <- assay(se)[is_protein_coding, group_sample_ids, drop = FALSE]
  if (any(is.na(counts_pc))) {
    counts_pc[is.na(counts_pc)] <- 0
  }
  
  # Per-group filterByExpr
  dge_initial <- DGEList(counts = counts_pc)
  keep <- filterByExpr(dge_initial)
  counts_filtered <- counts_pc[keep, , drop = FALSE]
  n_genes_used <- nrow(counts_filtered)
  
  # MUREN normalization
  muren_coeff <- muren_norm(
    counts_filtered,
    refs = "saturated",
    pairwise_method = config$MUREN_METHOD,
    single_param = TRUE,
    res_return = "scaling_coeff",
    workers = config$MUREN_WORKERS
  )
  
  # Create DGEList with MUREN normalization factors
  dge <- DGEList(counts = counts_filtered)
  dge$samples$norm.factors <- muren_coeff / mean(muren_coeff)
  
  # logCPM transformation
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
  
  # Return results (including SD/OD for plotting)
  list(
    group = group_name,
    n_samples = length(group_sample_ids),
    n_genes = n_genes_used,
    n_outliers = outlier_result$n_outliers,
    outlier_ids = outlier_ids,
    outlier_flags = outlier_result$outliers,
    K = pa_result$K,
    elapsed_time = elapsed,
    sample_ids = group_sample_ids,
    SD = outlier_result$SD,
    OD = outlier_result$OD,
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
          expr = process_group(gn, groups_ids[[gn]], se, is_protein_coding, CONFIG)
        )
      } else {
        process_group(gn, groups_ids[[gn]], se, is_protein_coding, CONFIG)
      }
    }, error = function(e) {
      list(
        group = gn,
        n_samples = length(groups_ids[[gn]]),
        n_genes = NA,
        n_outliers = 0,
        outlier_ids = character(0),
        outlier_flags = rep(FALSE, length(groups_ids[[gn]])),
        K = NA,
        elapsed_time = NA,
        sample_ids = groups_ids[[gn]],
        SD = rep(NA, length(groups_ids[[gn]])),
        OD = rep(NA, length(groups_ids[[gn]])),
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
          expr = process_group(gn, groups_ids[[gn]], se, is_protein_coding, CONFIG)
        )
      } else {
        results[[gn]] <- process_group(gn, groups_ids[[gn]], se, is_protein_coding, CONFIG)
      }
    }, error = function(e) {
      cat("  Error:", e$message, "\n")
      results[[gn]] <- list(
        group = gn,
        n_samples = length(groups_ids[[gn]]),
        n_genes = NA,
        n_outliers = 0,
        outlier_ids = character(0),
        outlier_flags = rep(FALSE, length(groups_ids[[gn]])),
        K = NA,
        elapsed_time = NA,
        sample_ids = groups_ids[[gn]],
        SD = rep(NA, length(groups_ids[[gn]])),
        OD = rep(NA, length(groups_ids[[gn]])),
        error = e$message
      )
    })
  }
}

overall_end <- Sys.time()
total_time <- difftime(overall_end, overall_start, units = "secs")

cat(sprintf("\nAll groups processed in %.1f seconds\n", as.numeric(total_time)))

# ============================================================================
# Create SD vs OD diagnostic plots
# ============================================================================

cat("\n--- Creating SD vs OD diagnostic plots ---\n")

create_sd_od_plot <- function(result) {
  if (is.null(result$SD) || is.null(result$OD) || 
      all(is.na(result$SD)) || all(is.na(result$OD))) {
    return(NULL)
  }
  
  plot_df <- data.frame(
    SD = result$SD,
    OD = result$OD,
    outlier = factor(result$outlier_flags, levels = c(FALSE, TRUE), 
                     labels = c("Normal", "Outlier")),
    sample_id = result$sample_ids
  )
  
  thr_SD <- result$diagnostics$thr_SD
  thr_OD <- result$diagnostics$thr_OD
  
  p <- ggplot(plot_df, aes(x = SD, y = OD, color = outlier)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("Normal" = "#4A90E2", "Outlier" = "#E94B3C"),
                       name = "Status") +
    geom_vline(xintercept = thr_SD, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = thr_OD, linetype = "dashed", color = "gray40") +
    labs(
      title = sprintf("%s (n=%d, K=%d)", result$group, result$n_samples, result$K),
      subtitle = sprintf("Genes: %d, Outliers: %d", result$n_genes, result$n_outliers),
      x = "Score Distance (SD)",
      y = "Orthogonal Distance (OD)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "bottom"
    )
  
  return(p)
}

# Create plots for all groups
all_plots <- lapply(results, create_sd_od_plot)
all_plots <- all_plots[!sapply(all_plots, is.null)]

# RET groups (4×2 layout)
ret_groups <- c("R0_tumor", "R0_normal", "R_Low_tumor", "R_Low_normal",
                "R_High_tumor", "R_High_normal", "R1_tumor", "R1_normal")
ret_plots <- all_plots[names(all_plots) %in% ret_groups]

if (length(ret_plots) > 0) {
  # Reorder plots
  ret_order <- intersect(ret_groups, names(ret_plots))
  ret_plots <- ret_plots[ret_order]
  
  cat("  Displaying RET groups (4×2)...\n")
  ret_combined <- gridExtra::grid.arrange(
    grobs = ret_plots,
    ncol = 2, nrow = 4,
    top = grid::textGrob("RET Groups: SD vs OD Diagnostic", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  print(ret_combined)
}

# BRAF groups (2×2 layout)
braf_groups <- c("B0_tumor", "B0_normal", "B1_tumor", "B1_normal")
braf_plots <- all_plots[names(all_plots) %in% braf_groups]

if (length(braf_plots) > 0) {
  # Reorder plots
  braf_order <- intersect(braf_groups, names(braf_plots))
  braf_plots <- braf_plots[braf_order]
  
  cat("  Displaying BRAF groups (2×2)...\n")
  braf_combined <- gridExtra::grid.arrange(
    grobs = braf_plots,
    ncol = 2, nrow = 2,
    top = grid::textGrob("BRAF Groups: SD vs OD Diagnostic", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  print(braf_combined)
}

# Save plots to PDF
pdf_file <- paste0(paths$output, "05_sd_od_diagnostic.pdf")
pdf(pdf_file, width = 10, height = 14)

if (length(ret_plots) > 0) {
  gridExtra::grid.arrange(
    grobs = ret_plots,
    ncol = 2, nrow = 4,
    top = grid::textGrob("RET Groups: SD vs OD Diagnostic", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
}

if (length(braf_plots) > 0) {
  gridExtra::grid.arrange(
    grobs = braf_plots,
    ncol = 2, nrow = 2,
    top = grid::textGrob("BRAF Groups: SD vs OD Diagnostic", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
}

dev.off()
cat("  Saved: 05_sd_od_diagnostic.pdf\n")

# ============================================================================
# Create outlier mapping by sample_id
# ============================================================================

cat("\n--- Creating outlier mapping ---\n")

# Collect all outlier sample IDs
outlier_samples <- unlist(lapply(results, function(r) r$outlier_ids))
outlier_samples <- unique(outlier_samples)

cat("Total unique outlier samples:", length(outlier_samples), "\n")

# Report by group
cat("\nOutliers by group:\n")
for (g in names(results)) {
  r <- results[[g]]
  if (!is.null(r$error)) {
    cat(sprintf("  %s: ERROR - %s\n", g, r$error))
  } else {
    cat(sprintf("  %s: %d outliers (n=%d, K=%d, genes=%d)\n", 
                g, r$n_outliers, r$n_samples, r$K, r$n_genes))
    if (r$n_outliers > 0) {
      cat(sprintf("    Samples: %s\n", paste(r$outlier_ids, collapse = ", ")))
    }
  }
}

# ============================================================================
# Update case_master with outlier flags
# ============================================================================

cat("\n--- Updating case_master with outlier flags ---\n")

# Initialize new columns with NA (for unpaired cases)
thyr_case_master_stage1_filtered$has_outlier_tumor <- NA_integer_
thyr_case_master_stage1_filtered$has_outlier_normal <- NA_integer_

# For paired cases, initialize with 0
paired_idx <- which(thyr_case_master_stage1_filtered$has_pair)
thyr_case_master_stage1_filtered$has_outlier_tumor[paired_idx] <- 0L
thyr_case_master_stage1_filtered$has_outlier_normal[paired_idx] <- 0L

# Groups that were QC'd
qc_target_groups <- c("R0", "R_Low", "R_High", "R1", "B0", "B1")

# Update flags based on outlier detection
for (i in seq_len(nrow(thyr_case_master_stage1_filtered))) {
  case_id <- thyr_case_master_stage1_filtered$case_submitter_id[i]
  group <- thyr_case_master_stage1_filtered$group[i]
  
  # Skip if not in QC target groups or not paired
  if (!(group %in% qc_target_groups) || !thyr_case_master_stage1_filtered$has_pair[i]) {
    next
  }
  
  # Find tumor and normal samples for this case
  tumor_samples <- sample_metadata$sample_submitter_id[
    sample_metadata$case_submitter_id == case_id &
      sample_metadata$sample_type == "Primary Tumor" &
      grepl("_merged", sample_metadata$sample_submitter_id)
  ]
  
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
    n_paired = sum(has_pair),
    n_unpaired = sum(!has_pair),
    n_tumor_outliers = sum(has_outlier_tumor == 1, na.rm = TRUE),
    n_normal_outliers = sum(has_outlier_normal == 1, na.rm = TRUE),
    n_both_outliers = sum(has_outlier_tumor == 1 & has_outlier_normal == 1, na.rm = TRUE),
    n_clean = sum(has_outlier_tumor == 0 & has_outlier_normal == 0, na.rm = TRUE),
    .groups = "drop"
  )

print(outlier_summary)

# Overall statistics (QC target groups only)
qc_target_cases <- thyr_case_master_stage1_filtered %>%
  filter(group %in% qc_target_groups & has_pair)
qc_target_count <- nrow(qc_target_cases)
clean_qc_target <- sum(qc_target_cases$has_outlier_tumor == 0 & 
                         qc_target_cases$has_outlier_normal == 0, na.rm = TRUE)

non_qc_cases <- sum(!(thyr_case_master_stage1_filtered$group %in% qc_target_groups) |
                      !thyr_case_master_stage1_filtered$has_pair)

cat(sprintf("\nQC target cases (R0/R_Low/R_High/R1/B0/B1 paired): %d - Clean: %d (%.1f%%)\n",
            qc_target_count, clean_qc_target, clean_qc_target/qc_target_count*100))
cat(sprintf("Non-QC cases (B_Low/B_High or unpaired): %d\n", non_qc_cases))

# ============================================================================
# Processing summary
# ============================================================================

cat("\n=== Processing Summary ===\n")
cat(sprintf("%-12s %6s %6s %4s %8s %6s %6s %8s\n", 
            "Group", "N_Samp", "Genes", "K", "Outliers", "SD_Out", "OD_Out", "Time(s)"))
cat(paste(rep("-", 70), collapse = ""), "\n")

for (g in names(results)) {
  r <- results[[g]]
  if (!is.null(r$error)) {
    cat(sprintf("%-12s %6d %6s %4s %8s %6s %6s %8s ERROR\n",
                g, r$n_samples, "NA", "NA", "NA", "NA", "NA", "NA"))
  } else {
    cat(sprintf("%-12s %6d %6d %4d %8d %6d %6d %8.1f\n",
                g, r$n_samples, r$n_genes, r$K, r$n_outliers,
                r$diagnostics$hit_SD, r$diagnostics$hit_OD,
                r$elapsed_time))
  }
}

cat(paste(rep("-", 70), collapse = ""), "\n")

# ============================================================================
# Save results
# ============================================================================

cat("\n--- Saving results ---\n")

# Save the updated case_master
output_file <- paste0(paths$processed, "thyr_case_master_stage1_filtered.rds")
saveRDS(thyr_case_master_stage1_filtered, output_file)
cat("  Updated case_master saved:", basename(output_file), "\n")
cat(sprintf("    Contains: %d cases (R0/R_Low/R_High/R1 + B0/B_Low/B_High/B1)\n", nrow(thyr_case_master_stage1_filtered)))
cat("    Columns: has_outlier_tumor, has_outlier_normal (0/1/NA)\n")
cat("    Note: NA = unpaired, PCA QC not performed\n")

# Save processing details for diagnostic purposes
processing_info <- list(
  date = Sys.Date(),
  version = "v7.12",
  config = CONFIG,
  results = results,
  outlier_samples = outlier_samples,
  processing_time = total_time,
  summary = outlier_summary,
  gene_filter = "protein_coding + per-group filterByExpr",
  normalization = "MUREN"
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

cat("\n=== PCA Outlier Detection Complete (v7.12) ===\n")
cat("Configuration:\n")
cat("  Gene filter: Protein coding + per-group filterByExpr\n")
cat("  Normalization: MUREN (LTS)\n")
cat("  Method: Adaptive thresholds (IQR/MAD)\n")
cat("  IQR multiplier:", CONFIG$IQR_MULTIPLIER, "\n")
cat("  MAD multiplier:", CONFIG$MAD_MULTIPLIER, "\n")
cat("  Max outliers:", sprintf("%.0f%%", CONFIG$MAX_OUTLIERS_PCT * 100), "\n")
cat("  Seed:", CONFIG$PA_SEED, "\n")
cat("\nProcessing:\n")
cat("  Groups analyzed:", length(groups_ids), "(paired samples only)\n")
cat("  Total time:", sprintf("%.1f seconds\n", as.numeric(total_time)))
cat("\nResults:\n")
cat("  Total outlier samples:", length(outlier_samples), "\n")
cat("  QC target cases:", qc_target_count, "- Clean:", clean_qc_target, 
    sprintf("(%.1f%%)\n", clean_qc_target/qc_target_count*100))
cat("  Non-QC cases:", non_qc_cases, "(B_Low/B_High or unpaired)\n")
cat("\nOutputs:\n")
cat("  Main: thyr_case_master_stage1_filtered.rds\n")
cat("  Info: stage1_pca_info.rds, stage1_outlier_summary.csv\n")
cat("  Plot: 05_sd_od_diagnostic.pdf\n")
cat("\nNext: Use thyr_case_master_stage1_filtered for downstream analysis\n")
cat("      For QC target: Filter with has_outlier_tumor == 0 & has_outlier_normal == 0\n")
cat("      For non-QC (B_Low/B_High, unpaired): has_outlier_tumor/normal is NA\n")

# Restore SE to original variable name
thyr_se_strand2_nonzero <- se

# Clean up
rm(list = setdiff(ls(), c("paths", "thyr_case_master", "thyr_case_master_stage1_filtered", "thyr_se_strand2_nonzero")))
gc()