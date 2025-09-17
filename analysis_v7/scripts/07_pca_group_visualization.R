# 07_pca_group_visualization.R - PCA Group Visualization with CDM
# Purpose: Visualize high-purity samples using CDM-based PCA for 4 groups
# Groups: RET_tumor, RET_normal, BRAF_tumor, BRAF_normal
# Input: thyr_case_master_stage2_filtered, thyr_se_strand2_nonzero
# Output: PCA plots colored by POC status (0/1)
# Version: v7.1 - Serial outer loop with parallel CDM (BLAS 3 threads)
# Date: 2025-09-16

source("analysis_v7/setup.R")

cat("\n=== PCA Group Visualization with CDM (v7.1) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Analysis: 4 separate PCAs for RET/BRAF × tumor/normal\n")
cat("Method: CDM with internal BLAS parallelization\n")
cat("Processing: Serial groups, parallel CDM (3 threads)\n")

# Load packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(ggplot2)
  library(dplyr)
  library(Rcpp)
  library(RcppArmadillo)
  library(gridExtra)
})

# ============================================================================
# Configuration
# ============================================================================

CONFIG <- list(
  PRIOR_COUNT = 0.5,           # Same as 05/06
  USE_HIGH_PURITY_ONLY = TRUE, # Use only low_purity == 0
  CDM_THREADS = 3,             # BLAS threads for CDM (same as 06)
  VERBOSE_CDM = FALSE,         # CDM verbosity
  PC_TO_PLOT = c(1, 2),       # Plot PC1 vs PC2
  POINT_SIZE = 3,              # Point size for plots
  ALPHA = 0.8                  # Point transparency
)

cat("\nConfiguration:\n")
cat("  Prior count:", CONFIG$PRIOR_COUNT, "\n")
cat("  High purity only:", CONFIG$USE_HIGH_PURITY_ONLY, "\n")
cat("  CDM BLAS threads:", CONFIG$CDM_THREADS, "\n")
cat("  PCs to plot: PC", CONFIG$PC_TO_PLOT[1], " vs PC", CONFIG$PC_TO_PLOT[2], "\n")

# Source BLAS thread control utility
if (file.exists("utils/with_openblas_threads.R")) {
  source("utils/with_openblas_threads.R")
  cat("  BLAS control: with_openblas_threads loaded\n")
  has_blas_control <- TRUE
} else {
  cat("  BLAS control: Manual fallback\n")
  has_blas_control <- FALSE
  # Set threads to 1 by default (will be overridden during CDM)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
}

# ============================================================================
# Compile CDM
# ============================================================================

cat("\nCompiling CDM implementation...\n")
sourceCpp("utils/CDM_fast3_arma_enhanced.cpp")
cat("  CDM compiled successfully\n")

# ============================================================================
# Load data
# ============================================================================

cat("\n--- Loading data ---\n")

# Load SummarizedExperiment
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
cat("SE dimensions:", format(nrow(se), big.mark=","), "genes × ", 
    format(ncol(se), big.mark=","), "samples\n")

# Load case_master_stage2_filtered
if (exists("thyr_case_master_stage2_filtered")) {
  cat("Using existing thyr_case_master_stage2_filtered\n")
  case_master <- thyr_case_master_stage2_filtered
} else {
  case_master_path <- paste0(paths$processed, "thyr_case_master_stage2_filtered.rds")
  if (!file.exists(case_master_path)) {
    stop("thyr_case_master_stage2_filtered.rds not found. Please run 06_purity_analysis.R first.")
  }
  cat("Loading case_master from file...\n")
  case_master <- readRDS(case_master_path)
}
cat("Case master loaded:", nrow(case_master), "clean cases\n")

# Get sample metadata
sample_metadata <- as.data.frame(colData(se))

# ============================================================================
# Filter to high purity cases
# ============================================================================

if (CONFIG$USE_HIGH_PURITY_ONLY) {
  cat("\n--- Filtering to high purity cases ---\n")
  high_purity_cases <- case_master %>%
    filter(low_purity == 0)
  cat("High purity cases:", nrow(high_purity_cases), "/", nrow(case_master), 
      sprintf("(%.1f%%)\n", nrow(high_purity_cases)/nrow(case_master)*100))
} else {
  high_purity_cases <- case_master
  cat("\n--- Using all clean cases ---\n")
}

# Group distribution
group_summary <- high_purity_cases %>%
  group_by(group) %>%
  summarise(n = n(), .groups = "drop")
cat("\nHigh purity cases by group:\n")
print(group_summary)

# ============================================================================
# Helper function: Extract samples for a driver/tissue combination (using lapply)
# ============================================================================

extract_samples_for_analysis <- function(driver_name, tissue_type, case_df, sample_meta) {
  # Get cases for this driver (R0+R1 or B0+B1)
  if (driver_name == "RET") {
    driver_cases <- case_df %>%
      filter(group %in% c("R0", "R1"))
  } else if (driver_name == "BRAF") {
    driver_cases <- case_df %>%
      filter(group %in% c("B0", "B1"))
  } else {
    stop("Unknown driver: ", driver_name)
  }
  
  if (nrow(driver_cases) == 0) {
    return(list(sample_ids = character(0), poc_status = numeric(0), cases = character(0)))
  }
  
  # Determine sample type filter
  if (tissue_type == "tumor") {
    sample_type_filter <- "Primary Tumor"
  } else if (tissue_type == "normal") {
    sample_type_filter <- "Solid Tissue Normal"
  } else {
    stop("Unknown tissue type: ", tissue_type)
  }
  
  # Use lapply instead of for loop for sample extraction
  result_list <- lapply(seq_len(nrow(driver_cases)), function(i) {
    case_id <- driver_cases$case_submitter_id[i]
    case_group <- driver_cases$group[i]
    
    # Find sample for this case
    sample_id <- sample_meta$sample_submitter_id[
      sample_meta$case_submitter_id == case_id &
        sample_meta$sample_type == sample_type_filter &
        grepl("_merged", sample_meta$sample_submitter_id)
    ]
    
    if (length(sample_id) > 0) {
      list(
        sample_id = sample_id[1],  # Take first if multiple
        poc_status = ifelse(case_group %in% c("R0", "B0"), 0, 1),
        case_id = case_id
      )
    } else {
      NULL
    }
  })
  
  # Remove NULL entries and combine results
  result_list <- result_list[!sapply(result_list, is.null)]
  
  if (length(result_list) == 0) {
    return(list(sample_ids = character(0), poc_status = numeric(0), cases = character(0)))
  }
  
  # Extract vectors from list
  sample_ids <- sapply(result_list, function(x) x$sample_id)
  poc_status <- sapply(result_list, function(x) x$poc_status)
  case_ids <- sapply(result_list, function(x) x$case_id)
  
  return(list(
    sample_ids = sample_ids,
    poc_status = poc_status,
    cases = case_ids,
    driver = driver_name,
    tissue = tissue_type
  ))
}

# ============================================================================
# Define analysis groups
# ============================================================================

cat("\n--- Defining analysis groups ---\n")

# Create 4 analysis groups using lapply
group_definitions <- list(
  RET_tumor = list(driver = "RET", tissue = "tumor"),
  RET_normal = list(driver = "RET", tissue = "normal"),
  BRAF_tumor = list(driver = "BRAF", tissue = "tumor"),
  BRAF_normal = list(driver = "BRAF", tissue = "normal")
)

analysis_groups <- lapply(names(group_definitions), function(g) {
  def <- group_definitions[[g]]
  extract_samples_for_analysis(def$driver, def$tissue, high_purity_cases, sample_metadata)
})
names(analysis_groups) <- names(group_definitions)

# Report group sizes
cat("\nAnalysis groups:\n")
for (g in names(analysis_groups)) {
  n_total <- length(analysis_groups[[g]]$sample_ids)
  n_poc0 <- sum(analysis_groups[[g]]$poc_status == 0)
  n_poc1 <- sum(analysis_groups[[g]]$poc_status == 1)
  cat(sprintf("  %s: %d samples (POC=0: %d, POC=1: %d)\n", 
              g, n_total, n_poc0, n_poc1))
}

# ============================================================================
# Unified gene filtering
# ============================================================================

cat("\n--- Applying unified gene filtering ---\n")

# Collect all sample IDs from 4 groups
all_analysis_sample_ids <- unlist(lapply(analysis_groups, function(x) x$sample_ids))
cat("Total samples for filtering:", length(all_analysis_sample_ids), "\n")

# Extract counts for these samples
all_analysis_counts <- assay(se)[, all_analysis_sample_ids, drop = FALSE]

# Handle NA values
if (any(is.na(all_analysis_counts))) {
  cat("  Warning: Found", sum(is.na(all_analysis_counts)), "NA values. Replacing with 0.\n")
  all_analysis_counts[is.na(all_analysis_counts)] <- 0
}

# Create group labels for filterByExpr (4 groups)
group_labels <- factor(rep(names(analysis_groups), 
                           sapply(analysis_groups, function(x) length(x$sample_ids))))

# Apply filterByExpr
keep_genes <- filterByExpr(all_analysis_counts, group = group_labels)

cat("Gene filtering results:\n")
cat("  Starting genes:", nrow(se), "\n")
cat("  After filterByExpr:", sum(keep_genes), "\n")
cat("  Retention rate:", sprintf("%.1f%%\n", sum(keep_genes) / nrow(se) * 100))

# ============================================================================
# CDM PCA function with BLAS thread control
# ============================================================================

run_cdm_pca <- function(group_info, se, keep_genes, prior_count = 0.5, 
                        cdm_threads = 3, verbose = FALSE) {
  # Extract counts
  counts <- assay(se)[keep_genes, group_info$sample_ids, drop = FALSE]
  
  if (any(is.na(counts))) {
    counts[is.na(counts)] <- 0
  }
  
  # logCPM transformation
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge)
  logcpm <- cpm(dge, log = TRUE, prior.count = prior_count)
  
  # Run CDM with appropriate BLAS threads
  if (exists("has_blas_control") && has_blas_control) {
    # Use with_openblas_threads for controlled parallelization
    cdm_result <- with_openblas_threads(
      threads = cdm_threads,
      expr = CDM_fast3_arma(logcpm, verbose = verbose)
    )
  } else {
    # Manual BLAS control
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      old_threads <- RhpcBLASctl::blas_get_num_procs()
      RhpcBLASctl::blas_set_num_threads(cdm_threads)
      cdm_result <- CDM_fast3_arma(logcpm, verbose = verbose)
      RhpcBLASctl::blas_set_num_threads(old_threads)
    } else {
      # Fallback: run without thread control
      cdm_result <- CDM_fast3_arma(logcpm, verbose = verbose)
    }
  }
  
  # Create result object
  result <- list(
    cdm = cdm_result,
    group_info = group_info,
    n_samples = length(group_info$sample_ids),
    n_genes = sum(keep_genes)
  )
  
  return(result)
}

# ============================================================================
# Run CDM for each group (serial processing, parallel CDM)
# ============================================================================

cat("\n--- Running CDM for each group (serial, with parallel BLAS) ---\n")

cdm_results <- list()

for (group_name in names(analysis_groups)) {
  cat(sprintf("[%s] Processing %s...\n", Sys.time(), group_name))
  
  if (length(analysis_groups[[group_name]]$sample_ids) < 3) {
    cat("  Skipping: insufficient samples\n")
    next
  }
  
  start_time <- Sys.time()
  
  cdm_results[[group_name]] <- run_cdm_pca(
    group_info = analysis_groups[[group_name]],
    se = se,
    keep_genes = keep_genes,
    prior_count = CONFIG$PRIOR_COUNT,
    cdm_threads = CONFIG$CDM_THREADS,
    verbose = CONFIG$VERBOSE_CDM
  )
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat(sprintf("  Completed: %d samples, %d genes (%.1f seconds)\n", 
              cdm_results[[group_name]]$n_samples,
              cdm_results[[group_name]]$n_genes,
              elapsed))
}

# ============================================================================
# Create PCA plots using lapply
# ============================================================================

cat("\n--- Creating PCA plots ---\n")

# Function to create a single PCA plot
create_pca_plot <- function(cdm_result, title, pc_x = 1, pc_y = 2) {
  # Extract scores
  scores <- cdm_result$cdm$scores
  
  # Check if we have enough PCs
  if (ncol(scores) < max(pc_x, pc_y)) {
    warning(sprintf("Not enough PCs for %s. Available: %d", title, ncol(scores)))
    return(NULL)
  }
  
  # Create data frame for plotting
  plot_df <- data.frame(
    PC_x = scores[, pc_x],
    PC_y = scores[, pc_y],
    POC = factor(cdm_result$group_info$poc_status,
                 levels = c(0, 1),
                 labels = c("POC=NA", "POC≥66.6")),
    sample_id = cdm_result$group_info$sample_ids
  )
  
  # Create plot (without eigenvalue percentages)
  p <- ggplot(plot_df, aes(x = PC_x, y = PC_y, color = POC)) +
    geom_point(size = CONFIG$POINT_SIZE, alpha = CONFIG$ALPHA) +
    scale_color_manual(values = c("POC=NA" = "#4A90E2", "POC≥66.6" = "#E94B3C"),
                       name = "POC Status") +
    labs(
      title = title,
      x = sprintf("PC%d", pc_x),
      y = sprintf("PC%d", pc_y),
      subtitle = sprintf("n = %d samples", nrow(plot_df))
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    coord_fixed()
  
  return(p)
}

# Create plots for each group using lapply
plot_list <- lapply(names(cdm_results), function(group_name) {
  # Create title
  parts <- strsplit(group_name, "_")[[1]]
  title <- sprintf("%s %s", parts[1], parts[2])
  
  create_pca_plot(
    cdm_result = cdm_results[[group_name]],
    title = title,
    pc_x = CONFIG$PC_TO_PLOT[1],
    pc_y = CONFIG$PC_TO_PLOT[2]
  )
})
names(plot_list) <- names(cdm_results)

# Remove NULL entries
plots <- plot_list[!sapply(plot_list, is.null)]

# ============================================================================
# Combine and save plots (stage3 naming)
# ============================================================================

if (length(plots) > 0) {
  cat("\n--- Saving plots ---\n")
  
  # Combine plots in 2x2 grid
  combined_plot <- gridExtra::grid.arrange(
    plots[["RET_tumor"]], plots[["RET_normal"]],
    plots[["BRAF_tumor"]], plots[["BRAF_normal"]],
    ncol = 2, nrow = 2,
    top = "PCA Analysis - High Purity Samples (Stage 3)"
  )
  
  # Save as PDF with stage3 naming
  pdf_file <- paste0(paths$output, "stage3_pca_plots.pdf")
  pdf(pdf_file, width = 12, height = 10)
  gridExtra::grid.arrange(
    plots[["RET_tumor"]], plots[["RET_normal"]],
    plots[["BRAF_tumor"]], plots[["BRAF_normal"]],
    ncol = 2, nrow = 2,
    top = "PCA Analysis - High Purity Samples (Stage 3)"
  )
  dev.off()
  cat("  Saved combined plot:", basename(pdf_file), "\n")
  
  # Save individual plots with stage3 naming
  for (group_name in names(plots)) {
    png_file <- paste0(paths$output, sprintf("stage3_pca_%s.png", group_name))
    ggsave(png_file, plots[[group_name]], width = 8, height = 6, dpi = 300)
    cat("  Saved individual plot:", basename(png_file), "\n")
  }
}

# ============================================================================
# Save CDM results (stage3 naming)
# ============================================================================

cat("\n--- Saving CDM results ---\n")

cdm_output <- list(
  date = Sys.Date(),
  config = CONFIG,
  groups = analysis_groups,
  cdm_results = cdm_results,
  keep_genes = which(keep_genes),
  n_genes_used = sum(keep_genes)
)

saveRDS(cdm_output, paste0(paths$output, "stage3_pca_results.rds"))
cat("  CDM results saved: stage3_pca_results.rds\n")

# ============================================================================
# Final report
# ============================================================================

cat("\n=== Stage 3 Visualization Complete ===\n")
cat("Configuration:\n")
cat("  Method: CDM with internal BLAS parallelization\n")
cat("  Processing: Serial groups, ", CONFIG$CDM_THREADS, " BLAS threads per CDM\n")
cat("  Prior count:", CONFIG$PRIOR_COUNT, "\n")
cat("  Genes used:", sum(keep_genes), sprintf("(%.1f%%)\n", sum(keep_genes)/nrow(se)*100))
cat("\nGroups analyzed:\n")
for (g in names(analysis_groups)) {
  n_total <- length(analysis_groups[[g]]$sample_ids)
  n_poc0 <- sum(analysis_groups[[g]]$poc_status == 0)
  n_poc1 <- sum(analysis_groups[[g]]$poc_status == 1)
  cat(sprintf("  %s: %d samples (POC=0: %d, POC=1: %d)\n", 
              g, n_total, n_poc0, n_poc1))
}
cat("\nOutputs:\n")
cat("  Plots: stage3_pca_plots.pdf (combined)\n")
cat("  Individual plots: stage3_pca_*.png\n")
cat("  Data: stage3_pca_results.rds\n")

# Clean up (keep only essential objects)
rm(list = setdiff(ls(), c("paths", "thyr_case_master_stage2_filtered", "thyr_se_strand2_nonzero")))
gc()