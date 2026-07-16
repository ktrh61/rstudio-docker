# ==============================================================================
# REBC-THYR Tumor Purity Analysis Script v6 - After PCA v6 (Simple logCPM)
# 06_purity_analysis_v6.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(edgeR)
library(dplyr)

# Source required normalization functions
source("./utils/utils_improved.R")
source("./utils/norm_improved.R")
source("./utils/contamde_purity_functions.R")

cat("Starting tumor purity analysis v6...\n")
cat("Using PCA v6 results (simple logCPM + CDM)\n")

# ==============================================================================
# 1. Load Data and PCA v6 Results
# ==============================================================================

cat("\nLoading PCA v6 results and data...\n")

# Load PCA v6 filtered sample lists
load("./data/processed/pca_filtered_sample_lists_v6.rda")

# Load SummarizedExperiment object for count data
if (!exists("se_thyr")) {
  load("./data/raw/thyr_data.rda")
  se_thyr <- data
  rm(data)
}

# Get stranded_second count data
count_data_full <- assay(se_thyr, "stranded_second")
gene_info <- rowData(se_thyr)

cat("Original count data dimensions:", dim(count_data_full), "\n")

# Display PCA v6 filtered sample counts
cat("\nPCA v6 filtered sample counts (input):\n")
for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(pca_filtered_sample_lists_v6)) {
    n_pairs <- length(pca_filtered_sample_lists_v6[[group]]$cases)
    cat(sprintf("  %s: %d pairs\n", group, n_pairs))
  }
}

# ==============================================================================
# 2. Filter to Protein-Coding Genes
# ==============================================================================

cat("\nFiltering to protein-coding genes...\n")

protein_coding_genes <- gene_info$gene_type == "protein_coding"
count_data_pc <- count_data_full[protein_coding_genes, ]
gene_info_pc <- gene_info[protein_coding_genes, ]

cat("After protein-coding filter:", dim(count_data_pc), "\n")

# ==============================================================================
# 3. Tumor Purity Estimation Function
# ==============================================================================

estimate_tumor_purity <- function(count_matrix, group_name, verbose = TRUE) {
  
  if (verbose) cat(sprintf("\nEstimating tumor purity for %s...\n", group_name))
  
  # Apply filterByExpr for low expression filtering
  dge <- DGEList(counts = count_matrix)
  n_pairs <- ncol(count_matrix) / 2
  group_factor <- factor(c(rep("normal", n_pairs), rep("tumor", n_pairs)))
  keep <- filterByExpr(dge, group = group_factor)
  dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
  
  if (verbose) {
    cat(sprintf("  After filterByExpr: %d genes (from %d)\n", 
                nrow(dge_filtered), nrow(count_matrix)))
  }
  
  # Run contamDE purity estimation with LTS options for stability
  tryCatch({
    contamde_result <- contamde_purity(
      counts = dge_filtered$counts,
      subtype = NULL,
      covariate = NULL,
      contaminated = TRUE,
      pairwise_method = "lts",  # Use LTS for robust estimation
      workers = "auto",
      voom_norm = "quantile",
      prior.count = 0,
      verbose = verbose
    )
    
    # Extract purity estimates
    purity_estimates <- contamde_result$proportion
    names(purity_estimates) <- paste0("Pair_", 1:length(purity_estimates))
    
    if (verbose) {
      cat(sprintf("  Purity range: [%.3f, %.3f], mean=%.3f\n",
                  min(purity_estimates), max(purity_estimates), mean(purity_estimates)))
    }
    
    return(list(
      purity_estimates = purity_estimates,
      contamde_result = contamde_result,
      filtered_genes = rownames(dge_filtered),
      n_pairs = n_pairs,
      success = TRUE,
      error_message = NULL
    ))
    
  }, error = function(e) {
    if (verbose) cat(sprintf("  Error: %s\n", e$message))
    return(list(
      purity_estimates = rep(NA, n_pairs),
      contamde_result = NULL,
      filtered_genes = character(0),
      n_pairs = n_pairs,
      success = FALSE,
      error_message = e$message
    ))
  })
}

# ==============================================================================
# 4. Process Each Group for Tumor Purity Estimation
# ==============================================================================

cat("\n=== Tumor Purity Estimation ===\n")

purity_results <- list()
purity_cutoff <- 0.6  # 60% purity threshold

for (group in c("R0", "R1", "B0", "B1")) {
  
  # Check if group has samples after PCA v6 filtering
  if (!group %in% names(pca_filtered_sample_lists_v6) || 
      length(pca_filtered_sample_lists_v6[[group]]$tumor) == 0) {
    cat(sprintf("%s: No samples available, skipping...\n", group))
    purity_results[[group]] <- NULL
    next
  }
  
  # Get PCA v6 filtered samples
  tumor_samples <- pca_filtered_sample_lists_v6[[group]]$tumor
  normal_samples <- pca_filtered_sample_lists_v6[[group]]$normal
  case_ids <- pca_filtered_sample_lists_v6[[group]]$cases
  
  # Verify sample availability
  tumor_available <- tumor_samples %in% colnames(count_data_pc)
  normal_available <- normal_samples %in% colnames(count_data_pc)
  
  if (!all(tumor_available) || !all(normal_available)) {
    missing_tumor <- sum(!tumor_available)
    missing_normal <- sum(!normal_available)
    cat(sprintf("%s: Warning - Missing samples (tumor: %d, normal: %d)\n", 
                group, missing_tumor, missing_normal))
  }
  
  # Use only available samples
  tumor_samples_final <- tumor_samples[tumor_available]
  normal_samples_final <- normal_samples[normal_available]
  min_pairs <- min(length(tumor_samples_final), length(normal_samples_final))
  
  if (min_pairs == 0) {
    cat(sprintf("%s: No valid pairs, skipping...\n", group))
    purity_results[[group]] <- NULL
    next
  }
  
  # Prepare count data: [normal1, normal2, ..., tumor1, tumor2, ...]
  group_samples <- c(normal_samples_final[1:min_pairs], tumor_samples_final[1:min_pairs])
  group_counts <- count_data_pc[, group_samples]
  
  cat(sprintf("%s: Processing %d pairs\n", group, min_pairs))
  
  # Run purity estimation
  purity_result <- estimate_tumor_purity(group_counts, group, verbose = TRUE)
  
  if (purity_result$success) {
    purity_results[[group]] <- list(
      group_name = group,
      purity_estimates = purity_result$purity_estimates,
      contamde_result = purity_result$contamde_result,
      tumor_samples = tumor_samples_final[1:min_pairs],
      normal_samples = normal_samples_final[1:min_pairs],
      case_ids = case_ids[1:min_pairs],
      filtered_genes = purity_result$filtered_genes,
      n_pairs = min_pairs,
      analysis_date = Sys.time(),
      success = TRUE
    )
  } else {
    cat(sprintf("%s: Purity estimation failed - %s\n", group, purity_result$error_message))
    purity_results[[group]] <- NULL
  }
}

# ==============================================================================
# 5. Purity Quality Assessment and High-Purity Sample Selection
# ==============================================================================

cat("\n=== Purity-Based Sample Selection ===\n")

# Initialize summary table
purity_summary <- data.frame(
  Group = c("R0", "R1", "B0", "B1"),
  Input_Pairs = integer(4),
  Mean_Purity = numeric(4),
  SD_Purity = numeric(4),
  Min_Purity = numeric(4),
  Max_Purity = numeric(4),
  High_Purity_Count = integer(4),
  Retention_Rate = numeric(4),
  stringsAsFactors = FALSE
)

# Create final high-purity filtered sample lists
final_high_purity_sample_lists_v6 <- list()

for (i in 1:4) {
  group <- purity_summary$Group[i]
  
  if (!is.null(purity_results[[group]])) {
    # Get purity estimates
    purity_est <- purity_results[[group]]$purity_estimates
    n_pairs <- purity_results[[group]]$n_pairs
    
    # Calculate summary statistics
    purity_summary$Input_Pairs[i] <- n_pairs
    purity_summary$Mean_Purity[i] <- round(mean(purity_est, na.rm = TRUE), 3)
    purity_summary$SD_Purity[i] <- round(sd(purity_est, na.rm = TRUE), 3)
    purity_summary$Min_Purity[i] <- round(min(purity_est, na.rm = TRUE), 3)
    purity_summary$Max_Purity[i] <- round(max(purity_est, na.rm = TRUE), 3)
    
    # Apply purity cutoff
    high_purity_indices <- which(purity_est >= purity_cutoff)
    n_high_purity <- length(high_purity_indices)
    
    purity_summary$High_Purity_Count[i] <- n_high_purity
    purity_summary$Retention_Rate[i] <- round(n_high_purity / n_pairs, 3)
    
    cat(sprintf("%s: %d/%d pairs ≥%.1f purity (%.1f%% retention)\n",
                group, n_high_purity, n_pairs, purity_cutoff, 
                purity_summary$Retention_Rate[i] * 100))
    
    # Create high-purity sample lists
    if (n_high_purity > 0) {
      final_high_purity_sample_lists_v6[[group]] <- list(
        tumor = purity_results[[group]]$tumor_samples[high_purity_indices],
        normal = purity_results[[group]]$normal_samples[high_purity_indices],
        cases = purity_results[[group]]$case_ids[high_purity_indices]
      )
    } else {
      final_high_purity_sample_lists_v6[[group]] <- list(
        tumor = character(0),
        normal = character(0), 
        cases = character(0)
      )
    }
    
  } else {
    # Group failed or had no data
    purity_summary[i, 2:8] <- c(0, NA, NA, NA, NA, 0, 0)
    
    final_high_purity_sample_lists_v6[[group]] <- list(
      tumor = character(0),
      normal = character(0),
      cases = character(0)
    )
  }
}

print(purity_summary)

# ==============================================================================
# 6. Summary and Quality Control
# ==============================================================================

cat("\n=== Complete Filtering Cascade Summary ===\n")

# Load original sample lists for comparison
load("./data/processed/sample_lists.rda")

cat("Filtering cascade (Original → PCA v6 → Purity-filtered):\n")
for (group in c("R0", "R1", "B0", "B1")) {
  orig_count <- length(sample_lists[[group]]$tumor)
  pca_count <- length(pca_filtered_sample_lists_v6[[group]]$tumor)
  purity_count <- length(final_high_purity_sample_lists_v6[[group]]$tumor)
  
  cat(sprintf("  %s: %d → %d → %d pairs\n", group, orig_count, pca_count, purity_count))
}

# Check minimum sample requirements
min_samples <- 5
viable_groups <- c()
cat("\nViability check (minimum 5 pairs per group):\n")
for (group in c("R0", "R1", "B0", "B1")) {
  final_count <- length(final_high_purity_sample_lists_v6[[group]]$tumor)
  if (final_count >= min_samples) {
    cat(sprintf("  %s: %d pairs ✓\n", group, final_count))
    viable_groups <- c(viable_groups, group)
  } else {
    cat(sprintf("  %s: %d pairs (insufficient)\n", group, final_count))
  }
}

total_final_pairs <- sum(sapply(final_high_purity_sample_lists_v6, function(x) length(x$tumor)))
cat(sprintf("\nTotal high-purity pairs: %d\n", total_final_pairs))
cat(sprintf("Viable groups for DEG analysis: %s\n", 
            ifelse(length(viable_groups) > 0, paste(viable_groups, collapse = ", "), "None")))

# ==============================================================================
# 7. Save Results
# ==============================================================================

cat("\nSaving purity analysis v6 results...\n")

# Create processed directory if needed
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save final high-purity sample lists (main output)
save(final_high_purity_sample_lists_v6, 
     file = "./data/processed/final_high_purity_sample_lists_v6.rda")

# Save comprehensive purity analysis results
purity_analysis_v6_results <- list(
  purity_results = purity_results,
  purity_summary = purity_summary,
  final_high_purity_sample_lists = final_high_purity_sample_lists_v6,
  purity_cutoff = purity_cutoff,
  viable_groups = viable_groups,
  protein_coding_genes = rownames(count_data_pc),
  analysis_date = Sys.time(),
  analysis_version = "v6_simple_logcpm_lts_muren"
)

save(purity_analysis_v6_results, 
     file = "./data/processed/purity_analysis_v6_results.rda")

# ==============================================================================
# 8. Define Viable Comparisons
# ==============================================================================

cat("\n=== Viable Comparisons for DEG Analysis ===\n")

possible_comparisons <- list()
if ("R0" %in% viable_groups && "R1" %in% viable_groups) {
  r0_pairs <- length(final_high_purity_sample_lists_v6[["R0"]]$tumor)
  r1_pairs <- length(final_high_purity_sample_lists_v6[["R1"]]$tumor)
  possible_comparisons$R0_vs_R1_tumor <- sprintf("RET Tumor: R0(%d) vs R1(%d)", r0_pairs, r1_pairs)
  possible_comparisons$R0_vs_R1_normal <- sprintf("RET Normal: R0(%d) vs R1(%d)", r0_pairs, r1_pairs)
}

if ("B0" %in% viable_groups && "B1" %in% viable_groups) {
  b0_pairs <- length(final_high_purity_sample_lists_v6[["B0"]]$tumor)
  b1_pairs <- length(final_high_purity_sample_lists_v6[["B1"]]$tumor)
  possible_comparisons$B0_vs_B1_tumor <- sprintf("BRAF Tumor: B0(%d) vs B1(%d)", b0_pairs, b1_pairs)
  possible_comparisons$B0_vs_B1_normal <- sprintf("BRAF Normal: B0(%d) vs B1(%d)", b0_pairs, b1_pairs)
}

if (length(possible_comparisons) > 0) {
  cat("Viable comparisons:\n")
  for (comp_name in names(possible_comparisons)) {
    cat(sprintf("  - %s: %s\n", comp_name, possible_comparisons[[comp_name]]))
  }
  
  # Save comparison info for next analysis
  comparison_info_v6 <- list(
    viable_groups = viable_groups,
    possible_comparisons = possible_comparisons,
    primary_comparison = names(possible_comparisons)[1],
    analysis_version = "v6"
  )
  save(comparison_info_v6, file = "./data/processed/comparison_info_v6.rda")
} else {
  cat("No viable comparisons available.\n")
  cat("Consider adjusting purity threshold or reviewing PCA filtering.\n")
}

# ==============================================================================
# 9. Final Summary
# ==============================================================================

cat("\n=== Purity Analysis v6 Complete ===\n")

if (total_final_pairs >= 10) {
  cat(sprintf("SUCCESS: %d high-purity pairs ready for downstream analysis\n", total_final_pairs))
  cat("Next steps:\n")
  cat("  1. Phase 2 PCA for group comparisons\n")
  cat("  2. DEG analysis\n")
} else if (total_final_pairs > 0) {
  cat(sprintf("LIMITED: Only %d high-purity pairs available\n", total_final_pairs))
  cat("Proceed with caution or adjust parameters.\n")
} else {
  cat("FAILED: No high-purity pairs available\n")
  cat("Review filtering parameters and data quality.\n")
}

cat("\nPurity analysis v6 completed!\n")