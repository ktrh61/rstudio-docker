# 08_pca_visualization.R - PCA Visualization of Normalized Data
# Purpose: Visualize normalized expression data using PCA with POC coloring
#          Quality checkpoint before DEG analysis
# Method: Standard PCA (prcomp) on DEGES-normalized logCPM
# Input: thyr_normalized_counts.rds (from 07), thyr_case_master_stage2_filtered.rds (from 06)
# Output: PCA plots for R0/R1 and B0/B1 comparisons
# Version: v7.3
# Date: 2025-12-08

source("analysis_v7/setup.R")

cat("\n=== PCA Visualization (v7.3) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Purpose: Quality checkpoint - visualize normalized data before DEG analysis\n")
cat("Method: PCA with POC coloring\n")

# Load packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})

# ============================================================================
# Configuration
# ============================================================================

CONFIG <- list(
  PRIOR_COUNT = 0.5,   # For logCPM transformation
  N_TOP_VAR = 2000     # Top variable genes for PCA (NULL = use all)
)

cat("\nConfiguration:\n")
cat("  Prior count:", CONFIG$PRIOR_COUNT, "\n")
cat("  Top variable genes:", ifelse(is.null(CONFIG$N_TOP_VAR), "All", CONFIG$N_TOP_VAR), "\n")

# ============================================================================
# Load data
# ============================================================================

cat("\n--- Loading data ---\n")

# Define comparisons to load
comparisons <- c("R0_vs_R1_tumor", "R0_vs_R1_normal", 
                 "B0_vs_B1_tumor", "B0_vs_B1_normal")

# Load normalized CPM files from 07
thyr_normalized_counts <- list()
missing_files <- c()

for (comp in comparisons) {
  cpm_path <- paste0(paths$processed, "analysis_cpm_", comp, ".rds")
  if (file.exists(cpm_path)) {
    thyr_normalized_counts[[comp]] <- readRDS(cpm_path)
    cat(sprintf("  Loaded: analysis_cpm_%s.rds\n", comp))
  } else {
    missing_files <- c(missing_files, comp)
  }
}

if (length(missing_files) == length(comparisons)) {
  stop("No normalized CPM files found. Please run 07_deges_normalization.R first.")
}

if (length(missing_files) > 0) {
  cat("  Warning: Missing files for:", paste(missing_files, collapse = ", "), "\n")
}

cat("Normalized counts loaded:", length(thyr_normalized_counts), "comparisons\n")

# Load case master from 06
case_master_path <- paste0(paths$processed, "thyr_case_master_stage2_filtered.rds")
if (!file.exists(case_master_path)) {
  stop("thyr_case_master_stage2_filtered.rds not found. Please run 06_purity_analysis.R first.")
}
thyr_case_master <- readRDS(case_master_path)
cat("Case master loaded:", nrow(thyr_case_master), "cases\n")

# Load SE for sample metadata
if (exists("thyr_se_strand2_nonzero")) {
  se <- thyr_se_strand2_nonzero
  cat("Using existing thyr_se_strand2_nonzero\n")
} else {
  se_path <- paste0(paths$processed, "thyr_se_strand2_nonzero.rds")
  se <- readRDS(se_path)
  cat("SE loaded from file\n")
}
sample_metadata <- as.data.frame(colData(se))

# ============================================================================
# Helper functions
# ============================================================================

# Get POC value for a sample
get_poc_value <- function(sample_id, sample_meta, case_master) {
  # sample_id is sample_submitter_id format (e.g., SC299536_merged)
  # Need to find corresponding case_submitter_id
  case_id <- sample_meta$case_submitter_id[sample_meta$sample_submitter_id == sample_id]
  if (length(case_id) == 0) return(NA)
  
  # POC column is uppercase in case_master
  poc <- case_master$POC[case_master$case_submitter_id == case_id[1]]
  if (length(poc) == 0) return(NA)
  return(poc)
}

# Run PCA and create plot
run_pca_visualization <- function(comparison_name, norm_counts, sample_meta, case_master, config) {
  
  cat(sprintf("\n--- Processing: %s ---\n", comparison_name))
  
  # Get normalized CPM matrix
  if (is.null(norm_counts[[comparison_name]])) {
    cat("  Skipping: No data available\n")
    return(NULL)
  }
  
  cpm_matrix <- norm_counts[[comparison_name]]
  cat(sprintf("  Samples: %d, Genes: %d\n", ncol(cpm_matrix), nrow(cpm_matrix)))
  
  # Convert to logCPM for PCA
  logcpm <- log2(cpm_matrix + config$PRIOR_COUNT)
  
  # Select top variable genes if specified
  if (!is.null(config$N_TOP_VAR) && config$N_TOP_VAR < nrow(logcpm)) {
    gene_vars <- apply(logcpm, 1, var)
    top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:config$N_TOP_VAR]
    logcpm_subset <- logcpm[top_genes, ]
    cat(sprintf("  Using top %d variable genes\n", config$N_TOP_VAR))
  } else {
    logcpm_subset <- logcpm
  }
  
  # Run PCA
  pca_result <- prcomp(t(logcpm_subset), center = TRUE, scale. = FALSE)
  
  # Calculate variance explained
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  cat(sprintf("  PC1: %.1f%%, PC2: %.1f%%\n", var_explained[1], var_explained[2]))
  
  # Extract group and tissue info from comparison name
  # Format: "R0_vs_R1_tumor" or "B0_vs_B1_normal"
  parts <- strsplit(comparison_name, "_")[[1]]
  group1 <- parts[1]  # R0 or B0
  group2 <- parts[3]  # R1 or B1
  tissue <- parts[4]  # tumor or normal
  
  # Create data frame for plotting
  sample_ids <- colnames(cpm_matrix)
  
  # Debug: show first few sample IDs
  cat(sprintf("  Sample IDs (first 3): %s\n", paste(head(sample_ids, 3), collapse = ", ")))
  
  # Determine group for each sample
  sample_groups <- sapply(sample_ids, function(sid) {
    case_id <- sample_meta$case_submitter_id[sample_meta$sample_submitter_id == sid]
    if (length(case_id) == 0) return(NA)
    grp <- case_master$group[case_master$case_submitter_id == case_id[1]]
    if (length(grp) == 0) return(NA)
    return(grp)
  })
  
  # Get POC values
  poc_values <- sapply(sample_ids, function(sid) {
    get_poc_value(sid, sample_meta, case_master)
  })
  
  # Ensure POC is numeric
  poc_values <- as.numeric(poc_values)
  
  plot_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    sample_id = sample_ids,
    group = sample_groups,
    POC = poc_values,
    stringsAsFactors = FALSE
  )
  
  # Debug: report POC stats
  cat(sprintf("  POC range: %.1f - %.1f (NA: %d)\n", 
              min(poc_values, na.rm = TRUE), 
              max(poc_values, na.rm = TRUE),
              sum(is.na(poc_values))))
  
  # Create plot
  p <- ggplot(plot_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = POC, shape = group), size = 3.5, alpha = 0.8) +
    scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                          midpoint = 50, limits = c(0, 100),
                          name = "POC", na.value = "gray50") +
    scale_shape_manual(values = c(16, 17, 15, 18),
                       name = "Group") +
    labs(
      title = sprintf("%s: %s vs %s", tools::toTitleCase(tissue), group1, group2),
      subtitle = sprintf("PC1: %.1f%%, PC2: %.1f%% | n=%d | %d genes",
                         var_explained[1], var_explained[2], 
                         ncol(cpm_matrix), nrow(logcpm_subset)),
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2])
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "right"
    )
  
  return(list(
    plot = p,
    pca = pca_result,
    var_explained = var_explained,
    plot_df = plot_df
  ))
}

# ============================================================================
# Generate PCA visualizations
# ============================================================================

cat("\n--- Generating PCA visualizations ---\n")

# Run PCA for each comparison (using comparisons defined during data loading)
pca_results <- list()
all_plots <- list()

for (comp in names(thyr_normalized_counts)) {
  result <- run_pca_visualization(comp, thyr_normalized_counts, 
                                  sample_metadata, thyr_case_master, CONFIG)
  if (!is.null(result)) {
    pca_results[[comp]] <- result
    all_plots[[comp]] <- result$plot
  }
}

# ============================================================================
# Display plots
# ============================================================================

cat("\n--- Displaying plots ---\n")

# RET comparisons (2×2)
ret_comps <- c("R0_vs_R1_tumor", "R0_vs_R1_normal")
ret_plots <- all_plots[names(all_plots) %in% ret_comps]

if (length(ret_plots) > 0) {
  ret_order <- intersect(ret_comps, names(ret_plots))
  ret_plots <- ret_plots[ret_order]
  
  cat("  Displaying RET comparisons (1×2)...\n")
  ret_combined <- gridExtra::grid.arrange(
    grobs = ret_plots,
    ncol = 2, nrow = 1,
    top = grid::textGrob("RET Fusion: R0 vs R1 Comparison", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  print(ret_combined)
}

# BRAF comparisons (2×2)
braf_comps <- c("B0_vs_B1_tumor", "B0_vs_B1_normal")
braf_plots <- all_plots[names(all_plots) %in% braf_comps]

if (length(braf_plots) > 0) {
  braf_order <- intersect(braf_comps, names(braf_plots))
  braf_plots <- braf_plots[braf_order]
  
  cat("  Displaying BRAF comparisons (1×2)...\n")
  braf_combined <- gridExtra::grid.arrange(
    grobs = braf_plots,
    ncol = 2, nrow = 1,
    top = grid::textGrob("BRAF V600E: B0 vs B1 Comparison", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  print(braf_combined)
}

# ============================================================================
# Save plots to PDF
# ============================================================================

cat("\n--- Saving plots ---\n")

pdf_file <- paste0(paths$output, "08_pca_visualization.pdf")
pdf(pdf_file, width = 12, height = 6)

if (length(ret_plots) > 0) {
  gridExtra::grid.arrange(
    grobs = ret_plots,
    ncol = 2, nrow = 1,
    top = grid::textGrob("RET Fusion: R0 vs R1 Comparison", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
}

if (length(braf_plots) > 0) {
  gridExtra::grid.arrange(
    grobs = braf_plots,
    ncol = 2, nrow = 1,
    top = grid::textGrob("BRAF V600E: B0 vs B1 Comparison", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
}

dev.off()
cat("  Saved: 08_pca_visualization.pdf (2 pages)\n")

# ============================================================================
# Save PCA results
# ============================================================================

pca_output <- list(
  date = Sys.Date(),
  version = "v7.3",
  config = CONFIG,
  results = pca_results
)

saveRDS(pca_output, paste0(paths$output, "08_pca_results.rds"))
cat("  Saved: 08_pca_results.rds\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n=== PCA Visualization Complete (v7.3) ===\n")
cat("\nSummary:\n")

for (comp in names(pca_results)) {
  res <- pca_results[[comp]]
  cat(sprintf("  %s: PC1=%.1f%%, PC2=%.1f%%, n=%d\n",
              comp, res$var_explained[1], res$var_explained[2],
              nrow(res$plot_df)))
}

cat("\nOutputs:\n")
cat("  Plot: 08_pca_visualization.pdf\n")
cat("  Data: 08_pca_results.rds\n")

cat("\nInterpretation guide:\n")
cat("  - POC coloring: Blue (0) → White (50) → Red (100)\n
")
cat("  - Shape indicates group (R0/R1 or B0/B1)\n")
cat("  - If groups cluster separately, there may be batch effects\n")
cat("  - If POC gradient visible, radiation signature may be detectable\n")

cat("\nNext: Proceed to 09_deg_analysis.R if visualization looks acceptable\n")

# ============================================================================
# Cleanup
# ============================================================================

rm(list = setdiff(ls(), c("paths", "thyr_case_master", "thyr_case_master_stage2_filtered", 
                          "thyr_se_strand2_nonzero", "thyr_normalized_counts")))
gc()

cat("\n=== Done ===\n")