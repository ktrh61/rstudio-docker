# ==============================================================================
# REBC-THYR DEG Overlap Analysis Script v6
# 10_deg_overlap_analysis_v6.R
# ==============================================================================

# Required libraries
library(UpSetR)
library(pheatmap)
library(dplyr)

cat("Starting DEG overlap analysis v6 with directional patterns...\n")

# ==============================================================================
# 1. Load DEG Analysis Results v6
# ==============================================================================

cat("Loading DEG analysis v6 results...\n")

# Load v6 DEG analysis results
load("./data/processed/final_deg_v6_results.rda")

cat("DEG analysis v6 results loaded\n")
cat("Available comparisons:", names(final_deg_v6_results$deg_analysis_results), "\n")

# Use v6 results
deg_analysis_results <- final_deg_v6_results$deg_analysis_results

# ==============================================================================
# 2. Prepare Directional DEG Lists
# ==============================================================================

# Initialize directional DEG lists
directional_deg_lists <- list()

# Extract up/down-regulated genes for each comparison
for (comp_name in names(deg_analysis_results)) {
  cat(sprintf("Processing %s for overlap analysis...\n", comp_name))
  
  results_df <- deg_analysis_results[[comp_name]]$deg_summary$results_df
  
  # Skip if no significant results
  if (nrow(results_df) == 0 || sum(results_df$significant, na.rm = TRUE) == 0) {
    cat(sprintf("No significant DEGs in %s, skipping...\n", comp_name))
    next
  }
  
  # Extract upregulated genes
  up_genes <- results_df$gene_id[results_df$significant & results_df$log2FC > 0 & !is.na(results_df$log2FC)]
  up_list_name <- paste0(comp_name, "_UP")
  directional_deg_lists[[up_list_name]] <- up_genes
  
  # Extract downregulated genes  
  down_genes <- results_df$gene_id[results_df$significant & results_df$log2FC < 0 & !is.na(results_df$log2FC)]
  down_list_name <- paste0(comp_name, "_DOWN")
  directional_deg_lists[[down_list_name]] <- down_genes
  
  cat(sprintf("  %s: %d upregulated, %d downregulated\n", 
              comp_name, length(up_genes), length(down_genes)))
}

# Display final lists summary
cat("\nDirectional DEG lists summary:\n")
for (list_name in names(directional_deg_lists)) {
  cat(sprintf("  %s: %d genes\n", list_name, length(directional_deg_lists[[list_name]])))
}

# ==============================================================================
# 3. Analyze Key Overlap Patterns
# ==============================================================================

cat("\n=== Analyzing key overlap patterns ===\n")

# Get lists with safe checking
get_deg_list <- function(name) {
  if (name %in% names(directional_deg_lists)) {
    return(directional_deg_lists[[name]])
  } else {
    return(character(0))
  }
}

# RET comparisons
ret_tumor_up <- get_deg_list("R0_vs_R1_tumor_UP")
ret_tumor_down <- get_deg_list("R0_vs_R1_tumor_DOWN")
ret_normal_up <- get_deg_list("R0_vs_R1_normal_UP")
ret_normal_down <- get_deg_list("R0_vs_R1_normal_DOWN")

# Calculate consistent patterns
ret_consistent_up <- character(0)
ret_consistent_down <- character(0)

if (length(ret_tumor_up) > 0 && length(ret_normal_up) > 0) {
  ret_consistent_up <- intersect(ret_tumor_up, ret_normal_up)
  cat(sprintf("RET consistent upregulation (tumor ∩ normal): %d genes\n", length(ret_consistent_up)))
}

if (length(ret_tumor_down) > 0 && length(ret_normal_down) > 0) {
  ret_consistent_down <- intersect(ret_tumor_down, ret_normal_down)
  cat(sprintf("RET consistent downregulation (tumor ∩ normal): %d genes\n", length(ret_consistent_down)))
}

# Total consistent genes
total_consistent <- length(ret_consistent_up) + length(ret_consistent_down)
cat(sprintf("\nTotal consistent DEGs (UP + DOWN): %d genes\n", total_consistent))

# Opposite directions
if (length(ret_tumor_up) > 0 && length(ret_normal_down) > 0) {
  ret_opposite_1 <- intersect(ret_tumor_up, ret_normal_down)
  cat(sprintf("RET opposite regulation (tumor↑, normal↓): %d genes\n", length(ret_opposite_1)))
}

if (length(ret_tumor_down) > 0 && length(ret_normal_up) > 0) {
  ret_opposite_2 <- intersect(ret_tumor_down, ret_normal_up)
  cat(sprintf("RET opposite regulation (tumor↓, normal↑): %d genes\n", length(ret_opposite_2)))
}

# ==============================================================================
# 4. Create Summary Table
# ==============================================================================

cat("\n=== Creating overlap summary ===\n")

# Summary of key patterns
overlap_summary <- data.frame(
  Pattern = c(
    "R0_vs_R1_tumor DEGs",
    "R0_vs_R1_normal DEGs", 
    "Consistent UP (tumor ∩ normal)",
    "Consistent DOWN (tumor ∩ normal)",
    "Total Consistent",
    "Tumor-specific UP",
    "Tumor-specific DOWN",
    "Normal-specific UP",
    "Normal-specific DOWN"
  ),
  Count = c(
    length(ret_tumor_up) + length(ret_tumor_down),
    length(ret_normal_up) + length(ret_normal_down),
    length(ret_consistent_up),
    length(ret_consistent_down),
    total_consistent,
    length(setdiff(ret_tumor_up, ret_normal_up)),
    length(setdiff(ret_tumor_down, ret_normal_down)),
    length(setdiff(ret_normal_up, ret_tumor_up)),
    length(setdiff(ret_normal_down, ret_tumor_down))
  ),
  stringsAsFactors = FALSE
)

print(overlap_summary)

# Save summary
write.csv(overlap_summary, "./output/reports/deg_overlap_summary_v6.csv", row.names = FALSE)

# ==============================================================================
# 5. Export Key Gene Lists
# ==============================================================================

cat("\n=== Exporting gene lists ===\n")

# Create directory
if (!dir.exists("./output/gene_lists")) {
  dir.create("./output/gene_lists", recursive = TRUE)
}

# Export consistent genes if they exist
if (length(ret_consistent_up) > 0) {
  write.table(ret_consistent_up, 
              file = "./output/gene_lists/RET_consistent_UP_v6.txt",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat(sprintf("Exported %d consistent UP genes\n", length(ret_consistent_up)))
}

if (length(ret_consistent_down) > 0) {
  write.table(ret_consistent_down,
              file = "./output/gene_lists/RET_consistent_DOWN_v6.txt", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat(sprintf("Exported %d consistent DOWN genes\n", length(ret_consistent_down)))
}

# ==============================================================================
# 6. Comparison with Previous Results
# ==============================================================================

cat("\n=== Comparison with previous analysis ===\n")

# Previous results (from documentation)
cat("Previous analysis (v3):\n")
cat("  Target: 321 consistent genes (247 UP + 74 DOWN)\n")
cat("  Status: Not achieved due to CDM modifications\n")

cat("\nCurrent analysis (v6):\n")
cat(sprintf("  R0_vs_R1_tumor: %d DEGs\n", length(ret_tumor_up) + length(ret_tumor_down)))
cat(sprintf("  R0_vs_R1_normal: %d DEGs\n", length(ret_normal_up) + length(ret_normal_down)))
cat(sprintf("  Consistent: %d genes (%d UP + %d DOWN)\n", 
            total_consistent, length(ret_consistent_up), length(ret_consistent_down)))

# Assessment
if (total_consistent > 100) {
  cat("\n✓ Sufficient consistent genes for feature selection\n")
} else if (total_consistent > 50) {
  cat("\n△ Moderate number of consistent genes\n")
} else {
  cat("\n✗ Limited consistent genes - consider tumor-only approach\n")
}

# ==============================================================================
# 7. Final Summary
# ==============================================================================

cat("\n==============================================\n")
cat("DEG Overlap Analysis v6 Summary\n")
cat("==============================================\n")

cat(sprintf("Total R0_vs_R1_tumor DEGs: %d\n", length(ret_tumor_up) + length(ret_tumor_down)))
cat(sprintf("Total R0_vs_R1_normal DEGs: %d\n", length(ret_normal_up) + length(ret_normal_down)))
cat(sprintf("Consistent directional genes: %d\n", total_consistent))

cat("\nRecommendations for feature selection:\n")
if (total_consistent >= 50) {
  cat("1. Use consistent genes for robust biomarkers\n")
  cat("2. Prioritize genes with largest fold changes\n")
} else {
  cat("1. Focus on R0_vs_R1_tumor DEGs (1,140 genes available)\n")
  cat("2. Use tumor-specific patterns for classification\n")
}

cat("\nFiles created:\n")
cat("  ./output/reports/deg_overlap_summary_v6.csv\n")
cat("  ./output/gene_lists/RET_consistent_*_v6.txt\n")

cat("\nOverlap analysis v6 completed!\n")