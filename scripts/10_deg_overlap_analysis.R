# ==============================================================================
# REBC-THYR DEG Overlap Analysis Script
# 10_deg_overlap_analysis.R
# ==============================================================================

# Required libraries
library(UpSetR)
library(pheatmap)
library(dplyr)

cat("Starting DEG overlap analysis with directional patterns...\n")

# ==============================================================================
# 1. Load DEG Analysis Results
# ==============================================================================

cat("Loading DEG analysis results from 09...\n")

# Load DEG analysis results
load("./data/processed/final_deg_results.rda")

cat("DEG analysis results loaded\n")
cat("Available comparisons:", names(final_deg_results$deg_analysis_results), "\n")

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
# 2. Create UpSet Plot
# ==============================================================================

if (length(directional_deg_lists) > 0) {
  cat("Creating UpSet plot for directional DEG overlap...\n")
  
  # Remove empty lists
  directional_deg_lists <- directional_deg_lists[sapply(directional_deg_lists, length) > 0]
  
  if (length(directional_deg_lists) >= 2) {
    # Create UpSet plot
    upset_plot <- upset(
      fromList(directional_deg_lists),
      order.by = "freq",
      decreasing = TRUE,
      nsets = length(directional_deg_lists),
      sets.bar.color = rep(c("red", "blue"), length.out = length(directional_deg_lists)),
      mainbar.y.label = "DEG Intersection Size",
      sets.x.label = "DEG Count per Comparison",
      text.scale = c(1.3, 1.3, 1, 1, 2, 0.75),
      point.size = 3.5,
      line.size = 2
    )
    
    # Save UpSet plot
    png("./output/plots/upset_directional_degs.png", width = 12, height = 8, units = "in", res = 300)
    print(upset_plot)
    dev.off()
    
    cat("UpSet plot saved to ./output/plots/upset_directional_degs.png\n")
    
  } else {
    cat("Insufficient comparisons for UpSet plot (need at least 2)\n")
  }
} else {
  cat("No directional DEG lists available for UpSet plot\n")
}

# ==============================================================================
# 3. Create Overlap Summary Table
# ==============================================================================

cat("Creating overlap summary table...\n")

# Calculate pairwise overlaps
overlap_summary <- data.frame(
  List1 = character(0),
  List2 = character(0),
  List1_Size = integer(0),
  List2_Size = integer(0),
  Overlap_Size = integer(0),
  Jaccard_Index = numeric(0),
  stringsAsFactors = FALSE
)

list_names <- names(directional_deg_lists)
for (i in 1:(length(list_names)-1)) {
  for (j in (i+1):length(list_names)) {
    list1_name <- list_names[i]
    list2_name <- list_names[j]
    
    list1_genes <- directional_deg_lists[[list1_name]]
    list2_genes <- directional_deg_lists[[list2_name]]
    
    # Calculate overlap
    overlap_genes <- intersect(list1_genes, list2_genes)
    union_genes <- union(list1_genes, list2_genes)
    
    # Jaccard index
    jaccard <- length(overlap_genes) / length(union_genes)
    
    # Add to summary
    overlap_summary <- rbind(overlap_summary, data.frame(
      List1 = list1_name,
      List2 = list2_name,
      List1_Size = length(list1_genes),
      List2_Size = length(list2_genes),
      Overlap_Size = length(overlap_genes),
      Jaccard_Index = round(jaccard, 3),
      stringsAsFactors = FALSE
    ))
  }
}

# Sort by overlap size
overlap_summary <- overlap_summary[order(overlap_summary$Overlap_Size, decreasing = TRUE), ]

# Display summary
cat("\nTop 10 DEG overlaps:\n")
print(head(overlap_summary, 10))

# Save overlap summary
write.csv(overlap_summary, "./output/reports/deg_overlap_summary.csv", row.names = FALSE)
cat("Overlap summary saved to ./output/reports/deg_overlap_summary.csv\n")

# ==============================================================================
# 4. Identify Key Overlap Patterns
# ==============================================================================

cat("Analyzing key overlap patterns...\n")

# Pattern 1: Same direction in different tissues (RET)
ret_tumor_up <- directional_deg_lists[["R0_vs_R1_tumor_UP"]]
ret_tumor_down <- directional_deg_lists[["R0_vs_R1_tumor_DOWN"]]
ret_normal_up <- directional_deg_lists[["R0_vs_R1_normal_UP"]]
ret_normal_down <- directional_deg_lists[["R0_vs_R1_normal_DOWN"]]

if (length(ret_tumor_up) > 0 && length(ret_normal_up) > 0) {
  ret_consistent_up <- intersect(ret_tumor_up, ret_normal_up)
  cat(sprintf("RET consistent upregulation (tumor + normal): %d genes\n", length(ret_consistent_up)))
}

if (length(ret_tumor_down) > 0 && length(ret_normal_down) > 0) {
  ret_consistent_down <- intersect(ret_tumor_down, ret_normal_down)
  cat(sprintf("RET consistent downregulation (tumor + normal): %d genes\n", length(ret_consistent_down)))
}

# Pattern 2: Opposite directions in different tissues (RET)
if (length(ret_tumor_up) > 0 && length(ret_normal_down) > 0) {
  ret_opposite_1 <- intersect(ret_tumor_up, ret_normal_down)
  cat(sprintf("RET opposite regulation (tumor↑, normal↓): %d genes\n", length(ret_opposite_1)))
}

if (length(ret_tumor_down) > 0 && length(ret_normal_up) > 0) {
  ret_opposite_2 <- intersect(ret_tumor_down, ret_normal_up)
  cat(sprintf("RET opposite regulation (tumor↓, normal↑): %d genes\n", length(ret_opposite_2)))
}

# Pattern 3: Driver-specific patterns (same tissue)
braf_normal_up <- directional_deg_lists[["B0_vs_B1_normal_UP"]]
braf_normal_down <- directional_deg_lists[["B0_vs_B1_normal_DOWN"]]

if (length(ret_normal_up) > 0 && length(braf_normal_up) > 0) {
  driver_specific_up <- setdiff(ret_normal_up, braf_normal_up)
  cat(sprintf("RET-specific upregulation (normal tissue): %d genes\n", length(driver_specific_up)))
}

if (length(ret_normal_down) > 0 && length(braf_normal_down) > 0) {
  driver_specific_down <- setdiff(ret_normal_down, braf_normal_down)
  cat(sprintf("RET-specific downregulation (normal tissue): %d genes\n", length(driver_specific_down)))
}

# ==============================================================================
# 5. Create Direction-Specific Heatmap
# ==============================================================================

cat("Creating direction-specific overlap heatmap...\n")

# Get all unique DEGs
all_unique_degs <- unique(unlist(directional_deg_lists))

if (length(all_unique_degs) > 0 && length(directional_deg_lists) > 0) {
  # Create binary matrix for heatmap
  direction_matrix <- matrix(0, 
                             nrow = length(all_unique_degs), 
                             ncol = length(directional_deg_lists))
  rownames(direction_matrix) <- all_unique_degs
  colnames(direction_matrix) <- names(directional_deg_lists)
  
  # Fill matrix
  for (list_name in names(directional_deg_lists)) {
    genes_in_list <- directional_deg_lists[[list_name]]
    direction_matrix[genes_in_list, list_name] <- 1
  }
  
  # Create heatmap only if reasonable size
  if (nrow(direction_matrix) <= 5000) {  # Limit for visualization
    library(pheatmap)
    
    # Create annotation for UP/DOWN
    col_annotation <- data.frame(
      Direction = ifelse(grepl("_UP$", colnames(direction_matrix)), "UP", "DOWN"),
      Comparison = gsub("_(UP|DOWN)$", "", colnames(direction_matrix))
    )
    rownames(col_annotation) <- colnames(direction_matrix)
    
    # Define colors
    ann_colors <- list(
      Direction = c("UP" = "red", "DOWN" = "blue"),
      Comparison = rainbow(length(unique(col_annotation$Comparison)))
    )
    names(ann_colors$Comparison) <- unique(col_annotation$Comparison)
    
    # Create and save heatmap
    pheatmap(direction_matrix,
             show_rownames = FALSE,
             annotation_col = col_annotation,
             annotation_colors = ann_colors,
             color = c("white", "black"),
             clustering_distance_cols = "binary",
             filename = "./output/plots/deg_direction_heatmap.png",
             width = 10, height = 8)
    
    cat("Direction-specific heatmap saved to ./output/plots/deg_direction_heatmap.png\n")
  } else {
    cat(sprintf("Too many DEGs (%d) for heatmap visualization, skipping...\n", nrow(direction_matrix)))
  }
}

# ==============================================================================
# 6. Export Key Gene Lists for Feature Selection
# ==============================================================================

cat("Exporting key gene lists for feature selection...\n")

# Create directory for gene lists
if (!dir.exists("./output/gene_lists")) {
  dir.create("./output/gene_lists", recursive = TRUE)
}

# Export all directional lists
for (list_name in names(directional_deg_lists)) {
  gene_list <- directional_deg_lists[[list_name]]
  if (length(gene_list) > 0) {
    write.table(gene_list, 
                file = paste0("./output/gene_lists/", list_name, ".txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

# Export specific pattern lists (if they exist)
pattern_lists <- list()

if (exists("ret_consistent_up") && length(ret_consistent_up) > 0) {
  pattern_lists[["RET_consistent_UP"]] <- ret_consistent_up
}
if (exists("ret_consistent_down") && length(ret_consistent_down) > 0) {
  pattern_lists[["RET_consistent_DOWN"]] <- ret_consistent_down
}
if (exists("ret_opposite_1") && length(ret_opposite_1) > 0) {
  pattern_lists[["RET_tumor_UP_normal_DOWN"]] <- ret_opposite_1
}
if (exists("ret_opposite_2") && length(ret_opposite_2) > 0) {
  pattern_lists[["RET_tumor_DOWN_normal_UP"]] <- ret_opposite_2
}

# Export pattern lists
for (pattern_name in names(pattern_lists)) {
  gene_list <- pattern_lists[[pattern_name]]
  write.table(gene_list,
              file = paste0("./output/gene_lists/pattern_", pattern_name, ".txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

cat("Gene lists exported to ./output/gene_lists/\n")

# ==============================================================================
# 7. Summary Report
# ==============================================================================

cat("\n==============================================\n")
cat("DEG Overlap Analysis Summary\n")
cat("==============================================\n")

total_unique_degs <- length(all_unique_degs)
total_directional_lists <- length(directional_deg_lists)

cat(sprintf("Total unique DEGs across all comparisons: %d\n", total_unique_degs))
cat(sprintf("Total directional lists analyzed: %d\n", total_directional_lists))

if (nrow(overlap_summary) > 0) {
  max_overlap <- max(overlap_summary$Overlap_Size)
  max_jaccard <- max(overlap_summary$Jaccard_Index)
  cat(sprintf("Largest overlap: %d genes\n", max_overlap))
  cat(sprintf("Highest Jaccard index: %.3f\n", max_jaccard))
  
  # Most overlapping pair
  top_overlap <- overlap_summary[1, ]
  cat(sprintf("Most overlapping pair: %s vs %s (%d genes, Jaccard=%.3f)\n",
              top_overlap$List1, top_overlap$List2, 
              top_overlap$Overlap_Size, top_overlap$Jaccard_Index))
}

cat("\nFiles created:\n")
cat("  ./output/plots/upset_directional_degs.png - UpSet plot\n")
cat("  ./output/plots/deg_direction_heatmap.png - Direction heatmap\n")
cat("  ./output/reports/deg_overlap_summary.csv - Overlap statistics\n")
cat("  ./output/gene_lists/*.txt - Individual gene lists\n")

cat("\nKey insights for feature selection:\n")
cat("  - Consistent patterns across tissues indicate robust targets\n")
cat("  - Opposite patterns suggest tissue-specific regulation\n")
cat("  - Driver-specific patterns show molecular specificity\n")
cat("  - Use overlap analysis to prioritize gene pairs for 12→1-2 reduction\n")

cat("DEG overlap visualization completed!\n")
cat("==============================================\n")

