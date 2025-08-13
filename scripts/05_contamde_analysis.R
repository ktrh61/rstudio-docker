# ==============================================================================
# REBC-THYR ContamDE Analysis Script
# 05_contamde_analysis.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(edgeR)
library(dplyr)

# Source the contamDE function
source("./utils/contamde_functions.R")

cat("Starting ContamDE analysis for tumor purity estimation...\n")

# ==============================================================================
# 1. Load Data and Sample Lists
# ==============================================================================

cat("Loading data and sample lists...\n")

# Load sample lists from previous step
load("./data/processed/sample_lists.rda")

# Load SummarizedExperiment object (assuming se_thyr is available)
# If not loaded, uncomment the following:
# load("./data/raw/thyr_data.rda")
# se_thyr <- data
# rm(data)

# Get stranded_second count data for ContamDE analysis
count_data_full <- assay(se_thyr, "stranded_second")
gene_info <- rowData(se_thyr)

cat("Original count data dimensions:", dim(count_data_full), "\n")

# ==============================================================================
# 2. Filter to Protein-Coding Genes
# ==============================================================================

cat("Filtering to protein-coding genes...\n")

# Filter to protein-coding genes only
protein_coding_genes <- gene_info$gene_type == "protein_coding"
count_data_pc <- count_data_full[protein_coding_genes, ]
gene_info_pc <- gene_info[protein_coding_genes, ]

cat("After protein-coding filter:", dim(count_data_pc), "\n")

# ==============================================================================
# 3. Process Each Group Separately
# ==============================================================================

cat("Processing each group for ContamDE analysis...\n")

# Initialize results storage
contamde_results <- list()
filtered_sample_lists <- list()
purity_cutoff <- 0.6  # Absolute purity threshold (60%)

for (group in c("R0", "R1", "B0", "B1")) {
  cat(sprintf("\n--- Processing group %s ---\n", group))
  
  # Check if group has samples
  if (length(sample_lists[[group]]$tumor) == 0) {
    cat(sprintf("Group %s has no samples, skipping...\n", group))
    contamde_results[[group]] <- NULL
    filtered_sample_lists[[group]] <- list(tumor = character(0), 
                                           normal = character(0), 
                                           cases = character(0))
    next
  }
  
  # Get tumor and normal samples for this group
  tumor_samples <- sample_lists[[group]]$tumor
  normal_samples <- sample_lists[[group]]$normal
  
  # Verify sample availability in count data
  tumor_available <- tumor_samples %in% colnames(count_data_pc)
  normal_available <- normal_samples %in% colnames(count_data_pc)
  
  if (!all(tumor_available) || !all(normal_available)) {
    missing_tumor <- sum(!tumor_available)
    missing_normal <- sum(!normal_available)
    cat(sprintf("Warning: Missing samples - tumor: %d, normal: %d\n", 
                missing_tumor, missing_normal))
  }
  
  # Use only available samples
  tumor_samples_avail <- tumor_samples[tumor_available]
  normal_samples_avail <- normal_samples[normal_available]
  
  # Ensure equal number of tumor and normal samples
  min_pairs <- min(length(tumor_samples_avail), length(normal_samples_avail))
  if (min_pairs == 0) {
    cat(sprintf("Group %s has no valid pairs, skipping...\n", group))
    contamde_results[[group]] <- NULL
    filtered_sample_lists[[group]] <- list(tumor = character(0), 
                                           normal = character(0), 
                                           cases = character(0))
    next
  }
  
  tumor_samples_final <- tumor_samples_avail[1:min_pairs]
  normal_samples_final <- normal_samples_avail[1:min_pairs]
  
  # Prepare count data in [normal1, normal2, ..., tumor1, tumor2, ...] format
  group_samples <- c(normal_samples_final, tumor_samples_final)
  group_counts <- count_data_pc[, group_samples]
  
  cat(sprintf("Group %s: %d pairs (%d samples total)\n", 
              group, min_pairs, ncol(group_counts)))
  
  # Apply filterByExpr for low expression filtering
  dge <- DGEList(counts = group_counts)
  # Create group factor for filterByExpr (normal=0, tumor=1)
  group_factor <- factor(c(rep("normal", min_pairs), rep("tumor", min_pairs)))
  keep <- filterByExpr(dge, group = group_factor, min.count = 1, min.total.count = 15)
  dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
  
  cat(sprintf("After filterByExpr: %d genes (from %d)\n", 
              nrow(dge_filtered), nrow(group_counts)))
  
  # Run ContamDE analysis
  cat("Running ContamDE analysis...\n")
  
  tryCatch({
    contamde_result <- contamde_lm(
      counts = dge_filtered$counts,
      subtype = NULL,
      covariate = NULL,
      contaminated = TRUE,
      robust = TRUE
    )
    
    # Extract purity estimates
    purity_estimates <- contamde_result$proportion
    names(purity_estimates) <- paste0("Pair_", 1:length(purity_estimates))
    
    cat("Purity estimates summary:\n")
    print(summary(purity_estimates))
    
    # Apply absolute purity cutoff (proportion > 0.6)
    purity_threshold <- purity_cutoff  # Use absolute threshold
    high_purity_pairs <- which(purity_estimates > purity_threshold)
    
    cat(sprintf("Purity threshold (absolute): %.1f\n", purity_threshold))
    cat(sprintf("High purity pairs (>%.1f): %d/%d\n", 
                purity_threshold, length(high_purity_pairs), length(purity_estimates)))
    
    # Update sample lists with high purity samples only
    if (length(high_purity_pairs) > 0) {
      filtered_tumor <- tumor_samples_final[high_purity_pairs]
      filtered_normal <- normal_samples_final[high_purity_pairs]
      filtered_cases <- sample_lists[[group]]$cases[high_purity_pairs]
      
      filtered_sample_lists[[group]] <- list(
        tumor = filtered_tumor,
        normal = filtered_normal,
        cases = filtered_cases
      )
    } else {
      cat("Warning: No samples passed purity filter!\n")
      filtered_sample_lists[[group]] <- list(
        tumor = character(0),
        normal = character(0),
        cases = character(0)
      )
    }
    
    # Store results
    contamde_results[[group]] <- list(
      result = contamde_result,
      purity_estimates = purity_estimates,
      purity_threshold = purity_threshold,
      high_purity_pairs = high_purity_pairs,
      filtered_genes = rownames(dge_filtered),
      original_pairs = min_pairs,
      final_pairs = length(high_purity_pairs)
    )
    
    cat(sprintf("Group %s processing completed successfully\n", group))
    
  }, error = function(e) {
    cat(sprintf("Error processing group %s: %s\n", group, e$message))
    contamde_results[[group]] <- NULL
    filtered_sample_lists[[group]] <- list(tumor = character(0), 
                                           normal = character(0), 
                                           cases = character(0))
  })
}

# ==============================================================================
# 4. Summary and Quality Control
# ==============================================================================

cat("\n==============================================\n")
cat("ContamDE Analysis Summary\n")
cat("==============================================\n")

# Create summary table
summary_table <- data.frame(
  Group = c("R0", "R1", "B0", "B1"),
  Original_Pairs = integer(4),
  Final_Pairs = integer(4),
  Reduction = character(4),
  Mean_Purity = numeric(4),
  Purity_Threshold = numeric(4),
  stringsAsFactors = FALSE
)

for (i in 1:4) {
  group <- summary_table$Group[i]
  
  if (!is.null(contamde_results[[group]])) {
    original <- contamde_results[[group]]$original_pairs
    final <- contamde_results[[group]]$final_pairs
    purity_est <- contamde_results[[group]]$purity_estimates
    threshold <- contamde_results[[group]]$purity_threshold
    
    summary_table$Original_Pairs[i] <- original
    summary_table$Final_Pairs[i] <- final
    summary_table$Reduction[i] <- sprintf("%.1f%%", (1 - final/original) * 100)
    summary_table$Mean_Purity[i] <- mean(purity_est, na.rm = TRUE)
    summary_table$Purity_Threshold[i] <- purity_cutoff  # Always 0.6
  } else {
    summary_table$Original_Pairs[i] <- 0
    summary_table$Final_Pairs[i] <- 0
    summary_table$Reduction[i] <- "N/A"
    summary_table$Mean_Purity[i] <- NA
    summary_table$Purity_Threshold[i] <- purity_cutoff  # Always 0.6
  }
}

print(summary_table)

# Check minimum sample requirements
min_samples <- 5
sufficient_groups <- summary_table$Final_Pairs >= min_samples
cat("\nSample size check after purity filtering (minimum 5 per group):\n")
for (i in 1:4) {
  group <- summary_table$Group[i]
  status <- ifelse(sufficient_groups[i], "✅ PASS", "❌ INSUFFICIENT")
  cat(sprintf("  %s: %d pairs %s\n", group, summary_table$Final_Pairs[i], status))
}

total_final_pairs <- sum(summary_table$Final_Pairs)
cat(sprintf("\nTotal pairs after purity filtering: %d\n", total_final_pairs))

if (all(sufficient_groups)) {
  cat("✅ All groups meet minimum sample requirements\n")
} else {
  insufficient_groups <- summary_table$Group[!sufficient_groups]
  cat(sprintf("⚠️  Insufficient samples in: %s\n", 
              paste(insufficient_groups, collapse = ", ")))
}

# ==============================================================================
# 5. Save Results
# ==============================================================================

cat("Saving ContamDE results...\n")

# Create processed directory if it doesn't exist
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save filtered sample lists
save(filtered_sample_lists, file = "./data/processed/filtered_sample_lists.rda")

# Save detailed ContamDE results
contamde_analysis_results <- list(
  contamde_results = contamde_results,
  filtered_sample_lists = filtered_sample_lists,
  summary_table = summary_table,
  purity_cutoff = purity_cutoff,
  protein_coding_genes = rownames(count_data_pc),
  analysis_date = Sys.time()
)

save(contamde_analysis_results, file = "./data/processed/contamde_analysis_results.rda")

cat("Results saved to ./data/processed/\n")

# ==============================================================================
# 6. Prepare for Next Analysis
# ==============================================================================

cat("\nPreparing for next analysis steps...\n")

# Check which groups are viable for DEG analysis
viable_groups <- summary_table$Group[sufficient_groups]
cat("Groups viable for DEG analysis:", paste(viable_groups, collapse = ", "), "\n")

# Define possible comparisons (tissue-specific)
possible_comparisons <- list()
if ("R0" %in% viable_groups && "R1" %in% viable_groups) {
  possible_comparisons$R0_vs_R1_tumor <- "RET Tumor: Unexposed vs RadHigh (PRIMARY)"
  possible_comparisons$R0_vs_R1_normal <- "RET Normal: Unexposed vs RadHigh"
}
if ("B0" %in% viable_groups && "B1" %in% viable_groups) {
  possible_comparisons$B0_vs_B1_tumor <- "BRAF Tumor: Unexposed vs RadHigh"
  possible_comparisons$B0_vs_B1_normal <- "BRAF Normal: Unexposed vs RadHigh"
}

cat("Possible DEG comparisons:\n")
for (comp_name in names(possible_comparisons)) {
  cat(sprintf("  %s: %s\n", comp_name, possible_comparisons[[comp_name]]))
}

# Store comparison info
comparison_info <- list(
  viable_groups = viable_groups,
  possible_comparisons = possible_comparisons,
  primary_comparison = ifelse("R0_vs_R1_tumor" %in% names(possible_comparisons), 
                              "R0_vs_R1_tumor", "None")
)

save(comparison_info, file = "./data/processed/comparison_info.rda")

cat("\nContamDE analysis completed successfully!\n")
cat("Next step: PCA analysis for outlier detection\n")
cat("==============================================\n")