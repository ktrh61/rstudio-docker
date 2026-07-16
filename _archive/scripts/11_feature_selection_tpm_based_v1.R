# ==============================================================================
# REBC-THYR Feature Selection Script - TPM Ratio Based v1 (Phase 1)
# 11_feature_selection_tpm_based_v1.R
# ==============================================================================
# Purpose: Basic data preparation and initial candidate count verification
# Strategy: Statistical performance absolute priority

# Required libraries
library(dplyr)

cat("Starting TPM ratio-based feature selection v1 (Phase 1: Basic verification)...\n")

# ==============================================================================
# 1. Load Data and DEG Results
# ==============================================================================

cat("Loading DEG v3 results and TPM data...\n")

# Load DEG analysis v3 results
if (!file.exists("./data/processed/final_deg_v3_results.rda")) {
  stop("DEG v3 results not found. Please run 09_deg_analysis_v3.R first.")
}
load("./data/processed/final_deg_v3_results.rda")

# Load TPM data
if (!file.exists("./data/processed/thyr_tpm.rda")) {
  stop("TPM data not found. Please run 02_tpm_calculation.R first.")
}
load("./data/processed/thyr_tpm.rda")

# Load high-purity sample lists
if (!file.exists("./data/processed/final_high_purity_sample_lists.rda")) {
  stop("High-purity sample lists not found.")
}
load("./data/processed/final_high_purity_sample_lists.rda")
high_purity_sample_lists <- final_high_purity_sample_lists

cat("Data loaded successfully\n")

# ==============================================================================
# 2. Extract R0_vs_R1_tumor DEG Results
# ==============================================================================

cat("Extracting R0_vs_R1_tumor DEG results...\n")

# Check if R0_vs_R1_tumor comparison is available
if (!"R0_vs_R1_tumor" %in% names(final_deg_v3_results$deg_analysis_results)) {
  stop("R0_vs_R1_tumor comparison not found in DEG results")
}

# Get tumor DEG results
tumor_deg_result <- final_deg_v3_results$deg_analysis_results[["R0_vs_R1_tumor"]]
tumor_degs <- tumor_deg_result$deg_summary$results_df

# Extract significant UP and DOWN genes
up_genes <- tumor_degs$gene_id[tumor_degs$significance_category == "Upregulated"]
down_genes <- tumor_degs$gene_id[tumor_degs$significance_category == "Downregulated"]

cat(sprintf("R0_vs_R1_tumor DEG summary:\n"))
cat(sprintf("  Total DEGs: %d\n", nrow(tumor_degs)))
cat(sprintf("  UP genes: %d\n", length(up_genes)))
cat(sprintf("  DOWN genes: %d\n", length(down_genes)))
cat(sprintf("  Total possible pairs: %d\n", length(up_genes) * length(down_genes)))

# ==============================================================================
# 3. Prepare TPM Data for Analysis
# ==============================================================================

cat("Preparing TPM data for analysis...\n")

# Get R0 and R1 tumor samples
r0_samples <- high_purity_sample_lists$R0$tumor
r1_samples <- high_purity_sample_lists$R1$tumor

cat(sprintf("Sample counts:\n"))
cat(sprintf("  R0 (ERR=0%%): %d samples\n", length(r0_samples)))
cat(sprintf("  R1 (ERR>=66.6%%): %d samples\n", length(r1_samples)))

# Verify samples are in TPM data
r0_available <- r0_samples[r0_samples %in% colnames(tpm)]
r1_available <- r1_samples[r1_samples %in% colnames(tpm)]

cat(sprintf("Available in TPM data:\n"))
cat(sprintf("  R0: %d/%d samples\n", length(r0_available), length(r0_samples)))
cat(sprintf("  R1: %d/%d samples\n", length(r1_available), length(r1_samples)))

if (length(r0_available) < 10 || length(r1_available) < 10) {
  stop("Insufficient samples available in TPM data")
}

# Create analysis TPM matrix
analysis_samples <- c(r0_available, r1_available)
analysis_tpm <- tpm[, analysis_samples]

# Create outcome vector (0=R0, 1=R1)
outcome <- c(rep(0, length(r0_available)), rep(1, length(r1_available)))
names(outcome) <- analysis_samples

cat(sprintf("Analysis dataset prepared:\n"))
cat(sprintf("  Total samples: %d (R0=%d, R1=%d)\n", 
            length(analysis_samples), length(r0_available), length(r1_available)))
cat(sprintf("  Total genes in TPM: %d\n", nrow(analysis_tpm)))

# ==============================================================================
# 4. Check Gene Availability in TPM Data
# ==============================================================================

cat("Checking gene availability in TPM data...\n")

# Function to clean gene IDs (remove version numbers)
clean_gene_id <- function(gene_ids) {
  return(gsub("\\.\\d+$", "", gene_ids))
}

# Clean gene IDs for matching
up_genes_clean <- clean_gene_id(up_genes)
down_genes_clean <- clean_gene_id(down_genes)
tpm_genes_clean <- clean_gene_id(rownames(analysis_tpm))

# Find available genes
up_genes_avail_clean <- up_genes_clean[up_genes_clean %in% tpm_genes_clean]
down_genes_avail_clean <- down_genes_clean[down_genes_clean %in% tpm_genes_clean]

# Map back to original IDs with versions
up_genes_avail <- character(0)
for (clean_id in up_genes_avail_clean) {
  # Find the first matching original gene ID
  matching_original <- up_genes[clean_gene_id(up_genes) == clean_id][1]
  matching_tpm <- rownames(analysis_tpm)[tpm_genes_clean == clean_id][1]
  if (!is.na(matching_original) && !is.na(matching_tpm)) {
    up_genes_avail <- c(up_genes_avail, matching_tpm)
  }
}

down_genes_avail <- character(0)
for (clean_id in down_genes_avail_clean) {
  matching_original <- down_genes[clean_gene_id(down_genes) == clean_id][1]
  matching_tpm <- rownames(analysis_tpm)[tpm_genes_clean == clean_id][1]
  if (!is.na(matching_original) && !is.na(matching_tpm)) {
    down_genes_avail <- c(down_genes_avail, matching_tpm)
  }
}

cat(sprintf("Gene availability check:\n"))
cat(sprintf("  UP genes: %d/%d available in TPM (%.1f%%)\n", 
            length(up_genes_avail), length(up_genes), 
            length(up_genes_avail)/length(up_genes)*100))
cat(sprintf("  DOWN genes: %d/%d available in TPM (%.1f%%)\n", 
            length(down_genes_avail), length(down_genes), 
            length(down_genes_avail)/length(down_genes)*100))
cat(sprintf("  Possible pairs with available genes: %d\n", 
            length(up_genes_avail) * length(down_genes_avail)))

# ==============================================================================
# 5. Basic Statistics Check
# ==============================================================================

cat("Performing basic statistics check...\n")

# Check for zero expression
up_zero_count <- sum(rowSums(analysis_tpm[up_genes_avail, ] == 0) > 0)
down_zero_count <- sum(rowSums(analysis_tpm[down_genes_avail, ] == 0) > 0)

cat(sprintf("Zero expression check:\n"))
cat(sprintf("  UP genes with any zero values: %d/%d\n", up_zero_count, length(up_genes_avail)))
cat(sprintf("  DOWN genes with any zero values: %d/%d\n", down_zero_count, length(down_genes_avail)))

# Sample one pair for basic calculation test
if (length(up_genes_avail) > 0 && length(down_genes_avail) > 0) {
  test_up_gene <- up_genes_avail[1]
  test_down_gene <- down_genes_avail[1]
  
  cat(sprintf("\nTesting basic calculation with sample pair:\n"))
  cat(sprintf("  UP gene: %s\n", test_up_gene))
  cat(sprintf("  DOWN gene: %s\n", test_down_gene))
  
  # Extract TPM values
  test_up_tpm <- as.numeric(analysis_tpm[test_up_gene, ])
  test_down_tpm <- as.numeric(analysis_tpm[test_down_gene, ])
  
  # Calculate ratio with pseudocount
  pseudocount <- 0.1
  test_ratios <- log2((test_up_tpm + pseudocount) / (test_down_tpm + pseudocount))
  
  # Basic statistics
  r0_indices <- which(outcome == 0)
  r1_indices <- which(outcome == 1)
  
  test_r0_ratios <- test_ratios[r0_indices]
  test_r1_ratios <- test_ratios[r1_indices]
  
  cat(sprintf("  R0 ratios: mean=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
              mean(test_r0_ratios), median(test_r0_ratios), 
              min(test_r0_ratios), max(test_r0_ratios)))
  cat(sprintf("  R1 ratios: mean=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
              mean(test_r1_ratios), median(test_r1_ratios), 
              min(test_r1_ratios), max(test_r1_ratios)))
  
  # Check direction
  r0_median <- median(test_r0_ratios)
  r1_median <- median(test_r1_ratios)
  direction_ok <- r0_median != r1_median
  separation <- abs(r1_median - r0_median)
  
  cat(sprintf("  Direction check: %s (separation=%.3f)\n", 
              ifelse(direction_ok, "OK", "FAILED"), separation))
  
  # Basic ERR=0 specificity check
  if (r1_median > r0_median) {
    threshold_est <- max(test_r0_ratios) + 0.1
    err0_spec <- sum(test_r0_ratios <= threshold_est) / length(test_r0_ratios)
  } else {
    threshold_est <- min(test_r0_ratios) - 0.1
    err0_spec <- sum(test_r0_ratios >= threshold_est) / length(test_r0_ratios)
  }
  
  cat(sprintf("  Estimated ERR=0 specificity: %.1f%%\n", err0_spec * 100))
  
} else {
  cat("⚠️  No genes available for testing!\n")
}

# ==============================================================================
# 6. Prepare Index Arrays
# ==============================================================================

cat("Preparing analysis indices...\n")

# Create index arrays for efficient computation
r0_indices <- which(outcome == 0)
r1_indices <- which(outcome == 1)

cat(sprintf("Analysis indices prepared:\n"))
cat(sprintf("  R0 indices: %d samples (positions %s)\n", 
            length(r0_indices), paste(range(r0_indices), collapse="-")))
cat(sprintf("  R1 indices: %d samples (positions %s)\n", 
            length(r1_indices), paste(range(r1_indices), collapse="-")))

# ==============================================================================
# 7. Save Phase 1 Results
# ==============================================================================

cat("Saving Phase 1 results...\n")

if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Phase 1 results
phase1_results <- list(
  # Data info
  total_degs = nrow(tumor_degs),
  up_genes_total = length(up_genes),
  down_genes_total = length(down_genes),
  up_genes_available = length(up_genes_avail),
  down_genes_available = length(down_genes_avail),
  total_possible_pairs = length(up_genes_avail) * length(down_genes_avail),
  
  # Sample info
  r0_samples_count = length(r0_available),
  r1_samples_count = length(r1_available),
  total_samples = length(analysis_samples),
  
  # Gene lists
  up_genes_avail = up_genes_avail,
  down_genes_avail = down_genes_avail,
  
  # Analysis data
  analysis_samples = analysis_samples,
  outcome = outcome,
  r0_indices = r0_indices,
  r1_indices = r1_indices,
  
  # Test results
  test_pair = if(exists("test_up_gene")) c(test_up_gene, test_down_gene) else NULL,
  test_separation = if(exists("separation")) separation else NULL,
  test_err0_specificity = if(exists("err0_spec")) err0_spec else NULL,
  
  analysis_date = Sys.time(),
  analysis_version = "phase1_basic_verification"
)

save(phase1_results, file = "./data/processed/feature_selection_phase1_results.rda")

# Also save the analysis TPM matrix for next phases
save(analysis_tpm, file = "./data/processed/analysis_tpm_matrix.rda")

cat("Phase 1 results saved successfully\n")

# ==============================================================================
# 8. Summary and Next Steps
# ==============================================================================

cat("\n==============================================\n")
cat("Feature Selection Phase 1 Summary\n")
cat("==============================================\n")

cat("Data Preparation:\n")
cat(sprintf("  ✅ DEG results loaded: %d total DEGs\n", nrow(tumor_degs)))
cat(sprintf("  ✅ UP genes: %d (%d available in TPM)\n", length(up_genes), length(up_genes_avail)))
cat(sprintf("  ✅ DOWN genes: %d (%d available in TPM)\n", length(down_genes), length(down_genes_avail)))
cat(sprintf("  ✅ Samples: R0=%d, R1=%d (total=%d)\n", 
            length(r0_available), length(r1_available), length(analysis_samples)))

cat("\nCandidate Pairs:\n")
cat(sprintf("  📊 Total possible pairs: %s\n", 
            format(length(up_genes_avail) * length(down_genes_avail), big.mark=",")))

if (exists("separation") && exists("err0_spec")) {
  cat("\nSample Pair Test:\n")
  cat(sprintf("  📈 Separation: %.3f log2 units\n", separation))
  cat(sprintf("  🎯 ERR=0 specificity: %.1f%%\n", err0_spec * 100))
  
  if (separation >= log2(1.5) && err0_spec >= 0.90) {
    cat("  ✅ Sample pair shows promising performance\n")
  } else {
    cat("  ⚠️  Sample pair shows modest performance\n")
  }
}

cat("\nNext Steps:\n")
if (length(up_genes_avail) > 0 && length(down_genes_avail) > 0) {
  cat("  1. ✅ Ready for Phase 2: FoldChange filtering tests\n")
  cat("  2. Run different FoldChange control patterns\n")
  cat("  3. Check candidate counts for each pattern\n")
  cat("  4. Select optimal filtering strategy\n")
} else {
  cat("  ❌ No genes available - check DEG results and TPM data\n")
}

cat("\nPhase 1 completed successfully!\n")
cat("==============================================\n")

