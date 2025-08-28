# ==============================================================================
# REBC-THYR TPM Calculation Script v6 - Simplified Gene-Level TPM
# 02_tpm_calculation_v6.R
# ==============================================================================
# Purpose: Calculate TPM using simple gene length normalization
# Changes from original: Removed effective length calculation (inappropriate for ENSG)

# Required libraries
library(SummarizedExperiment)
library(GenomicFeatures)
library(GenomicRanges)
library(AnnotationDbi)
library(txdbmaker)
library(dplyr)

cat("Starting TPM calculation v6 with simplified gene-length normalization...\n")
cat("Note: This version uses union exon length without effective length adjustment\n")

# ==============================================================================
# 1. Load Data
# ==============================================================================

cat("\nLoading REBC-THYR data...\n")
load("./data/raw/thyr_data.rda")

# Rename data object for clarity
se_thyr <- data
rm(data)

cat("Data loaded successfully!\n")
cat("Data class:", class(se_thyr), "\n")
cat("Dimensions:", dim(se_thyr), "\n")
cat("Available assays:", paste(assayNames(se_thyr), collapse = ", "), "\n")

# Extract count data (stranded_second based on reverse strand protocol)
count_data <- assay(se_thyr, "stranded_second")
cat("Count data dimensions:", dim(count_data), "\n")

# ==============================================================================
# 2. TxDb Creation or Loading
# ==============================================================================

cat("\nPreparing gene annotation...\n")

if (file.exists("./data/raw/gencode.v36.annotation.txdb.sqlite")) {
  cat("Loading existing TxDb...\n")
  txdb <- AnnotationDbi::loadDb("./data/raw/gencode.v36.annotation.txdb.sqlite")
} else {
  cat("Creating TxDb from GTF file...\n")
  if (!file.exists("./data/raw/gencode.v36.annotation.gtf")) {
    stop("GTF file not found in ./data/raw/gencode.v36.annotation.gtf")
  }
  txdb <- txdbmaker::makeTxDbFromGFF("./data/raw/gencode.v36.annotation.gtf")
  AnnotationDbi::saveDb(txdb, file = "./data/raw/gencode.v36.annotation.txdb.sqlite")
  cat("TxDb saved for future use\n")
}

# ==============================================================================
# 3. Gene Length Calculation
# ==============================================================================

cat("\nCalculating gene lengths (union exon model)...\n")

# Get exons grouped by gene
exons_by_gene <- GenomicFeatures::exonsBy(txdb, by = "gene")

# Calculate total exon length per gene (union of all exons)
# reduce() merges overlapping exons to avoid double-counting
gene_lengths <- sum(width(GenomicRanges::reduce(exons_by_gene)))
names(gene_lengths) <- names(exons_by_gene)

cat("Number of genes with length annotation:", length(gene_lengths), "\n")
cat("Gene length summary:\n")
print(summary(gene_lengths))

# Check for problematic gene lengths
zero_length_genes <- sum(gene_lengths == 0)
if (zero_length_genes > 0) {
  cat(sprintf("Warning: %d genes with zero length detected\n", zero_length_genes))
  gene_lengths[gene_lengths == 0] <- 1  # Set minimum length to 1
}

# ==============================================================================
# 4. Match Genes Between Counts and Annotation
# ==============================================================================

cat("\nMatching genes between count data and annotation...\n")

# Get common genes
gene_ids_counts <- rownames(count_data)
gene_ids_lengths <- names(gene_lengths)
common_genes <- intersect(gene_ids_counts, gene_ids_lengths)

cat(sprintf("Genes in count data: %d\n", length(gene_ids_counts)))
cat(sprintf("Genes with length annotation: %d\n", length(gene_ids_lengths)))
cat(sprintf("Common genes: %d\n", length(common_genes)))

if (length(common_genes) == 0) {
  stop("No common genes found between count data and gene lengths")
}

# Filter to common genes
count_data_filtered <- count_data[common_genes, ]
gene_lengths_filtered <- gene_lengths[common_genes]

# Ensure same order
count_data_filtered <- count_data_filtered[names(gene_lengths_filtered), ]

cat("Final dimensions for TPM calculation:", dim(count_data_filtered), "\n")

# ==============================================================================
# 5. TPM Calculation
# ==============================================================================

cat("\nCalculating TPM values...\n")

# Step 1: Calculate RPK (Reads Per Kilobase)
# Divide each gene's count by its length in kilobases
rpk <- count_data_filtered / (gene_lengths_filtered / 1000)

# Check for any infinite or NA values
if (any(!is.finite(rpk))) {
  problem_count <- sum(!is.finite(rpk))
  cat(sprintf("Warning: %d non-finite RPK values detected, setting to 0\n", problem_count))
  rpk[!is.finite(rpk)] <- 0
}

# Step 2: Calculate scaling factor (sum of RPK per sample)
scaling_factors <- colSums(rpk)

# Check scaling factors
if (any(scaling_factors == 0)) {
  zero_samples <- sum(scaling_factors == 0)
  stop(sprintf("Error: %d samples have zero total RPK", zero_samples))
}

cat("Scaling factors summary:\n")
print(summary(scaling_factors))

# Step 3: Calculate TPM
tpm <- t(t(rpk) / scaling_factors) * 1e6

# ==============================================================================
# 6. Quality Control
# ==============================================================================

cat("\nPerforming quality control checks...\n")

# Check TPM properties
tpm_sums <- colSums(tpm)
cat("TPM column sums (should be ~1,000,000):\n")
print(summary(tpm_sums))

# Verify all sums are close to 1 million
if (any(abs(tpm_sums - 1e6) > 1)) {
  cat("Warning: Some samples have TPM sums deviating from 1,000,000\n")
}

# Check for negative values
negative_count <- sum(tpm < 0, na.rm = TRUE)
cat(sprintf("Number of negative TPM values: %d\n", negative_count))

# Check for NA values
na_count <- sum(is.na(tpm))
cat(sprintf("Number of NA TPM values: %d\n", na_count))

# Overall statistics
cat("\nTPM matrix statistics:\n")
cat(sprintf("  Dimensions: %d genes × %d samples\n", nrow(tpm), ncol(tpm)))
cat(sprintf("  Non-zero values: %d (%.1f%%)\n", 
            sum(tpm > 0), sum(tpm > 0) / length(tpm) * 100))
cat(sprintf("  Median TPM (non-zero): %.2f\n", median(tpm[tpm > 0])))

# Gene expression distribution
gene_max_tpm <- apply(tpm, 1, max)
expressed_genes <- sum(gene_max_tpm > 1)
cat(sprintf("  Genes with max TPM > 1: %d (%.1f%%)\n", 
            expressed_genes, expressed_genes / nrow(tpm) * 100))

# ==============================================================================
# 7. Save Results
# ==============================================================================

cat("\nSaving TPM data...\n")

# Create processed directory if it doesn't exist
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save TPM matrix
tpm_v6 <- tpm
save(tpm_v6, file = "./data/processed/tpm_v6.rda")
cat("TPM data saved to ./data/processed/tpm_v6.rda\n")

# Save metadata
tpm_v6_info <- list(
  calculation_date = Sys.time(),
  method = "simple_gene_length",
  description = "TPM calculated using union exon length per gene without effective length correction",
  n_genes = nrow(tpm),
  n_samples = ncol(tpm),
  gene_ids = rownames(tpm),
  sample_ids = colnames(tpm),
  gene_lengths = gene_lengths_filtered,
  scaling_factors = scaling_factors,
  version = "v6"
)

save(tpm_v6_info, file = "./data/processed/tpm_v6_info.rda")
cat("TPM metadata saved to ./data/processed/tpm_v6_info.rda\n")

# Optional: Save as CSV for external use
write.csv(tpm, file = "./data/processed/tpm_v6.csv", row.names = TRUE)
cat("TPM matrix also saved as CSV: ./data/processed/tpm_v6.csv\n")

# ==============================================================================
# 8. Comparison with Original TPM (Optional)
# ==============================================================================

if (file.exists("./data/processed/thyr_tpm.rda")) {
  cat("\nComparing with original TPM calculation...\n")
  load("./data/processed/thyr_tpm.rda")
  
  # Compare dimensions
  cat(sprintf("Original TPM: %d genes × %d samples\n", nrow(tpm), ncol(tpm)))
  cat(sprintf("New TPM v6: %d genes × %d samples\n", nrow(tpm_v6), ncol(tpm_v6)))
  
  # Compare common genes
  common_compare <- intersect(rownames(tpm), rownames(tpm_v6))
  if (length(common_compare) > 0) {
    # Sample correlation
    if (ncol(tpm) == ncol(tpm_v6)) {
      sample_cors <- sapply(1:ncol(tpm), function(i) {
        cor(tpm[common_compare, i], tpm_v6[common_compare, i], method = "spearman")
      })
      cat("Sample-wise correlation between old and new TPM:\n")
      print(summary(sample_cors))
    }
  }
}

# ==============================================================================
# 9. Final Summary
# ==============================================================================

cat("\n==============================================\n")
cat("TPM Calculation v6 Summary\n")
cat("==============================================\n")

cat("Calculation complete using simplified gene-length normalization\n")
cat("Key features:\n")
cat("  - Union exon length per gene (no effective length)\n")
cat("  - Suitable for gene-level (ENSG) analysis\n")
cat("  - Compatible with PCR/Western blot validation\n")

cat(sprintf("\nFinal TPM matrix: %d genes × %d samples\n", nrow(tpm_v6), ncol(tpm_v6)))
cat(sprintf("Expressed genes (max TPM > 1): %d\n", expressed_genes))

cat("\nNext steps:\n")
cat("1. Use tpm_v6 for feature selection (11_feature_selection_v6.R)\n")
cat("2. Apply 3-filter strategy for gene pair identification\n")
cat("3. Focus on R0_vs_R1_tumor DEGs (1,140 candidates)\n")

cat("\nTPM calculation v6 completed successfully!\n")
cat("==============================================\n")