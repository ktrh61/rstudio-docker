# 03_prepare_counts.R - Extract stranded_second and remove zero-count genes
# Purpose: Extract stranded_second assay and remove completely unexpressed genes
# Input: thyr_se_raw.rds (full SummarizedExperiment)
# Output: thyr_se_strand2_nonzero.rds (SE with stranded_second only, zero genes removed, gene lengths added)
# Version: v1.1 - Added gene length from GENCODE v36

source("analysis_v7/setup.R")

cat("\n=== Extract stranded_second and Remove Zero-count Genes ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# ============================================================================
# Load SummarizedExperiment
# ============================================================================
cat("\n--- Loading SummarizedExperiment ---\n")

se_path <- paste0(paths$processed, "thyr_se_raw.rds")
if (!file.exists(se_path)) {
  stop("thyr_se_raw.rds not found. Please run 02_data_loading.R first.")
}

thyr_se <- readRDS(se_path)
cat("Loaded SE:", nrow(thyr_se), "genes x", ncol(thyr_se), "samples\n")

# Verify stranded_second is the appropriate assay
cat("\n--- Verifying assay selection ---\n")
total_counts <- c(
  unstranded = sum(assay(thyr_se, "unstranded")),
  stranded_first = sum(assay(thyr_se, "stranded_first")),
  stranded_second = sum(assay(thyr_se, "stranded_second"))
)
cat("Total counts by assay:\n")
print(total_counts)

if (names(which.max(total_counts)) != "stranded_second") {
  warning("stranded_second is not the highest count assay!")
}

# ============================================================================
# Extract stranded_second counts for processing
# ============================================================================
cat("\n--- Extracting stranded_second counts ---\n")

counts <- assay(thyr_se, "stranded_second")
cat("Count matrix dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")

# Basic statistics before filtering
cat("\nBasic statistics (before filtering):\n")
cat("  Total counts:", format(sum(counts), big.mark = ","), "\n")
cat("  Non-zero entries:", format(sum(counts > 0), big.mark = ","), "\n")
cat("  Sparsity:", sprintf("%.2f%%", sum(counts == 0) / length(counts) * 100), "\n")

# ============================================================================
# Gene filtering - Zero count removal only
# ============================================================================
cat("\n--- Gene filtering ---\n")

# Filter: Remove genes with zero counts across all samples
zero_genes <- rowSums(counts) == 0
cat("Genes with zero counts across all samples:", sum(zero_genes), "\n")
cat("Genes to keep:", sum(!zero_genes), "\n")

# Apply filter
keep_genes <- !zero_genes
cat("\n--- Filter summary ---\n")
cat("  Starting genes:", nrow(counts), "\n")
cat("  After zero removal:", sum(keep_genes), "\n")
cat("  Retention rate:", sprintf("%.1f%%", sum(keep_genes) / nrow(counts) * 100), "\n")

# ============================================================================
# Create filtered SummarizedExperiment
# ============================================================================
cat("\n--- Creating filtered SummarizedExperiment ---\n")

# Create new SE with only stranded_second and non-zero genes
thyr_se_strand2_nonzero <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = counts[keep_genes, , drop = FALSE]),
  colData = colData(thyr_se),
  rowData = rowData(thyr_se)[keep_genes, , drop = FALSE]
)

cat("Filtered SE created:\n")
cat("  Genes:", nrow(thyr_se_strand2_nonzero), "\n")
cat("  Samples:", ncol(thyr_se_strand2_nonzero), "\n")
cat("  Assays: counts (stranded_second only)\n")
cat("  Gene metadata preserved: gene_id, gene_name, gene_type\n")

# ============================================================================
# Save results
# ============================================================================
cat("\n--- Saving results ---\n")

# Save filtered SummarizedExperiment (first save)
output_file <- paste0(paths$processed, "thyr_se_strand2_nonzero.rds")
saveRDS(thyr_se_strand2_nonzero, output_file)
cat("Saved filtered SE to:", output_file, "\n")

# Save gene filter information for reference
filter_info <- data.frame(
  gene_id = rownames(thyr_se),
  gene_name = rowData(thyr_se)$gene_name,
  gene_type = rowData(thyr_se)$gene_type,
  zero_count = zero_genes,
  kept = keep_genes,
  stringsAsFactors = FALSE
)
filter_file <- paste0(paths$logs, "gene_filter_info_", 
                      format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
saveRDS(filter_info, filter_file)
cat("Saved filter information to logs/\n")

# ============================================================================
# Add gene lengths to rowData (for TPM calculation)
# ============================================================================

cat("\n--- Adding gene lengths from GENCODE v36 ---\n")

gene_length_added <- FALSE
n_genes_with_length <- 0

gene_lengths_path <- paste0(paths$processed, "gene_lengths.rds")
if (!file.exists(gene_lengths_path)) {
  cat("  Warning: gene_lengths.rds not found. Run 00_prepare_gene_lengths.R first.\n")
  cat("  Continuing without gene lengths.\n")
} else {
  gene_lengths <- readRDS(gene_lengths_path)
  
  # Match gene IDs (remove version numbers if present)
  gene_ids_clean <- gsub("\\..*", "", rownames(thyr_se_strand2_nonzero))
  matched_lengths <- gene_lengths[gene_ids_clean]
  
  # Report matching statistics
  n_genes_with_length <- sum(!is.na(matched_lengths))
  cat(sprintf("  Matched %d/%d genes (%.1f%%)\n", 
              n_genes_with_length, 
              nrow(thyr_se_strand2_nonzero),
              n_genes_with_length / nrow(thyr_se_strand2_nonzero) * 100))
  
  # Check for any unmatched genes
  if (n_genes_with_length < nrow(thyr_se_strand2_nonzero)) {
    unmatched <- gene_ids_clean[is.na(matched_lengths)]
    cat(sprintf("  Warning: %d genes without length information\n", 
                length(unmatched)))
    if (length(unmatched) <= 10) {
      cat("    Unmatched genes:", paste(head(unmatched, 10), collapse=", "), "\n")
    }
  }
  
  # Add to rowData
  rowData(thyr_se_strand2_nonzero)$gene_length <- matched_lengths
  gene_length_added <- TRUE
  
  # Re-save the SE with gene lengths
  saveRDS(thyr_se_strand2_nonzero, output_file)
  cat("  SE updated with gene lengths and re-saved\n")
  
  # Clean up
  rm(gene_lengths, gene_ids_clean, matched_lengths)
}

# Save summary statistics
summary_stats <- list(
  date = Sys.Date(),
  input_genes = nrow(thyr_se),
  input_samples = ncol(thyr_se),
  filtered_genes = nrow(thyr_se_strand2_nonzero),
  filtered_samples = ncol(thyr_se_strand2_nonzero),
  gene_length_added = gene_length_added,
  n_genes_with_length = n_genes_with_length,
  filter_criteria = list(
    zero_count_removal = TRUE,
    low_count_filter = FALSE,
    protein_coding_only = FALSE
  )
)
summary_file <- paste0(paths$logs, "counts_preparation_summary_",
                       format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
saveRDS(summary_stats, summary_file)

# ============================================================================
# Final summary
# ============================================================================
cat("\n=== Summary ===\n")
cat("Input:", nrow(thyr_se), "genes x", ncol(thyr_se), "samples\n")
cat("Output:", nrow(thyr_se_strand2_nonzero), "genes x", ncol(thyr_se_strand2_nonzero), "samples\n")
cat("Retention:", sprintf("%.1f%% genes", nrow(thyr_se_strand2_nonzero)/nrow(thyr_se)*100), "\n")

if (gene_length_added) {
  cat("\nGene lengths:\n")
  cat("  Source: GENCODE v36\n")
  cat("  Genes with length:", n_genes_with_length, 
      sprintf("(%.1f%%)", n_genes_with_length/nrow(thyr_se_strand2_nonzero)*100), "\n")
}

cat("\nFiles created:\n")
cat("  - thyr_se_strand2_nonzero.rds (main output)\n")
cat("  - gene_filter_info_*.rds (in logs/)\n")
cat("  - counts_preparation_summary_*.rds (in logs/)\n")

# Clean up environment
rm(list = setdiff(ls(), c("paths", "thyr_se_strand2_nonzero")))
gc()

cat("\n=== Preparation Complete ===\n")
cat("Next: Run 04_clinical_integration.R\n")