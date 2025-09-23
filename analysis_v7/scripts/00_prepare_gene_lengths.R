# 00_prepare_gene_lengths.R - Prepare Gene Lengths from GENCODE GTF
# Purpose: Calculate exon union lengths for all genes from GENCODE v36 GTF
# Output: gene_lengths.rds with ENSG ID to length mapping
# Version: v1.0 - For TPM calculation in downstream analyses
# Date: 2025-01-23

source("analysis_v7/setup.R")

cat("\n=== Gene Length Preparation from GENCODE v36 ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Purpose: Calculate exon union lengths for TPM normalization\n")

# Load required packages
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
})

# ============================================================================
# Configuration
# ============================================================================

CONFIG <- list(
  # GENCODE v36 URL (matches TCGA/GDC annotation)
  GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz",
  
  # Local file names
  GTF_GZ = "gencode.v36.annotation.gtf.gz",
  GTF_FILE = "gencode.v36.annotation.gtf",
  
  # Output
  OUTPUT_FILE = "gene_lengths.rds",
  
  # Processing
  KEEP_GTF = FALSE  # Keep GTF file after processing
)

# ============================================================================
# Step 1: Download GTF file
# ============================================================================

cat("\n--- Step 1: Downloading GENCODE v36 GTF ---\n")

gtf_gz_path <- paste0(paths$raw, CONFIG$GTF_GZ)

if (file.exists(gtf_gz_path)) {
  cat("  GTF file already exists, skipping download\n")
} else {
  cat("  Downloading from GENCODE...\n")
  cat("  URL:", CONFIG$GTF_URL, "\n")
  
  download_start <- Sys.time()
  
  tryCatch({
    download.file(
      url = CONFIG$GTF_URL,
      destfile = gtf_gz_path,
      method = "curl",
      quiet = FALSE
    )
    
    download_time <- difftime(Sys.time(), download_start, units = "secs")
    cat(sprintf("  Download completed in %.1f seconds\n", download_time))
    
    # Check file size
    file_size_mb <- file.info(gtf_gz_path)$size / 1024^2
    cat(sprintf("  File size: %.1f MB\n", file_size_mb))
    
  }, error = function(e) {
    stop("Failed to download GTF file: ", e$message)
  })
}

# ============================================================================
# Step 2: Decompress GTF file
# ============================================================================

cat("\n--- Step 2: Decompressing GTF file ---\n")

gtf_path <- paste0(paths$raw, CONFIG$GTF_FILE)

if (file.exists(gtf_path)) {
  cat("  Decompressed GTF already exists\n")
} else {
  cat("  Decompressing...\n")
  
  decompress_start <- Sys.time()
  
  # Use R.utils for cross-platform compatibility
  if (requireNamespace("R.utils", quietly = TRUE)) {
    R.utils::gunzip(gtf_gz_path, gtf_path, remove = FALSE, overwrite = TRUE)
  } else {
    # Fallback to system command
    system(sprintf("gunzip -c %s > %s", gtf_gz_path, gtf_path))
  }
  
  decompress_time <- difftime(Sys.time(), decompress_start, units = "secs")
  cat(sprintf("  Decompression completed in %.1f seconds\n", decompress_time))
  
  # Check decompressed file size
  file_size_mb <- file.info(gtf_path)$size / 1024^2
  cat(sprintf("  Decompressed size: %.1f MB\n", file_size_mb))
}

# ============================================================================
# Step 3: Load GTF and filter for exons
# ============================================================================

cat("\n--- Step 3: Loading GTF file ---\n")
cat("  This may take a few minutes...\n")

load_start <- Sys.time()

# Import GTF using rtracklayer
gtf <- import(gtf_path)

load_time <- difftime(Sys.time(), load_start, units = "mins")
cat(sprintf("  GTF loaded in %.1f minutes\n", load_time))
cat(sprintf("  Total features: %s\n", format(length(gtf), big.mark = ",")))

# Filter for exons only
exons <- gtf[gtf$type == "exon"]
cat(sprintf("  Exon features: %s\n", format(length(exons), big.mark = ",")))

# Extract gene IDs (remove version numbers)
exons$gene_id_clean <- gsub("\\..*", "", exons$gene_id)

# Check unique genes
unique_genes <- unique(exons$gene_id_clean)
cat(sprintf("  Unique genes: %s\n", format(length(unique_genes), big.mark = ",")))

# ============================================================================
# Step 4: Calculate exon union lengths
# ============================================================================

cat("\n--- Step 4: Calculating exon union lengths ---\n")
cat("  Computing merged exon lengths per gene...\n")

calc_start <- Sys.time()

# Split exons by gene
exons_by_gene <- split(exons, exons$gene_id_clean)

# Calculate union length for each gene
calculate_union_length <- function(gene_exons) {
  # Reduce overlapping exons to get union
  gene_ranges <- reduce(gene_exons)
  # Sum the widths
  sum(width(gene_ranges))
}

# Apply to all genes (with progress)
n_genes <- length(exons_by_gene)
cat(sprintf("  Processing %s genes...\n", format(n_genes, big.mark = ",")))

# Process in chunks for memory efficiency
chunk_size <- 1000
gene_lengths <- numeric(n_genes)
names(gene_lengths) <- names(exons_by_gene)

pb <- txtProgressBar(min = 0, max = n_genes, style = 3)

for (i in seq_len(n_genes)) {
  gene_lengths[i] <- calculate_union_length(exons_by_gene[[i]])
  
  if (i %% chunk_size == 0 || i == n_genes) {
    setTxtProgressBar(pb, i)
  }
}

close(pb)

calc_time <- difftime(Sys.time(), calc_start, units = "secs")
cat(sprintf("\n  Calculation completed in %.1f seconds\n", calc_time))

# ============================================================================
# Step 5: Quality control
# ============================================================================

cat("\n--- Step 5: Quality control ---\n")

# Check for zero-length genes
zero_length <- sum(gene_lengths == 0)
if (zero_length > 0) {
  cat(sprintf("  Warning: %d genes have zero length\n", zero_length))
  gene_lengths <- gene_lengths[gene_lengths > 0]
}

# Summary statistics
cat("\nGene length statistics:\n")
cat(sprintf("  Total genes: %s\n", format(length(gene_lengths), big.mark = ",")))
cat(sprintf("  Min length: %d bp\n", min(gene_lengths)))
cat(sprintf("  Max length: %s bp\n", format(max(gene_lengths), big.mark = ",")))
cat(sprintf("  Median length: %s bp\n", format(median(gene_lengths), big.mark = ",")))
cat(sprintf("  Mean length: %s bp\n", format(round(mean(gene_lengths)), big.mark = ",")))

# Distribution check
length_quantiles <- quantile(gene_lengths, probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
cat("\nLength distribution (percentiles):\n")
for (i in seq_along(length_quantiles)) {
  cat(sprintf("  %3s%%: %s bp\n", 
              names(length_quantiles)[i], 
              format(round(length_quantiles[i]), big.mark = ",")))
}

# Check for common genes
common_genes <- c("GAPDH", "ACTB", "TP53", "EGFR", "KRAS")
cat("\nSpot check of common genes:\n")
for (gene in common_genes) {
  # Find gene (may have version number in original)
  gene_pattern <- paste0("^", gene, "($|\\.)")
  matching <- grep(gene_pattern, names(gene_lengths))
  
  if (length(matching) > 0) {
    gene_id <- names(gene_lengths)[matching[1]]
    cat(sprintf("  %s (%s): %s bp\n", 
                gene, gene_id,
                format(gene_lengths[matching[1]], big.mark = ",")))
  } else {
    # Try searching in the original gene names
    matching_orig <- grep(gene, exons$gene_name)
    if (length(matching_orig) > 0) {
      gene_id_orig <- unique(exons$gene_id_clean[matching_orig])[1]
      if (gene_id_orig %in% names(gene_lengths)) {
        cat(sprintf("  %s (%s): %s bp\n", 
                    gene, gene_id_orig,
                    format(gene_lengths[gene_id_orig], big.mark = ",")))
      }
    }
  }
}

# ============================================================================
# Step 6: Save results
# ============================================================================

cat("\n--- Step 6: Saving results ---\n")

output_path <- paste0(paths$processed, CONFIG$OUTPUT_FILE)

# Save as RDS
saveRDS(gene_lengths, output_path)
cat(sprintf("  Gene lengths saved to: %s\n", basename(output_path)))

# Also save as CSV for reference
csv_path <- gsub("\\.rds$", ".csv", output_path)
gene_length_df <- data.frame(
  gene_id = names(gene_lengths),
  length = gene_lengths,
  stringsAsFactors = FALSE
)
write.csv(gene_length_df, csv_path, row.names = FALSE)
cat(sprintf("  CSV backup saved to: %s\n", basename(csv_path)))

# ============================================================================
# Step 7: Cleanup
# ============================================================================

cat("\n--- Step 7: Cleanup ---\n")

if (!CONFIG$KEEP_GTF) {
  cat("  Removing temporary GTF files...\n")
  
  if (file.exists(gtf_path)) {
    file.remove(gtf_path)
    cat("    Removed:", CONFIG$GTF_FILE, "\n")
  }
  
  # Optionally keep the compressed version
  # if (file.exists(gtf_gz_path)) {
  #   file.remove(gtf_gz_path)
  #   cat("    Removed:", CONFIG$GTF_GZ, "\n")
  # }
  
  cat("    Kept compressed GTF for potential reuse\n")
} else {
  cat("  GTF files retained for future use\n")
}

# ============================================================================
# Final summary
# ============================================================================

cat("\n=== Gene Length Preparation Complete ===\n")
cat("Summary:\n")
cat(sprintf("  GENCODE version: v36\n"))
cat(sprintf("  Total genes processed: %s\n", format(length(gene_lengths), big.mark = ",")))
cat(sprintf("  Output file: %s\n", CONFIG$OUTPUT_FILE))
cat(sprintf("  File location: %s\n", output_path))

cat("\nUsage in downstream scripts:\n")
cat("  gene_lengths <- readRDS('", output_path, "')\n", sep = "")
cat("  # For TPM calculation:\n")
cat("  # tpm <- calculate_tpm(counts, gene_lengths[rownames(counts)])\n")

cat("\nNext steps:\n")
cat("  1. Update 03_prepare_counts.R to add gene lengths to SE\n")
cat("  2. Use in 11_reo_pair_selection.R for TPM calculation\n")

# Clean up workspace
rm(gtf, exons, exons_by_gene)
gc()

cat("\n=== Script Complete ===\n")