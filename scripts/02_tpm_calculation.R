# ==============================================================================
# REBC-THYR TPM Calculation Script
# 02_tpm_calculation.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(GenomicFeatures)
library(GenomicRanges)
library(AnnotationDbi)
library(txdbmaker)
library(GenomicDataCommons)
library(dplyr)

# Load data (assuming se_thyr is already loaded from 01_data_loading.R)
# If not, uncomment the following lines:
# load("./data/raw/thyr_data.rda")
# se_thyr <- data
# rm(data)

# ==============================================================================
# 1. TxDb Creation
# ==============================================================================

if (!exists("txdb")) {
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
    cat("TxDb saved to ./data/raw/gencode.v36.annotation.txdb.sqlite\n")
  }
}

# ==============================================================================
# 2. Exon Length Calculation
# ==============================================================================

cat("Calculating exon lengths...\n")
exons_list <- GenomicFeatures::exonsBy(txdb, by = "gene")
exon_lengths <- sum(width(GenomicRanges::reduce(exons_list)))
names(exon_lengths) <- names(exons_list)

cat("Number of genes with exon length:", length(exon_lengths), "\n")
cat("Exon length summary:\n")
print(summary(exon_lengths))

# ==============================================================================
# 3. Read Length Acquisition
# ==============================================================================

cat("Fetching read length information from GDC...\n")
read_length_data <- GenomicDataCommons::files() %>%
  GenomicDataCommons::select(c("average_read_length", "cases.samples.submitter_id")) %>%
  GenomicDataCommons::filter(
    cases.project.project_id == "REBC-THYR" &
      data_type == "Aligned Reads" &
      experimental_strategy == "RNA-Seq" &
      data_format == "BAM"
  ) %>%
  GenomicDataCommons::response_all()

# Extract sample read lengths
valid_samples <- !is.na(read_length_data$results$average_read_length)
sample_read_lengths <- read_length_data$results$average_read_length[valid_samples]
sample_ids <- unlist(read_length_data$results$cases[valid_samples])
names(sample_read_lengths) <- sample_ids

cat("Number of samples with read length info:", length(sample_read_lengths), "\n")
cat("Read length summary:\n")
print(summary(sample_read_lengths))

# ==============================================================================
# 4. Data Preparation and Validation
# ==============================================================================

cat("Preparing count data...\n")
# Extract stranded_second count data
count_data <- assay(se_thyr, "stranded_second")
cat("Original count data dimensions:", dim(count_data), "\n")

# Gene ID matching and filtering
gene_ids_count <- rownames(count_data)
gene_ids_exon <- names(exon_lengths)
common_genes <- intersect(gene_ids_count, gene_ids_exon)

if (length(common_genes) == 0) {
  stop("No common genes found between count data and exon lengths")
}

count_data <- count_data[common_genes, ]
exon_lengths <- exon_lengths[common_genes]

cat("After gene filtering - Count data dimensions:", dim(count_data), "\n")
cat("Number of genes with exon lengths:", length(exon_lengths), "\n")

# Sample ID matching and filtering
sample_ids_count <- colnames(count_data)
sample_ids_readlen <- names(sample_read_lengths)
common_samples <- intersect(sample_ids_count, sample_ids_readlen)

if (length(common_samples) == 0) {
  stop("No common samples found between count data and read length data")
}

count_data <- count_data[, common_samples]
sample_read_lengths <- sample_read_lengths[common_samples]

cat("After sample filtering - Count data dimensions:", dim(count_data), "\n")
cat("Number of samples with read lengths:", length(sample_read_lengths), "\n")

# ==============================================================================
# 5. Effective Length Calculation
# ==============================================================================

cat("Calculating effective lengths...\n")
# Create effective length matrix for each gene-sample combination
effective_lengths <- outer(exon_lengths, sample_read_lengths, 
                           function(exon_len, read_len) {
                             pmax(exon_len - read_len + 1, 1)
                           })

# Verify dimensions
if (!all(dim(effective_lengths) == dim(count_data))) {
  stop("Dimension mismatch between effective lengths and count data")
}

cat("Effective length matrix dimensions:", dim(effective_lengths), "\n")
cat("Effective length summary:\n")
print(summary(as.vector(effective_lengths)))

# ==============================================================================
# 6. TPM Calculation
# ==============================================================================

cat("Calculating TPM values...\n")
# Calculate RPK (Reads Per Kilobase)
rpk <- count_data / effective_lengths

# Calculate scaling factors (per million mapped reads)
scaling_factors <- colSums(rpk)
cat("Scaling factors summary:\n")
print(summary(scaling_factors))

# Calculate TPM
tpm <- t(t(rpk) / scaling_factors) * 1e6

# Verify TPM properties
cat("TPM calculation complete!\n")
cat("TPM dimensions:", dim(tpm), "\n")
cat("TPM summary:\n")
print(summary(as.vector(tpm[1:1000, 1:10])))  # Sample subset for summary

# Check if TPM sums to ~1M per sample (should be exactly 1M)
tpm_sums <- colSums(tpm)
cat("TPM column sums (should be ~1,000,000):\n")
print(summary(tpm_sums))

# ==============================================================================
# 7. Quality Control
# ==============================================================================

cat("Performing quality control checks...\n")

# Check for negative values
negative_count <- sum(tpm < 0, na.rm = TRUE)
cat("Number of negative TPM values:", negative_count, "\n")

# Check for infinite values
infinite_count <- sum(is.infinite(tpm), na.rm = TRUE)
cat("Number of infinite TPM values:", infinite_count, "\n")

# Check for NA values
na_count <- sum(is.na(tpm))
cat("Number of NA TPM values:", na_count, "\n")

# Overall quality assessment
if (negative_count == 0 && infinite_count == 0 && na_count < nrow(tpm) * ncol(tpm) * 0.01) {
  cat("Quality assessment: EXCELLENT\n")
} else if (negative_count == 0 && infinite_count == 0) {
  cat("Quality assessment: GOOD (some NA values present)\n")
} else {
  cat("Quality assessment: NEEDS REVIEW (negative or infinite values present)\n")
}

# ==============================================================================
# 8. Save Results
# ==============================================================================

cat("Saving TPM data...\n")
# Create processed directory if it doesn't exist
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save TPM data
save(tpm, file = "./data/processed/thyr_tpm.rda")
cat("TPM data saved to ./data/processed/thyr_tpm.rda\n")

# Optionally save additional information
tpm_info <- list(
  dimensions = dim(tpm),
  gene_ids = rownames(tpm),
  sample_ids = colnames(tpm),
  calculation_date = Sys.time(),
  read_lengths = sample_read_lengths,
  exon_lengths = exon_lengths[rownames(tpm)]
)
save(tpm_info, file = "./data/processed/thyr_tpm_info.rda")
cat("TPM metadata saved to ./data/processed/thyr_tpm_info.rda\n")

# Clean up large objects to free memory
rm(txdb, exons_list, effective_lengths, rpk, read_length_data)
gc()

cat("TPM calculation completed successfully!\n")
cat("==============================================\n")

