# ==============================================================================
# REBC-THYR Data Loading Script
# 01_data_loading.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)

# Load thyroid data
cat("Loading REBC-THYR data...\n")
load("./data/raw/thyr_data.rda")

# Rename data object for clarity
se_thyr <- data
rm(data)

cat("Data loaded successfully!\n")
cat("Data class:", class(se_thyr), "\n")
cat("Dimensions:", dim(se_thyr), "\n")
cat("Available assays:", paste(assayNames(se_thyr), collapse = ", "), "\n")

# Verify strandedness (from previous analysis)
cat("Verifying strandedness...\n")
first_counts <- colSums(assay(se_thyr, "stranded_first"))
second_counts <- colSums(assay(se_thyr, "stranded_second"))
strand_ratios <- second_counts / first_counts

cat("Strand ratio summary (second/first):\n")
print(summary(strand_ratios))

# Check for potential outliers
outliers <- which(strand_ratios < 10)
if(length(outliers) > 0) {
  cat("Samples with low strand ratios:", names(outliers), "\n")
} else {
  cat("All samples show consistent reverse strand protocol\n")
}

cat("Data loading completed!\n")

