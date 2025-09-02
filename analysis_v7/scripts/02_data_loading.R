# 02_data_loading.R - REBC-THYR Data Loading (v7)
# Purpose: Load 906 STAR-Counts files and create SummarizedExperiment

source("analysis_v7/setup.R")

cat("\n=== REBC-THYR Data Loading ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# Get TSV file paths
cat("\nScanning for TSV files...\n")
tsv_files <- list.files(
  paths$gdc,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
cat("Found", length(tsv_files), "TSV files\n")

if(length(tsv_files) == 0) {
  stop("No TSV files found. Please run gdc-client download first.")
}

# Read first file to get gene structure
cat("\nReading first file to get gene structure...\n")
first_data <- data.table::fread(tsv_files[1])
cat("Columns:", paste(colnames(first_data), collapse = ", "), "\n")
cat("Number of rows:", nrow(first_data), "\n")

# Extract gene info (skip first 4 summary rows)
gene_info <- first_data[5:nrow(first_data), 1:3]  # gene_id, gene_name, gene_type
summary_stats <- first_data[1:4, ]
n_genes <- nrow(gene_info)
n_files <- length(tsv_files)

cat("Genes:", n_genes, "\n")
cat("Files:", n_files, "\n")

# Get file IDs first
file_ids <- sapply(tsv_files, function(f) basename(dirname(f)))

# Initialize count matrices with proper dimensions and names
cat("\nInitializing count matrices...\n")
counts_unstranded <- matrix(0, nrow = n_genes, ncol = n_files)
counts_stranded_first <- matrix(0, nrow = n_genes, ncol = n_files)
counts_stranded_second <- matrix(0, nrow = n_genes, ncol = n_files)

# Set row and column names upfront
rownames(counts_unstranded) <- gene_info$gene_id
rownames(counts_stranded_first) <- gene_info$gene_id
rownames(counts_stranded_second) <- gene_info$gene_id

colnames(counts_unstranded) <- file_ids
colnames(counts_stranded_first) <- file_ids
colnames(counts_stranded_second) <- file_ids

# Function to read one file
read_star_counts <- function(file_path) {
  tryCatch({
    data <- data.table::fread(file_path)
    # Skip first 4 summary rows and select count columns
    counts <- data[5:nrow(data), ]
    return(list(
      unstranded = as.numeric(counts$unstranded),
      stranded_first = as.numeric(counts$stranded_first),
      stranded_second = as.numeric(counts$stranded_second),
      success = TRUE
    ))
  }, error = function(e) {
    warning("Error reading: ", basename(dirname(file_path)), " - ", e$message)
    return(list(success = FALSE))
  })
}

# Read all files with progress
cat("\nReading count files...\n")
pb <- txtProgressBar(min = 0, max = n_files, style = 3)

for(i in seq_len(n_files)) {
  result <- read_star_counts(tsv_files[i])
  
  if(result$success) {
    counts_unstranded[, i] <- result$unstranded
    counts_stranded_first[, i] <- result$stranded_first
    counts_stranded_second[, i] <- result$stranded_second
  }
  setTxtProgressBar(pb, i)
}
close(pb)
cat("\n")

# Create column data
cat("\nPreparing sample information...\n")
colData_df <- data.frame(
  file_id = file_ids,
  stringsAsFactors = FALSE
)

# Load manifest for metadata
manifest_files <- list.files(
  paths$raw,
  pattern = "^manifest_gene_counts_.*\\.tsv$",
  full.names = TRUE
)
if(length(manifest_files) > 0) {
  manifest <- data.table::fread(manifest_files[length(manifest_files)])
  # Use correct column names: file_name instead of filename
  colData_df <- merge(
    colData_df,
    manifest[, .(id, file_name, md5sum)],
    by.x = "file_id", by.y = "id",
    all.x = TRUE
  )
  cat("Manifest info added\n")
}

# Create row data (use gene_info which has correct columns)
rowData_df <- data.frame(
  gene_id = gene_info$gene_id,
  gene_name = gene_info$gene_name,
  gene_type = gene_info$gene_type,
  stringsAsFactors = FALSE
)

# Create SummarizedExperiment
cat("\nCreating SummarizedExperiment...\n")
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    unstranded = counts_unstranded,
    stranded_first = counts_stranded_first,
    stranded_second = counts_stranded_second
  ),
  colData = colData_df,
  rowData = rowData_df
)

cat("SummarizedExperiment created:\n")
cat("  Genes:", nrow(se), "\n")
cat("  Samples:", ncol(se), "\n")
cat("  Assays:", paste(names(assays(se)), collapse = ", "), "\n")

# Basic QC
cat("\nBasic QC:\n")
lib_sizes <- colSums(assay(se, "unstranded"))
cat("  Library sizes (unstranded):\n")
cat("    Min:", format(min(lib_sizes), big.mark = ","), "\n")
cat("    Median:", format(median(lib_sizes), big.mark = ","), "\n")
cat("    Max:", format(max(lib_sizes), big.mark = ","), "\n")

gene_detection <- rowSums(assay(se, "unstranded") > 0)
cat("  Genes detected (>0 counts):\n")
cat("    In any sample:", sum(gene_detection > 0), "\n")
cat("    In >50% samples:", sum(gene_detection > ncol(se)/2), "\n")

# Save results
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_file <- paste0(paths$processed, "thyr_se_raw_", timestamp, ".rds")
saveRDS(se, output_file)
cat("\nSaved to:", output_file, "\n")

# Also save without timestamp for easier access
saveRDS(se, paste0(paths$processed, "thyr_se_raw.rds"))
save(se, file = paste0(paths$processed, "thyr_se_raw.RData"))

# Save metadata
metadata_file <- paste0(paths$processed, "loading_metadata_", timestamp, ".rds")
saveRDS(list(
  loading_date = Sys.Date(),
  n_files = n_files,
  n_genes = n_genes,
  file_ids = file_ids,
  lib_sizes = lib_sizes
), metadata_file)

cat("\n=== Data Loading Completed ===\n")
cat("Files created:\n")
cat("  - thyr_se_raw.rds\n")
cat("  - thyr_se_raw.RData\n")
cat("  - loading_metadata_*.rds\n")
