# 02_data_loading.R - REBC-THYR Data Loading with Sample Mapping (v7)
# Purpose: Create file-sample mapping and load STAR-Counts into SummarizedExperiment

source("analysis_v7/setup.R")

cat("\n=== REBC-THYR Data Loading with Sample Mapping ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# ============================================================================
# Part 1: Create file_id to sample_submitter_id mapping
# ============================================================================

cat("\n--- Part 1: Creating sample mapping ---\n")

# Load manifest
manifest_files <- list.files(
  paths$raw,
  pattern = "^manifest_gene_counts_.*\\.tsv$",
  full.names = TRUE
)

if(length(manifest_files) == 0) {
  stop("No manifest files found. Please run 01_data_download.R first.")
}

# Use most recent manifest
manifest_file <- tail(sort(manifest_files), 1)
cat("Using manifest:", basename(manifest_file), "\n")
manifest <- data.table::fread(manifest_file)
cat("Total files in manifest:", nrow(manifest), "\n")

# Get all file IDs
all_file_ids <- manifest$id
n_files <- length(all_file_ids)

# Query GDC API with expand for cases and samples
cat("\nQuerying GDC API for sample information...\n")
cat("This may take a few moments for", n_files, "files\n")

# Process in batches of 50
batch_size <- 50
n_batches <- ceiling(n_files / batch_size)
all_results <- list()

pb <- txtProgressBar(min = 0, max = n_batches, style = 3)
for(i in seq_len(n_batches)) {
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, n_files)
  batch_ids <- all_file_ids[start_idx:end_idx]
  
  # Query batch
  batch_result <- tryCatch({
    GenomicDataCommons::files() %>%
      GenomicDataCommons::filter(~ file_id %in% batch_ids) %>%
      GenomicDataCommons::expand(c("cases", "cases.samples")) %>%
      GenomicDataCommons::results_all()
  }, error = function(e) {
    warning("Error in batch ", i, ": ", e$message)
    return(NULL)
  })
  
  if(!is.null(batch_result)) {
    all_results[[i]] <- batch_result
  }
  
  setTxtProgressBar(pb, i)
  Sys.sleep(0.1)  # Small delay to be nice to API
}
close(pb)

# Combine results
cat("\nProcessing API results...\n")

# 完全版の実装（Metastaticサンプルにも対応）
mapping_list <- list()

for(batch_result in all_results) {
  if(is.null(batch_result)) next
  
  file_id_to_name <- setNames(batch_result$file_name, batch_result$id)
  
  for(j in seq_along(batch_result$cases)) {
    file_id <- names(batch_result$cases)[j]
    case_data <- batch_result$cases[[j]]
    
    # Extract case information
    case_submitter_id <- case_data$submitter_id[1]
    
    # Extract sample information
    if(!is.null(case_data$samples) && length(case_data$samples) > 0) {
      sample_data <- case_data$samples[[1]]
      
      # sample_typeの処理（存在しない場合は代替情報を使用）
      if(length(sample_data$sample_type) > 0) {
        sample_type_value <- as.character(sample_data$sample_type)[1]
      } else if(length(sample_data$tumor_descriptor) > 0 && length(sample_data$specimen_type) > 0) {
        # Metastaticサンプルの場合、tumor_descriptorとspecimen_typeを結合
        sample_type_value <- paste(sample_data$tumor_descriptor[1], 
                                   sample_data$specimen_type[1], 
                                   "Tumor", sep = " ")
      } else {
        sample_type_value <- "Unknown"
      }
      
      mapping_list[[file_id]] <- data.frame(
        file_id = file_id,
        file_name = as.character(file_id_to_name[file_id])[1],
        case_submitter_id = as.character(case_submitter_id)[1],
        sample_submitter_id = as.character(sample_data$submitter_id)[1],
        sample_type = sample_type_value,
        sample_id = as.character(sample_data$sample_id)[1],
        stringsAsFactors = FALSE
      )
    }
  }
}

# Combine into single data frame
file_sample_map <- data.table::rbindlist(mapping_list, fill = TRUE)
cat("Mapping created for", nrow(file_sample_map), "files\n")

# サンプルタイプの分布を確認
cat("\nSample type distribution:\n")
print(table(file_sample_map$sample_type))

# Metastaticサンプルの詳細確認
metastatic_samples <- file_sample_map[grep("Metastatic", file_sample_map$sample_type), ]
cat("\nMetastatic samples found:", nrow(metastatic_samples), "\n")
if(nrow(metastatic_samples) > 0) {
  cat("Example Metastatic sample IDs:\n")
  print(head(metastatic_samples$sample_submitter_id))
}

# 確認
head(file_sample_map[, c("file_id", "sample_submitter_id", "case_submitter_id")])
# Verify all manifest files are mapped
unmapped <- setdiff(manifest$id, file_sample_map$file_id)
if(length(unmapped) > 0) {
  warning("Unmapped file IDs: ", length(unmapped))
  cat("First few unmapped IDs:\n")
  print(head(unmapped))
}

# Save mapping
mapping_file <- paste0(paths$processed, "file_sample_mapping_", 
                       format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv")
data.table::fwrite(file_sample_map, mapping_file, sep = "\t")
cat("Mapping saved to:", basename(mapping_file), "\n")

# Summary
cat("\nMapping summary:\n")
cat("  Unique cases:", length(unique(file_sample_map$case_submitter_id)), "\n")
cat("  Unique samples:", length(unique(file_sample_map$sample_submitter_id)), "\n")
sample_type_table <- table(file_sample_map$sample_type)
cat("  Sample types:\n")
for(st in names(sample_type_table)) {
  cat("    ", st, ":", sample_type_table[st], "\n")
}

# ============================================================================
# Part 2: Load count data with correct sample IDs
# ============================================================================

cat("\n--- Part 2: Loading count data ---\n")

# Get TSV file paths
tsv_files <- list.files(
  paths$gdc,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
cat("Found", length(tsv_files), "TSV files\n")

if(length(tsv_files) == 0) {
  stop("No TSV files found. Please ensure gdc-client download completed.")
}

# Extract file_id from path (directory name)
extract_file_id <- function(path) {
  parts <- strsplit(path, "/")[[1]]
  # The file_id is the directory name above the TSV file
  parts[length(parts) - 1]
}

file_ids_from_path <- sapply(tsv_files, extract_file_id)

# Match with mapping
tsv_file_map <- data.frame(
  path = tsv_files,
  file_id = file_ids_from_path,
  stringsAsFactors = FALSE
)

# Merge with sample mapping
tsv_file_map <- merge(
  tsv_file_map,
  file_sample_map[, c("file_id", "sample_submitter_id", "case_submitter_id", "sample_type")],
  by = "file_id",
  all.x = TRUE
)

# Check for missing mappings
missing_mapping <- is.na(tsv_file_map$sample_submitter_id)
if(any(missing_mapping)) {
  warning("Files without sample mapping: ", sum(missing_mapping))
  cat("Files without mapping:\n")
  print(tsv_file_map$file_id[missing_mapping])
}

# Remove files without mapping
tsv_file_map <- tsv_file_map[!missing_mapping, ]
cat("\nFiles with valid mapping:", nrow(tsv_file_map), "\n")

# Sort by sample_submitter_id for consistent ordering
tsv_file_map <- tsv_file_map[order(tsv_file_map$sample_submitter_id), ]

# Read first file to get gene structure
cat("\nReading first file to get gene structure...\n")
first_data <- data.table::fread(tsv_file_map$path[1])

# Extract gene info (skip first 4 summary rows)
gene_info <- first_data[5:nrow(first_data), 1:3]
colnames(gene_info) <- c("gene_id", "gene_name", "gene_type")
n_genes <- nrow(gene_info)
n_samples <- nrow(tsv_file_map)

cat("Genes:", n_genes, "\n")
cat("Samples:", n_samples, "\n")

# Initialize count matrices
cat("\nInitializing count matrices...\n")
counts_unstranded <- matrix(0, nrow = n_genes, ncol = n_samples)
counts_stranded_first <- matrix(0, nrow = n_genes, ncol = n_samples)
counts_stranded_second <- matrix(0, nrow = n_genes, ncol = n_samples)

# Set row and column names
rownames(counts_unstranded) <- gene_info$gene_id
rownames(counts_stranded_first) <- gene_info$gene_id
rownames(counts_stranded_second) <- gene_info$gene_id

# Use sample_submitter_id as column names
colnames(counts_unstranded) <- tsv_file_map$sample_submitter_id
colnames(counts_stranded_first) <- tsv_file_map$sample_submitter_id
colnames(counts_stranded_second) <- tsv_file_map$sample_submitter_id

# Function to read one file
read_star_counts <- function(file_path) {
  tryCatch({
    data <- data.table::fread(file_path)
    # Skip first 4 summary rows
    counts <- data[5:nrow(data), ]
    return(list(
      unstranded = as.numeric(counts$unstranded),
      stranded_first = as.numeric(counts$stranded_first),
      stranded_second = as.numeric(counts$stranded_second),
      success = TRUE
    ))
  }, error = function(e) {
    warning("Error reading: ", basename(file_path), " - ", e$message)
    return(list(success = FALSE))
  })
}

# Read all files with progress
cat("\nReading count files...\n")
pb <- txtProgressBar(min = 0, max = n_samples, style = 3)

for(i in seq_len(n_samples)) {
  result <- read_star_counts(tsv_file_map$path[i])
  
  if(result$success) {
    counts_unstranded[, i] <- result$unstranded
    counts_stranded_first[, i] <- result$stranded_first
    counts_stranded_second[, i] <- result$stranded_second
  }
  setTxtProgressBar(pb, i)
}
close(pb)
cat("\n")

# Create column data with all metadata
cat("\nPreparing sample metadata...\n")
colData_df <- data.frame(
  sample_submitter_id = tsv_file_map$sample_submitter_id,
  case_submitter_id = tsv_file_map$case_submitter_id,
  sample_type = tsv_file_map$sample_type,
  file_id = tsv_file_map$file_id,
  stringsAsFactors = FALSE
)
rownames(colData_df) <- colData_df$sample_submitter_id

# Create row data
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

# Sample type distribution
cat("\n  Sample types in SE:\n")
sample_type_se <- table(colData(se)$sample_type)
for(st in names(sample_type_se)) {
  cat("    ", st, ":", sample_type_se[st], "\n")
}

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
  n_files_manifest = nrow(manifest),
  n_files_loaded = n_samples,
  n_genes = n_genes,
  file_sample_map = file_sample_map,
  lib_sizes = lib_sizes
), metadata_file)

cat("\n=== Data Loading Completed ===\n")
cat("Files created:\n")
cat("  - file_sample_mapping_*.tsv (mapping table)\n")
cat("  - thyr_se_raw.rds (SummarizedExperiment)\n")
cat("  - thyr_se_raw.RData\n")
cat("  - loading_metadata_*.rds\n")
cat("\nIMPORTANT: Sample IDs are now correctly mapped to sample_submitter_id\n")
cat("Example sample IDs:", paste(head(colnames(se), 3), collapse = ", "), "\n")