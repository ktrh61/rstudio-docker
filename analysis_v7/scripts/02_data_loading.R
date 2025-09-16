# 02_data_loading.R - REBC-THYR Data Loading with Sample Mapping (v7.2)
# Purpose: Create file-sample mapping and load STAR-Counts into SummarizedExperiment
# Output: SE with all assays (counts + expression values) and sample metadata
# v7.2: Fixed column structure (9 cols), added all expression assays to SE, memory management

source("analysis_v7/setup.R")

cat("\n=== REBC-THYR Data Loading with Sample Mapping (v7.2) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# Check if running on Unix-like system for parallel processing
if(.Platform$OS.type != "unix") {
  stop("This version requires Unix/Linux environment (WSL2, Mac, or Linux) for parallel processing")
}

library(parallel)
library(data.table)

# ============================================================================
# Part 1: Create file_id to sample_submitter_id mapping (with parallel result processing)
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
# NOTE: API queries remain sequential to respect rate limits
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

# Clean up progress bar variable
rm(pb)

# ============================================================================
# Process API results in parallel
# ============================================================================
cat("\nProcessing API results in parallel...\n")

# Function to process one batch result
process_batch_result <- function(batch_result) {
  if(is.null(batch_result)) return(NULL)
  
  # Set data.table threads to 1 in worker
  data.table::setDTthreads(1L)
  
  batch_mapping <- list()
  file_id_to_name <- setNames(batch_result$file_name, batch_result$id)
  
  for(j in seq_along(batch_result$cases)) {
    file_id <- names(batch_result$cases)[j]
    case_data <- batch_result$cases[[j]]
    
    # Extract case information
    case_submitter_id <- case_data$submitter_id[1]
    
    # Extract sample information
    if(!is.null(case_data$samples) && length(case_data$samples) > 0) {
      sample_data <- case_data$samples[[1]]
      
      # Handle sample_type (may not exist for some samples)
      if(length(sample_data$sample_type) > 0) {
        sample_type_value <- as.character(sample_data$sample_type)[1]
      } else if(length(sample_data$tumor_descriptor) > 0 && length(sample_data$specimen_type) > 0) {
        # Metastatic samples: combine tumor_descriptor and specimen_type
        sample_type_value <- paste(sample_data$tumor_descriptor[1], 
                                   sample_data$specimen_type[1], 
                                   "Tumor", sep = " ")
      } else {
        sample_type_value <- "Unknown"
      }
      
      batch_mapping[[file_id]] <- data.frame(
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
  
  return(batch_mapping)
}

# Determine number of cores for result processing
n_cores <- detectCores()
cores_to_use_part1 <- max(1, min(n_cores - 1, 4))  # Cap at 4 cores for this step
cat("Using", cores_to_use_part1, "cores for result processing\n")

# Process all batch results in parallel
mapping_lists <- mclapply(
  all_results,
  process_batch_result,
  mc.cores = cores_to_use_part1,
  mc.preschedule = FALSE  # Dynamic scheduling for variable processing time
)

# Flatten the nested list structure
mapping_list <- unlist(mapping_lists, recursive = FALSE)
mapping_list <- mapping_list[!sapply(mapping_list, is.null)]

# Combine into single data frame
file_sample_map <- data.table::rbindlist(mapping_list, fill = TRUE)
cat("Mapping created for", nrow(file_sample_map), "files\n")

# Clean up Part 1 intermediate objects
rm(all_results, mapping_lists, mapping_list, batch_result, 
   all_file_ids, batch_ids, batch_size, n_batches, i, 
   start_idx, end_idx, n_files)

# Sample type distribution
cat("\nSample type distribution:\n")
print(table(file_sample_map$sample_type))

# Check for Metastatic samples
metastatic_samples <- file_sample_map[grep("Metastatic", file_sample_map$sample_type), ]
cat("\nMetastatic samples found:", nrow(metastatic_samples), "\n")
if(nrow(metastatic_samples) > 0) {
  cat("Example Metastatic sample IDs:\n")
  print(head(metastatic_samples$sample_submitter_id))
}

# Verify all manifest files are mapped
unmapped <- setdiff(manifest$id, file_sample_map$file_id)
if(length(unmapped) > 0) {
  warning("Unmapped file IDs: ", length(unmapped))
  cat("First few unmapped IDs:\n")
  print(head(unmapped))
}

# Clean up verification variables
rm(unmapped, metastatic_samples)

# Save mapping
mapping_file <- paste0(paths$logs, "file_sample_mapping_", 
                       format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv")
data.table::fwrite(file_sample_map, mapping_file, sep = "\t")
cat("Mapping saved to:", basename(mapping_file), " (in logs/)\n")

# Summary
cat("\nMapping summary:\n")
cat("  Unique cases:", length(unique(file_sample_map$case_submitter_id)), "\n")
cat("  Unique samples:", length(unique(file_sample_map$sample_submitter_id)), "\n")
sample_type_table <- table(file_sample_map$sample_type)
cat("  Sample types:\n")
for(st in names(sample_type_table)) {
  cat("    ", st, ":", sample_type_table[st], "\n")
}

# Clean up summary variables
rm(sample_type_table, st)

# Clean up manifest (no longer needed)
rm(manifest, manifest_file, manifest_files)

# ============================================================================
# Part 2: Load count data with correct sample IDs (FULLY OPTIMIZED PARALLEL VERSION)
# ============================================================================

cat("\n--- Part 2: Loading count data (Fully Optimized Parallel) ---\n")

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
  parts[length(parts) - 1]
}

file_ids_from_path <- sapply(tsv_files, extract_file_id)

# Match with mapping
tsv_file_map <- data.frame(
  path = tsv_files,
  file_id = file_ids_from_path,
  stringsAsFactors = FALSE
)

# Clean up intermediate variables
rm(tsv_files, file_ids_from_path)

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

# Clean up missing_mapping check
rm(missing_mapping)

# IMPORTANT: Sort by sample_submitter_id for consistent ordering
tsv_file_map <- tsv_file_map[order(tsv_file_map$sample_submitter_id), ]

# Read first file to get gene structure (optimized)
cat("\nReading first file to get gene structure...\n")

# Read with proper handling of comment and summary rows
first_data <- data.table::fread(
  tsv_file_map$path[1],
  skip = 1L,  # Skip only comment line
  showProgress = FALSE,
  nThread = 1L
)

# Remove the first 4 rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
first_gene_info <- first_data[5:nrow(first_data), 1:3]
colnames(first_gene_info) <- c("gene_id", "gene_name", "gene_type")
n_genes <- nrow(first_gene_info)
n_samples <- nrow(tsv_file_map)

# Clean up first_data (no longer needed)
rm(first_data)

cat("Genes:", n_genes, "\n")
cat("Samples:", n_samples, "\n")

# ============================================================================
# OPTIMIZED: Prepare lightweight arguments (just vectors, not data.frames)
# ============================================================================
tsv_file_paths <- tsv_file_map$path  # FIXED: Changed from 'paths' to 'tsv_file_paths'
sample_ids <- tsv_file_map$sample_submitter_id
idx <- seq_len(nrow(tsv_file_map))

# ============================================================================
# OPTIMIZED: Reading function with correct column structure (9 columns total)
# ============================================================================
read_star_counts_optimized <- function(i, paths_vec, ids_vec) {
  # CRITICAL: Set data.table threads to 1 to avoid nested parallelism
  data.table::setDTthreads(1L)
  
  path <- paths_vec[i]
  sid <- ids_vec[i]
  
  tryCatch({
    # Read file skipping only the comment line
    DT <- data.table::fread(
      path,
      skip = 1L,  # Skip only the comment line, keep header
      showProgress = FALSE,
      nThread = 1L
    )
    
    # Remove the first 4 rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
    DT <- DT[5:nrow(DT), ]
    
    # Extract count and expression columns by position (9 columns total)
    list(
      idx = i,
      sample_id = sid,
      unstranded = as.numeric(DT[[4]]),       # Column 4: unstranded counts
      stranded_first = as.numeric(DT[[5]]),   # Column 5: stranded_first counts
      stranded_second = as.numeric(DT[[6]]),  # Column 6: stranded_second counts
      tpm_unstranded = as.numeric(DT[[7]]),   # Column 7: TPM unstranded
      fpkm_unstranded = as.numeric(DT[[8]]),  # Column 8: FPKM unstranded
      fpkm_uq_unstranded = as.numeric(DT[[9]]), # Column 9: FPKM-UQ unstranded
      success = TRUE
    )
  }, error = function(e) {
    warning("Error reading file ", i, " (", basename(path), "): ", e$message)
    list(
      idx = i,
      sample_id = sid,
      success = FALSE
    )
  })
}

# Determine number of cores to use
n_cores <- detectCores()
cores_to_use <- max(1, min(n_cores - 1, 8))
cat("\nUsing", cores_to_use, "cores for TSV file processing\n")

# Read all files in parallel with optimized settings
cat("\nReading count files in parallel (fully optimized)...\n")
start_time <- Sys.time()

# CRITICAL: Use lightweight arguments and dynamic scheduling
results <- mclapply(
  idx,
  read_star_counts_optimized,
  paths_vec = tsv_file_paths,      # FIXED: Using renamed variable
  ids_vec = sample_ids,   # Lightweight vector  
  mc.cores = cores_to_use,
  mc.preschedule = FALSE  # CRITICAL: Dynamic scheduling for variable I/O times
)

end_time <- Sys.time()
cat("Reading completed in", format(end_time - start_time), "\n")

# Clean up reading variables (no longer needed)
rm(tsv_file_paths, idx)

# Check for failed reads
failed_reads <- sapply(results, function(x) !isTRUE(x$success))
if(any(failed_reads)) {
  n_failed <- sum(failed_reads)
  warning("Failed to read ", n_failed, " files")
  if(n_failed <= 10) {
    failed_samples <- sapply(results[failed_reads], function(x) x$sample_id)
    cat("Failed samples:\n")
    print(failed_samples)
    rm(failed_samples)
  } else {
    cat("Too many failures to list (", n_failed, " files)\n")
  }
  rm(n_failed)
}

# Clean up failure check variables
rm(failed_reads)

# ============================================================================
# Efficient matrix assembly
# ============================================================================
cat("\nAssembling count and expression matrices...\n")

# Pre-allocate count matrices
counts_unstranded <- matrix(0, nrow = n_genes, ncol = n_samples)
counts_stranded_first <- matrix(0, nrow = n_genes, ncol = n_samples)
counts_stranded_second <- matrix(0, nrow = n_genes, ncol = n_samples)

# Pre-allocate expression matrices (TPM and FPKM)
tpm_unstranded <- matrix(0, nrow = n_genes, ncol = n_samples)
fpkm_unstranded <- matrix(0, nrow = n_genes, ncol = n_samples)
fpkm_uq_unstranded <- matrix(0, nrow = n_genes, ncol = n_samples)

# Set row names once
rownames(counts_unstranded) <- first_gene_info$gene_id
rownames(counts_stranded_first) <- first_gene_info$gene_id
rownames(counts_stranded_second) <- first_gene_info$gene_id
rownames(tpm_unstranded) <- first_gene_info$gene_id
rownames(fpkm_unstranded) <- first_gene_info$gene_id
rownames(fpkm_uq_unstranded) <- first_gene_info$gene_id

# Set column names once
colnames(counts_unstranded) <- sample_ids
colnames(counts_stranded_first) <- sample_ids
colnames(counts_stranded_second) <- sample_ids
colnames(tpm_unstranded) <- sample_ids
colnames(fpkm_unstranded) <- sample_ids
colnames(fpkm_uq_unstranded) <- sample_ids

# Fill matrices using direct indexing (faster than lookup)
cat("Populating count and expression matrices...\n")
pb <- txtProgressBar(min = 0, max = length(results), style = 3)

for(i in seq_along(results)) {
  if(results[[i]]$success) {
    col_idx <- results[[i]]$idx  # Direct index, no lookup needed
    # Count matrices
    counts_unstranded[, col_idx] <- results[[i]]$unstranded
    counts_stranded_first[, col_idx] <- results[[i]]$stranded_first
    counts_stranded_second[, col_idx] <- results[[i]]$stranded_second
    # Expression matrices
    tpm_unstranded[, col_idx] <- results[[i]]$tpm_unstranded
    fpkm_unstranded[, col_idx] <- results[[i]]$fpkm_unstranded
    fpkm_uq_unstranded[, col_idx] <- results[[i]]$fpkm_uq_unstranded
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# Clean up progress bar and loop variables
rm(pb, i, col_idx)

# CRITICAL: Clean up the large results object (several GB)
rm(results)
gc()  # Force garbage collection to free memory

# Verify non-zero counts
cat("\nVerification:\n")
cat("  Non-zero counts in unstranded:", sum(counts_unstranded > 0), "\n")
cat("  Non-zero counts in stranded_first:", sum(counts_stranded_first > 0), "\n")
cat("  Non-zero counts in stranded_second:", sum(counts_stranded_second > 0), "\n")
cat("  Non-zero expression in TPM:", sum(tpm_unstranded > 0), "\n")
cat("  Non-zero expression in FPKM:", sum(fpkm_unstranded > 0), "\n")
cat("  Non-zero expression in FPKM-UQ:", sum(fpkm_uq_unstranded > 0), "\n")

# Memory usage check
mem_used <- as.numeric(
  object.size(counts_unstranded) + 
    object.size(counts_stranded_first) + 
    object.size(counts_stranded_second) +
    object.size(tpm_unstranded) +
    object.size(fpkm_unstranded) +
    object.size(fpkm_uq_unstranded)
) / 1024^3
cat("  Memory used by all matrices: ", round(mem_used, 2), " GB\n")

# Clean up memory calculation variable
rm(mem_used)

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

# Create row data with gene_name and gene_type preserved
rowData_df <- data.frame(
  gene_id = first_gene_info$gene_id,
  gene_name = first_gene_info$gene_name,
  gene_type = first_gene_info$gene_type,
  stringsAsFactors = FALSE
)

# Clean up gene info (keeping rowData_df for SE)
rm(first_gene_info)

# Create SummarizedExperiment with all assays
cat("\nCreating SummarizedExperiment with all assays...\n")
thyr_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    # Count assays
    unstranded = counts_unstranded,
    stranded_first = counts_stranded_first,
    stranded_second = counts_stranded_second,
    # Expression assays
    tpm_unstranded = tpm_unstranded,
    fpkm_unstranded = fpkm_unstranded,
    fpkm_uq_unstranded = fpkm_uq_unstranded
  ),
  colData = colData_df,
  rowData = rowData_df
)

# Clean up individual matrices (now stored in thyr_se)
rm(counts_unstranded, counts_stranded_first, counts_stranded_second)
rm(tpm_unstranded, fpkm_unstranded, fpkm_uq_unstranded)

# Verify alignment
cat("\nAlignment verification:\n")
cat("  Assay colnames match colData rownames:", 
    identical(colnames(assay(thyr_se, "stranded_second")), rownames(colData(thyr_se))), "\n")
cat("  All names are sample_submitter_id:", 
    all(colnames(assay(thyr_se, "stranded_second")) == colData(thyr_se)$sample_submitter_id), "\n")

cat("\nSummarizedExperiment created:\n")
cat("  Genes:", nrow(thyr_se), "\n")
cat("  Samples:", ncol(thyr_se), "\n")
cat("  Assays:", paste(names(assays(thyr_se)), collapse = ", "), "\n")
cat("  Gene metadata: gene_id, gene_name, gene_type\n")

# Save SummarizedExperiment (fixed name only)
saveRDS(thyr_se, paste0(paths$processed, "thyr_se_raw.rds"))
cat("\nSaved to:", paste0(paths$processed, "thyr_se_raw.rds"), "\n")

# Save metadata to logs
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
metadata_file <- paste0(paths$logs, "loading_metadata_", timestamp, ".rds")
saveRDS(list(
  loading_date = Sys.Date(),
  n_files_manifest = nrow(file_sample_map),
  n_files_loaded = n_samples,
  n_genes = n_genes,
  file_sample_map = file_sample_map,
  parallel_cores_part1 = cores_to_use_part1,
  parallel_cores_part2 = cores_to_use,
  processing_time = end_time - start_time,
  optimization_flags = list(
    dynamic_scheduling = TRUE,
    dt_threads = 1,
    selective_columns = TRUE,
    direct_numeric = TRUE,
    all_assays_preserved = TRUE
  )
), metadata_file)
cat("Metadata saved to logs/\n")

# Clean up timing and file variables
rm(start_time, end_time, timestamp, metadata_file)

# NOTE: Individual matrices already cleaned up after assembly (line 439-441)
# No need to remove them again here

# Clean up data frames and metadata
rm(tsv_file_map, colData_df, rowData_df, file_sample_map)

# Clean up processing variables
rm(n_genes, n_samples, n_cores, cores_to_use, cores_to_use_part1)
rm(sample_ids, mapping_file)

# Clean up the functions we defined
rm(extract_file_id, process_batch_result, read_star_counts_optimized)

# Force garbage collection
gc()

# ============================================================================
# Strand-specific count comparison and final cleanup
# ============================================================================
cat("\n=== Strand-specific count comparison ===\n")

# Compare total counts across strand-specific assays
total_counts <- c(
  unstranded = sum(assay(thyr_se, "unstranded")),
  stranded_first = sum(assay(thyr_se, "stranded_first")),
  stranded_second = sum(assay(thyr_se, "stranded_second"))
)
cat("Total counts by assay:\n")
print(total_counts)
cat("Highest count assay:", names(which.max(total_counts)), "\n")

# Check non-zero gene counts
nonzero_genes <- c(
  unstranded = sum(rowSums(assay(thyr_se, "unstranded")) > 0),
  stranded_first = sum(rowSums(assay(thyr_se, "stranded_first")) > 0),
  stranded_second = sum(rowSums(assay(thyr_se, "stranded_second")) > 0)
)
cat("\nGenes with non-zero counts:\n")
print(nonzero_genes)

# Recommendation
if (names(which.max(total_counts)) == "stranded_second") {
  cat("\n[OK] stranded_second appears to be the appropriate assay (highest counts)\n")
} else {
  cat("\n[WARNING] Note:", names(which.max(total_counts)), "has highest counts, not stranded_second\n")
}

# Quick sample check for downstream verification
cat("\n=== Sample ID check for downstream ===\n")
cat("First 5 sample IDs:\n")
print(head(colnames(thyr_se), 5))
cat("Sample ID format: ", 
    ifelse(all(grepl("^SC", head(colnames(thyr_se), 5))), "SC format (correct)", "Unexpected format"), "\n")

# Store essential info before final cleanup
final_sample_count <- ncol(thyr_se)
final_gene_count <- nrow(thyr_se)
final_assay_names <- names(assays(thyr_se))

# Final cleanup - remove ALL objects except thyr_se and paths
rm(list = setdiff(ls(), c("thyr_se", "paths", "final_sample_count", "final_gene_count", "final_assay_names")))

# Force garbage collection
gc()

cat("\n=== Data Loading Completed ===\n")
cat("Retained in environment:\n")
cat("  - thyr_se: SummarizedExperiment object\n")
cat("    Dimensions:", final_gene_count, "genes x", final_sample_count, "samples\n")
cat("    Assays:", paste(final_assay_names, collapse = ", "), "\n")
cat("    Gene metadata: gene_id, gene_name, gene_type\n")
cat("  - paths: Directory structure from setup.R\n")
cat("\nAll other objects removed for clean environment.\n")
cat("\nFor downstream analysis:\n")
cat("  - Load SE: thyr_se <- readRDS('analysis_v7/data/processed/thyr_se_raw.rds')\n")
cat("  - Access counts: assay(thyr_se, 'stranded_second')\n")
cat("  - Access metadata: colData(thyr_se)\n")
cat("  - Access gene info: rowData(thyr_se)\n")

# Clean up the final info variables
rm(final_sample_count, final_gene_count, final_assay_names)