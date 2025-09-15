# 01_data_download.R - REBC-THYR Data Download (v7.4)
# Purpose: Create manifest and download Gene Expression Quantification files
# v7.4: Robust download with gdcdata return value verification

source("analysis_v7/setup.R")

cat("\n=== REBC-THYR Data Download (v7.4) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# ============================================================================
# Part 1: Create Manifest
# ============================================================================

cat("\n--- Part 1: Creating Manifest ---\n")

# Query - Gene Expression Quantification のみ (open access)
q <- GenomicDataCommons::files() %>%
  GenomicDataCommons::filter(
    ~ cases.project.project_id == "REBC-THYR" &
      data_category == "Transcriptome Profiling" &
      experimental_strategy == "RNA-Seq" &
      data_type == "Gene Expression Quantification" &
      analysis.workflow_type == "STAR - Counts" &
      access == "open"
  )

# Get manifest
cat("Creating manifest...\n")
manifest <- GenomicDataCommons::manifest(q)

# Validate manifest
if (nrow(manifest) == 0) stop("No files found in manifest")
if (!("id" %in% names(manifest))) stop("'id' column missing in manifest")

cat("Manifest entries:", nrow(manifest), "\n")

# Save manifest
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
manifest_file <- paste0(paths$raw, "manifest_gene_counts_", timestamp, ".tsv")
data.table::fwrite(manifest, manifest_file, sep = "\t", quote = FALSE)
cat("Manifest saved:", basename(manifest_file), "\n")

# ============================================================================
# Part 2: Download Files with Return Value Tracking
# ============================================================================

cat("\n--- Part 2: Downloading Files ---\n")

# Create GDC directory if needed
if (!dir.exists(paths$gdc)) {
  dir.create(paths$gdc, recursive = TRUE)
  cat("Created directory:", paths$gdc, "\n")
}

file_ids <- manifest$id
n_total <- length(file_ids)

# Check for existing files in project directory
existing_in_gdc <- list.dirs(paths$gdc, full.names = FALSE, recursive = FALSE)
files_to_download <- setdiff(file_ids, existing_in_gdc)
n_to_download <- length(files_to_download)

cat("\nDownload plan:\n")
cat("  Total files:", n_total, "\n")
cat("  Already in project:", length(existing_in_gdc), "\n")
cat("  To download:", n_to_download, "\n")

if (n_to_download > 0) {
  cat("\nStarting download...\n")
  cat("Note: Files will be downloaded to cache, then copied to project\n")
  cat("Estimated time:", round(n_to_download * 3 / 60, 1), "minutes\n\n")
  
  # Initialize tracking
  download_start_time <- Sys.time()
  successful <- character()
  failed <- character()
  download_paths <- list()
  
  # Batch processing for progress display
  batch_size <- 10
  n_batches <- ceiling(n_to_download / batch_size)
  
  for (batch_idx in seq_len(n_batches)) {
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, n_to_download)
    batch_ids <- files_to_download[start_idx:end_idx]
    
    cat(sprintf("Batch %d/%d (files %d-%d of %d)\n", 
                batch_idx, n_batches, start_idx, end_idx, n_to_download))
    
    for (file_id in batch_ids) {
      cat("  ", file_id, "... ")
      
      # Check if already in project directory
      target_dir <- file.path(paths$gdc, file_id)
      if (dir.exists(target_dir) && length(list.files(target_dir)) > 0) {
        cat("already in project (skipping)\n")
        successful <- c(successful, file_id)
        next
      }
      
      # Try to download
      result <- try(
        GenomicDataCommons::gdcdata(
          file_id,
          destination_dir = paths$gdc,  # Will be ignored, goes to cache
          progress = FALSE
        ),
        silent = TRUE
      )
      
      if (inherits(result, "try-error")) {
        cat("FAILED\n")
        failed <- c(failed, file_id)
      } else {
        # Store the returned path(s)
        download_paths[[file_id]] <- as.character(result)
        
        # Copy from cache to project directory
        cached_path <- as.character(result)[1]
        if (file.exists(cached_path)) {
          # Create target directory
          if (!dir.exists(target_dir)) {
            dir.create(target_dir, recursive = TRUE)
          }
          
          # Copy file
          target_file <- file.path(target_dir, basename(cached_path))
          if (file.copy(cached_path, target_file, overwrite = TRUE)) {
            cat("done (copied to project)\n")
            successful <- c(successful, file_id)
          } else {
            cat("FAILED (copy error)\n")
            failed <- c(failed, file_id)
          }
        } else {
          cat("FAILED (file not found)\n")
          failed <- c(failed, file_id)
        }
      }
    }
    
    # Progress report
    elapsed_time <- difftime(Sys.time(), download_start_time, units = "mins")
    completed <- length(successful) + length(failed)
    cat(sprintf("\nProgress: %d/%d completed (%.1f%%) - Elapsed: %.1f minutes\n\n",
                completed, n_to_download,
                completed / n_to_download * 100,
                as.numeric(elapsed_time)))
  }
  
  # Final report
  cat("=== Download Complete ===\n")
  cat("Successful:", length(successful), "\n")
  cat("Failed:", length(failed), "\n")
  
  # Save retry manifest if needed
  if (length(failed) > 0) {
    retry_manifest <- manifest[manifest$id %in% failed, ]
    retry_file <- paste0(paths$raw, "manifest_retry_", timestamp, ".tsv")
    data.table::fwrite(retry_manifest, retry_file, sep = "\t", quote = FALSE)
    cat("Retry manifest saved:", basename(retry_file), "\n")
  }
  
  total_time <- difftime(Sys.time(), download_start_time, units = "mins")
  cat(sprintf("Total time: %.1f minutes\n", as.numeric(total_time)))
  
} else {
  cat("\n✓ All files already present in project directory!\n")
}

# ============================================================================
# Part 3: Final Verification
# ============================================================================

cat("\n--- Part 3: Final Verification ---\n")

# Verify files in project directory
project_dirs <- list.dirs(paths$gdc, full.names = FALSE, recursive = FALSE)
project_files <- list.files(paths$gdc, pattern = "\\.tsv$", recursive = TRUE)

cat("Project directory status:\n")
cat("  Directories (file IDs):", length(project_dirs), "\n")
cat("  TSV files:", length(project_files), "\n")

# Check completeness
if (length(project_dirs) == n_total) {
  cat("\n✓ All", n_total, "files successfully downloaded and organized!\n")
} else {
  missing_count <- n_total - length(project_dirs)
  cat("\n⚠ Warning:", missing_count, "files may be missing\n")
  cat("Re-run this script to retry missing files\n")
}

# Show cache location if files were downloaded
if (exists("download_paths") && length(download_paths) > 0) {
  cache_dirs <- unique(dirname(unlist(download_paths)))
  if (length(cache_dirs) > 0) {
    cat("\nNote: Original downloads cached in:\n")
    for (d in head(cache_dirs, 3)) cat("  ", d, "\n")
  }
}

# ============================================================================
# Complete cleanup - Remove ALL variables except 'paths'
# ============================================================================

# Final summary before cleanup
cat(sprintf(
  "\n=== Final Summary ===\nManifest: %s\nTotal files: %d\nProject directory: %s\n",
  basename(manifest_file), n_total, paths$gdc
))

# Count and cleanup
cleanup_count <- length(setdiff(ls(), "paths"))
rm(list = setdiff(ls(), c("paths", "cleanup_count")))
cat(sprintf("\nEnvironment cleanup: %d objects removed, keeping only 'paths'\n", cleanup_count))
rm(cleanup_count)

cat("\nNext: Run 02_data_loading.R\n")