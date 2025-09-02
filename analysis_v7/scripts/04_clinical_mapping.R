# 04_clinical_mapping.R - Clinical Data Integration (v7)
# Purpose: Integrate Science paper supplemental data and map samples

source("analysis_v7/setup.R")

cat("\n=== Clinical Data Integration ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# ---- 1. Load existing data ----
cat("\nLoading existing data...\n")

# SummarizedExperiment
se <- readRDS(paste0(paths$processed, "thyr_se_raw.rds"))
cat("SE loaded:", ncol(se), "files\n")

# Sample lists from 03
sample_lists <- readRDS(paste0(paths$processed, "sample_lists_v7.rds"))
cat("Sample lists loaded:\n")
cat("  Total samples:", nrow(sample_lists$sample_data), "\n")
cat("  Tumor samples:", length(sample_lists$tumor_ids), "\n")
cat("  Normal samples:", length(sample_lists$normal_ids), "\n")

# ---- 2. Download and load Science supplemental data ----
cat("\nLoading Science supplemental data...\n")
science_file <- paste0(paths$raw, "abg2538-data-s1.txt")

if(!file.exists(science_file)) {
  cat("Downloading Science supplemental data...\n")
  download.file(
    url = "https://www.science.org/doi/suppl/10.1126/science.abg2538/suppl_file/abg2538-data-s1.txt",
    destfile = science_file,
    method = "curl"
  )
  cat("Download completed: ", science_file, "\n")
}

clinical_data <- data.table::fread(science_file)
cat("Clinical data loaded:", nrow(clinical_data), "cases\n")

if("ERK_REBC_DEgenes" %in% colnames(clinical_data)) {
  cat("ERR proxy column found: ERK_REBC_DEgenes\n")
}

# ---- 3. Map file IDs to sample/case IDs (VECTORIZED) ----
cat("\nMapping file IDs to samples and cases...\n")

# Load manifest
manifest_files <- list.files(
  paths$raw,
  pattern = "^manifest_gene_counts_.*\\.tsv$",
  full.names = TRUE
)
manifest <- data.table::fread(manifest_files[length(manifest_files)])
cat("Manifest loaded:", nrow(manifest), "files\n")

# Create file mapping
file_mapping <- data.frame(
  file_id = manifest$id,
  file_name = manifest$file_name,
  submitter_id = manifest$submitter_id,
  stringsAsFactors = FALSE
)

# ---- 4. VECTORIZED matching with sample data ----
cat("\nMatching files to samples (vectorized)...\n")

# Direct merge instead of loop - MUCH FASTER!
file_mapping_matched <- merge(
  file_mapping,
  sample_lists$sample_data,
  by.x = "submitter_id",
  by.y = "sample_id",
  all.x = TRUE
)

# Check matching results
n_matched <- sum(!is.na(file_mapping_matched$case_id))
cat("Matched", n_matched, "out of", nrow(file_mapping_matched), "files\n")

if(n_matched == 0) {
  cat("\nNo exact matches. Checking ID formats...\n")
  cat("First 5 submitter_ids from manifest:\n")
  print(head(file_mapping$submitter_id, 5))
  cat("\nFirst 5 sample IDs from sample_data:\n")
  print(head(sample_lists$sample_data$sample_id, 5))
}

# ---- 5. Integrate clinical data (VECTORIZED) ----
cat("\nIntegrating clinical data...\n")

# Merge with clinical data
integrated_data <- merge(
  file_mapping_matched,
  clinical_data,
  by.x = "case_id",
  by.y = "REBC_ID",
  all.x = TRUE
)

cat("Clinical data matched for", sum(!is.na(integrated_data$SEX)), "files\n")

# ---- 6. Identify tumor/normal pairs (VECTORIZED) ----
cat("\nIdentifying tumor/normal pairs...\n")

# Use data.table for speed
dt_mapping <- data.table::as.data.table(file_mapping_matched)

cases_with_pairs <- dt_mapping[!is.na(case_id), .(
  n_files = .N,
  has_tumor = any(grepl("Tumor", sample_type)),
  has_normal = any(grepl("Normal", sample_type))
), by = case_id][has_tumor == TRUE & has_normal == TRUE]

cat("Cases with tumor/normal pairs:", nrow(cases_with_pairs), "\n")

# ---- 7. Create paired sample list ----
cat("\nCreating paired sample list...\n")

# Get tumor and normal files for each paired case
paired_samples <- dt_mapping[case_id %in% cases_with_pairs$case_id]

tumor_files <- paired_samples[grepl("Tumor", sample_type), .(
  case_id,
  tumor_file_id = file_id,
  tumor_sample_id = submitter_id
)]

normal_files <- paired_samples[grepl("Normal", sample_type), .(
  case_id,
  normal_file_id = file_id,
  normal_sample_id = submitter_id
)]

# Join to create pairs
paired_files <- merge(tumor_files, normal_files, by = "case_id", allow.cartesian = TRUE)
cat("Total paired file combinations:", nrow(paired_files), "\n")

# ---- 8. Save results ----
cat("\nSaving results...\n")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

integrated_metadata <- list(
  file_mapping = file_mapping_matched,
  clinical_data = clinical_data,
  integrated_data = integrated_data,
  cases_with_pairs = cases_with_pairs,
  paired_files = paired_files,
  sample_lists = sample_lists
)

output_file <- paste0(paths$processed, "integrated_metadata_", timestamp, ".rds")
saveRDS(integrated_metadata, output_file)
saveRDS(integrated_metadata, paste0(paths$processed, "integrated_metadata.rds"))

cat("\nSaved to:", output_file, "\n")

# ---- 9. Summary report ----
cat("\n=== Summary ===\n")
cat("Total files:", nrow(file_mapping_matched), "\n")
cat("Files matched to samples:", n_matched, "\n")
cat("Files with clinical data:", sum(!is.na(integrated_data$SEX)), "\n")
cat("Cases with tumor/normal pairs:", nrow(cases_with_pairs), "\n")
cat("Paired file combinations:", nrow(paired_files), "\n")

cat("\n=== Clinical Data Integration Completed ===\n")
