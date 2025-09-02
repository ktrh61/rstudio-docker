# 04_link_files_samples.R - Link RNA-seq files to samples
# Purpose: Add sample and case metadata to SummarizedExperiment

source("analysis_v7/setup.R")

cat("\n=== Linking Files to Samples (v7) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# Load SE
se <- readRDS(paste0(paths$processed, "thyr_se_raw.rds"))
cat("RNA-seq files loaded:", ncol(se), "\n")

# Get all file IDs
all_file_ids <- colnames(se)

# Fetch file-sample-case mappings from GDC
cat("\nFetching file-sample-case mappings from GDC...\n")
file_sample_map <- GenomicDataCommons::files() %>%
  GenomicDataCommons::filter(~ file_id %in% all_file_ids) %>%
  GenomicDataCommons::expand(c("cases", "cases.samples")) %>%
  GenomicDataCommons::results_all()

cat("Retrieved metadata for", length(file_sample_map$file_id), "files\n")

# Extract sample information with handling for missing fields
cat("Extracting sample information...\n")
file_sample_list <- lapply(seq_along(file_sample_map$file_id), function(i) {
  # Get sample_type, handle missing field
  sample_type_val <- if(length(file_sample_map$cases[[i]]$samples[[1]]$sample_type) > 0) {
    file_sample_map$cases[[i]]$samples[[1]]$sample_type
  } else if(length(file_sample_map$cases[[i]]$samples[[1]]$specimen_type) > 0) {
    # Use specimen_type if sample_type is missing (e.g., Metastatic samples)
    paste(file_sample_map$cases[[i]]$samples[[1]]$specimen_type,
          file_sample_map$cases[[i]]$samples[[1]]$tissue_type)
  } else {
    NA_character_
  }
  
  data.frame(
    file_id = file_sample_map$file_id[i],
    file_name = file_sample_map$file_name[i],
    case_id = file_sample_map$cases[[i]]$case_id,
    case_submitter_id = file_sample_map$cases[[i]]$submitter_id,
    sample_id = file_sample_map$cases[[i]]$samples[[1]]$sample_id,
    sample_submitter_id = file_sample_map$cases[[i]]$samples[[1]]$submitter_id,
    sample_type = sample_type_val,
    stringsAsFactors = FALSE
  )
})

file_sample_df <- data.table::rbindlist(file_sample_list, fill = TRUE)

# Add sample type flags (修正版)
file_sample_df$is_primary_tumor <- grepl("Primary Tumor", file_sample_df$sample_type)
file_sample_df$is_solid_tumor <- grepl("Solid Tissue Tumor", file_sample_df$sample_type)
file_sample_df$is_tumor <- file_sample_df$is_primary_tumor | file_sample_df$is_solid_tumor
file_sample_df$is_normal <- grepl("Solid Tissue Normal", file_sample_df$sample_type)
file_sample_df$is_metastatic <- grepl("Lymphoid Tumor", file_sample_df$sample_type)

# Summary (修正版)
cat("\\nMapping summary:\\n")
cat("  Primary Tumor:", sum(file_sample_df$is_primary_tumor), "\\n")
cat("  Solid Tissue Tumor:", sum(file_sample_df$is_solid_tumor), "\\n")
cat("  Normal samples:", sum(file_sample_df$is_normal), "\\n")
cat("  Metastatic (Lymphoid):", sum(file_sample_df$is_metastatic), "\\n")
cat("  Total:", nrow(file_sample_df), "\\n")

# Case summary
case_summary <- file_sample_df[, .(
  n_files = .N,
  n_tumor = sum(is_tumor),
  n_normal = sum(is_normal),
  n_metastatic = sum(is_metastatic),
  has_pair = sum(is_tumor) > 0 & sum(is_normal) > 0
), by = case_submitter_id]

cat("\nCase summary:\n")
cat("  Total cases:", nrow(case_summary), "\n")
cat("  Paired cases (T/N):", sum(case_summary$has_pair), "\n")
cat("  Cases with metastatic:", sum(case_summary$n_metastatic > 0), "\n")

# Update SE colData
colData_df <- data.frame(file_id = all_file_ids, stringsAsFactors = FALSE)
colData_df <- merge(
  colData_df,
  file_sample_df,
  by = "file_id",
  all.x = TRUE,
  sort = FALSE
)

# Maintain original order
colData_df <- colData_df[match(all_file_ids, colData_df$file_id), ]
colData(se) <- DataFrame(colData_df)

# Save outputs
saveRDS(se, paste0(paths$processed, "thyr_se_annotated.rds"))
save(se, file = paste0(paths$processed, "thyr_se_annotated.RData"))

data.table::fwrite(
  file_sample_df,
  paste0(paths$processed, "file_sample_case_mapping.tsv"),
  sep = "\t"
)

data.table::fwrite(
  case_summary,
  paste0(paths$processed, "case_summary.tsv"),
  sep = "\t"
)

cat("\n=== Linking Completed ===\n")
cat("Files created:\n")
cat("  - thyr_se_annotated.rds\n")
cat("  - thyr_se_annotated.RData\n")
cat("  - file_sample_case_mapping.tsv\n")
cat("  - case_summary.tsv\n")
