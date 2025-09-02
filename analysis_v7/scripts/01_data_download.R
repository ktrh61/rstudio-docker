# 01_data_download.R - REBC-THYR Data Download (v7)
# Purpose: Download Gene Expression Quantification only (906 files)

source("analysis_v7/setup.R")

cat("\n=== REBC-THYR Data Download (Filtered) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# Query - Gene Expression Quantification のみ
q <- GenomicDataCommons::files() %>%
  GenomicDataCommons::filter(
    ~ cases.project.project_id == "REBC-THYR" &
      data_category == "Transcriptome Profiling" &
      experimental_strategy == "RNA-Seq" &
      data_type == "Gene Expression Quantification" &  # ← 追加
      analysis.workflow_type == "STAR - Counts"
  ) %>%
  GenomicDataCommons::select(c(
    "file_name",
    "access",
    "data_type",
    "cases.submitter_id",
    "cases.samples.submitter_id"
  ))

# Get manifest
cat("\nCreating filtered manifest...\n")
manifest <- GenomicDataCommons::manifest(q)

# アクセスレベル確認
access_table <- table(manifest$access)
cat("\nAccess levels:\n")
print(access_table)

# openデータのみにフィルタ（推奨）
if ("open" %in% names(access_table)) {
  manifest_open <- manifest[manifest$access == "open", ]
  cat("\nOpen access files:", nrow(manifest_open), "\n")
} else {
  manifest_open <- manifest
}

# Save manifest
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
manifest_file <- paste0(paths$raw, "manifest_gene_counts_", timestamp, ".tsv")
data.table::fwrite(manifest_open, manifest_file, sep = "\t", quote = FALSE)
cat("Manifest saved:", manifest_file, "\n")
cat("Total files to download:", nrow(manifest_open), "\n")

# Metadata保存
metadata_file <- paste0(paths$raw, "metadata_gene_counts_", timestamp, ".rds")
saveRDS(list(
  query_date = Sys.Date(),
  manifest = manifest_open,
  total_files = nrow(manifest_open),
  access_levels = access_table
), metadata_file)

cat("\n=== Download Instructions ===\n")
cat("Run in terminal:\n")
cat(sprintf("gdc-client download -m %s -d %s\n", manifest_file, paths$gdc))
cat("\nExpected: ~906 files (Gene Expression only)\n")
