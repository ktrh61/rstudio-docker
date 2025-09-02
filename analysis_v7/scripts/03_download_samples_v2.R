# 03_download_samples_v2.R - Download complete sample metadata
# Purpose: Get full sample metadata matching GDC Portal sample.tsv

source("analysis_v7/setup.R")

cat("\n=== Downloading Complete REBC-THYR Sample Data ===\n")

library(GenomicDataCommons)

# サンプル情報を詳細に取得
cat("Querying GDC for detailed sample information...\n")

# samples エンドポイントから直接取得
samples_query <- GenomicDataCommons::samples() %>%
  GenomicDataCommons::filter(~ project.project_id == "REBC-THYR") %>%
  GenomicDataCommons::expand(c("cases")) %>%
  GenomicDataCommons::results_all()

cat("Samples found:", length(samples_query$id), "\n")

# データフレームに整形
sample_data <- data.frame(
  sample_id = samples_query$id,
  sample_submitter_id = samples_query$submitter_id,
  sample_type = samples_query$sample_type,
  sample_type_id = samples_query$sample_type_id,
  stringsAsFactors = FALSE
)

# cases情報を追加
cat("\nExtracting case information...\n")
for(i in 1:nrow(sample_data)) {
  cases_info <- samples_query$cases[[i]]
  if(!is.null(cases_info) && nrow(cases_info) > 0) {
    sample_data$case_id[i] <- cases_info$id[1]
    sample_data$case_submitter_id[i] <- cases_info$submitter_id[1]
  } else {
    sample_data$case_id[i] <- NA
    sample_data$case_submitter_id[i] <- NA
  }
}

cat("Sample data prepared:", nrow(sample_data), "rows\n")

# サンプルタイプの分布
cat("\nSample type distribution:\n")
print(table(sample_data$sample_type))

# case_submitter_idの例を表示
cat("\nExample case_submitter_ids:\n")
print(head(unique(sample_data$case_submitter_id), 10))

# sample_submitter_idの例を表示
cat("\nExample sample_submitter_ids:\n")
print(head(sample_data$sample_submitter_id, 10))

# 保存
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# TSV形式で保存（GDC Portal sample.tsvと同等）
output_tsv <- paste0(paths$raw, "samples_complete_", timestamp, ".tsv")
data.table::fwrite(sample_data, output_tsv, sep = "\t")
cat("\nSaved TSV:", output_tsv, "\n")

# RDS形式でも保存
output_rds <- paste0(paths$processed, "samples_complete_", timestamp, ".rds")
saveRDS(sample_data, output_rds)
saveRDS(sample_data, paste0(paths$processed, "samples_complete.rds"))

cat("Saved RDS:", output_rds, "\n")

# リスト形式でも保存（後方互換性）
sample_lists <- list(
  sample_data = sample_data,
  tumor_ids = sample_data$sample_submitter_id[grepl("Tumor", sample_data$sample_type)],
  normal_ids = sample_data$sample_submitter_id[grepl("Normal", sample_data$sample_type)]
)

saveRDS(sample_lists, paste0(paths$processed, "sample_lists_v7_complete.rds"))

cat("\n=== Sample Download Completed ===\n")
cat("Key columns obtained:\n")
cat("  - sample_submitter_id (e.g., SC108351, REBC-ADM8-NB1-A)\n")
cat("  - case_submitter_id (e.g., REBC-ADM8)\n")
cat("These are needed to map files to samples\n")
