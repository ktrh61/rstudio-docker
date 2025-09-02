# 03_download_samples.R - Download sample metadata
# Purpose: Get complete sample and case information from GDC

source("analysis_v7/setup.R")

cat("\n=== Downloading REBC-THYR Sample Data ===\n")

library(GenomicDataCommons)

# cases経由でサンプル情報を取得
cat("Querying GDC for sample information via cases endpoint...\n")

cases_result <- GenomicDataCommons::cases() %>%
  GenomicDataCommons::filter(~ project.project_id == "REBC-THYR") %>%
  GenomicDataCommons::expand("samples") %>%
  GenomicDataCommons::results_all()

cat("Cases found:", length(cases_result$id), "\n")

# data.tableを使った高速処理
cat("Extracting sample data (vectorized)...\n")

# samplesリストをdata.tableに変換
sample_data <- data.table::rbindlist(
  lapply(seq_along(cases_result$id), function(i) {
    samples_df <- cases_result$samples[[i]]
    if(!is.null(samples_df) && nrow(samples_df) > 0) {
      # 必要な列のみ選択して、case情報を追加
      data.frame(
        sample_id = samples_df$sample_id,
        sample_submitter_id = samples_df$submitter_id,
        sample_type = samples_df$sample_type,
        sample_type_id = samples_df$sample_type_id,
        case_id = cases_result$id[i],
        case_submitter_id = cases_result$submitter_id[i],
        stringsAsFactors = FALSE
      )
    }
  }),
  fill = TRUE
)

cat("Total samples extracted:", nrow(sample_data), "\n")

cat("\nSample type distribution:\n")
print(table(sample_data$sample_type))

# 例を表示
cat("\nExample case_submitter_ids (REBC-XXX format):\n")
print(head(unique(sample_data$case_submitter_id), 5))

cat("\nExample sample_submitter_ids (e.g., REBC-XXX-XXX-A):\n")
print(head(sample_data$sample_submitter_id, 10))

# 保存
output_tsv <- paste0(paths$raw, "sample_metadata_", Sys.Date(), ".tsv")
data.table::fwrite(sample_data, output_tsv, sep = "\t")
cat("\nSaved TSV:", output_tsv, "\n")

# RDS形式でも保存
sample_lists <- list(
  sample_data = as.data.frame(sample_data),
  tumor_ids = sample_data$sample_submitter_id[grepl("Tumor", sample_data$sample_type)],
  normal_ids = sample_data$sample_submitter_id[grepl("Normal", sample_data$sample_type)]
)

saveRDS(sample_lists, paste0(paths$processed, "sample_lists_v7.rds"))
cat("Saved RDS: sample_lists_v7.rds\n")

cat("\n=== Download Completed ===\n")
cat("Summary:\n")
cat("  Cases:", length(cases_result$id), "\n")
cat("  Samples:", nrow(sample_data), "\n")
