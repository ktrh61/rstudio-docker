# 03_download_samples.R - Download sample metadata

source("analysis_v7/setup.R")

cat("\n=== Downloading REBC-THYR Sample Data ===\n")

library(GenomicDataCommons)

# NULLハンドリング
`%||%` <- function(x, y) if(is.null(x)) y else x

# Cases経由でサンプル情報を取得
cases_result <- GenomicDataCommons::cases() %>%
  GenomicDataCommons::filter(~ project.project_id == "REBC-THYR") %>%
  GenomicDataCommons::expand("samples") %>%
  GenomicDataCommons::results_all()

cat("Cases found:", length(cases_result$id), "\n")

# samplesはdata.frameのリスト
samples_list <- cases_result$samples
cat("Sample data frames:", length(samples_list), "\n")

# 全サンプルを結合
sample_data_list <- list()
for(i in seq_along(samples_list)) {
  df <- samples_list[[i]]
  if(!is.null(df) && nrow(df) > 0) {
    # case_idを追加
    df$case_id <- cases_result$submitter_id[i]
    sample_data_list[[i]] <- df[, c("submitter_id", "sample_type", "sample_type_id", "case_id")]
  }
}

# 結合
sample_data <- data.table::rbindlist(sample_data_list, fill = TRUE)
names(sample_data)[1] <- "sample_id"  # submitter_idをsample_idに

cat("Total samples:", nrow(sample_data), "\n")
cat("\nSample type distribution:\n")
print(table(sample_data$sample_type))

# 保存
output_file <- paste0(paths$raw, "sample_metadata_", Sys.Date(), ".tsv")
data.table::fwrite(sample_data, output_file, sep = "\t")
cat("\nSaved to:", output_file, "\n")

# リスト作成
tumor_samples <- sample_data[grepl("Tumor", sample_data$sample_type), ]$sample_id
normal_samples <- sample_data[grepl("Normal", sample_data$sample_type), ]$sample_id
cat("\nTumor samples:", length(tumor_samples), "\n")
cat("Normal samples:", length(normal_samples), "\n")

# RDS保存
saveRDS(list(
  sample_data = sample_data,
  tumor_ids = tumor_samples,
  normal_ids = normal_samples
), paste0(paths$processed, "sample_lists_v7.rds"))
