# 06_err_integration.R - Filter merged samples and integrate ERR
# Purpose: Filter to merged samples, add ERR, define final groups

source("analysis_v7/setup.R")

cat("\n=== ERR Integration and Sample Filtering (v7) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# Load data
se <- readRDS(paste0(paths$processed, "thyr_se_braf_ret.rds"))
metadata <- as.data.frame(colData(se))

cat("\nInput:", ncol(se), "samples\n")

# Filter to merged samples only
merged_mask <- grepl("merged", metadata$sample_id)
se_merged <- se[, merged_mask]
metadata_merged <- metadata[merged_mask, ]

# 列名が確実に設定されていることを確認
if(is.null(colnames(se_merged))) {
  colnames(se_merged) <- metadata_merged$sample_id
}

cat("After merged filter:", ncol(se_merged), "samples\n")

# Load ERR data
library(readxl)
err_braf <- read_excel(paste0(paths$raw, "thyr_poc_braf_mut.xlsx"))
err_ret <- read_excel(paste0(paths$raw, "thyr_poc_ret_fusion.xlsx"))

# Combine ERR data
err_combined <- rbind(err_braf, err_ret)
cat("\nERR data loaded:", nrow(err_combined), "cases\n")

# Merge ERR with metadata
metadata_err <- merge(
  metadata_merged,
  err_combined[, c("case_id", "ERR")],
  by = "case_id",
  all.x = TRUE
)

# Identify paired samples
paired_cases <- metadata_err %>%
  dplyr::group_by(case_id) %>%
  dplyr::summarise(
    has_tumor = any(is_tumor),
    has_normal = any(is_normal),
    has_pair = has_tumor & has_normal,
    .groups = "drop"
  )

metadata_err <- merge(metadata_err, paired_cases[, c("case_id", "has_pair")], by = "case_id")

# Define groups
metadata_err$group <- NA
metadata_err$group_type <- NA

# R0/R1 for RET
ret_paired <- metadata_err$case_driver == "RET" & metadata_err$has_pair
metadata_err$group[ret_paired & is.na(metadata_err$ERR)] <- "R0"
metadata_err$group[ret_paired & !is.na(metadata_err$ERR) & metadata_err$ERR >= 66.6] <- "R1"

# B0/B1 for BRAF
braf_paired <- metadata_err$case_driver == "BRAF" & metadata_err$has_pair
metadata_err$group[braf_paired & is.na(metadata_err$ERR)] <- "B0"
metadata_err$group[braf_paired & !is.na(metadata_err$ERR) & metadata_err$ERR >= 66.6] <- "B1"

# Marker test samples
intermediate_err <- !is.na(metadata_err$ERR) & metadata_err$ERR > 0 & metadata_err$ERR < 66.6
unpaired <- !metadata_err$has_pair
metadata_err$group_type[intermediate_err | unpaired] <- "marker_test"
metadata_err$group_type[!is.na(metadata_err$group)] <- "main_analysis"

# Summary
cat("\n=== Group Assignment ===\n")
group_summary <- table(metadata_err$group, metadata_err$sample_type, useNA = "ifany")
print(group_summary)

cat("\n=== Sample Categories ===\n")
cat("Main analysis (R0/R1/B0/B1 with pairs):", sum(metadata_err$group_type == "main_analysis"), "\n")
cat("Marker test (intermediate ERR or unpaired):", sum(metadata_err$group_type == "marker_test"), "\n")

# Update SE
colData(se_merged) <- DataFrame(metadata_err)

# 列名の最終確認（metadata_errの行順序に合わせる）
if(is.null(colnames(se_merged))) {
  colnames(se_merged) <- metadata_err$sample_id
}

# rownames も設定
rownames(colData(se_merged)) <- metadata_err$sample_id

# Save
saveRDS(se_merged, paste0(paths$processed, "thyr_se_merged_err.rds"))
save(se_merged, file = paste0(paths$processed, "thyr_se_merged_err.RData"))

data.table::fwrite(
  metadata_err,
  paste0(paths$processed, "merged_err_metadata.tsv"),
  sep = "\t"
)

cat("\n=== Completed ===\n")
cat("Files created:\n")
cat("  - thyr_se_merged_err.rds\n")
cat("  - merged_err_metadata.tsv\n")
