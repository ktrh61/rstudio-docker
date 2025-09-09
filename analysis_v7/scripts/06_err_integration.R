# 06_err_integration.R - Filter merged samples and integrate ERR
# Purpose: Filter to merged samples, add ERR, define 10 groups (4 training + 6 test)

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

# Define 10 groups with clear categorization
metadata_err$group <- NA
metadata_err$analysis_role <- NA
metadata_err$test_eligible <- FALSE
metadata_err$test_confidence <- NA

# ========== 1. Main analysis groups (Training) - 4 groups ==========
# R0/R1 for RET (paired only)
ret_paired <- metadata_err$case_driver == "RET" & metadata_err$has_pair
metadata_err$group[ret_paired & is.na(metadata_err$ERR)] <- "R0"
metadata_err$group[ret_paired & !is.na(metadata_err$ERR) & metadata_err$ERR >= 66.6] <- "R1"

# B0/B1 for BRAF (paired only)
braf_paired <- metadata_err$case_driver == "BRAF" & metadata_err$has_pair
metadata_err$group[braf_paired & is.na(metadata_err$ERR)] <- "B0"
metadata_err$group[braf_paired & !is.na(metadata_err$ERR) & metadata_err$ERR >= 66.6] <- "B1"

# Mark training samples
metadata_err$analysis_role[!is.na(metadata_err$group)] <- "training"

# ========== 2. Test groups - 6 groups ==========
# 2a. Intermediate ERR with pairs (highest confidence test) - 2 groups
intermediate_paired <- metadata_err$has_pair & 
  !is.na(metadata_err$ERR) & 
  metadata_err$ERR > 0 & 
  metadata_err$ERR < 66.6

metadata_err$group[intermediate_paired & metadata_err$case_driver == "RET"] <- "R_mid"
metadata_err$group[intermediate_paired & metadata_err$case_driver == "BRAF"] <- "B_mid"

# 2b. Unpaired tumor samples - 4 groups
unpaired_tumor <- !metadata_err$has_pair & metadata_err$is_tumor

# ERR=NA unpaired tumors
metadata_err$group[unpaired_tumor & metadata_err$case_driver == "RET" & is.na(metadata_err$ERR)] <- "R_unp_NA"
metadata_err$group[unpaired_tumor & metadata_err$case_driver == "BRAF" & is.na(metadata_err$ERR)] <- "B_unp_NA"

# ERR≥66.6 unpaired tumors
metadata_err$group[unpaired_tumor & metadata_err$case_driver == "RET" & 
                     !is.na(metadata_err$ERR) & metadata_err$ERR >= 66.6] <- "R_unp_high"
metadata_err$group[unpaired_tumor & metadata_err$case_driver == "BRAF" & 
                     !is.na(metadata_err$ERR) & metadata_err$ERR >= 66.6] <- "B_unp_high"

# Intermediate ERR unpaired tumors (not needed if no samples, but included for completeness)
metadata_err$group[unpaired_tumor & metadata_err$case_driver == "RET" & 
                     !is.na(metadata_err$ERR) & metadata_err$ERR > 0 & metadata_err$ERR < 66.6] <- "R_unp_mid"
metadata_err$group[unpaired_tumor & metadata_err$case_driver == "BRAF" & 
                     !is.na(metadata_err$ERR) & metadata_err$ERR > 0 & metadata_err$ERR < 66.6] <- "B_unp_mid"

# Mark test samples
test_groups <- c("R_mid", "B_mid", "R_unp_NA", "B_unp_NA", 
                 "R_unp_high", "B_unp_high", "R_unp_mid", "B_unp_mid")
metadata_err$analysis_role[metadata_err$group %in% test_groups] <- "test"

# Set test eligibility (tumor samples only)
metadata_err$test_eligible <- metadata_err$is_tumor & metadata_err$analysis_role == "test"

# Set test confidence levels
metadata_err$test_confidence[metadata_err$group %in% c("R_mid", "B_mid")] <- "high"
metadata_err$test_confidence[metadata_err$group %in% c("R_unp_mid", "B_unp_mid")] <- "medium"
metadata_err$test_confidence[metadata_err$group %in% c("R_unp_NA", "B_unp_NA", 
                                                       "R_unp_high", "B_unp_high")] <- "low"

# ========== Summary ==========
cat("\n=== Group Assignment (10 groups) ===\n")

# Training groups summary
training_summary <- metadata_err %>%
  dplyr::filter(analysis_role == "training") %>%
  dplyr::group_by(group, sample_type) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = sample_type, values_from = n, values_fill = 0)

cat("\nTraining Groups (4):\n")
print(as.data.frame(training_summary))

# Test groups summary
test_summary <- metadata_err %>%
  dplyr::filter(analysis_role == "test") %>%
  dplyr::group_by(group, has_pair, test_confidence) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    err_mean = mean(ERR, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nTest Groups (6):\n")
print(as.data.frame(test_summary))

# Overall summary
cat("\n=== Overall Summary ===\n")
cat("Total samples:", nrow(metadata_err), "\n")
cat("Training samples:", sum(metadata_err$analysis_role == "training", na.rm = TRUE), "\n")
cat("Test samples:", sum(metadata_err$analysis_role == "test", na.rm = TRUE), "\n")
cat("Unassigned samples:", sum(is.na(metadata_err$group)), "\n")

# Check for unassigned samples
if (sum(is.na(metadata_err$group)) > 0) {
  cat("\nWarning: Some samples were not assigned to any group!\n")
  unassigned <- metadata_err[is.na(metadata_err$group), 
                             c("sample_id", "case_driver", "sample_type", "has_pair", "ERR")]
  print(head(unassigned))
}

# Group distribution table
group_dist <- table(metadata_err$group, metadata_err$sample_type, useNA = "ifany")
cat("\n=== Detailed Group Distribution ===\n")
print(group_dist)

# Update SE
colData(se_merged) <- DataFrame(metadata_err)

# 列名の最終確認
if(is.null(colnames(se_merged))) {
  colnames(se_merged) <- metadata_err$sample_id
}

# rownames も設定
rownames(colData(se_merged)) <- metadata_err$sample_id

# Save main output
saveRDS(se_merged, paste0(paths$processed, "thyr_se_merged_err.rds"))
save(se_merged, file = paste0(paths$processed, "thyr_se_merged_err.RData"))

data.table::fwrite(
  metadata_err,
  paste0(paths$processed, "merged_err_metadata.tsv"),
  sep = "\t"
)

# Additional output: group summary for reference
group_summary <- metadata_err %>%
  dplyr::group_by(group, case_driver, analysis_role, has_pair, test_confidence) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_cases = dplyr::n_distinct(case_id),
    n_tumor = sum(is_tumor),
    n_normal = sum(is_normal),
    err_min = ifelse(any(!is.na(ERR)), min(ERR, na.rm = TRUE), NA),
    err_mean = mean(ERR, na.rm = TRUE),
    err_max = ifelse(any(!is.na(ERR)), max(ERR, na.rm = TRUE), NA),
    .groups = "drop"
  ) %>%
  dplyr::arrange(case_driver, analysis_role, group)

data.table::fwrite(
  group_summary,
  paste0(paths$processed, "group_summary_10groups.tsv"),
  sep = "\t"
)

# Helper function for downstream analysis
create_group_lists <- function(metadata) {
  list(
    training = list(
      R0 = metadata$sample_id[metadata$group == "R0"],
      R1 = metadata$sample_id[metadata$group == "R1"],
      B0 = metadata$sample_id[metadata$group == "B0"],
      B1 = metadata$sample_id[metadata$group == "B1"]
    ),
    test = list(
      R_mid = metadata$sample_id[metadata$group == "R_mid"],
      B_mid = metadata$sample_id[metadata$group == "B_mid"],
      R_unp_NA = metadata$sample_id[metadata$group == "R_unp_NA"],
      B_unp_NA = metadata$sample_id[metadata$group == "B_unp_NA"],
      R_unp_high = metadata$sample_id[metadata$group == "R_unp_high"],
      B_unp_high = metadata$sample_id[metadata$group == "B_unp_high"],
      R_unp_mid = metadata$sample_id[metadata$group == "R_unp_mid"],
      B_unp_mid = metadata$sample_id[metadata$group == "B_unp_mid"]
    )
  )
}

# Save group lists for easy access
group_lists <- create_group_lists(metadata_err)
saveRDS(group_lists, paste0(paths$processed, "group_lists_10groups.rds"))

cat("\n=== Completed ===\n")
cat("Files created:\n")
cat("  - thyr_se_merged_err.rds (main output)\n")
cat("  - merged_err_metadata.tsv\n")
cat("  - group_summary_10groups.tsv\n")
cat("  - group_lists_10groups.rds\n")
cat("\n12 groups defined:\n")
cat("  Training (4): R0, R1, B0, B1\n")
cat("  Test (6): R_mid, B_mid, R_unp_NA, B_unp_NA, R_unp_high, B_unp_high\n")
if (sum(metadata_err$group %in% c("R_unp_mid", "B_unp_mid")) > 0) {
  cat("  Additional test: R_unp_mid, B_unp_mid\n")
}
cat("  Excluded (2): R_unp_normal, B_unp_normal (unpaired normals)\n")