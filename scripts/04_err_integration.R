# ==============================================================================
# REBC-THYR ERR Integration and Sample Extraction Script (Corrected)
# 04_err_integration.R
# ==============================================================================

# Required libraries
library(readxl)
library(dplyr)
library(tidyr)

# Assuming previous data objects are loaded:
# - se_thyr (SummarizedExperiment object)

# Load required data if not already available
if (!exists("sample_data")) {
  cat("Loading sample metadata...\n")
  sample_data <- read.delim("./data/raw/sample.tsv", stringsAsFactors = FALSE)
}

if (!exists("abg_data")) {
  cat("Loading clinical data...\n")  
  abg_data <- read.delim("./data/raw/abg2538-data-s1.txt", stringsAsFactors = FALSE)
}

cat("Starting ERR integration and sample extraction (corrected)...\n")

# ==============================================================================
# 1. Identify RNA-seq Samples and Paired Cases
# ==============================================================================

cat("Identifying RNA-seq samples and paired cases...\n")

# Get actual RNA-seq sample IDs
rna_samples <- colnames(se_thyr)
cat("Total RNA-seq samples:", length(rna_samples), "\n")

# Filter sample_data to RNA-seq samples only
sample_in_rna <- sample_data$sample_submitter_id %in% rna_samples
rna_sample_data <- sample_data[sample_in_rna, ]
cat("RNA-seq sample records:", nrow(rna_sample_data), "\n")

# Identify cases with both tumor and normal samples
rna_case_summary <- rna_sample_data %>%
  group_by(case_submitter_id, sample_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = sample_type, values_from = count, values_fill = 0)

# Get paired cases (both Primary Tumor and Solid Tissue Normal)
paired_cases <- rna_case_summary[
  rna_case_summary$`Primary Tumor` >= 1 & 
    rna_case_summary$`Solid Tissue Normal` >= 1, ]

cat("Cases with tumor+normal pairs:", nrow(paired_cases), "\n")

# Create case-sample mapping for paired cases
paired_case_samples <- rna_sample_data[
  rna_sample_data$case_submitter_id %in% paired_cases$case_submitter_id, 
  c("case_submitter_id", "sample_submitter_id", "sample_type")
]

# Separate tumor and normal samples
tumor_map <- paired_case_samples[paired_case_samples$sample_type == "Primary Tumor", 
                                 c("case_submitter_id", "sample_submitter_id")]
names(tumor_map) <- c("REBC_ID", "tumor_sample_id")

normal_map <- paired_case_samples[paired_case_samples$sample_type == "Solid Tissue Normal", 
                                  c("case_submitter_id", "sample_submitter_id")]
names(normal_map) <- c("REBC_ID", "normal_sample_id")

# Merge tumor and normal mappings
case_sample_map <- merge(tumor_map, normal_map, by = "REBC_ID")

cat("Final paired cases with sample mapping:", nrow(case_sample_map), "\n")

# ==============================================================================
# 2. Load and Process ERR Data with Strict Criteria
# ==============================================================================

cat("Loading ERR data files...\n")

# Load RET and BRAF ERR data
ret_err <- readxl::read_xlsx("./data/raw/thyr_poc_ret_fusion.xlsx")
braf_err <- readxl::read_xlsx("./data/raw/thyr_poc_braf_mut.xlsx")

# Standardize column names
colnames(ret_err)[1] <- "REBC_ID"
colnames(braf_err)[1] <- "REBC_ID"

cat("Original ERR data:\n")
cat("  RET:", nrow(ret_err), "cases\n")
cat("  BRAF:", nrow(braf_err), "cases\n")

# Apply strict criteria for BRAF group
cat("Applying strict criteria for BRAF group...\n")

# First, load abg_data if not already loaded
if (!exists("abg_data")) {
  cat("Loading clinical data for BRAF filtering...\n")
  abg_data <- read.delim("./data/raw/abg2538-data-s1.txt", stringsAsFactors = FALSE)
}

# Define strict BRAF criteria
braf_strict_ids <- abg_data[
  (abg_data$WGS_CandidateDriverMutation == "BRAF" | 
     abg_data$WGS_CandidateDriverMutation == "") &
    (abg_data$RNA_CandidateDriverMutation == "BRAF" | 
       abg_data$RNA_CandidateDriverMutation == "") &
    abg_data$Designated_Driver == "BRAF.MutV600E", 
  "REBC_ID"]

cat("Strict BRAF criteria results:\n")
cat("  Total strict BRAF cases in abg_data:", length(braf_strict_ids), "\n")

# Filter BRAF ERR data to strict criteria only
braf_err_original <- braf_err
braf_err <- braf_err[braf_err$REBC_ID %in% braf_strict_ids, ]

cat("  BRAF ERR data after strict filtering:", nrow(braf_err), "/", nrow(braf_err_original), "\n")
cat("  Excluded BRAF cases:", nrow(braf_err_original) - nrow(braf_err), "\n")

# RET group uses original criteria (Fusion.RET is already specific)
cat("RET group (no additional filtering needed):\n")
cat("  RET ERR data:", nrow(ret_err), "cases\n")

# Verify ERR cases are within paired RNA-seq cases
ret_in_paired <- sum(ret_err$REBC_ID %in% case_sample_map$REBC_ID)
braf_in_paired <- sum(braf_err$REBC_ID %in% case_sample_map$REBC_ID)

cat("ERR cases with RNA-seq pairs:\n")
cat("  RET:", ret_in_paired, "/", nrow(ret_err), "\n")
cat("  BRAF (strict):", braf_in_paired, "/", nrow(braf_err), "\n")

# ==============================================================================
# 3. ERR-based Group Classification
# ==============================================================================

cat("Performing ERR-based group classification...\n")

# Define ERR cutoff (66.6%)
err_cutoff <- 66.6

# Classify RET samples
ret_err$group <- ifelse(
  is.na(ret_err$ERR), 
  "R0",  # RET_Unexposed
  ifelse(ret_err$ERR >= err_cutoff, "R1", "Exclude")  # RET_RadHigh or Exclude
)

# Classify BRAF samples  
braf_err$group <- ifelse(
  is.na(braf_err$ERR),
  "B0",  # BRAF_Unexposed
  ifelse(braf_err$ERR >= err_cutoff, "B1", "Exclude")  # BRAF_RadHigh or Exclude
)

# Filter to analysis groups only (exclude low ERR cases)
ret_analysis <- ret_err[ret_err$group != "Exclude", ]
braf_analysis <- braf_err[braf_err$group != "Exclude", ]

# Combine analysis data
analysis_data <- rbind(
  ret_analysis[, c("REBC_ID", "ERR", "group")],
  braf_analysis[, c("REBC_ID", "ERR", "group")]
)

# Summary of classification
cat("Group classification summary:\n")
group_summary <- table(analysis_data$group)
print(group_summary)

cat("Excluded cases (ERR < 66.6% but not NA):\n")
cat("  RET:", sum(ret_err$group == "Exclude"), "cases\n")
cat("  BRAF (strict criteria):", sum(braf_err$group == "Exclude"), "cases\n")

# Show detailed breakdown for BRAF strict filtering
cat("BRAF filtering summary:\n")
cat("  Original BRAF ERR cases:", nrow(braf_err_original), "\n")
cat("  Cases meeting strict criteria:", nrow(braf_err), "\n")
cat("  Cases removed by strict criteria:", nrow(braf_err_original) - nrow(braf_err), "\n")
cat("  B0 candidates (ERR=NA):", sum(is.na(braf_err$ERR)), "\n")
cat("  B1 candidates (ERR≥66.6%):", sum(braf_err$ERR >= err_cutoff, na.rm = TRUE), "\n")
cat("  Excluded (0%<ERR<66.6%):", sum(!is.na(braf_err$ERR) & braf_err$ERR < err_cutoff), "\n")

# ==============================================================================
# 4. Merge with Sample Mapping
# ==============================================================================

cat("Merging ERR groups with sample mapping...\n")

# Merge analysis data with case-sample mapping
final_data <- merge(analysis_data, case_sample_map, by = "REBC_ID", all.x = TRUE)

# Check for successful merging
successful_merge <- !is.na(final_data$tumor_sample_id) & !is.na(final_data$normal_sample_id)
cat("Cases successfully merged with sample pairs:", sum(successful_merge), "/", nrow(final_data), "\n")

# Filter to successfully merged cases
final_data_complete <- final_data[successful_merge, ]

# Final group summary
cat("Final groups with complete sample pairs:\n")
final_summary <- table(final_data_complete$group)
print(final_summary)

# ==============================================================================
# 5. Create Sample Lists by Group
# ==============================================================================

cat("Creating sample lists by group...\n")

# Initialize sample lists
sample_lists <- list()

for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(final_summary) && final_summary[group] > 0) {
    # Get cases for this group
    group_data <- final_data_complete[final_data_complete$group == group, ]
    
    # Extract tumor and normal sample IDs
    tumor_samples <- group_data$tumor_sample_id
    normal_samples <- group_data$normal_sample_id
    
    # Store as paired lists
    sample_lists[[group]] <- list(
      tumor = tumor_samples,
      normal = normal_samples,
      cases = group_data$REBC_ID
    )
    
    cat(sprintf("%s: %d cases (%d tumor + %d normal samples)\n", 
                group, length(tumor_samples), length(tumor_samples), length(normal_samples)))
  } else {
    sample_lists[[group]] <- list(
      tumor = character(0),
      normal = character(0),
      cases = character(0)
    )
    cat(sprintf("%s: 0 cases\n", group))
  }
}

# ==============================================================================
# 6. Verify Sample Availability in RNA-seq Data
# ==============================================================================

cat("Verifying sample availability in RNA-seq data...\n")

for (group in c("R0", "R1", "B0", "B1")) {
  if (length(sample_lists[[group]]$tumor) > 0) {
    tumor_in_rna <- sum(sample_lists[[group]]$tumor %in% rna_samples)
    normal_in_rna <- sum(sample_lists[[group]]$normal %in% rna_samples)
    
    cat(sprintf("%s: %d/%d tumor, %d/%d normal samples in RNA-seq data\n",
                group, tumor_in_rna, length(sample_lists[[group]]$tumor),
                normal_in_rna, length(sample_lists[[group]]$normal)))
  }
}

# ==============================================================================
# 7. Save Results
# ==============================================================================

cat("Saving results...\n")

# Create processed directory if it doesn't exist
if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Save sample lists
save(sample_lists, file = "./data/processed/sample_lists.rda")

# Save detailed data for reference
err_integration_data <- list(
  final_data_complete = final_data_complete,
  case_sample_map = case_sample_map,
  group_summary = final_summary,
  err_cutoff = err_cutoff,
  excluded_counts = list(
    RET = sum(ret_err$group == "Exclude"),
    BRAF_strict = sum(braf_err$group == "Exclude"),
    BRAF_removed_by_strict_criteria = nrow(braf_err_original) - nrow(braf_err)
  ),
  rna_paired_cases = nrow(paired_cases)
)
save(err_integration_data, file = "./data/processed/err_integration_data.rda")

cat("Results saved to ./data/processed/\n")

# ==============================================================================
# 8. Final Summary
# ==============================================================================

cat("==============================================\n")
cat("ERR Integration Summary (Corrected)\n")
cat("==============================================\n")

cat("RNA-seq Data Overview:\n")
cat(sprintf("  Total RNA-seq samples: %d\n", length(rna_samples)))
cat(sprintf("  Cases with tumor+normal pairs: %d\n", nrow(paired_cases)))

cat("Data Processing Summary:\n")
cat(sprintf("  Original BRAF ERR cases: %d\n", nrow(braf_err_original)))
cat(sprintf("  BRAF cases after strict criteria: %d\n", nrow(braf_err)))
cat(sprintf("  Cases removed by strict criteria: %d\n", nrow(braf_err_original) - nrow(braf_err)))

cat("\nGroup Definitions:\n")
cat("  R0 (RET_Unexposed): ERR = NA\n")
cat("  R1 (RET_RadHigh): ERR >= 66.6%\n")
cat("  B0 (BRAF_Unexposed): ERR = NA + Strict criteria\n")
cat("  B1 (BRAF_RadHigh): ERR >= 66.6% + Strict criteria\n")

cat("\nStrict BRAF Criteria:\n")
cat("  WGS_CandidateDriverMutation: 'BRAF' or empty\n")
cat("  RNA_CandidateDriverMutation: 'BRAF' or empty\n")
cat("  Designated_Driver: 'BRAF.MutV600E'\n")

cat("\nFinal Sample Counts:\n")
for (group in c("R0", "R1", "B0", "B1")) {
  count <- length(sample_lists[[group]]$tumor)
  driver <- ifelse(grepl("R", group), "RET", "BRAF")
  exposure <- ifelse(grepl("0", group), "Unexposed", "RadHigh")
  cat(sprintf("  %s (%s_%s): %d cases\n", group, driver, exposure, count))
}

total_cases <- sum(sapply(sample_lists, function(x) length(x$tumor)))
cat(sprintf("\nTotal cases for analysis: %d\n", total_cases))

# Check minimum sample requirements
min_samples <- 5
sufficient_groups <- sapply(sample_lists, function(x) length(x$tumor) >= min_samples)
cat("\nSample size check (minimum 5 per group):\n")
for (group in names(sufficient_groups)) {
  status <- ifelse(sufficient_groups[group], "✅ PASS", "❌ INSUFFICIENT")
  cat(sprintf("  %s: %s\n", group, status))
}

if (all(sufficient_groups)) {
  cat("\n✅ All groups meet minimum sample requirements for analysis\n")
} else {
  insufficient_groups <- names(sufficient_groups)[!sufficient_groups]
  cat(sprintf("\n⚠️  Insufficient samples in: %s\n", paste(insufficient_groups, collapse = ", ")))
}

cat("\nSample Lists Structure:\n")
cat("  Each group contains:\n")
cat("    $tumor: tumor sample IDs\n")
cat("    $normal: normal sample IDs\n") 
cat("    $cases: REBC_ID case identifiers\n")

cat("\nNext steps:\n")
cat("1. Tumor purity calculation and filtering (ContamDE)\n")
cat("2. PCA analysis for outlier detection (using TPM data)\n") 
cat("3. DEG analysis (primary comparison: R0 vs R1)\n")

cat("ERR integration completed successfully!\n")
cat("==============================================\n")

