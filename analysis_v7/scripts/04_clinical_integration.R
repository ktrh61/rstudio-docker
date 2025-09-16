# 04_clinical_integration.R - Integrate clinical data and define analysis groups
# Purpose: Create comprehensive case master table with all clinical/POC data
# Input: thyr_se_strand2_nonzero.rds (for pair checking), clinical data, POC data
# Output: thyr_case_master_full.rds (all cases), thyr_case_master.rds (grouped cases)
# Version: v7.7 - Fixed BRAF mask to strictly match BRAF.MutV600E
# Date: 2025-01-20

source("analysis_v7/setup.R")

cat("\n=== Clinical Data Integration and Group Definition (v7.7) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# Load required libraries
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(readxl)
  library(dplyr)
  library(data.table)
})

# ============================================================================
# Load SummarizedExperiment (for pair checking only)
# ============================================================================
cat("\n--- Loading SummarizedExperiment ---\n")

se_path <- paste0(paths$processed, "thyr_se_strand2_nonzero.rds")
if (!file.exists(se_path)) {
  stop("thyr_se_strand2_nonzero.rds not found. Please run 03_prepare_counts.R first.")
}

thyr_se_strand2_nonzero <- readRDS(se_path)
cat("Loaded SE:", format(nrow(thyr_se_strand2_nonzero), big.mark=","), "genes x", 
    format(ncol(thyr_se_strand2_nonzero), big.mark=","), "samples\n")

# Extract sample metadata from SE for pair checking
sample_metadata <- as.data.frame(colData(thyr_se_strand2_nonzero))

# ============================================================================
# Load comprehensive clinical data
# ============================================================================
cat("\n--- Loading clinical data ---\n")

clinical_file <- paste0(paths$raw, "abg2538-data-s1.txt")
if (!file.exists(clinical_file)) {
  stop("Clinical data file not found: ", clinical_file)
}

thyr_clinical_full <- data.table::fread(clinical_file)
cat("Clinical data loaded:", nrow(thyr_clinical_full), "cases x", 
    ncol(thyr_clinical_full), "variables\n")

# Rename REBC_ID to case_submitter_id for consistency
names(thyr_clinical_full)[names(thyr_clinical_full) == "REBC_ID"] <- "case_submitter_id"

# Check for columns with extreme NA rates (>90%)
na_rate <- colSums(is.na(thyr_clinical_full)) / nrow(thyr_clinical_full)
high_na_cols <- names(na_rate[na_rate > 0.9])
if (length(high_na_cols) > 0) {
  cat("Columns with >90% NA (will be kept but noted):", length(high_na_cols), "\n")
  if (length(high_na_cols) <= 10) {
    cat("  ", paste(high_na_cols, collapse=", "), "\n")
  }
}

# ============================================================================
# Load POC data
# ============================================================================
cat("\n--- Loading POC data ---\n")

poc_braf_file <- paste0(paths$raw, "thyr_poc_braf_mut.xlsx")
poc_ret_file <- paste0(paths$raw, "thyr_poc_ret_fusion.xlsx")

if (!file.exists(poc_braf_file) || !file.exists(poc_ret_file)) {
  stop("POC data files not found")
}

# Read POC data
poc_braf <- readxl::read_excel(poc_braf_file)
poc_ret <- readxl::read_excel(poc_ret_file)

cat("POC data columns:\n")
cat("  BRAF:", paste(names(poc_braf), collapse=", "), "\n")
cat("  RET:", paste(names(poc_ret), collapse=", "), "\n")

# Combine POC data
poc_combined <- rbind(poc_braf, poc_ret)
cat("Combined POC data:", nrow(poc_combined), "cases\n")

# POC distribution
poc_dist <- poc_combined %>%
  summarise(
    n_total = n(),
    poc_na = sum(is.na(POC)),
    poc_0 = sum(!is.na(POC) & POC == 0),
    poc_0_66 = sum(!is.na(POC) & POC > 0 & POC < 66.6),
    poc_66plus = sum(!is.na(POC) & POC >= 66.6)
  )
cat("POC distribution:\n")
print(poc_dist)

# ============================================================================
# Merge clinical and POC data
# ============================================================================
cat("\n--- Merging clinical and POC data ---\n")

# Merge POC into clinical data
thyr_case_master_full <- merge(
  thyr_clinical_full,
  poc_combined[, c("case_submitter_id", "POC")],
  by = "case_submitter_id",
  all.x = TRUE
)

cat("Master case table:", nrow(thyr_case_master_full), "cases\n")
cat("Variables retained:", ncol(thyr_case_master_full), "\n")

# ============================================================================
# Identify paired samples
# ============================================================================
cat("\n--- Identifying paired samples ---\n")

# Check which cases have both tumor and normal samples
paired_info <- sample_metadata %>%
  group_by(case_submitter_id) %>%
  summarise(
    n_samples = n(),
    n_tumor = sum(grepl("Tumor", sample_type)),
    n_normal = sum(grepl("Normal", sample_type)),
    has_pair = (n_tumor > 0) & (n_normal > 0),
    .groups = "drop"
  )

thyr_case_master_full <- merge(
  thyr_case_master_full,
  paired_info[, c("case_submitter_id", "has_pair")],
  by = "case_submitter_id",
  all.x = TRUE
)

# Set has_pair to FALSE for cases not in sample data
thyr_case_master_full$has_pair[is.na(thyr_case_master_full$has_pair)] <- FALSE

cat("Cases with paired samples:", sum(thyr_case_master_full$has_pair), "\n")

# ============================================================================
# Define driver groups based on multiple columns
# ============================================================================
cat("\n--- Defining driver groups ---\n")

# Create simplified driver classification
# Priority: RNA_CandidateDriverFusion > WGS_CandidateDriverFusion > Designated_Driver
thyr_case_master_full$driver <- NA_character_

# Identify BRAF cases (mutation-driven)
# CRITICAL: Must be BRAF.MutV600E specifically to avoid PTEN-PTX1 contamination
braf_mask <- (
  # WGS: either empty/NA or BRAF only (no co-mutations)
  (is.na(thyr_case_master_full$WGS_CandidateDriverMutation) | 
     thyr_case_master_full$WGS_CandidateDriverMutation == "" | 
     thyr_case_master_full$WGS_CandidateDriverMutation == "BRAF") &
    # RNA: either empty/NA or BRAF only (no co-mutations)
    (is.na(thyr_case_master_full$RNA_CandidateDriverMutation) | 
       thyr_case_master_full$RNA_CandidateDriverMutation == "" | 
       thyr_case_master_full$RNA_CandidateDriverMutation == "BRAF") &
    # Designated_Driver must be exactly BRAF.MutV600E
    thyr_case_master_full$Designated_Driver == "BRAF.MutV600E" &
    # Additional safety: exclude fusions
    !grepl("Fusion", thyr_case_master_full$Designated_DriverGroup)
)
thyr_case_master_full$driver[braf_mask] <- "BRAF"

# Identify RET cases (fusion-driven)
ret_mask <- grepl("RET", thyr_case_master_full$RNA_CandidateDriverFusion) |
  grepl("RET", thyr_case_master_full$WGS_CandidateDriverFusion) |
  grepl("RET", thyr_case_master_full$Designated_Driver)
thyr_case_master_full$driver[ret_mask] <- "RET"

cat("Driver assignment:\n")
cat("  BRAF:", sum(thyr_case_master_full$driver == "BRAF", na.rm = TRUE), "\n")
cat("  RET:", sum(thyr_case_master_full$driver == "RET", na.rm = TRUE), "\n")
cat("  Unassigned:", sum(is.na(thyr_case_master_full$driver)), "\n")

# Debug: Check if any non-BRAF.MutV600E cases are captured
if (any(braf_mask, na.rm = TRUE)) {
  braf_cases <- thyr_case_master_full[braf_mask, ]
  unique_designated <- unique(braf_cases$Designated_Driver)
  cat("\nBRAF cases - Designated_Driver values:\n")
  print(table(braf_cases$Designated_Driver))
  
  # Verify all are BRAF.MutV600E
  if (length(unique_designated) == 1 && unique_designated[1] == "BRAF.MutV600E") {
    cat("  ✓ All BRAF cases are correctly BRAF.MutV600E\n")
  } else {
    warning("Unexpected Designated_Driver values in BRAF cases!")
  }
}

# ============================================================================
# Define analysis groups
# ============================================================================
cat("\n--- Defining analysis groups ---\n")

# Initialize group assignment
thyr_case_master_full$group <- NA_character_
thyr_case_master_full$group_type <- NA_character_

# Main analysis groups (paired samples only)
# R0: RET + POC=NA + paired
idx_R0 <- thyr_case_master_full$driver == "RET" & 
  is.na(thyr_case_master_full$POC) & 
  thyr_case_master_full$has_pair
thyr_case_master_full$group[idx_R0] <- "R0"

# R1: RET + POC>=66.6 + paired
idx_R1 <- thyr_case_master_full$driver == "RET" & 
  !is.na(thyr_case_master_full$POC) & 
  thyr_case_master_full$POC >= 66.6 & 
  thyr_case_master_full$has_pair
thyr_case_master_full$group[idx_R1] <- "R1"

# B0: BRAF + POC=NA + paired
idx_B0 <- thyr_case_master_full$driver == "BRAF" & 
  is.na(thyr_case_master_full$POC) & 
  thyr_case_master_full$has_pair
thyr_case_master_full$group[idx_B0] <- "B0"

# B1: BRAF + POC>=66.6 + paired
idx_B1 <- thyr_case_master_full$driver == "BRAF" & 
  !is.na(thyr_case_master_full$POC) & 
  thyr_case_master_full$POC >= 66.6 & 
  thyr_case_master_full$has_pair
thyr_case_master_full$group[idx_B1] <- "B1"

# Mark main analysis groups
thyr_case_master_full$group_type[!is.na(thyr_case_master_full$group)] <- "main_analysis"

# Test groups (unpaired or intermediate POC)
# R_test: RET cases not in main analysis
idx_R_test <- thyr_case_master_full$driver == "RET" & is.na(thyr_case_master_full$group)
thyr_case_master_full$group[idx_R_test] <- "R_test"

# B_test: BRAF cases not in main analysis
idx_B_test <- thyr_case_master_full$driver == "BRAF" & is.na(thyr_case_master_full$group)
thyr_case_master_full$group[idx_B_test] <- "B_test"

# Mark test groups
thyr_case_master_full$group_type[thyr_case_master_full$group %in% c("R_test", "B_test")] <- "test"

# ============================================================================
# Create filtered version (grouped cases only)
# ============================================================================
cat("\n--- Creating filtered case master ---\n")

thyr_case_master <- thyr_case_master_full[!is.na(thyr_case_master_full$group), ]

cat(sprintf("Total cases: %d → Grouped cases: %d (%.1f%%)\n",
            nrow(thyr_case_master_full),
            nrow(thyr_case_master),
            nrow(thyr_case_master) / nrow(thyr_case_master_full) * 100))

# ============================================================================
# Reorganize column order for better usability
# ============================================================================
cat("\n--- Reorganizing column order ---\n")

# Define priority columns (in desired order)
priority_cols <- c(
  # Primary ID
  "case_submitter_id",
  # Analysis key variables
  "driver", "POC", "group", "group_type", "has_pair",
  # Basic demographics
  "SEX", "AGE_SURGERY", "AGE_EXPOSURE", "DOSE",
  # Driver details
  "Designated_Driver", "RNA_CandidateDriverFusion", "WGS_CandidateDriverFusion",
  "Designated_DriverGroup", "Designated_DriverPathway", "Designated_DriverType",
  "DriverGroup_CoMutFig",
  # Mutation details
  "RNA_CandidateDriverMutation", "WGS_CandidateDriverMutation"
)

# Get remaining columns
all_cols <- names(thyr_case_master_full)
remaining_cols <- setdiff(all_cols, priority_cols)

# Reorder columns for both datasets
new_col_order <- c(priority_cols[priority_cols %in% all_cols], remaining_cols)

thyr_case_master_full <- thyr_case_master_full[, new_col_order, with = FALSE]
thyr_case_master <- thyr_case_master[, new_col_order, with = FALSE]

cat("Columns reordered. Priority columns at front:\n")
cat("  ", paste(head(names(thyr_case_master), 10), collapse=", "), "...\n")

# ============================================================================
# Summary reports
# ============================================================================
cat("\n=== Group Summary ===\n")

# Main analysis groups
cat("Main analysis groups (paired samples required):\n")
for (g in c("R0", "R1", "B0", "B1")) {
  n_cases <- sum(thyr_case_master$group == g)
  if (n_cases > 0) {
    cat(sprintf("  %s: %d cases\n", g, n_cases))
  }
}

# Test groups
cat("\nTest groups:\n")
for (g in c("R_test", "B_test")) {
  n_cases <- sum(thyr_case_master$group == g)
  n_paired <- sum(thyr_case_master$group == g & thyr_case_master$has_pair)
  n_unpaired <- n_cases - n_paired
  if (n_cases > 0) {
    cat(sprintf("  %s: %d cases (%d paired, %d unpaired)\n", 
                g, n_cases, n_paired, n_unpaired))
  }
}

# RET fusion partner details
cat("\n=== RET Fusion Partner Distribution ===\n")
ret_cases <- thyr_case_master[thyr_case_master$driver == "RET", ]
if (nrow(ret_cases) > 0) {
  # RNA fusion partners (more reliable)
  rna_partners <- table(ret_cases$RNA_CandidateDriverFusion, useNA = "ifany")
  cat("RNA-based fusion partners:\n")
  print(sort(rna_partners[rna_partners > 0], decreasing = TRUE))
}

# BRAF mutation details (verification)
cat("\n=== BRAF Mutation Verification ===\n")
braf_cases <- thyr_case_master[thyr_case_master$driver == "BRAF", ]
if (nrow(braf_cases) > 0) {
  cat("BRAF cases - Designated_Driver distribution:\n")
  print(table(braf_cases$Designated_Driver, useNA = "ifany"))
  
  # Check for any non-BRAF.MutV600E
  non_v600e <- braf_cases$Designated_Driver != "BRAF.MutV600E"
  if (any(non_v600e, na.rm = TRUE)) {
    warning("Found non-BRAF.MutV600E cases in BRAF group!")
    print(unique(braf_cases$Designated_Driver[non_v600e]))
  } else {
    cat("  ✓ All BRAF cases are correctly BRAF.MutV600E\n")
  }
}

# POC distribution by group
cat("\n=== POC Distribution by Group ===\n")
poc_summary <- thyr_case_master %>%
  group_by(group) %>%
  summarise(
    n = n(),
    poc_na = sum(is.na(POC)),
    poc_0_66 = sum(!is.na(POC) & POC > 0 & POC < 66.6),
    poc_66plus = sum(!is.na(POC) & POC >= 66.6),
    .groups = "drop"
  ) %>%
  arrange(group)
print(as.data.frame(poc_summary))

# ============================================================================
# Save outputs
# ============================================================================
cat("\n--- Saving outputs ---\n")

# Save full case master (all cases with all clinical data)
saveRDS(thyr_case_master_full, paste0(paths$processed, "thyr_case_master_full.rds"))
cat("Full case master saved: thyr_case_master_full.rds\n")
cat("  Contains: All", nrow(thyr_case_master_full), "cases with", 
    ncol(thyr_case_master_full), "variables\n")

# Save filtered case master (main output - grouped cases only)
saveRDS(thyr_case_master, paste0(paths$processed, "thyr_case_master.rds"))
cat("Filtered case master saved: thyr_case_master.rds\n")
cat("  Contains:", nrow(thyr_case_master), "grouped cases\n")

# Export to CSV for review (subset of columns for readability)
key_cols <- c("case_submitter_id", "driver", "POC", "group", "group_type", "has_pair",
              "SEX", "AGE_SURGERY", "AGE_EXPOSURE", "DOSE",
              "Designated_Driver", "RNA_CandidateDriverFusion", "WGS_CandidateDriverFusion",
              "Designated_DriverGroup", "Designated_DriverPathway")
export_cols <- intersect(key_cols, names(thyr_case_master))

write.csv(thyr_case_master[, ..export_cols],
          paste0(paths$output, "thyr_case_master_key.csv"),
          row.names = FALSE)
cat("Key columns exported to: output/thyr_case_master_key.csv\n")

# Create detailed summary report
summary_report <- thyr_case_master %>%
  group_by(group, group_type) %>%
  summarise(
    n_cases = n(),
    n_paired = sum(has_pair),
    poc_na = sum(is.na(POC)),
    poc_mean = mean(POC, na.rm = TRUE),
    dose_available = sum(!is.na(DOSE)),
    age_surgery_mean = mean(AGE_SURGERY, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(group_type, group)

write.csv(summary_report,
          paste0(paths$output, "group_summary_report.csv"),
          row.names = FALSE)
cat("Summary report saved to: output/group_summary_report.csv\n")

# ============================================================================
# Final cleanup
# ============================================================================
# Keep only essential objects
rm(list = setdiff(ls(), c("paths", "thyr_case_master", "thyr_case_master_full", 
                          "thyr_se_strand2_nonzero")))
gc()

cat("\n=== Clinical Integration Complete ===\n")
cat("Main outputs:\n")
cat("  1. thyr_case_master_full: All", nrow(thyr_case_master_full), 
    "cases with complete clinical data\n")
cat("  2. thyr_case_master:", nrow(thyr_case_master), 
    "grouped cases (main working dataset)\n")
cat("\nGroups defined: R0, R1, B0, B1, R_test, B_test\n")
cat("Clinical variables preserved:", ncol(thyr_case_master), "\n")
cat("\nNext: Use thyr_case_master for downstream analysis\n")