# check_01_groups_revised.R
# Purpose: Check group classification and sample counts for manuscript (REVISED)
# Section 1: Dataset and Quality Control
# Date: 2025-01-26
# Fixed: Correct column names (group, has_outlier_*, tumor_purity, low_purity)

# Setup
source("analysis_v7/setup.R")

cat("=====================================\n")
cat("GROUP CLASSIFICATION CHECK (REVISED)\n")
cat("=====================================\n\n")

# ============================================================================
# Check for SummarizedExperiment (906 samples)
# ============================================================================

cat("--- Checking initial dataset (906 samples) ---\n")
se_file <- paste0(paths$processed, "thyr_se_raw.rds")
if (file.exists(se_file)) {
  se <- readRDS(se_file)
  cat("✓ SummarizedExperiment loaded\n")
  cat("□ Total samples in SE:", ncol(se), "\n")
  cat("□ Total genes in SE:", nrow(se), "\n")
  
  # Check sample types
  if ("sample_type" %in% names(colData(se))) {
    sample_types <- table(se$sample_type)
    cat("\n--- Sample types ---\n")
    print(sample_types)
  }
} else {
  cat("✗ SummarizedExperiment file not found\n")
}

# ============================================================================
# Load case master files
# ============================================================================

cat("\n--- Loading case master files ---\n")

# Full case master
if (file.exists(paste0(paths$processed, "thyr_case_master_full.rds"))) {
  case_master_full <- readRDS(paste0(paths$processed, "thyr_case_master_full.rds"))
  cat("✓ Full case master loaded:", nrow(case_master_full), "cases\n")
} else {
  cat("✗ Full case master not found\n")
}

# Case master with groups
if (file.exists(paste0(paths$processed, "thyr_case_master.rds"))) {
  case_master <- readRDS(paste0(paths$processed, "thyr_case_master.rds"))
  cat("✓ Grouped case master loaded:", nrow(case_master), "cases\n")
} else {
  cat("✗ Grouped case master not found\n")
}

# After Stage 1 (CDM)
if (file.exists(paste0(paths$processed, "thyr_case_master_stage1_filtered.rds"))) {
  case_master_s1 <- readRDS(paste0(paths$processed, "thyr_case_master_stage1_filtered.rds"))
  cat("✓ Stage 1 filtered loaded:", nrow(case_master_s1), "cases\n")
} else {
  cat("✗ Stage 1 filtered not found\n")
}

# After Stage 2 (Purity)
if (file.exists(paste0(paths$processed, "thyr_case_master_stage2_filtered.rds"))) {
  case_master_s2 <- readRDS(paste0(paths$processed, "thyr_case_master_stage2_filtered.rds"))
  cat("✓ Stage 2 filtered loaded:", nrow(case_master_s2), "cases\n")
} else {
  cat("✗ Stage 2 filtered not found\n")
}

# ============================================================================
# Initial Dataset Analysis
# ============================================================================

cat("\n=====================================\n")
cat("INITIAL DATASET\n")
cat("=====================================\n")

if (exists("case_master_full")) {
  cat("\n□ Total cases:", nrow(case_master_full), "\n")
  
  # Check driver mutations
  if ("driver" %in% colnames(case_master_full)) {
    cat("\n--- Driver mutations ---\n")
    driver_table <- table(case_master_full$driver, useNA = "ifany")
    print(driver_table)
    
    cat("\n□ RET fusion positive:", sum(case_master_full$driver == "RET", na.rm = TRUE), "\n")
    cat("□ BRAF V600E positive:", sum(case_master_full$driver == "BRAF", na.rm = TRUE), "\n")
    cat("□ Unknown/Other:", sum(is.na(case_master_full$driver) | 
                                  case_master_full$driver == "Unknown", na.rm = TRUE), "\n")
  }
}

# ============================================================================
# POC-based 4-group classification
# ============================================================================

cat("\n=====================================\n")
cat("POC-BASED GROUP CLASSIFICATION\n")
cat("=====================================\n")

if (exists("case_master")) {
  # Check group classification (should be R0, R1, B0, B1)
  if ("group" %in% colnames(case_master)) {
    cat("\n--- 4-group distribution ---\n")
    group_table <- table(case_master$group, useNA = "ifany")
    print(group_table)
    
    cat("\n□ R0 (RET, POC=0%):", sum(case_master$group == "R0", na.rm = TRUE), "\n")
    cat("□ R1 (RET, POC≥66.6%):", sum(case_master$group == "R1", na.rm = TRUE), "\n")
    cat("□ B0 (BRAF, POC=0%):", sum(case_master$group == "B0", na.rm = TRUE), "\n")
    cat("□ B1 (BRAF, POC≥66.6%):", sum(case_master$group == "B1", na.rm = TRUE), "\n")
    cat("□ Unassigned:", sum(is.na(case_master$group)), "\n")
  }
  
  # Check POC values
  if ("POC" %in% colnames(case_master)) {
    cat("\n--- POC distribution ---\n")
    poc_summary <- summary(case_master$POC)
    print(poc_summary)
    cat("□ POC = 0:", sum(case_master$POC == 0, na.rm = TRUE), "\n")
    cat("□ POC ≥ 66.6:", sum(case_master$POC >= 66.6, na.rm = TRUE), "\n")
    cat("□ POC NA:", sum(is.na(case_master$POC)), "\n")
  }
}

# ============================================================================
# Stage 1 Filtering (CDM Outliers)
# ============================================================================

cat("\n=====================================\n")
cat("STAGE 1 FILTERING (CDM OUTLIERS)\n")
cat("=====================================\n")

if (exists("case_master_s1")) {
  # Check outlier flags (has_outlier_tumor, has_outlier_normal)
  if ("has_outlier_tumor" %in% colnames(case_master_s1)) {
    cat("\n--- Tumor outliers ---\n")
    tumor_outlier_table <- table(case_master_s1$has_outlier_tumor, useNA = "ifany")
    print(tumor_outlier_table)
    cat("□ Cases with tumor outliers:", sum(case_master_s1$has_outlier_tumor == TRUE, na.rm = TRUE), "\n")
  }
  
  if ("has_outlier_normal" %in% colnames(case_master_s1)) {
    cat("\n--- Normal outliers ---\n")
    normal_outlier_table <- table(case_master_s1$has_outlier_normal, useNA = "ifany")
    print(normal_outlier_table)
    cat("□ Cases with normal outliers:", sum(case_master_s1$has_outlier_normal == TRUE, na.rm = TRUE), "\n")
  }
  
  # Cases without outliers
  no_outliers <- case_master_s1[!case_master_s1$has_outlier_tumor & !case_master_s1$has_outlier_normal, ]
  cat("\n□ Cases without any outliers:", nrow(no_outliers), "\n")
  
  if ("group" %in% colnames(no_outliers)) {
    cat("\n--- Groups after Stage 1 ---\n")
    s1_groups <- table(no_outliers$group)
    print(s1_groups)
  }
  
  # Compare before and after Stage 1
  if (exists("case_master") && "group" %in% colnames(case_master)) {
    cat("\n--- Stage 1 exclusions by group ---\n")
    for (grp in c("R0", "R1", "B0", "B1")) {
      before <- sum(case_master$group == grp, na.rm = TRUE)
      after <- sum(no_outliers$group == grp, na.rm = TRUE)
      excluded <- before - after
      cat(sprintf("□ %s: %d → %d (excluded: %d)\n", grp, before, after, excluded))
    }
  }
}

# ============================================================================
# Stage 2 Filtering (Tumor Purity)
# ============================================================================

cat("\n=====================================\n")
cat("STAGE 2 FILTERING (TUMOR PURITY)\n")
cat("=====================================\n")

if (exists("case_master_s2")) {
  # Check tumor purity information
  if ("tumor_purity" %in% colnames(case_master_s2)) {
    cat("\n--- Tumor purity statistics ---\n")
    purity_summary <- summary(case_master_s2$tumor_purity)
    print(purity_summary)
  }
  
  if ("low_purity" %in% colnames(case_master_s2)) {
    cat("\n--- Purity filtering ---\n")
    low_purity_table <- table(case_master_s2$low_purity, useNA = "ifany")
    print(low_purity_table)
    
    cat("\n□ Low purity cases (<60%):", sum(case_master_s2$low_purity == TRUE, na.rm = TRUE), "\n")
    cat("□ High purity cases (≥60%):", sum(case_master_s2$low_purity == FALSE, na.rm = TRUE), "\n")
    
    # By group
    if ("group" %in% colnames(case_master_s2)) {
      cat("\n--- Low purity by group ---\n")
      purity_by_group <- table(case_master_s2$group, 
                               case_master_s2$low_purity,
                               useNA = "ifany")
      print(purity_by_group)
    }
  }
  
  # Final high-purity, no-outlier cases
  final_cases <- case_master_s2
  
  # Apply both filters
  if ("has_outlier_tumor" %in% colnames(final_cases) && 
      "has_outlier_normal" %in% colnames(final_cases) &&
      "low_purity" %in% colnames(final_cases)) {
    
    final_cases <- final_cases[!final_cases$has_outlier_tumor & 
                                 !final_cases$has_outlier_normal & 
                                 !final_cases$low_purity, ]
    
    cat("\n--- Final high-purity, no-outlier cases ---\n")
    cat("□ Total cases remaining:", nrow(final_cases), "\n")
    
    if ("group" %in% colnames(final_cases)) {
      final_groups <- table(final_cases$group)
      print(final_groups)
      
      cat("\n--- Final counts by group ---\n")
      for (grp in c("R0", "R1", "B0", "B1")) {
        count <- sum(final_cases$group == grp, na.rm = TRUE)
        cat(sprintf("□ %s: %d cases\n", grp, count))
      }
    }
  }
  
  # Stage 2 exclusions
  if (exists("no_outliers") && "group" %in% colnames(no_outliers)) {
    cat("\n--- Stage 2 exclusions by group ---\n")
    for (grp in c("R0", "R1", "B0", "B1")) {
      before <- sum(no_outliers$group == grp, na.rm = TRUE)
      after <- sum(final_cases$group == grp, na.rm = TRUE)
      excluded <- before - after
      cat(sprintf("□ %s: %d → %d (excluded: %d)\n", grp, before, after, excluded))
    }
  }
}

# ============================================================================
# Paired Samples Analysis
# ============================================================================

cat("\n=====================================\n")
cat("PAIRED SAMPLES (TUMOR + NORMAL)\n")
cat("=====================================\n")

# Load sample lists
if (file.exists(paste0(paths$processed, "analysis_sample_lists.rds"))) {
  sample_lists <- readRDS(paste0(paths$processed, "analysis_sample_lists.rds"))
  
  cat("\n--- Paired sample counts ---\n")
  
  if ("R0_vs_R1" %in% names(sample_lists)) {
    r0r1 <- sample_lists$R0_vs_R1
    cat("□ R0 tumor samples:", length(r0r1$R0_tumor), "\n")
    cat("□ R0 normal samples:", length(r0r1$R0_normal), "\n")
    cat("□ R1 tumor samples:", length(r0r1$R1_tumor), "\n")
    cat("□ R1 normal samples:", length(r0r1$R1_normal), "\n")
  }
  
  if ("B0_vs_B1" %in% names(sample_lists)) {
    b0b1 <- sample_lists$B0_vs_B1
    cat("□ B0 tumor samples:", length(b0b1$B0_tumor), "\n")
    cat("□ B0 normal samples:", length(b0b1$B0_normal), "\n")
    cat("□ B1 tumor samples:", length(b0b1$B1_tumor), "\n")
    cat("□ B1 normal samples:", length(b0b1$B1_normal), "\n")
  }
}

# ============================================================================
# SUMMARY FOR MANUSCRIPT
# ============================================================================

cat("\n=====================================\n")
cat("SUMMARY FOR MANUSCRIPT\n")
cat("=====================================\n")

cat("\nKey numbers for Section 1:\n")
cat("----------------------------------\n")

# Initial data
if (exists("se")) {
  cat("Initial samples in SE:", ncol(se), "\n")
  cat("Initial genes in SE:", nrow(se), "\n")
}

if (exists("case_master_full")) {
  cat("Initial cases:", nrow(case_master_full), "\n")
}

# POC groups
if (exists("case_master") && "group" %in% colnames(case_master)) {
  initial_groups <- table(case_master$group)
  cat("\nInitial POC groups:\n")
  for (grp in c("R0", "R1", "B0", "B1")) {
    if (grp %in% names(initial_groups)) {
      cat(sprintf("  %s: %d\n", grp, initial_groups[grp]))
    }
  }
}

# After filtering
if (exists("final_cases") && "group" %in% colnames(final_cases)) {
  final_table <- table(final_cases$group)
  cat("\nFinal cases after QC:\n")
  for (grp in c("R0", "R1", "B0", "B1")) {
    if (grp %in% names(final_table)) {
      cat(sprintf("  %s: %d\n", grp, final_table[grp]))
    }
  }
}

cat("\n=== END OF GROUP CHECK (REVISED) ===\n")