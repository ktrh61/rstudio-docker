#!/usr/bin/env Rscript
# =============================================================================
# 14_reo_validation_driver.R
# REO Panel Validation - Phase 9: Driver Generalization (Exploratory)
# 
# Purpose:
#   Assess REO panel behavior in BRAF V600E background (B0/B1)
#   This is EXPLORATORY - generalization is NOT expected
#
# Scientific Rationale:
#   - RET fusion: strongly associated with radiation-induced cancer (60-84%)
#   - BRAF V600E: associated with sporadic thyroid cancer
#   - Panel was constructed on RET background (R0/R1)
#   - Testing on BRAF provides driver-specificity evidence
#
# Input:
#   - thyr_case_master_stage2_filtered.rds (case metadata with QC flags)
#   - reo_final_panel.rds (panel definition and threshold)
#   - thyr_se_strand2_nonzero.rds (expression data)
#
# Output:
#   - reo_validation_driver_results.rds (all results)
#   - reo_validation_driver_summary.csv (group-level summary)
#   - reo_validation_driver_samples.csv (sample-level results)
#
# Reference: REO_Panel_Implementation_Progress_v4.md
# Date: 2025-12-08
#
# IMPORTANT: TPM calculation must match training (11_reo_pair_selection.R)
#            Using log2(TPM) without +1 offset (data is zero-free)
# =============================================================================

source("analysis_v7/setup.R")

cat("\n=============================================================================\n")
cat("14_reo_validation_driver.R - REO Panel Driver Generalization (Phase 9)\n")
cat("=============================================================================\n\n")

cat("=== EXPLORATORY ANALYSIS ===\n")
cat("NOTE: BRAF V600E is associated with sporadic cancer.\n")
cat("      Generalization of RET-derived panel is NOT expected.\n")
cat("      This analysis assesses driver-specificity of the panel.\n\n")

# Load additional packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(dplyr)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
cat("=== Loading Data ===\n")

# Load case master with QC flags
case_master_path <- file.path(paths$processed, "thyr_case_master_stage2_filtered.rds")
if (!file.exists(case_master_path)) {
  stop("thyr_case_master_stage2_filtered.rds not found. Please run 06_purity_analysis.R first.")
}
case_master <- readRDS(case_master_path)
cat(sprintf("Case master loaded: %d cases\n", nrow(case_master)))

# Load REO panel
panel_path <- file.path(paths$processed, "reo_final_panel.rds")
if (!file.exists(panel_path)) {
  stop("reo_final_panel.rds not found. Please run 12_reo_panel_finalize.R first.")
}
panel_data <- readRDS(panel_path)
selected_pairs <- panel_data$selected_pairs
threshold_T <- panel_data$threshold_T
dead_zone_threshold <- panel_data$params$dead_zone_threshold

cat(sprintf("Panel loaded: %d pairs, threshold T=%d\n", nrow(selected_pairs), threshold_T))
cat(sprintf("Dead zone threshold: |r| < %.3f\n", dead_zone_threshold))

# Load SummarizedExperiment
se_path <- file.path(paths$processed, "thyr_se_strand2_nonzero.rds")
if (!file.exists(se_path)) {
  stop("thyr_se_strand2_nonzero.rds not found. Please run 03_prepare_counts.R first.")
}
se <- readRDS(se_path)
sample_metadata <- as.data.frame(colData(se))
cat(sprintf("SE loaded: %d genes x %d samples\n", nrow(se), ncol(se)))

cat("\n")

# -----------------------------------------------------------------------------
# Compute log2(TPM) matrix - MUST MATCH TRAINING
# -----------------------------------------------------------------------------
cat("=== Preparing log2(TPM) Matrix ===\n")
cat("NOTE: Using log2(TPM) without +1 offset to match training data\n")

# Get gene lengths
gene_info <- as.data.frame(rowData(se))
gene_lengths <- gene_info$gene_length
names(gene_lengths) <- rownames(se)

# Compute TPM (same method as 11_reo_pair_selection.R)
counts <- assay(se)
calc_tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
tpm_matrix <- apply(counts, 2, calc_tpm, lengths = gene_lengths)
rownames(tpm_matrix) <- rownames(counts)

# Log2 transform WITHOUT +1 offset (data is zero-free from thyr_se_strand2_nonzero.rds)
log2_tpm <- log2(tpm_matrix)

cat(sprintf("log2(TPM) matrix: %d genes x %d samples\n", nrow(log2_tpm), ncol(log2_tpm)))
cat("\n")

# -----------------------------------------------------------------------------
# Define BRAF Sample Categories
# -----------------------------------------------------------------------------
cat("=== Defining BRAF Sample Categories ===\n")

# Filter to BRAF groups only
braf_groups <- c("B0", "B1")
braf_cases <- case_master %>%
  filter(group %in% braf_groups)

cat(sprintf("BRAF cases total: %d\n", nrow(braf_cases)))

# Define QC status categories (same logic as Phase 8)
# QC_clear: has_outlier_tumor=0 AND has_outlier_normal=0 AND low_purity=0
# Outlier: has_outlier_tumor=1 OR has_outlier_normal=1
# Low_purity: low_purity=1 (but not outlier)
# Unpaired: is.na(has_outlier_tumor) - QC check was not possible

braf_cases <- braf_cases %>%
  mutate(
    qc_status = case_when(
      is.na(has_outlier_tumor) ~ "Unpaired",
      has_outlier_tumor == 1 | has_outlier_normal == 1 ~ "Outlier",
      low_purity == 1 ~ "Low_purity",
      has_outlier_tumor == 0 & has_outlier_normal == 0 & low_purity == 0 ~ "QC_clear",
      TRUE ~ "Unknown"
    )
  )

# Summary by group and QC status
cat("\n--- Sample Distribution ---\n")
qc_summary <- braf_cases %>%
  group_by(group, qc_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = qc_status, values_from = n, values_fill = 0)

print(as.data.frame(qc_summary))

# -----------------------------------------------------------------------------
# Helper Function: Compute Reversal Count (same as Phase 8)
# -----------------------------------------------------------------------------
compute_reversal_for_sample <- function(sample_id, selected_pairs, log2_tpm, 
                                        r0_majority_signs, dead_zone) {
  n_pairs <- nrow(selected_pairs)
  reversal_count <- 0
  reversal_details <- rep(NA, n_pairs)
  names(reversal_details) <- selected_pairs$pair_id
  
  for (k in 1:n_pairs) {
    up_g <- selected_pairs$up_gene[k]
    down_g <- selected_pairs$down_gene[k]
    pair_id <- selected_pairs$pair_id[k]
    
    # Check if sample exists
    if (!(sample_id %in% colnames(log2_tpm))) {
      next
    }
    
    # Compute r value
    r_val <- log2_tpm[up_g, sample_id] - log2_tpm[down_g, sample_id]
    
    # Handle infinite values
    if (is.infinite(r_val) || is.na(r_val)) {
      reversal_details[pair_id] <- NA
      next
    }
    
    # Dead zone check
    if (abs(r_val) < dead_zone) {
      reversal_details[pair_id] <- 0  # Non-reversal (conservative)
      next
    }
    
    # Check reversal against R0 majority sign
    majority_sign <- r0_majority_signs[pair_id]
    is_reversed <- sign(r_val) != majority_sign
    
    reversal_details[pair_id] <- as.integer(is_reversed)
    if (is_reversed) reversal_count <- reversal_count + 1
  }
  
  return(list(
    count = reversal_count,
    details = reversal_details
  ))
}

# -----------------------------------------------------------------------------
# Get R0 Majority Signs from Panel Data
# -----------------------------------------------------------------------------
cat("\n=== R0 Majority Signs (from RET panel construction) ===\n")

# Use the R0 samples from panel construction
r0_samples_panel <- panel_data$r0_samples
panel_log2_tpm <- panel_data$log2_tpm

# Compute majority signs for each pair
r0_majority_signs <- sapply(1:nrow(selected_pairs), function(k) {
  up_g <- selected_pairs$up_gene[k]
  down_g <- selected_pairs$down_gene[k]
  r_r0 <- panel_log2_tpm[up_g, r0_samples_panel] - panel_log2_tpm[down_g, r0_samples_panel]
  ifelse(sum(r_r0 > 0) > sum(r_r0 < 0), 1, -1)
})
names(r0_majority_signs) <- selected_pairs$pair_id

cat("R0 majority signs (reference from RET background):\n")
print(r0_majority_signs)
cat("\n")

# -----------------------------------------------------------------------------
# Process All BRAF Samples
# -----------------------------------------------------------------------------
cat("=== Processing BRAF Samples ===\n")

# Get tumor sample IDs for each case
get_tumor_sample <- function(case_id, sample_meta) {
  sample_meta$sample_submitter_id[
    sample_meta$case_submitter_id == case_id &
      sample_meta$sample_type == "Primary Tumor" &
      grepl("_merged", sample_meta$sample_submitter_id)
  ][1]
}

# Process each case
results_list <- list()

for (i in 1:nrow(braf_cases)) {
  case_id <- braf_cases$case_submitter_id[i]
  group <- braf_cases$group[i]
  qc_status <- braf_cases$qc_status[i]
  poc <- braf_cases$POC[i]
  purity <- braf_cases$tumor_purity[i]
  
  # Get tumor sample
  tumor_sample <- get_tumor_sample(case_id, sample_metadata)
  
  if (is.na(tumor_sample)) {
    # No tumor sample available
    results_list[[i]] <- data.frame(
      case_submitter_id = case_id,
      sample_id = NA_character_,
      group = group,
      qc_status = qc_status,
      POC = poc,
      tumor_purity = purity,
      reversal_count = NA_integer_,
      panel_result = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }
  
  # Compute reversal
  rev_result <- compute_reversal_for_sample(
    tumor_sample, selected_pairs, log2_tpm,
    r0_majority_signs, dead_zone_threshold
  )
  
  # Panel result
  panel_result <- ifelse(rev_result$count >= threshold_T, "Panel(+)", "Panel(-)")
  
  results_list[[i]] <- data.frame(
    case_submitter_id = case_id,
    sample_id = tumor_sample,
    group = group,
    qc_status = qc_status,
    POC = poc,
    tumor_purity = purity,
    reversal_count = rev_result$count,
    panel_result = panel_result,
    stringsAsFactors = FALSE
  )
}

results_df <- bind_rows(results_list)

cat(sprintf("Processed: %d cases\n", nrow(results_df)))
cat(sprintf("  With reversal data: %d\n", sum(!is.na(results_df$reversal_count))))
cat(sprintf("  Missing tumor sample: %d\n", sum(is.na(results_df$sample_id))))
cat("\n")

# -----------------------------------------------------------------------------
# Section 9a: QC-Cleared B0/B1 Results
# -----------------------------------------------------------------------------
cat("=== Section 9a: QC-Cleared B0/B1 Results ===\n")

qc_clear_braf <- results_df %>%
  filter(qc_status == "QC_clear", !is.na(reversal_count))

cat("\n--- B0 (QC_clear) ---\n")
b0_clear <- qc_clear_braf %>% filter(group == "B0")
if (nrow(b0_clear) > 0) {
  cat(sprintf("N = %d\n", nrow(b0_clear)))
  cat(sprintf("Reversal count: range=%d-%d, mean=%.1f, median=%.0f\n",
              min(b0_clear$reversal_count),
              max(b0_clear$reversal_count),
              mean(b0_clear$reversal_count),
              median(b0_clear$reversal_count)))
  cat(sprintf("Panel(+): %d/%d (%.1f%%)\n",
              sum(b0_clear$panel_result == "Panel(+)"),
              nrow(b0_clear),
              sum(b0_clear$panel_result == "Panel(+)") / nrow(b0_clear) * 100))
  cat("Distribution:\n")
  print(table(b0_clear$reversal_count, useNA = "ifany"))
} else {
  cat("No QC-cleared samples\n")
}

cat("\n--- B1 (QC_clear) ---\n")
b1_clear <- qc_clear_braf %>% filter(group == "B1")
if (nrow(b1_clear) > 0) {
  cat(sprintf("N = %d\n", nrow(b1_clear)))
  cat(sprintf("Reversal count: range=%d-%d, mean=%.1f, median=%.0f\n",
              min(b1_clear$reversal_count),
              max(b1_clear$reversal_count),
              mean(b1_clear$reversal_count),
              median(b1_clear$reversal_count)))
  cat(sprintf("Panel(+): %d/%d (%.1f%%)\n",
              sum(b1_clear$panel_result == "Panel(+)"),
              nrow(b1_clear),
              sum(b1_clear$panel_result == "Panel(+)") / nrow(b1_clear) * 100))
  cat("Distribution:\n")
  print(table(b1_clear$reversal_count, useNA = "ifany"))
} else {
  cat("No QC-cleared samples\n")
}

# -----------------------------------------------------------------------------
# Section 9b: QC Non-Cleared BRAF Samples
# -----------------------------------------------------------------------------
cat("\n=== Section 9b: QC Non-Cleared BRAF Samples ===\n")

non_clear_braf <- results_df %>%
  filter(qc_status != "QC_clear")

for (g in braf_groups) {
  cat(sprintf("\n--- %s (Non-QC-clear) ---\n", g))
  
  group_non_clear <- non_clear_braf %>% filter(group == g)
  
  if (nrow(group_non_clear) == 0) {
    cat("No non-QC-clear samples\n")
    next
  }
  
  # By QC status
  for (qc_stat in c("Outlier", "Low_purity", "Unpaired")) {
    subset <- group_non_clear %>% filter(qc_status == qc_stat)
    if (nrow(subset) == 0) next
    
    valid_subset <- subset %>% filter(!is.na(reversal_count))
    
    cat(sprintf("\n  [%s] N=%d (valid=%d)\n", qc_stat, nrow(subset), nrow(valid_subset)))
    
    if (nrow(valid_subset) > 0) {
      cat(sprintf("  Reversal: range=%d-%d, mean=%.1f\n",
                  min(valid_subset$reversal_count),
                  max(valid_subset$reversal_count),
                  mean(valid_subset$reversal_count)))
      cat(sprintf("  Panel(+): %d/%d (%.1f%%)\n",
                  sum(valid_subset$panel_result == "Panel(+)"),
                  nrow(valid_subset),
                  sum(valid_subset$panel_result == "Panel(+)") / nrow(valid_subset) * 100))
    }
  }
}

# -----------------------------------------------------------------------------
# Summary Statistics by Group and QC Status
# -----------------------------------------------------------------------------
cat("\n=== Summary by Group and QC Status ===\n")

summary_stats <- results_df %>%
  filter(!is.na(reversal_count)) %>%
  group_by(group, qc_status) %>%
  summarise(
    n = n(),
    rev_min = min(reversal_count),
    rev_max = max(reversal_count),
    rev_mean = round(mean(reversal_count), 1),
    rev_median = median(reversal_count),
    panel_positive = sum(panel_result == "Panel(+)"),
    panel_positive_pct = round(sum(panel_result == "Panel(+)") / n() * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(group, qc_status)

print(as.data.frame(summary_stats))

# -----------------------------------------------------------------------------
# Comparison with RET Results (Cross-Reference)
# -----------------------------------------------------------------------------
cat("\n=== Comparison with RET Background (Cross-Reference) ===\n")

# Load Phase 8 results if available
poc_results_path <- file.path(paths$processed, "reo_validation_poc_results.rds")
if (file.exists(poc_results_path)) {
  poc_results <- readRDS(poc_results_path)
  ret_summary <- poc_results$summary_by_group_qc
  
  cat("\n--- RET Background (from Phase 8) ---\n")
  ret_qc_clear <- ret_summary %>% filter(qc_status == "QC_clear")
  print(as.data.frame(ret_qc_clear))
  
  cat("\n--- BRAF Background (current Phase 9) ---\n")
  braf_qc_clear <- summary_stats %>% filter(qc_status == "QC_clear")
  print(as.data.frame(braf_qc_clear))
  
  # Direct comparison table
  cat("\n--- Direct Comparison (QC_clear only) ---\n")
  comparison <- data.frame(
    Background = c("RET (R0)", "RET (R1)", "BRAF (B0)", "BRAF (B1)"),
    N = c(
      ret_qc_clear$n[ret_qc_clear$group == "R0"],
      ret_qc_clear$n[ret_qc_clear$group == "R1"],
      ifelse(nrow(b0_clear) > 0, nrow(b0_clear), 0),
      ifelse(nrow(b1_clear) > 0, nrow(b1_clear), 0)
    ),
    Rev_mean = c(
      ret_qc_clear$rev_mean[ret_qc_clear$group == "R0"],
      ret_qc_clear$rev_mean[ret_qc_clear$group == "R1"],
      ifelse(nrow(b0_clear) > 0, round(mean(b0_clear$reversal_count), 1), NA),
      ifelse(nrow(b1_clear) > 0, round(mean(b1_clear$reversal_count), 1), NA)
    ),
    Panel_positive_pct = c(
      ret_qc_clear$panel_positive_pct[ret_qc_clear$group == "R0"],
      ret_qc_clear$panel_positive_pct[ret_qc_clear$group == "R1"],
      ifelse(nrow(b0_clear) > 0, 
             round(sum(b0_clear$panel_result == "Panel(+)") / nrow(b0_clear) * 100, 1), NA),
      ifelse(nrow(b1_clear) > 0, 
             round(sum(b1_clear$panel_result == "Panel(+)") / nrow(b1_clear) * 100, 1), NA)
    )
  )
  print(comparison)
} else {
  cat("Phase 8 results not found. Run 13_reo_validation_poc.R first for comparison.\n")
}

# -----------------------------------------------------------------------------
# POC Relationship in BRAF (Exploratory)
# -----------------------------------------------------------------------------
cat("\n=== POC vs Reversal in BRAF (Exploratory) ===\n")

poc_braf <- results_df %>%
  filter(qc_status == "QC_clear", !is.na(POC), !is.na(reversal_count))

if (nrow(poc_braf) >= 5) {
  cor_result <- cor.test(poc_braf$POC, poc_braf$reversal_count, method = "spearman")
  cat(sprintf("Spearman correlation (POC vs reversal_count):\n"))
  cat(sprintf("  rho = %.3f, p-value = %.4f\n", cor_result$estimate, cor_result$p.value))
  
  cat("\nInterpretation:\n")
  if (cor_result$p.value < 0.05 && cor_result$estimate > 0.3) {
    cat("  Significant positive correlation detected.\n")
    cat("  This suggests some POC-related signal even in BRAF background.\n")
  } else {
    cat("  No significant correlation or weak relationship.\n")
    cat("  This is consistent with BRAF being primarily sporadic.\n")
  }
} else {
  cat("Insufficient samples with POC data for correlation analysis\n")
}

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------
cat("\n=== Saving Results ===\n")

# Full results object
output <- list(
  # Parameters
  threshold_T = threshold_T,
  dead_zone_threshold = dead_zone_threshold,
  
  # Sample-level results
  results = results_df,
  
  # Summary statistics
  summary_by_group_qc = summary_stats,
  
  # R0 majority signs used (from RET)
  r0_majority_signs = r0_majority_signs,
  
  # Metadata
  analysis_type = "exploratory",
  note = "BRAF V600E is associated with sporadic cancer. Generalization not expected.",
  timestamp = Sys.time(),
  version = "v1.1"  # Updated version to reflect TPM fix
)

saveRDS(output, file.path(paths$processed, "reo_validation_driver_results.rds"))
cat(sprintf("Saved: %s\n", file.path(paths$processed, "reo_validation_driver_results.rds")))

# Summary CSV
write.csv(summary_stats, 
          file.path(paths$output, "reo_validation_driver_summary.csv"),
          row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_validation_driver_summary.csv")))

# Sample-level results CSV
write.csv(results_df,
          file.path(paths$output, "reo_validation_driver_samples.csv"),
          row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_validation_driver_samples.csv")))

# -----------------------------------------------------------------------------
# Final Summary
# -----------------------------------------------------------------------------
cat("\n=============================================================================\n")
cat("Phase 9 Exploratory Analysis Complete\n")
cat("=============================================================================\n")

cat("\n--- Key Findings ---\n")

# B0 summary
if (nrow(b0_clear) > 0) {
  cat(sprintf("\nB0 (Non-exposed BRAF, QC_clear): N=%d\n", nrow(b0_clear)))
  cat(sprintf("  Reversal: mean=%.1f, range=%d-%d\n",
              mean(b0_clear$reversal_count),
              min(b0_clear$reversal_count),
              max(b0_clear$reversal_count)))
  cat(sprintf("  Panel(+): %.1f%%\n",
              sum(b0_clear$panel_result == "Panel(+)") / nrow(b0_clear) * 100))
}

# B1 summary
if (nrow(b1_clear) > 0) {
  cat(sprintf("\nB1 (Exposed BRAF, QC_clear): N=%d\n", nrow(b1_clear)))
  cat(sprintf("  Reversal: mean=%.1f, range=%d-%d\n",
              mean(b1_clear$reversal_count),
              min(b1_clear$reversal_count),
              max(b1_clear$reversal_count)))
  cat(sprintf("  Panel(+): %.1f%%\n",
              sum(b1_clear$panel_result == "Panel(+)") / nrow(b1_clear) * 100))
}

cat("\n--- Interpretation Guidelines ---\n")
cat("Expected outcome (if panel is RET-specific):\n")
cat("  - B0 and B1 should show similar, low reversal patterns\n")
cat("  - Both groups should have low Panel(+) rates\n")
cat("  - No clear separation between B0 and B1\n")
cat("\nAlternative outcome (if cross-driver signal exists):\n")
cat("  - B1 shows elevated reversals compared to B0\n")
cat("  - This would suggest radiation signature transcends driver mutation\n")

cat("\n--- Next Steps ---\n")
cat("1. Compare B0/B1 results with R0/R1 baseline\n")
cat("2. Assess whether driver-specificity supports main findings\n")
cat("3. Document exploratory results for supplementary material\n")

cat("\n=== Done ===\n")