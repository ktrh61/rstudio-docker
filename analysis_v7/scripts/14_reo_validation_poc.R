#!/usr/bin/env Rscript
# =============================================================================
# 14_reo_validation_poc.R
# REO Panel Validation - Phase 8: POC Intermediate Value Validation
# 
# Purpose:
#   - 8a: Evaluate panel performance on R_Low/R_High (QC-cleared samples)
#   - 8b: Evaluate QC impact using non-cleared samples across all RET groups
#
# Input:
#   - thyr_case_master_stage2_filtered.rds (case metadata with QC flags)
#   - reo_final_panel.rds (panel definition and threshold)
#   - thyr_se_strand2_nonzero.rds (expression data)
#
# Output:
#   - reo_validation_poc_results.rds (all results)
#   - reo_validation_poc_summary.csv (group-level summary)
#   - reo_validation_poc_qc_impact.csv (QC impact assessment)
#   - reo_validation_poc_samples.csv (sample-level results)
#
# Reference: REO_Panel_Implementation_Progress_v4.md
# Date: 2025-12-07
# =============================================================================

source("analysis_v7/setup.R")

cat("\n=============================================================================\n")
cat("14_reo_validation_poc.R - REO Panel POC Validation (Phase 8)\n")
cat("=============================================================================\n\n")

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
  stop("reo_final_panel.rds not found. Please run 13_reo_panel_finalize.R first.")
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
cat(sprintf("SE loaded: %d genes × %d samples\n", nrow(se), ncol(se)))

cat("\n")

# -----------------------------------------------------------------------------
# Compute log2(TPM) matrix
# -----------------------------------------------------------------------------
cat("=== Preparing log2(TPM) Matrix ===\n")

# Get gene lengths
gene_info <- as.data.frame(rowData(se))
gene_lengths <- gene_info$gene_length
names(gene_lengths) <- rownames(se)

# Compute TPM
counts <- assay(se)
rpk <- sweep(counts, 1, gene_lengths / 1000, "/")
tpm <- sweep(rpk, 2, colSums(rpk) / 1e6, "/")
log2_tpm <- log2(tpm + 1)

cat(sprintf("log2(TPM+1) matrix: %d genes × %d samples\n", nrow(log2_tpm), ncol(log2_tpm)))
cat("\n")

# -----------------------------------------------------------------------------
# Define Sample Categories
# -----------------------------------------------------------------------------
cat("=== Defining Sample Categories ===\n")

# Filter to RET groups only
ret_groups <- c("R0", "R_Low", "R_High", "R1")
ret_cases <- case_master %>%
  filter(group %in% ret_groups)

cat(sprintf("RET cases total: %d\n", nrow(ret_cases)))

# Define QC status categories
# QC_clear: has_outlier_tumor=0 AND has_outlier_normal=0 AND low_purity=0
# Outlier: has_outlier_tumor=1 OR has_outlier_normal=1
# Low_purity: low_purity=1 (but not outlier)
# Unpaired: is.na(has_outlier_tumor) - QC check was not possible

ret_cases <- ret_cases %>%
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
qc_summary <- ret_cases %>%
  group_by(group, qc_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = qc_status, values_from = n, values_fill = 0)

print(as.data.frame(qc_summary))

# -----------------------------------------------------------------------------
# Helper Function: Compute Reversal Count
# -----------------------------------------------------------------------------
compute_reversal_for_sample <- function(sample_id, selected_pairs, log2_tpm, 
                                        r0_majority_signs, dead_zone) {
  # sample_id: sample_submitter_id of the tumor sample
  # r0_majority_signs: named vector of majority signs for each pair (from panel construction)
  
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
cat("\n=== Getting R0 Majority Signs ===\n")

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

cat("R0 majority signs (from panel construction):\n")
print(r0_majority_signs)
cat("\n")

# -----------------------------------------------------------------------------
# Process All Samples
# -----------------------------------------------------------------------------
cat("=== Processing Samples ===\n")

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

for (i in 1:nrow(ret_cases)) {
  case_id <- ret_cases$case_submitter_id[i]
  group <- ret_cases$group[i]
  qc_status <- ret_cases$qc_status[i]
  poc <- ret_cases$POC[i]
  purity <- ret_cases$tumor_purity[i]
  
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
# Section 8a: QC-Cleared R_Low/R_High Results
# -----------------------------------------------------------------------------
cat("=== Section 8a: QC-Cleared R_Low/R_High Results ===\n")

qc_clear_validation <- results_df %>%
  filter(group %in% c("R_Low", "R_High"), qc_status == "QC_clear")

if (nrow(qc_clear_validation) > 0) {
  cat("\n--- R_Low (QC_clear) ---\n")
  r_low_clear <- qc_clear_validation %>% filter(group == "R_Low")
  if (nrow(r_low_clear) > 0) {
    cat(sprintf("N = %d\n", nrow(r_low_clear)))
    cat(sprintf("Reversal count: range=%d-%d, mean=%.1f, median=%.0f\n",
                min(r_low_clear$reversal_count, na.rm = TRUE),
                max(r_low_clear$reversal_count, na.rm = TRUE),
                mean(r_low_clear$reversal_count, na.rm = TRUE),
                median(r_low_clear$reversal_count, na.rm = TRUE)))
    cat(sprintf("Panel(+): %d/%d (%.1f%%)\n",
                sum(r_low_clear$panel_result == "Panel(+)", na.rm = TRUE),
                nrow(r_low_clear),
                sum(r_low_clear$panel_result == "Panel(+)", na.rm = TRUE) / nrow(r_low_clear) * 100))
    cat("Distribution:\n")
    print(table(r_low_clear$reversal_count, useNA = "ifany"))
  } else {
    cat("No QC-cleared samples\n")
  }
  
  cat("\n--- R_High (QC_clear) ---\n")
  r_high_clear <- qc_clear_validation %>% filter(group == "R_High")
  if (nrow(r_high_clear) > 0) {
    cat(sprintf("N = %d\n", nrow(r_high_clear)))
    cat(sprintf("Reversal count: range=%d-%d, mean=%.1f, median=%.0f\n",
                min(r_high_clear$reversal_count, na.rm = TRUE),
                max(r_high_clear$reversal_count, na.rm = TRUE),
                mean(r_high_clear$reversal_count, na.rm = TRUE),
                median(r_high_clear$reversal_count, na.rm = TRUE)))
    cat(sprintf("Panel(+): %d/%d (%.1f%%)\n",
                sum(r_high_clear$panel_result == "Panel(+)", na.rm = TRUE),
                nrow(r_high_clear),
                sum(r_high_clear$panel_result == "Panel(+)", na.rm = TRUE) / nrow(r_high_clear) * 100))
    cat("Distribution:\n")
    print(table(r_high_clear$reversal_count, useNA = "ifany"))
  } else {
    cat("No QC-cleared samples\n")
  }
} else {
  cat("No QC-cleared samples in R_Low/R_High\n")
}

# -----------------------------------------------------------------------------
# Section 8b: QC Non-Cleared/Difficult Samples
# -----------------------------------------------------------------------------
cat("\n=== Section 8b: QC Non-Cleared/Difficult Samples ===\n")

non_clear_samples <- results_df %>%
  filter(qc_status != "QC_clear")

for (g in ret_groups) {
  cat(sprintf("\n--- %s (Non-QC-clear) ---\n", g))
  
  group_non_clear <- non_clear_samples %>% filter(group == g)
  
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
# QC Impact Assessment
# -----------------------------------------------------------------------------
cat("\n=== QC Impact Assessment ===\n")

# Compare QC_clear vs non-clear within each group
qc_impact <- results_df %>%
  filter(!is.na(reversal_count)) %>%
  mutate(qc_category = ifelse(qc_status == "QC_clear", "QC_clear", "Non_clear")) %>%
  group_by(group, qc_category) %>%
  summarise(
    n = n(),
    rev_mean = round(mean(reversal_count), 2),
    rev_sd = round(sd(reversal_count), 2),
    panel_positive_pct = round(sum(panel_result == "Panel(+)") / n() * 100, 1),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = qc_category,
    values_from = c(n, rev_mean, rev_sd, panel_positive_pct),
    names_sep = "_"
  )

print(as.data.frame(qc_impact))

# -----------------------------------------------------------------------------
# POC Monotonicity Check
# -----------------------------------------------------------------------------
cat("\n=== POC Monotonicity Check ===\n")

# For QC_clear samples only
poc_check <- results_df %>%
  filter(qc_status == "QC_clear", !is.na(POC), !is.na(reversal_count))

if (nrow(poc_check) > 0) {
  # Correlation between POC and reversal count
  cor_result <- cor.test(poc_check$POC, poc_check$reversal_count, method = "spearman")
  
  cat(sprintf("Spearman correlation (POC vs reversal_count):\n"))
  cat(sprintf("  rho = %.3f, p-value = %.4f\n", cor_result$estimate, cor_result$p.value))
  
  # Group means
  poc_group_means <- poc_check %>%
    group_by(group) %>%
    summarise(
      n = n(),
      poc_mean = mean(POC),
      rev_mean = mean(reversal_count),
      .groups = "drop"
    )
  
  cat("\nGroup means (QC_clear only):\n")
  print(as.data.frame(poc_group_means))
}

# -----------------------------------------------------------------------------
# Comparison with Panel Construction Data (R0/R1)
# -----------------------------------------------------------------------------
cat("\n=== Comparison with Panel Construction (Reference) ===\n")

# R0 from panel construction
cat("R0 (from panel construction):\n")
cat(sprintf("  N = %d\n", length(panel_data$r0_samples)))
cat(sprintf("  Reversal: range=%d-%d, mean=%.1f\n",
            min(panel_data$r0_reversals),
            max(panel_data$r0_reversals),
            mean(panel_data$r0_reversals)))
cat(sprintf("  Panel(+): %d/%d\n",
            sum(panel_data$r0_panel_positive),
            length(panel_data$r0_samples)))

# R1 from panel construction
cat("\nR1 (from panel construction):\n")
cat(sprintf("  N = %d\n", length(panel_data$r1_samples)))
cat(sprintf("  Reversal: range=%d-%d, mean=%.1f\n",
            min(panel_data$r1_reversals),
            max(panel_data$r1_reversals),
            mean(panel_data$r1_reversals)))
cat(sprintf("  Panel(+): %d/%d\n",
            sum(panel_data$r1_panel_positive),
            length(panel_data$r1_samples)))

# Current validation - R0 QC_clear
cat("\nR0 (current validation, QC_clear):\n")
r0_current <- results_df %>% filter(group == "R0", qc_status == "QC_clear", !is.na(reversal_count))
if (nrow(r0_current) > 0) {
  cat(sprintf("  N = %d\n", nrow(r0_current)))
  cat(sprintf("  Reversal: range=%d-%d, mean=%.1f\n",
              min(r0_current$reversal_count),
              max(r0_current$reversal_count),
              mean(r0_current$reversal_count)))
  cat(sprintf("  Panel(+): %d/%d\n",
              sum(r0_current$panel_result == "Panel(+)"),
              nrow(r0_current)))
}

# Current validation - R1 QC_clear
cat("\nR1 (current validation, QC_clear):\n")
r1_current <- results_df %>% filter(group == "R1", qc_status == "QC_clear", !is.na(reversal_count))
if (nrow(r1_current) > 0) {
  cat(sprintf("  N = %d\n", nrow(r1_current)))
  cat(sprintf("  Reversal: range=%d-%d, mean=%.1f\n",
              min(r1_current$reversal_count),
              max(r1_current$reversal_count),
              mean(r1_current$reversal_count)))
  cat(sprintf("  Panel(+): %d/%d\n",
              sum(r1_current$panel_result == "Panel(+)"),
              nrow(r1_current)))
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
  qc_impact = qc_impact,
  
  # R0 majority signs used
  r0_majority_signs = r0_majority_signs,
  
  # Metadata
  timestamp = Sys.time(),
  version = "v1.0"
)

saveRDS(output, file.path(paths$processed, "reo_validation_poc_results.rds"))
cat(sprintf("Saved: %s\n", file.path(paths$processed, "reo_validation_poc_results.rds")))

# Summary CSV
write.csv(summary_stats, 
          file.path(paths$output, "reo_validation_poc_summary.csv"),
          row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_validation_poc_summary.csv")))

# QC impact CSV
write.csv(qc_impact,
          file.path(paths$output, "reo_validation_poc_qc_impact.csv"),
          row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_validation_poc_qc_impact.csv")))

# Sample-level results CSV
write.csv(results_df,
          file.path(paths$output, "reo_validation_poc_samples.csv"),
          row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_validation_poc_samples.csv")))

# -----------------------------------------------------------------------------
# Final Summary
# -----------------------------------------------------------------------------
cat("\n=============================================================================\n")
cat("Phase 8 Validation Complete\n")
cat("=============================================================================\n")

cat("\n--- Key Findings ---\n")

# R_Low/R_High QC_clear summary
validation_summary <- results_df %>%
  filter(group %in% c("R_Low", "R_High"), qc_status == "QC_clear", !is.na(reversal_count))

if (nrow(validation_summary) > 0) {
  cat(sprintf("\nPOC Intermediate (R_Low + R_High, QC_clear): N=%d\n", nrow(validation_summary)))
  cat(sprintf("  Reversal range: %d-%d\n",
              min(validation_summary$reversal_count),
              max(validation_summary$reversal_count)))
  cat(sprintf("  Panel(+) rate: %.1f%%\n",
              sum(validation_summary$panel_result == "Panel(+)") / nrow(validation_summary) * 100))
}

# QC impact summary
cat("\nQC Impact:\n")
cat("  See reo_validation_poc_qc_impact.csv for detailed comparison\n")

cat("\n--- Next Steps ---\n")
cat("1. Review POC monotonicity (POC↑ → reversal↑ expected)\n")
cat("2. Decide on threshold T adjustment if needed\n")
cat("3. Proceed to Phase 9: Driver generalization (B0/B1)\n")

cat("\n=== Done ===\n")