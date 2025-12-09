#!/usr/bin/env Rscript
# =============================================================================
# 13_reo_evaluation_poc.R
# REO Panel Evaluation - Performance Assessment on R_Low/R_Mid
# 
# Purpose:
#   - Apply panel to R_Low and R_Mid (intermediate POC samples)
#   - Report performance for reference (NOT for threshold adjustment)
#   - This is an exploratory/descriptive analysis only
#
# Input:
#   - thyr_case_master_stage2_filtered.rds (case metadata with QC flags)
#   - reo_final_panel.rds (panel definition and boundary zone)
#   - thyr_se_strand2_nonzero.rds (expression data)
#
# Output:
#   - reo_evaluation_results.rds (all results)
#   - reo_evaluation_summary.csv (group-level summary)
#   - reo_evaluation_samples.csv (sample-level results)
#
# Version: v7.1
# Date: 2025-12-09
#
# IMPORTANT: 
#   - This script does NOT modify the panel or boundary zone
#   - Results are for exploratory/descriptive purposes only
#   - Design decisions were made in 11_ and 12_ using R0/R1 only
# =============================================================================

source("analysis_v7/setup.R")

cat("\n=============================================================================\n")
cat("13_reo_evaluation_poc.R - REO Panel Evaluation (v7.1)\n")
cat("Date: 2025-12-09\n")
cat("=============================================================================\n\n")

cat("*** NOTE: This is an exploratory/descriptive analysis only ***\n")
cat("*** Panel and boundary zone are NOT modified based on these results ***\n\n")

# Load additional packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
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
boundary <- panel_data$boundary
dead_zone <- panel_data$params$dead_zone
r0_majority_signs <- panel_data$r0_majority_signs

cat(sprintf("Panel loaded: %d pairs (version: %s)\n", 
            nrow(selected_pairs), panel_data$version))
cat(sprintf("Boundary: negative <= %d, positive >= %d\n", 
            boundary$negative_max, boundary$positive_min))

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
# Define Sample Categories (before TPM to identify required samples)
# -----------------------------------------------------------------------------
cat("=== Sample Categories ===\n")

# Filter to RET groups
ret_groups <- c("R0", "R_Low", "R_Mid", "R1")
ret_cases <- case_master %>%
  filter(group %in% ret_groups)

cat(sprintf("RET cases total: %d\n", nrow(ret_cases)))

# Define QC status
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

# Summary by group
cat("\n--- Sample Distribution by Group ---\n")
group_dist <- ret_cases %>%
  group_by(group) %>%
  summarise(
    total = n(),
    QC_clear = sum(qc_status == "QC_clear"),
    .groups = "drop"
  )
print(as.data.frame(group_dist))

# -----------------------------------------------------------------------------
# Helper Function: Get Tumor Sample ID
# -----------------------------------------------------------------------------
get_tumor_sample <- function(case_id, sample_meta) {
  matches <- sample_meta$sample_submitter_id[
    sample_meta$case_submitter_id == case_id &
      sample_meta$sample_type == "Primary Tumor" &
      grepl("_merged", sample_meta$sample_submitter_id)
  ]
  if (length(matches) == 0) return(NA_character_)
  return(matches[1])
}

# -----------------------------------------------------------------------------
# Compute log2(TPM) for Panel Genes x RET Samples Only
# -----------------------------------------------------------------------------
cat("\n=== Preparing log2(TPM) Matrix (filtered) ===\n")

# Get required samples (tumor samples for RET cases)
required_samples <- sapply(ret_cases$case_submitter_id, get_tumor_sample, 
                           sample_meta = sample_metadata)
required_samples <- required_samples[!is.na(required_samples)]
required_samples <- intersect(required_samples, colnames(se))

cat(sprintf("Required samples: %d\n", length(required_samples)))

# Get required genes (panel genes only)
panel_genes <- unique(c(selected_pairs$up_gene, selected_pairs$down_gene))
panel_genes <- intersect(panel_genes, rownames(se))

cat(sprintf("Required genes (panel): %d\n", length(panel_genes)))

# Extract subset of counts
counts_subset <- assay(se)[panel_genes, required_samples, drop = FALSE]
gene_lengths_subset <- rowData(se)$gene_length[match(panel_genes, rownames(se))]

calc_tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

tpm_matrix <- apply(counts_subset, 2, calc_tpm, lengths = gene_lengths_subset)
rownames(tpm_matrix) <- panel_genes

log2_tpm <- log2(tpm_matrix)

cat(sprintf("log2(TPM) matrix: %d genes x %d samples\n", nrow(log2_tpm), ncol(log2_tpm)))
cat("\n")

# -----------------------------------------------------------------------------
# Helper Function: Compute Reversal Count
# -----------------------------------------------------------------------------
compute_reversal_for_sample <- function(sample_id, selected_pairs, log2_tpm, 
                                        r0_majority_signs, dead_zone) {
  n_pairs <- nrow(selected_pairs)
  reversal_count <- 0
  
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
    
    # Handle infinite/NA values
    if (is.infinite(r_val) || is.na(r_val)) {
      next
    }
    
    # Dead zone check: if inside, treat as non-reversal
    if (abs(r_val) < dead_zone) {
      next
    }
    
    # Check reversal against R0 majority sign
    majority_sign <- r0_majority_signs[pair_id]
    if (sign(r_val) != majority_sign) {
      reversal_count <- reversal_count + 1
    }
  }
  
  return(reversal_count)
}

# -----------------------------------------------------------------------------
# Classification Function
# -----------------------------------------------------------------------------
classify_sample <- function(reversal_count, boundary) {
  if (is.na(reversal_count)) return(NA_character_)
  
  if (reversal_count <= boundary$negative_max && 
      (boundary$has_gap || reversal_count < boundary$positive_min)) {
    return("negative")
  } else if (reversal_count >= boundary$positive_min && 
             (boundary$has_gap || reversal_count > boundary$negative_max)) {
    return("positive")
  } else {
    return("undetermined")
  }
}

# -----------------------------------------------------------------------------
# Apply Panel to All Samples
# -----------------------------------------------------------------------------
cat("\n=== Applying Panel to All RET Samples ===\n")

results_list <- list()

for (i in 1:nrow(ret_cases)) {
  case_id <- ret_cases$case_submitter_id[i]
  group <- ret_cases$group[i]
  qc_status <- ret_cases$qc_status[i]
  poc <- ret_cases$POC[i]
  
  # Get tumor sample ID
  tumor_id <- get_tumor_sample(case_id, sample_metadata)
  
  # Compute reversal
  if (is.na(tumor_id) || !tumor_id %in% colnames(log2_tpm)) {
    reversal_count <- NA
    classification <- NA
  } else {
    reversal_count <- compute_reversal_for_sample(
      tumor_id, selected_pairs, log2_tpm, r0_majority_signs, dead_zone
    )
    classification <- classify_sample(reversal_count, boundary)
  }
  
  results_list[[i]] <- data.frame(
    case_id = case_id,
    tumor_id = ifelse(is.na(tumor_id), NA_character_, tumor_id),
    group = group,
    qc_status = qc_status,
    POC = poc,
    reversal_count = reversal_count,
    classification = classification,
    stringsAsFactors = FALSE
  )
}

results_df <- do.call(rbind, results_list)

cat(sprintf("Processed %d cases\n", nrow(results_df)))
cat(sprintf("Valid results: %d\n", sum(!is.na(results_df$reversal_count))))

# -----------------------------------------------------------------------------
# Summary by Group
# -----------------------------------------------------------------------------
cat("\n=== Results by Group ===\n")

summary_by_group <- results_df %>%
  filter(!is.na(reversal_count)) %>%
  group_by(group) %>%
  summarise(
    n = n(),
    rev_mean = round(mean(reversal_count), 2),
    rev_sd = round(sd(reversal_count), 2),
    rev_min = min(reversal_count),
    rev_max = max(reversal_count),
    n_negative = sum(classification == "negative"),
    n_undetermined = sum(classification == "undetermined"),
    n_positive = sum(classification == "positive"),
    pct_positive = round(sum(classification == "positive") / n() * 100, 1),
    .groups = "drop"
  )

print(as.data.frame(summary_by_group))

# QC_clear only
cat("\n=== Results by Group (QC_clear only) ===\n")

summary_qc_clear <- results_df %>%
  filter(!is.na(reversal_count), qc_status == "QC_clear") %>%
  group_by(group) %>%
  summarise(
    n = n(),
    rev_mean = round(mean(reversal_count), 2),
    rev_sd = round(sd(reversal_count), 2),
    rev_min = min(reversal_count),
    rev_max = max(reversal_count),
    n_negative = sum(classification == "negative"),
    n_undetermined = sum(classification == "undetermined"),
    n_positive = sum(classification == "positive"),
    pct_positive = round(sum(classification == "positive") / n() * 100, 1),
    .groups = "drop"
  )

print(as.data.frame(summary_qc_clear))

# -----------------------------------------------------------------------------
# Detailed View: Intermediate Groups (R_Low, R_Mid)
# -----------------------------------------------------------------------------
cat("\n=== Detailed: Intermediate POC Groups ===\n")

intermediate_samples <- results_df %>%
  filter(group %in% c("R_Low", "R_Mid"), !is.na(reversal_count)) %>%
  arrange(group, reversal_count)

if (nrow(intermediate_samples) > 0) {
  cat("\nR_Low samples:\n")
  r_low <- intermediate_samples %>% filter(group == "R_Low")
  if (nrow(r_low) > 0) {
    print(r_low[, c("tumor_id", "POC", "qc_status", "reversal_count", "classification")])
  } else {
    cat("  No R_Low samples\n")
  }
  
  cat("\nR_Mid samples:\n")
  r_mid <- intermediate_samples %>% filter(group == "R_Mid")
  if (nrow(r_mid) > 0) {
    print(r_mid[, c("tumor_id", "POC", "qc_status", "reversal_count", "classification")])
  } else {
    cat("  No R_Mid samples\n")
  }
}

# -----------------------------------------------------------------------------
# POC vs Reversal Correlation
# -----------------------------------------------------------------------------
cat("\n=== POC vs Reversal Analysis ===\n")

poc_analysis <- results_df %>%
  filter(!is.na(reversal_count), !is.na(POC))

if (nrow(poc_analysis) > 5) {
  cor_result <- cor.test(poc_analysis$POC, poc_analysis$reversal_count, 
                         method = "spearman")
  
  cat(sprintf("Spearman correlation (POC vs reversal count):\n"))
  cat(sprintf("  rho = %.3f\n", cor_result$estimate))
  cat(sprintf("  p-value = %.4f\n", cor_result$p.value))
  
  if (cor_result$estimate > 0 && cor_result$p.value < 0.05) {
    cat("  -> Positive correlation as expected (higher POC = more reversals)\n")
  } else if (cor_result$p.value >= 0.05) {
    cat("  -> No significant correlation detected\n")
  }
}

# -----------------------------------------------------------------------------
# Comparison with Training Data (R0/R1 from 12_)
# -----------------------------------------------------------------------------
cat("\n=== Comparison with Panel Construction ===\n")

cat("\nFrom panel construction (12_):\n")
cat(sprintf("  R0: n=%d, range=%d-%d, mean=%.1f\n",
            length(panel_data$r0_reversals),
            min(panel_data$r0_reversals),
            max(panel_data$r0_reversals),
            mean(panel_data$r0_reversals)))
cat(sprintf("  R1: n=%d, range=%d-%d, mean=%.1f\n",
            length(panel_data$r1_reversals),
            min(panel_data$r1_reversals),
            max(panel_data$r1_reversals),
            mean(panel_data$r1_reversals)))

cat("\nFrom this evaluation:\n")
for (g in c("R0", "R_Low", "R_Mid", "R1")) {
  g_data <- results_df %>% 
    filter(group == g, !is.na(reversal_count))
  if (nrow(g_data) > 0) {
    cat(sprintf("  %s: n=%d, range=%d-%d, mean=%.1f, positive=%.0f%%\n",
                g, nrow(g_data),
                min(g_data$reversal_count),
                max(g_data$reversal_count),
                mean(g_data$reversal_count),
                sum(g_data$classification == "positive") / nrow(g_data) * 100))
  }
}

# -----------------------------------------------------------------------------
# Expected Pattern Check
# -----------------------------------------------------------------------------
cat("\n=== Expected Pattern Check ===\n")

# Expected: R0 (low) < R_Low < R_Mid < R1 (high)
group_means <- results_df %>%
  filter(!is.na(reversal_count)) %>%
  group_by(group) %>%
  summarise(mean_rev = mean(reversal_count), .groups = "drop") %>%
  mutate(group = factor(group, levels = c("R0", "R_Low", "R_Mid", "R1"))) %>%
  arrange(group)

cat("\nMean reversal count by group (expected: increasing order):\n")
print(as.data.frame(group_means))

# Check monotonicity
if (nrow(group_means) >= 2) {
  means <- group_means$mean_rev
  is_monotonic <- all(diff(means) >= 0)
  if (is_monotonic) {
    cat("\n✓ Monotonically increasing (as expected)\n")
  } else {
    cat("\n⚠ Not strictly monotonic (may warrant investigation)\n")
  }
}

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------
cat("\n=== Saving Results ===\n")

output <- list(
  # Version info
  version = "v7.1",
  date = "2025-12-09",
  
  # Panel info used
  panel_version = panel_data$version,
  boundary = boundary,
  dead_zone = dead_zone,
  n_pairs = nrow(selected_pairs),
  
  # Results
  results = results_df,
  summary_by_group = summary_by_group,
  summary_qc_clear = summary_qc_clear,
  
  # Metadata
  timestamp = Sys.time()
)

saveRDS(output, file.path(paths$processed, "reo_evaluation_results.rds"))
cat(sprintf("Saved: %s\n", file.path(paths$processed, "reo_evaluation_results.rds")))

# Summary CSV
write.csv(summary_by_group, 
          file.path(paths$output, "reo_evaluation_summary.csv"),
          row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_evaluation_summary.csv")))

# Sample-level CSV
write.csv(results_df,
          file.path(paths$output, "reo_evaluation_samples.csv"),
          row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_evaluation_samples.csv")))

# -----------------------------------------------------------------------------
# Final Summary
# -----------------------------------------------------------------------------
cat("\n=============================================================================\n")
cat("Evaluation Complete\n")
cat("=============================================================================\n")

cat("\n*** REMINDER: This is exploratory/descriptive only ***\n")
cat("*** Panel and boundary zone remain unchanged ***\n")

cat(sprintf("\nPanel: %d pairs\n", nrow(selected_pairs)))
cat(sprintf("Boundary: negative <= %d, positive >= %d\n", 
            boundary$negative_max, boundary$positive_min))

cat("\nKey findings:\n")

# R_Low/R_Mid summary
for (g in c("R_Low", "R_Mid")) {
  g_data <- summary_by_group %>% filter(group == g)
  if (nrow(g_data) > 0) {
    cat(sprintf("  %s (n=%d): %.0f%% positive, mean reversal=%.1f\n",
                g, g_data$n, g_data$pct_positive, g_data$rev_mean))
  }
}

cat("\n=== Done ===\n")