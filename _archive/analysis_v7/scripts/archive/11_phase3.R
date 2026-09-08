# 11_phase3.R - REO Pair Evaluation and Filtering
# Purpose: Systematic evaluation of all REO pairs with proper matrix formatting
# Input: reo_phase1_data.rds, reo_phase2_data.rds
# Output: Clean REO matrices and filtered top pairs
# Version: v1.0
# Date: 2025-01-23

source("analysis_v7/setup.R")

cat("\n=== Phase 3: REO Pair Evaluation and Filtering ===\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
})

# ============================================================================
# Load Phase 1 and 2 data
# ============================================================================

cat("--- Loading previous phase data ---\n")

# Phase 1: TPM data
phase1_data <- readRDS(paste0(paths$processed, "reo_phase1_data.rds"))
cat(sprintf("  Phase 1: %d up genes, %d down genes\n", 
            length(phase1_data$up_genes), 
            length(phase1_data$down_genes)))

# Phase 2: Pair information
phase2_data <- readRDS(paste0(paths$processed, "reo_phase2_data.rds"))
cat(sprintf("  Phase 2: %s pairs\n", format(phase2_data$n_pairs, big.mark = ",")))

# ============================================================================
# Step 1: Create clean REO matrices with proper row names
# ============================================================================

cat("\n--- Step 1: Creating clean REO matrices ---\n")

# Extract pair information
pair_info <- phase2_data$pair_info

# Create proper row names (geneA/geneB format)
pair_names <- paste0(pair_info$gene_up, "/", pair_info$gene_down)
cat(sprintf("  Created %d pair names\n", length(pair_names)))

# Get REO matrices (remove incorrect row names)
reo_matrix_r0 <- phase2_data$reo_matrix_r0
reo_matrix_r1 <- phase2_data$reo_matrix_r1

# Set proper row and column names
rownames(reo_matrix_r0) <- pair_names
rownames(reo_matrix_r1) <- pair_names
colnames(reo_matrix_r0) <- phase2_data$r0_samples
colnames(reo_matrix_r1) <- phase2_data$r1_samples

cat(sprintf("  REO matrix R0: %d pairs × %d samples\n", 
            nrow(reo_matrix_r0), ncol(reo_matrix_r0)))
cat(sprintf("  REO matrix R1: %d pairs × %d samples\n", 
            nrow(reo_matrix_r1), ncol(reo_matrix_r1)))

# ============================================================================
# Step 2: Calculate consistency metrics
# ============================================================================

cat("\n--- Step 2: Calculating consistency metrics ---\n")

# R0 consistency (how many samples show the same direction)
r0_positive <- rowSums(reo_matrix_r0 > 0)
r0_negative <- rowSums(reo_matrix_r0 <= 0)
r0_consistency <- pmax(r0_positive, r0_negative)

# Expected direction in R0 (based on majority)
r0_expected_direction <- ifelse(r0_positive >= r0_negative, 1, -1)

# R1 reversal (how many samples show opposite direction to R0 expectation)
r1_reversal <- numeric(nrow(reo_matrix_r1))
for(i in 1:nrow(reo_matrix_r1)) {
  if(r0_expected_direction[i] == 1) {
    # R0 expects positive, count negatives in R1
    r1_reversal[i] <- sum(reo_matrix_r1[i,] <= 0)
  } else {
    # R0 expects negative, count positives in R1
    r1_reversal[i] <- sum(reo_matrix_r1[i,] > 0)
  }
}

# Add to pair information
pair_info$r0_consistency <- r0_consistency
pair_info$r0_expected_direction <- r0_expected_direction
pair_info$r1_reversal <- r1_reversal

cat(sprintf("  Calculated consistency for %s pairs\n", 
            format(nrow(pair_info), big.mark = ",")))

# ============================================================================
# Step 3: Apply minimum TPM filter
# ============================================================================

cat("\n--- Step 3: Applying minimum TPM filter ---\n")

# Minimum TPM threshold to avoid noise in low expression
MIN_TPM_THRESHOLD <- 1.0

# Get TPM matrices from Phase 1 data
tpm_r0 <- phase1_data$tpm_r0
tpm_r1 <- phase1_data$tpm_r1

cat(sprintf("  Checking minimum TPM threshold: %.1f\n", MIN_TPM_THRESHOLD))

# Vectorized approach for efficiency
# Check which genes meet the minimum TPM in all samples
genes_pass_r0 <- rownames(tpm_r0)[rowSums(tpm_r0 >= MIN_TPM_THRESHOLD) == ncol(tpm_r0)]
genes_pass_r1 <- rownames(tpm_r1)[rowSums(tpm_r1 >= MIN_TPM_THRESHOLD) == ncol(tpm_r1)]

# Genes must pass in both R0 and R1
genes_pass_both <- intersect(genes_pass_r0, genes_pass_r1)

cat(sprintf("  Genes passing TPM filter in all samples: %d/%d\n", 
            length(genes_pass_both), nrow(tpm_r0)))

# Check which pairs have both genes passing
min_tpm_pass <- (pair_info$gene_up %in% genes_pass_both) & 
  (pair_info$gene_down %in% genes_pass_both)

# Add to pair information
pair_info$min_tpm_pass <- min_tpm_pass

cat(sprintf("  Pairs passing minimum TPM filter: %s (%.1f%%)\n",
            format(sum(min_tpm_pass), big.mark = ","),
            sum(min_tpm_pass) / nrow(pair_info) * 100))

# ============================================================================
# Step 4: Progressive filtering with TPM constraint
# ============================================================================

cat("\n--- Step 4: Progressive filtering ---\n")

# Filter 1: Minimum TPM requirement
filter1_idx <- pair_info$min_tpm_pass
cat(sprintf("  Filter 1 (min TPM >= %.1f): %s pairs (%.1f%%)\n",
            MIN_TPM_THRESHOLD,
            format(sum(filter1_idx), big.mark = ","),
            sum(filter1_idx) / nrow(pair_info) * 100))

# Filter 2: Effect size threshold (increased to 2.0 for stronger signal)
filter2_idx <- filter1_idx & abs(pair_info$median_reo_r0) > 2.0
cat(sprintf("  Filter 2 (+ |median_r0| > 2.0): %s pairs (%.1f%%)\n", 
            format(sum(filter2_idx), big.mark = ","),
            sum(filter2_idx) / nrow(pair_info) * 100))

# Filter 3: R0 consistency
filter3_idx <- filter2_idx & (pair_info$r0_consistency >= 10)
cat(sprintf("  Filter 3 (+ R0 >= 10/11): %s pairs (%.1f%%)\n", 
            format(sum(filter3_idx), big.mark = ","),
            sum(filter3_idx) / nrow(pair_info) * 100))

# Filter 4: R1 reversal
filter4_idx <- filter3_idx & (pair_info$r1_reversal >= 10)
cat(sprintf("  Filter 4 (+ R1 >= 10/13): %s pairs (%.1f%%)\n", 
            format(sum(filter4_idx), big.mark = ","),
            sum(filter4_idx) / nrow(pair_info) * 100))

# ============================================================================
# Step 5: Extract and rank filtered pairs
# ============================================================================

cat("\n--- Step 5: Ranking filtered pairs ---\n")

# Extract filtered pairs
filtered_pairs <- pair_info[filter4_idx, ]
filtered_reo_r0 <- reo_matrix_r0[filter4_idx, ]
filtered_reo_r1 <- reo_matrix_r1[filter4_idx, ]

# Calculate additional metrics for ranking
filtered_pairs$r1_reversal_rate <- filtered_pairs$r1_reversal / 13
filtered_pairs$r0_consistency_rate <- filtered_pairs$r0_consistency / 11
filtered_pairs$abs_median_diff <- abs(filtered_pairs$median_diff)

# Ranking criteria (prioritize complete reversal)
filtered_pairs$rank_score <- (
  filtered_pairs$r1_reversal * 1000 +  # Primary: R1 reversal count
    filtered_pairs$r0_consistency * 10 +  # Secondary: R0 consistency
    filtered_pairs$abs_median_diff        # Tertiary: Effect size
)

# Sort by ranking score
filtered_pairs <- filtered_pairs[order(filtered_pairs$rank_score, decreasing = TRUE), ]
filtered_reo_r0 <- filtered_reo_r0[order(filtered_pairs$rank_score, decreasing = TRUE), ]
filtered_reo_r1 <- filtered_reo_r1[order(filtered_pairs$rank_score, decreasing = TRUE), ]

cat(sprintf("  Ranked %d candidate pairs\n", nrow(filtered_pairs)))

# ============================================================================
# Step 6: Display top pairs and select independent pairs
# ============================================================================

cat("\n--- Step 6: Top candidate pairs ---\n\n")

# Show top 20 pairs
n_display <- min(20, nrow(filtered_pairs))
for(i in 1:n_display) {
  cat(sprintf("%2d. %s/%s\n", i,
              filtered_pairs$gene_up[i],
              filtered_pairs$gene_down[i]))
  cat(sprintf("    R0: %d/11 (%.0f%%), R1: %d/13 reversed (%.0f%%)\n",
              filtered_pairs$r0_consistency[i],
              filtered_pairs$r0_consistency_rate[i] * 100,
              filtered_pairs$r1_reversal[i],
              filtered_pairs$r1_reversal_rate[i] * 100))
  cat(sprintf("    R0 median: %.2f, R1 median: %.2f, diff: %.2f\n",
              filtered_pairs$median_reo_r0[i],
              filtered_pairs$median_reo_r1[i],
              filtered_pairs$median_diff[i]))
  cat(sprintf("    DEG FC: up=%.2f, down=%.2f\n",
              filtered_pairs$gene_up_deg_fc[i],
              filtered_pairs$gene_down_deg_fc[i]))
}

# Perfect pairs (13/13 reversal)
perfect_pairs <- filtered_pairs[filtered_pairs$r1_reversal == 13, ]
cat(sprintf("\nPerfect reversal pairs (13/13): %d\n", nrow(perfect_pairs)))

# Select independent pairs (no gene overlap)
cat("\n--- Selecting independent pairs ---\n")
selected_pairs <- list()
used_genes <- character()

for(i in 1:nrow(filtered_pairs)) {
  gene_up <- filtered_pairs$gene_up[i]
  gene_down <- filtered_pairs$gene_down[i]
  
  # Select if neither gene has been used
  if(!gene_up %in% used_genes && !gene_down %in% used_genes) {
    selected_pairs[[length(selected_pairs) + 1]] <- filtered_pairs[i, ]
    used_genes <- c(used_genes, gene_up, gene_down)
    
    cat(sprintf("  Independent pair %d: %s/%s\n", 
                length(selected_pairs), gene_up, gene_down))
    cat(sprintf("    R0: %d/11, R1: %d/13, effect: %.2f\n",
                filtered_pairs$r0_consistency[i],
                filtered_pairs$r1_reversal[i],
                filtered_pairs$median_reo_r0[i]))
    
    # Stop after 5 independent pairs
    if(length(selected_pairs) >= 5) break
  }
}

cat(sprintf("\nSelected %d independent pairs for validation\n", length(selected_pairs)))

# ============================================================================
# Step 7: Create sample × pair matrix for downstream analysis
# ============================================================================

cat("\n--- Step 7: Creating sample × pair matrix ---\n")

# Combine R0 and R1 samples
sample_by_pair_matrix <- rbind(
  t(filtered_reo_r0),  # R0 samples (11 × n_pairs)
  t(filtered_reo_r1)   # R1 samples (13 × n_pairs)
)

# Set row names (samples) and column names (pairs)
rownames(sample_by_pair_matrix) <- c(phase2_data$r0_samples, phase2_data$r1_samples)
colnames(sample_by_pair_matrix) <- rownames(filtered_reo_r0)

# Create sample group labels
sample_groups <- c(rep("R0", 11), rep("R1", 13))
names(sample_groups) <- rownames(sample_by_pair_matrix)

cat(sprintf("  Created matrix: %d samples × %d pairs\n",
            nrow(sample_by_pair_matrix), ncol(sample_by_pair_matrix)))

# ============================================================================
# Step 9: Statistical feature selection using limma
# ============================================================================

cat("\n--- Step 9: Statistical feature selection ---\n")

# Load limma
suppressPackageStartupMessages(library(limma))

# Prepare data for limma analysis
# Matrix: samples (rows) x pairs (columns)
expr_matrix <- t(sample_by_pair_matrix)  # Now 330 pairs x 24 samples

# Design matrix
group <- factor(c(rep("R0", 11), rep("R1", 13)))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Apply array weights to handle R1 heterogeneity
cat("  Calculating sample weights to handle R1 heterogeneity...\n")
aw <- arrayWeights(expr_matrix, design)

cat("  Sample weights:\n")
cat(sprintf("    R0 samples: mean=%.2f, range=[%.2f, %.2f]\n", 
            mean(aw[1:11]), min(aw[1:11]), max(aw[1:11])))
cat(sprintf("    R1 samples: mean=%.2f, range=[%.2f, %.2f]\n", 
            mean(aw[12:24]), min(aw[12:24]), max(aw[12:24])))

# Identify potentially non-radiation-induced samples
low_weight_idx <- which(aw < 0.5)
if(length(low_weight_idx) > 0) {
  cat(sprintf("  Note: %d samples have low weights (<0.5), possibly natural carcinogenesis\n", 
              length(low_weight_idx)))
}

# Fit linear model with weights
cat("\n  Fitting weighted linear model...\n")
fit <- lmFit(expr_matrix, design, weights = aw)

# Define contrast (R1 - R0)
contrast.matrix <- makeContrasts(R1vsR0 = R1 - R0, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get top pairs by different criteria
cat("  Extracting top pairs by multiple criteria:\n")

# By B-statistic (log-odds of differential expression)
top_by_B <- topTable(fit2, number = 50, sort.by = "B")
cat(sprintf("    Top by B-statistic: %d pairs\n", nrow(top_by_B)))

# By adjusted p-value
top_by_p <- topTable(fit2, number = 50, sort.by = "p")
cat(sprintf("    Top by p-value: %d pairs\n", nrow(top_by_p)))

# Combine and rank
limma_selected <- unique(c(rownames(top_by_B), rownames(top_by_p)))
cat(sprintf("  Total unique pairs from limma: %d\n", length(limma_selected)))

# Get independent pairs from limma results
cat("\n  Selecting independent pairs from limma results:\n")
limma_pairs_info <- filtered_pairs[rownames(filtered_pairs) %in% limma_selected, ]
limma_pairs_info <- limma_pairs_info[order(top_by_B[rownames(limma_pairs_info), "B"], 
                                           decreasing = TRUE, na.last = TRUE), ]

# Select independent pairs
limma_independent <- list()
limma_used_genes <- character()

for(i in 1:nrow(limma_pairs_info)) {
  gene_up <- limma_pairs_info$gene_up[i]
  gene_down <- limma_pairs_info$gene_down[i]
  
  if(!gene_up %in% limma_used_genes && !gene_down %in% limma_used_genes) {
    limma_independent[[length(limma_independent) + 1]] <- limma_pairs_info[i, ]
    limma_used_genes <- c(limma_used_genes, gene_up, gene_down)
    
    # Get limma statistics
    pair_name <- rownames(limma_pairs_info)[i]
    if(pair_name %in% rownames(top_by_B)) {
      b_stat <- top_by_B[pair_name, "B"]
      p_adj <- top_by_B[pair_name, "adj.P.Val"]
    } else {
      b_stat <- NA
      p_adj <- NA
    }
    
    cat(sprintf("    Pair %d: %s/%s\n", 
                length(limma_independent), gene_up, gene_down))
    cat(sprintf("      B-stat: %.2f, adj.p: %.3e\n", b_stat, p_adj))
    
    if(length(limma_independent) >= 5) break
  }
}

cat(sprintf("\n  Selected %d independent pairs from limma\n", length(limma_independent)))

# ============================================================================
# Optional: Elastic Net + Stability Selection
# ============================================================================

cat("\n--- Optional: Elastic Net with Stability Selection ---\n")
cat("  (Can be run as an alternative or complementary approach)\n")

# Prepare for Elastic Net
suppressPackageStartupMessages(library(glmnet))

# Create response variable (0 for R0, 1 for R1)
y_binary <- c(rep(0, 11), rep(1, 13))

# Stability selection
cat("  Running stability selection (100 iterations)...\n")
stability_results <- replicate(100, {
  # Subsample 75% of data
  idx <- sample(1:24, 18)
  X_sub <- sample_by_pair_matrix[idx, ]
  y_sub <- y_binary[idx]
  
  # Fit elastic net (alpha=0.9, close to LASSO but more stable)
  cv_fit <- cv.glmnet(X_sub, y_sub, family = "binomial", alpha = 0.9, nfolds = 5)
  
  # Get selected features at lambda.min
  coef_vals <- coef(cv_fit, s = "lambda.min")
  selected <- which(coef_vals[-1] != 0)  # Exclude intercept
  
  return(selected)
})

# Count selection frequency
selection_freq <- table(unlist(stability_results))
stable_features <- as.numeric(names(selection_freq)[selection_freq >= 60])  # 60% threshold

cat(sprintf("  Features selected in >60%% of iterations: %d\n", length(stable_features)))

if(length(stable_features) > 0) {
  # Get pair names for stable features
  stable_pair_names <- colnames(sample_by_pair_matrix)[stable_features]
  
  # Compare with limma results
  overlap_with_limma <- sum(stable_pair_names %in% limma_selected)
  cat(sprintf("  Overlap with limma selection: %d pairs\n", overlap_with_limma))
  
  # Show top 5 most stable pairs
  top_stable <- head(sort(selection_freq[as.character(stable_features)], decreasing = TRUE), 5)
  cat("\n  Top 5 most stable pairs:\n")
  for(i in 1:length(top_stable)) {
    pair_idx <- as.numeric(names(top_stable)[i])
    pair_name <- colnames(sample_by_pair_matrix)[pair_idx]
    cat(sprintf("    %s (selected in %d%% of iterations)\n", 
                pair_name, top_stable[i]))
  }
}

# ============================================================================
# Step 10: Validation on R_test samples (updated to test limma-selected pairs)
# ============================================================================

cat("\n--- Step 10: Validation on R_test samples ---\n")

# Load R_test data
case_master_full <- readRDS(paste0(paths$processed, "thyr_case_master_full.rds"))
se <- readRDS(paste0(paths$processed, "thyr_se_strand2_nonzero.rds"))
sample_meta <- as.data.frame(colData(se))

# Get R_test cases
r_test_cases <- case_master_full[case_master_full$group == "R_test", ]
cat(sprintf("  R_test cases: %d\n", nrow(r_test_cases)))

# Get R_test tumor samples
r_test_tumor_samples <- sample_meta$sample_submitter_id[
  sample_meta$case_submitter_id %in% r_test_cases$case_submitter_id &
    sample_meta$sample_type == "Primary Tumor" &
    grepl("_merged", sample_meta$sample_submitter_id)
]
cat(sprintf("  R_test tumor samples: %d\n", length(r_test_tumor_samples)))

# Calculate TPM for R_test samples
r_test_counts <- assay(se)[, r_test_tumor_samples]
gene_lengths <- rowData(se)$gene_length
tpm_r_test <- t(t(r_test_counts / gene_lengths * 1000) / 
                  colSums(r_test_counts / gene_lengths * 1000) * 1e6)

# Get POC values for R_test samples
correct_case_ids <- sample_meta$case_submitter_id[
  match(r_test_tumor_samples, sample_meta$sample_submitter_id)
]
sample_poc <- r_test_cases$POC[
  match(correct_case_ids, r_test_cases$case_submitter_id)
]

# Categorize by POC
poc_groups <- cut(sample_poc, breaks = c(-Inf, 33.3, 66.6, Inf),
                  labels = c("Low(0-33%)", "Mid(33-67%)", "High(67%+)"))

cat("\n  POC distribution:\n")
print(table(poc_groups, useNA = "ifany"))

# Test independent pairs on R_test
if(length(selected_pairs) > 0) {
  cat("\n  Testing independent pairs:\n\n")
  
  for(j in 1:min(5, length(selected_pairs))) {
    pair <- selected_pairs[[j]]
    gene_up <- pair$gene_up
    gene_down <- pair$gene_down
    
    if(gene_up %in% rownames(tpm_r_test) & gene_down %in% rownames(tpm_r_test)) {
      reo_test <- log2(tpm_r_test[gene_up, ] / tpm_r_test[gene_down, ])
      expected_sign <- pair$r0_expected_direction
      
      cat(sprintf("  Pair %d: %s/%s\n", j, gene_up, gene_down))
      cat(sprintf("    Expected direction: %s\n", ifelse(expected_sign > 0, "positive", "negative")))
      
      # Test by POC group
      for(poc_level in levels(poc_groups)) {
        poc_idx <- which(!is.na(poc_groups) & poc_groups == poc_level)
        if(length(poc_idx) > 0) {
          reversal_rate <- sum(sign(reo_test[poc_idx]) != expected_sign) / length(poc_idx)
          cat(sprintf("    %s (n=%d): %.1f%% reversed\n",
                      poc_level, length(poc_idx), reversal_rate * 100))
        }
      }
    }
  }
  
  # Overall performance across all independent pairs
  cat("\n  Average reversal rates across independent pairs:\n")
  
  for(poc_level in levels(poc_groups)) {
    poc_idx <- which(!is.na(poc_groups) & poc_groups == poc_level)
    if(length(poc_idx) > 0) {
      all_reversals <- numeric()
      
      for(j in 1:min(5, length(selected_pairs))) {
        pair <- selected_pairs[[j]]
        gene_up <- pair$gene_up
        gene_down <- pair$gene_down
        
        if(gene_up %in% rownames(tpm_r_test) & gene_down %in% rownames(tpm_r_test)) {
          reo_test <- log2(tpm_r_test[gene_up, ] / tpm_r_test[gene_down, ])
          expected_sign <- pair$r0_expected_direction
          reversal_rate <- sum(sign(reo_test[poc_idx]) != expected_sign) / length(poc_idx)
          all_reversals <- c(all_reversals, reversal_rate)
        }
      }
      
      if(length(all_reversals) > 0) {
        avg_reversal <- mean(all_reversals) * 100
        cat(sprintf("    %s: %.1f%%\n", poc_level, avg_reversal))
      }
    }
  }
}

# ============================================================================
# Step 11: Save results
# ============================================================================

cat("\n--- Step 11: Saving results ---\n")

# Main results
phase3_results <- list(
  # Filtered pairs with all metrics
  filtered_pairs = filtered_pairs,
  
  # Independent pairs selected (original method)
  independent_pairs = do.call(rbind, selected_pairs),
  
  # Limma-selected independent pairs
  limma_independent_pairs = if(exists("limma_independent") && length(limma_independent) > 0) {
    do.call(rbind, limma_independent)
  } else {
    NULL
  },
  
  # Limma results
  limma_results = if(exists("top_by_B")) {
    list(
      top_by_B = top_by_B,
      sample_weights = aw,
      low_weight_samples = which(aw < 0.5)
    )
  } else {
    NULL
  },
  
  # REO matrices for filtered pairs
  reo_matrix_r0 = filtered_reo_r0,
  reo_matrix_r1 = filtered_reo_r1,
  
  # Sample × pair matrix
  sample_by_pair_matrix = sample_by_pair_matrix,
  sample_groups = sample_groups,
  
  # Metadata
  n_total_pairs = phase2_data$n_pairs,
  n_filtered_pairs = nrow(filtered_pairs),
  n_independent_pairs = length(selected_pairs),
  n_limma_pairs = if(exists("limma_independent")) length(limma_independent) else 0,
  filter_criteria = list(
    median_threshold = 2.0,
    r0_consistency_min = 10,
    r1_reversal_min = 10
  ),
  date_created = Sys.Date()
)

saveRDS(phase3_results, paste0(paths$processed, "reo_phase3_results.rds"))
cat(sprintf("  Saved to: reo_phase3_results.rds\n"))

# Export top pairs as CSV for review
top_100 <- head(filtered_pairs[, c("gene_up", "gene_down", 
                                   "r0_consistency", "r1_reversal",
                                   "median_diff", "gene_up_deg_fc", 
                                   "gene_down_deg_fc")], 100)

write.csv(top_100, paste0(paths$output, "reo_top100_pairs.csv"), row.names = FALSE)
cat("  Top 100 pairs exported to: output/reo_top100_pairs.csv\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n=== Phase 3 Complete ===\n")
cat("Summary:\n")
cat(sprintf("  Input: %s pairs\n", format(phase2_data$n_pairs, big.mark = ",")))
cat(sprintf("  After TPM filter (>= %.1f): %d pairs\n", MIN_TPM_THRESHOLD, sum(filter1_idx)))
cat(sprintf("  After effect filter (|median| > 2): %d pairs\n", sum(filter2_idx)))
cat(sprintf("  Final filtered: %d high-quality pairs\n", nrow(filtered_pairs)))
cat(sprintf("  Perfect reversals (13/13): %d pairs\n", nrow(perfect_pairs)))
cat(sprintf("  Independent pairs (original): %d\n", length(selected_pairs)))
if(exists("limma_independent")) {
  cat(sprintf("  Independent pairs (limma): %d\n", length(limma_independent)))
}

cat("\nFeature selection methods applied:\n")
cat("  1. Effect size and consistency filtering\n")
cat("  2. Weighted limma with arrayWeights (handles R1 heterogeneity)\n")

cat("\nNext steps:\n")
cat("  1. Compare limma-selected vs original pairs on R_test\n")
cat("  2. Validate on BRAF samples for specificity\n")
cat("  3. Select final 2-3 pairs based on:\n")
cat("     - POC dependence (Mid/Low ratio > 1.5)\n")
cat("     - Statistical significance (limma B-statistic)\n")
cat("     - Biological interpretability\n")

cat("\nNext steps:\n")
cat("  1. Test selected independent pairs on R_test samples\n")
cat("  2. Validate on BRAF samples for specificity\n")
cat("  3. Proceed with top 3-5 pairs for clinical validation\n")

# Clean up
rm(reo_matrix_r0, reo_matrix_r1)
gc()

cat("\n=== Script Complete ===\n")