#!/usr/bin/env Rscript
# =============================================================================
# 11_reo_pair_selection.R
# REO Panel Construction - Pair Generation and Selection
# 
# Purpose: Generate gene pairs from DEGs and apply selection criteria
# Input: 
#   - thyr_se_strand2_nonzero.rds (expression data)
#   - thyr_deg_results.rds (DEG results)
#   - analysis_sample_lists.rds (R0/R1 sample lists)
# Output:
#   - reo_candidate_pairs.rds (filtered candidate pairs)
#   - reo_candidate_pairs.csv
#
# Version: v7.1
# Date: 2025-12-09
#
# Selection Criteria (New Policy):
#   1. Dead zone: |r| < log2(1.2) → treated as "non-reversal"
#   2. R0 conditions:
#      - Sign consistency: all or max 1 exception (among dead zone外)
#      - Strength: q10(|r|) >= log2(1.5) (all samples)
#   3. R1 conditions:
#      - Majority reversal: >50% reversal (denominator = total samples)
#      - No full reversal: reversal < 100% (ban overfitting)
#   4. Ranking: |median_r_R0 - median_r_R1| descending
#   5. Constraints: no gene reuse, no high correlation pairs
# =============================================================================

source("analysis_v7/setup.R")

cat("\n=============================================================================\n")
cat("11_reo_pair_selection.R - REO Pair Selection (v7.1)\n")
cat("Date: 2025-12-09\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
PARAMS <- list(
  
  # Dead zone threshold
  dead_zone = log2(1.2),             # 0.263 - below this = non-reversal
  
  # R0 conditions
  r0_exception_max = 1,              # max exceptions for sign consistency
  r0_q10_threshold = log2(1.5),      # strength criterion (10th percentile)
  
  # R1 conditions
  r1_majority_required = TRUE,       # >50% reversal required
  r1_full_reversal_ban = TRUE,       # 100% reversal prohibited
  
  # Ranking method
  ranking_method = "median_diff",    # |median_R0 - median_R1|
  
  # For later use (Phase 5)
  correlation_threshold = 0.75       # Spearman correlation threshold
)

cat("=== Parameters ===\n")
cat(sprintf("Dead zone: |r| < log2(1.2) = %.4f\n", PARAMS$dead_zone))
cat(sprintf("R0 sign consistency: max %d exception(s) (among valid samples)\n", 
            PARAMS$r0_exception_max))
cat(sprintf("R0 strength (q10): >= log2(1.5) = %.4f\n", PARAMS$r0_q10_threshold))
cat(sprintf("R1 reversal: >50%% required, 100%% prohibited\n"))
cat(sprintf("Ranking: |median_R0 - median_R1| descending\n"))
cat("\n")

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
cat("=== Loading Data ===\n")

# Expression data
se <- readRDS(file.path(paths$processed, "thyr_se_strand2_nonzero.rds"))
cat(sprintf("Expression data: %d genes x %d samples\n", nrow(se), ncol(se)))

# DEG results
deg_results <- readRDS(file.path(paths$processed, "thyr_deg_results.rds"))
deg_df <- deg_results$deg_results$R0_vs_R1_tumor$deg_summary$results_df
cat(sprintf("DEG results: %d genes\n", nrow(deg_df)))

# Sample lists
sample_lists <- readRDS(file.path(paths$processed, "analysis_sample_lists.rds"))
r0_samples <- sample_lists$R0$tumor
r1_samples <- sample_lists$R1$tumor
n_r0 <- length(r0_samples)
n_r1 <- length(r1_samples)
cat(sprintf("R0 samples: %d\n", n_r0))
cat(sprintf("R1 samples: %d\n", n_r1))

cat("\n")

# -----------------------------------------------------------------------------
# Phase 2: Identify DEG UP and DOWN (before TPM to filter genes)
# -----------------------------------------------------------------------------
cat("=== Phase 2: Candidate Genes ===\n")

sig_deg <- deg_df[deg_df$significant == TRUE, ]
cat(sprintf("Significant DEGs: %d\n", nrow(sig_deg)))

deg_up <- sig_deg[sig_deg$PI > 0.5, ]    # R1 > R0
deg_down <- sig_deg[sig_deg$PI < 0.5, ]  # R1 < R0

cat(sprintf("DEG UP (R1 > R0): %d genes\n", nrow(deg_up)))
cat(sprintf("DEG DOWN (R1 < R0): %d genes\n", nrow(deg_down)))

# Get DEG gene IDs present in SE
deg_up_ids <- intersect(deg_up$gene_id, rownames(se))
deg_down_ids <- intersect(deg_down$gene_id, rownames(se))

cat(sprintf("DEG UP in SE: %d genes\n", length(deg_up_ids)))
cat(sprintf("DEG DOWN in SE: %d genes\n", length(deg_down_ids)))

# -----------------------------------------------------------------------------
# Prepare TPM Matrix (DEG genes x R0/R1 samples only)
# -----------------------------------------------------------------------------
cat("\n=== Preparing TPM Matrix (filtered) ===\n")

# Filter to required genes and samples only
required_genes <- unique(c(deg_up_ids, deg_down_ids))
required_samples <- c(r0_samples, r1_samples)

cat(sprintf("Required genes: %d\n", length(required_genes)))
cat(sprintf("Required samples: %d\n", length(required_samples)))

# Extract subset of counts
counts_subset <- assay(se, "counts")[required_genes, required_samples, drop = FALSE]
gene_lengths_subset <- rowData(se)$gene_length[match(required_genes, rownames(se))]

calc_tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

tpm_matrix <- apply(counts_subset, 2, calc_tpm, lengths = gene_lengths_subset)
rownames(tpm_matrix) <- required_genes

log2_tpm <- log2(tpm_matrix)

cat(sprintf("log2(TPM) matrix: %d genes x %d samples\n", nrow(log2_tpm), ncol(log2_tpm)))

# Remove genes with -Inf
has_inf_up <- sapply(deg_up_ids, function(g) any(is.infinite(log2_tpm[g, ])))
has_inf_down <- sapply(deg_down_ids, function(g) any(is.infinite(log2_tpm[g, ])))

n_inf_up <- sum(has_inf_up)
n_inf_down <- sum(has_inf_down)

if (n_inf_up > 0 || n_inf_down > 0) {
  cat(sprintf("\nRemoving genes with zero expression (log2 = -Inf):\n"))
  cat(sprintf("  DEG UP: %d genes removed\n", n_inf_up))
  cat(sprintf("  DEG DOWN: %d genes removed\n", n_inf_down))
  
  deg_up_ids <- deg_up_ids[!has_inf_up]
  deg_down_ids <- deg_down_ids[!has_inf_down]
  
  cat(sprintf("\nAfter removal:\n"))
  cat(sprintf("  DEG UP: %d genes\n", length(deg_up_ids)))
  cat(sprintf("  DEG DOWN: %d genes\n", length(deg_down_ids)))
}
cat("\n")

# -----------------------------------------------------------------------------
# Phase 3: Generate All UP x DOWN Pairs
# -----------------------------------------------------------------------------
cat("=== Phase 3: Pair Generation ===\n")

n_up <- length(deg_up_ids)
n_down <- length(deg_down_ids)
total_pairs <- n_up * n_down

cat(sprintf("Total possible pairs: %d x %d = %d\n", n_up, n_down, total_pairs))

pair_grid <- expand.grid(i = 1:n_up, j = 1:n_down)
cat(sprintf("Pair grid created: %d pairs\n", nrow(pair_grid)))
cat("\n")

# -----------------------------------------------------------------------------
# Phase 4: Apply Selection Criteria (New Policy)
# -----------------------------------------------------------------------------
cat("=== Phase 4: Selection Criteria (New Policy v7.1) ===\n")

cat("Evaluating pairs with new criteria...\n")

expr_up <- log2_tpm[deg_up_ids, , drop = FALSE]
expr_down <- log2_tpm[deg_down_ids, , drop = FALSE]

# Vectorized evaluation function
evaluate_pairs_v71 <- function(pair_grid, expr_up, expr_down, 
                               r0_samples, r1_samples, params) {
  
  n_pairs <- nrow(pair_grid)
  n_r0 <- length(r0_samples)
  n_r1 <- length(r1_samples)
  
  # q10 position for n_r0 samples (0.10 quantile)
  # For n=11: pos = 0.10 * 10 + 1 = 2 (2nd smallest)
  q10_pos <- max(1, ceiling(n_r0 * 0.10))
  
  # median position
  median_pos_r0 <- (n_r0 + 1) / 2
  median_pos_r1 <- (n_r1 + 1) / 2
  
  chunk_size <- 10000  # Conservative chunk size
  n_chunks <- ceiling(n_pairs / chunk_size)
  
  results_list <- vector("list", n_chunks)
  
  cat(sprintf("  Total chunks: %d (chunk size: %d)\n", n_chunks, chunk_size))
  
  for (chunk in 1:n_chunks) {
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, n_pairs)
    chunk_indices <- start_idx:end_idx
    n_chunk <- length(chunk_indices)
    
    i_chunk <- pair_grid$i[chunk_indices]
    j_chunk <- pair_grid$j[chunk_indices]
    
    # Compute r values: r = log2(UP) - log2(DOWN)
    r_r0 <- expr_up[i_chunk, r0_samples, drop = FALSE] - 
      expr_down[j_chunk, r0_samples, drop = FALSE]
    r_r1 <- expr_up[i_chunk, r1_samples, drop = FALSE] - 
      expr_down[j_chunk, r1_samples, drop = FALSE]
    
    # --- R0 Sign Consistency (among dead zone外 samples) ---
    valid_r0 <- abs(r_r0) >= params$dead_zone
    
    r0_pos_valid <- rowSums((r_r0 > 0) & valid_r0)
    r0_neg_valid <- rowSums((r_r0 < 0) & valid_r0)
    n_valid_r0 <- rowSums(valid_r0)
    
    r0_majority_sign <- ifelse(r0_pos_valid >= r0_neg_valid, 1, -1)
    r0_exceptions <- pmin(r0_pos_valid, r0_neg_valid)
    pass_consistency <- r0_exceptions <= params$r0_exception_max
    
    # --- R0 Strength: q10 criterion ---
    # Use partial sort for efficiency (only find q10_pos-th smallest)
    abs_r_r0 <- abs(r_r0)
    q10_values <- numeric(n_chunk)
    for (row in 1:n_chunk) {
      q10_values[row] <- sort(abs_r_r0[row, ], partial = q10_pos)[q10_pos]
    }
    pass_q10 <- q10_values >= params$r0_q10_threshold
    
    # --- R1 Reversal ---
    outside_dz_r1 <- abs(r_r1) >= params$dead_zone
    opposite_sign <- sign(r_r1) != r0_majority_sign
    reversal_count <- rowSums(outside_dz_r1 & opposite_sign)
    reversal_rate <- reversal_count / n_r1
    pass_majority <- reversal_rate > 0.5
    pass_not_all <- reversal_count < n_r1
    
    # --- Ranking: |median_R0 - median_R1| ---
    # Use partial sort for median (faster than apply + median)
    median_r0 <- numeric(n_chunk)
    median_r1 <- numeric(n_chunk)
    mid_r0 <- ceiling(n_r0 / 2)
    mid_r1 <- ceiling(n_r1 / 2)
    
    for (row in 1:n_chunk) {
      sorted_r0 <- sort(r_r0[row, ], partial = mid_r0)
      sorted_r1 <- sort(r_r1[row, ], partial = mid_r1)
      if (n_r0 %% 2 == 1) {
        median_r0[row] <- sorted_r0[mid_r0]
      } else {
        median_r0[row] <- (sorted_r0[mid_r0] + sorted_r0[mid_r0 + 1]) / 2
      }
      if (n_r1 %% 2 == 1) {
        median_r1[row] <- sorted_r1[mid_r1]
      } else {
        median_r1[row] <- (sorted_r1[mid_r1] + sorted_r1[mid_r1 + 1]) / 2
      }
    }
    median_diff <- abs(median_r0 - median_r1)
    
    # Store results
    results_list[[chunk]] <- data.frame(
      i = i_chunk,
      j = j_chunk,
      # R0 metrics
      n_valid_r0 = n_valid_r0,
      r0_majority_sign = r0_majority_sign,
      r0_exceptions = r0_exceptions,
      q10 = q10_values,
      # R1 metrics
      reversal_count = reversal_count,
      reversal_rate = reversal_rate,
      # Ranking
      median_r0 = median_r0,
      median_r1 = median_r1,
      median_diff = median_diff,
      # Pass flags
      pass_consistency = pass_consistency,
      pass_q10 = pass_q10,
      pass_majority = pass_majority,
      pass_not_all = pass_not_all
    )
    
    if (chunk %% 10 == 0 || chunk == n_chunks) {
      cat(sprintf("  Processed chunk %d/%d\n", chunk, n_chunks))
    }
  }
  
  return(do.call(rbind, results_list))
}

cat("Evaluating pairs...\n")
pair_results <- evaluate_pairs_v71(pair_grid, expr_up, expr_down, 
                                   r0_samples, r1_samples, PARAMS)

cat("\n=== Filtering Results ===\n")

cat(sprintf("Total pairs: %d\n", nrow(pair_results)))

# Filter summary
cat(sprintf("\n--- Filter Step by Step ---\n"))
cat(sprintf("Pass R0 consistency (<= %d exception): %d (%.1f%%)\n", 
            PARAMS$r0_exception_max,
            sum(pair_results$pass_consistency),
            sum(pair_results$pass_consistency) / nrow(pair_results) * 100))

cat(sprintf("Pass R0 q10 (>= %.3f): %d (%.1f%%)\n", 
            PARAMS$r0_q10_threshold,
            sum(pair_results$pass_q10),
            sum(pair_results$pass_q10) / nrow(pair_results) * 100))

cat(sprintf("Pass R1 majority (>50%%): %d (%.1f%%)\n", 
            sum(pair_results$pass_majority),
            sum(pair_results$pass_majority) / nrow(pair_results) * 100))

cat(sprintf("Pass R1 not-all (<100%%): %d (%.1f%%)\n", 
            sum(pair_results$pass_not_all),
            sum(pair_results$pass_not_all) / nrow(pair_results) * 100))

# Combined filter
pass_all <- pair_results$pass_consistency & 
  pair_results$pass_q10 & 
  pair_results$pass_majority & 
  pair_results$pass_not_all

cat(sprintf("\n--- Combined ---\n"))
cat(sprintf("Pass all criteria: %d (%.2f%%)\n", sum(pass_all), 
            sum(pass_all) / nrow(pair_results) * 100))

# Check if any candidates remain
if (sum(pass_all) == 0) {
  cat("\n!!! WARNING: No candidates passed all criteria !!!\n")
  cat("Consider relaxing parameters.\n")
  
  # Save diagnostic info
  saveRDS(list(
    params = PARAMS,
    pair_results = pair_results,
    pass_summary = data.frame(
      criterion = c("consistency", "q10", "majority", "not_all", "all"),
      n_pass = c(sum(pair_results$pass_consistency),
                 sum(pair_results$pass_q10),
                 sum(pair_results$pass_majority),
                 sum(pair_results$pass_not_all),
                 sum(pass_all))
    )
  ), file.path(paths$processed, "reo_diagnostic_no_candidates.rds"))
  
  stop("No candidates found. Diagnostic info saved.")
}

# Extract passing pairs
candidates <- pair_results[pass_all, ]
candidates$up_gene <- deg_up_ids[candidates$i]
candidates$down_gene <- deg_down_ids[candidates$j]

# Add gene names
candidates$up_name <- sapply(candidates$up_gene, function(g) {
  deg_df$gene_name[deg_df$gene_id == g]
})
candidates$down_name <- sapply(candidates$down_gene, function(g) {
  deg_df$gene_name[deg_df$gene_id == g]
})

# Sort by median_diff (descending) - primary ranking criterion
candidates <- candidates[order(-candidates$median_diff), ]

# Add rank
candidates$rank <- 1:nrow(candidates)

cat("\n=== Top Candidates (by |median_R0 - median_R1|) ===\n")
print(head(candidates[, c("rank", "up_name", "down_name", "median_diff", 
                          "reversal_count", "reversal_rate", "q10")], 20))

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------
cat("\n=== Candidate Statistics ===\n")

cat(sprintf("Median diff range: %.3f to %.3f\n", 
            min(candidates$median_diff), max(candidates$median_diff)))
cat(sprintf("Median diff mean: %.3f\n", mean(candidates$median_diff)))

cat(sprintf("\nReversal count distribution:\n"))
print(table(candidates$reversal_count))

cat(sprintf("\nq10 range: %.3f to %.3f\n", 
            min(candidates$q10), max(candidates$q10)))

cat(sprintf("\nR0 exceptions distribution:\n"))
print(table(candidates$r0_exceptions))

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------
cat("\n=== Saving Results ===\n")

output <- list(
  # Parameters
  params = PARAMS,
  version = "v7.1",
  date = "2025-12-09",
  
  # Gene lists
  deg_up_ids = deg_up_ids,
  deg_down_ids = deg_down_ids,
  
  # Results
  candidates = candidates,
  all_results = pair_results,
  
  # Expression data
  log2_tpm = log2_tpm,
  
  # Sample info
  r0_samples = r0_samples,
  r1_samples = r1_samples,
  n_r0 = n_r0,
  n_r1 = n_r1,
  
  # Metadata
  timestamp = Sys.time()
)

saveRDS(output, file.path(paths$processed, "reo_candidate_pairs.rds"))
cat(sprintf("Saved: %s\n", file.path(paths$processed, "reo_candidate_pairs.rds")))

# CSV output
csv_cols <- c("rank", "up_name", "down_name", "up_gene", "down_gene",
              "median_diff", "median_r0", "median_r1",
              "reversal_count", "reversal_rate", "q10",
              "r0_exceptions", "n_valid_r0")
write.csv(candidates[, csv_cols],
          file.path(paths$output, "reo_candidate_pairs.csv"),
          row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_candidate_pairs.csv")))

cat("\n=== Phase 2-4 Complete ===\n")
cat(sprintf("Candidate pairs for Phase 5: %d\n", nrow(candidates)))
cat(sprintf("Ranking criterion: |median_R0 - median_R1| (descending)\n"))