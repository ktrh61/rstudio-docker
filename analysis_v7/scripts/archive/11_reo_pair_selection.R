#!/usr/bin/env Rscript
# =============================================================================
# 11_reo_pair_selection.R
# REO Panel Construction - Phase 2-4: Pair Generation and Selection
# 
# Purpose: Generate gene pairs and apply selection criteria
# Input: 
#   - thyr_se_strand2_nonzero.rds (expression data)
#   - thyr_deg_results.rds (DEG results)
#   - analysis_sample_lists.rds (R0/R1 sample lists)
# Output:
#   - reo_candidate_pairs.rds (filtered candidate pairs)
#
# Reference: REO_Panel_Protocol_v2.md
# Date: 2025-12-08
# v1.3: Added separation score for ranking
# v1.4: Fixed separation calculation for both R0 majority signs (+1 and -1)
# =============================================================================

source("analysis_v7/setup.R")

cat("\n=============================================================================\n")
cat("11_reo_pair_selection.R - REO Pair Selection\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
# Selection criteria and rationale:
#
# 1. Dead zone threshold: log2(1.5) = 0.585
#    - 1.5-fold difference is commonly used as biological significance threshold
#
# 2. Dead zone min samples: 22/23 (~96%)
#    - Allows 1 sample exception for technical noise
#
# 3. R0 consistency: 10/11 (~91%)
#    - Allows 1 sample exception for biological variation
#
# 4. Wilson CI lower > 0.5
#    - Statistical guarantee that true reversal rate exceeds 50%
#
# 5. Separation score
#    - Measures gap between R0 and R1 distributions at boundary
#    - For R0 negative (R1 positive expected): r1_k - max(R0)
#    - For R0 positive (R1 negative expected): min(R0) - r1_k
#    - k = 3: ignores 2 extreme R1 values (aligned with 10/12 Wilson CI)
#    - Positive = distributions separated for 10/12 reversal
#    - Used for ranking candidates (higher = more robust)

PARAMS <- list(
  dead_zone_threshold = log2(1.5),  # 0.585
  dead_zone_min_samples = 22,       # out of 23 (~96%)
  r0_consistency_min = 10,          # out of 11 (~91%)
  wilson_alpha = 0.05,              # 95% confidence interval
  correlation_threshold = 0.75      # for pair independence (used in Phase 5)
)

cat("=== Parameters ===\n")
cat(sprintf("Dead zone threshold: log2(1.5) = %.3f\n", PARAMS$dead_zone_threshold))
cat(sprintf("Dead zone min samples: %d\n", PARAMS$dead_zone_min_samples))
cat(sprintf("R0 consistency min: %d\n", PARAMS$r0_consistency_min))
cat(sprintf("R1 reversal: Wilson %.0f%% CI lower > 0.5\n", (1-PARAMS$wilson_alpha)*100))
cat("Ranking: separation score (accounts for R0 majority sign)\n")
cat("  separation > 0: R0/R1 distributions separated for 10/12 reversal\n")
cat("\n")

# -----------------------------------------------------------------------------
# Wilson CI Function
# -----------------------------------------------------------------------------
wilson_ci_lower <- function(x, n, alpha = 0.05) {
  if (n == 0) return(NA_real_)
  p_hat <- x / n
  z <- qnorm(1 - alpha/2)
  denominator <- 1 + z^2/n
  center <- (p_hat + z^2/(2*n)) / denominator
  margin <- z * sqrt(p_hat*(1-p_hat)/n + z^2/(4*n^2)) / denominator
  return(center - margin)
}

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

# Show Wilson CI reference values
cat(sprintf("=== Wilson CI Reference (n=%d, alpha=0.05) ===\n", n_r1))
for (rev in max(1, n_r1-4):n_r1) {
  ci_lower <- wilson_ci_lower(rev, n_r1, PARAMS$wilson_alpha)
  status <- ifelse(ci_lower > 0.5, "PASS", "FAIL")
  cat(sprintf("  %2d/%d reversal: CI lower = %.3f (%s)\n", rev, n_r1, ci_lower, status))
}
cat("\n")

# -----------------------------------------------------------------------------
# Prepare TPM Matrix
# -----------------------------------------------------------------------------
cat("=== Preparing TPM Matrix ===\n")

counts <- assay(se, "counts")
gene_lengths <- rowData(se)$gene_length

calc_tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

tpm_matrix <- apply(counts, 2, calc_tpm, lengths = gene_lengths)
rownames(tpm_matrix) <- rownames(counts)

log2_tpm <- log2(tpm_matrix)

all_samples <- c(r0_samples, r1_samples)
log2_tpm <- log2_tpm[, all_samples]

cat(sprintf("TPM matrix: %d genes x %d samples\n", nrow(log2_tpm), ncol(log2_tpm)))
cat("\n")

# -----------------------------------------------------------------------------
# Phase 2: Identify DEG UP and DOWN
# -----------------------------------------------------------------------------
cat("=== Phase 2: Candidate Genes ===\n")

sig_deg <- deg_df[deg_df$significant == TRUE, ]
cat(sprintf("Significant DEGs: %d\n", nrow(sig_deg)))

deg_up <- sig_deg[sig_deg$PI > 0.5, ]    # R1 > R0
deg_down <- sig_deg[sig_deg$PI < 0.5, ]  # R1 < R0

cat(sprintf("DEG UP (R1 > R0): %d genes\n", nrow(deg_up)))
cat(sprintf("DEG DOWN (R1 < R0): %d genes\n", nrow(deg_down)))

deg_up_ids <- intersect(deg_up$gene_id, rownames(log2_tpm))
deg_down_ids <- intersect(deg_down$gene_id, rownames(log2_tpm))

cat(sprintf("DEG UP in expression matrix: %d genes\n", length(deg_up_ids)))
cat(sprintf("DEG DOWN in expression matrix: %d genes\n", length(deg_down_ids)))

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
# Phase 4: Apply Selection Criteria (Vectorized)
# -----------------------------------------------------------------------------
cat("=== Phase 4: Selection Criteria ===\n")

cat("Computing r values and separation scores for all pairs...\n")

expr_up <- log2_tpm[deg_up_ids, , drop = FALSE]
expr_down <- log2_tpm[deg_down_ids, , drop = FALSE]

evaluate_pairs <- function(pair_grid, expr_up, expr_down, r0_samples, r1_samples, params) {
  
  n_pairs <- nrow(pair_grid)
  n_r0 <- length(r0_samples)
  n_r1 <- length(r1_samples)
  
  chunk_size <- 10000
  n_chunks <- ceiling(n_pairs / chunk_size)
  
  results_list <- vector("list", n_chunks)
  
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
    
    # Dead zone check
    dz_r0 <- rowSums(abs(r_r0) >= params$dead_zone_threshold)
    dz_r1 <- rowSums(abs(r_r1) >= params$dead_zone_threshold)
    dz_total <- dz_r0 + dz_r1
    
    # R0 consistency
    r0_pos <- rowSums(r_r0 > 0)
    r0_neg <- rowSums(r_r0 < 0)
    r0_consistency <- pmax(r0_pos, r0_neg)
    r0_majority_sign <- ifelse(r0_pos > r0_neg, 1, -1)
    
    # R1 reversal
    r1_signs <- sign(r_r1)
    reversal_matrix <- sweep(r1_signs, 1, r0_majority_sign, FUN = function(x, y) x != y)
    reversal <- rowSums(reversal_matrix)
    
    # Separation score
    # Aligned with Wilson CI requirement:
    #   Wilson CI > 0.5 requires >= 10/12 reversal
    #   So 2 samples are allowed to not reverse
    #   k = (allowed non-reversal) + 1 = 3
    #
    # For R0 majority sign = -1 (R0 negative, R1 positive expected):
    #   r0_boundary = max(R0_r)  -- most positive (closest to boundary)
    #   r1_boundary = k-th smallest R1_r  -- ignoring 2 lowest
    #   separation = r1_boundary - r0_boundary
    #
    # For R0 majority sign = +1 (R0 positive, R1 negative expected):
    #   r0_boundary = min(R0_r)  -- most negative (closest to boundary)
    #   r1_boundary = k-th largest R1_r  -- ignoring 2 highest
    #   separation = r0_boundary - r1_boundary
    #
    # In both cases, separation > 0 means distributions separated for 10/12 reversal
    
    min_reversal_for_wilson <- 10
    allowed_non_reversal <- n_r1 - min_reversal_for_wilson  # = 2
    k <- allowed_non_reversal + 1  # = 3
    
    r0_max <- apply(r_r0, 1, max)
    r0_min <- apply(r_r0, 1, min)
    r1_k_small <- apply(r_r1, 1, function(x) sort(x, partial = k)[k])
    r1_k_large <- apply(r_r1, 1, function(x) -sort(-x, partial = k)[k])
    
    # Calculate separation based on R0 majority sign
    separation <- ifelse(
      r0_majority_sign == -1,
      r1_k_small - r0_max,  # R0 negative case
      r0_min - r1_k_large   # R0 positive case
    )
    
    # Store boundary values for reference
    r0_boundary <- ifelse(r0_majority_sign == -1, r0_max, r0_min)
    r1_boundary <- ifelse(r0_majority_sign == -1, r1_k_small, r1_k_large)
    
    results_list[[chunk]] <- data.frame(
      i = i_chunk,
      j = j_chunk,
      dz_total = dz_total,
      r0_consistency = r0_consistency,
      r0_majority_sign = r0_majority_sign,
      reversal = reversal,
      r0_boundary = r0_boundary,
      r1_boundary = r1_boundary,
      separation = separation
    )
    
    if (chunk %% 10 == 0) {
      cat(sprintf("  Processed chunk %d/%d\n", chunk, n_chunks))
    }
  }
  
  return(do.call(rbind, results_list))
}

cat("Evaluating pairs...\n")
pair_results <- evaluate_pairs(pair_grid, expr_up, expr_down, 
                               r0_samples, r1_samples, PARAMS)

# Calculate Wilson CI
cat("Computing Wilson CI...\n")
pair_results$wilson_lower <- sapply(pair_results$reversal, wilson_ci_lower, 
                                    n = n_r1, alpha = PARAMS$wilson_alpha)

cat("\n=== Filtering Results ===\n")

cat(sprintf("Total pairs: %d\n", nrow(pair_results)))

# Filter 1: Dead zone
pass_dz <- pair_results$dz_total >= PARAMS$dead_zone_min_samples
cat(sprintf("Passed dead zone (>= %d/%d): %d\n", 
            PARAMS$dead_zone_min_samples, n_r0 + n_r1, sum(pass_dz)))

# Filter 2: R0 consistency
pass_r0 <- pair_results$r0_consistency >= PARAMS$r0_consistency_min
cat(sprintf("Passed dead zone + R0 consistency (>= %d/%d): %d\n", 
            PARAMS$r0_consistency_min, n_r0, sum(pass_dz & pass_r0)))

# Filter 3: Wilson CI
pass_wilson <- pair_results$wilson_lower > 0.5
cat(sprintf("Passed all + Wilson CI lower > 0.5: %d\n", 
            sum(pass_dz & pass_r0 & pass_wilson)))

# Combined filter
pass_all <- pass_dz & pass_r0 & pass_wilson
cat(sprintf("\nTotal candidates: %d\n", sum(pass_all)))

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

# Sort by separation score (descending) - primary ranking criterion
candidates <- candidates[order(-candidates$separation), ]

cat("\n=== Top Candidates (by Separation Score) ===\n")
print(head(candidates[, c("up_name", "down_name", "separation", 
                          "reversal", "r0_consistency", "dz_total")], 20))

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------
cat("\n=== Separation Score Distribution ===\n")
cat("  (gap between R0 and R1 boundary values, k=3)\n")
cat(sprintf("  Range: %.2f to %.2f\n", min(candidates$separation), max(candidates$separation)))
cat(sprintf("  Mean: %.2f\n", mean(candidates$separation)))
cat(sprintf("  Median: %.2f\n", median(candidates$separation)))
cat(sprintf("  Positive (separated for 10/12 reversal): %d (%.1f%%)\n", 
            sum(candidates$separation > 0),
            sum(candidates$separation > 0) / nrow(candidates) * 100))

# R0 majority sign distribution
cat("\n=== R0 Majority Sign Distribution ===\n")
sign_table <- table(candidates$r0_majority_sign)
cat(sprintf("  Negative (UP < DOWN in R0): %d\n", 
            sum(candidates$r0_majority_sign == -1)))
cat(sprintf("  Positive (UP > DOWN in R0): %d\n", 
            sum(candidates$r0_majority_sign == 1)))

cat("\n=== Reversal Distribution ===\n")
rev_table <- table(candidates$reversal)
print(rev_table)

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------
cat("\n=== Saving Results ===\n")

output <- list(
  params = PARAMS,
  deg_up_ids = deg_up_ids,
  deg_down_ids = deg_down_ids,
  candidates = candidates,
  all_results = pair_results,
  log2_tpm = log2_tpm,
  r0_samples = r0_samples,
  r1_samples = r1_samples,
  n_r0 = length(r0_samples),
  n_r1 = length(r1_samples),
  timestamp = Sys.time()
)

saveRDS(output, file.path(paths$processed, "reo_candidate_pairs.rds"))
cat(sprintf("Saved: %s\n", file.path(paths$processed, "reo_candidate_pairs.rds")))

write.csv(candidates[, c("up_name", "down_name", "up_gene", "down_gene",
                         "separation", "r0_boundary", "r1_boundary", "reversal", 
                         "r0_consistency", "dz_total", "wilson_lower")],
          file.path(paths$output, "reo_candidate_pairs.csv"),
          row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_candidate_pairs.csv")))

cat("\n=== Phase 2-4 Complete ===\n")
cat(sprintf("Candidate pairs for Phase 5: %d\n", nrow(candidates)))