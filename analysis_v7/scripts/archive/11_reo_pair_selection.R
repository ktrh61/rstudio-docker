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
# Date: 2025-12-07
# =============================================================================

source("analysis_v7/setup.R")

cat("\n=============================================================================\n")
cat("11_reo_pair_selection.R - REO Pair Selection\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
# Selection criteria:
#   - Dead zone: >= 22/24 samples with |r| >= log2(1.5)
#     → "Technical measurement error tolerance (~8% samples)"
#   - R0 consistency: >= 10/11 same sign
#     → "Allow 1 sample exception as biological variation"
#   - R1 reversal: Wilson 95% CI lower bound > 0.5
#     → "Statistical guarantee of reversal rate > 50%"
#
# Rationale for relaxation from Protocol v2:
#   - Strict conditions (DZ=24, R0=11) yielded only 1 candidate
#   - DZ relaxation: Technical tolerance for measurement noise
#   - R0 relaxation: Biological variation tolerance, adjusted via voting threshold T

PARAMS <- list(
  dead_zone_threshold = log2(1.5),  # 0.585
  dead_zone_min_samples = 22,       # out of 24 (~92%)
  r0_consistency_min = 10,          # out of 11 (~91%)
  wilson_alpha = 0.05,              # 95% confidence interval
  correlation_threshold = 0.75      # for pair independence (used in Phase 5)
)

cat("=== Parameters ===\n")
cat(sprintf("Dead zone threshold: log2(1.5) = %.3f\n", PARAMS$dead_zone_threshold))
cat(sprintf("Dead zone min samples: %d/24 (%.1f%%)\n", 
            PARAMS$dead_zone_min_samples, PARAMS$dead_zone_min_samples/24*100))
cat(sprintf("R0 consistency min: %d/11 (%.1f%%)\n", 
            PARAMS$r0_consistency_min, PARAMS$r0_consistency_min/11*100))
cat(sprintf("R1 reversal: Wilson %.0f%% CI lower > 0.5\n", (1-PARAMS$wilson_alpha)*100))
cat(sprintf("Correlation threshold: %.2f (for Phase 5)\n", PARAMS$correlation_threshold))
cat("\n")

# -----------------------------------------------------------------------------
# Wilson CI Function
# -----------------------------------------------------------------------------
# Wilson score interval for binomial proportion
# Returns lower bound of confidence interval
wilson_ci_lower <- function(x, n, alpha = 0.05) {
  if (n == 0) return(NA_real_)
  p_hat <- x / n
  z <- qnorm(1 - alpha/2)
  denominator <- 1 + z^2/n
  center <- (p_hat + z^2/(2*n)) / denominator
  margin <- z * sqrt(p_hat*(1-p_hat)/n + z^2/(4*n^2)) / denominator
  return(center - margin)
}

# Show reference values
cat("=== Wilson CI Reference (n=13, alpha=0.05) ===\n")
for (rev in 8:13) {
  ci_lower <- wilson_ci_lower(rev, 13, PARAMS$wilson_alpha)
  status <- ifelse(ci_lower > 0.5, "PASS", "FAIL")
  cat(sprintf("  %2d/13 reversal: CI lower = %.3f (%s)\n", rev, ci_lower, status))
}
cat("\n")

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
cat("=== Loading Data ===\n")

# Expression data
se <- readRDS(file.path(paths$processed, "thyr_se_strand2_nonzero.rds"))
cat(sprintf("Expression data: %d genes × %d samples\n", nrow(se), ncol(se)))

# DEG results
deg_results <- readRDS(file.path(paths$processed, "thyr_deg_results.rds"))
deg_df <- deg_results$deg_results$R0_vs_R1_tumor$deg_summary$results_df
cat(sprintf("DEG results: %d genes\n", nrow(deg_df)))

# Sample lists
sample_lists <- readRDS(file.path(paths$processed, "analysis_sample_lists.rds"))
r0_samples <- sample_lists$R0$tumor
r1_samples <- sample_lists$R1$tumor
cat(sprintf("R0 samples: %d\n", length(r0_samples)))
cat(sprintf("R1 samples: %d\n", length(r1_samples)))

cat("\n")

# -----------------------------------------------------------------------------
# Prepare TPM Matrix
# -----------------------------------------------------------------------------
cat("=== Preparing TPM Matrix ===\n")

# Get counts and gene lengths
counts <- assay(se, "counts")
gene_lengths <- rowData(se)$gene_length

# Calculate TPM
calc_tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

tpm_matrix <- apply(counts, 2, calc_tpm, lengths = gene_lengths)
rownames(tpm_matrix) <- rownames(counts)

# Log2 transform (data is already zero-free)
log2_tpm <- log2(tpm_matrix)

# Subset to R0 + R1 samples
all_samples <- c(r0_samples, r1_samples)
log2_tpm <- log2_tpm[, all_samples]

cat(sprintf("TPM matrix: %d genes × %d samples\n", nrow(log2_tpm), ncol(log2_tpm)))
cat("\n")

# -----------------------------------------------------------------------------
# Phase 2: Identify DEG UP and DOWN
# -----------------------------------------------------------------------------
cat("=== Phase 2: Candidate Genes ===\n")

# Filter significant DEGs
sig_deg <- deg_df[deg_df$significant == TRUE, ]
cat(sprintf("Significant DEGs: %d\n", nrow(sig_deg)))

# Split by direction (based on PI - Probability Index from Brunner-Munzel)
deg_up <- sig_deg[sig_deg$PI > 0.5, ]    # R1 > R0
deg_down <- sig_deg[sig_deg$PI < 0.5, ]  # R1 < R0

cat(sprintf("DEG UP (R1 > R0): %d genes\n", nrow(deg_up)))
cat(sprintf("DEG DOWN (R1 < R0): %d genes\n", nrow(deg_down)))

# Filter to genes present in expression matrix
deg_up_ids <- intersect(deg_up$gene_id, rownames(log2_tpm))
deg_down_ids <- intersect(deg_down$gene_id, rownames(log2_tpm))

cat(sprintf("DEG UP in expression matrix: %d genes\n", length(deg_up_ids)))
cat(sprintf("DEG DOWN in expression matrix: %d genes\n", length(deg_down_ids)))

# Remove genes with -Inf (zero expression in any sample)
# This ensures all pairs have valid r values for all samples
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
# Phase 3: Generate All UP × DOWN Pairs
# -----------------------------------------------------------------------------
cat("=== Phase 3: Pair Generation ===\n")

n_up <- length(deg_up_ids)
n_down <- length(deg_down_ids)
total_pairs <- n_up * n_down

cat(sprintf("Total possible pairs: %d × %d = %d\n", n_up, n_down, total_pairs))

# Create pair indices
pair_grid <- expand.grid(i = 1:n_up, j = 1:n_down)
cat(sprintf("Pair grid created: %d pairs\n", nrow(pair_grid)))
cat("\n")

# -----------------------------------------------------------------------------
# Phase 4: Apply Selection Criteria (Vectorized)
# -----------------------------------------------------------------------------
cat("=== Phase 4: Selection Criteria ===\n")

cat("Computing r values for all pairs...\n")

# Extract expression matrices for UP and DOWN genes
expr_up <- log2_tpm[deg_up_ids, , drop = FALSE]
expr_down <- log2_tpm[deg_down_ids, , drop = FALSE]

# Function to evaluate all pairs (vectorized, chunked)
evaluate_pairs <- function(pair_grid, expr_up, expr_down, r0_samples, r1_samples, params) {
  
  n_pairs <- nrow(pair_grid)
  n_r0 <- length(r0_samples)
  n_r1 <- length(r1_samples)
  
  # Process in chunks to manage memory
  chunk_size <- 10000
  n_chunks <- ceiling(n_pairs / chunk_size)
  
  results_list <- vector("list", n_chunks)
  
  for (chunk in 1:n_chunks) {
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, n_pairs)
    chunk_indices <- start_idx:end_idx
    n_chunk <- length(chunk_indices)
    
    # Get pair indices for this chunk
    i_chunk <- pair_grid$i[chunk_indices]
    j_chunk <- pair_grid$j[chunk_indices]
    
    # Compute r values: r = log2(UP) - log2(DOWN)
    r_r0 <- expr_up[i_chunk, r0_samples, drop = FALSE] - 
            expr_down[j_chunk, r0_samples, drop = FALSE]
    r_r1 <- expr_up[i_chunk, r1_samples, drop = FALSE] - 
            expr_down[j_chunk, r1_samples, drop = FALSE]
    
    # Dead zone check: |r| >= threshold
    dz_r0 <- rowSums(abs(r_r0) >= params$dead_zone_threshold)
    dz_r1 <- rowSums(abs(r_r1) >= params$dead_zone_threshold)
    dz_total <- dz_r0 + dz_r1
    
    # R0 consistency: count majority sign
    r0_pos <- rowSums(r_r0 > 0)
    r0_neg <- rowSums(r_r0 < 0)
    r0_consistency <- pmax(r0_pos, r0_neg)
    
    # Determine R0 majority sign for each pair
    r0_majority_sign <- ifelse(r0_pos > r0_neg, 1, -1)
    
    # R1 reversal: count samples with opposite sign to R0 majority
    r1_signs <- sign(r_r1)
    reversal_matrix <- sweep(r1_signs, 1, r0_majority_sign, FUN = function(x, y) x != y)
    reversal <- rowSums(reversal_matrix)
    
    results_list[[chunk]] <- data.frame(
      i = i_chunk,
      j = j_chunk,
      dz_total = dz_total,
      r0_consistency = r0_consistency,
      reversal = reversal
    )
    
    if (chunk %% 10 == 0) {
      cat(sprintf("  Processed chunk %d/%d\n", chunk, n_chunks))
    }
  }
  
  return(do.call(rbind, results_list))
}

# Evaluate all pairs
cat("Evaluating pairs...\n")
pair_results <- evaluate_pairs(pair_grid, expr_up, expr_down, 
                               r0_samples, r1_samples, PARAMS)

# Calculate Wilson CI lower bound for each pair
cat("Computing Wilson CI...\n")
n_r1 <- length(r1_samples)
pair_results$wilson_lower <- sapply(pair_results$reversal, wilson_ci_lower, 
                                    n = n_r1, alpha = PARAMS$wilson_alpha)

cat("\n=== Filtering Results ===\n")

# Apply filters progressively
cat(sprintf("Total pairs: %d\n", nrow(pair_results)))

# Filter 1: Dead zone
pass_dz <- pair_results$dz_total >= PARAMS$dead_zone_min_samples
cat(sprintf("Passed dead zone (>= %d/24): %d\n", 
            PARAMS$dead_zone_min_samples, sum(pass_dz)))

# Filter 2: R0 consistency
pass_r0 <- pair_results$r0_consistency >= PARAMS$r0_consistency_min
cat(sprintf("Passed dead zone + R0 consistency (>= %d/11): %d\n", 
            PARAMS$r0_consistency_min, sum(pass_dz & pass_r0)))

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

# Sort by reversal (descending), then R0 consistency, then DZ
candidates <- candidates[order(-candidates$reversal, 
                               -candidates$r0_consistency, 
                               -candidates$dz_total), ]

cat("\n=== Top Candidates ===\n")
print(head(candidates[, c("up_name", "down_name", "dz_total", 
                          "r0_consistency", "reversal", "wilson_lower")], 20))

# -----------------------------------------------------------------------------
# Reversal Distribution
# -----------------------------------------------------------------------------
cat("\n=== Reversal Distribution ===\n")
rev_table <- table(candidates$reversal)
print(rev_table)

cat(sprintf("\nPairs with 13/13 reversal: %d\n", sum(candidates$reversal == 13)))
cat(sprintf("Pairs with 12/13 reversal: %d\n", sum(candidates$reversal == 12)))
cat(sprintf("Pairs with 11/13 reversal: %d\n", sum(candidates$reversal == 11)))

# ZNF560 analysis
cat("\n=== ZNF560 Dependency ===\n")
znf560_pairs <- sum(candidates$up_name == "ZNF560")
non_znf560_pairs <- sum(candidates$up_name != "ZNF560")
cat(sprintf("Pairs with ZNF560: %d\n", znf560_pairs))
cat(sprintf("Pairs without ZNF560: %d\n", non_znf560_pairs))

if (non_znf560_pairs > 0) {
  cat("\nNon-ZNF560 candidates:\n")
  print(candidates[candidates$up_name != "ZNF560", 
                   c("up_name", "down_name", "dz_total", "r0_consistency", 
                     "reversal", "wilson_lower")])
}

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------
cat("\n=== Saving Results ===\n")

output <- list(
  params = PARAMS,
  deg_up_ids = deg_up_ids,
  deg_down_ids = deg_down_ids,
  candidates = candidates,
  all_results = pair_results,  # Keep for sensitivity analysis
  log2_tpm = log2_tpm,
  r0_samples = r0_samples,
  r1_samples = r1_samples,
  n_r0 = length(r0_samples),
  n_r1 = length(r1_samples),
  timestamp = Sys.time()
)

saveRDS(output, file.path(paths$processed, "reo_candidate_pairs.rds"))
cat(sprintf("Saved: %s\n", file.path(paths$processed, "reo_candidate_pairs.rds")))

# Summary CSV
write.csv(candidates[, c("up_name", "down_name", "up_gene", "down_gene",
                         "dz_total", "r0_consistency", "reversal", "wilson_lower")],
          file.path(paths$output, "reo_candidate_pairs.csv"),
          row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(paths$output, "reo_candidate_pairs.csv")))

cat("\n=== Phase 2-4 Complete ===\n")
cat(sprintf("Candidate pairs for Phase 5: %d\n", nrow(candidates)))
cat(sprintf("Non-ZNF560 pairs available: %d\n", non_znf560_pairs))
