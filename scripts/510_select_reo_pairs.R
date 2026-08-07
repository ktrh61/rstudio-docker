# 510_select_reo_pairs.R
# REO (relative expression ordering) candidate-pair generation and selection,
# for the R cohort (RET/PTC) tumour comparison High vs Sporadic.
# Input : processed/thyr_expression_test.rds     (from 410; per-gene BM effect)
#         processed/thyr_normalized_counts.rds   (from 310; R_Tumor sample set)
#         processed/thyr_se_raw.rds              (from 120; single count assay)
#         processed/gene_lengths.rds             (from 020; unversioned lengths)
# Output: processed/thyr_reo_candidate_pairs.rds, output/reo_candidate_pairs.csv
#
# A REO pair (up-gene, down-gene) has a stable within-sample order in Sporadic
# (R0) and flips in High (R1). Because 410's exact permutation BM leaves no
# thresholded DEG list, the candidate gene pool is the top N genes by BM effect
# magnitude |effect - 0.5| (not a DEG cut), split into up (effect > 0.5, higher
# in High) and down (effect < 0.5). The pool only feeds pair generation; the
# reversal criteria below are the actual filter and are DE-independent.
#
# Expression is TPM (length-normalized), required because REO compares gene i vs
# gene j WITHIN a sample and raw/CPM counts confound that with gene length. TPM
# is computed here from the SE's single count assay (the strand chosen
# upstream in 110-120; the SE stores no TPM). REO is
# normalization-free (within-sample ranks), so DEGES normalization is irrelevant.
#
# Selection (following the v7 policy):
#   r = log2(TPM_up) - log2(TPM_down) per sample; dead zone |r| < log2(1.2).
#   R0 (Sporadic): sign consistent among non-dead-zone samples (<= 1 exception)
#     AND q10(|r|) >= log2(1.5) (a robust reference order).
#   R1 (High)   : > 50% of samples reverse the R0 reference sign, but not 100%
#     (a full reversal is banned as overfitting).
#   Rank by |median(r_R0) - median(r_R1)| descending.
# Redundancy removal (gene reuse, correlated pairs) is deferred to 520.

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "reo.R"))
source(file.path(paths$root, "lib", "annotation.R"))

# --- Configuration ---------------------------------------------------------
N_CANDIDATES <- 500L # top genes by |effect - 0.5| (~half up, half down)
PARAMS <- list(
  dead_zone = DEAD_ZONE,        # from config.R
  r0_exception_max = 1L,        # max sign exceptions for R0 consistency
  r0_q10_threshold = log2(1.5), # R0 strength: 10th percentile of |r|
  chunk_size = 10000L
)

# --- Load inputs -----------------------------------------------------------
test_path <- file.path(paths$processed, "thyr_expression_test.rds")
norm_path <- file.path(paths$processed, "thyr_normalized_counts.rds")
se_path <- file.path(paths$processed, "thyr_se_raw.rds")
len_path <- file.path(paths$processed, "gene_lengths.rds")
for (p in c(test_path, norm_path, se_path, len_path)) {
  if (!file.exists(p)) stop("missing input: ", p)
}
genes <- readRDS(test_path)$units$R_Tumor$genes
dge <- readRDS(norm_path)$units$R_Tumor$dgelist
se <- readRDS(se_path)
gene_lengths <- readRDS(len_path)

# R0 (Sporadic) and R1 (High) tumour sample IDs, from the analysed R_Tumor unit.
arms <- unit_arms(dge$samples$group, "R_Tumor")
r0_samples <- colnames(dge)[arms$sporadic]
r1_samples <- colnames(dge)[arms$high]
message("R0 (Sporadic) tumour: ", length(r0_samples),
  " ; R1 (High) tumour: ", length(r1_samples))

# --- log2 TPM from the single count assay (shared definition) ---------------
log2_tpm <- reo_log2_tpm(se, gene_lengths, c(r0_samples, r1_samples))
message("TPM matrix: ", nrow(log2_tpm), " length-annotated genes x ", ncol(log2_tpm), " samples")

# --- Candidate gene pool: top N by |effect - 0.5|, split up / down ---------
genes$ensembl <- strip_ensembl_version(genes$gene_id)
genes <- genes[genes$ensembl %in% rownames(log2_tpm), , drop = FALSE]
genes <- genes[order(abs(genes$effect - 0.5), decreasing = TRUE), , drop = FALSE]
pool <- head(genes, N_CANDIDATES)
up_ids <- pool$ensembl[pool$effect > 0.5]
down_ids <- pool$ensembl[pool$effect < 0.5]

# Drop candidates with any zero-expression (log2 TPM = -Inf) sample.
finite_gene <- function(ids) ids[apply(log2_tpm[ids, , drop = FALSE], 1, function(x) all(is.finite(x)))]
up_ids <- finite_gene(up_ids)
down_ids <- finite_gene(down_ids)
message(sprintf(
  "Candidate pool (top %d by |effect-0.5|): up %d, down %d",
  N_CANDIDATES, length(up_ids), length(down_ids)
))

# --- Evaluate all up x down pairs against the reversal criteria -------------
expr_up <- log2_tpm[up_ids, r0_samples, drop = FALSE]
expr_up_r1 <- log2_tpm[up_ids, r1_samples, drop = FALSE]
expr_down <- log2_tpm[down_ids, r0_samples, drop = FALSE]
expr_down_r1 <- log2_tpm[down_ids, r1_samples, drop = FALSE]
n_r0 <- length(r0_samples)
n_r1 <- length(r1_samples)
q10_pos <- max(1L, ceiling(n_r0 * 0.10))

pair_grid <- expand.grid(i = seq_along(up_ids), j = seq_along(down_ids))
n_pairs <- nrow(pair_grid)
message("Evaluating ", n_pairs, " up x down pairs")

row_median <- function(mat) {
  apply(mat, 1, stats::median)
}

n_chunks <- ceiling(n_pairs / PARAMS$chunk_size)
results <- vector("list", n_chunks)
for (chunk in seq_len(n_chunks)) {
  idx <- ((chunk - 1L) * PARAMS$chunk_size + 1L):min(chunk * PARAMS$chunk_size, n_pairs)
  ii <- pair_grid$i[idx]
  jj <- pair_grid$j[idx]

  r_r0 <- expr_up[ii, , drop = FALSE] - expr_down[jj, , drop = FALSE]
  r_r1 <- expr_up_r1[ii, , drop = FALSE] - expr_down_r1[jj, , drop = FALSE]

  # R0 sign consistency among non-dead-zone samples.
  valid_r0 <- abs(r_r0) >= PARAMS$dead_zone
  r0_pos <- rowSums((r_r0 > 0) & valid_r0)
  r0_neg <- rowSums((r_r0 < 0) & valid_r0)
  r0_sign <- ifelse(r0_pos >= r0_neg, 1, -1)
  pass_consistency <- pmin(r0_pos, r0_neg) <= PARAMS$r0_exception_max

  # R0 strength: 10th-percentile |r| clears the threshold.
  abs_r0 <- abs(r_r0)
  q10 <- apply(abs_r0, 1, function(x) sort(x, partial = q10_pos)[q10_pos])
  pass_q10 <- q10 >= PARAMS$r0_q10_threshold

  # R1 reversal: opposite sign to R0 among non-dead-zone; majority but not all.
  outside_r1 <- abs(r_r1) >= PARAMS$dead_zone
  reversal_count <- rowSums(outside_r1 & (sign(r_r1) != r0_sign))
  reversal_rate <- reversal_count / n_r1
  pass_reversal <- reversal_rate > 0.5 & reversal_count < n_r1

  keep <- pass_consistency & pass_q10 & pass_reversal
  if (!any(keep)) next
  k <- which(keep)
  results[[chunk]] <- data.frame(
    up = up_ids[ii[k]],
    down = down_ids[jj[k]],
    r0_sign = r0_sign[k],
    r0_q10 = q10[k],
    reversal_count = reversal_count[k],
    reversal_rate = reversal_rate[k],
    median_r0 = row_median(r_r0[k, , drop = FALSE]),
    median_r1 = row_median(r_r1[k, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
}

pairs <- do.call(rbind, results)
if (is.null(pairs) || nrow(pairs) == 0L) {
  stop("No candidate pair passed the reversal criteria.")
}
pairs$median_diff <- abs(pairs$median_r0 - pairs$median_r1)
pairs <- pairs[order(pairs$median_diff, decreasing = TRUE), , drop = FALSE]
rownames(pairs) <- NULL

message("Candidate pairs passing all criteria: ", nrow(pairs))
message(sprintf(
  "  median_diff range [%.3f, %.3f] ; reversal_rate range [%.2f, %.2f]",
  min(pairs$median_diff), max(pairs$median_diff),
  min(pairs$reversal_rate), max(pairs$reversal_rate)
))

# --- Assemble and save -----------------------------------------------------
name_of <- gene_name_map(se, strip_version = TRUE)
pairs$up_name <- unname(name_of[pairs$up])
pairs$down_name <- unname(name_of[pairs$down])

thyr_reo_candidate_pairs <- list(
  date = Sys.Date(),
  config = list(
    n_candidates = N_CANDIDATES,
    params = PARAMS,
    n_r0 = n_r0, n_r1 = n_r1,
    ranking_metric = "abs(effect - 0.5) from 410 BM",
    expression = "TPM from the selected-strand count assay"
  ),
  pairs = pairs
)
out_rds <- file.path(paths$processed, "thyr_reo_candidate_pairs.rds")
saveRDS(thyr_reo_candidate_pairs, out_rds)
message("Saved: ", out_rds, " (", nrow(pairs), " pairs)")

if (dir.exists(paths$output)) {
  out_csv <- file.path(paths$output, "reo_candidate_pairs.csv")
  utils::write.csv(pairs, out_csv, row.names = FALSE)
  message("Saved: ", out_csv)
}
