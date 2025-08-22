# ==============================================================================
# REBC-THYR | REO-based feature selection (DEG up x down, TPM, zero-free genes)
# 11_feature_selection_reo_clean_impl.R
# ==============================================================================
# Goal:
#  - Use DEG up/down lists
#  - Filter to genes with NO zero TPM across analysis samples
#  - REO with r = log2(p/q) (no pseudocount)
#  - R0: all (or all-1) samples one-direction (p>q or p<q)
#  - R1: median opposite sign; reverse-fraction within [alpha_low, alpha_high]
#  - Strength: R0 side q10 >= log2(1.5) (tunable)
#  - Efficient block processing; optional parallel over up-gene blocks
# ==============================================================================

suppressPackageStartupMessages({
  library(matrixStats)      # fast row/col stats
  library(future.apply)     # optional parallel (falls back to sequential if not planned)
})

cat("Starting REO clean implementation...\n")

# ------------------------------
# Parameters (tunable defaults)
# ------------------------------
allow_one_exception_R0 <- TRUE      # TRUE = "all or all-1" in R0
alpha_low  <- 0.60                  # required lower bound of reverse fraction in R1
alpha_high <- 0.95                  # upper bound to avoid 100% reverse in R1
strength_fc <- log2(1.5)            # ~0.585; required strength on R0 side
strength_mode <- "q10"              # "q10" or "min"
block_up_size <- 64                 # number of up-genes per block (tune to memory/CPU)
use_parallel <- TRUE                # require future::plan() prior to run if TRUE

# ------------------------------
# Load inputs
# ------------------------------
if (!file.exists("./data/processed/final_deg_v3_results.rda")) {
  stop("DEG results not found: ./data/processed/final_deg_v3_results.rda")
}
load("./data/processed/final_deg_v3_results.rda")

if (!file.exists("./data/processed/thyr_tpm.rda")) {
  stop("TPM data not found: ./data/processed/thyr_tpm.rda")
}
load("./data/processed/thyr_tpm.rda")   # expect object `tpm`

if (!file.exists("./data/processed/final_high_purity_sample_lists.rda")) {
  stop("High-purity sample lists not found.")
}
load("./data/processed/final_high_purity_sample_lists.rda")
high_purity_sample_lists <- final_high_purity_sample_lists

# ------------------------------
# DEG: use R0_vs_R1_tumor
# ------------------------------
if (!"R0_vs_R1_tumor" %in% names(final_deg_v3_results$deg_analysis_results)) {
  stop("R0_vs_R1_tumor not found in DEG results")
}
tumor_deg_result <- final_deg_v3_results$deg_analysis_results[["R0_vs_R1_tumor"]]
tumor_degs <- tumor_deg_result$deg_summary$results_df

up_ids   <- tumor_degs$gene_id[tumor_degs$significance_category == "Upregulated"]
down_ids <- tumor_degs$gene_id[tumor_degs$significance_category == "Downregulated"]

# Clean ENSG.x style IDs to match TPM rownames if needed
clean_id <- function(x) sub("\\.\\d+$", "", x)
tpm_clean <- clean_id(rownames(tpm))
up_clean <- clean_id(up_ids)
down_clean <- clean_id(down_ids)

# Map by first match (and warn on duplicates)
map_ids <- function(clean_vec, ref_clean, ref_names) {
  m <- match(clean_vec, ref_clean)
  idx <- which(!is.na(m))
  res <- ref_names[m[idx]]
  # Drop duplicates (same TPM row matched multiple times)
  res <- unique(res)
  res
}
up_genes   <- map_ids(up_clean, tpm_clean, rownames(tpm))
down_genes <- map_ids(down_clean, tpm_clean, rownames(tpm))

cat(sprintf("DEG available in TPM: up=%d, down=%d\n", length(up_genes), length(down_genes)))

# ------------------------------
# Samples and analysis TPM
# ------------------------------
r0_samples <- high_purity_sample_lists$R0$tumor
r1_samples <- high_purity_sample_lists$R1$tumor

r0_avail <- intersect(r0_samples, colnames(tpm))
r1_avail <- intersect(r1_samples, colnames(tpm))
stopifnot(length(r0_avail) > 0, length(r1_avail) > 0)

analysis_samples <- c(r0_avail, r1_avail)
analysis_tpm <- tpm[, analysis_samples, drop = FALSE]
analysis_tpm <- as.matrix(analysis_tpm)
storage.mode(analysis_tpm) <- "double"

idxR0 <- match(r0_avail, analysis_samples)
idxR1 <- match(r1_avail, analysis_samples)

cat(sprintf("Samples: R0=%d, R1=%d, total=%d\n",
            length(idxR0), length(idxR1), ncol(analysis_tpm)))

# ------------------------------
# Zero-free gene filter (strict)
# ------------------------------
keep <- rowSums(analysis_tpm == 0) == 0
analysis_tpm <- analysis_tpm[keep, , drop = FALSE]
cat(sprintf("Zero-free genes: %d / %d\n", nrow(analysis_tpm), length(keep)))

# remap up/down to zero-free set
up_genes   <- intersect(up_genes, rownames(analysis_tpm))
down_genes <- intersect(down_genes, rownames(analysis_tpm))
stopifnot(length(up_genes) > 0, length(down_genes) > 0)

# ------------------------------
# Log2 transform (no pseudocount)
# ------------------------------
y <- log2(analysis_tpm)               # gene x sample
yA <- y[, idxR0, drop = FALSE]        # R0
yB <- y[, idxR1, drop = FALSE]        # R1

# ------------------------------
# Helpers for per-up processing
# ------------------------------
require_strength <- function(rA_row) {
  if (strength_mode == "min") {
    min(rA_row) >= strength_fc
  } else { # q10
    # robust: 10th percentile across R0 samples
    as.numeric(quantile(rA_row, probs = 0.10, type = 7)) >= strength_fc
  }
}

one_direction_R0 <- function(rA_row) {
  # returns list(dir = "pos"/"neg"/NA, frac = fraction of majority sign)
  pos_frac <- mean(rA_row > 0)
  neg_frac <- mean(rA_row < 0)
  thr <- if (allow_one_exception_R0) 1 - 1/length(rA_row) else 1.0
  if (pos_frac >= thr) {
    list(dir = "pos", frac = pos_frac)
  } else if (neg_frac >= thr) {
    list(dir = "neg", frac = neg_frac)
  } else {
    list(dir = NA_character_, frac = max(pos_frac, neg_frac))
  }
}

passes_R1_mixed <- function(rB_row, dir) {
  # require median opposite sign AND reverse fraction in [alpha_low, alpha_high]
  med <- median(rB_row)
  if (dir == "pos") {
    frac_rev <- mean(rB_row < 0)
    median_ok <- (med < 0)
  } else {
    frac_rev <- mean(rB_row > 0)
    median_ok <- (med > 0)
  }
  list(ok = (median_ok && frac_rev >= alpha_low && frac_rev <= alpha_high),
       frac_rev = frac_rev,
       med = med)
}

# ------------------------------
# Blocked enumeration up x down
# ------------------------------
U <- length(up_genes); D <- length(down_genes)
cat(sprintf("Enumerating pairs: up=%d x down=%d = %s\n",
            U, D, format(U*D, big.mark=",")))

up_blocks <- split(seq_len(U), ceiling(seq_len(U) / block_up_size))

process_block <- function(block_idx) {
  ub <- up_blocks[[block_idx]]
  up_ids_block <- up_genes[ub]
  
  # pre-extract down matrices once per block
  YdA <- yA[down_genes, , drop = FALSE]  # D x nA
  YdB <- yB[down_genes, , drop = FALSE]  # D x nB
  
  out_list <- vector("list", length(up_ids_block))
  names(out_list) <- up_ids_block
  
  for (u_i in seq_along(up_ids_block)) {
    upg <- up_ids_block[u_i]
    yuA <- yA[upg, , drop = TRUE]        # length nA
    yuB <- yB[upg, , drop = TRUE]        # length nB
    
    # rA: (p - q) for all downs vs this up (D x nA)
    # rB: (p - q) for all downs vs this up (D x nB)
    rA <- matrix(yuA, nrow = D, ncol = length(yuA), byrow = TRUE) - YdA
    rB <- matrix(yuB, nrow = D, ncol = length(yuB), byrow = TRUE) - YdB
    
    # R0 direction (all or all-1)
    dir_info <- apply(rA, 1, one_direction_R0)
    dir_vec  <- vapply(dir_info, function(z) z$dir, character(1))
    fracR0   <- vapply(dir_info, function(z) z$frac, numeric(1))
    
    keep_rows <- which(!is.na(dir_vec))
    if (length(keep_rows) == 0L) {
      out_list[[u_i]] <- NULL
      next
    }
    
    rA_keep <- rA[keep_rows, , drop = FALSE]
    rB_keep <- rB[keep_rows, , drop = FALSE]
    dir_keep <- dir_vec[keep_rows]
    fracR0_keep <- fracR0[keep_rows]
    down_keep <- down_genes[keep_rows]
    
    # Strength on R0 side
    if (strength_mode == "min") {
      strong_ok <- rowMins(rA_keep) >= strength_fc
      strength_metric <- rowMins(rA_keep)
    } else {
      strength_metric <- rowQuantiles(rA_keep, probs = 0.10)
      strong_ok <- strength_metric >= strength_fc
    }
    
    if (!any(strong_ok)) {
      out_list[[u_i]] <- NULL
      next
    }
    
    # R1 mixed + median opposite
    # Row-wise loop (vectorized enough since filtered)
    idx_ok <- logical(length(keep_rows))
    frac_rev <- numeric(length(keep_rows))
    medR0 <- rowMedians(rA_keep)
    medR1 <- rowMedians(rB_keep)
    
    for (k in seq_along(keep_rows)) {
      di <- dir_keep[k]
      res <- passes_R1_mixed(rB_keep[k, ], di)
      idx_ok[k] <- res$ok
      frac_rev[k] <- res$frac_rev
    }
    
    final_idx <- which(strong_ok & idx_ok)
    if (length(final_idx) == 0L) {
      out_list[[u_i]] <- NULL
      next
    }
    
    df <- data.frame(
      up_gene = upg,
      down_gene = down_keep[final_idx],
      r0_direction = dir_keep[final_idx],                   # "pos" means R0: p>q
      r0_frac_major = round(fracR0_keep[final_idx], 4),
      r0_strength = round(strength_metric[final_idx], 4),   # q10 or min of rA
      r0_median = round(medR0[final_idx], 4),
      r1_median = round(medR1[final_idx], 4),
      r1_frac_reverse = round(frac_rev[final_idx], 4),
      stringsAsFactors = FALSE
    )
    out_list[[u_i]] <- df
  }
  
  do.call(rbind, out_list)
}

cat("Processing blocks...\n")
t_start <- proc.time()[3]

if (use_parallel) {
  # If user hasn't set plan(), provide a reasonable default
  if (is.null(future::plan("list")$strategy)) {
    os <- tolower(Sys.info()[["sysname"]])
    if (os == "windows") {
      future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1))
    } else {
      future::plan(future::multicore, workers = max(1, parallel::detectCores() - 1))
    }
  }
  res_list <- future_lapply(seq_along(up_blocks), process_block, future.seed = TRUE)
} else {
  res_list <- lapply(seq_along(up_blocks), process_block)
}

candidates <- do.call(rbind, res_list)
elapsed <- proc.time()[3] - t_start
cat(sprintf("Done. Time: %.1f sec. Candidates: %s\n",
            elapsed, ifelse(is.null(candidates), 0, nrow(candidates))))

if (is.null(candidates) || nrow(candidates) == 0L) {
  cat("No candidates satisfied the criteria. Consider relaxing:\n",
      "- strength_fc (currently log2(1.5))\n",
      "- alpha_low/high (currently 0.60–0.95)\n",
      "- allow_one_exception_R0 (TRUE→FALSE)\n", sep = "")
} else {
  # Sort by: stronger R0 strength, larger median gap, then R1 reverse fraction
  candidates$median_gap <- round(candidates$r0_median - candidates$r1_median, 4)
  candidates <- candidates[order(-candidates$r0_strength, -candidates$median_gap,
                                 -candidates$r1_frac_reverse), ]
}

# ------------------------------
# Save outputs
# ------------------------------
dir.create("./data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("./output/reports", showWarnings = FALSE, recursive = TRUE)

saveRDS(list(
  params = list(
    allow_one_exception_R0 = allow_one_exception_R0,
    alpha_low = alpha_low, alpha_high = alpha_high,
    strength_fc = strength_fc, strength_mode = strength_mode,
    block_up_size = block_up_size
  ),
  samples = list(R0 = r0_avail, R1 = r1_avail),
  up_genes = up_genes, down_genes = down_genes,
  zero_free_genes = rownames(analysis_tpm),
  candidates = candidates
), file = "./data/processed/reo_deg_ud_clean_results.rds")

if (!is.null(candidates) && nrow(candidates) > 0L) {
  write.csv(candidates, "./output/reports/reo_deg_ud_candidates.csv", row.names = FALSE)
  cat("Candidates written to ./output/reports/reo_deg_ud_candidates.csv\n")
}

cat("All done.\n")
