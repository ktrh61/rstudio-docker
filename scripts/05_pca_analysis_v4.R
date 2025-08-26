# ==============================================================================
# REBC-THYR PCA Analysis Script v4 - CDM (fast) on logCPM(TMM) + Pair Consistency
# 05_pca_analysis_v4.R
# ==============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(matrixStats)
  library(dplyr)
  library(Rcpp)
  library(RhpcBLASctl)
})

# ---- CDM core (as in v3) ----
source("./utils/with_openblas_threads.R")
sourceCpp("./utils/CDM_fast3_arma_enhanced.cpp")

cat("Starting PCA analysis v4: CDM on logCPM(TMM), protein-coding, MAD top-K...\n")

# ==============================================================================
# 0. Config
# ==============================================================================
PRIOR_COUNT <- 0.5
TOPK_GENES  <- 5000
VAR_DROP_FRAC <- 0.05   # drop bottom 5% variance before MAD ranking
MIN_SAMPLES <- 3
MAX_PC_OUT  <- 10

# ==============================================================================
# 1. CDM wrapper (unchanged I/O)
# ==============================================================================
CDM_fast_compatible <- function(X, by_sample = FALSE, center = TRUE, scale. = TRUE,
                                k = NULL, return_scores = TRUE, verbose = FALSE) {
  if (by_sample) {
    sample_names <- rownames(X)
    X <- t(X)  # genes x samples
  } else {
    sample_names <- colnames(X)
  }
  res <- with_openblas_threads("auto-2", {
    CDM_fast3_arma(X, verbose = verbose)
  })
  if (!is.null(res$scores) && !is.null(sample_names)) {
    rownames(res$scores) <- sample_names
    colnames(res$scores) <- paste0("PC", seq_len(ncol(res$scores)))
  }
  if (!is.null(k) && k < length(res$values)) {
    res$values  <- res$values[1:k]
    res$vectors <- res$vectors[, 1:k, drop = FALSE]
    if (!is.null(res$scores)) res$scores <- res$scores[, 1:k, drop = FALSE]
  }
  ve <- res$values^2 / sum(res$values^2)
  res$variance_explained  <- ve
  res$cumulative_variance <- cumsum(ve)
  res$n_components <- length(res$values)
  res
}

# ==============================================================================
# 2. Load Data & Sample Lists
# ==============================================================================
cat("Loading sample lists and SE counts...\n")
load("./data/processed/sample_lists.rda")  # -> sample_lists

if (!exists("se_thyr")) {
  load("./data/raw/thyr_data.rda")         # provides `data` or `se_thyr`
  if (exists("data")) { se_thyr <- data; rm(data) }
}
count_data_full <- assay(se_thyr, "stranded_second")
gene_info <- as.data.frame(rowData(se_thyr))

cat("Count matrix dims: ", paste(dim(count_data_full), collapse=" x "), "\n")
cat("Available groups: ", paste(names(sample_lists), collapse=", "), "\n")

# ==============================================================================
# 3. Helpers: biotype col, logCPM(TMM) + gene selection, scaling
# ==============================================================================
get_biotype_col <- function(df) {
  cands <- c("biotype","gene_type","type")
  hit <- cands[cands %in% colnames(df)]
  if (length(hit)) hit[1] else NA_character_
}

get_logcpm_for_pca <- function(count_mat, genes_df, samples,
                               K = TOPK_GENES, prior.count = PRIOR_COUNT,
                               var_drop_frac = VAR_DROP_FRAC) {
  smp <- samples[samples %in% colnames(count_mat)]
  if (length(smp) < MIN_SAMPLES) stop("Too few samples after intersect.")
  dge <- DGEList(counts = count_mat[, smp, drop=FALSE], genes = genes_df)
  dge <- calcNormFactors(dge, method = "TMM")
  logcpm <- cpm(dge, log=TRUE, prior.count=prior.count, normalized.lib.sizes=TRUE)
  
  biocol <- get_biotype_col(genes_df)
  if (!is.na(biocol)) {
    pc_mask <- genes_df[[biocol]] %in% c("protein_coding","protein-coding","protein coding")
    if (any(pc_mask, na.rm=TRUE)) logcpm <- logcpm[pc_mask, , drop=FALSE]
  }
  
  vars <- rowVars(logcpm)
  keep0 <- is.finite(vars) & vars > 0
  logcpm <- logcpm[keep0, , drop=FALSE]
  
  if (var_drop_frac > 0) {
    thr <- quantile(rowVars(logcpm), probs = var_drop_frac, na.rm = TRUE)
    logcpm <- logcpm[rowVars(logcpm) > thr, , drop=FALSE]
  }
  
  mads <- rowMads(logcpm)
  K_eff <- min(K, length(mads))
  logcpm[order(mads, decreasing=TRUE)[seq_len(K_eff)], , drop=FALSE]
}

prepare_for_cdm <- function(mat) {
  M <- scale(mat, center=TRUE, scale=TRUE)
  attr(M, "scaled:center") <- NULL
  attr(M, "scaled:scale")  <- NULL
  M
}

# ==============================================================================
# 4. Analysis groups (unchanged)
# ==============================================================================
analysis_groups <- list(
  "R0_normal" = list(group="R0", tissue="normal", description="RET Unexposed Normal"),
  "R0_tumor"  = list(group="R0", tissue="tumor",  description="RET Unexposed Tumor"),
  "R1_normal" = list(group="R1", tissue="normal", description="RET RadHigh Normal"),
  "R1_tumor"  = list(group="R1", tissue="tumor",  description="RET RadHigh Tumor"),
  "B0_normal" = list(group="B0", tissue="normal", description="BRAF Unexposed Normal"),
  "B0_tumor"  = list(group="B0", tissue="tumor",  description="BRAF Unexposed Tumor"),
  "B1_normal" = list(group="B1", tissue="normal", description="BRAF RadHigh Normal"),
  "B1_tumor"  = list(group="B1", tissue="tumor",  description="BRAF RadHigh Tumor")
)

# ==============================================================================
# 5. Run CDM per group on logCPM(TMM)
# ==============================================================================
cat("Running CDM on logCPM(TMM) per group...\n")
pca_results <- list()

for (analysis_name in names(analysis_groups)) {
  cat(sprintf("\n--- %s ---\n", analysis_name))
  gi <- analysis_groups[[analysis_name]]
  g  <- gi$group; tt <- gi$tissue; desc <- gi$description
  
  if (!g %in% names(sample_lists)) { cat("Group missing, skip\n"); next }
  
  sids <- if (tt=="tumor") sample_lists[[g]]$tumor else sample_lists[[g]]$normal
  sids <- sids[sids %in% colnames(count_data_full)]
  if (length(sids) < MIN_SAMPLES) {
    cat(sprintf("Insufficient samples (%d), skip\n", length(sids))); next
  }
  cat(sprintf("%s: %d samples\n", desc, length(sids)))
  
  # logCPM(TMM) + protein-coding + low-var drop + MAD topK
  group_logcpm <- tryCatch({
    get_logcpm_for_pca(count_data_full, gene_info, sids,
                       K = TOPK_GENES, prior.count = PRIOR_COUNT)
  }, error=function(e){
    cat("get_logcpm_for_pca failed: ", e$message, "\n"); return(NULL)
  })
  if (is.null(group_logcpm)) next
  cat(sprintf("logCPM genes selected: %d\n", nrow(group_logcpm)))
  
  X <- prepare_for_cdm(group_logcpm)
  if (!all(is.finite(X))) { cat("Non-finite values after scaling, skip\n"); next }
  
  pca_result <- tryCatch({
    CDM_fast_compatible(
      X = X, by_sample = FALSE, center = FALSE, scale. = FALSE,
      k = min(MAX_PC_OUT, ncol(X) - 1), return_scores = TRUE, verbose = FALSE
    )
  }, error=function(e){
    cat("CDM failed: ", e$message, "\n"); return(NULL)
  })
  if (is.null(pca_result)) next
  
  pca_results[[analysis_name]] <- list(
    pca = pca_result,
    samples = colnames(X),
    description = desc,
    n_genes = nrow(X),
    n_samples = ncol(X),
    analysis_name = analysis_name,
    group = g,
    tissue = tt
  )
  
  ve12 <- sum(pca_result$variance_explained[1:min(2, length(pca_result$variance_explained))]) * 100
  cat(sprintf("CDM done: %d comps, PC1-2 var≈ %.1f%%\n",
              pca_result$n_components, ve12))
}

# ==============================================================================
# 6. Outlier detection (as v3)
# ==============================================================================
enhanced_outlier_detection <- function(analysis_name, pca_data, verbose = TRUE) {
  if (verbose) cat(sprintf("\n--- Outlier detection: %s ---\n", analysis_name))
  scores <- pca_data$pca$scores
  if (is.null(scores) || nrow(scores) < 4) { if (verbose) cat("Too few samples\n"); return(NULL) }
  
  n_pcs <- min(3, ncol(scores))
  pcs   <- scores[, 1:n_pcs, drop=FALSE]
  ctr   <- colMeans(pcs)
  covm  <- cov(pcs)
  
  mahal_out <- integer(0); mahal_dist <- rep(NA_real_, nrow(scores))
  thr <- qchisq(0.95, df = n_pcs)
  if (det(covm) > .Machine$double.eps) {
    mahal_dist <- mahalanobis(pcs, ctr, covm)
    mahal_out  <- which(mahal_dist > thr)
  }
  
  pc1 <- scores[,1]
  q1 <- quantile(pc1, 0.25); q3 <- quantile(pc1, 0.75); iqr <- q3 - q1
  lf <- q1 - 1.5 * iqr; uf <- q3 + 1.5 * iqr
  iqr_out <- which(pc1 < lf | pc1 > uf)
  
  comb_idx <- unique(c(mahal_out, iqr_out))
  comb_smp <- rownames(scores)[comb_idx]
  
  if (verbose) {
    cat(sprintf("  Mahalanobis: %d | IQR: %d | Combined: %d\n",
                length(mahal_out), length(iqr_out), length(comb_idx)))
    if (length(comb_idx)) cat("  Outliers: ", paste(comb_smp, collapse=", "), "\n")
  }
  list(mahal_outliers = mahal_out,
       iqr_outliers   = iqr_out,
       combined_outliers = comb_idx,
       combined_samples  = comb_smp,
       total_samples = nrow(scores),
       analysis_name = analysis_name)
}

cat("\n==============================================\nOutlier Detection\n==============================================\n")
outlier_detection_results <- list()
for (nm in names(pca_results)) {
  od <- enhanced_outlier_detection(nm, pca_results[[nm]], verbose = TRUE)
  if (!is.null(od)) outlier_detection_results[[nm]] <- od
}

# ==============================================================================
# 7. Pair consistency filtering (unchanged)
# ==============================================================================
cat("\n==============================================\nPair Consistency Filtering\n==============================================\n")
pca_filtered_sample_lists <- sample_lists
outliers_by_group <- list()

for (nm in names(outlier_detection_results)) {
  od <- outlier_detection_results[[nm]]; if (is.null(od)) next
  if (!length(od$combined_samples)) next
  pr <- pca_results[[nm]]
  g  <- pr$group; tt <- pr$tissue
  if (!g %in% names(outliers_by_group)) outliers_by_group[[g]] <- list(tumor=character(0), normal=character(0))
  outliers_by_group[[g]][[tt]] <- od$combined_samples
}

for (g in c("R0","R1","B0","B1")) {
  if (!g %in% names(pca_filtered_sample_lists)) next
  orig_tumor  <- sample_lists[[g]]$tumor
  orig_normal <- sample_lists[[g]]$normal
  orig_cases  <- sample_lists[[g]]$cases
  if (!length(orig_cases)) { cat(g, ": no pairs\n"); next }
  
  idx <- integer(0)
  if (g %in% names(outliers_by_group)) {
    if (length(outliers_by_group[[g]]$tumor))
      idx <- c(idx, which(orig_tumor %in% outliers_by_group[[g]]$tumor))
    if (length(outliers_by_group[[g]]$normal))
      idx <- c(idx, which(orig_normal %in% outliers_by_group[[g]]$normal))
  }
  idx <- unique(idx)
  if (length(idx)) {
    exc_cases <- orig_cases[idx]
    cat(sprintf("%s: exclude %d pair(s): %s\n", g, length(idx), paste(exc_cases, collapse=", ")))
    keep <- setdiff(seq_along(orig_cases), idx)
    pca_filtered_sample_lists[[g]]$tumor  <- orig_tumor[keep]
    pca_filtered_sample_lists[[g]]$normal <- orig_normal[keep]
    pca_filtered_sample_lists[[g]]$cases  <- orig_cases[keep]
    cat(sprintf("%s: %d -> %d pairs\n", g, length(orig_cases), length(keep)))
  } else {
    cat(sprintf("%s: no outliers, keep all %d pairs\n", g, length(orig_cases)))
  }
}

# ==============================================================================
# 8. Summary table
# ==============================================================================
cat("\n==============================================\nSummary\n==============================================\n")
pca_summary <- data.frame(Group=character(), Tissue=character(),
                          Original_N=integer(), Outliers=integer(), Final_N=integer(),
                          PC1_Var=numeric(), PC2_Var=numeric(), stringsAsFactors=FALSE)

for (nm in names(analysis_groups)) {
  gi <- analysis_groups[[nm]]
  if (nm %in% names(pca_results) && !is.null(pca_results[[nm]])) {
    pr <- pca_results[[nm]]
    n_out <- if (nm %in% names(outlier_detection_results)) length(outlier_detection_results[[nm]]$combined_outliers) else 0
    pc1v <- pr$pca$variance_explained[1] * 100
    pc2v <- ifelse(length(pr$pca$variance_explained) > 1, pr$pca$variance_explained[2] * 100, 0)
    pca_summary <- rbind(pca_summary, data.frame(
      Group=gi$group, Tissue=gi$tissue,
      Original_N=pr$n_samples, Outliers=n_out, Final_N=pr$n_samples - n_out,
      PC1_Var=round(pc1v,1), PC2_Var=round(pc2v,1), stringsAsFactors=FALSE
    ))
  }
}
print(pca_summary)

cat("\nFinal paired sample counts after filtering:\n")
for (g in c("R0","R1","B0","B1")) {
  if (g %in% names(pca_filtered_sample_lists)) {
    cat(sprintf("  %s: %d pairs\n", g, length(pca_filtered_sample_lists[[g]]$cases)))
  }
}

# ==============================================================================
# 9. Save
# ==============================================================================
cat("\nSaving v4 results...\n")
if (!dir.exists("./data/processed")) dir.create("./data/processed", recursive = TRUE)

pca_phase1_results_v4 <- list(
  pca_results = pca_results,
  outlier_detection_results = outlier_detection_results,
  pca_filtered_sample_lists = pca_filtered_sample_lists,
  pca_summary = pca_summary,
  analysis_groups = analysis_groups,
  analysis_date = Sys.time(),
  analysis_version = "v4_cdm_logcpm_TMM_topK",
  cdm_implementation = "CDM_fast3_arma_enhanced",
  config = list(prior.count=PRIOR_COUNT, topK=TOPK_GENES, var_drop_frac=VAR_DROP_FRAC)
)

save(pca_phase1_results_v4, file = "./data/processed/pca_phase1_results_v4.rda")
save(pca_filtered_sample_lists, file = "./data/processed/pca_filtered_sample_lists_v4.rda")
cat("Saved to ./data/processed/ (v4)\n")

# ==============================================================================
# 10. Final notes
# ==============================================================================
tot_out <- sum(pca_summary$Outliers, na.rm=TRUE)
cat("\n==============================================\n")
cat("PCA Phase 1 v4 Complete\n")
cat("==============================================\n")
cat(sprintf("Total outliers detected: %d\n", tot_out))
viable <- character(0)
for (g in c("R0","R1","B0","B1")) {
  if (g %in% names(pca_filtered_sample_lists)) {
    n_pairs <- length(pca_filtered_sample_lists[[g]]$cases)
    if (n_pairs >= 5) viable <- c(viable, g)
  }
}
cat(sprintf("\nGroups viable for next steps (≥%d pairs): %s\n", 5, paste(viable, collapse=", ")))
cat("\nNext: 06_contamde / 08_deges_normalization_v3 downstream as-is.\n")
