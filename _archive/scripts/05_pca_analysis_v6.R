# ==============================================================================
# REBC-THYR PCA Analysis Script v6 - Simple logCPM + CDM
# 05_pca_analysis_v6.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(dplyr)
library(edgeR)

cat("Starting PCA analysis v6 with simple logCPM + CDM...\n")

# ==============================================================================
# 1. Load Data and Sample Lists
# ==============================================================================

cat("Loading data and sample lists...\n")

load("./data/processed/sample_lists.rda")
load("./data/raw/thyr_data.rda")
se_thyr <- data
count_data <- assay(se_thyr, "stranded_second")
rm(data, se_thyr)

cat("Count data dimensions:", dim(count_data), "\n")
cat("Available groups:", names(sample_lists), "\n")

# ==============================================================================
# 2. Gene Filtering
# ==============================================================================

cat("Applying gene filtering...\n")

dgelist_all <- DGEList(counts = count_data)

all_samples <- unique(unlist(lapply(sample_lists, function(x) c(x$tumor, x$normal))))
available_samples <- intersect(all_samples, colnames(count_data))

if (length(available_samples) < 10) {
  stop("Insufficient samples available in count data")
}

dgelist_filtered <- dgelist_all[, available_samples]
keep <- filterByExpr(dgelist_filtered)
dgelist_filtered <- dgelist_filtered[keep, , keep.lib.sizes = FALSE]

cat(sprintf("After filterByExpr: %d genes (from %d)\n", 
            nrow(dgelist_filtered), nrow(count_data)))

filtered_counts <- dgelist_filtered$counts

# ==============================================================================
# 3. Define Analysis Groups
# ==============================================================================

analysis_groups <- list(
  "R0_normal" = list(group = "R0", tissue = "normal", description = "RET Unexposed Normal"),
  "R0_tumor" = list(group = "R0", tissue = "tumor", description = "RET Unexposed Tumor"),
  "R1_normal" = list(group = "R1", tissue = "normal", description = "RET RadHigh Normal"),
  "R1_tumor" = list(group = "R1", tissue = "tumor", description = "RET RadHigh Tumor"),
  "B0_normal" = list(group = "B0", tissue = "normal", description = "BRAF Unexposed Normal"),
  "B0_tumor" = list(group = "B0", tissue = "tumor", description = "BRAF Unexposed Tumor"),
  "B1_normal" = list(group = "B1", tissue = "normal", description = "BRAF RadHigh Normal"),
  "B1_tumor" = list(group = "B1", tissue = "tumor", description = "BRAF RadHigh Tumor")
)

# ==============================================================================
# 4. Simple logCPM Conversion Function
# ==============================================================================

prepare_cdm_data_simple <- function(count_matrix, sample_ids, analysis_name) {
  cat(sprintf("\n--- Preparing CDM data for %s ---\n", analysis_name))
  
  analysis_counts <- count_matrix[, sample_ids, drop = FALSE]
  cat(sprintf("Analysis counts: %d genes x %d samples\n", 
              nrow(analysis_counts), ncol(analysis_counts)))
  
  # Simple CPM calculation
  cat("Converting to CPM...\n")
  library_sizes <- colSums(analysis_counts)
  cpm_data <- sweep(analysis_counts, 2, library_sizes, "/") * 1e6
  
  # logCPM transformation
  cat("Converting to logCPM...\n")
  prior_count <- 0.5
  log_cpm <- log2(cpm_data + prior_count)
  
  # Additional low expression filtering within group
  # Remove genes with very low expression (mean logCPM < 1)
  mean_logcpm <- rowMeans(log_cpm)
  expressed_genes <- mean_logcpm >= 1
  
  if (sum(expressed_genes) < nrow(log_cpm)) {
    log_cpm <- log_cpm[expressed_genes, , drop = FALSE]
    cat(sprintf("Removed %d low-expression genes (mean logCPM < 1)\n", 
                sum(!expressed_genes)))
  }
  
  cat(sprintf("Final logCPM data: %d genes x %d samples\n", 
              nrow(log_cpm), ncol(log_cpm)))
  
  return(list(
    log_cpm = log_cpm,
    sample_ids = sample_ids,
    n_genes_final = nrow(log_cpm),
    library_sizes = library_sizes
  ))
}

# ==============================================================================
# 5. CDM-based PCA Function (unchanged)
# ==============================================================================

perform_cdm_pca <- function(log_cpm_data, analysis_name, verbose = TRUE) {
  cat(sprintf("Performing CDM-based PCA for %s...\n", analysis_name))
  
  X <- log_cpm_data
  n_genes <- nrow(X)
  n_samples <- ncol(X)
  
  if (verbose) {
    cat(sprintf("CDM input: %d genes x %d samples\n", n_genes, n_samples))
  }
  
  X_centered <- X - rowMeans(X)
  gene_sds <- apply(X_centered, 1, sd)
  
  non_zero_var <- gene_sds > .Machine$double.eps
  if (sum(non_zero_var) < nrow(X)) {
    cat(sprintf("Removing %d zero-variance genes\n", sum(!non_zero_var)))
    X_centered <- X_centered[non_zero_var, , drop = FALSE]
    gene_sds <- gene_sds[non_zero_var]
  }
  
  X_scaled <- X_centered / gene_sds
  
  # CDM approach
  n1 <- floor(ncol(X_scaled) / 2)
  n2 <- ncol(X_scaled) - n1
  
  X1 <- X_scaled[, 1:n1, drop = FALSE]
  X2 <- X_scaled[, (n1+1):ncol(X_scaled), drop = FALSE]
  
  C <- crossprod(X1, X2) / sqrt((n1-1) * (n2-1))
  
  svd_c <- svd(C)
  
  eps <- .Machine$double.eps
  d_safe <- pmax(svd_c$d, eps)
  
  W1 <- X1 %*% sweep(svd_c$u, 2, d_safe, "/")
  W2 <- X2 %*% sweep(svd_c$v, 2, d_safe, "/")
  loadings <- (W1 + W2) / 2
  
  qr_result <- qr(loadings)
  loadings <- qr.Q(qr_result)
  
  scores <- t(loadings) %*% X_scaled
  scores <- t(scores)
  
  eigenvalues <- svd_c$d^2
  total_var <- sum(eigenvalues)
  variance_explained <- eigenvalues / total_var
  cumulative_variance <- cumsum(variance_explained)
  
  max_components <- min(10, ncol(X_scaled) - 1, nrow(X_scaled))
  if (ncol(scores) > max_components) {
    scores <- scores[, 1:max_components, drop = FALSE]
    loadings <- loadings[, 1:max_components, drop = FALSE]
    variance_explained <- variance_explained[1:max_components]
    cumulative_variance <- cumulative_variance[1:max_components]
    eigenvalues <- eigenvalues[1:max_components]
  }
  
  rownames(scores) <- colnames(X)
  colnames(scores) <- paste0("PC", 1:ncol(scores))
  
  if (verbose) {
    cat(sprintf("PCA completed: %d components\n", ncol(scores)))
    cat(sprintf("PC1-2 explain %.1f%% of variance\n", 
                sum(variance_explained[1:min(2, length(variance_explained))]) * 100))
  }
  
  return(list(
    scores = scores,
    loadings = loadings,
    values = eigenvalues,
    variance_explained = variance_explained,
    cumulative_variance = cumulative_variance,
    n_components = ncol(scores)
  ))
}

# ==============================================================================
# 6. Outlier Detection Function (unchanged)
# ==============================================================================

enhanced_outlier_detection <- function(analysis_name, pca_result, verbose = TRUE) {
  if (verbose) cat(sprintf("\n--- Outlier detection for %s ---\n", analysis_name))
  
  scores <- pca_result$scores
  if (is.null(scores) || nrow(scores) < 4) {
    if (verbose) cat("Insufficient data for outlier detection\n")
    return(NULL)
  }
  
  n_pcs_for_outlier <- min(3, ncol(scores))
  pc_subset <- scores[, 1:n_pcs_for_outlier, drop = FALSE]
  
  center_point <- colMeans(pc_subset)
  cov_matrix <- cov(pc_subset)
  
  mahal_outliers <- integer(0)
  mahal_dist <- rep(NA, nrow(scores))
  threshold <- qchisq(0.95, df = n_pcs_for_outlier)
  
  if (det(cov_matrix) > .Machine$double.eps) {
    mahal_dist <- mahalanobis(pc_subset, center_point, cov_matrix)
    mahal_outliers <- which(mahal_dist > threshold)
  }
  
  pc1_values <- scores[, 1]
  q1 <- quantile(pc1_values, 0.25)
  q3 <- quantile(pc1_values, 0.75)
  iqr <- q3 - q1
  
  lower_fence <- q1 - 1.5 * iqr
  upper_fence <- q3 + 1.5 * iqr
  
  iqr_outliers <- which(pc1_values < lower_fence | pc1_values > upper_fence)
  
  combined_outliers <- unique(c(mahal_outliers, iqr_outliers))
  combined_samples <- rownames(scores)[combined_outliers]
  
  if (verbose) {
    cat(sprintf("  Mahalanobis outliers: %d\n", length(mahal_outliers)))
    cat(sprintf("  IQR outliers: %d\n", length(iqr_outliers)))
    cat(sprintf("  Combined outliers: %d\n", length(combined_outliers)))
    if (length(combined_outliers) > 0) {
      cat(sprintf("  Outlier samples: %s\n", paste(combined_samples, collapse = ", ")))
    }
  }
  
  return(list(
    mahal_outliers = mahal_outliers,
    iqr_outliers = iqr_outliers,
    combined_outliers = combined_outliers,
    combined_samples = combined_samples,
    total_samples = nrow(scores),
    analysis_name = analysis_name
  ))
}

# ==============================================================================
# 7. Main Analysis Loop
# ==============================================================================

cat("\nPerforming PCA analysis with simple logCPM + CDM for each group...\n")

pca_results_v6 <- list()

for (analysis_name in names(analysis_groups)) {
  cat(sprintf("\n=== Analyzing %s ===\n", analysis_name))
  
  group_info <- analysis_groups[[analysis_name]]
  group_name <- group_info$group
  tissue_type <- group_info$tissue
  description <- group_info$description
  
  if (!group_name %in% names(sample_lists)) {
    cat(sprintf("Group %s not found, skipping...\n", group_name))
    next
  }
  
  if (tissue_type == "tumor") {
    sample_ids <- sample_lists[[group_name]]$tumor
  } else {
    sample_ids <- sample_lists[[group_name]]$normal
  }
  
  if (length(sample_ids) == 0) {
    cat(sprintf("No samples for %s, skipping...\n", analysis_name))
    next
  }
  
  available_samples <- sample_ids[sample_ids %in% colnames(filtered_counts)]
  
  if (length(available_samples) < 3) {
    cat(sprintf("Insufficient samples (%d) for %s, skipping...\n", 
                length(available_samples), analysis_name))
    next
  }
  
  cat(sprintf("%s: %d samples\n", description, length(available_samples)))
  
  tryCatch({
    cdm_data <- prepare_cdm_data_simple(filtered_counts, available_samples, analysis_name)
    pca_result <- perform_cdm_pca(cdm_data$log_cpm, analysis_name, verbose = TRUE)
    
    pca_results_v6[[analysis_name]] <- list(
      pca = pca_result,
      samples = available_samples,
      description = description,
      n_genes = cdm_data$n_genes_final,
      n_samples = length(available_samples),
      analysis_name = analysis_name,
      group = group_name,
      tissue = tissue_type,
      library_sizes = cdm_data$library_sizes
    )
    
    cat(sprintf("PCA completed: %d components, %.1f%% variance explained by PC1-2\n",
                pca_result$n_components,
                sum(pca_result$variance_explained[1:min(2, length(pca_result$variance_explained))]) * 100))
    
  }, error = function(e) {
    cat(sprintf("Error in analysis for %s: %s\n", analysis_name, e$message))
    pca_results_v6[[analysis_name]] <- NULL
  })
}

# ==============================================================================
# 8. Outlier Detection
# ==============================================================================

cat("\n=== Outlier Detection Analysis ===\n")

outlier_detection_results_v6 <- list()

for (analysis_name in names(pca_results_v6)) {
  if (!is.null(pca_results_v6[[analysis_name]])) {
    outlier_result <- enhanced_outlier_detection(
      analysis_name, 
      pca_results_v6[[analysis_name]]$pca, 
      verbose = TRUE
    )
    
    if (!is.null(outlier_result)) {
      outlier_detection_results_v6[[analysis_name]] <- outlier_result
    }
  }
}

# ==============================================================================
# 9. Create Filtered Sample Lists with Pair Consistency
# ==============================================================================

cat("\n=== Creating filtered sample lists with pair consistency ===\n")

pca_filtered_sample_lists_v6 <- sample_lists

outliers_by_group <- list()

for (analysis_name in names(outlier_detection_results_v6)) {
  if (is.null(outlier_detection_results_v6[[analysis_name]])) next
  
  outlier_data <- outlier_detection_results_v6[[analysis_name]]
  if (length(outlier_data$combined_samples) == 0) next
  
  pca_data <- pca_results_v6[[analysis_name]]
  group_name <- pca_data$group
  tissue_type <- pca_data$tissue
  
  if (!group_name %in% names(outliers_by_group)) {
    outliers_by_group[[group_name]] <- list(tumor = character(0), normal = character(0))
  }
  
  outliers_by_group[[group_name]][[tissue_type]] <- outlier_data$combined_samples
}

for (group_name in c("R0", "R1", "B0", "B1")) {
  if (!group_name %in% names(pca_filtered_sample_lists_v6)) next
  
  orig_tumor <- sample_lists[[group_name]]$tumor
  orig_normal <- sample_lists[[group_name]]$normal
  orig_cases <- sample_lists[[group_name]]$cases
  
  if (length(orig_cases) == 0 || length(orig_tumor) == 0 || length(orig_normal) == 0) {
    cat(sprintf("%s: No paired samples to process\n", group_name))
    next
  }
  
  outlier_indices <- c()
  
  if (group_name %in% names(outliers_by_group)) {
    if (length(outliers_by_group[[group_name]]$tumor) > 0) {
      tumor_outlier_idx <- which(orig_tumor %in% outliers_by_group[[group_name]]$tumor)
      outlier_indices <- c(outlier_indices, tumor_outlier_idx)
      cat(sprintf("%s: Found %d tumor outlier(s)\n", group_name, length(tumor_outlier_idx)))
    }
    
    if (length(outliers_by_group[[group_name]]$normal) > 0) {
      normal_outlier_idx <- which(orig_normal %in% outliers_by_group[[group_name]]$normal)
      outlier_indices <- c(outlier_indices, normal_outlier_idx)
      cat(sprintf("%s: Found %d normal outlier(s)\n", group_name, length(normal_outlier_idx)))
    }
  }
  
  outlier_indices <- unique(outlier_indices)
  
  if (length(outlier_indices) > 0) {
    exclude_cases <- orig_cases[outlier_indices]
    cat(sprintf("%s: Excluding cases for pair consistency: %s\n", 
                group_name, paste(exclude_cases, collapse = ", ")))
    
    keep_indices <- setdiff(1:length(orig_cases), outlier_indices)
    
    pca_filtered_sample_lists_v6[[group_name]]$tumor <- orig_tumor[keep_indices]
    pca_filtered_sample_lists_v6[[group_name]]$normal <- orig_normal[keep_indices]
    pca_filtered_sample_lists_v6[[group_name]]$cases <- orig_cases[keep_indices]
    
    cat(sprintf("%s: %d pairs → %d pairs after filtering\n",
                group_name, length(orig_cases), length(keep_indices)))
  } else {
    cat(sprintf("%s: No outliers detected, keeping all %d pairs\n", 
                group_name, length(orig_cases)))
  }
}

# ==============================================================================
# 10. Summary Statistics
# ==============================================================================

cat("\n=== PCA Phase 1 Summary (Simple logCPM + CDM) ===\n")

writeLines(capture.output(sessionInfo()), "session_info_pca_v6.txt")

pca_summary_v6 <- data.frame(
  Group = character(0),
  Tissue = character(0),
  Original_N = integer(0),
  Outliers = integer(0),
  Final_N = integer(0),
  PC1_Var = numeric(0),
  PC2_Var = numeric(0),
  Genes_Used = integer(0),
  stringsAsFactors = FALSE
)

for (analysis_name in names(analysis_groups)) {
  group_info <- analysis_groups[[analysis_name]]
  
  if (analysis_name %in% names(pca_results_v6) && !is.null(pca_results_v6[[analysis_name]])) {
    pca_data <- pca_results_v6[[analysis_name]]
    
    n_outliers <- 0
    if (analysis_name %in% names(outlier_detection_results_v6)) {
      n_outliers <- length(outlier_detection_results_v6[[analysis_name]]$combined_outliers)
    }
    
    pc1_var <- pca_data$pca$variance_explained[1] * 100
    pc2_var <- ifelse(length(pca_data$pca$variance_explained) > 1, 
                      pca_data$pca$variance_explained[2] * 100, 0)
    
    pca_summary_v6 <- rbind(pca_summary_v6, data.frame(
      Group = group_info$group,
      Tissue = group_info$tissue,
      Original_N = pca_data$n_samples,
      Outliers = n_outliers,
      Final_N = pca_data$n_samples - n_outliers,
      PC1_Var = round(pc1_var, 1),
      PC2_Var = round(pc2_var, 1),
      Genes_Used = pca_data$n_genes,
      stringsAsFactors = FALSE
    ))
  }
}

print(pca_summary_v6)

cat("\nFinal paired sample counts after outlier removal:\n")
for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(pca_filtered_sample_lists_v6)) {
    n_pairs <- length(pca_filtered_sample_lists_v6[[group]]$cases)
    cat(sprintf("  %s: %d pairs\n", group, n_pairs))
  }
}

# ==============================================================================
# 11. Save Results
# ==============================================================================

cat("\nSaving PCA Phase 1 results v6...\n")

pca_phase1_results_v6 <- list(
  pca_results = pca_results_v6,
  outlier_detection_results = outlier_detection_results_v6,
  pca_filtered_sample_lists = pca_filtered_sample_lists_v6,
  pca_summary = pca_summary_v6,
  analysis_groups = analysis_groups,
  analysis_date = Sys.time(),
  analysis_version = "v6_simple_logcpm_cdm"
)

save(pca_phase1_results_v6, file = "./data/processed/pca_phase1_results_v6.rda")
save(pca_filtered_sample_lists_v6, file = "./data/processed/pca_filtered_sample_lists_v6.rda")

cat("Results saved to ./data/processed/\n")

# ==============================================================================
# 12. Final Summary
# ==============================================================================

cat("\n=== PCA Phase 1 Analysis v6 Complete ===\n")

total_outliers <- sum(pca_summary_v6$Outliers)
cat(sprintf("Total outliers detected: %d\n", total_outliers))

min_samples <- 5
viable_groups <- character(0)
for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(pca_filtered_sample_lists_v6)) {
    n_pairs <- length(pca_filtered_sample_lists_v6[[group]]$cases)
    if (n_pairs >= min_samples) {
      viable_groups <- c(viable_groups, group)
    }
  }
}

cat(sprintf("\nGroups viable for downstream analysis (≥%d pairs): %s\n",
            min_samples, paste(viable_groups, collapse = ", ")))

cat("\nKey features of v6:\n")
cat("1. Simple CPM + logCPM transformation\n")
cat("2. Within-group low expression filtering (mean logCPM >= 1)\n")
cat("3. No complex normalization methods\n")
cat("4. Simplified and faster processing\n")

cat("\nPCA Phase 1 v6 completed successfully!\n")