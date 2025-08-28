# ==============================================================================
# REBC-THYR PCA Phase 2 Analysis v6 - Group Comparison with logCPM
# 07_pca_phase2_v6.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(dplyr)
library(cluster)
library(edgeR)

# Source CDM functions from v6
source("./utils/CDM_fast.R")

cat("Starting PCA Phase 2 v6: Group comparison with logCPM normalization...\n")
cat("Using high-purity samples from 06_purity_analysis_v6.R\n")

# ==============================================================================
# 1. Helper Functions for logCPM and CDM
# ==============================================================================

# Simple logCPM transformation (exactly as in 05_pca_analysis_v6)
prepare_logcpm_data <- function(count_matrix, sample_ids, analysis_name) {
  
  cat(sprintf("Preparing logCPM data for %s...\n", analysis_name))
  
  # Extract samples
  sample_data <- count_matrix[, sample_ids, drop = FALSE]
  
  # Create DGEList for library size calculation
  dge <- DGEList(counts = sample_data)
  
  # Calculate simple logCPM (NO TMM normalization, matching 05)
  log_cpm <- cpm(dge, log = TRUE, prior.count = 2)
  
  # Filter genes: keep genes with mean logCPM >= 1 (exactly as in 05)
  mean_expression <- rowMeans(log_cpm)
  keep_genes <- mean_expression >= 1
  
  log_cpm_filtered <- log_cpm[keep_genes, , drop = FALSE]
  
  cat(sprintf("  Genes after filtering (mean logCPM >= 1): %d (from %d)\n", 
              nrow(log_cpm_filtered), nrow(log_cpm)))
  
  return(list(
    log_cpm = log_cpm_filtered,
    n_genes_original = nrow(log_cpm),
    n_genes_final = nrow(log_cpm_filtered),
    library_sizes = dge$samples$lib.size
  ))
}

# CDM-based PCA (from v6)
perform_cdm_pca <- function(log_cpm_data, analysis_name, verbose = FALSE) {
  
  cat(sprintf("Performing CDM-based PCA for %s...\n", analysis_name))
  
  X <- log_cpm_data
  n_genes <- nrow(X)
  n_samples <- ncol(X)
  
  if (verbose) {
    cat(sprintf("  CDM input: %d genes x %d samples\n", n_genes, n_samples))
  }
  
  # Center data
  X_centered <- X - rowMeans(X)
  gene_sds <- apply(X_centered, 1, sd)
  
  # Remove zero-variance genes
  non_zero_var <- gene_sds > .Machine$double.eps
  if (sum(non_zero_var) < nrow(X)) {
    if (verbose) {
      cat(sprintf("  Removing %d zero-variance genes\n", sum(!non_zero_var)))
    }
    X_centered <- X_centered[non_zero_var, , drop = FALSE]
    gene_sds <- gene_sds[non_zero_var]
  }
  
  # Scale data
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
  
  # Limit components
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
    cat(sprintf("  PCA completed: %d components\n", ncol(scores)))
    cat(sprintf("  PC1-2 explain %.1f%% of variance\n", 
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
# 2. Load Data from v6 Results
# ==============================================================================

cat("\nLoading high-purity sample lists from v6 and count data...\n")

# Load v6 results
load("./data/processed/final_high_purity_sample_lists_v6.rda")
high_purity_sample_lists <- final_high_purity_sample_lists_v6

# Load count data
if (!exists("se_thyr")) {
  load("./data/raw/thyr_data.rda")
  se_thyr <- data
  rm(data)
}

# Get stranded_second count data
count_data_full <- assay(se_thyr, "stranded_second")
gene_info <- rowData(se_thyr)

# Filter to protein-coding genes
protein_coding_genes <- gene_info$gene_type == "protein_coding"
count_data_pc <- count_data_full[protein_coding_genes, ]

cat("Count data dimensions after protein-coding filter:", dim(count_data_pc), "\n")

# Display sample counts
cat("\nHigh-purity sample counts from v6:\n")
for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(high_purity_sample_lists)) {
    n_pairs <- length(high_purity_sample_lists[[group]]$tumor)
    cat(sprintf("  %s: %d pairs\n", group, n_pairs))
  }
}

# ==============================================================================
# 3. Group Comparison PCA Analysis
# ==============================================================================

cat("\n=== Group Comparison PCA Analysis ===\n")

comparison_groups <- list(
  "R0_R1_tumor" = list(groups = c("R0", "R1"), tissue = "tumor", 
                       title = "RET Tumor: R0 vs R1"),
  "R0_R1_normal" = list(groups = c("R0", "R1"), tissue = "normal",
                        title = "RET Normal: R0 vs R1"),
  "B0_B1_tumor" = list(groups = c("B0", "B1"), tissue = "tumor",
                       title = "BRAF Tumor: B0 vs B1"),
  "B0_B1_normal" = list(groups = c("B0", "B1"), tissue = "normal",
                        title = "BRAF Normal: B0 vs B1")
)

comparison_results <- list()

for (comp_name in names(comparison_groups)) {
  cat(sprintf("\n--- %s ---\n", comp_name))
  
  comp_info <- comparison_groups[[comp_name]]
  groups <- comp_info$groups
  tissue <- comp_info$tissue
  
  # Collect samples
  all_samples <- c()
  all_labels <- c()
  
  for (group in groups) {
    if (group %in% names(high_purity_sample_lists)) {
      group_samples <- high_purity_sample_lists[[group]][[tissue]]
      if (length(group_samples) > 0) {
        all_samples <- c(all_samples, group_samples)
        all_labels <- c(all_labels, rep(group, length(group_samples)))
      }
    }
  }
  
  if (length(all_samples) < 6) {
    cat(sprintf("Insufficient samples (%d), skipping...\n", length(all_samples)))
    comparison_results[[comp_name]] <- NULL
    next
  }
  
  cat(sprintf("Total samples: %d (%s: %d, %s: %d)\n", 
              length(all_samples),
              groups[1], sum(all_labels == groups[1]),
              groups[2], sum(all_labels == groups[2])))
  
  # Check sample availability
  available_samples <- all_samples[all_samples %in% colnames(count_data_pc)]
  available_labels <- all_labels[all_samples %in% colnames(count_data_pc)]
  
  if (length(available_samples) < length(all_samples)) {
    missing_count <- length(all_samples) - length(available_samples)
    cat(sprintf("Warning: %d samples missing from count data\n", missing_count))
  }
  
  if (length(available_samples) < 6) {
    cat("Insufficient available samples, skipping...\n")
    comparison_results[[comp_name]] <- NULL
    next
  }
  
  # Prepare logCPM data
  logcpm_data <- prepare_logcpm_data(count_data_pc, available_samples, comp_name)
  
  # Run PCA with CDM
  tryCatch({
    pca_result <- perform_cdm_pca(logcpm_data$log_cpm, comp_name, verbose = TRUE)
    
    # Calculate separation distance
    unique_groups <- unique(available_labels)
    if (length(unique_groups) == 2) {
      group1_scores <- pca_result$scores[available_labels == unique_groups[1], 1:2, drop = FALSE]
      group2_scores <- pca_result$scores[available_labels == unique_groups[2], 1:2, drop = FALSE]
      
      group1_centroid <- colMeans(group1_scores)
      group2_centroid <- colMeans(group2_scores)
      separation_distance <- sqrt(sum((group1_centroid - group2_centroid)^2))
      
      # Calculate separation quality
      if (separation_distance > 2) {
        separation_quality <- "Excellent"
      } else if (separation_distance > 1) {
        separation_quality <- "Good"
      } else if (separation_distance > 0.5) {
        separation_quality <- "Moderate"
      } else {
        separation_quality <- "Poor"
      }
    } else {
      separation_distance <- NA
      separation_quality <- NA
    }
    
    comparison_results[[comp_name]] <- list(
      pca = pca_result,
      samples = available_samples,
      labels = available_labels,
      title = comp_info$title,
      separation_distance = separation_distance,
      separation_quality = separation_quality,
      sample_counts = table(available_labels),
      n_genes = logcpm_data$n_genes_final
    )
    
    cat(sprintf("Results: PC1-2 variance=%.1f%%, Separation=%.3f (%s)\n",
                sum(pca_result$variance_explained[1:2]) * 100, 
                separation_distance,
                separation_quality))
    
  }, error = function(e) {
    cat(sprintf("Error in PCA: %s\n", e$message))
    comparison_results[[comp_name]] <- NULL
  })
}

# ==============================================================================
# 4. Clustering Analysis (Simplified for v6)
# ==============================================================================

cat("\n=== Clustering Analysis ===\n")

clustering_results <- list()

for (comp_name in names(comparison_results)) {
  if (is.null(comparison_results[[comp_name]])) {
    next
  }
  
  cat(sprintf("\n--- Clustering for %s ---\n", comp_name))
  
  result <- comparison_results[[comp_name]]
  pca_scores <- result$pca$scores
  
  # Use first 3 PCs for clustering
  n_pcs <- min(3, ncol(pca_scores))
  pca_scores_for_clustering <- pca_scores[, 1:n_pcs, drop = FALSE]
  
  # K-means clustering
  set.seed(123)
  kmeans_result <- kmeans(pca_scores_for_clustering, centers = 2, nstart = 25)
  
  # Hierarchical clustering
  dist_matrix <- dist(pca_scores_for_clustering)
  hclust_result <- hclust(dist_matrix, method = "ward.D2")
  hclust_clusters <- cutree(hclust_result, k = 2)
  
  # Cross-tabulation
  kmeans_table <- table(Original = result$labels, KMeans = kmeans_result$cluster)
  hclust_table <- table(Original = result$labels, HClust = hclust_clusters)
  
  cat("K-means clustering:\n")
  print(kmeans_table)
  cat("\nHierarchical clustering:\n")
  print(hclust_table)
  
  # Calculate purity
  calc_purity <- function(cluster_table) {
    max_per_cluster <- apply(cluster_table, 2, max)
    sum(max_per_cluster) / sum(cluster_table)
  }
  
  kmeans_purity <- calc_purity(kmeans_table)
  hclust_purity <- calc_purity(hclust_table)
  
  cat(sprintf("\nClustering purity: K-means=%.3f, Hierarchical=%.3f\n", 
              kmeans_purity, hclust_purity))
  
  clustering_results[[comp_name]] <- list(
    kmeans_clusters = kmeans_result$cluster,
    hclust_clusters = hclust_clusters,
    kmeans_table = kmeans_table,
    hclust_table = hclust_table,
    kmeans_purity = kmeans_purity,
    hclust_purity = hclust_purity,
    n_pcs_used = n_pcs
  )
}

# ==============================================================================
# 5. Summary Statistics
# ==============================================================================

cat("\n=== Phase 2 Summary Statistics ===\n")

# Create summary table
phase2_summary <- data.frame(
  Comparison = character(0),
  N_Samples = integer(0),
  N_Genes = integer(0),
  PC1_Var = numeric(0),
  PC2_Var = numeric(0),
  Separation = numeric(0),
  Quality = character(0),
  KMeans_Purity = numeric(0),
  HClust_Purity = numeric(0),
  stringsAsFactors = FALSE
)

for (comp_name in names(comparison_results)) {
  if (!is.null(comparison_results[[comp_name]])) {
    result <- comparison_results[[comp_name]]
    clust_result <- clustering_results[[comp_name]]
    
    phase2_summary <- rbind(phase2_summary, data.frame(
      Comparison = comp_name,
      N_Samples = length(result$samples),
      N_Genes = result$n_genes,
      PC1_Var = round(result$pca$variance_explained[1] * 100, 1),
      PC2_Var = round(result$pca$variance_explained[2] * 100, 1),
      Separation = round(result$separation_distance, 3),
      Quality = result$separation_quality,
      KMeans_Purity = round(clust_result$kmeans_purity, 3),
      HClust_Purity = round(clust_result$hclust_purity, 3),
      stringsAsFactors = FALSE
    ))
  }
}

print(phase2_summary)

# ==============================================================================
# 6. Save Results
# ==============================================================================

cat("\nSaving Phase 2 PCA v6 results...\n")

if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

phase2_results_v6 <- list(
  comparison_results = comparison_results,
  clustering_results = clustering_results,
  phase2_summary = phase2_summary,
  high_purity_sample_lists = high_purity_sample_lists,
  analysis_date = Sys.time(),
  analysis_version = "v6_logcpm_cdm"
)

save(phase2_results_v6, file = "./data/processed/phase2_results_v6.rda")

if (!dir.exists("./output/plots")) {
  dir.create("./output/plots", recursive = TRUE)
}

# ==============================================================================
# 7. Visualization Preparation
# ==============================================================================

cat("\n=== Preparing visualization data ===\n")

# Export PCA scores for plotting
for (comp_name in names(comparison_results)) {
  if (!is.null(comparison_results[[comp_name]])) {
    result <- comparison_results[[comp_name]]
    
    # Create data frame for plotting
    plot_data <- data.frame(
      Sample = result$samples,
      Group = result$labels,
      PC1 = result$pca$scores[, 1],
      PC2 = result$pca$scores[, 2],
      stringsAsFactors = FALSE
    )
    
    # Add clustering results if available
    if (comp_name %in% names(clustering_results)) {
      plot_data$KMeans_Cluster <- clustering_results[[comp_name]]$kmeans_clusters
      plot_data$HClust_Cluster <- clustering_results[[comp_name]]$hclust_clusters
    }
    
    # Save as CSV for external plotting
    write.csv(plot_data, 
              file = paste0("./output/plots/pca_scores_", comp_name, "_v6.csv"),
              row.names = FALSE)
    
    cat(sprintf("Saved PCA scores for %s\n", comp_name))
  }
}

# ==============================================================================
# 8. Final Assessment
# ==============================================================================

cat("\n=== Phase 2 Analysis Assessment ===\n")

# Count successful analyses
n_successful <- sum(!sapply(comparison_results, is.null))
n_total <- length(comparison_groups)

cat(sprintf("Successful comparisons: %d/%d\n", n_successful, n_total))

# Assess separation quality
quality_counts <- table(phase2_summary$Quality)
if (length(quality_counts) > 0) {
  cat("\nSeparation quality distribution:\n")
  for (q in names(quality_counts)) {
    cat(sprintf("  %s: %d\n", q, quality_counts[q]))
  }
}

# Assess clustering quality
high_purity_clustering <- sum(phase2_summary$KMeans_Purity > 0.7 | 
                                phase2_summary$HClust_Purity > 0.7)
if (high_purity_clustering > 0) {
  cat(sprintf("\nHigh-purity clustering (>0.7): %d comparisons\n", high_purity_clustering))
}

# Recommendations
cat("\n=== Recommendations ===\n")

if (n_successful == n_total) {
  cat("✓ All group comparisons successful\n")
  cat("✓ Proceed with DEG analysis using these groupings\n")
} else {
  cat("⚠ Some comparisons failed - review sample availability\n")
}

best_separation <- phase2_summary[which.max(phase2_summary$Separation), ]
if (nrow(best_separation) > 0) {
  cat(sprintf("\nBest separation: %s (distance=%.3f)\n", 
              best_separation$Comparison, best_separation$Separation))
  cat("Consider prioritizing this comparison for DEG analysis\n")
}

cat("\nPhase 2 PCA v6 analysis completed successfully!\n")

