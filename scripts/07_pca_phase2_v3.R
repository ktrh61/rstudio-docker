# ==============================================================================
# REBC-THYR PCA Phase 2 Analysis v3 - Group Comparison + R0+R1 Clustering
# 07_pca_phase2_v3.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(dplyr)
library(cluster)
library(Rcpp)
library(RhpcBLASctl)

# Source functions
source("./utils/with_openblas_threads.R")
sourceCpp("./utils/CDM_fast3_arma_enhanced.cpp")

cat("Starting PCA Phase 2 v3: Group comparison + R0+R1 clustering...\n")

# ==============================================================================
# 1. CDM_fast Compatible Wrapper
# ==============================================================================

CDM_fast_compatible <- function(X, verbose = FALSE) {
  sample_names <- colnames(X)
  
  result <- with_openblas_threads("auto-2", {
    CDM_fast3_arma(X, verbose = verbose)
  })
  
  if (!is.null(result$scores) && !is.null(sample_names)) {
    rownames(result$scores) <- sample_names
    colnames(result$scores) <- paste0("PC", 1:ncol(result$scores))
  }
  
  variance_explained <- result$values^2 / sum(result$values^2)
  result$variance_explained <- variance_explained
  result$cumulative_variance <- cumsum(variance_explained)
  
  return(result)
}

# ==============================================================================
# 2. Load Data
# ==============================================================================

cat("Loading purity-filtered sample lists and TPM data...\n")

load("./data/processed/purity_filtered_sample_lists.rda")
load("./data/processed/thyr_tpm.rda")

# Display final sample counts
cat("High-purity sample counts:\n")
for (group in c("R0", "R1", "B0", "B1")) {
  if (group %in% names(purity_filtered_sample_lists)) {
    n_pairs <- length(purity_filtered_sample_lists[[group]]$tumor)
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
    if (group %in% names(purity_filtered_sample_lists)) {
      group_samples <- purity_filtered_sample_lists[[group]][[tissue]]
      if (length(group_samples) > 0) {
        all_samples <- c(all_samples, group_samples)
        all_labels <- c(all_labels, rep(group, length(group_samples)))
      }
    }
  }
  
  if (length(all_samples) < 6) {
    cat(sprintf("Insufficient samples (%d), skipping...\n", length(all_samples)))
    next
  }
  
  cat(sprintf("Total samples: %d\n", length(all_samples)))
  
  # Extract TPM and filter
  available_samples <- all_samples[all_samples %in% colnames(tpm)]
  available_labels <- all_labels[all_samples %in% colnames(tpm)]
  
  combined_tpm <- tpm[, available_samples, drop = FALSE]
  non_zero_genes <- rowSums(combined_tpm > 0) > 0
  combined_tpm_filtered <- combined_tpm[non_zero_genes, ]
  
  cat(sprintf("Using %d genes for PCA\n", nrow(combined_tpm_filtered)))
  
  # Run PCA
  pca_result <- CDM_fast_compatible(combined_tpm_filtered, verbose = FALSE)
  
  # Calculate separation distance
  unique_groups <- unique(available_labels)
  if (length(unique_groups) == 2) {
    group1_scores <- pca_result$scores[available_labels == unique_groups[1], 1:2]
    group2_scores <- pca_result$scores[available_labels == unique_groups[2], 1:2]
    
    group1_centroid <- colMeans(group1_scores)
    group2_centroid <- colMeans(group2_scores)
    separation_distance <- sqrt(sum((group1_centroid - group2_centroid)^2))
  } else {
    separation_distance <- NA
  }
  
  comparison_results[[comp_name]] <- list(
    pca = pca_result,
    samples = available_samples,
    labels = available_labels,
    title = comp_info$title,
    separation_distance = separation_distance,
    sample_counts = table(available_labels)
  )
  
  cat(sprintf("PC1-2 variance: %.1f%%, Separation: %.3f\n",
              sum(pca_result$variance_explained[1:2]) * 100, separation_distance))
}

# ==============================================================================
# 4. R0+R1 Clustering Analysis (Main Focus)
# ==============================================================================

cat("\n=== R0+R1 Clustering Analysis ===\n")

# Prepare R0+R1 combined data
r0_tumor <- purity_filtered_sample_lists$R0$tumor
r0_normal <- purity_filtered_sample_lists$R0$normal
r1_tumor <- purity_filtered_sample_lists$R1$tumor  
r1_normal <- purity_filtered_sample_lists$R1$normal

# Combine tumor samples
tumor_samples <- c(r0_tumor, r1_tumor)
tumor_labels <- c(rep("R0", length(r0_tumor)), rep("R1", length(r1_tumor)))

cat(sprintf("R0+R1 tumor clustering: %d samples (R0:%d, R1:%d)\n",
            length(tumor_samples), length(r0_tumor), length(r1_tumor)))

# Extract and filter TPM
available_tumor <- tumor_samples[tumor_samples %in% colnames(tpm)]
available_tumor_labels <- tumor_labels[tumor_samples %in% colnames(tpm)]

r01_tpm <- tpm[, available_tumor, drop = FALSE]
non_zero_genes <- rowSums(r01_tpm > 0) > 0
r01_tpm_filtered <- r01_tpm[non_zero_genes, ]

cat(sprintf("Using %d genes for R0+R1 analysis\n", nrow(r01_tpm_filtered)))

# PCA
r01_pca <- CDM_fast_compatible(r01_tpm_filtered, verbose = FALSE)

cat(sprintf("R0+R1 PCA: PC1-2 explains %.1f%% variance\n",
            sum(r01_pca$variance_explained[1:2]) * 100))

# K-means clustering (k=2)
set.seed(123)
kmeans_result <- kmeans(r01_pca$scores[, 1:3], centers = 2, nstart = 25)

# Hierarchical clustering
dist_matrix <- dist(r01_pca$scores[, 1:3])
hclust_result <- hclust(dist_matrix, method = "ward.D2")
hclust_clusters <- cutree(hclust_result, k = 2)

# Silhouette analysis
sil_kmeans <- silhouette(kmeans_result$cluster, dist_matrix)
sil_hclust <- silhouette(hclust_clusters, dist_matrix)

# Analyze clustering vs original labels
r01_clustering_results <- list(
  original_labels = available_tumor_labels,
  kmeans_clusters = kmeans_result$cluster,
  hclust_clusters = hclust_clusters,
  scores = r01_pca$scores,
  silhouette_kmeans = mean(sil_kmeans[, 3]),
  silhouette_hclust = mean(sil_hclust[, 3])
)

# Cross-tabulation analysis
cat("\nClustering vs Original Labels:\n")
cat("K-means clustering:\n")
kmeans_table <- table(Original = available_tumor_labels, Cluster = kmeans_result$cluster)
print(kmeans_table)

cat("Hierarchical clustering:\n")
hclust_table <- table(Original = available_tumor_labels, Cluster = hclust_clusters)
print(hclust_table)

# Calculate clustering purity
calc_purity <- function(cluster_table) {
  max_per_cluster <- apply(cluster_table, 2, max)
  sum(max_per_cluster) / sum(cluster_table)
}

kmeans_purity <- calc_purity(kmeans_table)
hclust_purity <- calc_purity(hclust_table)

cat(sprintf("\nClustering purity: K-means=%.3f, Hierarchical=%.3f\n", 
            kmeans_purity, hclust_purity))
cat(sprintf("Silhouette scores: K-means=%.3f, Hierarchical=%.3f\n",
            r01_clustering_results$silhouette_kmeans, r01_clustering_results$silhouette_hclust))

# ==============================================================================
# 5. Clustering Interpretation
# ==============================================================================

cat("\n=== Clustering Interpretation ===\n")

# Determine best clustering method
best_method <- if (r01_clustering_results$silhouette_kmeans > r01_clustering_results$silhouette_hclust) {
  "kmeans"
} else {
  "hclust"
}

best_clusters <- if (best_method == "kmeans") {
  kmeans_result$cluster
} else {
  hclust_clusters
}

best_table <- if (best_method == "kmeans") kmeans_table else hclust_table
best_purity <- if (best_method == "kmeans") kmeans_purity else hclust_purity

cat(sprintf("Best method: %s (purity=%.3f)\n", best_method, best_purity))

# Analyze radiation exposure hypothesis
if (best_purity > 0.7) {
  cat("✅ High clustering purity suggests potential radiation exposure separation\n")
  
  # Identify potential radiation-exposed cluster
  r0_in_cluster1 <- sum(best_table["R0", 1])
  r0_in_cluster2 <- sum(best_table["R0", 2])
  
  if (r0_in_cluster1 > r0_in_cluster2) {
    # Cluster 1 is likely non-exposed (more R0)
    cat("Cluster 1: Likely non-exposed (radiation negative)\n")
    cat("Cluster 2: Likely exposed (radiation positive)\n")
    radiation_cluster <- 2
  } else {
    # Cluster 2 is likely non-exposed (more R0)
    cat("Cluster 1: Likely exposed (radiation positive)\n") 
    cat("Cluster 2: Likely non-exposed (radiation negative)\n")
    radiation_cluster <- 1
  }
  
  # Count mixed cases in R1
  r1_in_radiation_cluster <- best_table["R1", radiation_cluster]
  r1_in_non_radiation_cluster <- best_table["R1", 3 - radiation_cluster]
  
  cat(sprintf("R1 cases - Radiation cluster: %d, Non-radiation cluster: %d\n",
              r1_in_radiation_cluster, r1_in_non_radiation_cluster))
  
  estimated_radiation_rate <- r1_in_radiation_cluster / sum(best_table["R1", ])
  cat(sprintf("Estimated R1 radiation rate: %.1f%%\n", estimated_radiation_rate * 100))
  
} else {
  cat("⚠️ Low clustering purity suggests limited separation\n")
}

# ==============================================================================
# 4. Comprehensive Clustering Analysis (All Groups)
# ==============================================================================

cat("\n=== Comprehensive Clustering Analysis ===\n")

clustering_configs <- list(
  "R0_R1_tumor" = list(
    samples = list(r0 = purity_filtered_sample_lists$R0$tumor,
                   r1 = purity_filtered_sample_lists$R1$tumor),
    title = "R0+R1 Tumor Clustering"
  ),
  "R0_R1_normal" = list(
    samples = list(r0 = purity_filtered_sample_lists$R0$normal,
                   r1 = purity_filtered_sample_lists$R1$normal),
    title = "R0+R1 Normal Clustering"
  ),
  "B0_B1_tumor" = list(
    samples = list(b0 = purity_filtered_sample_lists$B0$tumor,
                   b1 = purity_filtered_sample_lists$B1$tumor),
    title = "B0+B1 Tumor Clustering"
  ),
  "B0_B1_normal" = list(
    samples = list(b0 = purity_filtered_sample_lists$B0$normal,
                   b1 = purity_filtered_sample_lists$B1$normal),
    title = "B0+B1 Normal Clustering"
  )
)

# Initialize clustering results
clustering_results <- list()

for (config_name in names(clustering_configs)) {
  cat(sprintf("\n--- %s ---\n", config_name))
  
  config <- clustering_configs[[config_name]]
  
  # Combine samples
  group1_samples <- config$samples[[1]]
  group2_samples <- config$samples[[2]]
  combined_samples <- c(group1_samples, group2_samples)
  combined_labels <- c(rep(names(config$samples)[1], length(group1_samples)),
                       rep(names(config$samples)[2], length(group2_samples)))
  
  cat(sprintf("%s: %d samples (%s:%d, %s:%d)\n",
              config$title, length(combined_samples),
              names(config$samples)[1], length(group1_samples),
              names(config$samples)[2], length(group2_samples)))
  
  # Check minimum requirements
  if (length(combined_samples) < 10 || length(group1_samples) < 3 || length(group2_samples) < 3) {
    cat("Insufficient samples for clustering, skipping...\n")
    clustering_results[[config_name]] <- NULL
    next
  }
  
  # Extract and filter TPM
  available_samples <- combined_samples[combined_samples %in% colnames(tpm)]
  available_labels <- combined_labels[combined_samples %in% colnames(tpm)]
  
  if (length(available_samples) < length(combined_samples)) {
    cat(sprintf("Warning: %d samples missing from TPM data\n", 
                length(combined_samples) - length(available_samples)))
  }
  
  combined_tpm <- tpm[, available_samples, drop = FALSE]
  non_zero_genes <- rowSums(combined_tpm > 0) > 0
  combined_tpm_filtered <- combined_tpm[non_zero_genes, ]
  
  cat(sprintf("Using %d genes for clustering analysis\n", nrow(combined_tpm_filtered)))
  
  # PCA
  pca_result <- CDM_fast_compatible(combined_tpm_filtered, verbose = FALSE)
  
  cat(sprintf("PCA: PC1-2 explains %.1f%% variance\n",
              sum(pca_result$variance_explained[1:2]) * 100))
  
  # K-means clustering (k=2)
  set.seed(123)
  tryCatch({
    kmeans_result <- kmeans(pca_result$scores[, 1:3], centers = 2, nstart = 25)
    
    # Hierarchical clustering
    dist_matrix <- dist(pca_result$scores[, 1:3])
    hclust_result <- hclust(dist_matrix, method = "ward.D2")
    hclust_clusters <- cutree(hclust_result, k = 2)
    
    # Silhouette analysis
    sil_kmeans <- silhouette(kmeans_result$cluster, dist_matrix)
    sil_hclust <- silhouette(hclust_clusters, dist_matrix)
    
    # Cross-tabulation analysis
    kmeans_table <- table(Original = available_labels, Cluster = kmeans_result$cluster)
    hclust_table <- table(Original = available_labels, Cluster = hclust_clusters)
    
    cat("K-means clustering:\n")
    print(kmeans_table)
    cat("Hierarchical clustering:\n")
    print(hclust_table)
    
    # Calculate clustering purity
    calc_purity <- function(cluster_table) {
      max_per_cluster <- apply(cluster_table, 2, max)
      sum(max_per_cluster) / sum(cluster_table)
    }
    
    kmeans_purity <- calc_purity(kmeans_table)
    hclust_purity <- calc_purity(hclust_table)
    
    cat(sprintf("Clustering purity: K-means=%.3f, Hierarchical=%.3f\n", 
                kmeans_purity, hclust_purity))
    cat(sprintf("Silhouette scores: K-means=%.3f, Hierarchical=%.3f\n",
                mean(sil_kmeans[, 3]), mean(sil_hclust[, 3])))
    
    # Store results
    clustering_results[[config_name]] <- list(
      pca = pca_result,
      samples = available_samples,
      labels = available_labels,
      kmeans_result = kmeans_result,
      hclust_clusters = hclust_clusters,
      kmeans_table = kmeans_table,
      hclust_table = hclust_table,
      kmeans_purity = kmeans_purity,
      hclust_purity = hclust_purity,
      silhouette_kmeans = mean(sil_kmeans[, 3]),
      silhouette_hclust = mean(sil_hclust[, 3]),
      title = config$title
    )
    
    # Interpretation
    best_purity <- max(kmeans_purity, hclust_purity)
    if (best_purity > 0.7) {
      cat("✅ High clustering purity - potential distinct subgroups\n")
    } else if (best_purity > 0.6) {
      cat("⚠️ Moderate clustering purity - some separation\n")
    } else {
      cat("❌ Low clustering purity - limited separation\n")
    }
    
  }, error = function(e) {
    cat(sprintf("Error in clustering analysis: %s\n", e$message))
    clustering_results[[config_name]] <- NULL
  })
}

# ==============================================================================
# 5. Save Results
# ==============================================================================

cat("Saving Phase 2 PCA and clustering results...\n")

if (!dir.exists("./data/processed")) {
  dir.create("./data/processed", recursive = TRUE)
}

# Ensure clustering_results exists (fallback for errors)
if (!exists("clustering_results")) {
  cat("Warning: clustering_results not found, creating empty list\n")
  clustering_results <- list()
}

phase2_clustering_results <- list(
  comparison_results = comparison_results,
  clustering_results = clustering_results,
  analysis_date = Sys.time(),
  analysis_version = "v3_comprehensive_clustering"
)

save(phase2_clustering_results, file = "./data/processed/phase2_clustering_results.rda")

# ==============================================================================
# 6. Summary
# ==============================================================================

cat("\n==============================================\n")
cat("Phase 2 PCA + Comprehensive Clustering Summary\n")
cat("==============================================\n")

cat("Group Comparisons:\n")
for (comp_name in names(comparison_results)) {
  if (!is.null(comparison_results[[comp_name]])) {
    result <- comparison_results[[comp_name]]
    cat(sprintf("  %s: %.1f%% variance, separation=%.3f\n",
                comp_name, 
                sum(result$pca$variance_explained[1:2]) * 100,
                result$separation_distance))
  }
}

cat("\nClustering Analysis Summary:\n")
# Check if clustering_results exists and has content
if (exists("clustering_results") && length(clustering_results) > 0) {
  for (config_name in names(clustering_results)) {
    if (!is.null(clustering_results[[config_name]])) {
      result <- clustering_results[[config_name]]
      best_purity <- max(result$kmeans_purity, result$hclust_purity)
      cat(sprintf("  %s: Best purity=%.3f\n", config_name, best_purity))
    }
  }
} else {
  cat("  No clustering results available\n")
}

# Overall assessment
high_purity_count <- 0
moderate_purity_count <- 0
low_purity_count <- 0

if (exists("clustering_results") && length(clustering_results) > 0) {
  for (config_name in names(clustering_results)) {
    if (!is.null(clustering_results[[config_name]])) {
      result <- clustering_results[[config_name]]
      best_purity <- max(result$kmeans_purity, result$hclust_purity)
      
      if (best_purity > 0.7) {
        high_purity_count <- high_purity_count + 1
      } else if (best_purity > 0.6) {
        moderate_purity_count <- moderate_purity_count + 1
      } else {
        low_purity_count <- low_purity_count + 1
      }
    }
  }
}

cat(sprintf("\nOverall clustering assessment:\n"))
cat(sprintf("  High purity (>70%%): %d comparisons\n", high_purity_count))
cat(sprintf("  Moderate purity (60-70%%): %d comparisons\n", moderate_purity_count))
cat(sprintf("  Low purity (<60%%): %d comparisons\n", low_purity_count))

if (high_purity_count > 0) {
  cat("\n✅ Some groups show distinct subpopulations - proceed to DEG analysis\n")
} else if (moderate_purity_count > 0) {
  cat("\n⚠️ Moderate separation observed - careful DEG analysis recommended\n")
} else {
  cat("\n❌ Limited separation - consider alternative approaches or skip DEG analysis\n")
}

cat("\nNext steps:\n")
if (high_purity_count > 0 || moderate_purity_count > 0) {
  cat("1. Proceed to DEGES normalization (08)\n")
  cat("2. DEG analysis focusing on high-separation groups\n")
  cat("3. Consider clustering-informed sample selection\n")
} else {
  cat("1. Review data quality and filtering thresholds\n")
  cat("2. Consider alternative stratification approaches\n")
  cat("3. May skip extensive DEG analysis due to limited separation\n")
}

cat("Phase 2 comprehensive clustering analysis completed!\n")
cat("==============================================\n")

