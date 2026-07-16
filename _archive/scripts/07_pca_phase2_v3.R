# ==============================================================================
# REBC-THYR PCA Phase 2 Analysis v3 - Group Comparison + Clustering (Corrected)
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

cat("Starting PCA Phase 2 v3: Group comparison + clustering (corrected)...\n")

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
# 2. Load Data (Corrected Input)
# ==============================================================================

cat("Loading high-purity sample lists and TPM data...\n")

# Load the correct output from 06_purity_analysis_v3.R
load("./data/processed/final_high_purity_sample_lists.rda")
load("./data/processed/thyr_tpm.rda")

# Use the correct variable name from 06
high_purity_sample_lists <- final_high_purity_sample_lists

# Display final sample counts
cat("High-purity sample counts (input from 06):\n")
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
  
  # Collect samples from high-purity lists
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
  
  cat(sprintf("Total samples: %d\n", length(all_samples)))
  
  # Extract TPM and filter
  available_samples <- all_samples[all_samples %in% colnames(tpm)]
  available_labels <- all_labels[all_samples %in% colnames(tpm)]
  
  if (length(available_samples) < length(all_samples)) {
    missing_count <- length(all_samples) - length(available_samples)
    cat(sprintf("Warning: %d samples missing from TPM data\n", missing_count))
  }
  
  combined_tpm <- tpm[, available_samples, drop = FALSE]
  non_zero_genes <- rowSums(combined_tpm > 0) > 0
  combined_tpm_filtered <- combined_tpm[non_zero_genes, ]
  
  cat(sprintf("Using %d genes for PCA\n", nrow(combined_tpm_filtered)))
  
  # Run PCA
  tryCatch({
    pca_result <- CDM_fast_compatible(combined_tpm_filtered, verbose = FALSE)
    
    # Calculate separation distance
    unique_groups <- unique(available_labels)
    if (length(unique_groups) == 2) {
      group1_scores <- pca_result$scores[available_labels == unique_groups[1], 1:2, drop = FALSE]
      group2_scores <- pca_result$scores[available_labels == unique_groups[2], 1:2, drop = FALSE]
      
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
      sample_counts = table(available_labels),
      n_genes = nrow(combined_tpm_filtered)
    )
    
    cat(sprintf("PC1-2 variance: %.1f%%, Separation: %.3f\n",
                sum(pca_result$variance_explained[1:2]) * 100, separation_distance))
    
  }, error = function(e) {
    cat(sprintf("Error in PCA for %s: %s\n", comp_name, e$message))
    comparison_results[[comp_name]] <- NULL
  })
}

# ==============================================================================
# 4. Comprehensive Clustering Analysis
# ==============================================================================

cat("\n=== Comprehensive Clustering Analysis ===\n")

# Define clustering configurations using high-purity sample lists
clustering_configs <- list(
  "R0_R1_tumor" = list(
    samples = list(R0 = high_purity_sample_lists$R0$tumor,
                   R1 = high_purity_sample_lists$R1$tumor),
    title = "R0+R1 Tumor Clustering"
  ),
  "R0_R1_normal" = list(
    samples = list(R0 = high_purity_sample_lists$R0$normal,
                   R1 = high_purity_sample_lists$R1$normal),
    title = "R0+R1 Normal Clustering"
  ),
  "B0_B1_tumor" = list(
    samples = list(B0 = high_purity_sample_lists$B0$tumor,
                   B1 = high_purity_sample_lists$B1$tumor),
    title = "B0+B1 Tumor Clustering"
  ),
  "B0_B1_normal" = list(
    samples = list(B0 = high_purity_sample_lists$B0$normal,
                   B1 = high_purity_sample_lists$B1$normal),
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
  
  # Check if samples exist
  if (length(group1_samples) == 0 || length(group2_samples) == 0) {
    cat("One or both groups have no samples, skipping...\n")
    clustering_results[[config_name]] <- NULL
    next
  }
  
  combined_samples <- c(group1_samples, group2_samples)
  combined_labels <- c(rep(names(config$samples)[1], length(group1_samples)),
                       rep(names(config$samples)[2], length(group2_samples)))
  
  cat(sprintf("%s: %d samples (%s:%d, %s:%d)\n",
              config$title, length(combined_samples),
              names(config$samples)[1], length(group1_samples),
              names(config$samples)[2], length(group2_samples)))
  
  # Check minimum requirements
  if (length(combined_samples) < 8 || length(group1_samples) < 3 || length(group2_samples) < 3) {
    cat("Insufficient samples for clustering, skipping...\n")
    clustering_results[[config_name]] <- NULL
    next
  }
  
  # Extract and filter TPM
  available_samples <- combined_samples[combined_samples %in% colnames(tpm)]
  available_labels <- combined_labels[combined_samples %in% colnames(tpm)]
  
  if (length(available_samples) < length(combined_samples)) {
    missing_count <- length(combined_samples) - length(available_samples)
    cat(sprintf("Warning: %d samples missing from TPM data\n", missing_count))
  }
  
  if (length(available_samples) < 8) {
    cat("Insufficient available samples for clustering, skipping...\n")
    clustering_results[[config_name]] <- NULL
    next
  }
  
  combined_tpm <- tpm[, available_samples, drop = FALSE]
  non_zero_genes <- rowSums(combined_tpm > 0) > 0
  combined_tpm_filtered <- combined_tpm[non_zero_genes, ]
  
  cat(sprintf("Using %d genes for clustering analysis\n", nrow(combined_tpm_filtered)))
  
  # PCA
  tryCatch({
    pca_result <- CDM_fast_compatible(combined_tpm_filtered, verbose = FALSE)
    
    cat(sprintf("PCA: PC1-2 explains %.1f%% variance\n",
                sum(pca_result$variance_explained[1:2]) * 100))
    
    # Use first 3 PCs for clustering (or fewer if not available)
    n_pcs <- min(3, ncol(pca_result$scores))
    pca_scores_for_clustering <- pca_result$scores[, 1:n_pcs, drop = FALSE]
    
    # K-means clustering (k=2)
    set.seed(123)
    kmeans_result <- kmeans(pca_scores_for_clustering, centers = 2, nstart = 25)
    
    # Hierarchical clustering
    dist_matrix <- dist(pca_scores_for_clustering)
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
      title = config$title,
      n_genes = nrow(combined_tpm_filtered),
      n_pcs_used = n_pcs
    )
    
    # Interpretation
    best_purity <- max(kmeans_purity, hclust_purity)
    if (best_purity > 0.7) {
      cat("✅ High clustering purity - potential distinct subgroups\n")
      
      # Special analysis for R0+R1 combinations
      if (grepl("R0.*R1", config_name)) {
        # Determine which method is better
        best_method <- ifelse(kmeans_purity > hclust_purity, "kmeans", "hclust")
        best_table <- ifelse(kmeans_purity > hclust_purity, 
                             list(kmeans_table), list(hclust_table))[[1]]
        
        cat(sprintf("Best method for %s: %s (purity=%.3f)\n", 
                    config_name, best_method, best_purity))
        
        # Analyze radiation exposure hypothesis
        r0_in_cluster1 <- best_table["R0", 1]
        r0_in_cluster2 <- best_table["R0", 2]
        r1_in_cluster1 <- best_table["R1", 1]
        r1_in_cluster2 <- best_table["R1", 2]
        
        if (r0_in_cluster1 > r0_in_cluster2) {
          # Cluster 1 is likely non-exposed (more R0)
          radiation_cluster <- 2
          cat("Cluster 1: Likely non-exposed (more R0)\n")
          cat("Cluster 2: Likely exposed (more R1)\n")
        } else {
          # Cluster 2 is likely non-exposed (more R0)
          radiation_cluster <- 1
          cat("Cluster 1: Likely exposed (more R1)\n")
          cat("Cluster 2: Likely non-exposed (more R0)\n")
        }
        
        # Calculate estimated radiation exposure rate in R1
        r1_total <- r1_in_cluster1 + r1_in_cluster2
        r1_exposed <- ifelse(radiation_cluster == 1, r1_in_cluster1, r1_in_cluster2)
        
        if (r1_total > 0) {
          estimated_radiation_rate <- r1_exposed / r1_total
          cat(sprintf("Estimated R1 radiation exposure rate: %.1f%% (%d/%d)\n",
                      estimated_radiation_rate * 100, r1_exposed, r1_total))
        }
      }
      
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

# Comprehensive Phase 2 results
phase2_clustering_results_v3 <- list(
  comparison_results = comparison_results,
  clustering_results = clustering_results,
  high_purity_sample_lists = high_purity_sample_lists,
  analysis_date = Sys.time(),
  analysis_version = "v3_corrected_input_comprehensive"
)

save(phase2_clustering_results_v3, file = "./data/processed/phase2_clustering_results_v3.rda")

cat("Results saved to ./data/processed/phase2_clustering_results_v3.rda\n")

# ==============================================================================
# 6. Summary and Assessment
# ==============================================================================

cat("\n==============================================\n")
cat("Phase 2 PCA + Clustering Summary (Corrected)\n")
cat("==============================================\n")

cat("Group Comparisons:\n")
valid_comparisons <- 0
for (comp_name in names(comparison_results)) {
  if (!is.null(comparison_results[[comp_name]])) {
    result <- comparison_results[[comp_name]]
    cat(sprintf("  %s: %.1f%% variance, separation=%.3f\n",
                comp_name, 
                sum(result$pca$variance_explained[1:2]) * 100,
                result$separation_distance))
    valid_comparisons <- valid_comparisons + 1
  } else {
    cat(sprintf("  %s: FAILED (insufficient samples)\n", comp_name))
  }
}

cat(sprintf("\nValid group comparisons: %d/4\n", valid_comparisons))

cat("\nClustering Analysis Summary:\n")
high_purity_count <- 0
moderate_purity_count <- 0
low_purity_count <- 0
valid_clustering <- 0

for (config_name in names(clustering_results)) {
  if (!is.null(clustering_results[[config_name]])) {
    result <- clustering_results[[config_name]]
    best_purity <- max(result$kmeans_purity, result$hclust_purity)
    cat(sprintf("  %s: Best purity=%.3f\n", config_name, best_purity))
    
    if (best_purity > 0.7) {
      high_purity_count <- high_purity_count + 1
    } else if (best_purity > 0.6) {
      moderate_purity_count <- moderate_purity_count + 1
    } else {
      low_purity_count <- low_purity_count + 1
    }
    valid_clustering <- valid_clustering + 1
  } else {
    cat(sprintf("  %s: FAILED (insufficient samples)\n", config_name))
  }
}

cat(sprintf("\nValid clustering analyses: %d/4\n", valid_clustering))
cat(sprintf("Clustering quality assessment:\n"))
cat(sprintf("  High purity (>70%%): %d\n", high_purity_count))
cat(sprintf("  Moderate purity (60-70%%): %d\n", moderate_purity_count))
cat(sprintf("  Low purity (<60%%): %d\n", low_purity_count))

# Overall recommendation
if (high_purity_count > 0) {
  cat("\n✅ EXCELLENT: Strong subgroup separation detected\n")
  cat("Recommendation: Proceed with DEG analysis focusing on high-separation groups\n")
} else if (moderate_purity_count > 0) {
  cat("\n⚠️ MODERATE: Some separation observed\n")
  cat("Recommendation: Careful DEG analysis with interpretation caveats\n")
} else if (valid_clustering > 0) {
  cat("\n❌ POOR: Limited subgroup separation\n")
  cat("Recommendation: Consider alternative approaches or skip clustering-based analysis\n")
} else {
  cat("\n❌ FAILED: No valid clustering analyses\n")
  cat("Recommendation: Review sample filtering and data quality\n")
}

# Next steps guidance
cat("\nNext steps:\n")
if (high_purity_count > 0 || moderate_purity_count > 0) {
  cat("1. Proceed to DEGES normalization (08)\n")
  cat("2. Focus DEG analysis on groups with high separation\n")
  cat("3. Consider cluster-informed sample stratification\n")
  cat("4. Use clustering results to interpret DEG findings\n")
} else {
  cat("1. Review high-purity filtering thresholds in 06\n")
  cat("2. Consider standard group comparison without clustering\n")
  cat("3. May proceed with DEG analysis using original R0/R1 labels\n")
}

# Special note about R0+R1 clustering
r01_tumor_available <- "R0_R1_tumor" %in% names(clustering_results) && 
  !is.null(clustering_results[["R0_R1_tumor"]])
r01_normal_available <- "R0_R1_normal" %in% names(clustering_results) && 
  !is.null(clustering_results[["R0_R1_normal"]])

if (r01_tumor_available || r01_normal_available) {
  cat("\n🎯 R0+R1 Clustering Results Available:\n")
  cat("These results can inform radiation exposure hypothesis testing\n")
  cat("Consider using clustering-based stratification in feature selection\n")
} else {
  cat("\n⚠️ R0+R1 Clustering Not Available:\n")
  cat("Limited samples for radiation exposure clustering analysis\n")
}

cat("\nPhase 2 comprehensive analysis (v3 corrected) completed!\n")
cat("==============================================\n")

