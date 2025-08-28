# ==============================================================================
# Gene Count Optimization for CDM Analysis
# 05b_gene_count_optimization.R
# ==============================================================================

# Required libraries
library(dplyr)
library(edgeR)

# Source improved utilities
source("./utils/utils_improved.R")
source("./utils/norm_improved.R")

cat("Starting gene count optimization analysis...\n")

# ==============================================================================
# 1. Load Previous Results and Data
# ==============================================================================

load("./data/processed/sample_lists.rda")
load("./data/raw/thyr_data.rda")
se_thyr <- data
count_data <- assay(se_thyr, "stranded_second")
rm(data, se_thyr)

# Apply same filtering as v5
dgelist_all <- DGEList(counts = count_data)
all_samples <- unique(unlist(lapply(sample_lists, function(x) c(x$tumor, x$normal))))
available_samples <- intersect(all_samples, colnames(count_data))
dgelist_filtered <- dgelist_all[, available_samples]
keep <- filterByExpr(dgelist_filtered)
dgelist_filtered <- dgelist_filtered[keep, , keep.lib.sizes = FALSE]
filtered_counts <- dgelist_filtered$counts

cat("Base filtered count data:", dim(filtered_counts), "\n")

# ==============================================================================
# 2. Test Gene Counts
# ==============================================================================

gene_counts <- c(5000, 8000, 10000, 12000)
test_groups <- c("R0_tumor", "R1_tumor", "B0_normal", "B1_normal")  # Representative groups

# ==============================================================================
# 3. Modified Functions for Testing
# ==============================================================================

prepare_cdm_data_test <- function(count_matrix, sample_ids, analysis_name, max_genes) {
  analysis_counts <- count_matrix[, sample_ids, drop = FALSE]
  
  initialize_muren(n_genes = nrow(analysis_counts), 
                   n_samples = ncol(analysis_counts), 
                   priority = "robustness")
  
  tryCatch({
    muren_coeff <- muren_norm(
      reads = analysis_counts,
      refs = 'saturated',
      pairwise_method = "lts",
      single_param = TRUE,
      res_return = 'scaling_coeff',
      filter_gene = FALSE,
      trim = 10,
      workers = "auto"
    )
    
    normalized_counts <- sweep(analysis_counts, 2, 1/muren_coeff, "*")
    
  }, error = function(e) {
    dgelist_temp <- DGEList(counts = analysis_counts)
    dgelist_temp <- calcNormFactors(dgelist_temp)
    normalized_counts <- cpm(dgelist_temp, normalized.lib.sizes = TRUE)
  })
  
  total_counts <- colSums(normalized_counts)
  cpm_data <- sweep(normalized_counts, 2, total_counts, "/") * 1e6
  prior_count <- 0.5
  log_cpm <- log2(cpm_data + prior_count)
  
  # Gene selection by variance
  if (nrow(log_cpm) > max_genes) {
    gene_vars <- apply(log_cpm, 1, var, na.rm = TRUE)
    top_genes <- order(gene_vars, decreasing = TRUE)[1:max_genes]
    log_cpm <- log_cpm[top_genes, , drop = FALSE]
  }
  
  return(log_cpm)
}

perform_cdm_pca_test <- function(log_cpm_data) {
  X <- log_cpm_data
  X_centered <- X - rowMeans(X)
  gene_sds <- apply(X_centered, 1, sd)
  
  non_zero_var <- gene_sds > .Machine$double.eps
  if (sum(non_zero_var) < nrow(X)) {
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
  variance_explained <- eigenvalues / sum(eigenvalues)
  
  max_components <- min(5, ncol(X_scaled) - 1, nrow(X_scaled))
  if (ncol(scores) > max_components) {
    scores <- scores[, 1:max_components, drop = FALSE]
    variance_explained <- variance_explained[1:max_components]
  }
  
  rownames(scores) <- colnames(X)
  colnames(scores) <- paste0("PC", 1:ncol(scores))
  
  return(list(
    scores = scores,
    variance_explained = variance_explained,
    n_components = ncol(scores)
  ))
}

detect_outliers_test <- function(pca_result) {
  scores <- pca_result$scores
  if (is.null(scores) || nrow(scores) < 4) return(0)
  
  n_pcs_for_outlier <- min(3, ncol(scores))
  pc_subset <- scores[, 1:n_pcs_for_outlier, drop = FALSE]
  
  center_point <- colMeans(pc_subset)
  cov_matrix <- cov(pc_subset)
  
  mahal_outliers <- integer(0)
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
  
  return(length(combined_outliers))
}

# ==============================================================================
# 4. Analysis Groups Setup
# ==============================================================================

analysis_groups <- list(
  "R0_normal" = list(group = "R0", tissue = "normal"),
  "R0_tumor" = list(group = "R0", tissue = "tumor"),
  "R1_normal" = list(group = "R1", tissue = "normal"),
  "R1_tumor" = list(group = "R1", tissue = "tumor"),
  "B0_normal" = list(group = "B0", tissue = "normal"),
  "B0_tumor" = list(group = "B0", tissue = "tumor"),
  "B1_normal" = list(group = "B1", tissue = "normal"),
  "B1_tumor" = list(group = "B1", tissue = "tumor")
)

# ==============================================================================
# 5. Run Optimization Tests
# ==============================================================================

cat("\nRunning gene count optimization tests with LTS normalization...\n")

results_summary <- data.frame(
  Gene_Count = integer(0),
  Group = character(0),
  Tissue = character(0),
  N_Samples = integer(0),
  PC1_Var = numeric(0),
  PC2_Var = numeric(0),
  PC1_PC2_Sum = numeric(0),
  N_Outliers = integer(0),
  Outlier_Rate = numeric(0),
  PC1_Category = character(0),
  stringsAsFactors = FALSE
)

for (max_genes in gene_counts) {
  cat(sprintf("\n=== Testing with %d genes ===\n", max_genes))
  
  for (analysis_name in names(analysis_groups)) {
    group_info <- analysis_groups[[analysis_name]]
    group_name <- group_info$group
    tissue_type <- group_info$tissue
    
    if (!group_name %in% names(sample_lists)) next
    
    if (tissue_type == "tumor") {
      sample_ids <- sample_lists[[group_name]]$tumor
    } else {
      sample_ids <- sample_lists[[group_name]]$normal
    }
    
    if (length(sample_ids) == 0) next
    available_samples <- sample_ids[sample_ids %in% colnames(filtered_counts)]
    if (length(available_samples) < 3) next
    
    tryCatch({
      cdm_data <- prepare_cdm_data_test(filtered_counts, available_samples, analysis_name, max_genes)
      pca_result <- perform_cdm_pca_test(cdm_data)
      n_outliers <- detect_outliers_test(pca_result)
      
      pc1_var <- pca_result$variance_explained[1] * 100
      pc2_var <- ifelse(length(pca_result$variance_explained) > 1, 
                        pca_result$variance_explained[2] * 100, 0)
      pc1_pc2_sum <- pc1_var + pc2_var
      outlier_rate <- n_outliers / length(available_samples) * 100
      
      # PC1 category
      pc1_category <- if (pc1_var >= 90) {
        "Very_High"
      } else if (pc1_var >= 80) {
        "High"
      } else if (pc1_var >= 50) {
        "Optimal"
      } else {
        "Low"
      }
      
      results_summary <- rbind(results_summary, data.frame(
        Gene_Count = max_genes,
        Group = group_name,
        Tissue = tissue_type,
        N_Samples = length(available_samples),
        PC1_Var = round(pc1_var, 1),
        PC2_Var = round(pc2_var, 1),
        PC1_PC2_Sum = round(pc1_pc2_sum, 1),
        N_Outliers = n_outliers,
        Outlier_Rate = round(outlier_rate, 1),
        PC1_Category = pc1_category,
        stringsAsFactors = FALSE
      ))
      
      cat(sprintf("%s %s: PC1=%.1f%%, Outliers=%d (%.1f%%)\n", 
                  group_name, tissue_type, pc1_var, n_outliers, outlier_rate))
      
    }, error = function(e) {
      cat(sprintf("Error in %s %s: %s\n", group_name, tissue_type, e$message))
    })
  }
}

# ==============================================================================
# 6. Analysis and Recommendations
# ==============================================================================

cat("\n=== Gene Count Optimization Results ===\n")

# Summary by gene count
gene_summary <- results_summary %>%
  group_by(Gene_Count) %>%
  summarise(
    Mean_PC1_Var = round(mean(PC1_Var), 1),
    Mean_Outlier_Rate = round(mean(Outlier_Rate), 1),
    N_Very_High_PC1 = sum(PC1_Category == "Very_High"),
    N_High_PC1 = sum(PC1_Category == "High"),
    N_Optimal_PC1 = sum(PC1_Category == "Optimal"),
    N_Low_PC1 = sum(PC1_Category == "Low"),
    .groups = "drop"
  )

print(gene_summary)

cat("\nDetailed results by gene count:\n")
for (gc in gene_counts) {
  cat(sprintf("\n--- %d Genes ---\n", gc))
  subset_data <- results_summary[results_summary$Gene_Count == gc, ]
  print(subset_data[, c("Group", "Tissue", "PC1_Var", "N_Outliers", "Outlier_Rate", "PC1_Category")])
}

# Recommendations
cat("\n=== Recommendations ===\n")

optimal_criteria <- gene_summary %>%
  mutate(
    Score = (N_Optimal_PC1 * 3) + (N_High_PC1 * 1) - (N_Very_High_PC1 * 2) - (N_Low_PC1 * 1),
    Avg_Outlier_Penalty = ifelse(Mean_Outlier_Rate > 20, -2, 
                                 ifelse(Mean_Outlier_Rate > 15, -1, 0)),
    Final_Score = Score + Avg_Outlier_Penalty
  ) %>%
  arrange(desc(Final_Score))

cat("Gene count ranking (higher score = better):\n")
print(optimal_criteria[, c("Gene_Count", "Mean_PC1_Var", "Mean_Outlier_Rate", 
                           "N_Optimal_PC1", "N_Very_High_PC1", "Final_Score")])

best_gene_count <- optimal_criteria$Gene_Count[1]
cat(sprintf("\nRecommended gene count: %d\n", best_gene_count))

# Quality flags
quality_flags <- results_summary %>%
  mutate(
    Quality_Flag = case_when(
      PC1_Var >= 90 ~ "PC1_Too_High",
      Outlier_Rate >= 25 ~ "Too_Many_Outliers",
      PC1_Var < 40 ~ "PC1_Too_Low",
      TRUE ~ "Acceptable"
    )
  ) %>%
  filter(Quality_Flag != "Acceptable")

if (nrow(quality_flags) > 0) {
  cat("\nQuality concerns:\n")
  print(quality_flags[, c("Gene_Count", "Group", "Tissue", "PC1_Var", "Outlier_Rate", "Quality_Flag")])
}

# ==============================================================================
# 7. Save Results
# ==============================================================================

cat("\nSaving optimization results...\n")

gene_count_optimization <- list(
  results_summary = results_summary,
  gene_summary = gene_summary,
  optimal_criteria = optimal_criteria,
  recommended_gene_count = best_gene_count,
  quality_flags = quality_flags,
  test_parameters = list(
    gene_counts = gene_counts,
    analysis_date = Sys.time(),
    criteria = "PC1 variance 50-80%, outlier rate <20%"
  )
)

save(gene_count_optimization, file = "./data/processed/gene_count_optimization.rda")
write.csv(results_summary, file = "./output/reports/gene_count_optimization_detailed.csv", row.names = FALSE)
write.csv(gene_summary, file = "./output/reports/gene_count_optimization_summary.csv", row.names = FALSE)

cat("Results saved to ./data/processed/ and ./output/reports/\n")

# ==============================================================================
# 8. Final Summary
# ==============================================================================

cat(sprintf("\n=== Gene Count Optimization Complete ===\n"))
cat(sprintf("Recommended gene count: %d genes\n", best_gene_count))
cat(sprintf("This provides the best balance of:\n"))
cat(sprintf("- PC1 variance in optimal range (50-80%%)\n"))
cat(sprintf("- Reasonable outlier detection rate (<20%%)\n"))
cat(sprintf("- Minimal overfitting indicators\n"))

best_row <- optimal_criteria[1, ]
cat(sprintf("\nAt %d genes:\n", best_gene_count))
cat(sprintf("- Average PC1 variance: %.1f%%\n", best_row$Mean_PC1_Var))
cat(sprintf("- Average outlier rate: %.1f%%\n", best_row$Mean_Outlier_Rate))
cat(sprintf("- Groups in optimal PC1 range: %d/8\n", best_row$N_Optimal_PC1))

cat("\nRerun 05_pca_analysis_v5.R with the recommended gene count for optimal results.\n")