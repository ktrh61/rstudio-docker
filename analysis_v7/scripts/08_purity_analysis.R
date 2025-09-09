# 08_purity_analysis.R - Stage 2 Tumor Purity Estimation and Filtering
# Purpose: Estimate tumor purity using ContamDE for R0/R1/B0/B1 groups
# Method: MUREN normalization with contamDE purity estimation

source("analysis_v7/setup.R")

cat("\n=== Stage 2: Tumor Purity Analysis ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# Load packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
  library(qvalue)
  library(dplyr)
})

# Load improved functions
source("utils/utils_improved.R")
source("utils/norm_improved.R")
source("utils/contamde_purity_functions.R")

# Configuration
CONFIG <- list(
  PURITY_THRESHOLD = 0.6,      # 60% minimum purity
  PAIRWISE_METHOD = "lts",     # MUREN method: lts, median, trim10
  WORKERS = "auto",             # Parallel workers
  PRIOR_COUNT = 0,              # No pseudocount by default
  VERBOSE = TRUE
)

# Thread control
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  cat("BLAS/OMP threads set to 1\n")
}

# Load Stage 1 filtered data
cat("\nLoading Stage 1 filtered data...\n")
se <- readRDS(paste0(paths$processed, "thyr_se_stage1_filtered.rds"))
metadata <- as.data.frame(colData(se))

cat("Input samples:", ncol(se), "\n")
cat("Input genes:", nrow(se), "\n")

# Extract paired sample groups for purity estimation
extract_paired_groups <- function(metadata, se) {
  groups <- list()
  
  for (grp in c("R0", "R1", "B0", "B1")) {
    # Find paired samples in this group
    grp_meta <- metadata[metadata$group == grp & !is.na(metadata$group), ]
    
    # Get cases with both tumor and normal
    cases_with_pairs <- grp_meta %>%
      dplyr::group_by(case_id) %>%
      dplyr::summarise(
        has_tumor = any(is_tumor),
        has_normal = any(is_normal),
        n_tumor = sum(is_tumor),
        n_normal = sum(is_normal),
        .groups = "drop"
      ) %>%
      dplyr::filter(has_tumor & has_normal & n_tumor == n_normal)
    
    if (nrow(cases_with_pairs) > 0) {
      # Sort cases for consistent ordering
      cases_ordered <- sort(cases_with_pairs$case_id)
      
      # Get sample indices maintaining case order
      tumor_idx <- numeric(length(cases_ordered))
      normal_idx <- numeric(length(cases_ordered))
      
      for (i in seq_along(cases_ordered)) {
        case <- cases_ordered[i]
        tumor_idx[i] <- which(metadata$case_id == case & 
                                metadata$group == grp & 
                                metadata$is_tumor)[1]  # Take first if multiple
        normal_idx[i] <- which(metadata$case_id == case & 
                                 metadata$group == grp & 
                                 metadata$is_normal)[1]  # Take first if multiple
      }
      
      # Remove any NA indices (shouldn't happen but safe guard)
      valid <- !is.na(tumor_idx) & !is.na(normal_idx)
      tumor_idx <- tumor_idx[valid]
      normal_idx <- normal_idx[valid]
      cases_ordered <- cases_ordered[valid]
      
      groups[[grp]] <- list(
        tumor_idx = tumor_idx,
        normal_idx = normal_idx,
        tumor_ids = metadata$sample_id[tumor_idx],
        normal_ids = metadata$sample_id[normal_idx],
        cases = cases_ordered  # Now guaranteed to match index order
      )
    }
  }
  
  groups
}

paired_groups <- extract_paired_groups(metadata, se)

cat("\n=== Paired Sample Groups ===\n")
for (grp in names(paired_groups)) {
  cat(sprintf("%s: %d pairs\n", grp, length(paired_groups[[grp]]$cases)))
}

# Prepare counts matrix for ContamDE (normal first, then tumor)
prepare_contamde_matrix <- function(se, normal_idx, tumor_idx) {
  # Extract counts using stranded_second assay
  counts <- cbind(
    assay(se, "stranded_second")[, normal_idx, drop = FALSE],
    assay(se, "stranded_second")[, tumor_idx, drop = FALSE]
  )
  
  # Set informative column names
  n_pairs <- length(normal_idx)
  colnames(counts) <- c(
    paste0("Normal_", seq_len(n_pairs)),
    paste0("Tumor_", seq_len(n_pairs))
  )
  
  return(counts)
}

# Run purity estimation for each group
purity_results <- list()
purity_summaries <- list()

for (grp in names(paired_groups)) {
  cat(sprintf("\n--- Processing %s ---\n", grp))
  
  # Skip if no pairs
  if (length(paired_groups[[grp]]$cases) == 0) {
    cat("  No paired samples available\n")
    next
  }
  
  # Prepare counts matrix
  counts <- prepare_contamde_matrix(
    se,
    paired_groups[[grp]]$normal_idx,
    paired_groups[[grp]]$tumor_idx
  )
  
  cat(sprintf("  Count matrix: %d genes x %d samples\n", 
              nrow(counts), ncol(counts)))
  
  # Run ContamDE purity estimation
  tryCatch({
    purity_result <- contamde_purity(
      counts = counts,
      subtype = NULL,
      covariate = NULL,
      contaminated = TRUE,
      pairwise_method = CONFIG$PAIRWISE_METHOD,
      workers = CONFIG$WORKERS,
      prior.count = CONFIG$PRIOR_COUNT,
      verbose = CONFIG$VERBOSE
    )
    
    purity_results[[grp]] <- purity_result
    
    # Quality assessment
    quality_stats <- assess_purity_quality(
      purity_result,
      threshold = CONFIG$PURITY_THRESHOLD,
      verbose = TRUE
    )
    
    purity_summaries[[grp]] <- quality_stats
    
  }, error = function(e) {
    cat("  Error in purity estimation:", e$message, "\n")
    purity_results[[grp]] <- NULL
    purity_summaries[[grp]] <- NULL
  })
}

# Create filtered sample lists based on purity
cat("\n=== Creating Purity-Filtered Sample Lists ===\n")

filtered_groups <- list()
for (grp in names(purity_results)) {
  if (is.null(purity_results[[grp]])) {
    cat(sprintf("  %s: Skipped (no results)\n", grp))
    next
  }
  
  purity <- purity_results[[grp]]$proportion
  high_purity <- purity >= CONFIG$PURITY_THRESHOLD
  high_purity_idx <- which(high_purity)
  
  if (length(high_purity_idx) == 0) {
    cat(sprintf("  %s: No high purity samples\n", grp))
    filtered_groups[[grp]] <- list(
      tumor_idx = integer(0),
      normal_idx = integer(0),
      tumor_ids = character(0),
      normal_ids = character(0),
      cases = character(0)
    )
  } else {
    # Apply filter to indices
    orig <- paired_groups[[grp]]
    filtered_groups[[grp]] <- list(
      tumor_idx = orig$tumor_idx[high_purity_idx],
      normal_idx = orig$normal_idx[high_purity_idx],
      tumor_ids = orig$tumor_ids[high_purity_idx],
      normal_ids = orig$normal_ids[high_purity_idx],
      cases = orig$cases[high_purity_idx],
      purity_scores = purity[high_purity_idx]
    )
    
    cat(sprintf("  %s: %d → %d pairs (%.1f%% retention)\n",
                grp,
                length(orig$cases),
                length(filtered_groups[[grp]]$cases),
                length(filtered_groups[[grp]]$cases) / length(orig$cases) * 100))
  }
}

# Summary statistics
cat("\n=== Overall Summary ===\n")

summary_df <- data.frame(
  Group = names(paired_groups),
  Original_Pairs = sapply(paired_groups, function(x) length(x$cases)),
  Retained_Pairs = sapply(names(paired_groups), function(g) {
    if (g %in% names(filtered_groups)) length(filtered_groups[[g]]$cases) else 0
  }),
  stringsAsFactors = FALSE
)

summary_df$Retention_Rate <- ifelse(
  summary_df$Original_Pairs > 0,
  summary_df$Retained_Pairs / summary_df$Original_Pairs * 100,
  NA
)
summary_df$Mean_Purity <- sapply(names(paired_groups), function(g) {
  if (g %in% names(purity_results) && !is.null(purity_results[[g]])) {
    mean(purity_results[[g]]$proportion)
  } else NA
})

print(summary_df)

total_original <- sum(summary_df$Original_Pairs)
total_retained <- sum(summary_df$Retained_Pairs)

if (total_original > 0) {
  overall_retention <- total_retained / total_original * 100
  cat(sprintf("\nTotal: %d → %d pairs (%.1f%% overall retention)\n",
              total_original, total_retained, overall_retention))
} else {
  cat("\nNo paired samples found in any group\n")
}

# Create filtered SummarizedExperiment
all_retained_idx <- sort(unique(c(
  unlist(lapply(filtered_groups, function(x) x$tumor_idx)),
  unlist(lapply(filtered_groups, function(x) x$normal_idx))
)))

if (length(all_retained_idx) > 0) {
  se_purity_filtered <- se[, all_retained_idx]
  cat(sprintf("\nFinal SE: %d genes x %d samples\n",
              nrow(se_purity_filtered), ncol(se_purity_filtered)))
} else {
  se_purity_filtered <- se[, integer(0)]
  cat("\nWarning: No samples retained after purity filtering\n")
}

# Save results
cat("\nSaving results...\n")

# Main outputs
saveRDS(purity_results, paste0(paths$output, "stage2_purity_results.rds"))
saveRDS(filtered_groups, paste0(paths$output, "stage2_filtered_groups.rds"))
saveRDS(se_purity_filtered, paste0(paths$processed, "thyr_se_stage2_filtered.rds"))

# Purity scores table
purity_scores_list <- list()
for (grp in names(purity_results)) {
  if (!is.null(purity_results[[grp]])) {
    purity_scores_list[[grp]] <- data.frame(
      group = grp,
      case_id = paired_groups[[grp]]$cases,
      tumor_id = paired_groups[[grp]]$tumor_ids,
      normal_id = paired_groups[[grp]]$normal_ids,
      purity = purity_results[[grp]]$proportion,
      retained = purity_results[[grp]]$proportion >= CONFIG$PURITY_THRESHOLD,
      stringsAsFactors = FALSE
    )
  }
}

if (length(purity_scores_list) > 0) {
  purity_scores_df <- do.call(rbind, purity_scores_list)
  write.csv(purity_scores_df, paste0(paths$output, "stage2_purity_scores.csv"), row.names = FALSE)
} else {
  cat("  No purity scores to save\n")
}

# Summary statistics
write.csv(summary_df, paste0(paths$output, "stage2_summary_stats.csv"), row.names = FALSE)

# List of excluded samples
excluded_cases <- list()
for (grp in names(paired_groups)) {
  if (grp %in% names(filtered_groups)) {
    excluded <- setdiff(paired_groups[[grp]]$cases, filtered_groups[[grp]]$cases)
    if (length(excluded) > 0) {
      excluded_cases[[grp]] <- excluded
    }
  }
}

if (length(excluded_cases) > 0) {
  excluded_df <- data.frame(
    group = rep(names(excluded_cases), sapply(excluded_cases, length)),
    case_id = unlist(excluded_cases),
    stage = "Stage2_Purity",
    stringsAsFactors = FALSE
  )
  write.csv(excluded_df, paste0(paths$output, "stage2_excluded_cases.csv"), row.names = FALSE)
}

cat("\n=== Stage 2 Complete ===\n")
cat("Files created:\n")
cat("  - stage2_purity_results.rds\n")
cat("  - stage2_filtered_groups.rds\n")
cat("  - thyr_se_stage2_filtered.rds\n")
if (length(purity_scores_list) > 0) {
  cat("  - stage2_purity_scores.csv\n")
}
cat("  - stage2_summary_stats.csv\n")
if (length(excluded_cases) > 0) {
  cat("  - stage2_excluded_cases.csv\n")
}
cat("\nConfiguration:\n")
cat(sprintf("  Purity threshold: %.1f%%\n", CONFIG$PURITY_THRESHOLD * 100))
cat(sprintf("  MUREN method: %s\n", CONFIG$PAIRWISE_METHOD))
cat(sprintf("  Prior count: %g\n", CONFIG$PRIOR_COUNT))
cat("\nNext: Run 09_cooks_distance_filtering.R for Stage 3\n")