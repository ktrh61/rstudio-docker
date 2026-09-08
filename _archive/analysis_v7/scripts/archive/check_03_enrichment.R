# check_03_enrichment.R
# Purpose: Extract enrichment analysis results for manuscript Results section
# Section 3: Enrichment Analysis with GSEA Cross-comparison
# Version: v7.3 - Updated for GSEA cross-comparison analysis
# Date: 2025-10-16

# Setup
source("analysis_v7/setup.R")

cat("=====================================\n")
cat("ENRICHMENT ANALYSIS RESULTS CHECK\n")
cat("=====================================\n\n")

# Load packages
suppressPackageStartupMessages({
  library(dplyr)
})

# ============================================================================
# Load enrichment results
# ============================================================================

cat("--- Loading enrichment results ---\n")

# Check for results file
enrichment_file <- paste0(paths$output, "enrichment_analysis/all_enrichment_results.rds")

if (!file.exists(enrichment_file)) {
  stop("Enrichment results not found. Please run 10_enrichment_analysis.R first.")
}

all_enrichment_results <- readRDS(enrichment_file)
cat("✓ Enrichment results loaded\n")

# Check for cross-driver annotation (now in separate file)
cross_driver_file <- paste0(paths$output, "cross_driver_analysis/cross_driver_genes_annotated.rds")
if (file.exists(cross_driver_file)) {
  cross_driver_annotated <- readRDS(cross_driver_file)
  cat("✓ Cross-driver annotation loaded\n")
} else {
  cat("✗ Cross-driver annotation not found (run separate cross-driver script if needed)\n")
}

# ============================================================================
# Extract ORA results (R0 vs R1 tumor only)
# ============================================================================

cat("\n=====================================\n")
cat("ORA ANALYSIS - R0 vs R1 TUMOR\n")
cat("=====================================\n")

if ("R0_vs_R1_tumor_ORA" %in% names(all_enrichment_results)) {
  ora_results <- all_enrichment_results$R0_vs_R1_tumor_ORA
  
  # Process each direction and database
  for (direction in c("up", "down", "all")) {
    cat(sprintf("\n--- %s-regulated genes ---\n", toupper(direction)))
    
    for (db in c("GO_BP", "KEGG", "Reactome")) {
      result_key <- paste0(db, "_", direction)
      
      if (result_key %in% names(ora_results) && !is.null(ora_results[[result_key]])) {
        result <- ora_results[[result_key]]
        
        if (nrow(result@result) > 0) {
          # Count significant terms (q < 0.05)
          sig_terms <- result@result[result@result$qvalue < 0.05, ]
          cat(sprintf("□ %s: %d enriched terms\n", db, nrow(sig_terms)))
          
          # Show top 3 terms
          if (nrow(sig_terms) > 0) {
            top_terms <- head(sig_terms[order(sig_terms$qvalue), ], 3)
            for (i in 1:nrow(top_terms)) {
              cat(sprintf("  %d. %s (q=%g, n=%s)\n",
                          i, 
                          substr(top_terms$Description[i], 1, 50),
                          top_terms$qvalue[i],
                          top_terms$Count[i]))
            }
          }
        } else {
          cat(sprintf("□ %s: No enriched terms\n", db))
        }
      }
    }
  }
}

# ============================================================================
# Extract GSEA results (all comparisons)
# ============================================================================

cat("\n=====================================\n")
cat("GSEA ANALYSIS - ALL COMPARISONS\n")
cat("=====================================\n")

if ("GSEA" %in% names(all_enrichment_results)) {
  gsea_results <- all_enrichment_results$GSEA
  
  for (comp_name in names(gsea_results)) {
    cat(sprintf("\n--- %s ---\n", comp_name))
    
    comp_results <- gsea_results[[comp_name]]
    
    for (db in c("GO_BP", "KEGG", "Reactome")) {
      if (db %in% names(comp_results) && !is.null(comp_results[[db]])) {
        result <- comp_results[[db]]
        
        if (nrow(result@result) > 0) {
          # Count significant terms
          sig_terms <- result@result[result@result$qvalue < 0.05, ]
          cat(sprintf("□ %s: %d enriched terms\n", db, nrow(sig_terms)))
          
          # Count up/down enriched
          if (nrow(sig_terms) > 0) {
            n_up <- sum(sig_terms$NES > 0)
            n_down <- sum(sig_terms$NES < 0)
            cat(sprintf("  Up-enriched: %d, Down-enriched: %d\n", n_up, n_down))
            
            # Top enriched pathway
            top_term <- sig_terms[order(abs(sig_terms$NES), decreasing = TRUE)[1], ]
            direction_sym <- ifelse(top_term$NES > 0, "↑", "↓")
            cat(sprintf("  Top: %s %s (NES=%g, q=%g)\n", 
                        direction_sym,
                        substr(top_term$Description, 1, 40),
                        top_term$NES,
                        top_term$qvalue))
          }
        } else {
          cat(sprintf("□ %s: No enriched terms\n", db))
        }
      }
    }
  }
}

# ============================================================================
# GSEA Cross-comparison Analysis
# ============================================================================

cat("\n=====================================\n")
cat("GSEA CROSS-COMPARISON ANALYSIS\n")
cat("=====================================\n")

if ("cross_comparison" %in% names(all_enrichment_results)) {
  cross_comp <- all_enrichment_results$cross_comparison
  
  cat("\n--- Summary ---\n")
  
  # Count consistent pathways for each database
  databases <- c("GO_BP", "KEGG", "Reactome")
  
  for (db in databases) {
    if (!is.null(cross_comp[[db]])) {
      n_total <- nrow(cross_comp[[db]]$all)
      n_consistent <- nrow(cross_comp[[db]]$consistent)
      n_inconsistent <- nrow(cross_comp[[db]]$inconsistent)
      
      cat(sprintf("□ %s:\n", db))
      cat(sprintf("  Common pathways: %d\n", n_total))
      cat(sprintf("  NES-consistent: %d (%.1f%%)\n", 
                  n_consistent, 
                  100 * n_consistent / max(n_total, 1)))
      cat(sprintf("  NES-inconsistent: %d\n", n_inconsistent))
      
      # Show top 3 consistent pathways
      if (n_consistent > 0) {
        cat("  Top consistent pathways:\n")
        top3 <- head(cross_comp[[db]]$consistent, 3)
        for (i in 1:nrow(top3)) {
          cat(sprintf("    %d. %s\n", i, substr(top3$Description[i], 1, 50)))
          cat(sprintf("       NES_RET=%g, NES_BRAF=%g, Score=%g\n",
                      top3$NES_RET[i], top3$NES_BRAF[i], 
                      top3$consistency_score[i]))
        }
      }
    } else {
      cat(sprintf("□ %s: No cross-comparison results\n", db))
    }
  }
  
  # Overall statistics
  cat("\n--- Overall Statistics ---\n")
  
  total_consistent <- sum(sapply(databases, function(db) {
    if (!is.null(cross_comp[[db]])) {
      nrow(cross_comp[[db]]$consistent)
    } else 0
  }))
  
  cat(sprintf("Total consistent pathways across all databases: %d\n", total_consistent))
  
  # Statistical details for consistent pathways
  if (total_consistent > 0) {
    all_consistent_for_stats <- do.call(rbind, lapply(databases, function(db) {
      if (!is.null(cross_comp[[db]]) && nrow(cross_comp[[db]]$consistent) > 0) {
        cross_comp[[db]]$consistent
      } else NULL
    }))
    
    if (!is.null(all_consistent_for_stats)) {
      cat("\n--- Statistical Details ---\n")
      cat(sprintf("  Median q-value (RET): %.2e\n", median(all_consistent_for_stats$qvalue_RET)))
      cat(sprintf("  Median q-value (BRAF): %.2e\n", median(all_consistent_for_stats$qvalue_BRAF)))
      cat(sprintf("  Mean |NES|: %.2f\n", mean(abs(c(all_consistent_for_stats$NES_RET, 
                                                     all_consistent_for_stats$NES_BRAF)))))
      
      # NES correlation
      if (nrow(all_consistent_for_stats) > 1) {
        cor_value <- cor(all_consistent_for_stats$NES_RET, all_consistent_for_stats$NES_BRAF)
        cat(sprintf("  NES correlation: r = %.3f\n", cor_value))
      }
      
      # Check for pathways with large NES differences
      all_consistent_for_stats$NES_diff <- abs(all_consistent_for_stats$NES_RET - 
                                                 all_consistent_for_stats$NES_BRAF)
      high_diff <- all_consistent_for_stats[all_consistent_for_stats$NES_diff > 0.5, ]
      if (nrow(high_diff) > 0) {
        cat(sprintf("\n  Pathways with NES difference >0.5: %d\n", nrow(high_diff)))
        if (nrow(high_diff) <= 3) {
          for (i in 1:nrow(high_diff)) {
            cat(sprintf("    - %s (diff=%.2f)\n", 
                        substr(high_diff$Description[i], 1, 50),
                        high_diff$NES_diff[i]))
          }
        }
      }
    }
  }
  
  # Find most consistent pathways across all databases
  all_consistent_pathways <- do.call(rbind, lapply(databases, function(db) {
    if (!is.null(cross_comp[[db]]) && nrow(cross_comp[[db]]$consistent) > 0) {
      cross_comp[[db]]$consistent
    } else NULL
  }))
  
  if (!is.null(all_consistent_pathways) && nrow(all_consistent_pathways) > 0) {
    cat("\n--- Top 5 Most Consistent Pathways (All DBs) ---\n")
    all_consistent_pathways <- all_consistent_pathways[order(
      all_consistent_pathways$consistency_score, decreasing = TRUE), ]
    top5 <- head(all_consistent_pathways, 5)
    
    for (i in 1:nrow(top5)) {
      cat(sprintf("□ [%s] %s\n", 
                  top5$database[i],
                  substr(top5$Description[i], 1, 60)))
      cat(sprintf("  NES: RET=%g, BRAF=%g | Score=%g\n",
                  top5$NES_RET[i], top5$NES_BRAF[i], 
                  top5$consistency_score[i]))
    }
  }
  
} else {
  cat("No cross-comparison results found in enrichment results.\n")
  cat("This analysis compares GSEA results between R0_vs_R1_tumor and B0_vs_B1_tumor.\n")
}

# ============================================================================
# Cross-driver Gene Annotation (if available)
# ============================================================================

if (exists("cross_driver_annotated") && !is.null(cross_driver_annotated)) {
  cat("\n=====================================\n")
  cat("CROSS-DRIVER GENE ANNOTATION\n")
  cat("=====================================\n")
  
  cat(sprintf("Total cross-driver genes: %d\n", nrow(cross_driver_annotated)))
  
  # Direction analysis
  if ("direction_match" %in% colnames(cross_driver_annotated)) {
    n_consistent <- sum(cross_driver_annotated$direction_match)
    n_inconsistent <- sum(!cross_driver_annotated$direction_match)
    cat(sprintf("Direction consistent: %d (%.1f%%)\n", 
                n_consistent, 100 * n_consistent / nrow(cross_driver_annotated)))
    cat(sprintf("Direction inconsistent: %d\n", n_inconsistent))
  }
  
  # Gene type distribution
  if ("GENETYPE" %in% colnames(cross_driver_annotated)) {
    cat("\n--- Gene type distribution ---\n")
    type_table <- table(cross_driver_annotated$GENETYPE)
    type_df <- data.frame(
      Type = names(type_table),
      Count = as.vector(type_table),
      Percent = round(100 * as.vector(type_table) / sum(type_table), 1)
    )
    type_df <- type_df[order(type_df$Count, decreasing = TRUE), ]
    for (i in 1:min(5, nrow(type_df))) {
      cat(sprintf("□ %s: %d genes (%.1f%%)\n",
                  type_df$Type[i], type_df$Count[i], type_df$Percent[i]))
    }
  }
  
  # Top genes by consistency
  cat("\n--- Top 5 genes by consistency score ---\n")
  if ("consistency_score" %in% colnames(cross_driver_annotated)) {
    top5 <- head(cross_driver_annotated[order(
      cross_driver_annotated$consistency_score, decreasing = TRUE), ], 5)
    for (i in 1:nrow(top5)) {
      cat(sprintf("□ %s: FC_RET=%.2f, FC_BRAF=%.2f, Score=%.3f\n",
                  top5$SYMBOL[i],
                  top5$FC_RET[i],
                  top5$FC_BRAF[i],
                  top5$consistency_score[i]))
    }
  }
}

# ============================================================================
# Summary for manuscript
# ============================================================================

cat("\n=====================================\n")
cat("SUMMARY FOR MANUSCRIPT\n")
cat("=====================================\n")

cat("\nSection 3 key findings:\n")
cat("----------------------------------\n")

# Count total enriched pathways
total_ora <- 0
total_gsea <- 0

if ("R0_vs_R1_tumor_ORA" %in% names(all_enrichment_results)) {
  for (result in all_enrichment_results$R0_vs_R1_tumor_ORA) {
    if (!is.null(result) && nrow(result@result) > 0) {
      total_ora <- total_ora + sum(result@result$qvalue < 0.05)
    }
  }
}

if ("GSEA" %in% names(all_enrichment_results)) {
  for (comp in all_enrichment_results$GSEA) {
    for (result in comp) {
      if (!is.null(result) && nrow(result@result) > 0) {
        total_gsea <- total_gsea + sum(result@result$qvalue < 0.05)
      }
    }
  }
}

cat(sprintf("\nTotal enriched pathways:\n"))
cat(sprintf("  ORA (R0 vs R1 tumor): %d\n", total_ora))
cat(sprintf("  GSEA (all comparisons): %d\n", total_gsea))

# Cross-comparison summary
if ("cross_comparison" %in% names(all_enrichment_results)) {
  total_consistent <- sum(sapply(c("GO_BP", "KEGG", "Reactome"), function(db) {
    if (!is.null(all_enrichment_results$cross_comparison[[db]])) {
      nrow(all_enrichment_results$cross_comparison[[db]]$consistent)
    } else 0
  }))
  cat(sprintf("  GSEA cross-comparison (consistent): %d\n", total_consistent))
}

cat("\nMain biological themes:\n")

# Extract main themes from top pathways
themes_identified <- FALSE

# Check ORA results for themes
if ("R0_vs_R1_tumor_ORA" %in% names(all_enrichment_results)) {
  if (!is.null(all_enrichment_results$R0_vs_R1_tumor_ORA$GO_BP_down)) {
    top_go <- head(all_enrichment_results$R0_vs_R1_tumor_ORA$GO_BP_down@result, 10)
    if (nrow(top_go) > 0) {
      # Check for common themes
      if (any(grepl("cell cycle|mitosis|division", top_go$Description, ignore.case = TRUE))) {
        cat("  1. Cell cycle and proliferation (down-regulated)\n")
        themes_identified <- TRUE
      }
      if (any(grepl("immune|inflammatory", top_go$Description, ignore.case = TRUE))) {
        cat("  2. Immune response modulation\n")
        themes_identified <- TRUE
      }
      if (any(grepl("DNA|repair|damage", top_go$Description, ignore.case = TRUE))) {
        cat("  3. DNA damage response\n")
        themes_identified <- TRUE
      }
    }
  }
}

if (!themes_identified) {
  cat("  1. [TO BE DETERMINED from pathway inspection]\n")
  cat("  2. [TO BE DETERMINED from pathway inspection]\n")
  cat("  3. [TO BE DETERMINED from pathway inspection]\n")
}

# Cross-driver genes summary
if (exists("cross_driver_annotated") && !is.null(cross_driver_annotated)) {
  cat("\nCross-driver genes:\n")
  cat(sprintf("  Total: %d genes\n", nrow(cross_driver_annotated)))
  
  # Direction summary
  if (all(c("direction_RET", "direction_BRAF") %in% colnames(cross_driver_annotated))) {
    n_up_both <- sum(cross_driver_annotated$direction_RET == "UP" & 
                       cross_driver_annotated$direction_BRAF == "UP")
    n_down_both <- sum(cross_driver_annotated$direction_RET == "DOWN" & 
                         cross_driver_annotated$direction_BRAF == "DOWN")
    cat(sprintf("  Consistently up: %d genes\n", n_up_both))
    cat(sprintf("  Consistently down: %d genes\n", n_down_both))
  }
  
  # Top functional category
  if ("GENETYPE" %in% colnames(cross_driver_annotated)) {
    top_category <- names(sort(table(cross_driver_annotated$GENETYPE), 
                               decreasing = TRUE))[1]
    cat(sprintf("  Top category: %s\n", top_category))
  }
} else {
  cat("\nCross-driver genes: Analysis pending (run separate script)\n")
}

# Key insights from cross-comparison
if ("cross_comparison" %in% names(all_enrichment_results)) {
  cat("\nDriver-independent pathways:\n")
  
  # Count pathways enriched in same direction
  n_common_up <- 0
  n_common_down <- 0
  
  for (db in c("GO_BP", "KEGG", "Reactome")) {
    if (!is.null(all_enrichment_results$cross_comparison[[db]])) {
      consistent <- all_enrichment_results$cross_comparison[[db]]$consistent
      if (nrow(consistent) > 0) {
        n_common_up <- n_common_up + sum(consistent$NES_RET > 0)
        n_common_down <- n_common_down + sum(consistent$NES_RET < 0)
      }
    }
  }
  
  cat(sprintf("  Commonly up-regulated: %d pathways\n", n_common_up))
  cat(sprintf("  Commonly down-regulated: %d pathways\n", n_common_down))
  
  # Highlight if immune/inflammation pathways are common
  all_consistent <- do.call(rbind, lapply(c("GO_BP", "KEGG", "Reactome"), function(db) {
    if (!is.null(all_enrichment_results$cross_comparison[[db]])) {
      all_enrichment_results$cross_comparison[[db]]$consistent
    } else NULL
  }))
  
  if (!is.null(all_consistent) && nrow(all_consistent) > 0) {
    # Categorize pathways
    cat("\n  Pathway categories:\n")
    
    # Immune/inflammation
    immune_idx <- grep("immune|inflammatory|cytokine|chemokine|interferon|interleukin", 
                       all_consistent$Description, ignore.case = TRUE)
    cat(sprintf("  - Immune/inflammation: %d (%.1f%%)\n", 
                length(immune_idx), 100 * length(immune_idx) / nrow(all_consistent)))
    
    # Cell cycle/proliferation
    cycle_idx <- grep("cell cycle|mitosis|division|proliferation|chromosome", 
                      all_consistent$Description, ignore.case = TRUE)
    cat(sprintf("  - Cell cycle/proliferation: %d (%.1f%%)\n", 
                length(cycle_idx), 100 * length(cycle_idx) / nrow(all_consistent)))
    
    # DNA damage/repair
    dna_idx <- grep("DNA|repair|damage response|checkpoint", 
                    all_consistent$Description, ignore.case = TRUE)
    cat(sprintf("  - DNA damage/repair: %d (%.1f%%)\n", 
                length(dna_idx), 100 * length(dna_idx) / nrow(all_consistent)))
    
    # Cell death
    death_idx <- grep("apoptosis|death|killing|necrosis", 
                      all_consistent$Description, ignore.case = TRUE)
    cat(sprintf("  - Cell death: %d (%.1f%%)\n", 
                length(death_idx), 100 * length(death_idx) / nrow(all_consistent)))
    
    if (length(immune_idx) > 0) {
      cat("\n  * Immune/inflammatory pathways are driver-independent\n")
    }
  }
}

cat("\n=== END OF ENRICHMENT CHECK ===\n")

# ============================================================================
# Data Quality Check
# ============================================================================

cat("\n--- Data Quality Check ---\n")

# Check for unexpected patterns
if ("cross_comparison" %in% names(all_enrichment_results)) {
  # Check if all pathways are in same direction
  all_down <- TRUE
  all_up <- TRUE
  
  for (db in c("GO_BP", "KEGG", "Reactome")) {
    if (!is.null(all_enrichment_results$cross_comparison[[db]])) {
      consistent <- all_enrichment_results$cross_comparison[[db]]$consistent
      if (nrow(consistent) > 0) {
        if (any(consistent$NES_RET > 0)) all_down <- FALSE
        if (any(consistent$NES_RET < 0)) all_up <- FALSE
      }
    }
  }
  
  if (all_down) {
    cat("⚠ Note: All consistent pathways are down-regulated\n")
    cat("  This suggests strong suppressive effects of radiation exposure\n")
  } else if (all_up) {
    cat("⚠ Note: All consistent pathways are up-regulated\n")
    cat("  This is unusual and may warrant further investigation\n")
  }
  
  # Check NES consistency rate
  total_common <- sum(sapply(c("GO_BP", "KEGG", "Reactome"), function(db) {
    if (!is.null(all_enrichment_results$cross_comparison[[db]])) {
      nrow(all_enrichment_results$cross_comparison[[db]]$all)
    } else 0
  }))
  
  if (total_common > 0 && total_consistent == total_common) {
    cat("⚠ Note: 100% NES direction consistency between RET and BRAF\n")
    cat("  This indicates very strong driver-independent effects\n")
  }
}