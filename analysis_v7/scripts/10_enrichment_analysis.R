# 10_enrichment_analysis.R - Enrichment Analysis for DEG Results
# Purpose: Perform GO/KEGG/Reactome enrichment analysis on DEG results
# Method: ORA for R0_vs_R1_tumor, GSEA for all comparisons, Cross-comparison analysis
# Input: thyr_deg_results.rds (from 09_deg_analysis.R)
# Output: Enrichment results for each comparison and GSEA cross-comparison
# Version: v7.3 - Added GSEA cross-comparison analysis
# Date: 2025-10-16

source("analysis_v7/setup.R")

cat("\n=== Enrichment Analysis for DEG Results (v7.3) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Strategy:\n")
cat("  - R0_vs_R1_tumor: ORA (Up/Down/All) + GSEA for 3 databases\n")
cat("  - Other comparisons: GSEA only for 3 databases\n")
cat("  - Cross-comparison: Identify pathways consistent between RET and BRAF\n")

# Load packages
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(DOSE)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
})

# Record database versions
cat("\n--- Database Versions ---\n")
cat("  org.Hs.eg.db:", as.character(packageVersion("org.Hs.eg.db")), "\n")
cat("  ReactomePA:", as.character(packageVersion("ReactomePA")), "\n")
cat("  clusterProfiler:", as.character(packageVersion("clusterProfiler")), "\n")
cat("  KEGG data date: Current online version\n")

# Set seed for reproducibility
set.seed(1986)
cat("\n  Random seed set to: 1986\n")

# ============================================================================
# Configuration
# ============================================================================

CONFIG <- list(
  # Statistical threshold
  FDR_CUTOFF = 0.05,          # FDR threshold for Storey method
  
  # Gene set size limits (for both ORA and GSEA)
  MIN_GENESET_SIZE = 10,      # Minimum genes in a set
  MAX_GENESET_SIZE = 500,     # Maximum genes in a set
  
  # Output control
  VERBOSE = TRUE              # Verbose output
)

cat("\nConfiguration:\n")
cat("  FDR threshold:", CONFIG$FDR_CUTOFF, "\n")
cat("  Gene set size:", CONFIG$MIN_GENESET_SIZE, "-", CONFIG$MAX_GENESET_SIZE, "\n")
cat("  GSEA: Using fgseaMultilevel (adaptive)\n")

# ============================================================================
# Load DEG results
# ============================================================================

cat("\n--- Loading DEG results ---\n")

deg_results_path <- paste0(paths$processed, "thyr_deg_results.rds")
if (!file.exists(deg_results_path)) {
  stop("DEG results not found. Please run 09_deg_analysis.R first.")
}

deg_output <- readRDS(deg_results_path)
thyr_deg_results <- deg_output$comparison_results
consistency_results <- deg_output$consistency_results

cat("  Loaded DEG results for", length(thyr_deg_results), "comparisons\n")
cat("  Consistency results:", if(is.null(consistency_results)) "Not available" else "Loaded", "\n")

# ============================================================================
# Phase 1: Prepare Gene Mapping
# ============================================================================

cat("\n=== PHASE 1: GENE ID MAPPING ===\n")
cat("\n--- Preparing gene ID mappings ---\n")

# Extract all unique gene IDs from all comparisons
all_gene_ids <- unique(unlist(lapply(thyr_deg_results, function(comp) {
  comp$deg_summary$results_df$gene_id
})))

# Clean Ensembl IDs (remove version)
all_gene_ids_clean <- sub("\\..*", "", all_gene_ids)
cat("  Total unique genes:", length(all_gene_ids_clean), "\n")

# Map Ensembl to Entrez
gene_mapping <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = all_gene_ids_clean,
                                      columns = c("ENTREZID", "SYMBOL"),
                                      keytype = "ENSEMBL")

# Remove NA mappings
gene_mapping <- gene_mapping[!is.na(gene_mapping$ENTREZID), ]

# Check for 1:many mappings before removing duplicates
n_multi <- sum(duplicated(gene_mapping$ENSEMBL))
if (n_multi > 0) {
  cat(sprintf("  Note: %d Ensembl IDs map to multiple Entrez IDs (using first match)\n", n_multi))
}

gene_mapping <- gene_mapping[!duplicated(gene_mapping$ENSEMBL), ]

cat("  Successfully mapped:", nrow(gene_mapping), "genes\n")
cat("  Unmapped genes:", length(all_gene_ids_clean) - nrow(gene_mapping), "\n")

# Create output directory structure
output_dir <- paste0(paths$output, "enrichment_analysis/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\nPhase 1 setup complete!\n")

# ============================================================================
# Phase 2: ORA Analysis Functions
# ============================================================================

cat("\n=== PHASE 2: ORA FUNCTIONS ===\n")
cat("\n--- Setting up ORA analysis functions ---\n")

# Function to prepare gene lists for ORA
prepare_ora_genelist <- function(deg_results, direction = "all", gene_mapping) {
  # Get significant genes
  sig_genes <- deg_results$deg_summary$results_df[deg_results$deg_summary$results_df$significant, ]
  
  # Filter by direction if specified
  if (direction == "up") {
    sig_genes <- sig_genes[sig_genes$direction == "UP", ]
  } else if (direction == "down") {
    sig_genes <- sig_genes[sig_genes$direction == "DOWN", ]
  }
  
  # Clean Ensembl IDs
  sig_genes$ensembl_clean <- sub("\\..*", "", sig_genes$gene_id)
  
  # Map to Entrez IDs
  sig_genes_mapped <- merge(sig_genes, gene_mapping, 
                            by.x = "ensembl_clean", by.y = "ENSEMBL", 
                            all.x = FALSE)
  
  # Return unique Entrez IDs
  unique(sig_genes_mapped$ENTREZID)
}

# Function to perform ORA for GO
perform_go_ora <- function(gene_list, universe_list, ont = "BP") {
  enrichGO(
    gene = gene_list,
    universe = universe_list,
    OrgDb = org.Hs.eg.db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,      # Apply cutoff after Storey
    qvalueCutoff = 1,      # Apply cutoff after Storey
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE
  )
}

# Function to perform ORA for KEGG
perform_kegg_ora <- function(gene_list, universe_list) {
  enrichKEGG(
    gene = gene_list,
    universe = universe_list,
    organism = 'hsa',
    pAdjustMethod = "BH",
    pvalueCutoff = 1,      # Apply cutoff after Storey
    qvalueCutoff = 1,      # Apply cutoff after Storey
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE
  )
}

# Function to perform ORA for Reactome
perform_reactome_ora <- function(gene_list, universe_list) {
  enrichPathway(
    gene = gene_list,
    universe = universe_list,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,      # Apply cutoff after Storey
    qvalueCutoff = 1,      # Apply cutoff after Storey
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE
  )
}

# Function to apply Storey method for FDR control
apply_storey_to_enrichment <- function(enrich_result, alpha = CONFIG$FDR_CUTOFF) {
  if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
    return(enrich_result)
  }
  
  # Apply Storey's q-value
  library(qvalue)
  p_values <- enrich_result@result$pvalue
  
  # Handle edge cases
  if (length(p_values) < 2 || all(is.na(p_values))) {
    return(enrich_result)
  }
  
  # Calculate q-values
  q_obj <- tryCatch({
    qvalue(p_values, lambda = seq(0, 0.95, 0.05))
  }, error = function(e) {
    # Fallback to simple BH
    list(qvalues = p.adjust(p_values, method = "BH"))
  })
  
  # Update results
  enrich_result@result$qvalue <- q_obj$qvalues
  enrich_result@result$p.adjust <- q_obj$qvalues  # Update both columns
  
  # Filter by FDR threshold
  enrich_result@result <- enrich_result@result[enrich_result@result$qvalue < alpha, ]
  
  return(enrich_result)
}

# Function to save ORA results
save_ora_results <- function(result, comparison, database, direction, output_dir) {
  if (is.null(result) || nrow(result@result) == 0) {
    cat(sprintf("    No enriched terms for %s %s (%s)\n", comparison, database, direction))
    return(NULL)
  }
  
  # Create output subdirectory
  comp_dir <- paste0(output_dir, comparison, "/")
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save as CSV
  filename <- sprintf("%sORA_%s_%s.csv", comp_dir, database, direction)
  result_df <- as.data.frame(result)
  write.csv(result_df, filename, row.names = FALSE)
  
  cat(sprintf("    Saved %d enriched terms for %s (%s): %s\n", 
              nrow(result_df), database, direction, basename(filename)))
  
  return(result)
}

# Master function to run all ORA analyses for a comparison
run_ora_analysis <- function(comparison_name, thyr_deg_results, gene_mapping, output_dir) {
  cat(sprintf("\n  Running ORA for %s\n", comparison_name))
  
  # Get the comparison data
  deg_data <- thyr_deg_results[[comparison_name]]
  
  # Prepare universe (all tested genes)
  all_genes_clean <- sub("\\..*", "", deg_data$deg_summary$results_df$gene_id)
  universe_df <- merge(data.frame(ensembl_clean = all_genes_clean), 
                       gene_mapping,
                       by.x = "ensembl_clean", by.y = "ENSEMBL",
                       all.x = FALSE)
  universe_list <- unique(universe_df$ENTREZID)
  
  cat(sprintf("    Universe size: %d genes\n", length(universe_list)))
  
  # Store all results
  ora_results <- list()
  
  # Run ORA for each direction
  for (direction in c("up", "down", "all")) {
    # Prepare gene list
    gene_list <- prepare_ora_genelist(deg_data, direction, gene_mapping)
    cat(sprintf("    %s-regulated genes: %d\n", 
                if(direction == "all") "All" else toupper(direction),
                length(gene_list)))
    
    if (length(gene_list) < 5) {
      cat(sprintf("      Skipping %s (too few genes)\n", direction))
      next
    }
    
    # GO Biological Process
    cat("      Running GO BP...\n")
    ora_results[[paste0("GO_BP_", direction)]] <- tryCatch({
      result <- perform_go_ora(gene_list, universe_list, "BP")
      result <- apply_storey_to_enrichment(result)  # Apply Storey method
      save_ora_results(result, comparison_name, "GO_BP", direction, output_dir)
    }, error = function(e) {
      cat(sprintf("        Error in GO BP: %s\n", e$message))
      NULL
    })
    
    # KEGG
    cat("      Running KEGG...\n")
    ora_results[[paste0("KEGG_", direction)]] <- tryCatch({
      result <- perform_kegg_ora(gene_list, universe_list)
      result <- apply_storey_to_enrichment(result)  # Apply Storey method
      save_ora_results(result, comparison_name, "KEGG", direction, output_dir)
    }, error = function(e) {
      cat(sprintf("        Error in KEGG: %s\n", e$message))
      NULL
    })
    
    # Reactome
    cat("      Running Reactome...\n")
    ora_results[[paste0("Reactome_", direction)]] <- tryCatch({
      result <- perform_reactome_ora(gene_list, universe_list)
      result <- apply_storey_to_enrichment(result)  # Apply Storey method
      save_ora_results(result, comparison_name, "Reactome", direction, output_dir)
    }, error = function(e) {
      cat(sprintf("        Error in Reactome: %s\n", e$message))
      NULL
    })
  }
  
  return(ora_results)
}

cat("ORA analysis functions ready!\n")
cat("\nPhase 2 complete! Ready for GSEA setup.\n")

# ============================================================================
# Phase 3: GSEA Analysis Functions
# ============================================================================

cat("\n=== PHASE 3: GSEA FUNCTIONS ===\n")
cat("\n--- Setting up GSEA analysis functions ---\n")

# Function to prepare ranked gene list for GSEA
prepare_gsea_genelist <- function(deg_results, gene_mapping) {
  # Get all genes with statistics
  all_genes <- deg_results$deg_summary$results_df
  
  # Clean Ensembl IDs
  all_genes$ensembl_clean <- sub("\\..*", "", all_genes$gene_id)
  
  # Map to Entrez IDs
  all_genes_mapped <- merge(all_genes, gene_mapping,
                            by.x = "ensembl_clean", by.y = "ENSEMBL",
                            all.x = FALSE)
  
  # Create ranked list (by fold change)
  gene_list <- all_genes_mapped$FC
  names(gene_list) <- all_genes_mapped$ENTREZID
  
  # Sort by fold change (descending)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Remove duplicates if any
  gene_list <- gene_list[!duplicated(names(gene_list))]
  
  return(gene_list)
}

# Function to perform GSEA for GO
perform_go_gsea <- function(gene_list, ont = "BP") {
  gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,      # Apply cutoff after Storey
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE,
    eps = 0,
    seed = 1986
  )
}

# Function to perform GSEA for KEGG
perform_kegg_gsea <- function(gene_list) {
  gseKEGG(
    geneList = gene_list,
    organism = 'hsa',
    pAdjustMethod = "BH",
    pvalueCutoff = 1,      # Apply cutoff after Storey
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE,
    eps = 0,
    seed = 1986
  )
}

# Function to perform GSEA for Reactome
perform_reactome_gsea <- function(gene_list) {
  gsePathway(
    geneList = gene_list,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,      # Apply cutoff after Storey
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE,
    eps = 0,
    seed = 1986
  )
}

# Function to save GSEA results
save_gsea_results <- function(result, comparison, database, output_dir) {
  if (is.null(result) || nrow(result@result) == 0) {
    cat(sprintf("      No enriched terms for %s %s\n", comparison, database))
    return(NULL)
  }
  
  # Create output subdirectory
  comp_dir <- paste0(output_dir, comparison, "/")
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save as CSV
  filename <- sprintf("%sGSEA_%s.csv", comp_dir, database)
  result_df <- as.data.frame(result)
  write.csv(result_df, filename, row.names = FALSE)
  
  cat(sprintf("      Saved %d enriched terms for %s: %s\n", 
              nrow(result_df), database, basename(filename)))
  
  return(result)
}

# Master function to run all GSEA analyses for a comparison
run_gsea_analysis <- function(comparison_name, thyr_deg_results, gene_mapping, output_dir) {
  cat(sprintf("\n  Running GSEA for %s\n", comparison_name))
  
  # Get the comparison data
  deg_data <- thyr_deg_results[[comparison_name]]
  
  # Prepare ranked gene list
  gene_list <- prepare_gsea_genelist(deg_data, gene_mapping)
  cat(sprintf("    Ranked gene list: %d genes\n", length(gene_list)))
  cat(sprintf("    Top gene: %s (FC=%.2f)\n", 
              names(gene_list)[1], gene_list[1]))
  cat(sprintf("    Bottom gene: %s (FC=%.2f)\n", 
              names(gene_list)[length(gene_list)], gene_list[length(gene_list)]))
  
  # Store all results
  gsea_results <- list()
  
  # GO Biological Process
  cat("    Running GO BP GSEA...\n")
  gsea_results[["GO_BP"]] <- tryCatch({
    result <- perform_go_gsea(gene_list, "BP")
    result <- apply_storey_to_enrichment(result)  # Apply Storey method
    save_gsea_results(result, comparison_name, "GO_BP", output_dir)
  }, error = function(e) {
    cat(sprintf("      Error in GO BP GSEA: %s\n", e$message))
    NULL
  })
  
  # KEGG
  cat("    Running KEGG GSEA...\n")
  gsea_results[["KEGG"]] <- tryCatch({
    result <- perform_kegg_gsea(gene_list)
    result <- apply_storey_to_enrichment(result)  # Apply Storey method
    save_gsea_results(result, comparison_name, "KEGG", output_dir)
  }, error = function(e) {
    cat(sprintf("      Error in KEGG GSEA: %s\n", e$message))
    NULL
  })
  
  # Reactome
  cat("    Running Reactome GSEA...\n")
  gsea_results[["Reactome"]] <- tryCatch({
    result <- perform_reactome_gsea(gene_list)
    result <- apply_storey_to_enrichment(result)  # Apply Storey method
    save_gsea_results(result, comparison_name, "Reactome", output_dir)
  }, error = function(e) {
    cat(sprintf("      Error in Reactome GSEA: %s\n", e$message))
    NULL
  })
  
  return(gsea_results)
}

cat("GSEA analysis functions ready!\n")
cat("\nPhase 3 complete! Ready to run enrichment analyses.\n")

# ============================================================================
# Phase 4: Execute Enrichment Analyses
# ============================================================================

cat("\n=== PHASE 4: EXECUTE ENRICHMENT ANALYSES ===\n")

# Check if results already exist
results_file <- paste0(output_dir, "all_enrichment_results.rds")
if (file.exists(results_file)) {
  cat("Previous results detected. Loading from file...\n")
  all_enrichment_results <- readRDS(results_file)
  cat("Loaded enrichment results successfully.\n")
  
  # Count results
  if ("R0_vs_R1_tumor_ORA" %in% names(all_enrichment_results)) {
    ora_count <- sum(!sapply(all_enrichment_results[["R0_vs_R1_tumor_ORA"]], is.null))
    cat(sprintf("  ORA analyses: %d\n", ora_count))
  }
  
  if ("GSEA" %in% names(all_enrichment_results)) {
    gsea_count <- 0
    for (comp in names(all_enrichment_results[["GSEA"]])) {
      gsea_count <- gsea_count + sum(!sapply(all_enrichment_results[["GSEA"]][[comp]], is.null))
    }
    cat(sprintf("  GSEA analyses: %d\n", gsea_count))
  }
  
  cat("\nSkipping re-analysis. To force re-analysis, delete:\n")
  cat("  ", results_file, "\n")
  
} else {
  # Execute new analysis
  cat("No previous results found. Starting new analysis...\n")
  
  # Initialize results storage
  all_enrichment_results <- list()
  
  # ----------------------------------------------------------------------------
  # 1. ORA for R0_vs_R1_tumor
  # ----------------------------------------------------------------------------
  
  if ("R0_vs_R1_tumor" %in% names(thyr_deg_results)) {
    cat("\n[ORA Analysis] R0_vs_R1_tumor\n")
    cat("======================================\n")
    
    ora_results <- run_ora_analysis(
      comparison_name = "R0_vs_R1_tumor",
      thyr_deg_results = thyr_deg_results,
      gene_mapping = gene_mapping,
      output_dir = output_dir
    )
    
    all_enrichment_results[["R0_vs_R1_tumor_ORA"]] <- ora_results
    
    # Summary
    successful_ora <- sum(!sapply(ora_results, is.null))
    cat(sprintf("\n  ORA Complete: %d/%d analyses successful\n", 
                successful_ora, length(ora_results)))
  }
  
  # ----------------------------------------------------------------------------
  # 2. GSEA for all comparisons
  # ----------------------------------------------------------------------------
  
  cat("\n[GSEA Analysis] All Comparisons\n")
  cat("======================================\n")
  
  comparisons_to_analyze <- names(thyr_deg_results)
  gsea_results_all <- list()
  
  for (comp_name in comparisons_to_analyze) {
    cat(sprintf("\nProcessing: %s\n", comp_name))
    cat(paste(rep("-", 40), collapse = ""), "\n")
    
    gsea_results <- run_gsea_analysis(
      comparison_name = comp_name,
      thyr_deg_results = thyr_deg_results,
      gene_mapping = gene_mapping,
      output_dir = output_dir
    )
    
    gsea_results_all[[comp_name]] <- gsea_results
    
    # Summary for this comparison
    successful_gsea <- sum(!sapply(gsea_results, is.null))
    cat(sprintf("  GSEA Complete: %d/3 databases enriched\n", successful_gsea))
  }
  
  all_enrichment_results[["GSEA"]] <- gsea_results_all
  
  # ----------------------------------------------------------------------------
  # 3. Overall Summary
  # ----------------------------------------------------------------------------
  
  cat("\n=== ENRICHMENT ANALYSIS SUMMARY ===\n")
  cat("=====================================\n")
  
  # Count ORA results
  if ("R0_vs_R1_tumor_ORA" %in% names(all_enrichment_results)) {
    ora_count <- sum(!sapply(all_enrichment_results[["R0_vs_R1_tumor_ORA"]], is.null))
    cat(sprintf("ORA analyses completed: %d\n", ora_count))
  }
  
  # Count GSEA results
  gsea_count <- 0
  for (comp in names(all_enrichment_results[["GSEA"]])) {
    gsea_count <- gsea_count + sum(!sapply(all_enrichment_results[["GSEA"]][[comp]], is.null))
  }
  cat(sprintf("GSEA analyses completed: %d\n", gsea_count))
  
  # Total analyses
  total_analyses <- ifelse(exists("ora_count"), ora_count, 0) + gsea_count
  cat(sprintf("\nTotal enrichment analyses: %d\n", total_analyses))
  
  # Save consolidated results object
  saveRDS(all_enrichment_results, results_file)
  cat("\nResults saved: all_enrichment_results.rds\n")
  
  # List output directory contents
  output_files <- list.files(output_dir, recursive = TRUE, pattern = "\\.csv$")
  cat(sprintf("\nGenerated %d result files in %s\n", 
              length(output_files), output_dir))
}

cat("\nPhase 4 complete! Main enrichment analyses finished.\n")
cat("Next: GSEA cross-comparison analysis (Phase 5)\n")

# ============================================================================
# Phase 5: GSEA Cross-Comparison Analysis
# ============================================================================

cat("\n=== PHASE 5: GSEA CROSS-COMPARISON ANALYSIS ===\n")
cat("Finding pathways consistently enriched across RET and BRAF comparisons\n\n")

# Check if GSEA results exist for both comparisons
if ("GSEA" %in% names(all_enrichment_results) &&
    "R0_vs_R1_tumor" %in% names(all_enrichment_results[["GSEA"]]) &&
    "B0_vs_B1_tumor" %in% names(all_enrichment_results[["GSEA"]])) {
  
  # Create output directory for cross-comparison
  cross_comp_dir <- paste0(output_dir, "cross_comparison/")
  dir.create(cross_comp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize storage for all cross-comparison results
  cross_comparison_results <- list()
  
  # Function to find consistent pathways between two comparisons
  find_consistent_pathways <- function(result1, result2, db_name) {
    # Extract significant pathways from both comparisons
    if (is.null(result1) || is.null(result2)) {
      return(NULL)
    }
    
    pathways1 <- result1@result[result1@result$qvalue < CONFIG$FDR_CUTOFF, ]
    pathways2 <- result2@result[result2@result$qvalue < CONFIG$FDR_CUTOFF, ]
    
    if (nrow(pathways1) == 0 || nrow(pathways2) == 0) {
      return(NULL)
    }
    
    # Merge by pathway description
    common_pathways <- merge(
      pathways1[, c("ID", "Description", "NES", "qvalue", "setSize")],
      pathways2[, c("ID", "Description", "NES", "qvalue", "setSize")],
      by = "Description",
      suffixes = c("_RET", "_BRAF")
    )
    
    if (nrow(common_pathways) == 0) {
      return(NULL)
    }
    
    # Check for NES direction consistency
    common_pathways$NES_consistent <- sign(common_pathways$NES_RET) == sign(common_pathways$NES_BRAF)
    
    # Calculate consistency metrics
    common_pathways$avg_abs_NES <- (abs(common_pathways$NES_RET) + abs(common_pathways$NES_BRAF)) / 2
    common_pathways$min_qvalue <- pmin(common_pathways$qvalue_RET, common_pathways$qvalue_BRAF)
    common_pathways$consistency_score <- common_pathways$avg_abs_NES * ifelse(common_pathways$NES_consistent, 1, -0.5)
    
    # Sort by consistency score
    common_pathways <- common_pathways[order(common_pathways$consistency_score, decreasing = TRUE), ]
    
    # Add database name
    common_pathways$database <- db_name
    
    return(common_pathways)
  }
  
  # Analyze each database
  databases <- c("GO_BP", "KEGG", "Reactome")
  
  for (db in databases) {
    cat(sprintf("\n--- %s Cross-comparison ---\n", db))
    
    # Get results from both comparisons
    ret_result <- all_enrichment_results[["GSEA"]][["R0_vs_R1_tumor"]][[db]]
    braf_result <- all_enrichment_results[["GSEA"]][["B0_vs_B1_tumor"]][[db]]
    
    # Find consistent pathways
    consistent <- find_consistent_pathways(ret_result, braf_result, db)
    
    if (!is.null(consistent) && nrow(consistent) > 0) {
      # Filter for NES-consistent pathways
      consistent_only <- consistent[consistent$NES_consistent, ]
      inconsistent <- consistent[!consistent$NES_consistent, ]
      
      cat(sprintf("  Total common pathways: %d\n", nrow(consistent)))
      cat(sprintf("  NES-consistent pathways: %d\n", nrow(consistent_only)))
      cat(sprintf("  NES-inconsistent pathways: %d\n", nrow(inconsistent)))
      
      # Show top consistent pathways
      if (nrow(consistent_only) > 0) {
        cat("\n  Top 5 consistent pathways:\n")
        top5 <- head(consistent_only, 5)
        for (i in 1:nrow(top5)) {
          cat(sprintf("    %d. %s\n", i, top5$Description[i]))
          cat(sprintf("       NES_RET=%.2f, NES_BRAF=%.2f, Score=%.2f\n",
                      top5$NES_RET[i], top5$NES_BRAF[i], top5$consistency_score[i]))
        }
        
        # Save results
        filename <- paste0(cross_comp_dir, "cross_comparison_", db, "_consistent.csv")
        write.csv(consistent_only, filename, row.names = FALSE)
        cat(sprintf("\n  Saved %d consistent pathways to: %s\n", 
                    nrow(consistent_only), basename(filename)))
      }
      
      # Store in results
      cross_comparison_results[[db]] <- list(
        all = consistent,
        consistent = consistent_only,
        inconsistent = inconsistent
      )
      
    } else {
      cat("  No common significant pathways found\n")
      cross_comparison_results[[db]] <- NULL
    }
  }
  
  # ============================================================================
  # Create consolidated Excel report
  # ============================================================================
  
  cat("\n--- Creating consolidated cross-comparison report ---\n")
  
  wb <- createWorkbook()
  
  # Add overview sheet
  addWorksheet(wb, "Overview")
  overview_data <- data.frame(
    Database = databases,
    Common_Pathways = sapply(databases, function(db) {
      if (!is.null(cross_comparison_results[[db]])) nrow(cross_comparison_results[[db]]$all) else 0
    }),
    Consistent_NES = sapply(databases, function(db) {
      if (!is.null(cross_comparison_results[[db]])) nrow(cross_comparison_results[[db]]$consistent) else 0
    }),
    Inconsistent_NES = sapply(databases, function(db) {
      if (!is.null(cross_comparison_results[[db]])) nrow(cross_comparison_results[[db]]$inconsistent) else 0
    })
  )
  writeData(wb, "Overview", overview_data)
  
  # Add header style
  header_style <- createStyle(
    fontColour = "white", fgFill = "#4F81BD",
    halign = "center", textDecoration = "bold"
  )
  addStyle(wb, "Overview", header_style, rows = 1, cols = 1:4, gridExpand = TRUE)
  
  # Add sheets for each database with consistent pathways
  for (db in databases) {
    if (!is.null(cross_comparison_results[[db]]) && 
        nrow(cross_comparison_results[[db]]$consistent) > 0) {
      
      sheet_name <- paste0(db, "_Consistent")
      addWorksheet(wb, sheet_name)
      
      # Select columns to display
      display_data <- cross_comparison_results[[db]]$consistent[, c(
        "Description", "NES_RET", "qvalue_RET", "NES_BRAF", "qvalue_BRAF",
        "avg_abs_NES", "consistency_score"
      )]
      
      # Round numeric columns
      numeric_cols <- c("NES_RET", "qvalue_RET", "NES_BRAF", "qvalue_BRAF", 
                        "avg_abs_NES", "consistency_score")
      display_data[numeric_cols] <- lapply(display_data[numeric_cols], function(x) round(x, 4))
      
      writeData(wb, sheet_name, display_data)
      addStyle(wb, sheet_name, header_style, rows = 1, 
               cols = 1:ncol(display_data), gridExpand = TRUE)
      
      # Auto-size columns
      setColWidths(wb, sheet_name, cols = 1:ncol(display_data), widths = "auto")
      
      # Add conditional formatting for high consistency scores
      if (nrow(display_data) > 1) {
        high_consistency <- which(display_data$consistency_score > 2)
        if (length(high_consistency) > 0) {
          # Apply style row by row with error handling
          tryCatch({
            for (row_idx in high_consistency) {
              addStyle(wb, sheet_name,
                       createStyle(fgFill = "#E8F5E8"),
                       rows = row_idx + 1,
                       cols = 1:ncol(display_data),
                       gridExpand = TRUE)
            }
          }, error = function(e) {
            cat("  Note: Could not apply conditional formatting (non-critical).\n")
          })
        }
      }
    }
  }
  
  # Save Excel file
  excel_filename <- paste0(cross_comp_dir, "GSEA_cross_comparison_results.xlsx")
  tryCatch({
    saveWorkbook(wb, excel_filename, overwrite = TRUE)
    cat("  Saved Excel report to:", basename(excel_filename), "\n")
  }, error = function(e) {
    cat("  Warning: Could not save Excel file:", e$message, "\n")
    cat("  CSV files are still available in the output directory.\n")
  })
  
  # Add to all_enrichment_results
  all_enrichment_results[["cross_comparison"]] <- cross_comparison_results
  
  # Re-save consolidated results
  saveRDS(all_enrichment_results, results_file)
  cat("\nUpdated results saved with cross-comparison analysis\n")
  
  # ============================================================================
  # Summary statistics
  # ============================================================================
  
  cat("\n=== CROSS-COMPARISON SUMMARY ===\n")
  
  total_consistent <- sum(sapply(databases, function(db) {
    if (!is.null(cross_comparison_results[[db]])) {
      nrow(cross_comparison_results[[db]]$consistent)
    } else 0
  }))
  
  cat(sprintf("Total consistent pathways across all databases: %d\n", total_consistent))
  
  # Find the most consistent pathways across databases
  all_consistent <- do.call(rbind, lapply(databases, function(db) {
    if (!is.null(cross_comparison_results[[db]]) && 
        nrow(cross_comparison_results[[db]]$consistent) > 0) {
      cross_comparison_results[[db]]$consistent
    } else NULL
  }))
  
  if (!is.null(all_consistent) && nrow(all_consistent) > 0) {
    cat("\nTop 10 most consistent pathways (all databases):\n")
    all_consistent <- all_consistent[order(all_consistent$consistency_score, decreasing = TRUE), ]
    top10 <- head(all_consistent, 10)
    for (i in 1:nrow(top10)) {
      cat(sprintf("  %2d. [%s] %s (Score: %.2f)\n",
                  i, top10$database[i], 
                  substr(top10$Description[i], 1, 60),
                  top10$consistency_score[i]))
    }
  }
  
} else {
  cat("GSEA results not available for cross-comparison\n")
  cat("Please ensure GSEA analysis has been completed for both R0_vs_R1_tumor and B0_vs_B1_tumor\n")
}

# ============================================================================
# Cleanup and Final Summary
# ============================================================================

cat("\n=== ENRICHMENT ANALYSIS COMPLETE ===\n")
cat("Results available in:\n")
cat("  RDS: ", results_file, "\n")
cat("  CSV files: ", output_dir, "\n")

# Keep only essential objects
rm(list = setdiff(ls(), c("paths", "all_enrichment_results", "thyr_deg_results")))
gc()

cat("\nNext steps:\n")
cat("  1. Review enrichment results in output/enrichment_analysis/\n")
cat("  2. Review cross-comparison results in output/enrichment_analysis/cross_comparison/\n")
cat("  3. Proceed to cross-driver gene annotation if needed (separate script)\n")
cat("  4. Proceed to REO panel construction (20-26_reo_*.R)\n")

cat("\n=== Script Complete ===\n")