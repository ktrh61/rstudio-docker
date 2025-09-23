# 10_enrichment_analysis.R - Enrichment Analysis for DEG Results
# Purpose: Perform GO/KEGG/Reactome enrichment analysis on DEG results
# Method: ORA for R0_vs_R1_tumor, GSEA for all comparisons
# Input: thyr_deg_results.rds (from 09_deg_analysis.R)
# Output: Enrichment results for each comparison
# Version: v7.0 - Comprehensive enrichment with cross-driver annotation
# Date: 2025-01-22

source("analysis_v7/setup.R")

cat("\n=== Enrichment Analysis for DEG Results (v7.0) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Strategy:\n")
cat("  - R0_vs_R1_tumor: ORA (Up/Down/All) + GSEA for 3 databases\n")
cat("  - Other comparisons: GSEA only for 3 databases\n")
cat("  - Cross-driver: Individual gene annotation\n")

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
thyr_deg_results <- deg_output$deg_results

# Summary of available comparisons
cat("Available comparisons:\n")
for (comp_name in names(thyr_deg_results)) {
  n_degs <- thyr_deg_results[[comp_name]]$deg_summary$summary_stats$significant_genes
  n_total <- thyr_deg_results[[comp_name]]$deg_summary$summary_stats$total_genes_tested
  cat(sprintf("  %s: %d DEGs from %d genes (%.1f%%)\n", 
              comp_name, n_degs, n_total, n_degs/n_total*100))
}

# ============================================================================
# Prepare gene ID conversion
# ============================================================================

cat("\n--- Preparing gene ID conversion ---\n")

# Extract all unique gene IDs from results
all_gene_ids <- unique(unlist(lapply(thyr_deg_results, function(x) {
  x$deg_summary$results_df$gene_id
})))

# Remove Ensembl version numbers (e.g., ENSG00000001.5 -> ENSG00000001)
# This is standard practice as version numbers change between releases
all_gene_ids_clean <- sub("\\..*", "", all_gene_ids)

cat("Total unique genes:", length(all_gene_ids_clean), "\n")

# Create mapping from Ensembl to Entrez
cat("Creating Ensembl to Entrez ID mapping...\n")
gene_mapping <- tryCatch({
  suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db,
                          keys = all_gene_ids_clean,
                          columns = c("ENTREZID", "SYMBOL", "GENENAME"),
                          keytype = "ENSEMBL")
  )
}, error = function(e) {
  cat("  Warning: Some genes could not be mapped\n")
  suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db,
                          keys = all_gene_ids_clean,
                          columns = c("ENTREZID", "SYMBOL", "GENENAME"),
                          keytype = "ENSEMBL")
  )
})

# Remove duplicates and NAs
# For 1:many mappings (one Ensembl -> multiple Entrez), keep the first
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
  enrichGO(gene = gene_list,
           universe = universe_list,
           OrgDb = org.Hs.eg.db,
           ont = ont,
           pAdjustMethod = "BH",
           pvalueCutoff = 1,  # Get all results
           qvalueCutoff = 1,  # Get all results
           minGSSize = CONFIG$MIN_GENESET_SIZE,
           maxGSSize = CONFIG$MAX_GENESET_SIZE,
           readable = TRUE)
}

# Function to perform ORA for KEGG
perform_kegg_ora <- function(gene_list, universe_list) {
  enrichKEGG(gene = gene_list,
             universe = universe_list,
             organism = 'hsa',
             pAdjustMethod = "BH",
             pvalueCutoff = 1,  # Get all results
             qvalueCutoff = 1,  # Get all results
             minGSSize = CONFIG$MIN_GENESET_SIZE,
             maxGSSize = CONFIG$MAX_GENESET_SIZE)
}

# Function to perform ORA for Reactome
perform_reactome_ora <- function(gene_list, universe_list) {
  enrichPathway(gene = gene_list,
                universe = universe_list,
                organism = "human",
                pAdjustMethod = "BH",
                pvalueCutoff = 1,  # Get all results
                qvalueCutoff = 1,  # Get all results
                minGSSize = CONFIG$MIN_GENESET_SIZE,
                maxGSSize = CONFIG$MAX_GENESET_SIZE,
                readable = TRUE)
}

# Function to apply Storey method to enrichment results
apply_storey_to_enrichment <- function(enrich_result, alpha = CONFIG$FDR_CUTOFF) {
  if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
    return(enrich_result)
  }
  
  # Get p-values
  pvalues <- enrich_result@result$pvalue
  
  # Apply Storey method with bootstrap (same as DEG analysis)
  library(qvalue)
  q_obj <- qvalue(pvalues, pi0.method = "bootstrap")
  
  # Update qvalues
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
cat("\nPhase 2 complete! Ready to perform enrichment analyses.\n")

# ============================================================================
# Phase 3: GSEA Analysis Functions
# ============================================================================

cat("\n--- Setting up GSEA analysis functions ---\n")

# Function to prepare ranked gene list for GSEA
prepare_gsea_genelist <- function(deg_results, gene_mapping) {
  # Get all genes with their fold changes
  all_genes <- deg_results$deg_summary$results_df
  
  # Clean Ensembl IDs
  all_genes$ensembl_clean <- sub("\\..*", "", all_genes$gene_id)
  
  # Map to Entrez IDs
  all_genes_mapped <- merge(all_genes, gene_mapping,
                            by.x = "ensembl_clean", by.y = "ENSEMBL",
                            all.x = FALSE)
  
  # Remove genes with NA log2FC
  all_genes_mapped <- all_genes_mapped[!is.na(all_genes_mapped$log2FC), ]
  
  # Create named vector of log2FC values
  gene_list <- all_genes_mapped$log2FC
  names(gene_list) <- all_genes_mapped$ENTREZID
  
  # Handle duplicate Entrez IDs (different Ensembl IDs mapping to same Entrez)
  # Keep the one with largest absolute fold change
  if (any(duplicated(names(gene_list)))) {
    n_dup <- sum(duplicated(names(gene_list)))
    cat(sprintf("    Note: %d duplicate Entrez IDs found (using max |log2FC|)\n", n_dup))
    
    gene_list_dedup <- tapply(gene_list, names(gene_list), function(x) x[which.max(abs(x))])
    # Convert array to numeric vector
    gene_list <- as.numeric(gene_list_dedup)
    names(gene_list) <- names(gene_list_dedup)
  }
  
  # Sort by decreasing fold change
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}

# Function to perform GSEA for GO
perform_go_gsea <- function(gene_list, ont = "BP") {
  gseGO(geneList = gene_list,
        OrgDb = org.Hs.eg.db,
        ont = ont,
        minGSSize = CONFIG$MIN_GENESET_SIZE,
        maxGSSize = CONFIG$MAX_GENESET_SIZE,
        pvalueCutoff = 1,  # Get all results
        pAdjustMethod = "BH",
        verbose = FALSE,
        seed = 123)
}

# Function to perform GSEA for KEGG
perform_kegg_gsea <- function(gene_list) {
  gseKEGG(geneList = gene_list,
          organism = 'hsa',
          minGSSize = CONFIG$MIN_GENESET_SIZE,
          maxGSSize = CONFIG$MAX_GENESET_SIZE,
          pvalueCutoff = 1,  # Get all results
          pAdjustMethod = "BH",
          verbose = FALSE,
          seed = 123)
}

# Function to perform GSEA for Reactome
perform_reactome_gsea <- function(gene_list) {
  gsePathway(geneList = gene_list,
             organism = "human",
             minGSSize = CONFIG$MIN_GENESET_SIZE,
             maxGSSize = CONFIG$MAX_GENESET_SIZE,
             pvalueCutoff = 1,  # Get all results
             pAdjustMethod = "BH",
             verbose = FALSE,
             seed = 123)
}

# Function to save GSEA results
save_gsea_results <- function(result, comparison, database, output_dir) {
  if (is.null(result) || nrow(result@result) == 0) {
    cat(sprintf("    No enriched gene sets for %s %s\n", comparison, database))
    return(NULL)
  }
  
  # Create output subdirectory
  comp_dir <- paste0(output_dir, comparison, "/")
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save as CSV
  filename <- sprintf("%sGSEA_%s.csv", comp_dir, database)
  result_df <- as.data.frame(result)
  write.csv(result_df, filename, row.names = FALSE)
  
  cat(sprintf("    Saved %d enriched gene sets for %s: %s\n", 
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

cat("\n=== EXECUTING ENRICHMENT ANALYSES ===\n")

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
cat("Next: Cross-driver gene annotation (Phase 5)\n")