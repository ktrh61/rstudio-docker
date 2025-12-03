# 10_enrichment_analysis.R - Enrichment Analysis for DEG Results
# Purpose: Perform GO/KEGG/Reactome/Hallmark enrichment analysis on DEG results
# Method: GSEA for all comparisons (primary), ORA for R0_vs_R1_tumor (supplementary)
# Input: thyr_deg_results.rds (from 09_deg_analysis.R)
# Output: Enrichment results for each comparison and cross-comparison analysis
# Version: v7.6 - Refined GSEA strategy with selective comparisons
# Date: 2025-12-03

source("analysis_v7/setup.R")

cat("\n=== Enrichment Analysis for DEG Results (v7.6) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Strategy:\n")
cat("  - GSEA for 3 comparisons (excluding B0_vs_B1_normal)\n")
cat("  - Primary result: R0_vs_R1_tumor\n")
cat("  - Cross-comparison 1: RET × BRAF (tumor vs tumor)\n")
cat("  - Cross-comparison 2: Tumor × Normal (within RET)\n")
cat("  - Supplementary: ORA for R0_vs_R1_tumor\n")
cat("\nChanges from v7.5:\n")
cat("  - GSEA ranking: delta × min(-log10(q), 10) (restored)\n")
cat("  - Excluded: B0_vs_B1_normal (no biological signal)\n")
cat("  - Added: Tumor × Normal cross-comparison\n")

# ============================================================================
# Load packages
# ============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(msigdbr)
  library(enrichplot)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
})

# Record database versions
cat("\n--- Database Versions ---\n")
cat("  org.Hs.eg.db:", as.character(packageVersion("org.Hs.eg.db")), "\n")
cat("  ReactomePA:", as.character(packageVersion("ReactomePA")), "\n")
cat("  clusterProfiler:", as.character(packageVersion("clusterProfiler")), "\n")
cat("  msigdbr:", as.character(packageVersion("msigdbr")), "\n")
cat("  KEGG data: Current online version\n")

# Set seed for reproducibility
set.seed(1986)
cat("\n  Random seed set to: 1986\n")

# ============================================================================
# Configuration
# ============================================================================

CONFIG <- list(
  # Statistical threshold
  FDR_CUTOFF = 0.05,          # FDR threshold (BH method)
  
  # GSEA ranking cap
  Q_CAP = 10,                 # Cap for -log10(q) in ranking
  
  # Gene set size limits
  MIN_GENESET_SIZE = 10,      # Minimum genes in a set
  MAX_GENESET_SIZE = 500,     # Maximum genes in a set
  
  # Comparisons to analyze (exclude B0_vs_B1_normal)
  GSEA_COMPARISONS = c("R0_vs_R1_tumor", "R0_vs_R1_normal", "B0_vs_B1_tumor"),
  
  # Output control
  VERBOSE = TRUE              # Verbose output
)

cat("\nConfiguration:\n")
cat("  FDR correction: BH method\n")
cat("  FDR threshold:", CONFIG$FDR_CUTOFF, "\n")
cat("  GSEA ranking: delta × min(-log10(q),", CONFIG$Q_CAP, ")\n")
cat("  GSEA comparisons:", paste(CONFIG$GSEA_COMPARISONS, collapse = ", "), "\n")
cat("  Gene set size:", CONFIG$MIN_GENESET_SIZE, "-", CONFIG$MAX_GENESET_SIZE, "\n")

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
consistency_results <- deg_output$consistency_results

cat("  Loaded DEG results for", length(thyr_deg_results), "comparisons\n")
for (comp in names(thyr_deg_results)) {
  n_degs <- thyr_deg_results[[comp]]$deg_summary$summary_stats$significant_genes
  cat(sprintf("    %s: %d DEGs\n", comp, n_degs))
}

# ============================================================================
# Phase 1: Gene ID Mapping
# ============================================================================

cat("\n=== PHASE 1: GENE ID MAPPING ===\n")

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

# Check for 1:many mappings
n_multi <- sum(duplicated(gene_mapping$ENSEMBL))
if (n_multi > 0) {
  cat(sprintf("  Note: %d Ensembl IDs map to multiple Entrez IDs (using first match)\n", n_multi))
}

gene_mapping <- gene_mapping[!duplicated(gene_mapping$ENSEMBL), ]

cat("  Successfully mapped:", nrow(gene_mapping), "genes\n")
cat("  Unmapped genes:", length(all_gene_ids_clean) - nrow(gene_mapping), "\n")

# Create output directory
output_dir <- paste0(paths$output, "enrichment_analysis/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\nPhase 1 complete!\n")

# ============================================================================
# Phase 2: Prepare MSigDB Hallmark Gene Sets
# ============================================================================

cat("\n=== PHASE 2: PREPARE HALLMARK GENE SETS ===\n")

# Get Hallmark gene sets
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
cat("  Loaded MSigDB Hallmark:", length(unique(hallmark_df$gs_name)), "gene sets\n")

# Convert to list format for GSEA
hallmark_list <- split(hallmark_df$ncbi_gene, hallmark_df$gs_name)

# Clean names (remove "HALLMARK_" prefix for readability)
names(hallmark_list) <- gsub("^HALLMARK_", "", names(hallmark_list))

cat("  Gene sets prepared for GSEA\n")
cat("  Example sets:", paste(head(names(hallmark_list), 5), collapse = ", "), "...\n")

# Create TERM2GENE format for clusterProfiler
hallmark_t2g <- hallmark_df %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  dplyr::mutate(gs_name = gsub("^HALLMARK_", "", gs_name))

cat("\nPhase 2 complete!\n")

# ============================================================================
# Phase 3: GSEA Functions
# ============================================================================

cat("\n=== PHASE 3: GSEA FUNCTIONS ===\n")

# Function to prepare ranked gene list for GSEA
# Uses Cliff's delta × min(-log10(q), Q_CAP)
prepare_gsea_genelist <- function(deg_results, gene_mapping, q_cap = CONFIG$Q_CAP) {
  # Get all genes with statistics
  all_genes <- deg_results$deg_summary$results_df
  
  # Clean Ensembl IDs
  all_genes$ensembl_clean <- sub("\\..*", "", all_genes$gene_id)
  
  # Map to Entrez IDs
  all_genes_mapped <- merge(all_genes, gene_mapping,
                            by.x = "ensembl_clean", by.y = "ENSEMBL",
                            all.x = FALSE)
  
  # Calculate ranking metric: delta × min(-log10(q), q_cap)
  # This combines effect size (direction + magnitude) with statistical confidence
  all_genes_mapped$neg_log10_q <- pmin(-log10(pmax(all_genes_mapped$qvalue, 1e-300)), q_cap)
  all_genes_mapped$rank_metric <- all_genes_mapped$cliffs_delta * all_genes_mapped$neg_log10_q
  
  # Create named vector
  gene_list <- all_genes_mapped$rank_metric
  names(gene_list) <- all_genes_mapped$ENTREZID
  
  # Sort by ranking metric (descending)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Remove duplicates
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
    pvalueCutoff = 1,      # Filter later
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
    pvalueCutoff = 1,      # Filter later
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
    pvalueCutoff = 1,      # Filter later
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE,
    eps = 0,
    seed = 1986
  )
}

# Function to perform GSEA for Hallmark
perform_hallmark_gsea <- function(gene_list, hallmark_t2g) {
  GSEA(
    geneList = gene_list,
    TERM2GENE = hallmark_t2g,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,      # Filter later
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE,
    eps = 0,
    seed = 1986
  )
}

# Function to filter and save GSEA results
save_gsea_results <- function(result, comparison, database, output_dir, fdr_cutoff = CONFIG$FDR_CUTOFF) {
  if (is.null(result) || nrow(result@result) == 0) {
    cat(sprintf("      No enriched terms for %s %s\n", comparison, database))
    return(NULL)
  }
  
  # Filter by FDR
  result_filtered <- result
  result_filtered@result <- result@result[result@result$p.adjust < fdr_cutoff, ]
  
  if (nrow(result_filtered@result) == 0) {
    cat(sprintf("      No significant terms (FDR < %g) for %s %s\n", fdr_cutoff, comparison, database))
    return(NULL)
  }
  
  # Create output subdirectory
  comp_dir <- paste0(output_dir, comparison, "/")
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save as CSV
  filename <- sprintf("%sGSEA_%s.csv", comp_dir, database)
  result_df <- as.data.frame(result_filtered)
  write.csv(result_df, filename, row.names = FALSE)
  
  cat(sprintf("      Saved %d enriched terms for %s: %s\n", 
              nrow(result_df), database, basename(filename)))
  
  return(result_filtered)
}

# Master function to run GSEA for a comparison
run_gsea_analysis <- function(comparison_name, thyr_deg_results, gene_mapping, hallmark_t2g, output_dir) {
  cat(sprintf("\n  Running GSEA for %s\n", comparison_name))
  cat("  ----------------------------------------\n")
  
  # Get the comparison data
  deg_data <- thyr_deg_results[[comparison_name]]
  
  # Prepare ranked gene list
  gene_list <- prepare_gsea_genelist(deg_data, gene_mapping)
  cat(sprintf("    Ranked gene list: %d genes\n", length(gene_list)))
  cat(sprintf("    Ranking metric: delta × min(-log10(q), %d)\n", CONFIG$Q_CAP))
  cat(sprintf("    Top gene: %s (score=%.2f)\n", names(gene_list)[1], gene_list[1]))
  cat(sprintf("    Bottom gene: %s (score=%.2f)\n", 
              names(gene_list)[length(gene_list)], gene_list[length(gene_list)]))
  
  # Store results
  
  gsea_results <- list()
  
  # 1. GO Biological Process
  cat("    Running GO BP GSEA...\n")
  gsea_results[["GO_BP"]] <- tryCatch({
    result <- perform_go_gsea(gene_list, "BP")
    save_gsea_results(result, comparison_name, "GO_BP", output_dir)
  }, error = function(e) {
    cat(sprintf("      Error in GO BP: %s\n", e$message))
    NULL
  })
  
  # 2. KEGG
  cat("    Running KEGG GSEA...\n")
  gsea_results[["KEGG"]] <- tryCatch({
    result <- perform_kegg_gsea(gene_list)
    save_gsea_results(result, comparison_name, "KEGG", output_dir)
  }, error = function(e) {
    cat(sprintf("      Error in KEGG: %s\n", e$message))
    NULL
  })
  
  # 3. Reactome
  cat("    Running Reactome GSEA...\n")
  gsea_results[["Reactome"]] <- tryCatch({
    result <- perform_reactome_gsea(gene_list)
    save_gsea_results(result, comparison_name, "Reactome", output_dir)
  }, error = function(e) {
    cat(sprintf("      Error in Reactome: %s\n", e$message))
    NULL
  })
  
  # 4. Hallmark
  cat("    Running Hallmark GSEA...\n")
  gsea_results[["Hallmark"]] <- tryCatch({
    result <- perform_hallmark_gsea(gene_list, hallmark_t2g)
    save_gsea_results(result, comparison_name, "Hallmark", output_dir)
  }, error = function(e) {
    cat(sprintf("      Error in Hallmark: %s\n", e$message))
    NULL
  })
  
  # Summary
  successful <- sum(!sapply(gsea_results, is.null))
  cat(sprintf("    GSEA Complete: %d/4 databases with significant results\n", successful))
  
  return(gsea_results)
}

cat("GSEA functions ready!\n")
cat("\nPhase 3 complete!\n")

# ============================================================================
# Phase 4: ORA Functions (Supplementary)
# ============================================================================

cat("\n=== PHASE 4: ORA FUNCTIONS (SUPPLEMENTARY) ===\n")

# Function to prepare gene lists for ORA
prepare_ora_genelist <- function(deg_results, direction = "all", gene_mapping) {
  # Get significant genes
  sig_genes <- deg_results$deg_summary$results_df[deg_results$deg_summary$results_df$significant, ]
  
  # Filter by direction
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
  
  unique(sig_genes_mapped$ENTREZID)
}

# ORA functions for each database
perform_go_ora <- function(gene_list, universe_list, ont = "BP") {
  enrichGO(
    gene = gene_list,
    universe = universe_list,
    OrgDb = org.Hs.eg.db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE
  )
}

perform_kegg_ora <- function(gene_list, universe_list) {
  enrichKEGG(
    gene = gene_list,
    universe = universe_list,
    organism = 'hsa',
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE
  )
}

perform_reactome_ora <- function(gene_list, universe_list) {
  enrichPathway(
    gene = gene_list,
    universe = universe_list,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE
  )
}

perform_hallmark_ora <- function(gene_list, universe_list, hallmark_t2g) {
  enricher(
    gene = gene_list,
    universe = universe_list,
    TERM2GENE = hallmark_t2g,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE
  )
}

# Save ORA results
save_ora_results <- function(result, comparison, database, direction, output_dir, fdr_cutoff = CONFIG$FDR_CUTOFF) {
  if (is.null(result) || nrow(result@result) == 0) {
    cat(sprintf("        No enriched terms for %s %s (%s)\n", comparison, database, direction))
    return(NULL)
  }
  
  # Filter by FDR
  result_filtered <- result
  result_filtered@result <- result@result[result@result$p.adjust < fdr_cutoff, ]
  
  if (nrow(result_filtered@result) == 0) {
    cat(sprintf("        No significant terms (FDR < %g) for %s %s (%s)\n", 
                fdr_cutoff, comparison, database, direction))
    return(NULL)
  }
  
  # Create output subdirectory
  comp_dir <- paste0(output_dir, comparison, "/")
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save as CSV
  filename <- sprintf("%sORA_%s_%s.csv", comp_dir, database, direction)
  result_df <- as.data.frame(result_filtered)
  write.csv(result_df, filename, row.names = FALSE)
  
  cat(sprintf("        Saved %d enriched terms for %s (%s): %s\n", 
              nrow(result_df), database, direction, basename(filename)))
  
  return(result_filtered)
}

# Master function for ORA
run_ora_analysis <- function(comparison_name, thyr_deg_results, gene_mapping, hallmark_t2g, output_dir) {
  cat(sprintf("\n  Running ORA for %s (Supplementary)\n", comparison_name))
  cat("  ----------------------------------------\n")
  
  deg_data <- thyr_deg_results[[comparison_name]]
  
  # Prepare universe
  all_genes_clean <- sub("\\..*", "", deg_data$deg_summary$results_df$gene_id)
  universe_df <- merge(data.frame(ensembl_clean = all_genes_clean), 
                       gene_mapping,
                       by.x = "ensembl_clean", by.y = "ENSEMBL",
                       all.x = FALSE)
  universe_list <- unique(universe_df$ENTREZID)
  
  cat(sprintf("    Universe size: %d genes\n", length(universe_list)))
  
  ora_results <- list()
  
  for (direction in c("up", "down", "all")) {
    gene_list <- prepare_ora_genelist(deg_data, direction, gene_mapping)
    dir_label <- if(direction == "all") "All" else toupper(direction)
    cat(sprintf("    %s-regulated genes: %d\n", dir_label, length(gene_list)))
    
    if (length(gene_list) < 5) {
      cat(sprintf("      Skipping %s (too few genes)\n", direction))
      next
    }
    
    # GO BP
    cat(sprintf("      Running GO BP (%s)...\n", direction))
    ora_results[[paste0("GO_BP_", direction)]] <- tryCatch({
      result <- perform_go_ora(gene_list, universe_list, "BP")
      save_ora_results(result, comparison_name, "GO_BP", direction, output_dir)
    }, error = function(e) {
      cat(sprintf("        Error: %s\n", e$message))
      NULL
    })
    
    # KEGG
    cat(sprintf("      Running KEGG (%s)...\n", direction))
    ora_results[[paste0("KEGG_", direction)]] <- tryCatch({
      result <- perform_kegg_ora(gene_list, universe_list)
      save_ora_results(result, comparison_name, "KEGG", direction, output_dir)
    }, error = function(e) {
      cat(sprintf("        Error: %s\n", e$message))
      NULL
    })
    
    # Reactome
    cat(sprintf("      Running Reactome (%s)...\n", direction))
    ora_results[[paste0("Reactome_", direction)]] <- tryCatch({
      result <- perform_reactome_ora(gene_list, universe_list)
      save_ora_results(result, comparison_name, "Reactome", direction, output_dir)
    }, error = function(e) {
      cat(sprintf("        Error: %s\n", e$message))
      NULL
    })
    
    # Hallmark
    cat(sprintf("      Running Hallmark (%s)...\n", direction))
    ora_results[[paste0("Hallmark_", direction)]] <- tryCatch({
      result <- perform_hallmark_ora(gene_list, universe_list, hallmark_t2g)
      save_ora_results(result, comparison_name, "Hallmark", direction, output_dir)
    }, error = function(e) {
      cat(sprintf("        Error: %s\n", e$message))
      NULL
    })
  }
  
  successful <- sum(!sapply(ora_results, is.null))
  cat(sprintf("    ORA Complete: %d analyses with significant results\n", successful))
  
  return(ora_results)
}

cat("ORA functions ready!\n")
cat("\nPhase 4 complete!\n")

# ============================================================================
# Phase 5: Execute Enrichment Analyses
# ============================================================================

cat("\n=== PHASE 5: EXECUTE ENRICHMENT ANALYSES ===\n")

# Check for existing results
results_file <- paste0(output_dir, "all_enrichment_results_v7.6.rds")
if (file.exists(results_file)) {
  cat("Previous results detected. Loading...\n")
  all_enrichment_results <- readRDS(results_file)
  cat("Loaded existing results.\n")
  cat("To force re-analysis, delete:", results_file, "\n")
  
} else {
  cat("Starting new analysis...\n\n")
  
  all_enrichment_results <- list()
  
  # ----------------------------------------------------------------------------
  # 1. GSEA for all comparisons (Primary Analysis)
  # ----------------------------------------------------------------------------
  
  cat("[PRIMARY] GSEA Analysis - All Comparisons\n")
  cat("==========================================\n")
  
  gsea_results_all <- list()
  
  # Only analyze configured comparisons (excluding B0_vs_B1_normal)
  for (comp_name in CONFIG$GSEA_COMPARISONS) {
    if (!comp_name %in% names(thyr_deg_results)) {
      cat(sprintf("  Warning: %s not found in DEG results, skipping\n", comp_name))
      next
    }
    gsea_results <- run_gsea_analysis(
      comparison_name = comp_name,
      thyr_deg_results = thyr_deg_results,
      gene_mapping = gene_mapping,
      hallmark_t2g = hallmark_t2g,
      output_dir = output_dir
    )
    gsea_results_all[[comp_name]] <- gsea_results
  }
  
  all_enrichment_results[["GSEA"]] <- gsea_results_all
  
  # ----------------------------------------------------------------------------
  # 2. ORA for R0_vs_R1_tumor only (Supplementary)
  # ----------------------------------------------------------------------------
  
  if ("R0_vs_R1_tumor" %in% names(thyr_deg_results)) {
    cat("\n[SUPPLEMENTARY] ORA Analysis - R0_vs_R1_tumor\n")
    cat("==============================================\n")
    
    ora_results <- run_ora_analysis(
      comparison_name = "R0_vs_R1_tumor",
      thyr_deg_results = thyr_deg_results,
      gene_mapping = gene_mapping,
      hallmark_t2g = hallmark_t2g,
      output_dir = output_dir
    )
    
    all_enrichment_results[["R0_vs_R1_tumor_ORA"]] <- ora_results
  }
  
  # Save results
  saveRDS(all_enrichment_results, results_file)
  cat("\nResults saved:", results_file, "\n")
}

cat("\nPhase 5 complete!\n")

# ============================================================================
# Phase 6: Cross-Comparison Analysis
# ============================================================================

cat("\n=== PHASE 6: CROSS-COMPARISON ANALYSIS ===\n")

# Check if required GSEA results exist
if ("GSEA" %in% names(all_enrichment_results)) {
  
  cross_comp_dir <- paste0(output_dir, "cross_comparison/")
  dir.create(cross_comp_dir, showWarnings = FALSE, recursive = TRUE)
  
  all_cross_comparison_results <- list()
  
  # Function to find consistent pathways
  find_consistent_pathways <- function(result1, result2, label1, label2, 
                                       db_name, fdr_cutoff = CONFIG$FDR_CUTOFF) {
    if (is.null(result1) || is.null(result2)) return(NULL)
    
    pathways1 <- result1@result[result1@result$p.adjust < fdr_cutoff, ]
    pathways2 <- result2@result[result2@result$p.adjust < fdr_cutoff, ]
    
    if (nrow(pathways1) == 0 || nrow(pathways2) == 0) return(NULL)
    
    # Merge by description
    common <- merge(
      pathways1[, c("ID", "Description", "NES", "p.adjust", "setSize")],
      pathways2[, c("ID", "Description", "NES", "p.adjust", "setSize")],
      by = "Description",
      suffixes = c(paste0("_", label1), paste0("_", label2))
    )
    
    if (nrow(common) == 0) return(NULL)
    
    # Check NES direction consistency
    nes_col1 <- paste0("NES_", label1)
    nes_col2 <- paste0("NES_", label2)
    common$NES_consistent <- sign(common[[nes_col1]]) == sign(common[[nes_col2]])
    common$avg_abs_NES <- (abs(common[[nes_col1]]) + abs(common[[nes_col2]])) / 2
    common$consistency_score <- common$avg_abs_NES * ifelse(common$NES_consistent, 1, -0.5)
    common$database <- db_name
    
    common[order(common$consistency_score, decreasing = TRUE), ]
  }
  
  databases <- c("GO_BP", "KEGG", "Reactome", "Hallmark")
  
  # ----------------------------------------------------------------------------
  # Cross-comparison 1: RET × BRAF (tumor vs tumor)
  # ----------------------------------------------------------------------------
  
  cat("\n--- Cross-comparison 1: RET × BRAF (tumor vs tumor) ---\n")
  cat("Purpose: Identify pathways affected by radiation in both driver backgrounds\n\n")
  
  if ("R0_vs_R1_tumor" %in% names(all_enrichment_results[["GSEA"]]) &&
      "B0_vs_B1_tumor" %in% names(all_enrichment_results[["GSEA"]])) {
    
    ret_braf_results <- list()
    
    for (db in databases) {
      cat(sprintf("  [%s] ", db))
      
      ret_result <- all_enrichment_results[["GSEA"]][["R0_vs_R1_tumor"]][[db]]
      braf_result <- all_enrichment_results[["GSEA"]][["B0_vs_B1_tumor"]][[db]]
      
      consistent <- find_consistent_pathways(ret_result, braf_result, "RET", "BRAF", db)
      
      if (!is.null(consistent) && nrow(consistent) > 0) {
        consistent_only <- consistent[consistent$NES_consistent, ]
        
        cat(sprintf("Common: %d, NES-consistent: %d\n", 
                    nrow(consistent), nrow(consistent_only)))
        
        if (nrow(consistent_only) > 0) {
          # Show top 3
          top3 <- head(consistent_only, 3)
          for (i in 1:nrow(top3)) {
            cat(sprintf("      %d. %s (NES: RET=%.2f, BRAF=%.2f)\n",
                        i, substr(top3$Description[i], 1, 45),
                        top3$NES_RET[i], top3$NES_BRAF[i]))
          }
          
          # Save
          filename <- paste0(cross_comp_dir, "RET_vs_BRAF_", db, "_consistent.csv")
          write.csv(consistent_only, filename, row.names = FALSE)
        }
        
        ret_braf_results[[db]] <- list(all = consistent, consistent = consistent_only)
      } else {
        cat("No common significant pathways\n")
      }
    }
    
    all_cross_comparison_results[["RET_vs_BRAF"]] <- ret_braf_results
    
  } else {
    cat("  Required comparisons not available\n")
  }
  
  # ----------------------------------------------------------------------------
  # Cross-comparison 2: Tumor × Normal (within RET)
  # ----------------------------------------------------------------------------
  
  cat("\n--- Cross-comparison 2: Tumor × Normal (within RET) ---\n")
  cat("Purpose: Identify tumor-specific vs systemic radiation effects\n\n")
  
  if ("R0_vs_R1_tumor" %in% names(all_enrichment_results[["GSEA"]]) &&
      "R0_vs_R1_normal" %in% names(all_enrichment_results[["GSEA"]])) {
    
    tumor_normal_results <- list()
    
    for (db in databases) {
      cat(sprintf("  [%s] ", db))
      
      tumor_result <- all_enrichment_results[["GSEA"]][["R0_vs_R1_tumor"]][[db]]
      normal_result <- all_enrichment_results[["GSEA"]][["R0_vs_R1_normal"]][[db]]
      
      consistent <- find_consistent_pathways(tumor_result, normal_result, "Tumor", "Normal", db)
      
      if (!is.null(consistent) && nrow(consistent) > 0) {
        consistent_only <- consistent[consistent$NES_consistent, ]
        
        cat(sprintf("Common: %d, NES-consistent: %d\n", 
                    nrow(consistent), nrow(consistent_only)))
        
        if (nrow(consistent_only) > 0) {
          # Show top 3
          top3 <- head(consistent_only, 3)
          for (i in 1:nrow(top3)) {
            cat(sprintf("      %d. %s (NES: Tumor=%.2f, Normal=%.2f)\n",
                        i, substr(top3$Description[i], 1, 45),
                        top3$NES_Tumor[i], top3$NES_Normal[i]))
          }
          
          # Save
          filename <- paste0(cross_comp_dir, "Tumor_vs_Normal_", db, "_consistent.csv")
          write.csv(consistent_only, filename, row.names = FALSE)
        }
        
        tumor_normal_results[[db]] <- list(all = consistent, consistent = consistent_only)
      } else {
        cat("No common significant pathways\n")
      }
    }
    
    all_cross_comparison_results[["Tumor_vs_Normal"]] <- tumor_normal_results
    
  } else {
    cat("  Required comparisons not available\n")
  }
  
  # ----------------------------------------------------------------------------
  # Create comprehensive Excel report
  # ----------------------------------------------------------------------------
  
  cat("\n--- Creating Excel report ---\n")
  
  wb <- createWorkbook()
  
  header_style <- createStyle(fontColour = "white", fgFill = "#4F81BD",
                              halign = "center", textDecoration = "bold")
  
  # Overview sheet
  addWorksheet(wb, "Overview")
  
  overview_data <- data.frame(
    Cross_Comparison = character(),
    Database = character(),
    Total_Common = integer(),
    NES_Consistent = integer(),
    Top_Pathway = character(),
    stringsAsFactors = FALSE
  )
  
  for (cc_name in names(all_cross_comparison_results)) {
    for (db in databases) {
      cc_result <- all_cross_comparison_results[[cc_name]][[db]]
      if (!is.null(cc_result)) {
        top_path <- if (nrow(cc_result$consistent) > 0) {
          substr(cc_result$consistent$Description[1], 1, 50)
        } else {
          "-"
        }
        overview_data <- rbind(overview_data, data.frame(
          Cross_Comparison = cc_name,
          Database = db,
          Total_Common = nrow(cc_result$all),
          NES_Consistent = nrow(cc_result$consistent),
          Top_Pathway = top_path,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  writeData(wb, "Overview", overview_data)
  addStyle(wb, "Overview", header_style, rows = 1, cols = 1:5, gridExpand = TRUE)
  setColWidths(wb, "Overview", cols = 1:5, widths = c(20, 12, 15, 15, 50))
  
  # Detail sheets for each cross-comparison
  for (cc_name in names(all_cross_comparison_results)) {
    for (db in databases) {
      cc_result <- all_cross_comparison_results[[cc_name]][[db]]
      if (!is.null(cc_result) && nrow(cc_result$consistent) > 0) {
        sheet_name <- paste0(substr(cc_name, 1, 10), "_", db)
        if (nchar(sheet_name) > 31) sheet_name <- substr(sheet_name, 1, 31)
        
        addWorksheet(wb, sheet_name)
        
        # Select columns based on cross-comparison type
        if (cc_name == "RET_vs_BRAF") {
          display_cols <- c("Description", "NES_RET", "p.adjust_RET", 
                            "NES_BRAF", "p.adjust_BRAF", "avg_abs_NES")
        } else {
          display_cols <- c("Description", "NES_Tumor", "p.adjust_Tumor", 
                            "NES_Normal", "p.adjust_Normal", "avg_abs_NES")
        }
        
        display_data <- cc_result$consistent[, display_cols]
        
        # Round numeric columns
        numeric_cols <- names(display_data)[sapply(display_data, is.numeric)]
        display_data[numeric_cols] <- lapply(display_data[numeric_cols], round, 4)
        
        writeData(wb, sheet_name, display_data)
        addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:ncol(display_data), gridExpand = TRUE)
        setColWidths(wb, sheet_name, cols = 1:ncol(display_data), widths = "auto")
      }
    }
  }
  
  excel_file <- paste0(cross_comp_dir, "GSEA_cross_comparison_results.xlsx")
  tryCatch({
    saveWorkbook(wb, excel_file, overwrite = TRUE)
    cat("  Saved:", basename(excel_file), "\n")
  }, error = function(e) {
    cat("  Warning: Could not save Excel file -", e$message, "\n")
  })
  
  # Update results
  all_enrichment_results[["cross_comparison"]] <- all_cross_comparison_results
  saveRDS(all_enrichment_results, results_file)
  
} else {
  cat("GSEA results not available for cross-comparison\n")
}

cat("\nPhase 6 complete!\n")

# ============================================================================
# Final Summary
# ============================================================================

cat("\n=== ENRICHMENT ANALYSIS COMPLETE (v7.6) ===\n")
cat("============================================\n")

cat("\nAnalysis scope:\n")
cat("  GSEA comparisons:", paste(CONFIG$GSEA_COMPARISONS, collapse = ", "), "\n")
cat("  Excluded: B0_vs_B1_normal (no biological signal)\n")

cat("\nMethodology:\n")
cat("  GSEA ranking: delta × min(-log10(q),", CONFIG$Q_CAP, ")\n")
cat("  FDR correction: BH method\n")
cat("  FDR threshold:", CONFIG$FDR_CUTOFF, "\n")

cat("\nDatabases analyzed:\n")
cat("  - GO Biological Process\n")
cat("  - KEGG Pathways\n")
cat("  - Reactome Pathways\n")
cat("  - MSigDB Hallmark (50 gene sets)\n")

cat("\nCross-comparisons:\n")
cat("  1. RET × BRAF: Common radiation effects across drivers\n")
cat("  2. Tumor × Normal: Tissue-specific radiation effects\n")

cat("\nOutput location:", output_dir, "\n")

# Cleanup
rm(list = setdiff(ls(), c("paths", "all_enrichment_results", "thyr_deg_results")))
gc()

cat("\n=== Script Complete ===\n")