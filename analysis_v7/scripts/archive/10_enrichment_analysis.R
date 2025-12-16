# 10_enrichment_analysis.R - Enrichment Analysis for R0_vs_R1_tumor
# Purpose: Perform GO/KEGG/Reactome/Hallmark enrichment analysis on R0_vs_R1_tumor DEGs
# Method: GSEA (4 databases) + ORA (3 databases × 3 patterns)
# Input: thyr_deg_results.rds (from 09_deg_analysis.R)
# Output: Enrichment results with Excel reports and visualizations
# Version: v7.9 - R0_vs_R1_tumor specialized
# Date: 2025-12-10

source("analysis_v7/setup.R")

cat("\n=== Enrichment Analysis for R0_vs_R1_tumor (v7.9) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Target: R0_vs_R1_tumor only\n")
cat("GSEA: GO_BP, KEGG, Reactome, Hallmark\n")
cat("ORA: GO_BP, KEGG, Reactome (Up/Down/All)\n")

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
  
  # Output control
  VERBOSE = TRUE              # Verbose output
)

cat("\nConfiguration:\n")
cat("  FDR correction: BH method\n")
cat("  FDR threshold:", CONFIG$FDR_CUTOFF, "\n")
cat("  GSEA ranking: delta x min(-log10(q),", CONFIG$Q_CAP, ")\n")
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

# Extract R0_vs_R1_tumor only
if (!"R0_vs_R1_tumor" %in% names(thyr_deg_results)) {
  stop("R0_vs_R1_tumor not found in DEG results.")
}

r0r1_tumor <- thyr_deg_results[["R0_vs_R1_tumor"]]
results_df <- r0r1_tumor$deg_summary$results_df
summary_stats <- r0r1_tumor$deg_summary$summary_stats

cat("  R0_vs_R1_tumor loaded\n")
cat(sprintf("    Total genes: %d\n", summary_stats$total_genes_tested))
cat(sprintf("    Significant DEGs: %d\n", summary_stats$significant_genes))
cat(sprintf("    Upregulated: %d\n", summary_stats$upregulated))
cat(sprintf("    Downregulated: %d\n", summary_stats$downregulated))

# Create output directory
output_dir <- paste0(paths$output, "enrichment_analysis/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Phase 1: Gene ID Mapping
# ============================================================================

cat("\n=== PHASE 1: GENE ID MAPPING ===\n")

# Clean Ensembl IDs (remove version)
all_gene_ids_clean <- sub("\\..*", "", results_df$gene_id)
cat("  Total genes:", length(all_gene_ids_clean), "\n")

# Map Ensembl to Entrez
gene_mapping <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = all_gene_ids_clean,
                                      columns = c("ENTREZID", "SYMBOL"),
                                      keytype = "ENSEMBL")

# Remove NA mappings
gene_mapping <- gene_mapping[!is.na(gene_mapping$ENTREZID), ]

# Handle 1:many mappings (use first match)
n_multi <- sum(duplicated(gene_mapping$ENSEMBL))
if (n_multi > 0) {
  cat(sprintf("  Note: %d Ensembl IDs map to multiple Entrez IDs (using first match)\n", n_multi))
}
gene_mapping <- gene_mapping[!duplicated(gene_mapping$ENSEMBL), ]

cat("  Successfully mapped:", nrow(gene_mapping), "genes\n")
cat("  Unmapped genes:", length(all_gene_ids_clean) - nrow(gene_mapping), "\n")

# Add clean Ensembl ID to results_df
results_df$ensembl_clean <- sub("\\..*", "", results_df$gene_id)

# Merge mapping
results_df_mapped <- merge(results_df, gene_mapping,
                           by.x = "ensembl_clean", by.y = "ENSEMBL",
                           all.x = FALSE)

cat(sprintf("  Genes with Entrez mapping: %d\n", nrow(results_df_mapped)))

cat("\nPhase 1 complete!\n")

# ============================================================================
# Phase 2: Prepare Hallmark Gene Sets
# ============================================================================

cat("\n=== PHASE 2: PREPARE HALLMARK GENE SETS ===\n")

# Get Hallmark gene sets
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
cat("  Loaded MSigDB Hallmark:", length(unique(hallmark_df$gs_name)), "gene sets\n")

# Create TERM2GENE format for clusterProfiler
hallmark_t2g <- hallmark_df %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  dplyr::mutate(gs_name = gsub("^HALLMARK_", "", gs_name))

cat("\nPhase 2 complete!\n")

# ============================================================================
# Phase 3: GSEA Analysis
# ============================================================================

cat("\n=== PHASE 3: GSEA ANALYSIS ===\n")

# Prepare ranked gene list for GSEA
# Ranking: cliffs_delta × min(-log10(qvalue), Q_CAP)
results_df_mapped$neg_log10_q <- pmin(-log10(pmax(results_df_mapped$qvalue, 1e-300)), CONFIG$Q_CAP)
results_df_mapped$rank_metric <- results_df_mapped$cliffs_delta * results_df_mapped$neg_log10_q

# Create named vector
gene_list <- results_df_mapped$rank_metric
names(gene_list) <- results_df_mapped$ENTREZID

# Sort by ranking metric (descending)
gene_list <- sort(gene_list, decreasing = TRUE)

# Remove duplicates
gene_list <- gene_list[!duplicated(names(gene_list))]

cat(sprintf("  Prepared gene list: %d genes\n", length(gene_list)))
cat(sprintf("  Rank range: %.2f to %.2f\n", min(gene_list), max(gene_list)))

# Storage for GSEA results
gsea_results <- list()

# --- GO Biological Process GSEA ---
cat("\n  Running GO BP GSEA...\n")
gsea_go_bp <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  minGSSize = CONFIG$MIN_GENESET_SIZE,
  maxGSSize = CONFIG$MAX_GENESET_SIZE,
  eps = 0,
  seed = 1986
)

# Filter by FDR
gsea_go_bp_sig <- gsea_go_bp
if (nrow(gsea_go_bp@result) > 0) {
  gsea_go_bp_sig@result <- gsea_go_bp@result[gsea_go_bp@result$p.adjust < CONFIG$FDR_CUTOFF, ]
}
gsea_results[["GO_BP"]] <- gsea_go_bp_sig
cat(sprintf("    Significant pathways: %d\n", nrow(gsea_go_bp_sig@result)))

# --- KEGG GSEA ---
cat("\n  Running KEGG GSEA...\n")
gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  minGSSize = CONFIG$MIN_GENESET_SIZE,
  maxGSSize = CONFIG$MAX_GENESET_SIZE,
  eps = 0,
  seed = 1986
)

gsea_kegg_sig <- gsea_kegg
if (nrow(gsea_kegg@result) > 0) {
  gsea_kegg_sig@result <- gsea_kegg@result[gsea_kegg@result$p.adjust < CONFIG$FDR_CUTOFF, ]
}
gsea_results[["KEGG"]] <- gsea_kegg_sig
cat(sprintf("    Significant pathways: %d\n", nrow(gsea_kegg_sig@result)))

# --- Reactome GSEA ---
cat("\n  Running Reactome GSEA...\n")
gsea_reactome <- gsePathway(
  geneList = gene_list,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  minGSSize = CONFIG$MIN_GENESET_SIZE,
  maxGSSize = CONFIG$MAX_GENESET_SIZE,
  eps = 0,
  seed = 1986
)

gsea_reactome_sig <- gsea_reactome
if (nrow(gsea_reactome@result) > 0) {
  gsea_reactome_sig@result <- gsea_reactome@result[gsea_reactome@result$p.adjust < CONFIG$FDR_CUTOFF, ]
}
gsea_results[["Reactome"]] <- gsea_reactome_sig
cat(sprintf("    Significant pathways: %d\n", nrow(gsea_reactome_sig@result)))

# --- Hallmark GSEA ---
cat("\n  Running Hallmark GSEA...\n")
gsea_hallmark <- GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark_t2g,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  minGSSize = CONFIG$MIN_GENESET_SIZE,
  maxGSSize = CONFIG$MAX_GENESET_SIZE,
  eps = 0,
  seed = 1986
)

gsea_hallmark_sig <- gsea_hallmark
if (nrow(gsea_hallmark@result) > 0) {
  gsea_hallmark_sig@result <- gsea_hallmark@result[gsea_hallmark@result$p.adjust < CONFIG$FDR_CUTOFF, ]
}
gsea_results[["Hallmark"]] <- gsea_hallmark_sig
cat(sprintf("    Significant pathways: %d\n", nrow(gsea_hallmark_sig@result)))

cat("\nPhase 3 complete!\n")

# ============================================================================
# Phase 4: ORA Analysis
# ============================================================================

cat("\n=== PHASE 4: ORA ANALYSIS ===\n")

# Prepare DEG lists
sig_genes <- results_df_mapped[results_df_mapped$significant, ]
up_genes <- sig_genes[sig_genes$direction == "UP", "ENTREZID"]
down_genes <- sig_genes[sig_genes$direction == "DOWN", "ENTREZID"]
all_degs <- sig_genes$ENTREZID

# Background: all tested genes with Entrez mapping
universe <- results_df_mapped$ENTREZID

cat(sprintf("  DEGs with Entrez mapping:\n"))
cat(sprintf("    Up: %d\n", length(up_genes)))
cat(sprintf("    Down: %d\n", length(down_genes)))
cat(sprintf("    All: %d\n", length(all_degs)))
cat(sprintf("    Universe: %d\n", length(universe)))

# Storage for ORA results
ora_results <- list(Up = list(), Down = list(), All = list())

# Function to perform ORA for all databases
perform_ora <- function(gene_list, universe, label) {
  results <- list()
  
  # GO BP
  cat(sprintf("    %s: GO BP...\n", label))
  ora_go <- enrichGO(
    gene = gene_list,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE
  )
  if (!is.null(ora_go) && nrow(ora_go@result) > 0) {
    ora_go@result <- ora_go@result[ora_go@result$p.adjust < CONFIG$FDR_CUTOFF, ]
  }
  results[["GO_BP"]] <- ora_go
  
  # KEGG
  cat(sprintf("    %s: KEGG...\n", label))
  ora_kegg <- enrichKEGG(
    gene = gene_list,
    universe = universe,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE
  )
  if (!is.null(ora_kegg) && nrow(ora_kegg@result) > 0) {
    ora_kegg@result <- ora_kegg@result[ora_kegg@result$p.adjust < CONFIG$FDR_CUTOFF, ]
  }
  results[["KEGG"]] <- ora_kegg
  
  # Reactome
  cat(sprintf("    %s: Reactome...\n", label))
  ora_reactome <- enrichPathway(
    gene = gene_list,
    universe = universe,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = CONFIG$MIN_GENESET_SIZE,
    maxGSSize = CONFIG$MAX_GENESET_SIZE
  )
  if (!is.null(ora_reactome) && nrow(ora_reactome@result) > 0) {
    ora_reactome@result <- ora_reactome@result[ora_reactome@result$p.adjust < CONFIG$FDR_CUTOFF, ]
  }
  results[["Reactome"]] <- ora_reactome
  
  return(results)
}

# Perform ORA for each direction
cat("\n  Running ORA for Up-regulated genes...\n")
ora_results[["Up"]] <- perform_ora(up_genes, universe, "Up")

cat("\n  Running ORA for Down-regulated genes...\n")
ora_results[["Down"]] <- perform_ora(down_genes, universe, "Down")

cat("\n  Running ORA for All DEGs...\n")
ora_results[["All"]] <- perform_ora(all_degs, universe, "All")

# Summary
cat("\n  ORA Results Summary:\n")
for (direction in c("Up", "Down", "All")) {
  for (db in c("GO_BP", "KEGG", "Reactome")) {
    n_sig <- 0
    if (!is.null(ora_results[[direction]][[db]])) {
      n_sig <- nrow(ora_results[[direction]][[db]]@result)
    }
    cat(sprintf("    %s-%s: %d\n", direction, db, n_sig))
  }
}

cat("\nPhase 4 complete!\n")

# ============================================================================
# Phase 5: Save Results and Create Excel Reports
# ============================================================================

cat("\n=== PHASE 5: SAVE RESULTS ===\n")

# Save RDS
all_enrichment_results <- list(
  GSEA = gsea_results,
  ORA = ora_results,
  gene_mapping = gene_mapping,
  config = CONFIG
)

results_file <- paste0(paths$processed, "thyr_enrichment_results.rds")
saveRDS(all_enrichment_results, results_file)
cat("  Saved:", basename(results_file), "\n")

# --- Create GSEA Excel Report ---
cat("\n  Creating GSEA Excel report...\n")

wb_gsea <- createWorkbook()
header_style <- createStyle(fontColour = "white", fgFill = "#4F81BD",
                            halign = "center", textDecoration = "bold")

# Overview sheet
addWorksheet(wb_gsea, "Overview")
overview_data <- data.frame(
  Database = c("GO_BP", "KEGG", "Reactome", "Hallmark"),
  Significant_Pathways = sapply(c("GO_BP", "KEGG", "Reactome", "Hallmark"), function(db) {
    if (!is.null(gsea_results[[db]])) nrow(gsea_results[[db]]@result) else 0
  }),
  stringsAsFactors = FALSE
)
writeData(wb_gsea, "Overview", overview_data)
addStyle(wb_gsea, "Overview", header_style, rows = 1, cols = 1:2, gridExpand = TRUE)

# Detail sheets for each database
for (db in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
  if (!is.null(gsea_results[[db]]) && nrow(gsea_results[[db]]@result) > 0) {
    addWorksheet(wb_gsea, db)
    result_df <- as.data.frame(gsea_results[[db]]@result)
    
    # Select key columns
    cols_to_keep <- intersect(
      c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "core_enrichment"),
      colnames(result_df)
    )
    result_df <- result_df[, cols_to_keep]
    
    # Round numeric columns
    numeric_cols <- names(result_df)[sapply(result_df, is.numeric)]
    result_df[numeric_cols] <- lapply(result_df[numeric_cols], function(x) round(x, 4))
    
    writeData(wb_gsea, db, result_df)
    addStyle(wb_gsea, db, header_style, rows = 1, cols = 1:ncol(result_df), gridExpand = TRUE)
    setColWidths(wb_gsea, db, cols = 1:ncol(result_df), widths = "auto")
  }
}

gsea_excel_file <- paste0(output_dir, "GSEA_R0_vs_R1_tumor.xlsx")
saveWorkbook(wb_gsea, gsea_excel_file, overwrite = TRUE)
cat("  Saved:", basename(gsea_excel_file), "\n")

# --- Create ORA Excel Report ---
cat("\n  Creating ORA Excel report...\n")

wb_ora <- createWorkbook()

# Overview sheet
addWorksheet(wb_ora, "Overview")
ora_overview <- data.frame(
  Direction = rep(c("Up", "Down", "All"), each = 3),
  Database = rep(c("GO_BP", "KEGG", "Reactome"), 3),
  Significant_Terms = sapply(1:9, function(i) {
    dir <- c("Up", "Down", "All")[ceiling(i/3)]
    db <- c("GO_BP", "KEGG", "Reactome")[((i-1) %% 3) + 1]
    if (!is.null(ora_results[[dir]][[db]])) nrow(ora_results[[dir]][[db]]@result) else 0
  }),
  stringsAsFactors = FALSE
)
writeData(wb_ora, "Overview", ora_overview)
addStyle(wb_ora, "Overview", header_style, rows = 1, cols = 1:3, gridExpand = TRUE)

# Detail sheets
for (direction in c("Up", "Down", "All")) {
  for (db in c("GO_BP", "KEGG", "Reactome")) {
    ora_result <- ora_results[[direction]][[db]]
    if (!is.null(ora_result) && nrow(ora_result@result) > 0) {
      sheet_name <- paste0(direction, "_", db)
      addWorksheet(wb_ora, sheet_name)
      
      result_df <- as.data.frame(ora_result@result)
      
      # Select key columns
      cols_to_keep <- intersect(
        c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"),
        colnames(result_df)
      )
      result_df <- result_df[, cols_to_keep]
      
      # Round numeric columns
      numeric_cols <- names(result_df)[sapply(result_df, is.numeric)]
      result_df[numeric_cols] <- lapply(result_df[numeric_cols], function(x) round(x, 4))
      
      writeData(wb_ora, sheet_name, result_df)
      addStyle(wb_ora, sheet_name, header_style, rows = 1, cols = 1:ncol(result_df), gridExpand = TRUE)
      setColWidths(wb_ora, sheet_name, cols = 1:ncol(result_df), widths = "auto")
    }
  }
}

ora_excel_file <- paste0(output_dir, "ORA_R0_vs_R1_tumor.xlsx")
saveWorkbook(wb_ora, ora_excel_file, overwrite = TRUE)
cat("  Saved:", basename(ora_excel_file), "\n")

cat("\nPhase 5 complete!\n")

# ============================================================================
# Phase 6: Visualization
# ============================================================================

cat("\n=== PHASE 6: VISUALIZATION ===\n")

# --- GSEA Dotplots ---
cat("\n  Creating GSEA dotplots...\n")

for (db in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
  gsea_result <- gsea_results[[db]]
  if (!is.null(gsea_result) && nrow(gsea_result@result) > 0) {
    n_show <- min(20, nrow(gsea_result@result))
    
    p <- dotplot(gsea_result, showCategory = n_show, split = ".sign") +
      facet_grid(.~.sign) +
      ggtitle(paste0("GSEA: ", db, " (R0_vs_R1_tumor)")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(
      filename = paste0(output_dir, "GSEA_dotplot_", db, ".pdf"),
      plot = p, width = 12, height = 8
    )
    cat(sprintf("    Saved: GSEA_dotplot_%s.pdf\n", db))
  }
}

# --- ORA Dotplots ---
cat("\n  Creating ORA dotplots...\n")

for (direction in c("Up", "Down", "All")) {
  for (db in c("GO_BP", "KEGG", "Reactome")) {
    ora_result <- ora_results[[direction]][[db]]
    if (!is.null(ora_result) && nrow(ora_result@result) > 0) {
      n_show <- min(20, nrow(ora_result@result))
      
      p <- dotplot(ora_result, showCategory = n_show) +
        ggtitle(paste0("ORA: ", db, " - ", direction, " DEGs (R0_vs_R1_tumor)")) +
        theme(plot.title = element_text(hjust = 0.5))
      
      ggsave(
        filename = paste0(output_dir, "ORA_dotplot_", direction, "_", db, ".pdf"),
        plot = p, width = 10, height = 8
      )
      cat(sprintf("    Saved: ORA_dotplot_%s_%s.pdf\n", direction, db))
    }
  }
}

# --- GSEA Running Score Plots (Top 5 per database) ---
cat("\n  Creating GSEA running score plots...\n")

for (db in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
  gsea_result <- gsea_results[[db]]
  if (!is.null(gsea_result) && nrow(gsea_result@result) > 0) {
    # Get top 5 by absolute NES
    result_df <- gsea_result@result
    result_df <- result_df[order(abs(result_df$NES), decreasing = TRUE), ]
    top_ids <- head(result_df$ID, 5)
    
    for (i in seq_along(top_ids)) {
      pathway_id <- top_ids[i]
      p <- gseaplot2(gsea_result, geneSetID = pathway_id, title = pathway_id)
      
      ggsave(
        filename = paste0(output_dir, "GSEA_running_", db, "_top", i, ".pdf"),
        plot = p, width = 10, height = 6
      )
    }
    cat(sprintf("    Saved: GSEA_running_%s_top1-5.pdf\n", db))
  }
}

cat("\nPhase 6 complete!\n")

# ============================================================================
# Final Summary
# ============================================================================

cat("\n=== ENRICHMENT ANALYSIS COMPLETE (v7.9) ===\n")
cat("============================================\n")

cat("\nTarget: R0_vs_R1_tumor\n")
cat(sprintf("  Total DEGs: %d (Up: %d, Down: %d)\n",
            summary_stats$significant_genes,
            summary_stats$upregulated,
            summary_stats$downregulated))

cat("\n--- GSEA Results ---\n")
cat(sprintf("%-12s %8s\n", "Database", "Pathways"))
cat(paste(rep("-", 22), collapse = ""), "\n")
for (db in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
  n_sig <- if (!is.null(gsea_results[[db]])) nrow(gsea_results[[db]]@result) else 0
  cat(sprintf("%-12s %8d\n", db, n_sig))
}

cat("\n--- ORA Results ---\n")
cat(sprintf("%-8s %-12s %8s\n", "Direction", "Database", "Terms"))
cat(paste(rep("-", 32), collapse = ""), "\n")
for (direction in c("Up", "Down", "All")) {
  for (db in c("GO_BP", "KEGG", "Reactome")) {
    n_sig <- if (!is.null(ora_results[[direction]][[db]])) nrow(ora_results[[direction]][[db]]@result) else 0
    cat(sprintf("%-8s %-12s %8d\n", direction, db, n_sig))
  }
}

cat("\nMethodology:\n")
cat("  GSEA ranking: delta × min(-log10(q),", CONFIG$Q_CAP, ")\n")
cat("  FDR correction: BH method\n")
cat("  FDR threshold:", CONFIG$FDR_CUTOFF, "\n")

cat("\nDatabases:\n")
cat("  - GO Biological Process\n")
cat("  - KEGG Pathways\n")
cat("  - Reactome Pathways\n")
cat("  - MSigDB Hallmark (GSEA only)\n")

cat("\nOutput files:\n")
cat("  - thyr_enrichment_results.rds (processed/)\n")
cat("  - GSEA_R0_vs_R1_tumor.xlsx\n")
cat("  - ORA_R0_vs_R1_tumor.xlsx\n")
cat("  - GSEA_dotplot_*.pdf\n")
cat("  - ORA_dotplot_*.pdf\n")
cat("  - GSEA_running_*.pdf\n")

cat("\nOutput location:", output_dir, "\n")

# Cleanup
rm(list = setdiff(ls(), c("paths", "all_enrichment_results")))
gc()

cat("\n=== Script Complete ===\n")