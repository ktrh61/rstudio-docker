# 10_enrichment_analysis.R - GSEA Enrichment Analysis for R0_vs_R1_tumor
# Purpose: Perform GO/KEGG/Reactome/Hallmark GSEA on R0_vs_R1_tumor DEGs
# Method: GSEA with ranked gene list (Cliff's delta × -log10(q))
# Input: thyr_deg_results.rds (from 09_deg_analysis.R)
# Output: GSEA results with Excel report and visualizations
# Version: v7.10 - GSEA only, visualization can be re-generated independently
# Date: 2025-12-15

source("analysis_v7/setup.R")

cat("\n=== GSEA Enrichment Analysis (v7.10) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Target: R0_vs_R1_tumor\n")
cat("Databases: GO_BP, KEGG, Reactome, Hallmark\n")

# Load packages
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

set.seed(1986)
cat("  Random seed: 1986\n")

# ============================================================================
# Configuration
# ============================================================================

CONFIG <- list(
  FDR_CUTOFF = 0.05,
  Q_CAP = 10,
  MIN_GENESET_SIZE = 10,
  MAX_GENESET_SIZE = 500
)

cat("\nConfiguration:\n")
cat("  FDR correction: BH method\n")
cat("  FDR threshold:", CONFIG$FDR_CUTOFF, "\n")
cat("  GSEA ranking: delta × min(-log10(q),", CONFIG$Q_CAP, ")\n")
cat("  Gene set size:", CONFIG$MIN_GENESET_SIZE, "-", CONFIG$MAX_GENESET_SIZE, "\n")

# ============================================================================
# Loading DEG results
# ============================================================================

cat("\n--- Loading DEG results ---\n")

deg_results_path <- paste0(paths$processed, "thyr_deg_results.rds")
if (!file.exists(deg_results_path)) {
  stop("DEG results not found. Please run 09_deg_analysis.R first.")
}

deg_output <- readRDS(deg_results_path)
thyr_deg_results <- deg_output$deg_results

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
# Gene ID Mapping
# ============================================================================

cat("\n--- Gene ID Mapping ---\n")

all_gene_ids_clean <- sub("\\..*", "", results_df$gene_id)
cat("  Total genes:", length(all_gene_ids_clean), "\n")

gene_mapping <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = all_gene_ids_clean,
                                      columns = c("ENTREZID", "SYMBOL"),
                                      keytype = "ENSEMBL")

gene_mapping <- gene_mapping[!is.na(gene_mapping$ENTREZID), ]

n_multi <- sum(duplicated(gene_mapping$ENSEMBL))
if (n_multi > 0) {
  cat(sprintf("  Note: %d Ensembl IDs map to multiple Entrez IDs (using first)\n", n_multi))
}
gene_mapping <- gene_mapping[!duplicated(gene_mapping$ENSEMBL), ]

cat("  Successfully mapped:", nrow(gene_mapping), "genes\n")

results_df$ensembl_clean <- sub("\\..*", "", results_df$gene_id)
results_df_mapped <- merge(results_df, gene_mapping,
                           by.x = "ensembl_clean", by.y = "ENSEMBL",
                           all.x = FALSE)

cat(sprintf("  Genes with Entrez mapping: %d\n", nrow(results_df_mapped)))

# ============================================================================
# Prepare Hallmark Gene Sets
# ============================================================================

cat("\n--- Prepare Hallmark Gene Sets ---\n")

hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
cat("  Loaded MSigDB Hallmark:", length(unique(hallmark_df$gs_name)), "gene sets\n")

hallmark_t2g <- hallmark_df %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  dplyr::mutate(gs_name = gsub("^HALLMARK_", "", gs_name))

# ============================================================================
# GSEA Analysis
# ============================================================================

cat("\n--- GSEA Analysis ---\n")

# Prepare ranked gene list
results_df_mapped$neg_log10_q <- pmin(-log10(pmax(results_df_mapped$qvalue, 1e-300)), CONFIG$Q_CAP)
results_df_mapped$rank_metric <- results_df_mapped$cliffs_delta * results_df_mapped$neg_log10_q

gene_list <- results_df_mapped$rank_metric
names(gene_list) <- results_df_mapped$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]

cat(sprintf("  Prepared gene list: %d genes\n", length(gene_list)))
cat(sprintf("  Rank range: %.2f to %.2f\n", min(gene_list), max(gene_list)))

gsea_results <- list()

# GO Biological Process
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

gsea_go_bp_sig <- gsea_go_bp
if (nrow(gsea_go_bp@result) > 0) {
  gsea_go_bp_sig@result <- gsea_go_bp@result[gsea_go_bp@result$p.adjust < CONFIG$FDR_CUTOFF, ]
}
gsea_results[["GO_BP"]] <- gsea_go_bp_sig
cat(sprintf("    Significant pathways: %d\n", nrow(gsea_go_bp_sig@result)))

# KEGG
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

# Reactome
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

# Hallmark
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

# ============================================================================
# Saving results
# ============================================================================

cat("\n--- Saving results ---\n")

# Save RDS
thyr_enrichment_results <- list(
  GSEA = gsea_results,
  gene_list = gene_list,
  gene_mapping = gene_mapping,
  config = CONFIG,
  date = Sys.Date(),
  version = "v7.10"
)

results_file <- paste0(paths$processed, "thyr_enrichment_results.rds")
saveRDS(thyr_enrichment_results, results_file)
cat("  Saved:", basename(results_file), "\n")

# Create Excel report
cat("  Creating Excel report...\n")

wb <- createWorkbook()
header_style <- createStyle(fontColour = "white", fgFill = "#4F81BD",
                            halign = "center", textDecoration = "bold")

# Overview sheet
addWorksheet(wb, "Overview")
overview_data <- data.frame(
  Database = c("GO_BP", "KEGG", "Reactome", "Hallmark"),
  Significant_Pathways = sapply(c("GO_BP", "KEGG", "Reactome", "Hallmark"), function(db) {
    if (!is.null(gsea_results[[db]])) nrow(gsea_results[[db]]@result) else 0
  }),
  stringsAsFactors = FALSE
)
writeData(wb, "Overview", overview_data)
addStyle(wb, "Overview", header_style, rows = 1, cols = 1:2, gridExpand = TRUE)

# Detail sheets
for (db in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
  if (!is.null(gsea_results[[db]]) && nrow(gsea_results[[db]]@result) > 0) {
    addWorksheet(wb, db)
    result_df <- as.data.frame(gsea_results[[db]]@result)
    
    cols_to_keep <- intersect(
      c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "core_enrichment"),
      colnames(result_df)
    )
    result_df <- result_df[, cols_to_keep]
    
    numeric_cols <- names(result_df)[sapply(result_df, is.numeric)]
    result_df[numeric_cols] <- lapply(result_df[numeric_cols], function(x) round(x, 4))
    
    writeData(wb, db, result_df)
    addStyle(wb, db, header_style, rows = 1, cols = 1:ncol(result_df), gridExpand = TRUE)
    setColWidths(wb, db, cols = 1:ncol(result_df), widths = "auto")
  }
}

excel_file <- paste0(output_dir, "GSEA_R0_vs_R1_tumor.xlsx")
saveWorkbook(wb, excel_file, overwrite = TRUE)
cat("  Saved:", basename(excel_file), "\n")

# ============================================================================
# GSEA Processing Complete - Clean up before visualization
# ============================================================================

cat("\n=== GSEA Processing Complete ===\n")
cat("All results saved. Cleaning up intermediate variables...\n")

rm(list = setdiff(ls(), c("paths", "CONFIG")))
gc()

cat("Memory cleaned. Visualization can be run independently from here.\n")

# ============================================================================
# Generating Plots
# (This section can be re-run independently after GSEA processing)
# ============================================================================

cat("\n--- Generating Plots ---\n")

if (!exists("paths")) {
  source("analysis_v7/setup.R")
}

# Load saved results
cat("Loading saved GSEA results...\n")
enrichment_output <- readRDS(paste0(paths$processed, "thyr_enrichment_results.rds"))
gsea_results <- enrichment_output$GSEA

cat(sprintf("  Found %d databases\n", length(gsea_results)))

library(ggplot2)
library(enrichplot)

output_dir <- paste0(paths$output, "enrichment_analysis/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Dotplots
cat("\n  Creating dotplots...\n")

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

# Running score plots (top 5 per database)
cat("\n  Creating running score plots...\n")

for (db in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
  gsea_result <- gsea_results[[db]]
  if (!is.null(gsea_result) && nrow(gsea_result@result) > 0) {
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

# ============================================================================
# Final Summary
# ============================================================================

cat("\n=== GSEA Enrichment Analysis Complete (v7.10) ===\n")

cat("\nTarget: R0_vs_R1_tumor\n")

cat("\nGSEA Results:\n")
cat(sprintf("%-12s %8s\n", "Database", "Pathways"))
cat(paste(rep("-", 22), collapse = ""), "\n")
for (db in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
  n_sig <- if (!is.null(gsea_results[[db]])) nrow(gsea_results[[db]]@result) else 0
  cat(sprintf("%-12s %8d\n", db, n_sig))
}

cat("\nMethodology:\n")
cat("  Ranking: Cliff's delta × min(-log10(q), 10)\n")
cat("  FDR correction: BH method\n")
cat("  FDR threshold: 0.05\n")

cat("\nOutputs:\n")
cat("  Results: thyr_enrichment_results.rds\n")
cat("  Excel: GSEA_R0_vs_R1_tumor.xlsx\n")
cat("  Plots: GSEA_dotplot_*.pdf, GSEA_running_*.pdf\n")

cat("\nOutput location:", output_dir, "\n")

rm(list = setdiff(ls(), c("paths")))
gc()

cat("\n=== Script Complete ===\n")