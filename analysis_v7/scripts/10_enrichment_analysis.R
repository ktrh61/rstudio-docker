# 10_enrichment_analysis.R - GSEA Enrichment Analysis for R0_vs_R1_tumor
# Purpose: Perform GO/KEGG/Reactome/Hallmark GSEA on R0_vs_R1_tumor DEGs
# Method: GSEA with Cliff's delta ranking (p-value tiebreaker for reproducibility)
# Input: thyr_deg_results.rds (from 09_deg_analysis.R)
# Output: GSEA results with Excel report and visualizations
# Version: v7.11 - Delta-only ranking with p-value tiebreaker
# Date: 2025-12-15

source("analysis_v7/setup.R")

cat("\n=== GSEA Enrichment Analysis (v7.11) ===\n")
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
  MIN_GENESET_SIZE = 10,
  MAX_GENESET_SIZE = 500,
  TIEBREAKER_EPSILON = 1e-10
)

cat("\nConfiguration:\n")
cat("  FDR correction: BH method\n")
cat("  FDR threshold:", CONFIG$FDR_CUTOFF, "\n")
cat("  GSEA ranking: Cliff's delta\n")
cat("  Tiebreaker: p-value -> ENTREZID (deterministic)\n")
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
# Sorting priority: delta (desc) -> p-value (asc) -> ENTREZID (asc)
cat("  Sorting genes: delta -> p-value -> ENTREZID\n")

results_df_sorted <- results_df_mapped[order(
  -results_df_mapped$cliffs_delta,
  results_df_mapped$pvalue,
  results_df_mapped$ENTREZID
), ]

n <- nrow(results_df_sorted)

# Add negligible perturbation to preserve sorting order in GSEA
# Higher rank gets larger perturbation value
perturbation <- seq(CONFIG$TIEBREAKER_EPSILON, 0, length.out = n)

rank_metric <- results_df_sorted$cliffs_delta + perturbation

# Check for ties
n_base_ties <- sum(duplicated(results_df_sorted$cliffs_delta))
n_final_ties <- sum(duplicated(rank_metric))
cat(sprintf("  Base delta ties: %d\n", n_base_ties))
cat(sprintf("  Final ranking ties: %d\n", n_final_ties))

gene_list <- rank_metric
names(gene_list) <- results_df_sorted$ENTREZID
gene_list <- gene_list[!duplicated(names(gene_list))]

cat(sprintf("  Prepared gene list: %d genes\n", length(gene_list)))
cat(sprintf("  Rank range: %.6f to %.6f\n", min(gene_list), max(gene_list)))

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

# Apply simplify() to GO_BP and save
go_bp_simplified <- NULL
if (!is.null(gsea_results[["GO_BP"]]) && nrow(gsea_results[["GO_BP"]]@result) > 0) {
  cat("  Applying simplify() to GO_BP (cutoff=0.7, by=p.adjust)...\n")
  go_bp_simplified <- clusterProfiler::simplify(
    gsea_results[["GO_BP"]],
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  n_before <- nrow(gsea_results[["GO_BP"]]@result)
  n_after <- nrow(go_bp_simplified@result)
  cat(sprintf("    Reduced: %d -> %d pathways (%.1f%% reduction)\n", 
              n_before, n_after, (1 - n_after/n_before) * 100))
}

# Save RDS
thyr_enrichment_results <- list(
  GSEA = gsea_results,
  GO_BP_simplified = go_bp_simplified,
  gene_list = gene_list,
  gene_mapping = gene_mapping,
  config = CONFIG,
  date = Sys.Date(),
  version = "v7.11"
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

rm(list = setdiff(ls(), c("paths")))
gc()

cat("Memory cleaned. Visualization and summary can be run independently from here.\n")

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

cat("\nPlot generation complete.\n")

# ============================================================================
# Summary
# (This section can be re-run independently)
# ============================================================================

cat("\n--- Summary ---\n")

if (!exists("paths")) {
  source("analysis_v7/setup.R")
}

# Load saved results if not already loaded
if (!exists("enrichment_output")) {
  enrichment_output <- readRDS(paste0(paths$processed, "thyr_enrichment_results.rds"))
  gsea_results <- enrichment_output$GSEA
}

cat("\nTarget: R0_vs_R1_tumor\n")

# GSEA results with NES breakdown
cat("\nGSEA Results:\n")
cat(sprintf("%-12s %8s %8s %8s\n", "Database", "Total", "NES>0", "NES<0"))
cat(paste(rep("-", 40), collapse = ""), "\n")

for (db in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
  if (!is.null(gsea_results[[db]]) && nrow(gsea_results[[db]]@result) > 0) {
    result_df <- gsea_results[[db]]@result
    n_total <- nrow(result_df)
    n_pos <- sum(result_df$NES > 0, na.rm = TRUE)
    n_neg <- sum(result_df$NES < 0, na.rm = TRUE)
    cat(sprintf("%-12s %8d %8d %8d\n", db, n_total, n_pos, n_neg))
  } else {
    cat(sprintf("%-12s %8d %8d %8d\n", db, 0, 0, 0))
  }
}

# Top pathways per database
cat("\nTop 5 pathways by |NES| per database:\n")

for (db in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
  if (!is.null(gsea_results[[db]]) && nrow(gsea_results[[db]]@result) > 0) {
    cat(sprintf("\n  [%s]\n", db))
    result_df <- gsea_results[[db]]@result
    result_df <- result_df[order(abs(result_df$NES), decreasing = TRUE), ]
    top5 <- head(result_df, 5)
    
    for (i in 1:nrow(top5)) {
      direction <- ifelse(top5$NES[i] > 0, "+", "-")
      cat(sprintf("    %d. %s (NES=%s%.2f, q=%.2e)\n", 
                  i, 
                  substr(top5$Description[i], 1, 50),
                  direction,
                  abs(top5$NES[i]),
                  top5$p.adjust[i]))
    }
  }
}

# ============================================================================
# Detailed Results by NES Direction
# ============================================================================

cat("\n--- Detailed Results by NES Direction ---\n")

# Helper function to display top N pathways
display_top_pathways <- function(result_df, n = 10, direction_label = "") {
  if (nrow(result_df) == 0) {
    cat(sprintf("    No pathways found\n"))
    return()
  }
  
  n_show <- min(n, nrow(result_df))
  top_df <- head(result_df[order(result_df$p.adjust), ], n_show)
  
  for (i in 1:nrow(top_df)) {
    cat(sprintf("    %2d. %s\n", i, top_df$Description[i]))
    cat(sprintf("        NES=%.2f, p.adjust=%.2e, size=%d\n",
                top_df$NES[i], top_df$p.adjust[i], top_df$setSize[i]))
  }
}

# --- GO_BP with simplify() ---
cat("\n[GO_BP - with simplify()]\n")

go_bp_result <- gsea_results[["GO_BP"]]
go_bp_simplified <- enrichment_output$GO_BP_simplified

if (!is.null(go_bp_simplified) && nrow(go_bp_simplified@result) > 0) {
  
  simplified_df <- go_bp_simplified@result
  n_before <- nrow(go_bp_result@result)
  n_after <- nrow(simplified_df)
  cat(sprintf("  Using saved simplify() result (cutoff=0.7, by=p.adjust)\n"))
  cat(sprintf("  Reduced: %d -> %d pathways (%.1f%% reduction)\n", 
              n_before, n_after, (1 - n_after/n_before) * 100))
  
  # Split by NES direction
  pos_df <- simplified_df[simplified_df$NES > 0, ]
  neg_df <- simplified_df[simplified_df$NES < 0, ]
  
  cat(sprintf("  After simplify: NES>0: %d, NES<0: %d\n", nrow(pos_df), nrow(neg_df)))
  
  cat("\n  NES > 0 (Up in R1) - Top 10:\n")
  display_top_pathways(pos_df, 10)
  
  cat("\n  NES < 0 (Down in R1) - Top 10:\n")
  display_top_pathways(neg_df, 10)
  
} else {
  cat("  No GO_BP simplified results available\n")
}

# --- KEGG ---
cat("\n[KEGG]\n")

kegg_result <- gsea_results[["KEGG"]]
if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
  kegg_df <- kegg_result@result
  pos_df <- kegg_df[kegg_df$NES > 0, ]
  neg_df <- kegg_df[kegg_df$NES < 0, ]
  
  cat(sprintf("  NES>0: %d, NES<0: %d\n", nrow(pos_df), nrow(neg_df)))
  
  cat("\n  NES > 0 (Up in R1) - Top 10:\n")
  display_top_pathways(pos_df, 10)
  
  cat("\n  NES < 0 (Down in R1) - Top 10:\n")
  display_top_pathways(neg_df, 10)
  
} else {
  cat("  No KEGG results available\n")
}

# --- Reactome ---
cat("\n[Reactome]\n")

reactome_result <- gsea_results[["Reactome"]]
if (!is.null(reactome_result) && nrow(reactome_result@result) > 0) {
  reactome_df <- reactome_result@result
  pos_df <- reactome_df[reactome_df$NES > 0, ]
  neg_df <- reactome_df[reactome_df$NES < 0, ]
  
  cat(sprintf("  NES>0: %d, NES<0: %d\n", nrow(pos_df), nrow(neg_df)))
  
  cat("\n  NES > 0 (Up in R1) - Top 10:\n")
  display_top_pathways(pos_df, 10)
  
  cat("\n  NES < 0 (Down in R1) - Top 10:\n")
  display_top_pathways(neg_df, 10)
  
} else {
  cat("  No Reactome results available\n")
}

# --- Hallmark ---
cat("\n[Hallmark]\n")

hallmark_result <- gsea_results[["Hallmark"]]
if (!is.null(hallmark_result) && nrow(hallmark_result@result) > 0) {
  hallmark_df <- hallmark_result@result
  pos_df <- hallmark_df[hallmark_df$NES > 0, ]
  neg_df <- hallmark_df[hallmark_df$NES < 0, ]
  
  cat(sprintf("  NES>0: %d, NES<0: %d\n", nrow(pos_df), nrow(neg_df)))
  
  cat("\n  NES > 0 (Up in R1) - Top 10:\n")
  display_top_pathways(pos_df, 10)
  
  cat("\n  NES < 0 (Down in R1) - Top 10:\n")
  display_top_pathways(neg_df, 10)
  
} else {
  cat("  No Hallmark results available\n")
}

# ============================================================================
# NES Direction Dotplots (Per Database)
# ============================================================================

cat("\n--- NES Direction Dotplots ---\n")

# Ensure output_dir is defined
output_dir <- paste0(paths$output, "enrichment_analysis/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Helper function to extract top N pathways by NES direction
extract_top_pathways <- function(result_df, n_each = 5) {
  if (is.null(result_df) || nrow(result_df) == 0) {
    return(NULL)
  }
  
  # NES > 0: top n by p.adjust
  pos_df <- result_df[result_df$NES > 0, ]
  if (nrow(pos_df) > 0) {
    pos_df <- head(pos_df[order(pos_df$p.adjust), ], n_each)
  }
  
  # NES < 0: top n by p.adjust
  neg_df <- result_df[result_df$NES < 0, ]
  if (nrow(neg_df) > 0) {
    neg_df <- head(neg_df[order(neg_df$p.adjust), ], n_each)
  }
  
  rbind(pos_df, neg_df)
}

# Helper function to create dotplot for a single database
create_nes_dotplot <- function(plot_df, db_name, output_dir) {
  if (is.null(plot_df) || nrow(plot_df) == 0) {
    cat(sprintf("    %s: No data available\n", db_name))
    return(NULL)
  }
  
  # Wrap description at 50 characters
  plot_df$Description_wrapped <- stringr::str_wrap(plot_df$Description, width = 50)
  
  # Create the plot
  p <- ggplot(plot_df, aes(x = NES, 
                           y = reorder(Description_wrapped, NES),
                           size = setSize,
                           color = p.adjust)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_gradient(low = "red", high = "blue", 
                         trans = "log10",
                         name = "FDR\n(p.adjust)") +
    scale_size_continuous(name = "Gene Set\nSize", range = c(3, 10)) +
    labs(
      title = paste0("GSEA: ", db_name, " (R0 vs R1 Tumor)"),
      subtitle = "Top 5 pathways by FDR for each NES direction",
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  # Display in RStudio
  print(p)
  
  # Save PDF
  plot_file <- paste0(output_dir, "GSEA_NES_dotplot_", db_name, ".pdf")
  ggsave(plot_file, p, width = 10, height = 6)
  
  n_pos <- sum(plot_df$NES > 0)
  n_neg <- sum(plot_df$NES < 0)
  cat(sprintf("    %s: %d pathways (NES>0: %d, NES<0: %d) -> Saved\n", 
              db_name, nrow(plot_df), n_pos, n_neg))
  
  return(invisible(p))
}

cat("  Creating individual dotplots per database...\n\n")

# GO_BP (use simplify() result)
go_bp_top <- NULL
if (!is.null(go_bp_simplified) && nrow(go_bp_simplified@result) > 0) {
  go_bp_top <- extract_top_pathways(go_bp_simplified@result, 5)
  create_nes_dotplot(go_bp_top, "GO_BP", output_dir)
}

# KEGG
kegg_top <- NULL
if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
  kegg_top <- extract_top_pathways(kegg_result@result, 5)
  create_nes_dotplot(kegg_top, "KEGG", output_dir)
}

# Reactome
reactome_top <- NULL
if (!is.null(reactome_result) && nrow(reactome_result@result) > 0) {
  reactome_top <- extract_top_pathways(reactome_result@result, 5)
  create_nes_dotplot(reactome_top, "Reactome", output_dir)
}

# Hallmark
hallmark_top <- NULL
if (!is.null(hallmark_result) && nrow(hallmark_result@result) > 0) {
  hallmark_top <- extract_top_pathways(hallmark_result@result, 5)
  create_nes_dotplot(hallmark_top, "Hallmark", output_dir)
}

cat("\n  Dotplot generation complete.\n")

# ============================================================================
# GO_BP Emapplot (Semantic Similarity Network)
# ============================================================================

cat("\n--- GO_BP Emapplot ---\n")

if (!is.null(go_bp_simplified) && nrow(go_bp_simplified@result) > 0) {
  
  result_df <- go_bp_simplified@result
  
  # --- NES > 0 (Up in R1) ---
  pos_df <- result_df[result_df$NES > 0, ]
  if (nrow(pos_df) >= 5) {
    cat(sprintf("  NES > 0: %d pathways available\n", nrow(pos_df)))
    
    # Select top 50 by p.adjust
    n_show_pos <- min(50, nrow(pos_df))
    top_pos_ids <- head(pos_df[order(pos_df$p.adjust), "ID"], n_show_pos)
    
    # Create subset object FIRST, then calculate pairwise similarity
    go_bp_pos <- go_bp_simplified
    go_bp_pos@result <- result_df[result_df$ID %in% top_pos_ids, ]
    
    cat(sprintf("    Calculating pairwise similarity for %d pathways...\n", n_show_pos))
    go_bp_pos <- pairwise_termsim(go_bp_pos)
    
    cat(sprintf("    Plotting...\n"))
    
    p_pos <- emapplot(go_bp_pos, 
                      showCategory = n_show_pos,
                      color = "NES",
                      layout = "kk") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                           midpoint = 0, name = "NES") +
      labs(title = "GO_BP: NES > 0 (Up in R1)",
           subtitle = sprintf("Top %d pathways by FDR", n_show_pos)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5))
    
    print(p_pos)
    
    emapplot_pos_file <- paste0(output_dir, "GSEA_emapplot_GO_BP_positive.pdf")
    ggsave(emapplot_pos_file, p_pos, width = 10, height = 10)
    cat(sprintf("    Saved: %s\n", basename(emapplot_pos_file)))
    
  } else {
    cat("  NES > 0: Not enough pathways for emapplot\n")
  }
  
  # --- NES < 0 (Down in R1) ---
  neg_df <- result_df[result_df$NES < 0, ]
  if (nrow(neg_df) >= 5) {
    cat(sprintf("  NES < 0: %d pathways available\n", nrow(neg_df)))
    
    # Select top 50 by p.adjust
    n_show_neg <- min(50, nrow(neg_df))
    top_neg_ids <- head(neg_df[order(neg_df$p.adjust), "ID"], n_show_neg)
    
    # Create subset object FIRST, then calculate pairwise similarity
    go_bp_neg <- go_bp_simplified
    go_bp_neg@result <- result_df[result_df$ID %in% top_neg_ids, ]
    
    cat(sprintf("    Calculating pairwise similarity for %d pathways...\n", n_show_neg))
    go_bp_neg <- pairwise_termsim(go_bp_neg)
    
    cat(sprintf("    Plotting...\n"))
    
    p_neg <- emapplot(go_bp_neg,
                      showCategory = n_show_neg,
                      color = "NES",
                      layout = "kk") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = 0, name = "NES") +
      labs(title = "GO_BP: NES < 0 (Down in R1)",
           subtitle = sprintf("Top %d pathways by FDR", n_show_neg)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5))
    
    print(p_neg)
    
    emapplot_neg_file <- paste0(output_dir, "GSEA_emapplot_GO_BP_negative.pdf")
    ggsave(emapplot_neg_file, p_neg, width = 10, height = 10)
    cat(sprintf("    Saved: %s\n", basename(emapplot_neg_file)))
    
  } else {
    cat("  NES < 0: Not enough pathways for emapplot\n")
  }
  
  cat("\n  Emapplot generation complete.\n")
  
} else {
  cat("  No GO_BP simplified results available for emapplot\n")
}

# ============================================================================
# Leading-Edge Heatmap (Immune/Inflammatory Pathways)
# ============================================================================

cat("\n--- Leading-Edge Heatmap ---\n")

# Target pathways (immune/inflammatory theme, NES < 0)
target_pathways <- c(
  "GO:0032609",  # type II interferon production
  "GO:0032640",  # tumor necrosis factor production
  "GO:0071621",  # granulocyte chemotaxis
  "GO:0001909"   # leukocyte mediated cytotoxicity
)

# Get GO_BP result (use original, not simplified, for core_enrichment)
go_bp_full <- gsea_results[["GO_BP"]]

if (!is.null(go_bp_full) && nrow(go_bp_full@result) > 0) {
  
  go_bp_df <- go_bp_full@result
  
  # Check which pathways are available
  available_pathways <- target_pathways[target_pathways %in% go_bp_df$ID]
  missing_pathways <- target_pathways[!target_pathways %in% go_bp_df$ID]
  
  cat(sprintf("  Target pathways found: %d/%d\n", 
              length(available_pathways), length(target_pathways)))
  
  if (length(missing_pathways) > 0) {
    cat("  Missing pathways:\n")
    for (mp in missing_pathways) {
      cat(sprintf("    - %s\n", mp))
    }
  }
  
  if (length(available_pathways) >= 2) {
    
    # Step 1: Extract leading-edge genes from each pathway
    cat("  Extracting leading-edge genes...\n")
    
    pathway_genes <- list()
    for (pathway_id in available_pathways) {
      pathway_row <- go_bp_df[go_bp_df$ID == pathway_id, ]
      pathway_name <- pathway_row$Description
      
      # core_enrichment contains "/" separated ENTREZID
      core_genes <- unlist(strsplit(pathway_row$core_enrichment, "/"))
      pathway_genes[[pathway_id]] <- core_genes
      
      cat(sprintf("    %s: %d genes\n", pathway_name, length(core_genes)))
    }
    
    # Step 2: Count gene occurrences across pathways
    all_genes <- unlist(pathway_genes)
    gene_counts <- table(all_genes)
    gene_counts_sorted <- sort(gene_counts, decreasing = TRUE)
    
    cat(sprintf("  Total unique genes: %d\n", length(gene_counts)))
    cat(sprintf("  Genes in 2+ pathways: %d\n", sum(gene_counts >= 2)))
    cat(sprintf("  Genes in 3+ pathways: %d\n", sum(gene_counts >= 3)))
    
    # Step 3: Select genes (prioritize multi-pathway, limit to 30)
    # First take genes in 2+ pathways, then fill with single-pathway genes
    multi_pathway_genes <- names(gene_counts_sorted[gene_counts_sorted >= 2])
    single_pathway_genes <- names(gene_counts_sorted[gene_counts_sorted == 1])
    
    n_target <- 30
    if (length(multi_pathway_genes) >= n_target) {
      selected_genes <- head(multi_pathway_genes, n_target)
    } else {
      n_fill <- n_target - length(multi_pathway_genes)
      selected_genes <- c(multi_pathway_genes, head(single_pathway_genes, n_fill))
    }
    
    cat(sprintf("  Selected genes: %d\n", length(selected_genes)))
    
    # Step 4: Convert ENTREZID to SYMBOL
    gene_symbols <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = selected_genes,
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )
    
    # Remove NA mappings
    valid_idx <- !is.na(gene_symbols)
    selected_genes <- selected_genes[valid_idx]
    gene_symbols <- gene_symbols[valid_idx]
    
    cat(sprintf("  Genes with valid symbols: %d\n", length(gene_symbols)))
    
    # Step 5: Convert ENTREZID to ENSEMBL using gene_mapping from enrichment_output
    cat("  Converting ENTREZID to ENSEMBL...\n")
    
    gene_mapping <- enrichment_output$gene_mapping
    
    # Create ENTREZID -> ENSEMBL lookup
    entrez_to_ensembl <- gene_mapping$ENSEMBL
    names(entrez_to_ensembl) <- gene_mapping$ENTREZID
    
    # Convert selected genes
    selected_ensembl <- entrez_to_ensembl[selected_genes]
    valid_ensembl <- !is.na(selected_ensembl)
    
    selected_genes <- selected_genes[valid_ensembl]
    selected_ensembl <- selected_ensembl[valid_ensembl]
    gene_symbols <- gene_symbols[valid_ensembl]
    
    cat(sprintf("  Genes with valid ENSEMBL: %d\n", length(selected_ensembl)))
    
    # Step 6: Load normalized expression data
    cat("  Loading normalized expression data...\n")
    
    dgelist_file <- paste0(paths$processed, "analysis_dgelist_R0_vs_R1_tumor.rds")
    
    if (file.exists(dgelist_file)) {
      dgelist <- readRDS(dgelist_file)
      
      # Calculate normalized CPM
      norm_cpm <- edgeR::cpm(dgelist, normalized.lib.sizes = TRUE, log = TRUE)
      
      # DGEList rownames have version numbers (e.g., ENSG00000000003.15)
      # Create lookup: base ENSEMBL -> versioned ENSEMBL
      dgelist_ensembl_base <- sub("\\.\\d+$", "", rownames(norm_cpm))
      names(dgelist_ensembl_base) <- rownames(norm_cpm)
      versioned_lookup <- setNames(names(dgelist_ensembl_base), dgelist_ensembl_base)
      
      # Match genes using base ENSEMBL IDs
      matched_versioned <- versioned_lookup[selected_ensembl]
      valid_match <- !is.na(matched_versioned)
      
      available_versioned <- matched_versioned[valid_match]
      available_symbols <- gene_symbols[valid_match]
      available_entrez <- selected_genes[valid_match]
      
      cat(sprintf("  Genes matched in expression data: %d\n", length(available_versioned)))
      
      if (length(available_versioned) >= 10) {
        
        # Subset expression matrix
        expr_matrix <- norm_cpm[available_versioned, , drop = FALSE]
        
        # Use gene symbols as rownames
        rownames(expr_matrix) <- available_symbols
        
        cat(sprintf("  Expression matrix: %d genes x %d samples\n", 
                    nrow(expr_matrix), ncol(expr_matrix)))
        
        # Step 7: Gene-wise Z-score transformation
        expr_zscore <- t(scale(t(expr_matrix)))
        
        # Step 8: Prepare sample annotation
        sample_groups <- as.character(dgelist$samples$group)
        names(sample_groups) <- rownames(dgelist$samples)
        
        # Order samples: R0 first, then R1
        r0_samples <- names(sample_groups)[sample_groups == "R0"]
        r1_samples <- names(sample_groups)[sample_groups == "R1"]
        sample_order <- c(r0_samples, r1_samples)
        
        expr_zscore_ordered <- expr_zscore[, sample_order, drop = FALSE]
        
        cat(sprintf("  Sample order: R0 (%d) -> R1 (%d)\n", 
                    length(r0_samples), length(r1_samples)))
        
        # Step 9: Create heatmap with ComplexHeatmap
        cat("  Creating heatmap...\n")
        
        library(ComplexHeatmap)
        library(circlize)
        
        # Column annotation
        col_annotation <- HeatmapAnnotation(
          Group = factor(sample_groups[sample_order], levels = c("R0", "R1")),
          col = list(Group = c("R0" = "#3498db", "R1" = "#e74c3c")),
          annotation_name_side = "left"
        )
        
        # Color scale for Z-score (-3 to +3 to cover full data range)
        col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
        
        # Row annotation: pathway membership count
        gene_membership <- gene_counts[available_entrez]
        names(gene_membership) <- available_symbols
        gene_membership_ordered <- gene_membership[rownames(expr_zscore_ordered)]
        
        row_annotation <- rowAnnotation(
          Pathways = anno_barplot(
            as.numeric(gene_membership_ordered),
            width = unit(1.5, "cm"),
            gp = gpar(fill = "#95a5a6")
          ),
          annotation_name_rot = 0
        )
        
        # Create heatmap
        ht <- Heatmap(
          expr_zscore_ordered,
          name = "Z-score",
          col = col_fun,
          
          # Clustering
          cluster_rows = TRUE,
          cluster_columns = FALSE,
          
          # Annotations
          top_annotation = col_annotation,
          right_annotation = row_annotation,
          
          # Labels
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          
          # Title
          column_title = "Leading-edge genes from immune/inflammatory pathways\ndownregulated in R1",
          column_title_gp = gpar(fontsize = 11, fontface = "bold"),
          
          # Legend
          heatmap_legend_param = list(
            title = "Z-score",
            at = c(-3, -2, -1, 0, 1, 2, 3),
            labels = c("-3", "-2", "-1", "0", "1", "2", "3")
          )
        )
        
        # Display
        draw(ht)
        
        # Save PDF
        heatmap_file <- paste0(output_dir, "GSEA_leading_edge_heatmap.pdf")
        pdf(heatmap_file, width = 10, height = 8)
        draw(ht)
        dev.off()
        cat(sprintf("  Saved: %s\n", basename(heatmap_file)))
        
        # Report pathway info
        cat("\n  Pathways included:\n")
        for (pathway_id in available_pathways) {
          pathway_row <- go_bp_df[go_bp_df$ID == pathway_id, ]
          cat(sprintf("    - %s (NES=%.2f)\n", 
                      pathway_row$Description, pathway_row$NES))
        }
        
        cat("\n  Leading-edge heatmap complete.\n")
        
      } else {
        cat("  Not enough genes found in expression data\n")
      }
      
    } else {
      cat(sprintf("  DGEList file not found: %s\n", basename(dgelist_file)))
    }
    
  } else {
    cat("  Not enough target pathways found\n")
  }
  
} else {
  cat("  No GO_BP results available\n")
}

# Configuration from saved data
cat("\nMethodology:\n")
cat("  Ranking: Cliff's delta\n")
cat("  Tiebreaker: p-value (asc) -> ENTREZID (asc)\n")
cat("  Perturbation: 1e-10 scale to preserve order\n")
cat("  FDR correction: BH method\n")
cat("  FDR threshold:", enrichment_output$config$FDR_CUTOFF, "\n")
cat("  Gene set size:", enrichment_output$config$MIN_GENESET_SIZE, "-", 
    enrichment_output$config$MAX_GENESET_SIZE, "\n")
cat("  Version:", enrichment_output$version, "\n")
cat("  Date:", as.character(enrichment_output$date), "\n")

cat("\nOutputs:\n")
cat("  Results: thyr_enrichment_results.rds\n")
cat("  Excel: GSEA_R0_vs_R1_tumor.xlsx\n")
cat("  Plots: GSEA_NES_dotplot_*.pdf, GSEA_running_*.pdf, GSEA_emapplot_*.pdf\n")
cat("         GSEA_leading_edge_heatmap.pdf\n")

cat("\nOutput location:", output_dir, "\n")

# ============================================================================
# Script Complete
# ============================================================================

rm(list = setdiff(ls(), c("paths")))
gc()

cat("\n=== Script Complete ===\n")