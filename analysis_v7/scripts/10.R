# 10_cross_driver_analysis.R - Cross-driver Gene Analysis
# Purpose: Annotate and analyze 43 cross-driver consistent genes
# Input: thyr_deg_results.rds, gene_mapping from enrichment analysis
# Output: Annotated gene list with functional categories
# Version: v1.0 - Standalone analysis for cross-driver genes
# Date: 2025-01-22

source("analysis_v7/setup.R")

cat("\n=== Cross-driver Gene Analysis (43 genes) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# Load packages
suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(dplyr)
  library(openxlsx)
})

# ============================================================================
# Load data
# ============================================================================

cat("\n--- Loading DEG results ---\n")

# Load DEG results
deg_results_path <- paste0(paths$processed, "thyr_deg_results.rds")
if (!file.exists(deg_results_path)) {
  stop("DEG results not found. Please run 09_deg_analysis.R first.")
}

deg_output <- readRDS(deg_results_path)
consistency_results <- deg_output$consistency_results

# Extract cross-driver genes
if (!"cross_driver_tumor" %in% names(consistency_results)) {
  stop("Cross-driver analysis not found in consistency results.")
}

cross_driver_genes <- consistency_results$cross_driver_tumor$consistent_genes
cat("Cross-driver genes found:", length(cross_driver_genes), "\n")

# ============================================================================
# Get gene information from DEG results
# ============================================================================

cat("\n--- Extracting gene information ---\n")

# Get R0_vs_R1 and B0_vs_B1 tumor results
r0r1_tumor <- deg_output$deg_results$R0_vs_R1_tumor$deg_summary$results_df
b0b1_tumor <- deg_output$deg_results$B0_vs_B1_tumor$deg_summary$results_df

# Filter for cross-driver genes
r0r1_cross <- r0r1_tumor[r0r1_tumor$gene_id %in% cross_driver_genes, ]
b0b1_cross <- b0b1_tumor[b0b1_tumor$gene_id %in% cross_driver_genes, ]

# Merge information
cross_driver_df <- merge(
  r0r1_cross[, c("gene_id", "gene_name", "gene_type", "log2FC", "qvalue")],
  b0b1_cross[, c("gene_id", "log2FC", "qvalue")],
  by = "gene_id",
  suffixes = c("_RET", "_BRAF")
)

cat("Merged data for", nrow(cross_driver_df), "genes\n")

# ============================================================================
# Gene annotation
# ============================================================================

cat("\n--- Annotating genes ---\n")

# Clean Ensembl IDs for mapping
cross_driver_df$ensembl_clean <- sub("\\..*", "", cross_driver_df$gene_id)

# Get Entrez IDs and additional annotation
annotation <- tryCatch({
  suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db,
                          keys = cross_driver_df$ensembl_clean,
                          columns = c("ENTREZID", "SYMBOL", "GENENAME", "GO", "PATH"),
                          keytype = "ENSEMBL")
  )
}, error = function(e) {
  cat("Warning: Some annotations failed\n")
  suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db,
                          keys = cross_driver_df$ensembl_clean,
                          columns = c("ENTREZID", "SYMBOL", "GENENAME"),
                          keytype = "ENSEMBL")
  )
})

# Merge with cross-driver data
cross_driver_annotated <- merge(cross_driver_df, annotation,
                                by.x = "ensembl_clean", by.y = "ENSEMBL",
                                all.x = TRUE)

# Remove duplicate rows from multiple GO/PATH annotations
cross_driver_annotated <- cross_driver_annotated[!duplicated(cross_driver_annotated$gene_id), ]

# ============================================================================
# Functional categorization
# ============================================================================

cat("\n--- Categorizing genes by function ---\n")

# Define categories based on gene names and known functions
categorize_gene <- function(gene_name, genename) {
  gene_lower <- tolower(paste(gene_name, genename))
  
  if (grepl("immune|interferon|interleukin|cytokine|chemokine|hla|cd[0-9]", gene_lower)) {
    return("Immune response")
  } else if (grepl("dna|repair|damage|brca|rad|xrcc|parp", gene_lower)) {
    return("DNA damage/repair")
  } else if (grepl("cell cycle|cyclin|cdk|aurora|plk|chk", gene_lower)) {
    return("Cell cycle")
  } else if (grepl("apoptosis|casp|bcl|bax|fas|trail", gene_lower)) {
    return("Apoptosis")
  } else if (grepl("transcription|zinc finger|znf|sox|fox|myc", gene_lower)) {
    return("Transcription regulation")
  } else if (grepl("kinase|phosphatase|signal|mapk|akt|pi3k", gene_lower)) {
    return("Signal transduction")
  } else if (grepl("metabol|enzyme|synthase|dehydrogenase", gene_lower)) {
    return("Metabolism")
  } else if (grepl("transport|channel|pump|carrier", gene_lower)) {
    return("Transport")
  } else if (grepl("ribosom|translation|eif|rpl|rps", gene_lower)) {
    return("Translation")
  } else if (grepl("mitochondri|oxidative|respiratory", gene_lower)) {
    return("Mitochondrial")
  } else {
    return("Other/Unknown")
  }
}

# Apply categorization
cross_driver_annotated$functional_category <- mapply(
  categorize_gene,
  cross_driver_annotated$gene_name,
  cross_driver_annotated$GENENAME
)

# ============================================================================
# Calculate biomarker potential score
# ============================================================================

cat("\n--- Evaluating biomarker potential ---\n")

# Simple scoring based on:
# 1. Consistency of fold change direction
# 2. Magnitude of fold change
# 3. Statistical significance

cross_driver_annotated$fc_consistency <- with(cross_driver_annotated,
                                              ifelse(sign(log2FC_RET) == sign(log2FC_BRAF), 1, 0))

cross_driver_annotated$avg_abs_fc <- with(cross_driver_annotated,
                                          (abs(log2FC_RET) + abs(log2FC_BRAF)) / 2)

cross_driver_annotated$combined_significance <- with(cross_driver_annotated,
                                                     -log10(qvalue_RET * qvalue_BRAF))

# Normalize scores
normalize_score <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

cross_driver_annotated$biomarker_score <- 
  normalize_score(cross_driver_annotated$fc_consistency) * 0.3 +
  normalize_score(cross_driver_annotated$avg_abs_fc) * 0.4 +
  normalize_score(cross_driver_annotated$combined_significance) * 0.3

# Sort by biomarker score
cross_driver_annotated <- cross_driver_annotated[order(cross_driver_annotated$biomarker_score, decreasing = TRUE), ]

# ============================================================================
# Summary statistics
# ============================================================================

cat("\n--- Summary ---\n")

# Functional category distribution
category_table <- table(cross_driver_annotated$functional_category)
cat("\nFunctional categories:\n")
print(sort(category_table, decreasing = TRUE))

# Direction consistency
cat("\nFold change direction:\n")
cat("  Consistent:", sum(cross_driver_annotated$fc_consistency), "\n")
cat("  Inconsistent:", sum(!cross_driver_annotated$fc_consistency), "\n")

# Top biomarker candidates
cat("\nTop 10 biomarker candidates:\n")
top10 <- head(cross_driver_annotated[, c("gene_name", "log2FC_RET", "log2FC_BRAF", "biomarker_score", "functional_category")], 10)
print(top10)

# ============================================================================
# Save results
# ============================================================================

cat("\n--- Saving results ---\n")

# Create output directory
output_dir <- paste0(paths$output, "cross_driver_analysis/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save as CSV
write.csv(cross_driver_annotated,
          paste0(output_dir, "cross_driver_43genes_annotated.csv"),
          row.names = FALSE)

# Save as Excel with formatting
wb <- createWorkbook()
addWorksheet(wb, "Cross-driver Genes")
writeData(wb, "Cross-driver Genes", cross_driver_annotated)

# Add formatting
headerStyle <- createStyle(textDecoration = "bold", halign = "center")
addStyle(wb, "Cross-driver Genes", headerStyle, rows = 1, cols = 1:ncol(cross_driver_annotated))

saveWorkbook(wb, paste0(output_dir, "cross_driver_43genes_annotated.xlsx"), overwrite = TRUE)

cat("  Saved: cross_driver_43genes_annotated.csv\n")
cat("  Saved: cross_driver_43genes_annotated.xlsx\n")

# Save R object
saveRDS(cross_driver_annotated, paste0(output_dir, "cross_driver_annotated.rds"))
cat("  Saved: cross_driver_annotated.rds\n")

cat("\nCross-driver analysis complete!\n")