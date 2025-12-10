# check_enrichment_results.R - Review Enrichment Analysis Results
# Purpose: Summarize and inspect GSEA/ORA results from 10_enrichment_analysis.R

source("analysis_v7/setup.R")

cat("\n=== Enrichment Results Overview ===\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# Load results
results_path <- paste0(paths$processed, "thyr_enrichment_results.rds")
if (!file.exists(results_path)) {
  stop("Results not found. Run 10_enrichment_analysis.R first.")
}

enrichment <- readRDS(results_path)

# ============================================================================
# GSEA Results Summary
# ============================================================================

cat("=== GSEA RESULTS ===\n\n")

gsea <- enrichment$GSEA

for (db in names(gsea)) {
  result <- gsea[[db]]
  if (!is.null(result) && nrow(result@result) > 0) {
    df <- as.data.frame(result@result)
    
    cat(sprintf("--- %s (%d pathways) ---\n", db, nrow(df)))
    
    # Separate by NES direction
    up_pathways <- df[df$NES > 0, ]
    down_pathways <- df[df$NES < 0, ]
    
    cat(sprintf("  Activated (NES > 0): %d\n", nrow(up_pathways)))
    cat(sprintf("  Suppressed (NES < 0): %d\n", nrow(down_pathways)))
    
    # All pathways sorted by absolute NES
    df_sorted <- df[order(abs(df$NES), decreasing = TRUE), ]
    cat("\n  All pathways (sorted by |NES|):\n")
    for (i in 1:nrow(df_sorted)) {
      desc <- substr(df_sorted$Description[i], 1, 55)
      if (nchar(df_sorted$Description[i]) > 55) desc <- paste0(desc, "...")
      direction <- ifelse(df_sorted$NES[i] > 0, "UP", "DN")
      cat(sprintf("    %2d. [%s] %s\n", i, direction, desc))
      cat(sprintf("        NES=%.2f, p.adj=%.2e, size=%d\n",
                  df_sorted$NES[i], df_sorted$p.adjust[i], df_sorted$setSize[i]))
    }
    cat("\n")
  } else {
    cat(sprintf("--- %s: No significant pathways ---\n\n", db))
  }
}

# ============================================================================
# ORA Results Summary
# ============================================================================

cat("\n=== ORA RESULTS ===\n\n")

ora <- enrichment$ORA

for (direction in c("Up", "Down", "All")) {
  cat(sprintf("--- %s-regulated DEGs ---\n", direction))
  
  for (db in c("GO_BP", "KEGG", "Reactome")) {
    result <- ora[[direction]][[db]]
    if (!is.null(result) && nrow(result@result) > 0) {
      df <- as.data.frame(result@result)
      cat(sprintf("  %s: %d terms\n", db, nrow(df)))
      
      # All terms sorted by p.adjust
      df_sorted <- df[order(df$p.adjust), ]
      for (i in 1:nrow(df_sorted)) {
        desc <- substr(df_sorted$Description[i], 1, 50)
        if (nchar(df_sorted$Description[i]) > 50) desc <- paste0(desc, "...")
        cat(sprintf("    %2d. %s\n", i, desc))
        cat(sprintf("        p.adj=%.2e, count=%d, GeneRatio=%s\n",
                    df_sorted$p.adjust[i], df_sorted$Count[i], df_sorted$GeneRatio[i]))
      }
    } else {
      cat(sprintf("  %s: 0 terms\n", db))
    }
  }
  cat("\n")
}

# ============================================================================
# Biological Interpretation Hints
# ============================================================================

cat("\n=== BIOLOGICAL THEMES ===\n\n")

# Extract key terms from GSEA results
if (nrow(gsea[["GO_BP"]]@result) > 0) {
  go_df <- as.data.frame(gsea[["GO_BP"]]@result)
  
  # Keywords of interest for radiation biology
  keywords <- c("DNA", "repair", "damage", "immune", "inflamm", "apoptosis", 
                "cell cycle", "proliferat", "differenti", "metabol", "signal",
                "response", "stress", "oxidat", "thyroid", "hormone")
  
  cat("Key pathway categories in GO BP:\n")
  for (kw in keywords) {
    matches <- go_df[grepl(kw, go_df$Description, ignore.case = TRUE), ]
    if (nrow(matches) > 0) {
      cat(sprintf("  '%s': %d pathways\n", kw, nrow(matches)))
      # Show all matches
      for (i in 1:nrow(matches)) {
        direction <- ifelse(matches$NES[i] > 0, "UP", "DOWN")
        cat(sprintf("    - %s (%s, NES=%.2f)\n", 
                    substr(matches$Description[i], 1, 50), direction, matches$NES[i]))
      }
    }
  }
}

# Hallmark interpretation
cat("\n\nHallmark pathways detected:\n")
if (nrow(gsea[["Hallmark"]]@result) > 0) {
  hm_df <- as.data.frame(gsea[["Hallmark"]]@result)
  for (i in 1:nrow(hm_df)) {
    direction <- ifelse(hm_df$NES[i] > 0, "ACTIVATED", "SUPPRESSED")
    cat(sprintf("  - %s: %s (NES=%.2f, p.adj=%.2e)\n",
                hm_df$Description[i], direction, hm_df$NES[i], hm_df$p.adjust[i]))
  }
} else {
  cat("  None\n")
}

cat("\n=== Check Complete ===\n")