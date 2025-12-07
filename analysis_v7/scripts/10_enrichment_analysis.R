# ============================================================================
# GSEA結果の詳細確認スクリプト
# 10_enrichment_analysis.R v7.8 の結果を確認するための直打ち用スクリプト
# 更新日: 2024-12-07
# ============================================================================

library(org.Hs.eg.db)

# --- Entrez ID → Gene Symbol 変換関数 ---
entrez_to_symbol <- function(entrez_ids) {
  if (is.character(entrez_ids) && length(entrez_ids) == 1 && grepl("/", entrez_ids)) {
    entrez_ids <- strsplit(entrez_ids, "/")[[1]]
  }
  symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, 
                    column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  symbols[is.na(symbols)] <- entrez_ids[is.na(symbols)]  # NAの場合は元のIDを保持
  return(symbols)
}

# --- 結果の読み込み ---
results <- readRDS("analysis_v7/output/enrichment_analysis/all_enrichment_results_v7.8.rds")

# ============================================================================
# NESの方向別集計（全体像の把握）
# ============================================================================
cat("\n")
cat("=======================================================\n")
cat("  NES Direction Summary: R0_vs_R1_tumor\n")
cat("  (NES > 0: R1でUP, NES < 0: R1でDOWN)\n")
cat("=======================================================\n\n")

for (db in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
  res <- results$GSEA$R0_vs_R1_tumor[[db]]@result
  sig <- res[res$p.adjust < 0.05, ]
  up <- sum(sig$NES > 0)
  down <- sum(sig$NES < 0)
  cat(sprintf("  %-12s: UP(R1 high)=%3d, DOWN(R1 low)=%3d (total=%d)\n", 
              db, up, down, up + down))
}

# ============================================================================
# R0_vs_R1_tumor: Hallmark（全結果 - UP/DOWN分離）
# ============================================================================
cat("\n")
cat("=======================================================\n")
cat("  R0_vs_R1_tumor: Hallmark - All Significant\n")
cat("=======================================================\n")

hallmark_ret <- results$GSEA$R0_vs_R1_tumor$Hallmark@result
hallmark_ret_sig <- hallmark_ret[hallmark_ret$p.adjust < 0.05, ]

# UP (R1 > R0)
hallmark_up <- hallmark_ret_sig[hallmark_ret_sig$NES > 0, ]
hallmark_up <- hallmark_up[order(-hallmark_up$NES), ]
cat("\n--- UP in R1 (high radiation) ---\n")
if (nrow(hallmark_up) > 0) {
  print(hallmark_up[, c("Description", "NES", "p.adjust", "setSize")], row.names = FALSE)
} else {
  cat("  (None)\n")
}

# DOWN (R1 < R0)
hallmark_down <- hallmark_ret_sig[hallmark_ret_sig$NES < 0, ]
hallmark_down <- hallmark_down[order(hallmark_down$NES), ]
cat("\n--- DOWN in R1 (high radiation) ---\n")
if (nrow(hallmark_down) > 0) {
  print(hallmark_down[, c("Description", "NES", "p.adjust", "setSize")], row.names = FALSE)
} else {
  cat("  (None)\n")
}

# ============================================================================
# R0_vs_R1_tumor: KEGG（全結果 - UP/DOWN分離）
# ============================================================================
cat("\n")
cat("=======================================================\n")
cat("  R0_vs_R1_tumor: KEGG - All Significant\n")
cat("=======================================================\n")

kegg_ret <- results$GSEA$R0_vs_R1_tumor$KEGG@result
kegg_ret_sig <- kegg_ret[kegg_ret$p.adjust < 0.05, ]

# UP
kegg_up <- kegg_ret_sig[kegg_ret_sig$NES > 0, ]
kegg_up <- kegg_up[order(-kegg_up$NES), ]
cat("\n--- UP in R1 (high radiation) ---\n")
if (nrow(kegg_up) > 0) {
  print(kegg_up[, c("Description", "NES", "p.adjust", "setSize")], row.names = FALSE)
} else {
  cat("  (None)\n")
}

# DOWN
kegg_down <- kegg_ret_sig[kegg_ret_sig$NES < 0, ]
kegg_down <- kegg_down[order(kegg_down$NES), ]
cat("\n--- DOWN in R1 (high radiation) ---\n")
if (nrow(kegg_down) > 0) {
  print(kegg_down[, c("Description", "NES", "p.adjust", "setSize")], row.names = FALSE)
} else {
  cat("  (None)\n")
}

# ============================================================================
# R0_vs_R1_tumor: Reactome（上位30 - UP/DOWN分離）
# ============================================================================
cat("\n")
cat("=======================================================\n")
cat("  R0_vs_R1_tumor: Reactome - Top by |NES|\n")
cat("=======================================================\n")

reactome_ret <- results$GSEA$R0_vs_R1_tumor$Reactome@result
reactome_ret_sig <- reactome_ret[reactome_ret$p.adjust < 0.05, ]

# UP
reactome_up <- reactome_ret_sig[reactome_ret_sig$NES > 0, ]
reactome_up <- reactome_up[order(-reactome_up$NES), ]
cat("\n--- UP in R1 (all", nrow(reactome_up), ") ---\n")
if (nrow(reactome_up) > 0) {
  print(reactome_up[, c("Description", "NES", "p.adjust", "setSize")], row.names = FALSE)
} else {
  cat("  (None)\n")
}

# DOWN
reactome_down <- reactome_ret_sig[reactome_ret_sig$NES < 0, ]
reactome_down <- reactome_down[order(reactome_down$NES), ]
cat("\n--- DOWN in R1 (top 30 of", nrow(reactome_down), ") ---\n")
if (nrow(reactome_down) > 0) {
  print(head(reactome_down[, c("Description", "NES", "p.adjust", "setSize")], 30), row.names = FALSE)
} else {
  cat("  (None)\n")
}

# ============================================================================
# R0_vs_R1_tumor: GO_BP（上位30 - UP/DOWN分離）
# ============================================================================
cat("\n")
cat("=======================================================\n")
cat("  R0_vs_R1_tumor: GO_BP - Top by |NES|\n")
cat("=======================================================\n")

gobp_ret <- results$GSEA$R0_vs_R1_tumor$GO_BP@result
gobp_ret_sig <- gobp_ret[gobp_ret$p.adjust < 0.05, ]

# UP
gobp_up <- gobp_ret_sig[gobp_ret_sig$NES > 0, ]
gobp_up <- gobp_up[order(-gobp_up$NES), ]
cat("\n--- UP in R1 (all", nrow(gobp_up), ") ---\n")
if (nrow(gobp_up) > 0) {
  print(gobp_up[, c("Description", "NES", "p.adjust", "setSize")], row.names = FALSE)
} else {
  cat("  (None)\n")
}

# DOWN
gobp_down <- gobp_ret_sig[gobp_ret_sig$NES < 0, ]
gobp_down <- gobp_down[order(gobp_down$NES), ]
cat("\n--- DOWN in R1 (top 30 of", nrow(gobp_down), ") ---\n")
if (nrow(gobp_down) > 0) {
  print(head(gobp_down[, c("Description", "NES", "p.adjust", "setSize")], 30), row.names = FALSE)
} else {
  cat("  (None)\n")
}

# ============================================================================
# DNA修復関連パスウェイの詳細
# ============================================================================
cat("\n")
cat("=======================================================\n")
cat("  DNA Repair Related Pathways (KEGG + Reactome)\n")
cat("=======================================================\n")

# KEGG DNA修復関連
dna_keywords <- c("repair", "recombination", "Fanconi", "mismatch", "excision", 
                  "damage", "checkpoint", "cycle", "replication")
kegg_all <- results$GSEA$R0_vs_R1_tumor$KEGG@result
kegg_dna <- kegg_all[grepl(paste(dna_keywords, collapse = "|"), 
                           kegg_all$Description, ignore.case = TRUE), ]
kegg_dna <- kegg_dna[order(kegg_dna$NES), ]

cat("\n--- KEGG: DNA repair / Cell cycle related ---\n")
if (nrow(kegg_dna) > 0) {
  print(kegg_dna[, c("Description", "NES", "p.adjust", "setSize")], row.names = FALSE)
}

# Reactome DNA修復関連（全件）
reactome_all <- results$GSEA$R0_vs_R1_tumor$Reactome@result
reactome_dna <- reactome_all[grepl(paste(dna_keywords, collapse = "|"), 
                                   reactome_all$Description, ignore.case = TRUE), ]
reactome_dna <- reactome_dna[order(reactome_dna$NES), ]

cat("\n--- Reactome: DNA repair / Cell cycle related (all", nrow(reactome_dna), ") ---\n")
if (nrow(reactome_dna) > 0) {
  print(reactome_dna[, c("Description", "NES", "p.adjust", "setSize")], row.names = FALSE)
}

# ============================================================================
# Cross-comparison詳細
# ============================================================================
cat("\n")
cat("=======================================================\n")
cat("  Cross-comparison: RET x BRAF (Hallmark)\n")
cat("=======================================================\n")

cc_file <- "analysis_v7/output/enrichment_analysis/cross_comparison/RET_vs_BRAF_Hallmark_consistent.csv"
if (file.exists(cc_file)) {
  cc_hallmark <- read.csv(cc_file)
  print(cc_hallmark[, c("Description", "NES_RET", "NES_BRAF", "avg_abs_NES")])
} else {
  cat("  File not found\n")
}

cat("\n")
cat("=======================================================\n")
cat("  Cross-comparison: Tumor x Normal (Hallmark)\n")
cat("=======================================================\n")

cc_tn_file <- "analysis_v7/output/enrichment_analysis/cross_comparison/Tumor_vs_Normal_Hallmark_consistent.csv"
if (file.exists(cc_tn_file)) {
  cc_tn_hallmark <- read.csv(cc_tn_file)
  print(cc_tn_hallmark[, c("Description", "NES_Tumor", "NES_Normal", "avg_abs_NES")])
} else {
  cat("  File not found\n")
}

# ============================================================================
# Leading edge genes（主要パスウェイ）- Gene Symbol表示
# ============================================================================
cat("\n")
cat("=======================================================\n")
cat("  Leading Edge Genes: Key Pathways (with Gene Symbols)\n")
cat("=======================================================\n")

# E2F_TARGETS
cat("\n--- E2F_TARGETS ---\n")
e2f_row <- hallmark_ret[hallmark_ret$Description == "E2F_TARGETS", ]
if (nrow(e2f_row) > 0) {
  cat("NES:", round(e2f_row$NES, 3), ", q-value:", format(e2f_row$p.adjust, digits = 3), "\n")
  le_genes <- strsplit(e2f_row$core_enrichment, "/")[[1]]
  le_symbols <- entrez_to_symbol(le_genes)
  cat("Leading edge genes (", length(le_genes), "):\n", sep = "")
  cat(paste(head(le_symbols, 30), collapse = ", "), "\n")
  if (length(le_genes) > 30) cat("  ... and", length(le_genes) - 30, "more\n")
}

# G2M_CHECKPOINT
cat("\n--- G2M_CHECKPOINT ---\n")
g2m_row <- hallmark_ret[hallmark_ret$Description == "G2M_CHECKPOINT", ]
if (nrow(g2m_row) > 0) {
  cat("NES:", round(g2m_row$NES, 3), ", q-value:", format(g2m_row$p.adjust, digits = 3), "\n")
  le_genes <- strsplit(g2m_row$core_enrichment, "/")[[1]]
  le_symbols <- entrez_to_symbol(le_genes)
  cat("Leading edge genes (", length(le_genes), "):\n", sep = "")
  cat(paste(head(le_symbols, 30), collapse = ", "), "\n")
  if (length(le_genes) > 30) cat("  ... and", length(le_genes) - 30, "more\n")
}

# Mismatch repair (KEGG)
cat("\n--- Mismatch repair (KEGG) ---\n")
mmr_row <- kegg_ret[grepl("Mismatch repair", kegg_ret$Description), ]
if (nrow(mmr_row) > 0) {
  cat("NES:", round(mmr_row$NES[1], 3), ", q-value:", format(mmr_row$p.adjust[1], digits = 3), "\n")
  le_genes <- strsplit(mmr_row$core_enrichment[1], "/")[[1]]
  le_symbols <- entrez_to_symbol(le_genes)
  cat("Leading edge genes (", length(le_genes), "): ", sep = "")
  cat(paste(le_symbols, collapse = ", "), "\n")
}

# Homologous recombination (KEGG)
cat("\n--- Homologous recombination (KEGG) ---\n")
hr_row <- kegg_ret[grepl("Homologous recombination", kegg_ret$Description), ]
if (nrow(hr_row) > 0) {
  cat("NES:", round(hr_row$NES[1], 3), ", q-value:", format(hr_row$p.adjust[1], digits = 3), "\n")
  le_genes <- strsplit(hr_row$core_enrichment[1], "/")[[1]]
  le_symbols <- entrez_to_symbol(le_genes)
  cat("Leading edge genes (", length(le_genes), "):\n", sep = "")
  cat(paste(le_symbols, collapse = ", "), "\n")
}

# Fanconi anemia pathway (KEGG)
cat("\n--- Fanconi anemia pathway (KEGG) ---\n")
fa_row <- kegg_ret[grepl("Fanconi", kegg_ret$Description), ]
if (nrow(fa_row) > 0) {
  cat("NES:", round(fa_row$NES[1], 3), ", q-value:", format(fa_row$p.adjust[1], digits = 3), "\n")
  le_genes <- strsplit(fa_row$core_enrichment[1], "/")[[1]]
  le_symbols <- entrez_to_symbol(le_genes)
  cat("Leading edge genes (", length(le_genes), "):\n", sep = "")
  cat(paste(le_symbols, collapse = ", "), "\n")
}

# Cell Cycle Checkpoints (Reactome)
cat("\n--- Cell Cycle Checkpoints (Reactome) ---\n")
ccc_row <- reactome_ret[grepl("^Cell Cycle Checkpoints$", reactome_ret$Description), ]
if (nrow(ccc_row) > 0) {
  cat("NES:", round(ccc_row$NES[1], 3), ", q-value:", format(ccc_row$p.adjust[1], digits = 3), "\n")
  le_genes <- strsplit(ccc_row$core_enrichment[1], "/")[[1]]
  le_symbols <- entrez_to_symbol(le_genes)
  cat("Leading edge genes (", length(le_genes), "):\n", sep = "")
  cat(paste(head(le_symbols, 40), collapse = ", "), "\n")
  if (length(le_genes) > 40) cat("  ... and", length(le_genes) - 40, "more\n")
}

# HDR through HRR (Reactome)
cat("\n--- HDR through Homologous Recombination (Reactome) ---\n")
hdr_row <- reactome_ret[grepl("HDR through Homologous Recombination \\(HRR\\)$", reactome_ret$Description), ]
if (nrow(hdr_row) > 0) {
  cat("NES:", round(hdr_row$NES[1], 3), ", q-value:", format(hdr_row$p.adjust[1], digits = 3), "\n")
  le_genes <- strsplit(hdr_row$core_enrichment[1], "/")[[1]]
  le_symbols <- entrez_to_symbol(le_genes)
  cat("Leading edge genes (", length(le_genes), "):\n", sep = "")
  cat(paste(le_symbols, collapse = ", "), "\n")
}

# DNA Double-Strand Break Repair (Reactome)
cat("\n--- DNA Double-Strand Break Repair (Reactome) ---\n")
dsb_row <- reactome_ret[grepl("^DNA Double-Strand Break Repair$", reactome_ret$Description), ]
if (nrow(dsb_row) > 0) {
  cat("NES:", round(dsb_row$NES[1], 3), ", q-value:", format(dsb_row$p.adjust[1], digits = 3), "\n")
  le_genes <- strsplit(dsb_row$core_enrichment[1], "/")[[1]]
  le_symbols <- entrez_to_symbol(le_genes)
  cat("Leading edge genes (", length(le_genes), "):\n", sep = "")
  cat(paste(head(le_symbols, 40), collapse = ", "), "\n")
  if (length(le_genes) > 40) cat("  ... and", length(le_genes) - 40, "more\n")
}

# ============================================================================
# B0_vs_B1_tumor との比較
# ============================================================================
cat("\n")
cat("=======================================================\n")
cat("  Comparison: B0_vs_B1_tumor (BRAF background)\n")
cat("=======================================================\n")

cat("\n--- Hallmark ---\n")
hallmark_braf <- results$GSEA$B0_vs_B1_tumor$Hallmark@result
hallmark_braf_sig <- hallmark_braf[hallmark_braf$p.adjust < 0.05, ]
hallmark_braf_sig <- hallmark_braf_sig[order(hallmark_braf_sig$NES), ]
if (nrow(hallmark_braf_sig) > 0) {
  print(hallmark_braf_sig[, c("Description", "NES", "p.adjust", "setSize")], row.names = FALSE)
} else {
  cat("  (None significant)\n")
}

# NES方向集計
cat("\n--- NES Direction Summary (B0_vs_B1_tumor) ---\n")
for (db in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
  res <- results$GSEA$B0_vs_B1_tumor[[db]]@result
  sig <- res[res$p.adjust < 0.05, ]
  up <- sum(sig$NES > 0)
  down <- sum(sig$NES < 0)
  cat(sprintf("  %-12s: UP(B1 high)=%3d, DOWN(B1 low)=%3d (total=%d)\n", 
              db, up, down, up + down))
}

# ============================================================================
# R1でUPしているパスウェイの詳細（Leading edge）
# ============================================================================
cat("\n")
cat("=======================================================\n")
cat("  Leading Edge: Pathways UP in R1 (high radiation)\n")
cat("=======================================================\n")

# MYOGENESIS (Hallmark)
cat("\n--- MYOGENESIS (Hallmark) ---\n")
myo_row <- hallmark_ret[hallmark_ret$Description == "MYOGENESIS", ]
if (nrow(myo_row) > 0 && myo_row$p.adjust < 0.05) {
  cat("NES:", round(myo_row$NES, 3), ", q-value:", format(myo_row$p.adjust, digits = 3), "\n")
  le_genes <- strsplit(myo_row$core_enrichment, "/")[[1]]
  le_symbols <- entrez_to_symbol(le_genes)
  cat("Leading edge genes (", length(le_genes), "):\n", sep = "")
  cat(paste(head(le_symbols, 30), collapse = ", "), "\n")
  if (length(le_genes) > 30) cat("  ... and", length(le_genes) - 30, "more\n")
}

# Protein processing in ER (KEGG)
cat("\n--- Protein processing in ER (KEGG) ---\n")
er_row <- kegg_ret[grepl("Protein processing in endoplasmic reticulum", kegg_ret$Description), ]
if (nrow(er_row) > 0 && er_row$NES[1] > 0) {
  cat("NES:", round(er_row$NES[1], 3), ", q-value:", format(er_row$p.adjust[1], digits = 3), "\n")
  le_genes <- strsplit(er_row$core_enrichment[1], "/")[[1]]
  le_symbols <- entrez_to_symbol(le_genes)
  cat("Leading edge genes (", length(le_genes), "):\n", sep = "")
  cat(paste(head(le_symbols, 30), collapse = ", "), "\n")
  if (length(le_genes) > 30) cat("  ... and", length(le_genes) - 30, "more\n")
}

# AKT activation (GO_BP)
cat("\n--- activation of protein kinase B activity (GO_BP) ---\n")
akt_row <- gobp_ret[grepl("activation of protein kinase B activity", gobp_ret$Description), ]
if (nrow(akt_row) > 0 && akt_row$NES[1] > 0) {
  cat("NES:", round(akt_row$NES[1], 3), ", q-value:", format(akt_row$p.adjust[1], digits = 3), "\n")
  le_genes <- strsplit(akt_row$core_enrichment[1], "/")[[1]]
  le_symbols <- entrez_to_symbol(le_genes)
  cat("Leading edge genes (", length(le_genes), "): ", sep = "")
  cat(paste(le_symbols, collapse = ", "), "\n")
}

# ============================================================================
# 重要遺伝子のサマリー
# ============================================================================
cat("\n")
cat("=======================================================\n")
cat("  Key DNA Repair Genes Summary\n")
cat("=======================================================\n")

# 複数のパスウェイで共通して出現する遺伝子を特定
cat("\n--- Genes appearing in multiple DNA repair pathways ---\n")

# 各パスウェイからleading edge遺伝子を収集
all_dna_repair_genes <- list()

if (nrow(mmr_row) > 0) {
  all_dna_repair_genes[["MMR"]] <- strsplit(mmr_row$core_enrichment[1], "/")[[1]]
}
if (nrow(hr_row) > 0) {
  all_dna_repair_genes[["HR"]] <- strsplit(hr_row$core_enrichment[1], "/")[[1]]
}
if (nrow(fa_row) > 0) {
  all_dna_repair_genes[["FA"]] <- strsplit(fa_row$core_enrichment[1], "/")[[1]]
}

if (length(all_dna_repair_genes) > 1) {
  all_genes <- unlist(all_dna_repair_genes)
  gene_counts <- table(all_genes)
  multi_genes <- names(gene_counts[gene_counts > 1])
  
  if (length(multi_genes) > 0) {
    multi_symbols <- entrez_to_symbol(multi_genes)
    cat("Genes in 2+ pathways:\n")
    for (i in seq_along(multi_genes)) {
      pathways <- names(all_dna_repair_genes)[sapply(all_dna_repair_genes, function(x) multi_genes[i] %in% x)]
      cat(sprintf("  %s: %s\n", multi_symbols[i], paste(pathways, collapse = ", ")))
    }
  }
}

cat("\n=== スクリプト実行完了 ===\n")