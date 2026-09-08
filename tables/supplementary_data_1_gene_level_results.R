# supplementary_data_1_gene_level_results.R  (Supplementary Data 1)
# 遺伝子レベル全結果の完全開示。
#
# 整形のみ・計算なし。確定済み 410 出力(processed/thyr_expression_test.rds)から、
# 全検定遺伝子 × 4対比の遺伝子別統計量(効果 θ・exact 置換 p・Storey q)を一方向で
# CSV 化する。セットレベルの Supplementary Data 2(supplementary_data_2_set_level_results.R)の遺伝子レベル
# 対応物 — 単一閾値に依存させない開示方針(R(α) 保持・非2値化)の実体。
# 照合対象: 対比別の行数 = N-15(15,621 / 15,466 / 16,002 / 15,863)、
# q<0.10 員数 = N-16(1,765 / 0 / 0 / 1)、B_Normal の1件 = BHLHB9(N-22)。
# Input : processed/thyr_expression_test.rds (from 410), processed/thyr_se_raw.rds
# Output: output/tables/supplementary_data_1.csv
#

source("setup.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
})
source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "annotation.R"))

test <- readRDS(file.path(paths$processed, "thyr_expression_test.rds"))
se <- readRDS(file.path(paths$processed, "thyr_se_raw.rds"))
name_of <- gene_name_map(se)

out <- do.call(rbind, lapply(UNIT_ORDER, function(u) {
  g <- test$units[[u]]$genes
  data.frame(
    contrast = u,
    gene_id = g$gene_id,
    gene_symbol = unname(name_of[g$gene_id]),
    effect = g$effect,
    p_exact = g$p_exact,
    q_storey = g$q_storey,
    stringsAsFactors = FALSE
  )
}))
out <- out[order(match(out$contrast, UNIT_ORDER), out$p_exact, out$gene_id), ]

dir.create(file.path(paths$output, "tables"), recursive = TRUE, showWarnings = FALSE)
path <- file.path(paths$output, "tables", "supplementary_data_1.csv")
write.csv(out, path, row.names = FALSE)

for (u in UNIT_ORDER) {
  d <- out[out$contrast == u, ]
  message(sprintf("  %-9s genes %5d | q<%.2f %4d | min p %.3g",
    u, nrow(d), FDR_CUT, sum(d$q_storey < FDR_CUT), min(d$p_exact)))
}
bh <- out[out$contrast == "B_Normal" & out$q_storey < FDR_CUT, ]
message("  B_Normal hit: ", paste(bh$gene_id, bh$gene_symbol, collapse = " "))
cat("Saved:", path, "rows:", nrow(out), "\n")
