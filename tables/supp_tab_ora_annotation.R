# supp_tab_ora_annotation.R — Table S3(R_Tumor ORA 注釈の全表)の整形出力
#
# 整形のみ・計算なし。確定済み 430 出力(processed/thyr_deg_ora_annotation.rds の $table、
# 18,576 行 = 6,192 セット × up/down/combined)を一方向で CSV 化する。
# 照合対象: family×list の q_bh<0.10 員数が N-59(up 0/0/0/0、down 12/105/205/7、
# combined 6/37/71/2)と一致すること(初回実行 2026-08-21 で全一致を確認済み)。
#
# Input:  processed/thyr_deg_ora_annotation.rds
# Output: output/tables/supp_tab_ora_annotation.csv

ora <- readRDS(file.path("processed", "thyr_deg_ora_annotation.rds"))
tb <- ora$table
tb <- tb[order(tb$family, tb$list, tb$p_hyper),
         c("family", "list", "pathway", "set_size", "overlap", "expected", "p_hyper", "q_bh")]
tb$expected <- round(tb$expected, 4)

dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
out <- file.path("output", "tables", "supp_tab_ora_annotation.csv")
write.csv(tb, out, row.names = FALSE)

n_hit <- with(tb[tb$q_bh < 0.10, ], table(family, list))
cat("rows:", nrow(tb), "\n")
print(n_hit)
cat("Saved:", out, "\n")
