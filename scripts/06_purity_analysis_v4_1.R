# ==============================================================================
# 06_purity_analysis_v4_1.R  — contamDE.lm を用いた腫瘍相対純度QC（全群）
# ==============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(contamDE)
  library(matrixStats)
  library(dplyr)
})

cat("\n>>> [P0] start QC v4.1 (contamDE.lm) <<<\n")

## 設定（最初にまとめる）
PCA_FILTER_FILE    <- "./data/processed/pca_filtered_sample_lists_v4.rda"  # ←①のファイル
ASSAY_NAME         <- "stranded_second"
MIN_TUMOR_PURITY   <- 0.60
FILTER_BY_EXPR     <- TRUE
GROUPS             <- c("R0","R1","B0","B1")
OUT_TAG            <- "min060_pca"
MIN_SAMPLES_PER_ARM<- 3

## 入力読み込み
stopifnot(file.exists(PCA_FILTER_FILE))
load(PCA_FILTER_FILE)                                # => pca_filtered_sample_lists が入ってる
base_sample_lists <- pca_filtered_sample_lists       # ←ここを明示
cat(sprintf("[P0] using PCA file: %s\n", PCA_FILTER_FILE))

## カウント等
stopifnot(exists("se_thyr"))
counts_full <- assay(se_thyr, ASSAY_NAME)
gene_info   <- rowData(se_thyr)

## protein-coding（NAセーフ）
pc_mask <- if ("gene_type" %in% colnames(gene_info)) {
  !is.na(gene_info$gene_type) & gene_info$gene_type == "protein_coding"
} else {
  rep(TRUE, nrow(counts_full))
}
cat(sprintf("[P0] assay=%s | genes=%d (pc=%d) | samples=%d\n",
            ASSAY_NAME, nrow(counts_full), sum(pc_mask), ncol(counts_full)))

# ---------------------------
# 1) 1グループ処理関数（lm 専用・QCのみ）
# ---------------------------
run_contamde_lm_one <- function(g) {
  cat(sprintf("\n>>> [G:%s] start <<<\n", g))
  
  # ペア取得
  tum_ids <- intersect(base_sample_lists[[g]]$tumor,  colnames(counts_full))
  nor_ids <- intersect(base_sample_lists[[g]]$normal, colnames(counts_full))
  N <- min(length(tum_ids), length(nor_ids))
  cat(sprintf("[G:%s] pairs=%d (tum=%d, nor=%d)\n", g, N, length(tum_ids), length(nor_ids)))
  if (N < MIN_SAMPLES_PER_ARM) {
    return(list(ok=FALSE, reason=sprintf("too few pairs (%d)", N)))
  }
  
  # 並び: normal→tumor（(i, N+i)が同一ペア）
  cnt_nor <- counts_full[pc_mask, nor_ids, drop=FALSE][, seq_len(N), drop=FALSE]
  cnt_tum <- counts_full[pc_mask, tum_ids, drop=FALSE][, seq_len(N), drop=FALSE]
  counts_lm <- cbind(cnt_nor, cnt_tum)
  
  # edgeR 前処理（TMM）
  dge <- DGEList(counts_lm,
                 group=factor(c(rep("normal", N), rep("tumor", N)), levels=c("normal","tumor")))
  if (FILTER_BY_EXPR) {
    keep <- filterByExpr(dge, group=dge$samples$group)
    dge  <- dge[keep, , keep.lib.sizes=FALSE]
  }
  dge <- calcNormFactors(dge, "TMM")
  cat(sprintf("[G:%s] dim after TMM/filter=%d x %d\n", g, nrow(dge), ncol(dge)))
  
  # 低変動（腫瘍-正常の log2差が全くゼロ）を除外して mgcv::gam の非有限を防止
  Np <- ncol(dge)/2
  cnt_n <- dge$counts[, 1:Np, drop=FALSE]
  cnt_t <- dge$counts[, (Np+1):(2*Np), drop=FALSE]
  y_pre <- log2(cnt_t + 1) - log2(cnt_n + 1)
  v_pre <- rowVars(as.matrix(y_pre))
  keepV <- is.finite(v_pre) & v_pre > 0
  dge2  <- dge[keepV, , keep.lib.sizes=FALSE]
  cat(sprintf("[G:%s] kept genes (var>0): %d (%.1f%%)\n", g, sum(keepV), 100*sum(keepV)/length(keepV)))
  
  # contamDE.lm 本体（DE出力は使わず proportion のみ活用）
  res <- contamDE::contamDE.lm(counts=dge2$counts,
                               subtype=NULL, covariate=NULL,
                               is.contaminated=TRUE, robust=TRUE)
  
  prop <- as.numeric(res$proportion)  # 相対純度（max=1）: 長さは N
  if (length(prop) != Np) {
    return(list(ok=FALSE, reason=sprintf("proportion length mismatch: %d vs %d", length(prop), Np)))
  }
  tumor_names <- colnames(dge2)[(Np+1):(2*Np)]
  names(prop) <- tumor_names
  
  cat(sprintf("[G:%s] purity stats: median=%.3f, IQR=[%.3f, %.3f], max=%.3f\n",
              g, median(prop), quantile(prop,.25), quantile(prop,.75), max(prop)))
  
  # しきい値適用（腫瘍のみで判定）、pair consistency
  low_tum <- names(prop)[prop < MIN_TUMOR_PURITY]
  orig_cases  <- base_sample_lists[[g]]$cases
  orig_tumor  <- base_sample_lists[[g]]$tumor
  orig_normal <- base_sample_lists[[g]]$normal
  drop_idx <- which(orig_tumor %in% low_tum)
  keep_idx <- setdiff(seq_along(orig_cases), drop_idx)
  
  out_lists <- list(
    tumor  = orig_tumor[keep_idx],
    normal = orig_normal[keep_idx],
    cases  = orig_cases[keep_idx]
  )
  cat(sprintf("[G:%s] dropped %d pairs @thr=%.2f (max=1 scale) | remain=%d\n",
              g, length(drop_idx), MIN_TUMOR_PURITY, length(keep_idx)))
  
  list(
    ok=TRUE,
    group=g,
    N_pairs=N,
    purity_tumor=prop,
    dropped_pairs=length(drop_idx),
    remain_pairs=length(keep_idx),
    filtered_lists=out_lists
  )
}

# ---------------------------
# 2) 全群ループ
# ---------------------------
purity_results <- list()
final_high_purity_base_sample_lists <- base_sample_lists  # 上書き

for (g in GROUPS) {
  if (!g %in% names(base_sample_lists)) next
  rr <- run_contamde_lm_one(g)
  if (!isTRUE(rr$ok)) {
    cat(sprintf("[G:%s] SKIP: %s\n", g, rr$reason))
    next
  }
  purity_results[[g]] <- rr
  final_high_purity_base_sample_lists[[g]]$tumor  <- rr$filtered_lists$tumor
  final_high_purity_base_sample_lists[[g]]$normal <- rr$filtered_lists$normal
  final_high_purity_base_sample_lists[[g]]$cases  <- rr$filtered_lists$cases
}

# ---------------------------
# 3) サマリ出力
# ---------------------------
purity_summary <- bind_rows(lapply(names(purity_results), function(g){
  rr <- purity_results[[g]]
  data.frame(
    Group = g,
    Tumor_N_before   = rr$N_pairs,
    Tumor_low_purity = rr$dropped_pairs,
    Pairs_after      = rr$remain_pairs,
    Tumor_purity_median = median(rr$purity_tumor, na.rm=TRUE),
    stringsAsFactors = FALSE
  )
}))

print(purity_summary)

# ---------------------------
# 4) 保存（閾値タグ付き）
# ---------------------------
if (!dir.exists("./data/processed")) dir.create("./data/processed", recursive=TRUE)

obj_name <- paste0("final_high_purity_sample_lists_v4_main_", OUT_TAG)
assign(obj_name, final_high_purity_base_sample_lists)   # ← base を保存
save(list=obj_name, file=paste0("./data/processed/", obj_name, ".rda"))

write.csv(purity_summary,
          file=paste0("./data/processed/purity_summary_v4_1_", OUT_TAG, ".csv"),
          row.names=FALSE)

cat("\n>>> [P0] done | saved: ",
    paste0(obj_name, ".rda, purity_summary_v4_1_", OUT_TAG, ".csv"),
    " <<<\n")

