# ==============================================================================
# REBC-THYR Purity Analysis v4 - contamDE-based purity after PCA v4
# 06_purity_analysis_v4.R
# ==============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(contamDE)   # CRAN/GitHub: 'contamDE'
  library(matrixStats)
  library(dplyr)
})

cat("Starting purity analysis v4 (contamDE after PCA v4)...\n")

# ------------------------------------------------------------------------------
# 0) Config
# ------------------------------------------------------------------------------
USE_SENS <- FALSE        # 感度ライン（B1_normal復帰）を使うなら TRUE
ASSAY_NAME <- "stranded_second"

# contamDE しきい値（必要に応じて調整）
MIN_TUMOR_PURITY   <- 0.60   # 腫瘍サンプル: これ未満は低純度として除外
MAX_NORMAL_TUMORLIKE <- 0.20 # 正常サンプル: 腫瘍成分がこの閾値を越えたら混入疑いとして除外

# フィルタと前処理
FILTER_BY_EXPR <- TRUE
MIN_SAMPLES_PER_ARM <- 3     # contamDEを走らせる最低サンプル数/群
PRIOR_COUNT <- 0.5           # logCPMの擬似カウント（チェック用可視化時）

# 出力タグ
qc_tag <- if (USE_SENS) "v4_sens" else "v4_main"

# ------------------------------------------------------------------------------
# 1) Inputs: sample lists (from PCA v4) & counts
# ------------------------------------------------------------------------------
cat("Loading PCA-filtered sample lists...\n")
if (USE_SENS) {
  load("./data/processed/pca_filtered_sample_lists_v4_sens.rda")   # -> pca_filtered_sample_lists_sens
  sample_lists <- pca_filtered_sample_lists_sens
} else {
  load("./data/processed/pca_filtered_sample_lists_v4.rda")        # -> pca_filtered_sample_lists
  sample_lists <- pca_filtered_sample_lists
}

cat("Loading SummarizedExperiment counts...\n")
if (!exists("se_thyr")) {
  load("./data/raw/thyr_data.rda")   # -> se_thyr もしくは data
  if (exists("data")) { se_thyr <- data; rm(data) }
}
counts_full <- assay(se_thyr, ASSAY_NAME)
gene_info   <- as.data.frame(rowData(se_thyr))

cat("Counts dims: ", paste(dim(counts_full), collapse=" x "), "\n")

# biotype列の候補
bio_col <- {
  cands <- c("biotype","gene_type","type")
  hit <- cands[cands %in% colnames(gene_info)]
  if (length(hit)) hit[1] else NA_character_
}
if (!is.na(bio_col)) {
  pc_mask <- gene_info[[bio_col]] %in% c("protein_coding","protein-coding","protein coding")
} else {
  warning("biotype/gene_type列が見つからないため、protein-codingフィルタはスキップします")
  pc_mask <- rep(TRUE, nrow(gene_info))
}

# ------------------------------------------------------------------------------
# 2) Helper: run contamDE for a given (group, tissue)  【API自動検出版】
# ------------------------------------------------------------------------------
run_contamde_one <- function(group_name, tissue, samples_vec) {
  smp <- samples_vec[samples_vec %in% colnames(counts_full)]
  if (length(smp) < MIN_SAMPLES_PER_ARM) {
    return(list(ok = FALSE, reason = sprintf("too few samples (%d)", length(smp))))
  }
  
  # 同一 group 内の tumor / normal
  tum_ids <- intersect(sample_lists[[group_name]]$tumor, colnames(counts_full))
  nor_ids <- intersect(sample_lists[[group_name]]$normal, colnames(counts_full))
  if (length(tum_ids) < MIN_SAMPLES_PER_ARM || length(nor_ids) < MIN_SAMPLES_PER_ARM) {
    return(list(ok = FALSE, reason = "paired arms too small for contamDE"))
  }
  
  # protein-coding 抽出
  cnt_all <- cbind(counts_full[pc_mask, tum_ids, drop=FALSE],
                   counts_full[pc_mask, nor_ids, drop=FALSE])
  grp_all <- factor(c(rep("tumor", length(tum_ids)), rep("normal", length(nor_ids))),
                    levels = c("normal","tumor"))
  
  # edgeR 前処理
  dge <- edgeR::DGEList(counts = cnt_all, group = grp_all)
  if (FILTER_BY_EXPR) {
    keep <- edgeR::filterByExpr(dge, group = dge$samples$group)
    dge  <- dge[keep, , keep.lib.sizes = FALSE]
  }
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  
  # ---- contamDE 呼び出し（API自動検出）----
  run_cd <- function(fun) {
    fml <- names(formals(fun))
    args <- list()
    
    # カウント引数の候補
    if ("counts" %in% fml)       args$counts       <- dge$counts
    if ("X" %in% fml)            args$X            <- dge$counts
    if ("data" %in% fml)         args$data         <- dge$counts
    
    # グループ/ラベルの候補（0/1 or factor どちらも許容）
    g01 <- as.integer(dge$samples$group == "tumor")
    if ("group" %in% fml)        args$group        <- g01
    if ("Y" %in% fml)            args$Y            <- g01
    if ("label" %in% fml)        args$label        <- g01
    if ("Group" %in% fml)        args$Group        <- dge$samples$group
    
    # 正規化係数（もし受け取る実装なら渡す）
    if ("norm_factors" %in% fml) args$norm_factors <- dge$samples$norm.factors
    if ("norm.factors" %in% fml) args$norm.factors <- dge$samples$norm.factors
    
    # 冗長引数は渡さない（unused回避）
    do.call(fun, args)
  }
  
  res <- tryCatch({
    # まず contamDE()
    run_cd(contamDE::contamDE)
  }, error = function(e1) {
    tryCatch({
      # ダメなら contamDE.lm() も試す
      run_cd(contamDE::contamDE.lm)
    }, error = function(e2) {
      list(`.error` = paste("contamDE unused-args:", e1$message,
                            "| contamDE.lm unused-args:", e2$message))
    })
  })
  
  if (!is.null(res$.error)) {
    return(list(ok = FALSE, reason = res$.error))
  }
  
  # ---- 戻り値の正規化（purity / prop / 他）----
  # よくあるフィールド名から候補を拾う
  cand <- NULL
  for (nm in c("purity","Purity","pureness","prob_tumor","prop","proportion","contam","contamination")) {
    if (!is.null(res[[nm]])) { cand <- res[[nm]]; break }
  }
  # 汚染率なら純度に反転
  if (!is.null(cand) && any(names(res) %in% c("prop","proportion","contam","contamination"))) {
    cand <- 1 - cand
  }
  
  # インデックス
  # ---- 戻り値の正規化（prop/purity の“向き”と“対応”を自動判定）----
  idx_tum <- which(dge$samples$group == "tumor")
  idx_nor <- which(dge$samples$group == "normal")
  n_all <- ncol(dge); n_tum <- length(idx_tum); n_nor <- length(idx_nor)
  
  # 候補フィールドから拾う
  prop <- NULL; pur  <- NULL
  for (nm in c("prop","proportion","contam","contamination","alpha","pi")) if (!is.null(res[[nm]])) { prop <- as.numeric(res[[nm]]); break }
  for (nm in c("purity","Purity","pureness","prob_tumor")) if (!is.null(res[[nm]])) { pur <- as.numeric(res[[nm]]); break }
  
  # contamDE は“腫瘍だけ”返す実装が多い想定：
  # 1) purity がそのまま腫瘍長で返っている
  if (!is.null(pur) && length(pur) == n_tum) {
    purity_tumor <- pur
    names(purity_tumor) <- colnames(dge)[idx_tum]
    tumorlike_normal <- setNames(rep(0, n_nor), colnames(dge)[idx_nor])
    
    # 2) proportion（汚染率 or 腫瘍率）が腫瘍長で返っている
  } else if (!is.null(prop) && length(prop) == n_tum) {
    # 向きを自動決定： 1-prop と prop のどちらが“落ち過ぎないか”
    cand1 <- 1 - prop    # 仮: 汚染率→純度
    cand2 <- prop        # 仮: これ自体が純度
    drop1 <- sum(cand1 < MIN_TUMOR_PURITY, na.rm = TRUE)
    drop2 <- sum(cand2 < MIN_TUMOR_PURITY, na.rm = TRUE)
    
    # ドロップが少ない方を純度と解釈（極端な全落ちは避ける）
    if (drop2 < drop1) {
      purity_tumor <- cand2
      chosen <- "prop(as purity)"
    } else {
      purity_tumor <- cand1
      chosen <- "1-prop(as purity)"
    }
    names(purity_tumor) <- colnames(dge)[idx_tum]
    tumorlike_normal <- setNames(rep(0, n_nor), colnames(dge)[idx_nor])
    
    cat(sprintf("  [auto] purity orientation: %s | median=%.2f | drop=%d/%d\n",
                chosen, median(purity_tumor, na.rm=TRUE), sum(purity_tumor < MIN_TUMOR_PURITY), n_tum))
    
    # 3) 名前付きで返ってくる（順序不明）→ 列名でマッピング
  } else if (!is.null(prop) && length(names(prop)) > 0 && any(names(prop) %in% colnames(dge))) {
    purity_tumor <- setNames(rep(NA_real_, n_tum), colnames(dge)[idx_tum])
    # tumor/normal どちらに対応するかを重なりで決める
    over_t <- sum(names(prop) %in% colnames(dge)[idx_tum])
    over_n <- sum(names(prop) %in% colnames(dge)[idx_nor])
    
    vec <- prop
    # 向き判定（prop→pure or 1-prop）を重なり側で行う
    if (over_t >= over_n) {
      # 腫瘍名が多く含まれる → 腫瘍へマップ
      hit <- intersect(names(vec), colnames(dge)[idx_tum])
      # 両向きを試す（腫瘍側の値のみで判定）
      vec1 <- rep(NA_real_, n_tum); names(vec1) <- colnames(dge)[idx_tum]; vec1[hit] <- 1 - vec[hit]
      vec2 <- rep(NA_real_, n_tum); names(vec2) <- colnames(dge)[idx_tum]; vec2[hit] <- vec[hit]
      drop1 <- sum(vec1 < MIN_TUMOR_PURITY, na.rm = TRUE)
      drop2 <- sum(vec2 < MIN_TUMOR_PURITY, na.rm = TRUE)
      purity_tumor <- if (drop2 < drop1) vec2 else vec1
    } else {
      # normal 名が多い → normal 側の腫瘍様を返している可能性。腫瘍側は proxy で補う
      warning("prop appears to map to NORMAL samples; filling tumor purity via proxy.")
      logcpm <- edgeR::cpm(dge, log=TRUE, prior.count=PRIOR_COUNT, normalized.lib.sizes=TRUE)
      pc1 <- stats::prcomp(t(scale(t(logcpm))))$x[,1]
      r <- rank(pc1[idx_tum], ties.method="average")
      purity_tumor <- (r - min(r)) / (max(r) - min(r) + 1e-8)
      names(purity_tumor) <- colnames(dge)[idx_tum]
    }
    tumorlike_normal <- setNames(rep(0, n_nor), colnames(dge)[idx_nor])
    
    # 4) それ以外 → フォールバック（PCA proxy; 腫瘍のみ）
  } else {
    warning("contamDE returned no usable purity/prop; using PCA proxy for tumor only.")
    logcpm <- edgeR::cpm(dge, log=TRUE, prior.count=PRIOR_COUNT, normalized.lib.sizes=TRUE)
    pc1 <- stats::prcomp(t(scale(t(logcpm))))$x[,1]
    r <- rank(pc1[idx_tum], ties.method="average")
    purity_tumor <- (r - min(r)) / (max(r) - min(r) + 1e-8)
    names(purity_tumor) <- colnames(dge)[idx_tum]
    tumorlike_normal <- setNames(rep(0, n_nor), colnames(dge)[idx_nor])
  }
  
  # 最終スプリット
  list(
    ok = TRUE,
    group = group_name,
    tumor_ids = colnames(dge)[idx_tum],
    normal_ids = colnames(dge)[idx_nor],
    purity_tumor = purity_tumor,
    tumorlike_normal = tumorlike_normal,
    raw = res
  )
  
  
}


# ------------------------------------------------------------------------------
# 3) Run contamDE per group (paired tumor/normal within group)
# ------------------------------------------------------------------------------
groups <- intersect(names(sample_lists), c("R0","R1","B0","B1"))
purity_results <- list()
for (g in groups) {
  cat(sprintf("\n--- Running contamDE for %s (paired tumor/normal) ---\n", g))
  s_any <- unique(c(sample_lists[[g]]$tumor, sample_lists[[g]]$normal))
  rr <- run_contamde_one(g, "paired", s_any)
  if (!isTRUE(rr$ok)) {
    cat("  Skipped: ", rr$reason, "\n", sep="")
    next
  }
  purity_results[[g]] <- rr
  cat(sprintf("  Tumor samples: %d | Normal samples: %d\n",
              length(rr$tumor_ids), length(rr$normal_ids)))
}

if (!length(purity_results)) {
  stop("No groups successfully processed by contamDE. Check inputs.")
}

# ------------------------------------------------------------------------------
# 4) Thresholding: build final_high_purity_sample_lists
# ------------------------------------------------------------------------------
final_high_purity_sample_lists <- sample_lists  # ベースはPCA後のものを引き継ぐ

purity_summary_rows <- list()

for (g in names(purity_results)) {
  rr <- purity_results[[g]]
  
  # 腫瘍の低純度だけで除外を決定（正常は表示用のみ）
  low_tum <- names(rr$purity_tumor)[which(rr$purity_tumor < MIN_TUMOR_PURITY)]
  
  orig_cases  <- sample_lists[[g]]$cases
  orig_tumor  <- sample_lists[[g]]$tumor
  orig_normal <- sample_lists[[g]]$normal
  
  drop_idx <- which(orig_tumor %in% low_tum)
  
  if (length(drop_idx)) {
    keep <- setdiff(seq_along(orig_cases), drop_idx)
    final_high_purity_sample_lists[[g]]$tumor  <- orig_tumor[keep]
    final_high_purity_sample_lists[[g]]$normal <- orig_normal[keep]
    final_high_purity_sample_lists[[g]]$cases  <- orig_cases[keep]
    cat(sprintf("%s: dropped %d pairs by tumor purity threshold\n", g, length(drop_idx)))
  } else {
    cat(sprintf("%s: no pairs dropped by tumor purity threshold\n", g))
  }
  
  # サマリは参考として normal 側の “tumor-like=0” を表示
  purity_summary_rows[[g]] <- data.frame(
    Group = g,
    Tumor_N_before  = length(orig_cases),
    Tumor_low_purity = length(low_tum),
    Pairs_after     = length(final_high_purity_sample_lists[[g]]$cases),
    Tumor_purity_median = median(rr$purity_tumor, na.rm=TRUE),
    stringsAsFactors = FALSE
  )
  
}

purity_summary <- do.call(rbind, purity_summary_rows)
print(purity_summary)

# ------------------------------------------------------------------------------
# 5) Save
# ------------------------------------------------------------------------------
if (!dir.exists("./data/processed")) dir.create("./data/processed", recursive = TRUE)

save(purity_results,
     file = paste0("./data/processed/purity_results_", qc_tag, ".rda"))
final_high_purity_tagged <- paste0("final_high_purity_sample_lists_", qc_tag)
assign(final_high_purity_tagged, final_high_purity_sample_lists)
save(list = final_high_purity_tagged,
     file = paste0("./data/processed/final_high_purity_sample_lists_", qc_tag, ".rda"))

write.csv(purity_summary,
          file = paste0("./data/processed/purity_summary_", qc_tag, ".csv"),
          row.names = FALSE)

cat("\nSaved:\n")
cat("  - purity_results_", qc_tag, ".rda\n", sep="")
cat("  - final_high_purity_sample_lists_", qc_tag, ".rda\n", sep="")
cat("  - purity_summary_", qc_tag, ".csv\n", sep="")

# ------------------------------------------------------------------------------
# 6) Notes
# ------------------------------------------------------------------------------
cat("\nNotes:\n")
cat("* Thresholds: MIN_TUMOR_PURITY=", MIN_TUMOR_PURITY,
    ", MAX_NORMAL_TUMORLIKE=", MAX_NORMAL_TUMORLIKE, "\n", sep="")
cat("* Base sample lists: ", if (USE_SENS) "v4_sens" else "v4_main", "\n", sep="")
cat("* Downstream: 08_deges_normalization_v3.R を v3' として、上流読込を本スクリプトの出力に切替。\n")

