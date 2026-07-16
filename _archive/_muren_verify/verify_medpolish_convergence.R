#!/usr/bin/env Rscript
# verify_medpolish_convergence.R
# 目的: MUREN統合ステップ medpolish が maxiter=70 で収束するかの検証(③)。
#   本番データ構造に合わせず、無選別サンプル(906 tsv)からランダム抽出した
#   複数の行列で、収束の頑健性を確認する。
#
# 検証方法:
#   - 是正版 utils_improved.R / norm_improved.R を source し、本体と同一の
#     reg_sp・filter_gene_l・rs_mx構成ロジックを用いる。
#   - medpolish のみ観察用に trace.iter=TRUE で直接呼び、反復回数を捕捉。
#   - maxiter=70 と maxiter=500 で係数(overall+col)が一致するかを併せて確認。
#   - N(サンプル数)=20/60/150 × シード複数回で頑健性を見る。
#
# 論文基準: 収束は「70で収束済みか」の二値。70未満で収束し、かつ70/500で
#   係数が数値的に一致すれば、maxiter=70 は論文のLAD推定に到達しているとみなす。
#
# 前提:
#   - コンテナ rebc-r453-dev 内、/workspace にマウント。
#   - 生カウント: /workspace/analysis_v7/data/raw/gdc/<UUID>/*.augmented_star_gene_counts.tsv
#   - 是正版 MUREN: /workspace/utils/norm_improved.R, /workspace/utils/utils_improved.R
#   - BLASスレッド1固定を維持。worker=16。

suppressWarnings(suppressMessages({
  library(data.table)
}))

## ============================================================
## 設定
## ============================================================
WORKSPACE   <- "/workspace"
GDC_DIR     <- file.path(WORKSPACE, "analysis_v7/data/raw/gdc")
UTILS_R     <- file.path(WORKSPACE, "utils/utils_improved.R")
NORM_R      <- file.path(WORKSPACE, "utils/norm_improved.R")

COUNT_COL   <- 4L         # unstranded(第4列)。本番の stranded_second とは別列(本番非使用データ)。
STAT_ROWS   <- 4L         # 先頭統計4行(N_unmapped等)を除去
N_LEVELS    <- c(20L, 60L, 150L)   # サンプル数水準
SEEDS       <- c(101L, 202L, 303L) # 各水準の抽出シード
CPM_THRESH  <- 1          # 群非依存フィルタ: CPM>1 を
PROP_SAMPLE <- 0.5        # 過半数サンプルで満たす遺伝子を残す
MAXITER_MAIN<- 70L        # 本番の maxiter
MAXITER_REF <- 500L       # 収束裏取り用の十分大きい maxiter
WORKERS     <- 16L
EPS_MATCH   <- 1e-8       # 70/500 係数一致とみなす許容(数値誤差レベル)

## BLASスレッド1固定(本番方針)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
} else {
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
}

## 是正版 MUREN を source(reg_sp, filter_gene_l, lg, is.wholenumber 等)
stopifnot(file.exists(UTILS_R), file.exists(NORM_R))
source(UTILS_R)
# norm_improved.R は muren_norm を定義するが、本検証では内部関数のみ使用。
# source して reg_sp/filter_gene_l 等がグローバルに存在することを利用する。
source(NORM_R)

## ============================================================
## 段1: 無選別サンプルの列挙と、tsv読み込み関数
## ============================================================
tsv_files <- list.files(GDC_DIR, pattern = "\\.augmented_star_gene_counts\\.tsv$",
                        recursive = TRUE, full.names = TRUE)
n_total <- length(tsv_files)
cat(sprintf("無選別サンプル(tsv)総数: %d\n", n_total))
if (n_total < max(N_LEVELS)) stop("利用可能サンプル数が最大水準に満たない。")

# 1ファイル読み込み: コメント行(1行目)をskip、ヘッダ保持、先頭4統計行を除去、指定カウント列を返す
read_counts_one <- function(path, count_col = COUNT_COL, stat_rows = STAT_ROWS) {
  dt <- data.table::fread(path, skip = 1L, header = TRUE, sep = "\t",
                          showProgress = FALSE)
  dt <- dt[(stat_rows + 1L):nrow(dt), ]   # 統計行除去 -> 60660遺伝子
  as.numeric(dt[[count_col]])
}

# 遺伝子IDは段階1報告で全ファイル同一・同順(md5一致)と確認済み。
# 先頭ファイルからID取得し、以降は position 対応で cbind。
get_gene_ids <- function(path, stat_rows = STAT_ROWS) {
  dt <- data.table::fread(path, skip = 1L, header = TRUE, sep = "\t",
                          showProgress = FALSE)
  dt <- dt[(stat_rows + 1L):nrow(dt), ]
  as.character(dt[[1L]])   # gene_id 列
}

## ============================================================
## 段2: 群非依存フィルタ(CPM>1 を過半数サンプルで)
## ============================================================
filter_by_cpm <- function(mat, cpm_thresh = CPM_THRESH, prop = PROP_SAMPLE) {
  lib <- colSums(mat)
  cpm <- sweep(mat, 2, lib, `/`) * 1e6
  keep <- rowSums(cpm > cpm_thresh) >= (prop * ncol(mat))
  mat[keep, , drop = FALSE]
}

## ============================================================
## 段3: 本体と同一ロジックで rs_mx を構成し、medpolish を観察
## ============================================================
# 本体 norm_improved.R の rs_mx 構成を忠実に再現:
#   - saturated, include_self=TRUE(是正後): 各サンプルの参照は全サンプル(自己含む)
#   - pairs: (i=sample, j=ref) 全組
#   - locations = pairs[,2] + (pairs[,1]-1)*n_exp
#   - rs_mx[locations] <- res_pairwise (n_exp x n_exp)
#   - unused_refs は saturated なら空
# reg_sp は source済みの是正版を使用(options で lts)。

build_rs_mx <- function(count_mat, workers = WORKERS) {
  options(muren_pair_method = "lts")
  # 本体同様: 数値行列化 -> filter_gene_l(trim=10) -> log2(1+x)
  reads <- as.matrix(count_mat)
  i_gene <- filter_gene_l(reads, 10)      # 是正版: rowMaxs>=trim のみ
  reads <- reads[i_gene, , drop = FALSE]
  log_mx <- lg(reads)                     # log2(1+x)
  n_exp <- ncol(reads)

  # saturated, include_self=TRUE: 全ペア(自己含む)
  refs_list <- lapply(seq_len(n_exp), function(k) seq_len(n_exp))
  pairs <- do.call(rbind, lapply(seq_len(n_exp), function(i) cbind(i, refs_list[[i]])))
  locations <- pairs[, 2] + (pairs[, 1] - 1L) * n_exp

  # ペアワイズ reg_sp(並列)
  cl <- parallel::makeCluster(workers, type = "PSOCK")
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  doSNOW::registerDoSNOW(cl)
  parallel::clusterSetRNGStream(cl, 12345L)   # 本体と同じRNG固定
  invisible(parallel::clusterCall(cl, function() {
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1L); RhpcBLASctl::omp_set_num_threads(1L)
    } else Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    NULL
  }))
  invisible(parallel::clusterCall(cl, function(m) { options(muren_pair_method = m); NULL }, "lts"))
  parallel::clusterExport(cl, c("reg_sp", ".reg_backend"), envir = .GlobalEnv)

  n_pairs <- nrow(pairs)
  n_parts <- min(workers, n_pairs)
  split_idx <- split(seq_len(n_pairs),
                     rep(seq_len(n_parts), each = ceiling(n_pairs/n_parts), length.out = n_pairs))
  `%dopar%` <- foreach::`%dopar%`
  res_chunks <- foreach::foreach(
    idx = iterators::iter(split_idx),
    .packages = c("MASS","robustbase"),
    .export = ".reg_backend",
    .combine = "c"
  ) %dopar% {
    out <- numeric(length(idx))
    for (ii in seq_along(idx)) {
      p <- idx[ii]; i <- pairs[p,1]; j <- pairs[p,2]
      out[ii] <- reg_sp(log_mx[, i], log_mx[, j])
    }
    out
  }
  res_pairwise <- as.numeric(res_chunks)

  # rs_mx 構成(unused_refs は saturated で空)
  rs_mx <- matrix(NA_real_, nrow = n_exp, ncol = n_exp)
  rs_mx[locations] <- res_pairwise
  list(rs_mx = rs_mx, n_gene_kept = nrow(reads), n_exp = n_exp)
}

# medpolish を trace付きで呼び、反復回数を捕捉。
# medpolish は収束時のiter数を戻り値に持たない。trace.iter=TRUE の出力は
# 反復行("1: ...","2: ...")と最終行("Final: ...")から成る。初回確認で
# 実際の出力形式を確認済み。反復回数は "^[0-9]+:" にマッチする行のみ数え、
# "Final:" 行は除外する。
run_medpolish_traced <- function(rs_mx, maxiter) {
  tc <- textConnection("mp_trace", "w", local = TRUE)
  sink(tc)
  m <- stats::medpolish(rs_mx, na.rm = TRUE, trace.iter = TRUE, maxiter = maxiter)
  sink()
  close(tc)
  # 反復行のみカウント(Final行を除外)
  n_iter <- sum(grepl("^[0-9]+:", mp_trace))
  hit_max <- any(grepl("^Final:", mp_trace)) == FALSE  # Final無し=上限打ち切りの疑い
  list(coef = as.numeric(m$overall + m$col), n_iter = n_iter, no_final = hit_max)
}

## ============================================================
## メインループ
## ============================================================
cat("\n================ ③ medpolish 収束検証 ================\n")
cat(sprintf("水準 N = %s / シード = %s\n",
            paste(N_LEVELS, collapse=","), paste(SEEDS, collapse=",")))
cat(sprintf("フィルタ: CPM>%g を過半数(%.0f%%)サンプル / worker=%d\n\n",
            CPM_THRESH, PROP_SAMPLE*100, WORKERS))

results <- list()
row_i <- 0L
for (N in N_LEVELS) {
  for (sd in SEEDS) {
    row_i <- row_i + 1L
    set.seed(sd)
    pick <- sample.int(n_total, N)
    sel_files <- tsv_files[pick]

    # 行列組み立て(position cbind)
    gene_ids <- get_gene_ids(sel_files[1])
    mat <- vapply(sel_files, read_counts_one, numeric(length(gene_ids)))
    rownames(mat) <- gene_ids

    # 段2フィルタ
    mat_f <- filter_by_cpm(mat)

    # 段3 rs_mx 構成
    built <- build_rs_mx(mat_f)

    # medpolish: 70 と 500
    r70  <- run_medpolish_traced(built$rs_mx, MAXITER_MAIN)
    r500 <- run_medpolish_traced(built$rs_mx, MAXITER_REF)

    max_abs_diff <- max(abs(r70$coef - r500$coef))
    # 収束判定: Final行が出た(=medpolishが収束基準に達した)かつ上限未満
    converged_70 <- (!r70$no_final) && (r70$n_iter < MAXITER_MAIN)
    coef_match   <- (max_abs_diff < EPS_MATCH)

    results[[row_i]] <- data.frame(
      N = N, seed = sd,
      n_gene_kept = built$n_gene_kept,
      iter_70 = r70$n_iter, iter_500 = r500$n_iter,
      converged_70 = converged_70,
      max_abs_diff_70_500 = max_abs_diff,
      coef_match = coef_match
    )
    cat(sprintf("N=%3d seed=%3d | genes=%5d | iter@70=%2d iter@500=%3d | conv@70=%-5s | |Δ|max=%.2e | match=%s\n",
                N, sd, built$n_gene_kept, r70$n_iter, r500$n_iter,
                converged_70, max_abs_diff, coef_match))
  }
}

res_df <- do.call(rbind, results)
cat("\n================ 要約 ================\n")
print(res_df, row.names = FALSE)

all_conv  <- all(res_df$converged_70)
all_match <- all(res_df$coef_match)
cat(sprintf("\n全ケース 70未満収束: %s\n", all_conv))
cat(sprintf("全ケース 70/500係数一致: %s\n", all_match))
if (all_conv && all_match) {
  cat("=> maxiter=70 は全水準・全シードで論文のLAD推定に到達(収束済み)。\n")
} else {
  cat("=> 一部ケースで未収束または係数不一致。maxiter見直し or LP-LAD代替の検討対象。\n")
}
cat("======================================\n")