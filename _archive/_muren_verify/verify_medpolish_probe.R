#!/usr/bin/env Rscript
# verify_medpolish_probe.R
# 目的: ③本検証の前の初回確認(1ケースのみ)。以下の2点を検める。
#   (A) medpolish(trace.iter=TRUE) の出力形式を目視確認し、
#       「1反復=1行」で反復回数を数えてよいかを検証する。
#   (B) 検証スクリプトが再現する rs_mx が本体 muren_norm と一致するか、
#       本体 muren_norm(res_return="scaling_coeff") の係数と、
#       検証側 maxiter=70 で得る係数を照合して裏取りする。
#
#   N=20・シード1つのみ。計算コスト最小。ここで捕捉方法とrs_mx再現の
#   正しさを確認してから、全水準版(verify_medpolish_convergence.R)を回す。

suppressWarnings(suppressMessages({
  library(data.table)
}))

## ---- 設定 ----
WORKSPACE   <- "/workspace"
GDC_DIR     <- file.path(WORKSPACE, "analysis_v7/data/raw/gdc")
UTILS_R     <- file.path(WORKSPACE, "utils/utils_improved.R")
NORM_R      <- file.path(WORKSPACE, "utils/norm_improved.R")

COUNT_COL   <- 4L
STAT_ROWS   <- 4L
N_PROBE     <- 20L
SEED        <- 101L
CPM_THRESH  <- 1
PROP_SAMPLE <- 0.5
MAXITER_MAIN<- 70L
WORKERS     <- 16L
EPS_MATCH   <- 1e-8

## BLASスレッド1固定
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
} else {
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
}

stopifnot(file.exists(UTILS_R), file.exists(NORM_R))
source(UTILS_R)
source(NORM_R)

## ---- データ読み込み関数 ----
read_counts_one <- function(path, count_col = COUNT_COL, stat_rows = STAT_ROWS) {
  dt <- data.table::fread(path, skip = 1L, header = TRUE, sep = "\t", showProgress = FALSE)
  dt <- dt[(stat_rows + 1L):nrow(dt), ]
  as.numeric(dt[[count_col]])
}
get_gene_ids <- function(path, stat_rows = STAT_ROWS) {
  dt <- data.table::fread(path, skip = 1L, header = TRUE, sep = "\t", showProgress = FALSE)
  dt <- dt[(stat_rows + 1L):nrow(dt), ]
  as.character(dt[[1L]])
}
filter_by_cpm <- function(mat, cpm_thresh = CPM_THRESH, prop = PROP_SAMPLE) {
  lib <- colSums(mat)
  cpm <- sweep(mat, 2, lib, `/`) * 1e6
  keep <- rowSums(cpm > cpm_thresh) >= (prop * ncol(mat))
  mat[keep, , drop = FALSE]
}

## ---- 抽出 ----
tsv_files <- list.files(GDC_DIR, pattern = "\\.augmented_star_gene_counts\\.tsv$",
                        recursive = TRUE, full.names = TRUE)
n_total <- length(tsv_files)
cat(sprintf("無選別サンプル総数: %d\n", n_total))
set.seed(SEED)
pick <- sample.int(n_total, N_PROBE)
sel_files <- tsv_files[pick]

gene_ids <- get_gene_ids(sel_files[1])
mat <- vapply(sel_files, read_counts_one, numeric(length(gene_ids)))
rownames(mat) <- gene_ids
mat_f <- filter_by_cpm(mat)
cat(sprintf("抽出 N=%d, フィルタ後遺伝子=%d\n\n", N_PROBE, nrow(mat_f)))

## ============================================================
## (A) medpolish trace 出力形式の目視確認
## ============================================================
cat("======== (A) medpolish trace.iter=TRUE の生出力を表示 ========\n")
# 本体同様に rs_mx を構成(段3ロジック)
options(muren_pair_method = "lts")
reads <- as.matrix(mat_f)
i_gene <- filter_gene_l(reads, 10)
reads <- reads[i_gene, , drop = FALSE]
log_mx <- lg(reads)
n_exp <- ncol(reads)

refs_list <- lapply(seq_len(n_exp), function(k) seq_len(n_exp))
pairs <- do.call(rbind, lapply(seq_len(n_exp), function(i) cbind(i, refs_list[[i]])))
locations <- pairs[, 2] + (pairs[, 1] - 1L) * n_exp

cl <- parallel::makeCluster(WORKERS, type = "PSOCK")
doSNOW::registerDoSNOW(cl)
parallel::clusterSetRNGStream(cl, 12345L)
parallel::clusterCall(cl, function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L); RhpcBLASctl::omp_set_num_threads(1L)
  } else Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
  NULL
})
parallel::clusterCall(cl, function(m) { options(muren_pair_method = m); NULL }, "lts")
parallel::clusterExport(cl, c("reg_sp", ".reg_backend"), envir = .GlobalEnv)

n_pairs <- nrow(pairs)
n_parts <- min(WORKERS, n_pairs)
split_idx <- split(seq_len(n_pairs),
                   rep(seq_len(n_parts), each = ceiling(n_pairs/n_parts), length.out = n_pairs))
`%dopar%` <- foreach::`%dopar%`
res_chunks <- foreach::foreach(
  idx = iterators::iter(split_idx),
  .packages = c("MASS","robustbase"),
  .export = c("reg_sp",".reg_backend"),
  .combine = "c"
) %dopar% {
  out <- numeric(length(idx))
  for (ii in seq_along(idx)) {
    p <- idx[ii]; i <- pairs[p,1]; j <- pairs[p,2]
    out[ii] <- reg_sp(log_mx[, i], log_mx[, j])
  }
  out
}
try(parallel::stopCluster(cl), silent = TRUE)
res_pairwise <- as.numeric(res_chunks)

rs_mx <- matrix(NA_real_, nrow = n_exp, ncol = n_exp)
rs_mx[locations] <- res_pairwise

# trace出力を「そのまま画面に」出す(sinkせず)。反復ごとの行の形を目視する。
cat("---- medpolish の生 trace 出力(以下) ----\n")
m70 <- stats::medpolish(rs_mx, na.rm = TRUE, trace.iter = TRUE, maxiter = MAXITER_MAIN)
cat("---- 生 trace 出力(以上) ----\n\n")

# 同じ呼び出しを sink で捕捉し、行数を数えてみる(本検証の捕捉方法の妥当性確認)
tc <- textConnection("mp_trace", "w", local = TRUE)
sink(tc)
m70b <- stats::medpolish(rs_mx, na.rm = TRUE, trace.iter = TRUE, maxiter = MAXITER_MAIN)
sink(); close(tc)
cat(sprintf("sink捕捉した trace 行数(=反復回数とみなす候補): %d\n", length(mp_trace)))
cat("上の生出力と行数が一致し、各行が1反復に対応しているかを目視確認すること。\n")
cat(sprintf("捕捉行の内容(先頭数行):\n"))
print(utils::head(mp_trace, 5))
cat("\n")

coef_probe <- as.numeric(m70$overall + m70$col)

## ============================================================
## (B) 本体 muren_norm との係数照合(rs_mx再現の裏取り)
## ============================================================
cat("======== (B) 本体 muren_norm(scaling_coeff) との照合 ========\n")
# 本体 muren_norm を本番同等引数で呼ぶ。返り値は 1/coef_sp。
# 検証側 coef_probe は overall+col(log2)。本体の内部変換は:
#   coef_sp = 2^(overall+col); return 1/coef_sp
# よって本体返り値 == 1/(2^coef_probe) のはず(rs_mx再現が正しければ)。
main_coeff <- muren_norm(
  mat_f,
  refs = "saturated",
  pairwise_method = "lts",
  single_param = TRUE,
  res_return = "scaling_coeff",
  workers = WORKERS
)
probe_as_scaling <- 1 / (2^coef_probe)

# 並び(サンプル順)を揃えて比較
if (length(main_coeff) != length(probe_as_scaling)) {
  cat(sprintf("長さ不一致: 本体=%d, 検証=%d\n", length(main_coeff), length(probe_as_scaling)))
} else {
  max_abs_diff <- max(abs(as.numeric(main_coeff) - probe_as_scaling))
  cat(sprintf("本体 scaling_coeff と検証(1/2^(overall+col)) の最大絶対差: %.3e\n", max_abs_diff))
  if (max_abs_diff < EPS_MATCH) {
    cat("=> 一致。検証スクリプトの rs_mx 再現は本体と等価と確認。\n")
  } else {
    cat("=> 不一致。rs_mx 再現ロジックを本体と突き合わせて要修正。\n")
    cat("   (本体先頭5値と検証先頭5値)\n")
    print(utils::head(as.numeric(main_coeff), 5))
    print(utils::head(probe_as_scaling, 5))
  }
}

cat("\n======== 初回確認まとめ ========\n")
cat("(A) trace 行数と生出力を目視し、行数=反復回数が成立するか判断。\n")
cat("(B) 本体との最大絶対差が閾値未満なら rs_mx 再現は妥当。\n")
cat("両者が確認できたら全水準版 verify_medpolish_convergence.R を実行する。\n")
cat("================================\n")