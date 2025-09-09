# utils/norm.R — MUREN core (consolidated fast path)
# 依存: foreach, doSNOW, iterators, parallel, MASS（任意: robustbase）
# 事前に utils/utils.R を source（.reg_backend, reg_sp, mode_sp, reg_dp 等）

if (!requireNamespace("assertthat", quietly = TRUE)) stop("Package 'assertthat' is required.")
if (!requireNamespace("foreach", quietly = TRUE))    stop("Package 'foreach' is required.")
if (!requireNamespace("doSNOW", quietly = TRUE))     stop("Package 'doSNOW' is required.")
if (!requireNamespace("iterators", quietly = TRUE))  stop("Package 'iterators' is required.")
if (!requireNamespace("MASS", quietly = TRUE))       stop("Package 'MASS' is required.")
if (!requireNamespace("parallel", quietly = TRUE))   stop("Package 'parallel' is required.")
`%dopar%` <- foreach::`%dopar%`

muren_norm <- function(reads,
                       refs = 'saturated',
                       pairwise_method = "lts",   # "lts","mode","median","trim10","huber"
                       refs_cap = Inf,            # 参照本数の上限（中央値近傍から抽出）
                       single_param = TRUE,
                       res_return = 'counts',     # "counts" | "log_counts" | "scaling_coeff"(spのみ)
                       filter_gene = TRUE,
                       trim = 10,
                       maxiter = 70,
                       workers = 2,               # 数値 or "auto" or 既存cluster
                       include_self = FALSE,
                       ...) {
  
  # reg_sp が読むスイッチを options で渡す
  old_m <- getOption("muren_pair_method", NULL)
  on.exit({ options(muren_pair_method = old_m) }, add = TRUE)
  options(muren_pair_method = pairwise_method)
  
  # BLAS は 1 スレ（外側並列と二重化しない）
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
  }
  
  # ---- 引数チェック ----
  assertthat::assert_that(is.logical(single_param))
  ok_methods <- c("lts","mode","median","trim10","huber")
  assertthat::assert_that(is.character(pairwise_method), pairwise_method %in% ok_methods)
  assertthat::assert_that(maxiter > 0 & is.wholenumber(maxiter))
  assertthat::assert_that((is.numeric(workers) && workers >= 0 && is.wholenumber(workers)) ||
                            (is.character(workers) && workers %in% c("auto")) ||
                            inherits(workers, "cluster"))
  assertthat::assert_that(is.data.frame(reads) | is.matrix(reads))
  assertthat::assert_that(trim >= 0)
  assertthat::assert_that(nrow(reads) > 1, msg = "Bad input data !")
  
  raw_reads <- reads
  
  # 数値列のみ抽出
  i_sample <- rep(TRUE, ncol(raw_reads))
  if (is.data.frame(raw_reads)) i_sample <- sapply(raw_reads, is.numeric)
  reads <- as.matrix(raw_reads[, i_sample, drop = FALSE])
  
  # 遺伝子フィルタ（群非依存の軽い条件）
  i_gene <- filter_gene_l(reads, trim)
  if (filter_gene) reads <- reads[i_gene, , drop = FALSE]
  
  # log2(1+x)
  log_raw_reads_mx <- lg(reads)
  
  n_exp  <- ncol(reads)
  n_gene <- nrow(reads)
  REFSERR <- "Bad specification of references !"
  
  # ---- refs の解釈（自己参照を除外）----
  get_refs <- NULL
  if (is.character(refs)) {
    if (length(refs) == 1) {
      if (refs == "saturated") {
        get_refs <- function(k) if (include_self) seq_len(n_exp) else setdiff(seq_len(n_exp), k)
      } else if (refs %in% colnames(reads)) {
        ref_idx <- which(colnames(reads) %in% refs)
        get_refs <- function(k) { v <- setdiff(ref_idx, k); if (length(v)==0) ref_idx else v }
      } else stop(REFSERR)
    } else {
      if (all(refs %in% colnames(reads))) {
        ref_idx <- which(colnames(reads) %in% refs)
        get_refs <- function(k) { v <- setdiff(ref_idx, k); if (length(v)==0) ref_idx else v }
      } else stop(REFSERR)
    }
  } else if (all(is.wholenumber(refs))) {
    if (max(refs) <= n_exp) {
      if (length(refs) > 1) {
        get_refs <- function(k) { v <- setdiff(refs, k); if (length(v)==0) refs else v }
      } else {
        single <- as.integer(refs)
        get_refs <- function(k) single
      }
    } else stop(REFSERR)
  } else stop(REFSERR)
  
  # 回帰メソッドのラッパー名（実処理は reg_sp が options を読む）
  reg_wapper <- if (single_param) { if (pairwise_method == 'mode') 'mode_sp' else 'reg_sp' } else 'reg_dp'
  
  # ---- 各サンプルの参照セットを作る ----
  refs_list <- lapply(seq_len(n_exp), get_refs)
  
  # ---- refs_cap: ライブラリサイズ中央値近傍から上限本数に絞る ----
  if (is.finite(refs_cap)) {
    lib <- colSums(reads)
    pick_k <- function(r, k) {
      if (!length(r)) return(r)
      m <- stats::median(lib[r])
      r[order(abs(lib[r] - m))][seq_len(min(length(r), k))]
    }
    refs_list <- lapply(refs_list, pick_k, k = as.integer(refs_cap))
  }
  
  # cap 後に used/unused を確定
  used_refs   <- sort(unique(unlist(refs_list)))
  unused_refs <- setdiff(seq_len(n_exp), used_refs)
  
  # ---- ペアインデックス ----
  pairs <- do.call(rbind, lapply(seq_len(n_exp), function(i) {
    if (length(refs_list[[i]]) == 0) return(NULL)
    cbind(i, refs_list[[i]])
  }))
  if (is.null(pairs)) stop("No valid (sample, ref) pairs were generated.")
  locations <- pairs[,2] + (pairs[,1] - 1L) * n_exp
  
  # ---- 並列設定（auto/cluster対応 & 再現性）----
  own_cluster <- TRUE
  if ((is.character(workers) && workers == "auto") || (is.numeric(workers) && workers == 0L)) {
    workers <- max(1L, parallel::detectCores() - 1L)
  }
  if (inherits(workers, "cluster")) {
    cl <- workers; own_cluster <- FALSE; workers <- length(cl)
  } else {
    workers <- min(as.integer(workers), max(1L, n_exp))
    cl <- parallel::makeCluster(workers, type = "PSOCK")
  }
  doSNOW::registerDoSNOW(cl)
  parallel::clusterSetRNGStream(cl, 12345L)
  on.exit(if (isTRUE(own_cluster)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  
  # ★ 各ワーカーでBLAS/OMPを1スレッドに固定（外側並列との二重化防止）
  parallel::clusterCall(cl, function() {
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1L)
      RhpcBLASctl::omp_set_num_threads(1L)
    } else {
      Sys.setenv(OMP_NUM_THREADS = "1",
                 OPENBLAS_NUM_THREADS = "1",
                 MKL_NUM_THREADS = "1",
                 BLIS_NUM_THREADS = "1",
                 VECLIB_MAXIMUM_THREADS = "1")
    }
    NULL
  })
  
  # ★ これを追加：各ワーカーに pairwise_method を伝える
  parallel::clusterCall(cl, function(m) { options(muren_pair_method = m); NULL }, pairwise_method)
  
  # ワーカーに必要な関数を配布
  needed <- c("reg_sp","mode_sp","reg_dp",".reg_backend",
              "polish_one_gene","polish_coeff","lg","ep","TOL")
  missing <- needed[!vapply(needed, exists, logical(1), inherits = FALSE)]
  if (length(missing) > 0) {
    stop("Required MUREN helpers not found: ", paste(missing, collapse = ", "))
  }
  parallel::clusterExport(cl, varlist = needed, envir = environment())
  
  pkgs <- c("MASS")
  if (requireNamespace("robustbase", quietly = TRUE)) pkgs <- c(pkgs, "robustbase")
  
  # ---- チャンク配布で並列実行 ----
  n_pairs <- nrow(pairs)
  n_parts <- min(workers, n_pairs)
  split_idx <- split(seq_len(n_pairs),
                     rep(seq_len(n_parts), each = ceiling(n_pairs / n_parts), length.out = n_pairs))
  
  if (single_param) {
    # reg_sp / mode_sp：各ペア 1 スカラー → ベクトル連結
    res_chunks <- foreach::foreach(
      idx = iterators::iter(split_idx),
      .packages = pkgs,
      # ← 関数だけを明示。データ系は foreach の自動捕捉に任せる
      .export   = c("reg_sp","mode_sp",".reg_backend"),
      .combine  = "c"
    ) %dopar% {
      out <- numeric(length(idx))
      for (ii in seq_along(idx)) {
        p <- idx[ii]; i <- pairs[p, 1]; j <- pairs[p, 2]
        if (reg_wapper == 'reg_sp')      out[ii] <- reg_sp (log_raw_reads_mx[, i], log_raw_reads_mx[, j], ...)
        else                              out[ii] <- mode_sp(log_raw_reads_mx[, i], log_raw_reads_mx[, j], ...)
      }
      out
    }
    res_pairwise <- as.numeric(res_chunks)
    
  } else {
    # reg_dp：各ペア n_gene ベクトル → 行列cbind
    res_chunks <- foreach::foreach(
      idx = iterators::iter(split_idx),
      .packages = pkgs,
      .export   = c("reg_dp",".reg_backend",
                    "log_raw_reads_mx","pairs","n_gene"),
      .combine  = "cbind"
    ) %dopar% {
      out <- matrix(NA_real_, n_gene, length(idx))
      for (ii in seq_along(idx)) {
        p <- idx[ii]; i <- pairs[p, 1]; j <- pairs[p, 2]
        out[, ii] <- reg_dp(log_raw_reads_mx[, i], log_raw_reads_mx[, j], ...)
      }
      out
    }
    res_pairwise <- as.matrix(res_chunks)
  }
  
  # ---- まとめ ----
  if (single_param) {
    if (length(res_pairwise) != length(locations))
      stop("Pairwise result length mismatch (single_param).")
    
    # サンプル効果（log2係数）を推定
    coef_sp <- polish_coeff(
      fitted_n    = res_pairwise,
      n_exp       = n_exp,
      locations   = locations,
      unused_refs = unused_refs,
      maxiter     = maxiter
    )
    coef_sp <- 2^(as.vector(coef_sp))
    names(coef_sp) <- colnames(reads)
    
    if (res_return == 'scaling_coeff') return(1/coef_sp)
    
    # counts をスケール（diag生成を避ける）
    polished_mx <- sweep(as.matrix(raw_reads[, i_sample, drop = FALSE]), 2, coef_sp, `*`)
    if (res_return == 'log_counts') polished_mx <- lg(polished_mx)
    
  } else {
    if (!is.matrix(res_pairwise) || nrow(res_pairwise) != n_gene)
      stop("Pairwise result shape mismatch (double_param).")
    
    # gene-wise polish
    polished_mx <- foreach::foreach(
      n = iterators::iter(seq_len(n_gene)),
      .combine = rbind,
      .export  = c("polish_one_gene")
    ) %dopar% {
      polish_one_gene(
        fitted_n    = res_pairwise[n, ],
        n_exp       = n_exp,
        locations   = locations,
        unused_refs = unused_refs,
        maxiter     = maxiter
      )
    }
    
    # 0→0 と負の丸め
    polished_mx[log_raw_reads_mx < TOL | polished_mx < 0] <- 0
    
    if (res_return == 'counts') {
      polished_mx <- ep(polished_mx)  # 2^x - 1
    } else if (res_return == 'log_counts') {
      # そのまま
    } else if (res_return == 'scaling_coeff') {
      stop("res_return='scaling_coeff' is only valid when single_param=TRUE.")
    }
  }
  
  # ---- 付帯情報/型 ----
  if (!is.null(rownames(raw_reads))) {
    if (single_param) rownames(polished_mx) <- rownames(raw_reads)
    else              rownames(polished_mx) <- rownames(reads)
  }
  colnames(polished_mx) <- colnames(reads)
  
  if (is.data.frame(raw_reads)) {
    res_df <- raw_reads
    if (single_param) {
      res_df[, i_sample] <- polished_mx
    } else {
      res_df <- res_df[i_gene, , drop = FALSE]
      res_df[, i_sample] <- polished_mx
    }
    return(res_df)
  } else {
    return(as.matrix(polished_mx))
  }
}