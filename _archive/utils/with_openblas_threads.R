# with_openblas_threads.R
# OpenBLASスレッド数制御関数（中間版）

with_openblas_threads <- function(threads = "auto-2", expr) {
  # RhpcBLASctlパッケージの読み込み
  if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    stop("RhpcBLASctl package is required. Install with: install.packages('RhpcBLASctl')")
  }
  library(RhpcBLASctl)
  
  # 現在の設定を保存
  old <- blas_get_num_procs()
  max_threads <- get_num_procs()
  
  # スレッド数の決定
  if (threads == "auto-1") {
    new_threads <- max(1, max_threads - 1)
    cat(sprintf("Auto-setting BLAS threads: max-1 = %d (max: %d)\n", new_threads, max_threads))
  } else if (threads == "auto-2") {
    new_threads <- max(1, max_threads - 2)
    cat(sprintf("Auto-setting BLAS threads: max-2 = %d (max: %d)\n", new_threads, max_threads))
  } else if (is.numeric(threads)) {
    new_threads <- as.integer(threads)
    if (new_threads < 1) {
      stop("Thread count must be at least 1")
    }
    if (new_threads > max_threads) {
      warning(sprintf("Requested threads (%d) exceeds maximum (%d), using maximum", 
                      new_threads, max_threads))
      new_threads <- max_threads
    }
    cat(sprintf("Setting BLAS threads to: %d (max: %d)\n", new_threads, max_threads))
  } else {
    stop("threads must be 'auto-1', 'auto-2', or a numeric value")
  }
  
  # スレッド数変更
  blas_set_num_threads(new_threads)
  cat(sprintf("BLAS threads changed: %d -> %d\n", old, new_threads))
  
  # 終了時に必ず復元
  on.exit({
    blas_set_num_threads(old)
    cat(sprintf("BLAS threads restored: %d\n", old))
  }, add = TRUE)
  
  # 式を評価
  force(expr)
}

# 使用例とテスト関数
test_with_openblas_threads <- function() {
  cat("=== Testing with_openblas_threads ===\n")
  
  if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    cat("RhpcBLASctl not available. Skipping test.\n")
    return(invisible(NULL))
  }
  
  library(RhpcBLASctl)
  
  cat("System information:\n")
  cat("  Max cores:", get_num_cores(), "\n")
  cat("  Max processors:", get_num_procs(), "\n")
  cat("  Current BLAS threads:", blas_get_num_procs(), "\n")
  
  # テスト1: auto-2
  cat("\nTest 1: auto-2\n")
  result1 <- with_openblas_threads("auto-2", {
    Sys.sleep(0.1)  # 短い処理のシミュレーション
    cat("  Inside function - current threads:", blas_get_num_procs(), "\n")
    "test_result_auto2"
  })
  cat("  Returned:", result1, "\n")
  
  # テスト2: auto-1
  cat("\nTest 2: auto-1\n")
  result2 <- with_openblas_threads("auto-1", {
    Sys.sleep(0.1)
    cat("  Inside function - current threads:", blas_get_num_procs(), "\n")
    "test_result_auto1"
  })
  cat("  Returned:", result2, "\n")
  
  # テスト3: 明示的数値
  cat("\nTest 3: explicit number (4)\n")
  result3 <- with_openblas_threads(4, {
    Sys.sleep(0.1)
    cat("  Inside function - current threads:", blas_get_num_procs(), "\n")
    "test_result_4"
  })
  cat("  Returned:", result3, "\n")
  
  cat("\nFinal BLAS threads:", blas_get_num_procs(), "\n")
  cat("=== Test completed ===\n")
  
  invisible(list(result1, result2, result3))
}