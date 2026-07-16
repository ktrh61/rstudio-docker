# 09_bm_test_v4.R — logCPMに対するBrunner–Munzel検定 + Storey qvalue
suppressPackageStartupMessages({
  library(brunnermunzel)  # なければ: install.packages("brunnermunzel")
  library(qvalue); library(matrixStats)
})

IN_DIR  <- "./data/processed/de_inputs_v4"
OUT_DIR <- "./data/processed/de_results_v4"
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

files <- list.files(IN_DIR, pattern="^de_input_.*\\.rds$", full.names=TRUE)
bm_one <- function(path){
  dat <- readRDS(path)
  X <- dat$logCPM; meta <- dat$meta
  g <- factor(meta$group)
  if (length(unique(g))!=2) { cat("skip:", basename(path), "\n"); return(NULL) }
  idx1 <- which(g==levels(g)[1]); idx2 <- which(g==levels(g)[2])
  if (length(idx1)<3 || length(idx2)<3) { cat("too few per arm:", basename(path), "\n"); return(NULL) }
  
  # geneごとにBM（両側）
  pvec <- numeric(nrow(X)); est <- numeric(nrow(X))
  for (i in seq_len(nrow(X))){
    x1 <- as.numeric(X[i, idx1]); x2 <- as.numeric(X[i, idx2])
    # 変な値の防御
    if (!all(is.finite(x1)) || !all(is.finite(x2))) { pvec[i] <- NA; est[i] <- NA; next }
    res <- try(brunnermunzel.test(x1, x2, alternative="two.sided"), silent=TRUE)
    if (inherits(res,"try-error") || is.na(res$p.value)) { pvec[i] <- NA; est[i] <- NA
    } else { pvec[i] <- res$p.value; est[i] <- res$estimate } # estimate ~ P(X>Y) の推定
  }
  
  ok <- is.finite(pvec)
  qv <- rep(NA_real_, length(pvec))
  if (sum(ok) > 10) {
    qv[ok] <- tryCatch(qvalue(pvec[ok])$qvalues, error=function(e) p.adjust(pvec[ok], method="BH"))
  }
  
  res_df <- data.frame(
    gene = rownames(X),
    pval = pvec,
    qval = qv,
    bm_est = est,            # “効果方向”の参考（>0.5 で group2 優位、という解釈に注意）
    stringsAsFactors=FALSE
  )
  out_csv <- file.path(OUT_DIR, sub("^de_input_","de_bm_", sub("\\.rds$",".csv", basename(path))))
  write.csv(res_df, out_csv, row.names=FALSE)
  cat("Saved:", out_csv, "\n")
}

invisible(lapply(files, bm_one))
