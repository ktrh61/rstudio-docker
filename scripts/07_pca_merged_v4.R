# 07_pca_merged_v4.R — R0+R1 / B0+B1 の統合PCA（CDM, TMM→logCPM, protein_coding, ゼロ残す）
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(edgeR); library(dplyr); library(Rcpp)
})
source("./utils/with_openblas_threads.R")
sourceCpp("./utils/CDM_fast3_arma_enhanced.cpp")  # v4で使ってるCDM

ASSAY_NAME <- "stranded_second"
PRIOR <- 0.5
TOPK  <- NULL   # 例: 5000 にすると上位5k遺伝子に絞る。NULLで全遺伝子
IN_RDA <- "./data/processed/final_high_purity_sample_lists_v4_main_min060_pca.rda"
OUT_RDA <- "./data/processed/pca_merged_v4_results.rda"

stopifnot(file.exists(IN_RDA))
load(IN_RDA)  # => final_high_purity_sample_lists_v4_main_min060_pca
hp_lists <- final_high_purity_sample_lists_v4_main_min060_pca

stopifnot(exists("se_thyr"))
cts <- assay(se_thyr, ASSAY_NAME); gi <- rowData(se_thyr)
pc_mask <- if ("gene_type" %in% colnames(gi)) !is.na(gi$gene_type) & gi$gene_type=="protein_coding" else rep(TRUE,nrow(cts))
cts <- cts[pc_mask, , drop=FALSE]

make_pca <- function(groups, label){
  smp <- unlist(lapply(groups, function(g) c(hp_lists[[g]]$normal, hp_lists[[g]]$tumor)))
  smp <- intersect(smp, colnames(cts))
  stopifnot(length(smp) >= 6)
  
  dge <- DGEList(cts[, smp, drop=FALSE])
  keep <- filterByExpr(dge); dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge, "TMM")
  
  logcpm <- cpm(dge, log=TRUE, prior.count=PRIOR, normalized.lib.sizes=TRUE)
  if (!is.null(TOPK)) {
    rv <- matrixStats::rowVars(as.matrix(logcpm))
    idx <- order(rv, decreasing=TRUE)[seq_len(min(TOPK, length(rv)))]
    logcpm <- logcpm[idx, , drop=FALSE]
  }
  # CDM（genes x samples を渡す）
  res <- with_openblas_threads("auto-2", { CDM_fast3_arma(logcpm, verbose=FALSE) })
  # 名付け
  colnames(res$scores) <- paste0("PC", seq_len(ncol(res$scores)))
  rownames(res$scores) <- colnames(logcpm)
  ve <- res$values^2; ve <- ve/sum(ve)
  
  # 付帯情報
  grp <- rep(groups, times = sapply(groups, function(g) length(c(hp_lists[[g]]$normal, hp_lists[[g]]$tumor))))
  grp <- setNames(grp, unlist(lapply(groups, function(g) c(hp_lists[[g]]$normal, hp_lists[[g]]$tumor))))
  grp <- grp[rownames(res$scores)]
  tis <- ifelse(grepl("_merged$", names(grp)), NA, NA) # 任意。必要なら別メタをここで付ける
  
  list(
    label = label,
    samples = rownames(res$scores),
    scores = res$scores,
    variance_explained = ve,
    k = length(res$values)
  )
}

r01 <- make_pca(c("R0","R1"), "R0+R1")
b01 <- make_pca(c("B0","B1"), "B0+B1")

dir.create("./data/processed", showWarnings=FALSE, recursive=TRUE)
pca_merged_v4_results <- list(r01=r01, b01=b01, assay=ASSAY_NAME, prior=PRIOR, topk=TOPK)
save(pca_merged_v4_results, file=OUT_RDA)
cat("Saved:", OUT_RDA, "\n")
