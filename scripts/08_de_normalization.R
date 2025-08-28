# 08_de_normalization_v4.R — DE用に logCPM を比較ごとに保存（TMM, filterByExpr, protein_coding）
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(edgeR); library(dplyr)
})
ASSAY_NAME <- "stranded_second"
PRIOR <- 0.5
IN_HP <- "./data/processed/final_high_purity_sample_lists_v4_main_min060_pca.rda"
OUT_DIR <- "./data/processed/de_inputs_v4"
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

stopifnot(file.exists(IN_HP)); load(IN_HP)
hp <- final_high_purity_sample_lists_v4_main_min060_pca

stopifnot(exists("se_thyr"))
cts <- assay(se_thyr, ASSAY_NAME); gi <- rowData(se_thyr)
pc_mask <- if ("gene_type" %in% colnames(gi)) !is.na(gi$gene_type) & gi$gene_type=="protein_coding" else rep(TRUE,nrow(cts))
cts <- cts[pc_mask, , drop=FALSE]

comparisons <- list(
  R0_vs_R1_tumor  = list(g1="R0", g2="R1", tissue="tumor"),
  R0_vs_R1_normal = list(g1="R0", g2="R1", tissue="normal"),
  B0_vs_B1_tumor  = list(g1="B0", g2="B1", tissue="tumor"),
  B0_vs_B1_normal = list(g1="B0", g2="B1", tissue="normal")
)

run_one <- function(cmp_name, info){
  s1 <- hp[[info$g1]][[info$tissue]]
  s2 <- hp[[info$g2]][[info$tissue]]
  smp <- c(s1, s2); smp <- intersect(smp, colnames(cts))
  if (length(s1)<3 || length(s2)<3) { cat(cmp_name, ": too few samples\n"); return(NULL) }
  
  dge <- DGEList(cts[, smp, drop=FALSE], group = factor(c(rep(info$g1,length(s1)), rep(info$g2,length(s2)))))
  keep <- filterByExpr(dge, group=dge$samples$group); dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge, "TMM")
  logcpm <- cpm(dge, log=TRUE, prior.count=PRIOR, normalized.lib.sizes=TRUE)
  
  meta <- data.frame(sample=colnames(logcpm), group=as.character(dge$samples$group), stringsAsFactors=FALSE)
  out <- list(logCPM=logcpm, meta=meta, prior=PRIOR, assay=ASSAY_NAME)
  saveRDS(out, file=file.path(OUT_DIR, paste0("de_input_", cmp_name, ".rds")))
  cat("Saved:", cmp_name, "\n")
}

invisible(lapply(names(comparisons), function(nm) run_one(nm, comparisons[[nm]])))
