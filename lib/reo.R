# log2 TPM from a SummarizedExperiment's stranded_second counts, for REO.
#
# REO compares gene i vs gene j WITHIN a sample, so expression must be length-
# normalized (TPM); raw or CPM counts confound the within-sample order with gene
# length. The per-sample TPM scaling (dividing by the sample total) cancels in a
# pair ratio r = log2(TPM_i / TPM_j) = log2(rate_i / rate_j), so REO is free of
# library-size / composition normalization -- only gene length matters. TPM is
# kept here (harmless, matches the original panel) over the stranded_second
# assay chosen upstream; the unstranded TPM must not be used.
#
# All three REO scripts (510/520/530) call this so a panel built on one sample
# set is applied to another with an identical expression definition.
reo_log2_tpm <- function(se, gene_lengths, samples) {
  counts <- SummarizedExperiment::assay(se, "stranded_second")[, samples, drop = FALSE]
  ensembl <- sub("\\..*$", "", rownames(counts))
  has_len <- ensembl %in% names(gene_lengths)
  counts <- counts[has_len, , drop = FALSE]
  ensembl <- ensembl[has_len]
  len_kb <- gene_lengths[ensembl] / 1000
  rate <- counts / len_kb
  tpm <- sweep(rate, 2, colSums(rate), "/") * 1e6
  rownames(tpm) <- ensembl
  log2(tpm)
}

# Per-sample panel reversal score: the number of panel pairs whose within-
# sample order (|r| >= dead_zone) is opposite to the R0 reference sign.
# dead_zone is an explicit argument; callers pass the config value (520) or the
# value recorded in the saved panel artifact (530).
reversal_score <- function(l2tpm, samples, panel, dead_zone) {
  sc <- integer(length(samples))
  names(sc) <- samples
  for (k in seq_len(nrow(panel))) {
    r <- l2tpm[panel$up[k], samples] - l2tpm[panel$down[k], samples]
    sc <- sc + as.integer(abs(r) >= dead_zone & sign(r) != panel$r0_sign[k])
  }
  sc
}

# Read-A classification against the R0-based boundary.
classify_reversal <- function(score, boundary) {
  ifelse(score > boundary$negative_max, "positive", "negative")
}
