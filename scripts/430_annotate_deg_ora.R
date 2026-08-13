# 430_annotate_deg_ora.R
# Descriptive over-representation annotation of the R_Tumor DEG list --
# the composition of the already-discovered gene list, reported as
# hypothesis-generating annotation (claim map C-14; level and reading fixed
# there before the first run).
#
# Provenance: ORA was part of this study's original protocol (the Thyroid
# submission of the study's first form). It was removed during the pipeline
# rebuild on AI-era advice that the researcher does not accept as a prior
# commitment, and restored to the numbered stream as a designed analysis on
# 2026-08-13 (claim map C-14 revision notes). The computation was first
# executed 2026-08-12 as a writing-phase supplementary calculation
# (diagnostics/, superseded); algorithm, inputs and outputs are unchanged
# here, and the restored run was verified equal to the superseded one.
#
# Position: this is a characterization of the discovered list, not a second
# discovery machinery, and no claim of the paper moves with its outcome.
# Experiment-level set inference is 420 (label-permutation null,
# D6-calibrated). The weighting between the two engines is not a preference
# but a sampling-model fact (Goeman & Buhlmann 2007,
# doi:10.1093/bioinformatics/btm051): the two p-values refer to different
# randomness, and only the label-permutation null refers to the experiment
# actually performed -- gene-sampling p-values are based on a model that
# does not resemble the experiment and can be wildly anti-conservative
# under co-expression. The hypergeometric p used here is exact under its
# stated gene-sampling null; what that null ignores (co-expression) is
# disclosed wherever the results appear. Flags are read as "supported only
# under the gene-sampling null" when 420 is null for the same family; all
# results are reported including zero. The one forbidden reading is using
# these q-values as stand-alone experiment-level enrichment claims
# (objections ledger Q-15 (3b)).
#
# Design (fixed before first run): only R_Tumor (the one unit with DEGs);
# three lists -- up (effect > 0.5, high in the High arm), down, combined;
# universe = the unit's 15,621 tested genes (never the whole genome); set
# universe identical to 420 by construction (same collections loader, same
# gsea_pathway_index size window 15-500 against the same gene ids);
# one-sided hypergeometric upper tail; BH within family x list; flag at
# q_bh < 0.10 (the protocol-wide threshold -- no new number).
#
# Input : processed/thyr_expression_test.rds   (from 410; R_Tumor genes)
#         lib/gsea_collections.R               (families; size window)
#         lib/gsea_permutation.R               (gsea_pathway_index)
#         lib/annotation.R                     (strip_ensembl_version)
# Output: processed/thyr_deg_ora_annotation.rds

source("setup.R")

suppressPackageStartupMessages({
  library(msigdbr)
})

source(file.path(paths$root, "lib", "gsea_collections.R"))
source(file.path(paths$root, "lib", "gsea_permutation.R"))
source(file.path(paths$root, "lib", "annotation.R"))

SPECIES <- "Homo sapiens"
FDR_CUT <- 0.10

expression_test <- readRDS(
  file.path(paths$processed, "thyr_expression_test.rds")
)
genes <- expression_test$units$R_Tumor$genes
universe <- strip_ensembl_version(genes$gene_id)
deg <- genes$q_storey < FDR_CUT

gene_lists <- list(
  up = universe[deg & genes$effect > 0.5],
  down = universe[deg & genes$effect < 0.5],
  combined = universe[deg]
)
cat(sprintf(
  "R_Tumor universe %d tested genes; DEG lists: up %d / down %d / combined %d\n",
  length(universe), length(gene_lists$up), length(gene_lists$down),
  length(gene_lists$combined)
))

collections <- load_gene_set_collections(SPECIES)

rows <- list()
for (family in names(collections$gene_sets)) {
  # Same indexing helper and size window as 420: sets keep their members
  # present in the universe and survive only at intersected size 15-500.
  index <- gsea_pathway_index(
    collections$gene_sets[[family]], universe,
    min_size = GSEA_MIN_SET_SIZE, max_size = GSEA_MAX_SET_SIZE
  )
  cat(sprintf("Family %-18s %d sets after universe intersection\n",
              family, length(index)))
  for (list_name in names(gene_lists)) {
    hit_idx <- which(universe %in% gene_lists[[list_name]])
    k <- vapply(index, function(ix) sum(ix %in% hit_idx), integer(1))
    K <- lengths(index)
    n <- length(hit_idx)
    N <- length(universe)
    p <- phyper(k - 1L, K, N - K, n, lower.tail = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      family = family, list = list_name, pathway = names(index),
      set_size = K, overlap = k, expected = n * K / N,
      p_hyper = p, q_bh = p.adjust(p, method = "BH")
    )
  }
}
ora <- do.call(rbind, rows)
rownames(ora) <- NULL

cat("\nSummary (flag = q_bh <", FDR_CUT, "within family x list):\n")
for (family in unique(ora$family)) {
  for (list_name in names(gene_lists)) {
    blk <- ora[ora$family == family & ora$list == list_name, ]
    cat(sprintf(
      "  %-18s %-8s sets %4d | flagged %3d | min q %.3f\n",
      family, list_name, nrow(blk), sum(blk$q_bh < FDR_CUT), min(blk$q_bh)
    ))
  }
}

flagged <- ora[ora$q_bh < FDR_CUT, ]
flagged <- flagged[order(flagged$family, flagged$list, flagged$q_bh), ]
cat("\nFlagged sets (all, ordered; empty section = none):\n")
if (nrow(flagged) == 0) {
  cat("  (none)\n")
} else {
  for (i in seq_len(nrow(flagged))) {
    f <- flagged[i, ]
    cat(sprintf(
      "  %-18s %-8s %-55s k=%d/%d (exp %.1f) p %.2e q %.4f\n",
      f$family, f$list, substr(f$pathway, 1, 55), f$overlap, f$set_size,
      f$expected, f$p_hyper, f$q_bh
    ))
  }
}

result <- list(
  date = Sys.time(),
  config = list(
    unit = "R_Tumor", fdr_cut = FDR_CUT,
    universe_size = length(universe),
    list_sizes = lengths(gene_lists),
    test = "one-sided hypergeometric, BH within family x list",
    size_window = c(GSEA_MIN_SET_SIZE, GSEA_MAX_SET_SIZE)
  ),
  table = ora
)
saveRDS(result, file.path(paths$processed, "thyr_deg_ora_annotation.rds"))
cat("\nSaved:", file.path(paths$processed, "thyr_deg_ora_annotation.rds"), "\n")
