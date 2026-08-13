# external_gene_anchors.R
# Descriptive cross-reference of externally published radiation-associated
# thyroid gene lists against this study's per-unit DEG sets (claim map C-13;
# reading fixed there before this script's first run).
#
# Anchor corpus (tier A = qRT-PCR-validated core lists; curated in
# external_gene_anchors.csv with source DOIs): Abend 2013 BJC (8 normal /
# 6 tumour, dose-related within exposed cases), Abend 2012 PLoS One (11
# genes, tumour-minus-normal pair difference vs dose -- tagged "indirect"
# because no unit of this study measures a pair difference), Dom 2012 BJC
# (7 validated, exposed-vs-unexposed normal tissue), and CLIP2 (Hess 2011
# PNAS; the most heavily validated single tumour marker).
#
# This is a membership count, not a test. For every anchor set x unit we
# report k (list genes among the unit's Storey q<0.10 DEGs) over n (list
# genes among the unit's tested genes), plus per-gene effect/p/q detail for
# the tissue-matched units and for any unit where k > 0 (so a nonzero cell
# is always self-documenting in the log). No enrichment statistic is computed: the lists
# are small, the platforms and contrasts differ across sources, and the
# pre-fixed reading is symmetric -- any k is reported as description, and
# no claim of the paper moves with the outcome. Non-overlap is the norm
# between prior studies themselves (Dom 2012 found none of the earlier
# signatures enriched in its own normal-tissue contrast).
#
# Symbol handling: the CSV carries the published symbol and the current
# GENCODE symbol (renames: C1orf9->SUCO, FAM38A->PIEZO1); resolution to
# Ensembl gene ids happens here against the SE's own rowData, so the
# mapping is reproducible from repository contents. Unresolved or
# multi-mapped symbols are printed, not silently dropped.
#
# Input : diagnostics/external_gene_anchors.csv  (curated anchor lists)
#         processed/thyr_expression_test.rds     (from 410; genes table)
#         processed/thyr_se_raw.rds              (rowData: gene_id, gene_name)
# Output: diagnostics/output/external_gene_anchors.rds

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

anchors <- read.csv(
  file.path(paths$root, "diagnostics", "external_gene_anchors.csv"),
  stringsAsFactors = FALSE
)
expression_test <- readRDS(
  file.path(paths$processed, "thyr_expression_test.rds")
)
se <- readRDS(file.path(paths$processed, "thyr_se_raw.rds"))
annot <- as.data.frame(rowData(se))[, c("gene_id", "gene_name")]

cat(
  "Anchor lists:", nrow(anchors), "entries,",
  length(unique(anchors$anchor_set)), "sets;",
  length(unique(anchors$symbol_current)), "unique current symbols\n"
)

# Resolve current symbols to Ensembl gene ids via the SE annotation.
hits <- merge(
  anchors, annot,
  by.x = "symbol_current", by.y = "gene_name",
  all.x = TRUE, sort = FALSE
)
unresolved <- unique(hits$symbol_current[is.na(hits$gene_id)])
multi <- names(which(table(hits$symbol_current) >
  table(anchors$symbol_current)[names(table(hits$symbol_current))]))
if (length(unresolved) > 0) {
  cat("UNRESOLVED symbols (no gene_name match):",
      paste(unresolved, collapse = ", "), "\n")
}
if (length(multi) > 0) {
  cat("MULTI-MAPPED symbols (>1 gene_id):",
      paste(multi, collapse = ", "), "\n")
}
if (length(unresolved) == 0 && length(multi) == 0) {
  cat("All symbols resolved 1:1 against SE rowData\n")
}

units <- names(expression_test$units)
tissue_match <- list(
  normal = c("R_Normal", "B_Normal"),
  tumor = c("R_Tumor", "B_Tumor"),
  pair_difference = units # indirect: shown for all units
)

# Membership summary: one row per anchor_set x unit.
summary_rows <- list()
detail_rows <- list()
for (set_name in unique(hits$anchor_set)) {
  set_hits <- hits[hits$anchor_set == set_name & !is.na(hits$gene_id), ]
  set_tissue <- unique(hits$tissue[hits$anchor_set == set_name])
  for (u in units) {
    genes <- expression_test$units[[u]]$genes
    in_tested <- set_hits[set_hits$gene_id %in% genes$gene_id, ]
    g <- genes[match(in_tested$gene_id, genes$gene_id), ]
    k <- sum(g$q_storey < 0.10)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      anchor_set = set_name, tissue = set_tissue,
      matched = u %in% tissue_match[[set_tissue]],
      unit = u,
      n_list = sum(hits$anchor_set == set_name),
      n_tested = nrow(in_tested), k_deg = k
    )
    if ((u %in% tissue_match[[set_tissue]] || k > 0) && nrow(in_tested) > 0) {
      detail_rows[[length(detail_rows) + 1]] <- data.frame(
        anchor_set = set_name, unit = u,
        symbol = in_tested$symbol_current, gene_id = in_tested$gene_id,
        effect = g$effect, p_exact = g$p_exact,
        q_storey = g$q_storey, rank = g$rank
      )
    }
  }
}
summary_tab <- do.call(rbind, summary_rows)
detail_tab <- do.call(rbind, detail_rows)

# Cell labels keep two distinctions apart: [pair-diff] marks the Abend 2012
# list whose source estimand (tumour-minus-normal difference) has no
# corresponding unit here (CSV mapping = indirect); [tissue-matched] /
# [cross-tissue] describe whether the unit's tissue side agrees with the
# list's source tissue.
cat("\nMembership summary (k = list genes with q<0.10 in unit; n = tested):\n")
for (i in seq_len(nrow(summary_tab))) {
  s <- summary_tab[i, ]
  lab <- if (s$tissue == "pair_difference") {
    "[pair-diff]"
  } else if (s$matched) {
    "[tissue-matched]"
  } else {
    "[cross-tissue]"
  }
  cat(sprintf(
    "  %-18s %-9s %-16s k=%d / n=%d (of %d listed)\n",
    s$anchor_set, s$unit, lab, s$k_deg, s$n_tested, s$n_list
  ))
}

cat("\nPer-gene detail (tissue-matched units):\n")
for (i in seq_len(nrow(detail_tab))) {
  d <- detail_tab[i, ]
  cat(sprintf(
    "  %-18s %-9s %-9s effect %.3f  p %.3g  q %.3f  rank %d\n",
    d$anchor_set, d$unit, d$symbol, d$effect, d$p_exact, d$q_storey, d$rank
  ))
}

result <- list(
  date = Sys.time(),
  anchors = anchors,
  resolution = hits,
  summary = summary_tab,
  detail = detail_tab
)
out_dir <- file.path(paths$root, "diagnostics", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(result, file.path(out_dir, "external_gene_anchors.rds"))
cat("\nSaved:", file.path(out_dir, "external_gene_anchors.rds"), "\n")
