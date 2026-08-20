# external_gene_anchors.R
# Descriptive cross-reference of externally published radiation-associated
# thyroid gene lists against this study's per-unit DEG sets (claim map C-13).
#
# The corpus has two evidence classes. The validated-anchor class, curated in
# external_gene_anchors.csv, contains the qRT-PCR-validated Abend 2013 and Dom
# 2012 lists, the qRT-PCR-confirmed Abend 2012 tumour-minus-normal list, and
# CLIP2. The multivariate-signature class, curated separately in
# ory2026_gene_signatures.csv, contains the 50-gene shared-tissue, 45-gene
# normal-tissue, and 64-gene tumour-tissue transcriptomic signatures reported
# in Ory 2026 Supplementary Tables S2, S4, and S6. The classes remain explicit
# because the Ory lists were selected by sPLS-DA/DIABLO rather than gene-wise
# testing followed by qRT-PCR confirmation.
#
# This is a membership count, not a test. For every anchor set x unit we
# report k (list genes among the unit's Storey q<0.10 DEGs) over n (list genes
# among the unit's tested genes), plus per-gene effect/p/q detail for the
# tissue-matched units and for any unit where k > 0 (so a nonzero cell is
# always self-documenting in the log). No enrichment statistic is computed:
# the platforms, contrasts, and list-construction procedures differ across
# sources, and the reading is symmetric -- any k is reported descriptively,
# and no claim of the paper moves with the outcome.
#
# Symbol handling: the CSV carries the published symbol and the current
# GENCODE symbol (renames: C1orf9->SUCO, FAM38A->PIEZO1); resolution to
# Ensembl gene ids happens here against the SE's own rowData, so the
# mapping is reproducible from repository contents. Unresolved or
# multi-mapped symbols are printed, not silently dropped.
#
# Input : diagnostics/external_gene_anchors.csv  (validated anchor lists)
#         diagnostics/ory2026_gene_signatures.csv (multivariate signatures)
#         processed/thyr_expression_test.rds     (from 410; genes table)
#         processed/thyr_se_raw.rds              (rowData: gene_id, gene_name)
# Output: diagnostics/output/external_gene_anchors.rds

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

validated <- read.csv(
  file.path(paths$root, "diagnostics", "external_gene_anchors.csv"),
  stringsAsFactors = FALSE
)
validated$evidence_class <- "validated_anchor"
validated$source_table <- ""
validated$ensembl_id <- ""

ory <- read.csv(
  file.path(paths$root, "diagnostics", "ory2026_gene_signatures.csv"),
  stringsAsFactors = FALSE
)
names(ory)[names(ory) == "signature"] <- "anchor_set"
ory$evidence_class <- "multivariate_signature"

anchor_columns <- c(
  "anchor_set", "evidence_class", "tissue", "mapping", "doi",
  "source_table", "symbol_published", "symbol_current", "ensembl_id", "note"
)
anchors <- rbind(validated[, anchor_columns], ory[, anchor_columns])
stopifnot(!anyDuplicated(anchors[, c("anchor_set", "symbol_published")]))

expression_test <- readRDS(
  file.path(paths$processed, "thyr_expression_test.rds")
)
se <- readRDS(file.path(paths$processed, "thyr_se_raw.rds"))
annot <- as.data.frame(rowData(se))[, c("gene_id", "gene_name")]
annot$ensembl_base <- sub("\\..*$", "", annot$gene_id)

cat(
  "External gene lists:", nrow(anchors), "entries,",
  length(unique(anchors$anchor_set)), "sets;",
  length(unique(anchors$symbol_current)), "unique current symbols\n"
)
cat("Evidence classes:", paste(names(table(anchors$evidence_class)),
  as.integer(table(anchors$evidence_class)), sep = "=", collapse = "; "), "\n")

# Resolve current symbols to Ensembl gene ids via the SE annotation. An
# explicit Ensembl id is used when a GENCODE symbol is duplicated.
resolution_rows <- lapply(seq_len(nrow(anchors)), function(i) {
  entry <- anchors[i, , drop = FALSE]
  explicit_id <- !is.na(entry$ensembl_id) && nzchar(entry$ensembl_id)
  idx <- if (explicit_id) {
    which(annot$ensembl_base == entry$ensembl_id)
  } else {
    which(annot$gene_name == entry$symbol_current)
  }
  if (length(idx) == 0L) {
    resolved <- data.frame(
      gene_id = NA_character_, matched_symbol = NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    resolved <- data.frame(
      gene_id = annot$gene_id[idx], matched_symbol = annot$gene_name[idx],
      stringsAsFactors = FALSE
    )
  }
  cbind(entry[rep(1L, nrow(resolved)), , drop = FALSE], resolved)
})
hits <- do.call(rbind, resolution_rows)
row.names(hits) <- NULL

unresolved <- unique(hits$symbol_current[is.na(hits$gene_id)])
resolved_key <- interaction(
  hits$anchor_set[!is.na(hits$gene_id)],
  hits$symbol_current[!is.na(hits$gene_id)], drop = TRUE
)
multi_key <- names(which(table(resolved_key) > 1L))
if (length(unresolved) > 0) {
  cat("UNRESOLVED symbols (not present in GENCODE v36; counted as not tested):",
    paste(unresolved, collapse = ", "), "\n")
}
if (length(multi_key) > 0) {
  cat("MULTI-MAPPED list entries (>1 gene_id):",
    paste(multi_key, collapse = ", "), "\n")
}
if (length(unresolved) == 0 && length(multi_key) == 0) {
  cat("All list entries resolved 1:1 against SE rowData\n")
}

units <- names(expression_test$units)
tissue_match <- list(
  normal = c("R_Normal", "B_Normal"),
  tumor = c("R_Tumor", "B_Tumor"),
  shared = units,
  pair_difference = units # indirect: shown for all units
)

# Membership summary: one row per anchor_set x unit.
summary_rows <- list()
detail_rows <- list()
for (set_name in unique(hits$anchor_set)) {
  set_hits <- hits[hits$anchor_set == set_name & !is.na(hits$gene_id), ]
  set_tissue <- unique(hits$tissue[hits$anchor_set == set_name])
  set_class <- unique(hits$evidence_class[hits$anchor_set == set_name])
  stopifnot(length(set_tissue) == 1L, length(set_class) == 1L)
  for (u in units) {
    genes <- expression_test$units[[u]]$genes
    in_tested <- set_hits[set_hits$gene_id %in% genes$gene_id, ]
    g <- genes[match(in_tested$gene_id, genes$gene_id), ]
    k <- sum(g$q_storey < 0.10)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      anchor_set = set_name, evidence_class = set_class, tissue = set_tissue,
      matched = u %in% tissue_match[[set_tissue]],
      unit = u,
      n_list = sum(hits$anchor_set == set_name),
      n_tested = nrow(in_tested), k_deg = k
    )
    if ((u %in% tissue_match[[set_tissue]] || k > 0) && nrow(in_tested) > 0) {
      detail_rows[[length(detail_rows) + 1]] <- data.frame(
        anchor_set = set_name, evidence_class = set_class, unit = u,
        symbol = in_tested$symbol_current, gene_id = in_tested$gene_id,
        effect = g$effect, p_exact = g$p_exact,
        q_storey = g$q_storey, rank = g$rank
      )
    }
  }
}
summary_tab <- do.call(rbind, summary_rows)
detail_tab <- do.call(rbind, detail_rows)

# Cell labels keep three distinctions apart: [pair-diff] marks the Abend 2012
# list whose source estimand (tumour-minus-normal difference) has no
# corresponding unit here (CSV mapping = indirect); [shared-tissue] marks the
# Ory list selected jointly from normal and tumour tissues; [tissue-matched] /
# [cross-tissue] otherwise describe whether the unit's tissue side agrees
# with the list's source tissue.
cat("\nMembership summary (k = list genes with q<0.10 in unit; n = tested):\n")
for (i in seq_len(nrow(summary_tab))) {
  s <- summary_tab[i, ]
  lab <- if (s$tissue == "pair_difference") {
    "[pair-diff]"
  } else if (s$tissue == "shared") {
    "[shared-tissue]"
  } else if (s$matched) {
    "[tissue-matched]"
  } else {
    "[cross-tissue]"
  }
  cat(sprintf(
    "  %-18s %-22s %-9s %-16s k=%d / n=%d (of %d listed)\n",
    s$anchor_set, s$evidence_class, s$unit, lab,
    s$k_deg, s$n_tested, s$n_list
  ))
}

cat("\nPer-gene detail (matched/shared units and nonzero cross-tissue cells):\n")
for (i in seq_len(nrow(detail_tab))) {
  d <- detail_tab[i, ]
  cat(sprintf(
    "  %-18s %-22s %-9s %-9s effect %.3f  p %.3g  q %.3f  rank %d\n",
    d$anchor_set, d$evidence_class, d$unit, d$symbol,
    d$effect, d$p_exact, d$q_storey, d$rank
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
