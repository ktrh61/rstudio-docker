# 420_test_gene_sets.R
# Gene-set enrichment for each analysis unit, on the ranking produced by the
# Brunner-Munzel contrast, under the spec-B protocol (reorg plan v2 D2/D3/D7,
# decided 2026-08-07). The null comes from the same label shuffles 410 used
# (its saved perm_index, consumed by reference), so the gene-level and
# set-level results rest on one permutation null.
# Input : processed/thyr_normalized_counts.rds  (from 310; per-unit DGEList)
#         processed/thyr_expression_test.rds    (from 410; perm_index)
#         lib/stat_brunnermunzel.R, lib/gsea_permutation.R,
#         lib/gsea_collections.R
# Output: processed/thyr_enrichment_test.rds
#           list(date, config, collections, units)
#           units = { R_Tumor, R_Normal, B_Tumor, B_Normal }; each a data frame
#             sets : collection, pathway, size, ES, NES, pval, q_bh,
#                    redundant_with, leading_edge
#
# GSEA consumes the whole ranking. That choice is not a fallback: it is
# threshold-free (no DEG-list cut whose membership sits on the q cliff decides
# what is tested) and the label-shuffle null carries the inter-gene dependence
# exactly, which is what a diffuse-signal hypothesis requires.
#
# Ranking metric: tie-averaged normal scores of the signed Brunner-Munzel
# statistic (lib/gsea_permutation.R header). Enrichment score: weighted
# running sum at gseaParam = 1, evaluated at tie-block boundaries only, so no
# arbitrary tie-break injects order the data does not contain; on tie-free
# input it is exactly the standard GSEA statistic (verified in
# tests/testthat/test-gsea-block-es.R).
#
# Inference: per-set permutation p (each set against its own null row,
# sign-conditional) with BH inside each of the four families, q_bh < FDR_CUT
# (0.10). BH is the Storey q machinery with pi0 held at its conservative
# bound 1: at this level the pi0 plug-in is variance-dominated (dependent
# p-values riding one collective-drift mode) and errs anti-conservative in
# exactly the realizations where the procedure is most fragile, so it is not
# estimated. This is the procedure the held-out null calibration passed
# (global-null P(>=1) 0.045 pooled against nominal 0.10); the pooled
# tail-ratio FDR and a restandardized variant were measured miscalibrated
# (0.14 / 0.22) and are not used -- measurements in diagnostics/output/,
# decision record in reorg plan v2 appendix B. No cross-family claim is
# made; q < 0.25 is not used, not even as an exploratory bar.
#
# Families: H, C2:CP, C5:GO:BP, and C2:CGP:radiation -- an exploratory family
# whose curation regex was fixed before this inference touched real data
# (lib/gsea_collections.R, with the construct caveat recorded there). C6/C7
# stay excluded on relevance grounds. Everything is computed within
# collection, so adding a family cannot affect another family's q-values.
#
# The change from the previous spec (gseaParam = 0, Westfall-Young FWER as
# primary inference) is a protocol amendment, not a bug fix; appendix B keeps
# the before/after specs and what had been seen when. The spike-in positive
# control is re-run under this inference (diagnostics/gsea_spikein_control.R).
#
# The size window drops both noise-prone tiny sets and the vague giant ones;
# redundant_with flags a set whose leading edge is largely contained in a
# better-ranked set (block-granular edges), gated at q_bh < FDR_CUT.

source("setup.R")

suppressPackageStartupMessages({
  library(edgeR)
  library(msigdbr)
  library(parallel)
})

source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))
source(file.path(paths$root, "lib", "gsea_permutation.R"))
source(file.path(paths$root, "lib", "gsea_collections.R"))
source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "annotation.R"))

# --- Configuration ---------------------------------------------------------
SPECIES <- "Homo sapiens"
# The ES weight (gseaParam = 1, weight = |normal score|) is intrinsic to the
# block ES in lib/gsea_permutation.R, not a parameter of this script; the
# saved config records it as the literal it is.
MIN_SET_SIZE <- GSEA_MIN_SET_SIZE # lib/gsea_collections.R; shared with the
MAX_SET_SIZE <- GSEA_MAX_SET_SIZE # null calibration diagnostic
REDUNDANT_JACCARD <- 0.5
# FDR_CUT / WORKERS come from config.R via setup.R. 420 computes BM
# statistics only (no p-value path), so the exact-enumeration options that
# 310/410 set are not set here -- they would be inert knobs.

pin_blas_threads()

# --- Load inputs -----------------------------------------------------------
norm_path <- file.path(paths$processed, "thyr_normalized_counts.rds")
test_path <- file.path(paths$processed, "thyr_expression_test.rds")
if (!file.exists(norm_path)) stop("thyr_normalized_counts.rds not found (310)")
if (!file.exists(test_path)) stop("thyr_expression_test.rds not found (410)")
normalized <- readRDS(norm_path)
expression_test <- readRDS(test_path)
message(
  "Reusing 410's permutation index: ",
  expression_test$config$n_perm, " shuffles per unit"
)

# --- Gene sets -------------------------------------------------------------
collections <- load_gene_set_collections(SPECIES)
gene_sets <- collections$gene_sets
message("Gene sets before filtering: ", paste(
  names(gene_sets), lengths(gene_sets), sep = "=", collapse = " "
))
message(
  "C2:CGP:radiation curation: ",
  length(collections$curation$radiation_sets), " sets"
)

# --- Per-unit enrichment ---------------------------------------------------
test_unit <- function(dgelist, unit) {
  arms <- unit_arms(dgelist$samples$group, unit)
  sporadic <- arms$sporadic
  high <- arms$high
  cpm_matrix <- edgeR::cpm(
    dgelist,
    normalized.lib.sizes = TRUE, prior.count = 0, log = FALSE
  )[, c(sporadic, high), drop = FALSE]
  storage.mode(cpm_matrix) <- "double"
  # MSigDB carries unversioned Ensembl identifiers.
  ensembl <- strip_ensembl_version(rownames(cpm_matrix))
  cpm_matrix <- cpm_matrix[!duplicated(ensembl), , drop = FALSE]
  rownames(cpm_matrix) <- ensembl[!duplicated(ensembl)]
  nx <- length(sporadic)

  perm_index <- expression_test$units[[unit]]$perm_index
  n_perm <- ncol(perm_index)
  if (nrow(perm_index) != ncol(cpm_matrix)) {
    stop(unit, ": perm_index rows do not match the unit's sample count.")
  }

  index <- lapply(gene_sets, gsea_pathway_index,
    gene_ids = rownames(cpm_matrix),
    min_size = MIN_SET_SIZE, max_size = MAX_SET_SIZE
  )
  collection_of <- rep(names(index), lengths(index))
  flat <- unlist(unname(index), recursive = FALSE)

  scores <- gsea_normal_scores(brunnermunzel_statistics(cpm_matrix, nx))
  observed <- gsea_block_scores(scores, flat)

  null_columns <- mclapply(seq_len(n_perm), function(i) {
    gsea_block_scores(
      gsea_normal_scores(brunnermunzel_statistics(
        cpm_matrix[, perm_index[, i], drop = FALSE], nx
      )),
      flat
    )
  }, mc.cores = WORKERS)
  null <- gsea_bind_null_columns(null_columns, length(flat))
  rownames(null) <- names(flat)

  leading_edge <- gsea_block_leading_edge(scores, flat, rownames(cpm_matrix))
  per_collection <- lapply(names(index), function(collection) {
    keep <- collection_of == collection
    standardized <- gsea_nes(observed[keep], null[keep, , drop = FALSE])
    pval <- gsea_pathway_pvalues(standardized$nes, standardized$nes_null)
    result <- data.frame(
      collection = collection,
      pathway = names(observed)[keep],
      size = lengths(flat[keep]),
      ES = unname(observed[keep]),
      NES = unname(standardized$nes),
      pval = pval,
      q_bh = p.adjust(pval, method = "BH"),
      stringsAsFactors = FALSE
    )
    result$redundant_with <- gsea_redundancy(
      result, leading_edge, REDUNDANT_JACCARD,
      candidate = !is.na(result$q_bh) & result$q_bh < FDR_CUT
    )
    result$leading_edge <- I(unname(leading_edge[result$pathway]))
    result[order(result$q_bh, result$pval, -abs(result$NES)), , drop = FALSE]
  })

  sets <- do.call(rbind, per_collection)
  rownames(sets) <- NULL

  per_family <- vapply(names(index), function(collection) {
    rows <- sets[sets$collection == collection, , drop = FALSE]
    sprintf(
      "%s q<%.2f %d (min %.3f)", collection, FDR_CUT,
      sum(rows$q_bh < FDR_CUT, na.rm = TRUE),
      min(rows$q_bh, na.rm = TRUE)
    )
  }, character(1))
  message(sprintf(
    "  %-9s %d sets | %s", unit, nrow(sets),
    paste(per_family, collapse = " | ")
  ))
  sets
}

message("Testing units:")
units <- lapply(names(normalized$units), function(unit) {
  test_unit(normalized$units[[unit]]$dgelist, unit)
})
names(units) <- names(normalized$units)

# --- Assemble and save -----------------------------------------------------
thyr_enrichment_test <- list(
  date = Sys.Date(),
  config = list(
    ranking = "tie-averaged normal scores of the signed BM statistic",
    es = "tie-block invariant weighted running sum, gseaParam = 1",
    permutation = "sample labels; 410's saved perm_index, by reference",
    n_perm = expression_test$config$n_perm,
    perm_seed = expression_test$config$perm_seed,
    perm_index_hash = vapply(
      expression_test$units, function(u) u$perm_index_hash, character(1)
    ),
    inference = paste(
      "per-set permutation p + BH within collection, q_bh <", FDR_CUT,
      "(Storey q with pi0 at its conservative bound 1;",
      "reorg plan v2 appendix B)"
    ),
    q_threshold = FDR_CUT,
    excluded_collections = c("C2:CGP (except radiation subset)", "C6", "C7"),
    radiation_curation = collections$curation,
    min_set_size = MIN_SET_SIZE,
    max_set_size = MAX_SET_SIZE,
    redundant_jaccard = REDUNDANT_JACCARD,
    gsea_param = 1, # hardcoded in the block ES (weight = |normal score|)
    msigdbr_version = as.character(packageVersion("msigdbr"))
  ),
  collections = names(gene_sets),
  units = units
)

out <- file.path(paths$processed, "thyr_enrichment_test.rds")
saveRDS(thyr_enrichment_test, out)
message("Saved: ", out, " (", length(units), " units)")
