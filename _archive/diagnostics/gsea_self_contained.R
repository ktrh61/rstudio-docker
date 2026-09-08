# gsea_self_contained.R
# Set-level self-contained diagnostic (reorg plan v2 appendix B.6; reading
# rules fixed there before this first ran). The competitive test in 420 asks
# whether signal CONCENTRATES inside a set relative to outside; this
# diagnostic asks, per set, whether signal EXISTS among its members at all
# (Goeman & Buehlmann's competitive / self-contained distinction). Under the
# diffuse-signal hypothesis the expected picture is broad weak firing here
# with no competitive concentration -- a positive measurement of infiltration
# rather than an absence claim. Reporting both is not a contradiction: the
# two layers answer different questions (same non-contradiction form as
# reorg plan v2 s3.1b).
#
# Harm-free form: DISTRIBUTION-ONLY reporting. Per unit x family, the share
# of sets with self-contained permutation p below 0.05 / 0.01 is reported
# against the uniform expectation. Individual set p-values and names are
# deliberately not reported and carry no claims: under diffuse infiltration
# the sets are nearly interchangeable and naming would mislead. A side line
# records DEG coverage (share of q < 0.10 genes inside the tested-set union
# vs the same share for all genes) to separate "infiltrates the tested space"
# from "lives outside it". Descriptive only; never promoted to a claim.
#
# Statistic: the 410 omnibus count restricted to the set -- the number of
# member genes with |BM statistic| at or beyond the unit's pooled-null
# 95% / 99% quantile -- with the permutation p taken over the same label
# shuffles 410 saved (perm_index, by reference).
# Input : processed/thyr_normalized_counts.rds  (from 310)
#         processed/thyr_expression_test.rds    (from 410; perm_index, genes)
#         lib/stat_brunnermunzel.R, lib/gsea_permutation.R,
#         lib/gsea_collections.R
# Output: diagnostics/output/gsea_self_contained.rds

source("setup.R")

suppressPackageStartupMessages({
  library(edgeR)
  library(msigdbr)
  library(parallel)
  library(Matrix)
})

source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))
source(file.path(paths$root, "lib", "gsea_permutation.R"))
source(file.path(paths$root, "lib", "gsea_collections.R"))
source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "annotation.R"))

# --- Configuration ---------------------------------------------------------
SPECIES <- "Homo sapiens"
COUNT_QUANTILES <- c(0.95, 0.99) # pooled-null cuts for the per-set count
P_BARS <- c(0.05, 0.01) # reported ECDF points of the per-set p

pin_blas_threads()

# --- Load inputs -----------------------------------------------------------
norm_path <- file.path(paths$processed, "thyr_normalized_counts.rds")
test_path <- file.path(paths$processed, "thyr_expression_test.rds")
if (!file.exists(norm_path)) stop("thyr_normalized_counts.rds not found (310)")
if (!file.exists(test_path)) stop("thyr_expression_test.rds not found (410)")
normalized <- readRDS(norm_path)
expression_test <- readRDS(test_path)

collections <- load_gene_set_collections(SPECIES)
gene_sets <- collections$gene_sets

out_dir <- file.path(paths$root, "diagnostics", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Per-unit measurement --------------------------------------------------
measure_unit <- function(dgelist, unit) {
  arms <- unit_arms(dgelist$samples$group, unit)
  cpm_matrix <- edgeR::cpm(
    dgelist,
    normalized.lib.sizes = TRUE, prior.count = 0, log = FALSE
  )[, c(arms$sporadic, arms$high), drop = FALSE]
  storage.mode(cpm_matrix) <- "double"
  ensembl <- strip_ensembl_version(rownames(cpm_matrix))
  keep_row <- !duplicated(ensembl)
  cpm_matrix <- cpm_matrix[keep_row, , drop = FALSE]
  rownames(cpm_matrix) <- ensembl[keep_row]
  nx <- length(arms$sporadic)
  n_genes <- nrow(cpm_matrix)

  perm_index <- expression_test$units[[unit]]$perm_index
  n_perm <- ncol(perm_index)

  index <- lapply(gene_sets, gsea_pathway_index,
    gene_ids = rownames(cpm_matrix),
    min_size = GSEA_MIN_SET_SIZE, max_size = GSEA_MAX_SET_SIZE
  )
  collection_of <- rep(names(index), lengths(index))
  flat <- unlist(unname(index), recursive = FALSE)

  statistic <- abs(brunnermunzel_statistics(cpm_matrix, nx))
  null_statistic <- gsea_bind_null_columns(
    mclapply(seq_len(n_perm), function(i) {
      abs(brunnermunzel_statistics(
        cpm_matrix[, perm_index[, i], drop = FALSE], nx
      ))
    }, mc.cores = WORKERS),
    nrow(cpm_matrix)
  )

  # Sparse membership matrix: sets x genes.
  membership <- sparseMatrix(
    i = rep(seq_along(flat), lengths(flat)),
    j = unlist(flat),
    x = 1,
    dims = c(length(flat), n_genes)
  )

  # Same cut convention as 410's omnibus (order-statistic pick on the sorted
  # pooled null), so "the 410 count restricted to the set" is exact, not
  # approximate-by-quantile-interpolation.
  pooled <- sort(as.numeric(null_statistic))
  rows <- lapply(COUNT_QUANTILES, function(qq) {
    cut <- pooled[max(1L, round(length(pooled) * qq))]
    observed_counts <- as.numeric(membership %*% (statistic >= cut))
    null_counts <- as.matrix(membership %*% (null_statistic >= cut))
    p_set <- (rowSums(null_counts >= observed_counts) + 1) / (n_perm + 1)
    do.call(rbind, lapply(names(index), function(collection) {
      keep <- collection_of == collection
      data.frame(
        unit = unit, collection = collection, m_sets = sum(keep),
        count_quantile = qq,
        share_p_lt_05 = mean(p_set[keep] < P_BARS[1]),
        share_p_lt_01 = mean(p_set[keep] < P_BARS[2]),
        expected_05 = P_BARS[1], expected_01 = P_BARS[2],
        median_p = stats::median(p_set[keep]),
        stringsAsFactors = FALSE
      )
    }))
  })

  # DEG coverage side line: q < 0.10 genes inside the tested-set union.
  genes <- expression_test$units[[unit]]$genes
  deg_ids <- strip_ensembl_version(genes$gene_id[genes$q_storey < FDR_CUT])
  union_members <- rownames(cpm_matrix)[
    unique(unlist(flat))
  ]
  coverage <- data.frame(
    unit = unit,
    n_deg = length(deg_ids),
    deg_in_sets = if (length(deg_ids)) {
      mean(deg_ids %in% union_members)
    } else {
      NA_real_
    },
    background_in_sets = mean(rownames(cpm_matrix) %in% union_members),
    stringsAsFactors = FALSE
  )

  list(shares = do.call(rbind, rows), coverage = coverage)
}

message("Self-contained set diagnostic (distribution-only reporting):")
results <- lapply(names(normalized$units), function(unit) {
  message("  Unit ", unit)
  measure_unit(normalized$units[[unit]]$dgelist, unit)
})
names(results) <- names(normalized$units)

shares <- do.call(rbind, lapply(results, `[[`, "shares"))
coverage <- do.call(rbind, lapply(results, `[[`, "coverage"))
rownames(shares) <- rownames(coverage) <- NULL
print(shares, digits = 3, row.names = FALSE)
message("DEG coverage by the tested-set union (q < ", FDR_CUT, "):")
print(coverage, digits = 3, row.names = FALSE)

# --- Save ------------------------------------------------------------------
gsea_self_contained <- list(
  date = Sys.Date(),
  config = list(
    statistic = "per-set count of |BM| >= pooled-null quantile",
    count_quantiles = COUNT_QUANTILES,
    permutation = "410's saved perm_index, by reference",
    n_perm = expression_test$config$n_perm,
    reporting = "distribution-only; no per-set values or names (appendix B.6)",
    min_set_size = GSEA_MIN_SET_SIZE,
    max_set_size = GSEA_MAX_SET_SIZE,
    q_threshold = FDR_CUT
  ),
  shares = shares,
  coverage = coverage
)
out <- file.path(out_dir, "gsea_self_contained.rds")
saveRDS(gsea_self_contained, out)
message("Saved: ", out)
