# Shared gene-filter helpers. Two deliberately different recipes exist in the
# pipeline and are preserved as-is (they are NOT interchangeable):
#   - filter_pc_expr_mask(): protein-coding & filterByExpr computed on the
#     UNSUBSET count matrix (library sizes from all genes) -- the 220/ContamDE
#     recipe, also used by the purity diagnostic.
#   - 310 keeps its own inline subset-first recipe (protein-coding first, then
#     filterByExpr on the reduced DGEList) -- numerically different via the
#     library-size denominators; single consumer, left in place.
#   - logcpm_for_qc(): the PC-OD input recipe (no protein-coding filter,
#     filterByExpr without group, unnormalized log-CPM, prior.count = 2).

# Mask over ALL genes: protein-coding AND expressed (filterByExpr grouped by
# Normal/Tumor) for the given paired sample set.
filter_pc_expr_mask <- function(counts, protein_coding, normal_ids, tumor_ids) {
  cg <- counts[, c(normal_ids, tumor_ids), drop = FALSE]
  grp <- factor(c(
    rep("Normal", length(normal_ids)),
    rep("Tumor", length(tumor_ids))
  ))
  protein_coding & edgeR::filterByExpr(cg, group = grp)
}

# Unnormalized log-CPM for sample-outlier detection (PC-OD input). Raw library
# sizes preserve composition outliers that normalization would mask.
# keep_lib_sizes is a caller decision: 210 uses TRUE, the Low/Mid diagnostic
# has historically used FALSE (changes the CPM denominator after subsetting);
# the divergence is preserved and its unification is a recorded open item.
logcpm_for_qc <- function(counts, ids, keep_lib_sizes) {
  y <- edgeR::DGEList(counts = counts[, ids, drop = FALSE])
  keep <- edgeR::filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = keep_lib_sizes]
  edgeR::cpm(y, log = TRUE, normalized.lib.sizes = FALSE, prior.count = 2)
}
