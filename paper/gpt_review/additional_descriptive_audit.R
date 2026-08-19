#!/usr/bin/env Rscript

# Descriptive audit for the R_Tumor result.
#
# HISTORICAL NON-PRODUCTION SCRIPT. The selection-conditioned effect summaries
# and stats::prcomp PCA produced here were withdrawn from the manuscript. The
# PCA is not an HDLSS-specific method, is not PC-OD, and is not approved for
# production use. Retain this script only to document what was run.
#
# This script reads the finalized pipeline artifacts but does not rerun QC,
# normalization, differential-expression testing, or multiple-testing
# adjustment. It creates no new p-values and no alternative DEG list. Its
# original purpose was to add magnitude and covariate context to the existing
# result; that use was subsequently withdrawn.
#
# Run in the canonical container with the repository mounted at /workspace
# and a writable output directory mounted at /review:
#   Rscript /workspace/paper/gpt_review/additional_descriptive_audit.R /review

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1L]] else "/review"
if (!dir.exists(out_dir)) stop("Output directory does not exist: ", out_dir)

required <- c(
  "processed/thyr_expression_test.rds",
  "processed/thyr_normalized_counts.rds",
  "processed/thyr_analysis_cohorts.rds",
  "processed/thyr_clinical.rds",
  "processed/thyr_case_outliers.rds"
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing input(s): ", paste(missing, collapse = ", "))

test <- readRDS(required[[1L]])
norm <- readRDS(required[[2L]])
cohort <- readRDS(required[[3L]])
clinical <- readRDS(required[[4L]])
outliers <- readRDS(required[[5L]])

unit <- "R_Tumor"
genes <- test$units[[unit]]$genes
dge <- norm$units[[unit]]$dgelist
if (is.null(genes) || is.null(dge)) stop("R_Tumor input is absent")

group <- as.character(dge$samples$group)
sporadic <- which(group == "R_Sporadic")
high <- which(group == "R_High")
if (!length(sporadic) || !length(high)) stop("R_Tumor groups are incomplete")

lcpm <- edgeR::cpm(
  dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1
)
mean_log2_fc <- rowMeans(lcpm[, high, drop = FALSE]) -
  rowMeans(lcpm[, sporadic, drop = FALSE])

gene_audit <- genes
gene_audit$signed_rank_effect <- 2 * gene_audit$effect - 1
gene_audit$abs_signed_rank_effect <- abs(gene_audit$signed_rank_effect)
gene_audit$mean_log2_fc <- unname(mean_log2_fc[gene_audit$gene_id])
gene_audit$abs_mean_log2_fc <- abs(gene_audit$mean_log2_fc)
gene_audit$existing_q_lt_0_10 <- gene_audit$q_storey < 0.10
gene_audit$direction <- ifelse(
  gene_audit$signed_rank_effect > 0,
  "higher_in_High",
  ifelse(gene_audit$signed_rank_effect < 0, "lower_in_High", "tie")
)

annotation <- dge$genes
if (!is.null(annotation)) {
  ann <- as.data.frame(annotation, stringsAsFactors = FALSE)
  ann$gene_id <- rownames(dge)
  label_col <- intersect(c("gene_name", "gene_symbol", "symbol"), names(ann))
  if (length(label_col)) {
    label <- setNames(as.character(ann[[label_col[[1L]]]]), ann$gene_id)
    gene_audit$gene_name <- unname(label[gene_audit$gene_id])
  }
}
if (!"gene_name" %in% names(gene_audit)) gene_audit$gene_name <- NA_character_

qstats <- function(x) {
  x <- x[is.finite(x)]
  c(
    minimum = min(x),
    q25 = unname(stats::quantile(x, 0.25)),
    median = stats::median(x),
    q75 = unname(stats::quantile(x, 0.75)),
    maximum = max(x)
  )
}

summarize_gene_subset <- function(d, subset_name) {
  measures <- list(
    abs_signed_rank_effect = d$abs_signed_rank_effect,
    signed_rank_effect = d$signed_rank_effect,
    abs_mean_log2_fc = d$abs_mean_log2_fc,
    mean_log2_fc = d$mean_log2_fc
  )
  do.call(rbind, lapply(names(measures), function(measure) {
    s <- qstats(measures[[measure]])
    data.frame(
      subset = subset_name,
      n_genes = nrow(d),
      measure = measure,
      minimum = s[["minimum"]],
      q25 = s[["q25"]],
      median = s[["median"]],
      q75 = s[["q75"]],
      maximum = s[["maximum"]],
      stringsAsFactors = FALSE
    )
  }))
}

selected <- gene_audit[gene_audit$existing_q_lt_0_10, , drop = FALSE]
effect_summary <- rbind(
  summarize_gene_subset(gene_audit, "all_tested"),
  summarize_gene_subset(selected, "existing_q_lt_0.10"),
  summarize_gene_subset(
    selected[selected$direction == "higher_in_High", , drop = FALSE],
    "existing_q_lt_0.10_higher_in_High"
  ),
  summarize_gene_subset(
    selected[selected$direction == "lower_in_High", , drop = FALSE],
    "existing_q_lt_0.10_lower_in_High"
  )
)

top_genes <- gene_audit[order(gene_audit$p_exact, gene_audit$gene_id), , drop = FALSE]
top_genes <- head(top_genes, 20L)
top_genes <- top_genes[, c(
  "gene_id", "gene_name", "effect", "signed_rank_effect", "mean_log2_fc",
  "p_exact", "q_storey", "rank", "existing_q_lt_0_10", "direction"
)]

# Map tumor sample IDs back to cases and clinical covariates while preserving
# the DGEList column order used by PCA.
main_ret <- cohort[cohort$include_main_bm & cohort$driver == "RET", , drop = FALSE]
case_index <- match(colnames(dge), main_ret$tumor_id)
if (anyNA(case_index)) stop("Could not map all R_Tumor samples to cases")
sample_audit <- main_ret[case_index, c(
  "case_submitter_id", "tumor_id", "band", "tumor_purity"
)]
clinical_index <- match(sample_audit$case_submitter_id, clinical$REBC_ID)
if (anyNA(clinical_index)) stop("Could not map all R_Tumor cases to clinical data")
sample_audit$sex <- tolower(as.character(clinical$SEX[clinical_index]))
sample_audit$age_surgery <- as.numeric(clinical$AGE_SURGERY[clinical_index])
sample_audit$age_exposure <- as.numeric(clinical$AGE_EXPOSURE[clinical_index])
sample_audit$ret_fusion_partner <- as.character(
  clinical$Designated_Driver[clinical_index]
)
sample_audit$group <- factor(
  paste0("R_", sample_audit$band), levels = c("R_Sporadic", "R_High")
)

# PCA is unsupervised: all 15,621 already-tested genes enter, with no
# q-value-based or group-dependent feature selection. Gene-wise centering is
# performed by prcomp; genes are not variance-scaled.
pca <- stats::prcomp(t(lcpm), center = TRUE, scale. = FALSE)
variance <- pca$sdev^2 / sum(pca$sdev^2)
sample_audit$PC1 <- pca$x[, 1L]
sample_audit$PC2 <- pca$x[, 2L]
sample_audit$PC3 <- pca$x[, 3L]

fmt_range <- function(x, digits = 3L) {
  x <- x[is.finite(x)]
  sprintf(
    paste0("%.", digits, "f [%.", digits, "f-%.", digits, "f]"),
    stats::median(x), min(x), max(x)
  )
}

group_rows <- lapply(split(sample_audit, sample_audit$group), function(d) {
  partner <- sort(table(d$ret_fusion_partner), decreasing = TRUE)
  data.frame(
    group = as.character(d$group[[1L]]),
    n = nrow(d),
    female = sum(d$sex == "female"),
    male = sum(d$sex == "male"),
    age_surgery_median_range = fmt_range(d$age_surgery, 1L),
    relative_purity_median_range = fmt_range(d$tumor_purity, 3L),
    ret_fusion_partner_counts = paste(
      sprintf("%s %d", names(partner), as.integer(partner)), collapse = "; "
    ),
    PC1_median_range = fmt_range(d$PC1, 3L),
    PC2_median_range = fmt_range(d$PC2, 3L),
    stringsAsFactors = FALSE
  )
})
covariate_summary <- do.call(rbind, group_rows)
rownames(covariate_summary) <- NULL

rho_row <- function(pc_name, covariate_name) {
  x <- sample_audit[[pc_name]]
  y <- sample_audit[[covariate_name]]
  ok <- is.finite(x) & is.finite(y)
  data.frame(
    component = pc_name,
    covariate = covariate_name,
    n_complete = sum(ok),
    spearman_rho = unname(stats::cor(x[ok], y[ok], method = "spearman")),
    stringsAsFactors = FALSE
  )
}
pca_covariate_summary <- rbind(
  rho_row("PC1", "age_surgery"),
  rho_row("PC1", "tumor_purity"),
  rho_row("PC2", "age_surgery"),
  rho_row("PC2", "tumor_purity")
)

pca_variance <- data.frame(
  component = paste0("PC", seq_along(variance)),
  proportion = variance,
  cumulative = cumsum(variance),
  stringsAsFactors = FALSE
)

# Feasibility only: assess whether a prespecified covariate-adjusted linear
# design would be full rank. No adjusted model is fitted to any gene here.
design_data <- data.frame(
  group = stats::relevel(factor(sample_audit$group), ref = "R_Sporadic"),
  age_z = as.numeric(scale(sample_audit$age_surgery)),
  purity_z = as.numeric(scale(sample_audit$tumor_purity)),
  sex = factor(sample_audit$sex),
  partner = factor(sample_audit$ret_fusion_partner)
)
complete_design <- stats::complete.cases(design_data)
design_matrix <- stats::model.matrix(
  ~ group + age_z + purity_z + sex + partner,
  data = design_data[complete_design, , drop = FALSE]
)
design_feasibility <- data.frame(
  n_total = nrow(design_data),
  n_complete = sum(complete_design),
  n_columns = ncol(design_matrix),
  rank = qr(design_matrix)$rank,
  residual_df = nrow(design_matrix) - qr(design_matrix)$rank,
  full_column_rank = qr(design_matrix)$rank == ncol(design_matrix),
  condition_number = kappa(design_matrix),
  stringsAsFactors = FALSE
)

pcod_counts <- stats::aggregate(
  cbind(has_outlier_tumor, has_outlier_normal) ~ group,
  data = outliers,
  FUN = sum
)
pcod_counts <- pcod_counts[order(pcod_counts$group), , drop = FALSE]

utils::write.csv(
  effect_summary,
  file.path(out_dir, "additional_effect_size_summary.csv"), row.names = FALSE
)
utils::write.csv(
  top_genes,
  file.path(out_dir, "additional_top_gene_effects.csv"), row.names = FALSE
)
utils::write.csv(
  sample_audit,
  file.path(out_dir, "additional_sample_pca_covariates.csv"), row.names = FALSE
)
utils::write.csv(
  covariate_summary,
  file.path(out_dir, "additional_covariate_summary.csv"), row.names = FALSE
)
utils::write.csv(
  pca_covariate_summary,
  file.path(out_dir, "additional_pca_covariate_summary.csv"), row.names = FALSE
)
utils::write.csv(
  pca_variance,
  file.path(out_dir, "additional_pca_variance.csv"), row.names = FALSE
)
utils::write.csv(
  design_feasibility,
  file.path(out_dir, "additional_adjusted_design_feasibility.csv"),
  row.names = FALSE
)
utils::write.csv(
  pcod_counts,
  file.path(out_dir, "additional_pcod_flag_counts.csv"), row.names = FALSE
)

png(
  file.path(out_dir, "additional_rtumor_pca.png"),
  width = 2400, height = 2000, res = 200, type = "cairo"
)
par(mfrow = c(2, 2), mar = c(4.3, 4.3, 2.5, 1.0), cex = 0.9)
xlab <- sprintf("PC1 (%.1f%%)", 100 * variance[[1L]])
ylab <- sprintf("PC2 (%.1f%%)", 100 * variance[[2L]])
plot_panel <- function(col, pch, title) {
  plot(
    sample_audit$PC1, sample_audit$PC2,
    col = col, pch = pch, cex = 1.25,
    xlab = xlab, ylab = ylab, main = title
  )
}

group_cols <- ifelse(sample_audit$group == "R_High", "#D55E00", "#0072B2")
sex_pch <- ifelse(sample_audit$sex == "female", 16, 17)
plot_panel(group_cols, sex_pch, "Group and sex")
legend(
  "topright",
  legend = c("R_Sporadic", "R_High", "female", "male"),
  col = c("#0072B2", "#D55E00", "grey30", "grey30"),
  pch = c(15, 15, 16, 17), bty = "n", cex = 0.85
)

continuous_cols <- grDevices::colorRampPalette(
  c("#440154", "#21908C", "#FDE725")
)(100L)
map_continuous <- function(x) {
  if (diff(range(x)) == 0) return(rep(continuous_cols[[50L]], length(x)))
  index <- 1L + round(99 * (x - min(x)) / diff(range(x)))
  continuous_cols[index]
}

purity_cols <- map_continuous(sample_audit$tumor_purity)
plot_panel(purity_cols, 16, "Relative tumor-purity score")
purity_ticks <- seq(
  min(sample_audit$tumor_purity), max(sample_audit$tumor_purity), length.out = 3L
)
legend(
  "topright", legend = sprintf("%.2f", purity_ticks),
  col = map_continuous(purity_ticks), pch = 16, bty = "n", cex = 0.85
)

age_cols <- map_continuous(sample_audit$age_surgery)
plot_panel(age_cols, 16, "Age at surgery (years)")
age_ticks <- seq(
  min(sample_audit$age_surgery), max(sample_audit$age_surgery), length.out = 3L
)
legend(
  "topright", legend = sprintf("%.0f", age_ticks),
  col = map_continuous(age_ticks), pch = 16, bty = "n", cex = 0.85
)

partner_levels <- sort(unique(sample_audit$ret_fusion_partner))
partner_palette <- setNames(c("#009E73", "#CC79A7", "#E69F00"), partner_levels)
plot_panel(
  unname(partner_palette[sample_audit$ret_fusion_partner]), 16,
  "RET fusion partner"
)
legend(
  "topright", legend = partner_levels,
  col = unname(partner_palette[partner_levels]), pch = 16,
  bty = "n", cex = 0.85
)
dev.off()

sel_abs_effect <- qstats(selected$abs_signed_rank_effect)
sel_abs_lfc <- qstats(selected$abs_mean_log2_fc)
up <- selected[selected$direction == "higher_in_High", , drop = FALSE]
down <- selected[selected$direction == "lower_in_High", , drop = FALSE]

summary_lines <- c(
  "# Additional descriptive audit of R_Tumor",
  "",
  "This audit reads finalized analysis artifacts. It does not rerun QC, normalization, hypothesis testing, or multiple-testing adjustment. It creates no new p-values and no alternative gene list.",
  "",
  sprintf("- Existing q<0.10 list: %d genes (%d higher and %d lower in High).", nrow(selected), nrow(up), nrow(down)),
  sprintf(
    "- Among those genes, median absolute signed Brunner–Munzel effect was %.3f (IQR %.3f–%.3f).",
    sel_abs_effect[["median"]], sel_abs_effect[["q25"]], sel_abs_effect[["q75"]]
  ),
  sprintf(
    "- Their median absolute display-only mean log2 fold change was %.3f (IQR %.3f–%.3f).",
    sel_abs_lfc[["median"]], sel_abs_lfc[["q25"]], sel_abs_lfc[["q75"]]
  ),
  sprintf(
    "- PCA used all %d tested genes; PC1 and PC2 explained %.1f%% and %.1f%% of variance.",
    nrow(gene_audit), 100 * variance[[1L]], 100 * variance[[2L]]
  ),
  sprintf(
    "- In the finalized main cohort, median relative purity was %.3f in R_Sporadic and %.3f in R_High; these values are from the main-cohort analysis object and are distinct from the later all-band REO purity diagnostic.",
    stats::median(sample_audit$tumor_purity[sample_audit$group == "R_Sporadic"]),
    stats::median(sample_audit$tumor_purity[sample_audit$group == "R_High"])
  ),
  sprintf(
    "- A candidate adjusted design with group, age at surgery, relative purity, sex, and RET fusion partner used %d complete cases, had %d columns, rank %d, and %d residual degrees of freedom. This establishes numerical estimability only; no adjusted gene model was fitted.",
    design_feasibility$n_complete, design_feasibility$n_columns,
    design_feasibility$rank, design_feasibility$residual_df
  ),
  sprintf(
    "- The existing PC-OD output flagged %d RET tumor and %d RET normal samples; the only main-cohort target flag was one B_High tumor. Thus PC-OD did not alter the realized R_Tumor sample set.",
    sum(pcod_counts$has_outlier_tumor[grepl("^R_", pcod_counts$group)]),
    sum(pcod_counts$has_outlier_normal[grepl("^R_", pcod_counts$group)])
  ),
  "",
  "See the CSV files for complete summaries and sample-level values."
)
writeLines(summary_lines, file.path(out_dir, "additional_descriptive_audit.md"))

message("Wrote descriptive audit to: ", normalizePath(out_dir))
