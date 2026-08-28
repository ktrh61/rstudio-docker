# fig_ma_gene_bm.R
# MA plots companion to fig_gene_bm_evidence.R, one facet per analysis unit. The A axis
# is average abundance (mean log2 CPM across all samples in the unit) and the M
# axis is the log2 fold change High - Sporadic, computed from the 310 DEGES-
# normalized CPM for display only (the DE call is the rank-based Brunner-Munzel
# of 410, which has no fold change). Points are coloured by the 410 Storey q
# inference (q_storey < FDR_CUT). An MA plot also checks that the fold change
# does not trend with abundance -- a normalization sanity view. Status is
# tracked in figures/manifest.csv.
# Input : processed/thyr_normalized_counts.rds  (from 310; per-unit DGEList)
#         processed/thyr_expression_test.rds      (from 410; q_storey, effect)
#         processed/thyr_se_raw.rds               (gene_id -> gene_name)
# Output: output/figures/fig_ma_gene_bm.png (+ .tif 600 dpi, .pdf vector)
# Drawn at final width 175 mm, text 5.5-7 pt, no in-figure title -- artwork-
# guide alignment 2026-08-28. Colour-key wording follows the manuscript bands
# (High-AS), not the code-era 'exposed'.

source("setup.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
})

source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "annotation.R"))
source(file.path(paths$root, "lib", "plot_theme.R"))

# FDR_CUT comes from config.R via setup.R.
N_LABEL <- 8L # strongest genes to label per unit

norm <- readRDS(file.path(paths$processed, "thyr_normalized_counts.rds"))
test <- readRDS(file.path(paths$processed, "thyr_expression_test.rds"))
se <- readRDS(file.path(paths$processed, "thyr_se_raw.rds"))
name_of <- gene_name_map(se)

unit_order <- UNIT_ORDER
df <- do.call(rbind, lapply(unit_order, function(u) {
  dge <- norm$units[[u]]$dgelist
  arms <- unit_arms(dge$samples$group, u)
  hi <- arms$high
  sp <- arms$sporadic
  lcpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  A <- rowMeans(lcpm)
  M <- rowMeans(lcpm[, hi, drop = FALSE]) - rowMeans(lcpm[, sp, drop = FALSE])

  g <- test$units[[u]]$genes
  q <- setNames(g$q_storey, g$gene_id)[rownames(dge)]
  p <- setNames(g$p_exact, g$gene_id)[rownames(dge)]
  data.frame(
    unit = u,
    gene = gene_label(rownames(dge), name_of),
    A = A, M = M, q = unname(q), p_exact = unname(p),
    stringsAsFactors = FALSE
  )
}))
df$unit <- factor(df$unit, levels = unit_order)
lab_up <- sprintf("higher in High-AS (q<%.2f)", FDR_CUT)
lab_down <- sprintf("lower in High-AS (q<%.2f)", FDR_CUT)
df$sig <- factor(ifelse(df$q < FDR_CUT,
  ifelse(df$M >= 0, lab_up, lab_down), "n.s."),
  levels = c(lab_up, lab_down, "n.s."))

lab <- do.call(rbind, lapply(split(df, df$unit), function(d) head(d[order(d$p_exact), ], N_LABEL)))
pal <- setNames(c(COL_UP, COL_DOWN, COL_NS), c(lab_up, lab_down, "n.s."))

p <- ggplot(df, aes(x = A, y = M)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_point(aes(colour = sig), size = 0.6, alpha = 0.6) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), size = label_size(),
    family = FONT_FAMILY, segment.size = 0.2,
    max.overlaps = 20, min.segment.length = 0, colour = "grey20", seed = 1L) +
  scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
  facet_wrap(~unit, ncol = 2) +
  labs(
    x = expression("A:  average abundance  " * (log[2] * " CPM)")),
    y = expression("M:  " * log[2] * " fold change  (High - Sporadic)")
  ) +
  theme_thyr()

save_figure(p, "fig_ma_gene_bm.png", width = 175, height = 150)
for (u in unit_order) {
  d <- df[df$unit == u, ]
  message(sprintf("  %-9s |M| median %.3f | M range [%.2f, %.2f] | q<%.2f %d",
    u, stats::median(abs(d$M)), min(d$M), max(d$M), FDR_CUT, sum(d$q < FDR_CUT, na.rm = TRUE)))
}
