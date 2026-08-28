# fig_gene_bm_evidence.R
# Gene-level Brunner-Munzel evidence plot, one facet per analysis unit (not a
# volcano: a rank test has no fold change). The x axis is the signed BM effect
# 2*theta - 1 = Pr(X<Y) - Pr(X>Y) (Cliff's delta; > 0 = higher in the High
# arm) and the y axis is -log10(exact permutation p). Points are coloured by
# the Storey q inference (q_storey < FDR_CUT) and the strongest genes per unit
# are labelled. Status is tracked in figures/manifest.csv.
# Input : processed/thyr_expression_test.rds  (from 410)
#         processed/thyr_se_raw.rds            (gene_id -> gene_name)
# Output: output/figures/fig_gene_bm_evidence.png (+ .tif 600 dpi, .pdf vector)
# Drawn at final width 175 mm (BJC double column), text 5.5-7 pt, no in-figure
# title (the legend carries it) -- artwork-guide alignment 2026-08-28.

source("setup.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ggplot2)
  library(ggrepel)
})

source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "annotation.R"))
source(file.path(paths$root, "lib", "plot_theme.R"))

# FDR_CUT comes from config.R via setup.R.
N_LABEL <- 8L # strongest genes to label per unit

test <- readRDS(file.path(paths$processed, "thyr_expression_test.rds"))
se <- readRDS(file.path(paths$processed, "thyr_se_raw.rds"))
name_of <- gene_name_map(se)

unit_order <- UNIT_ORDER
df <- do.call(rbind, lapply(unit_order, function(u) {
  g <- test$units[[u]]$genes
  data.frame(
    unit = u,
    gene = gene_label(g$gene_id, name_of),
    x = 2 * g$effect - 1,
    y = -log10(pmax(g$p_exact, .Machine$double.xmin)),
    q = g$q_storey,
    p_exact = g$p_exact,
    stringsAsFactors = FALSE
  )
}))
df$unit <- factor(df$unit, levels = unit_order)
lab_up <- sprintf("higher in High-AS (q<%.2f)", FDR_CUT)
lab_down <- sprintf("lower in High-AS (q<%.2f)", FDR_CUT)
df$sig <- ifelse(df$q < FDR_CUT,
  ifelse(df$x >= 0, lab_up, lab_down),
  "n.s.")
df$sig <- factor(df$sig, levels = c(lab_up, lab_down, "n.s."))

# strongest genes per unit (by p) to label
lab <- do.call(rbind, lapply(split(df, df$unit), function(d) {
  head(d[order(d$p_exact), ], N_LABEL)
}))

pal <- setNames(c(COL_UP, COL_DOWN, COL_NS), c(lab_up, lab_down, "n.s."))

p <- ggplot(df, aes(x = x, y = y)) +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(colour = sig), size = 0.6, alpha = 0.7) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), size = label_size(),
    family = FONT_FAMILY, segment.size = 0.2,
    max.overlaps = 20, min.segment.length = 0, colour = "grey20", seed = 1L) +
  scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
  facet_wrap(~unit, ncol = 2) +
  labs(
    x = "signed Brunner–Munzel effect  2θ−1 = Pr(X<Y)−Pr(X>Y)  → higher in High-AS",
    y = expression(-log[10] * "(exact permutation p)")
  ) +
  theme_thyr()

save_figure(p, "fig_gene_bm_evidence.png", width = 175, height = 150)
for (u in unit_order) {
  g <- test$units[[u]]$genes
  message(sprintf("  %-9s genes %5d | min p_exact %.2e | q<%.2f %d",
    u, nrow(g), min(g$p_exact), FDR_CUT, sum(g$q_storey < FDR_CUT)))
}
