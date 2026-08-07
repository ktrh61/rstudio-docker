# fig_gene_bm_evidence.R  (PROVISIONAL / 仮置き)
# Volcano plots of the 410 gene-level Brunner-Munzel results, one facet per
# analysis unit. The x axis is the signed BM effect 2*effect - 1 = P(X<Y) -
# P(X>Y) (Cliff's delta; > 0 = higher in the High arm), since a rank test has no
# fold change. The y axis is -log10(exact permutation p). Points are coloured by
# the permutation-calibrated FDR (fdr_perm < FDR_CUT), and the strongest genes
# per unit are labelled. Provisional: not wired into the pipeline.
# Input : processed/thyr_expression_test.rds  (from 410)
#         processed/thyr_se_raw.rds            (gene_id -> gene_name)
# Output: output/volcano_expression.png

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
    fdr = g$fdr_perm,
    p_exact = g$p_exact,
    stringsAsFactors = FALSE
  )
}))
df$unit <- factor(df$unit, levels = unit_order)
lab_up <- sprintf("up in High (FDR<%.2f)", FDR_CUT)
lab_down <- sprintf("up in Sporadic (FDR<%.2f)", FDR_CUT)
df$sig <- ifelse(df$fdr < FDR_CUT,
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
  geom_point(aes(colour = sig), size = 1.0, alpha = 0.7) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.6,
    max.overlaps = 20, min.segment.length = 0, colour = "grey20", seed = 1L) +
  scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
  facet_wrap(~unit, ncol = 2) +
  labs(
    x = expression("signed Brunner-Munzel effect  " * (P(X < Y) - P(X > Y))
      * "   " %->% "  higher in High"),
    y = expression(-log[10] * "(exact permutation p)"),
    title = "Gene-level Brunner-Munzel volcano (410), per analysis unit",
    subtitle = paste0("Rank-based effect (no fold change); coloured by permutation FDR < ",
      FDR_CUT, ". R_Normal/B_Tumor are negative controls.")
  ) +
  theme_thyr()

save_figure(p, "volcano_expression.png", width = 9, height = 8)
for (u in unit_order) {
  g <- test$units[[u]]$genes
  message(sprintf("  %-9s genes %5d | min p_exact %.2e | fdr<%.2f %d",
    u, nrow(g), min(g$p_exact), FDR_CUT, sum(g$fdr_perm < FDR_CUT)))
}
