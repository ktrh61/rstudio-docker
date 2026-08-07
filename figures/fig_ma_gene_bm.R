# 075_plot_ma_provisional.R  (PROVISIONAL / 仮置き)
# MA plots companion to the 085 volcano, one facet per analysis unit. The A axis
# is average abundance (mean log2 CPM across all samples in the unit) and the M
# axis is the log2 fold change High - Sporadic, computed from the 070 DEGES-
# normalized CPM for display only (the DE call is the rank-based Brunner-Munzel
# of 080, which has no fold change). Points are coloured by the 080 permutation
# FDR (fdr_perm < FDR_CUT). An MA plot also checks that the fold change does not
# trend with abundance -- a normalization sanity view. Provisional: not wired
# into the pipeline (it reads 080's output despite the lower number).
# Input : processed/thyr_normalized_counts.rds  (from 070; per-unit DGEList)
#         processed/thyr_expression_test.rds      (from 080; fdr_perm, effect)
#         processed/thyr_se_raw.rds               (gene_id -> gene_name)
# Output: output/ma_expression.png

source("setup.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
})

FDR_CUT <- 0.10
N_LABEL <- 8L # strongest genes to label per unit

norm <- readRDS(file.path(paths$processed, "thyr_normalized_counts.rds"))
test <- readRDS(file.path(paths$processed, "thyr_expression_test.rds"))
se <- readRDS(file.path(paths$processed, "thyr_se_raw.rds"))
name_of <- setNames(as.data.frame(rowData(se))$gene_name, rownames(se))

unit_order <- c("R_Tumor", "R_Normal", "B_Tumor", "B_Normal")
df <- do.call(rbind, lapply(unit_order, function(u) {
  dge <- norm$units[[u]]$dgelist
  grp <- as.character(dge$samples$group)
  hi <- grepl("High", grp, fixed = TRUE)
  sp <- grepl("Sporadic", grp, fixed = TRUE)
  lcpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  A <- rowMeans(lcpm)
  M <- rowMeans(lcpm[, hi, drop = FALSE]) - rowMeans(lcpm[, sp, drop = FALSE])

  g <- test$units[[u]]$genes
  fdr <- setNames(g$fdr_perm, g$gene_id)[rownames(dge)]
  p <- setNames(g$p_exact, g$gene_id)[rownames(dge)]
  data.frame(
    unit = u,
    gene = ifelse(is.na(name_of[rownames(dge)]) | name_of[rownames(dge)] == "",
      sub("\\..*$", "", rownames(dge)), unname(name_of[rownames(dge)])),
    A = A, M = M, fdr = unname(fdr), p_exact = unname(p),
    stringsAsFactors = FALSE
  )
}))
df$unit <- factor(df$unit, levels = unit_order)
df$sig <- factor(ifelse(df$fdr < FDR_CUT,
  ifelse(df$M >= 0, "up in High (FDR<0.10)", "up in Sporadic (FDR<0.10)"), "n.s."),
  levels = c("up in High (FDR<0.10)", "up in Sporadic (FDR<0.10)", "n.s."))

lab <- do.call(rbind, lapply(split(df, df$unit), function(d) head(d[order(d$p_exact), ], N_LABEL)))
pal <- c("up in High (FDR<0.10)" = "#eb6834",
  "up in Sporadic (FDR<0.10)" = "#2a78d6", "n.s." = "grey75")

p <- ggplot(df, aes(x = A, y = M)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_point(aes(colour = sig), size = 0.9, alpha = 0.6) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.6,
    max.overlaps = 20, min.segment.length = 0, colour = "grey20", seed = 1L) +
  scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
  facet_wrap(~unit, ncol = 2) +
  labs(
    x = expression("A:  average abundance  " * (log[2] * " CPM)")),
    y = expression("M:  " * log[2] * " fold change  (High - Sporadic)"),
    title = "Gene-level MA plot (070 CPM, 080 Brunner-Munzel FDR), per unit",
    subtitle = paste0("Fold change for display only (DE call is rank-based BM); ",
      "coloured by permutation FDR < ", FDR_CUT, ". R_Normal/B_Tumor negative controls.")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

if (!dir.exists(paths$output)) dir.create(paths$output, recursive = TRUE)
out_png <- file.path(paths$output, "ma_expression.png")
ggsave(out_png, p, width = 9, height = 8, dpi = 160, type = "cairo")
message("Saved: ", out_png)
for (u in unit_order) {
  d <- df[df$unit == u, ]
  message(sprintf("  %-9s |M| median %.3f | M range [%.2f, %.2f] | fdr<%.2f %d",
    u, stats::median(abs(d$M)), min(d$M), max(d$M), FDR_CUT, sum(d$fdr < FDR_CUT, na.rm = TRUE)))
}
