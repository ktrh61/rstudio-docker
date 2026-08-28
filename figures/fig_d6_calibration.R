# fig_d6_calibration.R
# Held-out null calibration of the set-level inference, as a forest plot:
# one row per contrast x collection cell (16), x = share of null replicates
# with at least one discovery at q_bh < 0.10, with exact binomial 95% CI,
# against the nominal 0.10 level. Formatting only -- draws the frozen
# calibration summary from the D6 run; values must match N-24 (range),
# N-25 (B_Normal/H excess) and N-26 (boundary cells) on first run.
# The single cell whose CI sits above the nominal level (B_Normal/Hallmark)
# is emphasised; it is the disclosed calibration excess.
# Input : diagnostics/output/gsea_null_calibration.rds
# Output: output/figures/fig_d6_calibration.png (+ .tif 600 dpi, .pdf vector)
# Drawn at final width 175 mm, text 5.6-7 pt, no in-figure title/subtitle;
# key and axis wording follow the manuscript (Clopper-Pearson interval,
# Benjamini-Hochberg q) -- artwork-guide alignment 2026-08-28.

source("setup.R")
suppressPackageStartupMessages({
  library(ggplot2)
})

source(file.path(paths$root, "lib", "plot_theme.R"))

cal <- readRDS(file.path(paths$root, "diagnostics", "output",
                         "gsea_null_calibration.rds"))
s <- cal$summary
stopifnot(nrow(s) == 16L, cal$config$nominal == 0.1)

contrast_order <- c("R_Tumor", "R_Normal", "B_Tumor", "B_Normal")
coll_label <- c(H = "Hallmark", `C2:CP` = "C2:CP", `C5:GO:BP` = "C5:GO:BP",
                `C2:CGP:radiation` = "radiation (C2:CGP)")
s$contrast <- factor(s$unit, levels = contrast_order)
s$coll <- coll_label[s$collection]
s$row <- factor(paste(s$contrast, s$coll, sep = " · "),
                levels = rev(paste(rep(contrast_order, each = 4),
                                   coll_label[c("H", "C2:CP", "C5:GO:BP",
                                                "C2:CGP:radiation")],
                                   sep = " · ")))
# The disclosed excess: the cell whose CI lower bound exceeds the nominal level.
s$excess <- s$ci_lo > cal$config$nominal
lab_ok <- "interval overlaps nominal"
lab_ex <- "interval entirely above nominal"
s$flag <- factor(ifelse(s$excess, lab_ex, lab_ok), levels = c(lab_ex, lab_ok))

p <- ggplot(s, aes(x = p_any, y = row)) +
  geom_vline(xintercept = cal$config$nominal, linetype = "dashed",
             colour = "grey40") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi, colour = flag),
                orientation = "y", width = 0.25, linewidth = 0.4) +
  geom_point(aes(colour = flag), size = 1.4) +
  scale_colour_manual(values = setNames(c(COL_UP, "grey30"),
                                        c(lab_ex, lab_ok)), name = NULL) +
  scale_x_continuous(limits = c(0, 0.30), breaks = seq(0, 0.3, 0.05)) +
  labs(
    x = "proportion of held-out pseudo-observations with ≥1 discovery (Benjamini–Hochberg q<0.10)",
    y = NULL
  ) +
  theme_thyr()

save_figure(p, "fig_d6_calibration.png", width = 175, height = 120)

# First-run verification against the frozen ledger values.
message(sprintf("  p_any range %.2f-%.2f (N-24 expects 0.01-0.18)",
                min(s$p_any), max(s$p_any)))
ex <- s[s$excess, ]
message(sprintf("  excess cells: %d | %s %s p_any %.2f CI %.3f-%.3f (N-25)",
                nrow(ex), ex$unit, ex$collection, ex$p_any, ex$ci_lo, ex$ci_hi))
for (r in which(s$collection == "C2:CGP:radiation" & s$unit %in% c("B_Tumor", "B_Normal"))) {
  message(sprintf("  boundary: %s radiation p_any %.2f CI %.3f-%.3f (N-26)",
                  s$unit[r], s$p_any[r], s$ci_lo[r], s$ci_hi[r]))
}
