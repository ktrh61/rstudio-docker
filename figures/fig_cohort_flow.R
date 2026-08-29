# fig_cohort_flow.R  (フロー図(仮))
# Cohort flow figure, both driver strata side by side (researcher decision
# 2026-08-14: per-driver attrition is carried by this figure, not by text).
# Formatting only -- draws the frozen step counts from 230; values must
# match N-08 on first run. Base graphics; layout is a draft for 判断点5.
# Box wording fixed 2026-08-29: the application set requires a tumor sample only
# (230: driver RET & Low/Mid & has_tumor_reo); 31 of the 36 also have a normal.
# Extended 2026-08-15 (researcher Go): (i) step-1 box shows the total only
# (driver counts appear from the classification step onward), (ii) per-step
# exclusion annotations (differences of the frozen counts; REMARK-style),
# (iii) side branch at the band step for the REO evaluation set (counts from
# include_reo_evaluation in thyr_analysis_cohorts.rds; must match N-10).
# Input : processed/thyr_cohort_flow.rds, processed/thyr_analysis_cohorts.rds
# Output: output/figures/fig_cohort_flow.png (+ .tif 600 dpi, .pdf vector)
# Drawn at final width 175 mm with 7 pt base text (box titles 7 pt, counts
# 6.5 pt, annotations 5.5 pt; Liberation Sans) -- artwork-guide alignment
# 2026-08-28.

source("setup.R")
source(file.path(paths$root, "lib", "plot_theme.R")) # FONT_FAMILY

flow <- readRDS(file.path(paths$processed, "thyr_cohort_flow.rds"))
print(flow)
stopifnot(all(c("step", "n_total", "n_RET", "n_BRAF") %in% names(flow)))

co <- readRDS(file.path(paths$processed, "thyr_analysis_cohorts.rds"))
mb <- co[co$include_main_bm, ]
grp <- table(paste0(ifelse(mb$driver == "RET", "R", "B"), "_", mb$band))
cat("main BM groups:", paste(names(grp), grp, collapse = " | "), "\n")
ev <- co[co$include_reo_evaluation, ]
n_low <- sum(ev$band == "Low")
n_mid <- sum(ev$band == "Mid")
cat(sprintf("REO evaluation set: %d (Low %d | Mid %d)\n",
            nrow(ev), n_low, n_mid))

labels <- c(
  all_cases = "All cases",
  driver_classified = "Driver classified (RET / BRAF)",
  band_sporadic_or_high = "Dose-zero or High-AS band",
  paired = "Tumor/normal pair available",
  pcod_clean = "Outlier screen passed (both tissues)",
  purity_pass = "Relative purity ≥ 0.6"
)
reasons <- c(
  driver_classified = "no single classified driver",
  band_sporadic_or_high = "Low-AS / Mid-AS band or no reference",
  paired = "no tumor/normal pair",
  pcod_clean = "outlier-flagged tissue",
  purity_pass = "relative purity < 0.6"
)

n <- nrow(flow)
out_dir <- file.path(paths$output, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
draw_flow <- function() {
op <- par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(0, 0, type = "n", xlim = c(0, 13.5), ylim = c(0, n * 2), axes = FALSE,
     xlab = "", ylab = "")
xc <- 4.0 # main-column centre
for (i in seq_len(n)) {
  y <- (n - i) * 2 + 1
  lab <- labels[flow$step[i]]
  if (is.na(lab)) lab <- flow$step[i]
  rect(0.5, y - 0.55, 7.5, y + 0.55)
  text(xc, y + ifelse(i == n, 0.32, 0.22), lab, cex = 1.0, font = 2)
  if (i == 1) {
    text(xc, y - 0.25, sprintf("total %d", flow$n_total[i]), cex = 0.93)
  } else if (i == n) {
    text(xc, y - 0.02, sprintf("total %d   (RET %d | BRAF %d)",
         flow$n_total[i], flow$n_RET[i], flow$n_BRAF[i]), cex = 0.93)
    text(xc, y - 0.36, sprintf(
      "RET: dose-zero %d | High-AS %d    BRAF: dose-zero %d | High-AS %d",
      grp["R_Sporadic"], grp["R_High"], grp["B_Sporadic"], grp["B_High"]),
      cex = 0.79)
  } else {
    text(xc, y - 0.25, sprintf("total %d   (RET %d | BRAF %d)",
         flow$n_total[i], flow$n_RET[i], flow$n_BRAF[i]), cex = 0.93)
  }
  if (i < n) {
    arrows(xc, y - 0.58, xc, y - 1.42, length = 0.08)
    # exclusion annotation (REMARK-style), derived from the frozen counts
    d_tot <- flow$n_total[i] - flow$n_total[i + 1]
    d_ret <- flow$n_RET[i] - flow$n_RET[i + 1]
    d_braf <- flow$n_BRAF[i] - flow$n_BRAF[i + 1]
    step_to <- flow$step[i + 1]
    if (step_to == "driver_classified") {
      text(xc + 0.3, y - 1.0, sprintf("excluded %d: %s", d_tot,
           reasons[step_to]), cex = 0.79, adj = 0)
    } else if (step_to == "band_sporadic_or_high") {
      # two short lines so the branch arrow below stays clear
      text(xc + 0.3, y - 0.72, sprintf("excluded %d (RET %d | BRAF %d):",
           d_tot, d_ret, d_braf), cex = 0.79, adj = 0)
      text(xc + 0.3, y - 0.97, reasons[step_to], cex = 0.79, adj = 0)
    } else {
      text(xc + 0.3, y - 1.0, sprintf("excluded %d (RET %d | BRAF %d): %s",
           d_tot, d_ret, d_braf, reasons[step_to]), cex = 0.79, adj = 0)
    }
  }
}
# --- REO evaluation branch (leaves the flow at the band step) ---------------
y2 <- (n - 2) * 2 + 1 # y of the driver_classified box centre
n_lowmid_ret <- flow$n_RET[2] - flow$n_RET[3] # RET cases leaving at the band step
rect(9.0, y2 - 2.45, 13.3, y2 - 0.45)
text(11.15, y2 - 0.75, "REO application set", cex = 1.0, font = 2)
text(11.15, y2 - 1.12, sprintf("of %d RET Low-/Mid-AS cases,", n_lowmid_ret),
     cex = 0.86)
text(11.15, y2 - 1.44, sprintf("%d with a tumor RNA-seq sample", nrow(ev)),
     cex = 0.86)
text(11.15, y2 - 1.76, sprintf("(Low-AS %d | Mid-AS %d)", n_low, n_mid),
     cex = 0.86)
text(11.15, y2 - 2.12, "unfiltered by outlier / purity screens", cex = 0.79)
arrows(xc, y2 - 1.25, 8.95, y2 - 1.25, length = 0.08)
par(op)
}

W_IN <- 175 / 25.4 # final width: BJC double column
H_IN <- W_IN * (230 * n) / 2200 # aspect of the original layout
PT <- 7 # base point size at final width (artwork guide: text 5-7 pt)
png(file.path(out_dir, "fig_cohort_flow.png"), type = "cairo",
    width = W_IN, height = H_IN, units = "in", res = 160,
    pointsize = PT, family = FONT_FAMILY)
draw_flow()
dev.off()
cat("Saved:", file.path(out_dir, "fig_cohort_flow.png"), "\n")

# Submission copies: 600 dpi LZW TIFF and vector PDF (fonts embedded).
tiff(file.path(out_dir, "fig_cohort_flow.tif"), type = "cairo",
    width = W_IN, height = H_IN, units = "in", res = 600,
    pointsize = PT, family = FONT_FAMILY, compression = "lzw")
draw_flow()
dev.off()
cat("Saved:", file.path(out_dir, "fig_cohort_flow.tif"), "\n")

grDevices::cairo_pdf(file.path(out_dir, "fig_cohort_flow.pdf"),
    width = W_IN, height = H_IN, pointsize = PT, family = FONT_FAMILY)
draw_flow()
dev.off()
cat("Saved:", file.path(out_dir, "fig_cohort_flow.pdf"), "\n")
