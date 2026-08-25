# fig_cohort_flow.R  (フロー図(仮))
# Cohort flow figure, both driver strata side by side (researcher decision
# 2026-08-14: per-driver attrition is carried by this figure, not by text).
# Formatting only -- draws the frozen step counts from 230; values must
# match N-08 on first run. Base graphics; layout is a draft for 判断点5.
# Extended 2026-08-15 (researcher Go): (i) step-1 box shows the total only
# (driver counts appear from the classification step onward), (ii) per-step
# exclusion annotations (differences of the frozen counts; REMARK-style),
# (iii) side branch at the band step for the REO evaluation set (counts from
# include_reo_evaluation in thyr_analysis_cohorts.rds; must match N-10).
# Input : processed/thyr_cohort_flow.rds, processed/thyr_analysis_cohorts.rds
# Output: output/figures/fig_cohort_flow.png (+ .tif submission copy, 600 dpi)

source("setup.R")

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
  purity_pass = "Relative purity >= 0.6"
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
  text(xc, y + ifelse(i == n, 0.32, 0.22), lab, cex = 0.95, font = 2)
  if (i == 1) {
    text(xc, y - 0.25, sprintf("total %d", flow$n_total[i]), cex = 0.9)
  } else if (i == n) {
    text(xc, y - 0.02, sprintf("total %d   (RET %d | BRAF %d)",
         flow$n_total[i], flow$n_RET[i], flow$n_BRAF[i]), cex = 0.9)
    text(xc, y - 0.36, sprintf(
      "RET: dose-zero %d | High-AS %d    BRAF: dose-zero %d | High-AS %d",
      grp["R_Sporadic"], grp["R_High"], grp["B_Sporadic"], grp["B_High"]),
      cex = 0.75)
  } else {
    text(xc, y - 0.25, sprintf("total %d   (RET %d | BRAF %d)",
         flow$n_total[i], flow$n_RET[i], flow$n_BRAF[i]), cex = 0.9)
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
           reasons[step_to]), cex = 0.75, adj = 0)
    } else if (step_to == "band_sporadic_or_high") {
      # two short lines so the branch arrow below stays clear
      text(xc + 0.3, y - 0.72, sprintf("excluded %d (RET %d | BRAF %d):",
           d_tot, d_ret, d_braf), cex = 0.75, adj = 0)
      text(xc + 0.3, y - 0.97, reasons[step_to], cex = 0.75, adj = 0)
    } else {
      text(xc + 0.3, y - 1.0, sprintf("excluded %d (RET %d | BRAF %d): %s",
           d_tot, d_ret, d_braf, reasons[step_to]), cex = 0.75, adj = 0)
    }
  }
}
# --- REO evaluation branch (leaves the flow at the band step) ---------------
y2 <- (n - 2) * 2 + 1 # y of the driver_classified box centre
n_lowmid_ret <- flow$n_RET[2] - flow$n_RET[3] # RET cases leaving at the band step
rect(9.0, y2 - 2.45, 13.3, y2 - 0.45)
text(11.15, y2 - 0.75, "REO application set", cex = 0.95, font = 2)
text(11.15, y2 - 1.12, sprintf("of %d RET Low-/Mid-AS cases,", n_lowmid_ret),
     cex = 0.85)
text(11.15, y2 - 1.44, sprintf("%d with tumor/normal pairs", nrow(ev)),
     cex = 0.85)
text(11.15, y2 - 1.76, sprintf("(Low-AS %d | Mid-AS %d)", n_low, n_mid),
     cex = 0.85)
text(11.15, y2 - 2.12, "unfiltered by outlier / purity screens", cex = 0.75)
arrows(xc, y2 - 1.25, 8.95, y2 - 1.25, length = 0.08)
par(op)
}

png(file.path(out_dir, "fig_cohort_flow.png"), type = "cairo",
    width = 2200, height = 230 * n, res = 200)
draw_flow()
dev.off()
cat("Saved:", file.path(out_dir, "fig_cohort_flow.png"), "\n")

# Submission copy (BJC artwork: minimum 300 dpi) -- 600 dpi LZW TIFF,
# same physical size as the PNG preview.
tiff(file.path(out_dir, "fig_cohort_flow.tif"), type = "cairo",
    width = 2200 / 200, height = 230 * n / 200, units = "in",
    res = 600, compression = "lzw")
draw_flow()
dev.off()
cat("Saved:", file.path(out_dir, "fig_cohort_flow.tif"), "\n")
