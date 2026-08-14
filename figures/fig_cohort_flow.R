# fig_cohort_flow.R  (フロー図(仮))
# Cohort flow figure, both driver strata side by side (researcher decision
# 2026-08-14: per-driver attrition is carried by this figure, not by text).
# Formatting only -- draws the frozen step counts from 230; values must
# match N-08 on first run. Base graphics; layout is a draft for 判断点5.
# Input : processed/thyr_cohort_flow.rds
# Output: output/figures/fig_cohort_flow.png

source("setup.R")

flow <- readRDS(file.path(paths$processed, "thyr_cohort_flow.rds"))
print(flow)
stopifnot(all(c("step", "n_total", "n_RET", "n_BRAF") %in% names(flow)))

labels <- c(
  all_cases = "All cases",
  driver_classified = "Driver classified (RET / BRAF)",
  band_sporadic_or_high = "Sporadic or High band",
  paired = "Tumor/normal pair available",
  pcod_clean = "Outlier screen passed (both tissues)",
  purity_pass = "Relative purity >= 0.6"
)
n <- nrow(flow)
out_dir <- file.path(paths$output, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
png(file.path(out_dir, "fig_cohort_flow.png"),
    width = 1600, height = 220 * n, res = 200)
op <- par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, n * 2), axes = FALSE,
     xlab = "", ylab = "")
for (i in seq_len(n)) {
  y <- (n - i) * 2 + 1
  lab <- labels[flow$step[i]]
  if (is.na(lab)) lab <- flow$step[i]
  rect(0.5, y - 0.55, 9.5, y + 0.55)
  text(5, y + 0.22, lab, cex = 0.95, font = 2)
  text(5, y - 0.25, sprintf("total %d   (RET %d | BRAF %d)",
       flow$n_total[i], flow$n_RET[i], flow$n_BRAF[i]), cex = 0.9)
  if (i < n) arrows(5, y - 0.58, 5, y - 1.42, length = 0.08)
}
par(op)
dev.off()
cat("Saved:", file.path(out_dir, "fig_cohort_flow.png"), "\n")
