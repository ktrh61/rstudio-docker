# 280_plot_reo_evaluation_provisional.R  (PROVISIONAL / 仮置き)
# Figure for the REO out-of-sample validation, same concept as the v7
# 13_reo_evaluation_poc plot: REO reversal score (x) against assigned share
# (radiation attribution %, y), coloured by exposure band, training arms open /
# evaluation arms filled, with the AS band boundaries (33.3/66.6) and the A
# classification threshold marked. Shows the graded trend (band medians rise
# with assigned share). Provisional: not wired into the pipeline; regenerates
# output/reo_poc_vs_reversal.png from the saved panel/evaluation objects.
# Input : processed/thyr_reo_panel.rds        (from 110; panel + boundary + training)
#         processed/thyr_reo_evaluation.rds   (from 120; R_Low/Mid scores)
#         processed/thyr_se_raw.rds           (sample -> case for training AS)
#         processed/thyr_case_assigned_share.rds (assigned share)
# Output: output/reo_poc_vs_reversal.png

source("setup.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ggplot2)
})

reo <- readRDS(file.path(paths$processed, "thyr_reo_panel.rds"))
ev <- readRDS(file.path(paths$processed, "thyr_reo_evaluation.rds"))
se <- readRDS(file.path(paths$processed, "thyr_se_raw.rds"))
as_df <- as.data.frame(readRDS(file.path(paths$processed, "thyr_case_assigned_share.rds")))

# sample_id -> case -> assigned share (for the training arms)
cd <- as.data.frame(colData(se))
samp2case <- setNames(as.character(cd$case_submitter_id), as.character(cd$sample_submitter_id))
as_of <- setNames(as.numeric(as_df$assigned_share_approx), as.character(as_df$REBC_ID))
as_for_sample <- function(sid) unname(as_of[samp2case[sid]])

tr <- reo$training
spo <- data.frame(band = "R_Sporadic", score = as.integer(tr$r0_score),
  as = 0, stringsAsFactors = FALSE) # dose 0 by definition
hi_as <- as_for_sample(tr$r1_samples)
hi_as[is.na(hi_as)] <- 66.6
hi <- data.frame(band = "R_High", score = as.integer(tr$r1_score), as = hi_as)
mid <- data.frame(band = ev$samples$band, score = ev$samples$score,
  as = ev$samples$assigned_share, stringsAsFactors = FALSE)

d <- rbind(spo, hi, mid)
d$band <- factor(d$band, levels = c("R_Sporadic", "R_Low", "R_Mid", "R_High"))
d$set <- ifelse(d$band %in% c("R_Sporadic", "R_High"), "training", "evaluation")

thr <- reo$boundary$negative_max # A: positive if score > thr
# colorblind-safe, ordered by assigned share (blue -> aqua -> yellow -> orange)
pal <- c(R_Sporadic = "#2a78d6", R_Low = "#1baf7a", R_Mid = "#eda100", R_High = "#eb6834")
n_pairs <- nrow(reo$panel)

set.seed(1L)
p <- ggplot(d, aes(x = score, y = as)) +
  geom_hline(yintercept = c(33.3, 66.6), linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = thr + 0.5, linetype = "dashed", colour = "grey40") +
  annotate("text", x = thr + 0.5, y = 103, label = paste0("positive: score > ", thr),
    hjust = -0.03, vjust = 1, size = 3, colour = "grey30") +
  geom_jitter(aes(colour = band, shape = set), width = 0.15, height = 1.6,
    size = 2.6, alpha = 0.9, stroke = 0.8) +
  scale_colour_manual(values = pal, name = "Exposure band") +
  scale_shape_manual(values = c(training = 1L, evaluation = 16L), name = "Set") +
  scale_x_continuous(breaks = 0:n_pairs, limits = c(-0.5, n_pairs + 0.5)) +
  scale_y_continuous(breaks = c(0, 33.3, 66.6, 100), limits = c(-5, 105)) +
  labs(x = paste0("REO reversal score (panel of ", n_pairs, " pairs)"),
    y = "Assigned share  (radiation attribution, %)",
    title = "REO reversal grades with radiation attribution (RET/PTC tumours)",
    subtitle = "Out-of-sample R_Low / R_Mid (filled) vs training R_Sporadic / R_High (open); trend suggestive") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, colour = "grey30"),
    panel.grid.minor = element_blank(),
    legend.position = "right")

if (!dir.exists(paths$output)) dir.create(paths$output, recursive = TRUE)
out_png <- file.path(paths$output, "reo_poc_vs_reversal.png")
ggsave(out_png, p, width = 8.2, height = 5.6, dpi = 160, type = "cairo")
message("Saved: ", out_png)

for (b in levels(d$band)) {
  s <- d$score[d$band == b]
  message(sprintf("  %-11s n=%2d score median %.1f | AS median %.1f",
    b, length(s), stats::median(s), stats::median(d$as[d$band == b])))
}
