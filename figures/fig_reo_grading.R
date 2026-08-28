# fig_reo_grading.R
# Figure for the REO out-of-sample evaluation, same concept as the v7
# 13_reo_evaluation_poc plot: REO reversal score (x) against assigned share
# (assigned share %, y), coloured by exposure band, training groups open /
# evaluation arms filled, with the AS band boundaries (33.3/66.6) and the A
# classification threshold marked. Shows the graded pattern (band medians rise
# with assigned share) at the descriptive-observation level fixed in the plan
# (v2 s0.5). Status is tracked in figures/manifest.csv.
# Input : processed/thyr_reo_panel.rds        (from 520; panel + boundary + training)
#         processed/thyr_reo_evaluation.rds   (from 530; R_Low/Mid scores)
#         processed/thyr_se_raw.rds           (sample -> case for training AS)
#         processed/thyr_case_assigned_share.rds (assigned share)
# Output: output/figures/fig_reo_grading.png (+ .tif 600 dpi, .pdf vector)
# Drawn at final width 175 mm, text 5.5-7 pt, no in-figure title/subtitle --
# artwork-guide alignment 2026-08-28 (the retired subtitle's 'out-of-sample'
# wording also conflicted with the manuscript's non-validation stance).

source("setup.R")
source(file.path(paths$root, "lib", "plot_theme.R"))
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
as_of <- setNames(as.numeric(as_df$assigned_share), as.character(as_df$REBC_ID))
as_for_sample <- function(sid) unname(as_of[samp2case[sid]])

tr <- reo$training

# R_High: assigned share must resolve for every training sample -- no silent
# imputation (the earlier `is.na -> 66.6` fallback is retired; it never fired).
hi_as <- as_for_sample(tr$r1_samples)
stopifnot(!anyNA(hi_as))
hi <- data.frame(band = "R_High", score = as.integer(tr$r1_score), y = hi_as,
  stringsAsFactors = FALSE)
mid <- data.frame(band = ev$samples$band, score = ev$samples$score,
  y = ev$samples$assigned_share, stringsAsFactors = FALSE)

# R_Sporadic: AS is UNDEFINED for unexposed cases (not zero). The band is drawn
# in a separate strip below the AS axis; within the strip, identical scores are
# stacked, so vertical position is a case counter, not a value.
STRIP_BASE <- -10
STRIP_STEP <- -3.4
spo_scores <- as.integer(tr$r0_score)
spo <- do.call(rbind, lapply(split(spo_scores, spo_scores), function(v)
  data.frame(band = "R_Sporadic", score = v,
    y = STRIP_BASE + (seq_along(v) - 1) * STRIP_STEP, stringsAsFactors = FALSE)))

d <- rbind(spo, hi, mid)
# display labels follow the manuscript nomenclature (2026-08-21):
# analysis-file bands R_Sporadic/R_Low/R_Mid/R_High -> dose-zero/Low-AS/Mid-AS/High-AS
band_labels <- c(R_Sporadic = "dose-zero", R_Low = "Low-AS",
                 R_Mid = "Mid-AS", R_High = "High-AS")
d$set <- ifelse(d$band %in% c("R_Sporadic", "R_High"), "construction", "application")
d$band <- factor(unname(band_labels[d$band]), levels = unname(band_labels))

thr <- reo$boundary$negative_max # A: positive if score > thr
pal <- setNames(PAL_BANDS[names(band_labels)], unname(band_labels)) # lib/plot_theme.R
n_pairs <- nrow(reo$panel)
y_min <- min(d$y) - 4

# No jitter anywhere: scores are integers (exact) and AS separates the exposed
# points vertically on its own; the Sporadic strip stacks instead of jittering.
p <- ggplot(d, aes(x = score, y = y)) +
  geom_hline(yintercept = c(33.3, 66.6), linetype = "dashed", colour = "grey70") +
  geom_hline(yintercept = -5, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = thr + 0.5, linetype = "dashed", colour = "grey40") +
  annotate("text", x = thr + 0.5, y = 103, label = paste0("positive: score > ", thr),
    hjust = -0.03, vjust = 1, size = label_size(), family = FONT_FAMILY,
    colour = "grey30") +
  annotate("text", x = n_pairs + 0.5, y = STRIP_BASE, hjust = 1, vjust = 1,
    size = label_size(), family = FONT_FAMILY, colour = "grey30",
    label = "dose-zero strip: AS undefined (unexposed);\nstacked points count cases") +
  geom_point(aes(colour = band, shape = set), size = 1.8, alpha = 0.9, stroke = 0.6) +
  scale_colour_manual(values = pal, name = "AS band") +
  scale_shape_manual(values = c(construction = 1L, application = 16L), name = "Set") +
  scale_x_continuous(breaks = 0:n_pairs, limits = c(-0.5, n_pairs + 0.5)) +
  scale_y_continuous(breaks = c(0, 33.3, 66.6, 100), limits = c(y_min, 105)) +
  labs(x = paste0("REO reversal score (panel of ", n_pairs, " pairs)"),
    y = "Assigned share  (radiation attributability, %)") +
  theme_thyr(legend_position = "right")

save_figure(p, "fig_reo_grading.png", width = 175, height = 120)

for (b in levels(d$band)) {
  s <- d$score[d$band == b]
  as_msg <- if (b == "dose-zero") "AS undefined (strip)" else
    sprintf("AS median %.1f", stats::median(d$y[d$band == b]))
  message(sprintf("  %-11s n=%2d score median %.1f | %s", b, length(s),
    stats::median(s), as_msg))
}
