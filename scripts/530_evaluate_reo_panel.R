# 530_evaluate_reo_panel.R
# Apply the finalized REO panel (520) to the intermediate-exposure RET tumours
# it was never trained on -- R_Low (assigned share 0-33.3%) and R_Mid (33.3-
# 66.6%) -- as a graded, out-of-sample check. If the panel captures a radiation-
# attributable signal, reversal scores should grade with assigned share:
# Sporadic (train R0) < R_Low < R_Mid < High (train R1). This is exploratory and
# descriptive; it does not alter the panel or its boundary.
# Input : processed/thyr_reo_panel.rds             (from 520; panel + boundary)
#         processed/thyr_analysis_cohorts.rds      (from 230; include_reo_evaluation)
#         processed/thyr_se_raw.rds                (from 120; stranded_second)
#         processed/gene_lengths.rds               (from 020)
#         lib/reo.R
# Output: processed/thyr_reo_evaluation.rds, output/reo_evaluation_samples.csv
#
# DESIGN NOTES (deferred, revisit): the R_Low/R_Mid tumours here are NOT filtered
# by ContamDE purity nor by PC-OD outlier detection; all paired RET tumours in
# the two AS bands are used. Whether to purity-match them to the training set is
# an open decision flagged for later.

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

source(file.path(paths$root, "lib", "reo.R"))
source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))

# AS bands (percent). High (>= 66.6) is the training arm; these are below it.
# AS band boundaries come from config.R (AS_LOW_MAX / AS_HIGH_MIN); the band
# assignment itself is fixed upstream in 140/230.

# --- Load inputs -----------------------------------------------------------
panel_path <- file.path(paths$processed, "thyr_reo_panel.rds")
cohorts_path <- file.path(paths$processed, "thyr_analysis_cohorts.rds")
se_path <- file.path(paths$processed, "thyr_se_raw.rds")
len_path <- file.path(paths$processed, "gene_lengths.rds")
for (p in c(panel_path, cohorts_path, se_path, len_path)) {
  if (!file.exists(p)) stop("missing input: ", p)
}
reo <- readRDS(panel_path)
panel <- reo$panel
boundary <- reo$boundary
dead_zone <- reo$config$dead_zone
cohorts <- readRDS(cohorts_path)
se <- readRDS(se_path)
gene_lengths <- readRDS(len_path)

# --- R_Low / R_Mid RET tumours (adoption fixed in 230) ---------------------
eval_cases <- cohorts[cohorts$include_reo_evaluation, , drop = FALSE]
as_tbl <- data.frame(
  case_submitter_id = eval_cases$case_submitter_id,
  dose_mgy = eval_cases$dose_mgy,
  assigned_share = eval_cases$assigned_share_approx,
  band = paste0("R_", eval_cases$band),
  tumor_id = eval_cases$tumor_id,
  stringsAsFactors = FALSE
)
message("Intermediate RET tumours: ",
  paste(names(table(as_tbl$band)), table(as_tbl$band), sep = "=", collapse = " "))

# --- Apply the panel -------------------------------------------------------
log2_tpm <- reo_log2_tpm(se, gene_lengths, as_tbl$tumor_id)
reversal_score <- function(l2tpm, samples, panel, dead_zone) {
  sc <- integer(length(samples))
  names(sc) <- samples
  for (k in seq_len(nrow(panel))) {
    r <- l2tpm[panel$up[k], samples] - l2tpm[panel$down[k], samples]
    sc <- sc + as.integer(abs(r) >= dead_zone & sign(r) != panel$r0_sign[k])
  }
  sc
}
as_tbl$score <- reversal_score(log2_tpm, as_tbl$tumor_id, panel, dead_zone)

# --- Read A: classification vs the R0-based threshold -----------------------
classify <- function(score, b) ifelse(score > b$negative_max, "positive", "negative")
as_tbl$class <- classify(as_tbl$score, boundary)

# --- Read B: graded increase Low -> Mid (out-of-sample, permutation BM) -----
# Two ordered bands, so the ordered-alternative (Jonckheere-Terpstra) test is a
# one-sided Brunner-Munzel: is the R_Mid reversal score stochastically greater
# than R_Low? BM is used (not Wilcoxon) for the same reason as the main
# analysis: the mixture makes the arms unequally dispersed. Training arms
# (Sporadic/High) are NOT in this test -- they are shown only for the figure.
low_score <- as_tbl$score[as_tbl$band == "R_Low"]
mid_score <- as_tbl$score[as_tbl$band == "R_Mid"]
bm_low_mid <- brunnermunzel_mc_test(
  low_score, mid_score, alternative = "less", method = "auto", seed = SEED
)
message(sprintf(
  "\nRead B (out-of-sample): Mid > Low reversal score, one-sided BM p = %.4f (%s), effect P(Low<Mid)=%.3f",
  bm_low_mid$p.value, attr(bm_low_mid, "mc")$method, unname(bm_low_mid$estimate)
))

# --- Report: graded validation ---------------------------------------------
n_pairs <- nrow(panel)
message(sprintf("\nPanel size %d ; R0-based positive threshold: score > %d",
  n_pairs, boundary$negative_max))
message("Reversal score by exposure band (training arms shown for reference):")
tr <- reo$training
ref <- data.frame(
  band = c("R_Sporadic(train)", "R_High(train)"),
  n = c(length(tr$r0_score), length(tr$r1_score)),
  median = c(stats::median(tr$r0_score), stats::median(tr$r1_score)),
  min = c(min(tr$r0_score), min(tr$r1_score)),
  max = c(max(tr$r0_score), max(tr$r1_score))
)
obs <- do.call(rbind, lapply(c("R_Low", "R_Mid"), function(b) {
  s <- as_tbl$score[as_tbl$band == b]
  data.frame(band = b, n = length(s), median = stats::median(s),
    min = min(s), max = max(s))
}))
summary_tbl <- rbind(ref[1, ], obs, ref[2, ])
print(summary_tbl, row.names = FALSE)
message("\nClassification by band:")
print(table(band = as_tbl$band, class = as_tbl$class))

# --- Assemble and save -----------------------------------------------------
thyr_reo_evaluation <- list(
  date = Sys.Date(),
  config = list(as_low_max = AS_LOW_MAX, as_mid_max = AS_HIGH_MIN,
    panel_size = n_pairs, note = "R_Low/R_Mid unfiltered for purity/outliers"),
  samples = as_tbl[, c("case_submitter_id", "band", "assigned_share",
    "dose_mgy", "tumor_id", "score", "class")],
  summary = summary_tbl,
  # Read A: R0-based threshold classification. Read B: out-of-sample ordered
  # test on Low vs Mid only (training arms excluded). Sporadic/High are for the
  # graded figure, not the inference.
  read_A_threshold = boundary$negative_max,
  read_B = list(
    test = "one-sided Brunner-Munzel, Mid > Low reversal score",
    method = attr(bm_low_mid, "mc")$method,
    p_value = bm_low_mid$p.value,
    effect_P_low_lt_mid = unname(bm_low_mid$estimate)
  )
)
out_rds <- file.path(paths$processed, "thyr_reo_evaluation.rds")
saveRDS(thyr_reo_evaluation, out_rds)
message("Saved: ", out_rds)

if (dir.exists(paths$output)) {
  out_csv <- file.path(paths$output, "reo_evaluation_samples.csv")
  utils::write.csv(thyr_reo_evaluation$samples, out_csv, row.names = FALSE)
  message("Saved: ", out_csv)
}
