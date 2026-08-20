# 530_evaluate_reo_panel.R
# Apply the finalized REO panel (520), without refitting, to the intermediate-
# exposure RET tumours not used in its construction: R_Low (assigned share
# 0-33.3%) and R_Mid (33.3-66.6%). The graded scores are examined against the
# hypothesized band ordering. This application is descriptive, does not alter
# the panel or its cutoff, and is not independent external validation.
# Input : processed/thyr_reo_panel.rds             (from 520; panel + boundary)
#         processed/thyr_analysis_cohorts.rds      (from 230; include_reo_evaluation)
#         processed/thyr_se_raw.rds                (from 120; single count assay)
#         processed/gene_lengths.rds               (from 020)
#         lib/reo.R
# Output: processed/thyr_reo_evaluation.rds, output/reo_evaluation_samples.csv
#
# The R_Low/R_Mid tumours are not filtered by ContamDE purity or PC-OD. These
# quantities are examined separately as non-exclusion diagnostics and do not
# alter the intermediate-band application set.

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "reo.R"))
source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))

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
  assigned_share = eval_cases$assigned_share,
  band = paste0("R_", eval_cases$band),
  tumor_id = eval_cases$tumor_id,
  stringsAsFactors = FALSE
)
message("Intermediate RET tumours: ",
  paste(names(table(as_tbl$band)), table(as_tbl$band), sep = "=", collapse = " "))

# --- Apply the panel -------------------------------------------------------
log2_tpm <- reo_log2_tpm(se, gene_lengths, as_tbl$tumor_id)
as_tbl$score <- reversal_score(log2_tpm, as_tbl$tumor_id, panel, dead_zone)

# --- Read A: classification vs the R0-based threshold -----------------------
as_tbl$class <- classify_reversal(as_tbl$score, boundary)

# --- Read B: graded Mid vs Low application (permutation BM) -----------------
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
  "\nRead B (intermediate-band application): Mid > Low reversal score, one-sided BM p = %.4f (%s), effect Pr(Low<Mid)=%.3f",
  bm_low_mid$p.value, attr(bm_low_mid, "mc")$method, unname(bm_low_mid$estimate)
))

# --- Report: graded band summary -------------------------------------------
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
  # Read A: R0-based threshold classification. Read B: intermediate-band
  # ordered test on Low vs Mid only (construction arms excluded).
  # Sporadic/High are shown in the graded figure but do not enter this test.
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
