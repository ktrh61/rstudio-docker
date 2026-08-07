# reo_lowmid_confound.R  (PROVISIONAL / 仮置き)
# Layer-2 QC for the REO out-of-sample validation: reo_lowmid_purity.R showed the reversal score
# correlates with tumour purity within R_Low/Mid, and purity itself grades with
# assigned share (AS). So the graded Low<Mid score increase could be AS-driven
# (radiation) OR purity-driven. This script separates them:
#   (D) Partial Spearman: association of band (Low/Mid) with score, controlling
#       for purity. If it survives, the gradient is not merely purity.
#   (A) Purity-stratified ordered test: within each purity stratum, is Mid > Low
#       score (one-sided BM)? Purity does not vary within a stratum, so a
#       surviving Mid>Low is AS-attributable.
# Provisional: not wired into the pipeline; reads reo_lowmid_purity.R's provisional purity.
# Input : processed/thyr_reo_lowmid_purity_provisional.rds (from 260)
#         lib/stat_brunnermunzel.R
# Output: processed/thyr_reo_confound_provisional.rds

source("setup.R")
source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))

d <- readRDS(file.path(paths$processed, "thyr_reo_lowmid_purity_provisional.rds"))
d <- d[d$band %in% c("R_Low", "R_Mid") & !is.na(d$score) & !is.na(d$tumor_purity), , drop = FALSE]
d$band_num <- ifelse(d$band == "R_Mid", 1L, 0L) # ordered: Mid higher AS
message("R_Low/Mid with score+purity: ",
  paste(names(table(d$band)), table(d$band), sep = "=", collapse = " "))

# --- (D) Partial Spearman correlation, controlling for purity ---------------
# Partial rho of (band, score) given purity = correlation of the residuals of
# each on purity, using rank (Spearman) variables.
partial_spearman <- function(x, y, z) {
  rx <- rank(x); ry <- rank(y); rz <- rank(z)
  ex <- residuals(lm(rx ~ rz))
  ey <- residuals(lm(ry ~ rz))
  unname(cor(ex, ey))
}
raw_bs <- cor(d$band_num, d$score, method = "spearman")
raw_bp <- cor(d$band_num, d$tumor_purity, method = "spearman")
raw_sp <- cor(d$score, d$tumor_purity, method = "spearman")
par_bs <- partial_spearman(d$band_num, d$score, d$tumor_purity)

# Permutation p for the partial correlation: permute band labels, keep score &
# purity paired, so the null is "band unrelated to score given purity".
set.seed(SEED)
perm <- replicate(9999, {
  b <- sample(d$band_num)
  partial_spearman(b, d$score, d$tumor_purity)
})
p_partial <- (sum(perm >= par_bs) + 1) / (length(perm) + 1) # one-sided (Mid higher)

message("\n(D) Rank correlations:")
message(sprintf("  band-score  %+.3f | band-purity %+.3f | score-purity %+.3f",
  raw_bs, raw_bp, raw_sp))
message(sprintf("  PARTIAL band-score | purity = %+.3f ; permutation p(one-sided) = %.4f",
  par_bs, p_partial))
message("  (survives = gradient not merely purity ; collapses = purity-confounded)")

# --- (A) Purity-stratified ordered test (Mid > Low within stratum) ----------
message("\n(A) Purity-stratified Mid > Low reversal score (one-sided BM):")
cut <- stats::median(d$tumor_purity)
d$stratum <- ifelse(d$tumor_purity >= cut, "hi_purity", "lo_purity")
strata <- lapply(c("lo_purity", "hi_purity"), function(s) {
  ds <- d[d$stratum == s, ]
  lo <- ds$score[ds$band == "R_Low"]; mi <- ds$score[ds$band == "R_Mid"]
  if (length(lo) < 2 || length(mi) < 2) {
    message(sprintf("  %-9s n(Low=%d,Mid=%d): too few for a test", s, length(lo), length(mi)))
    return(NULL)
  }
  bm <- brunnermunzel_mc_test(lo, mi, alternative = "less", method = "auto", seed = SEED)
  message(sprintf(
    "  %-9s Low median %.1f (n=%d) | Mid median %.1f (n=%d) | Mid>Low p=%.4f, P(Low<Mid)=%.3f",
    s, stats::median(lo), length(lo), stats::median(mi), length(mi),
    bm$p.value, unname(bm$estimate)))
  data.frame(stratum = s, low_median = stats::median(lo), mid_median = stats::median(mi),
    p = bm$p.value, effect = unname(bm$estimate))
})
strata_tbl <- do.call(rbind, strata)

# --- Save ------------------------------------------------------------------
out <- list(
  partial = list(band_score = raw_bs, band_purity = raw_bp, score_purity = raw_sp,
    partial_band_score_given_purity = par_bs, p_partial_one_sided = p_partial),
  strata = strata_tbl,
  data = d[, c("case_submitter_id", "band", "assigned_share", "tumor_purity", "score", "stratum")]
)
out_rds <- file.path(paths$processed, "thyr_reo_confound_provisional.rds")
saveRDS(out, out_rds)
message("\nSaved: ", out_rds)
