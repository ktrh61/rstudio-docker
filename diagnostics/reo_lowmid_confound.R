# reo_lowmid_confound.R
# Ancillary QC for the REO intermediate-band application. This script provides
# two descriptive views of the band-score pattern in cases with relative-purity
# estimates:
#   (D) partial Spearman coefficient after rank-scale adjustment for purity;
#   (A) one-sided Mid-vs-Low Brunner-Munzel comparisons within two strata formed
#       at median purity.
# Neither analysis establishes an AS-band association independent of purity.
# This diagnostic runs outside the numbered stream and reads the output from
# reo_lowmid_purity.R.
# Input : diagnostics/output/reo_lowmid_purity.rds
#         lib/stat_brunnermunzel.R
# Output: diagnostics/output/reo_confound.rds

source("setup.R")
source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))

d <- readRDS(file.path(paths$root, "diagnostics", "output", "reo_lowmid_purity.rds"))
d <- d[d$band %in% c("R_Low", "R_Mid") & !is.na(d$score) & !is.na(d$tumor_purity), , drop = FALSE]
d$band_num <- ifelse(d$band == "R_Mid", 1L, 0L) # ordered: Mid higher AS
message("R_Low/Mid with score+purity: ",
  paste(names(table(d$band)), table(d$band), sep = "=", collapse = " "))

# --- (D) Partial Spearman after rank-scale adjustment for purity ------------
# The coefficient is the correlation between residuals from separate linear
# projections of band rank and score rank on purity rank.
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

# Descriptive permutation reference: permute band labels while retaining each
# score-purity pair. This breaks the associations of band with both variables;
# it is not conditional randomization at fixed purity.
set.seed(SEED)
perm <- replicate(9999, {
  b <- sample(d$band_num)
  partial_spearman(b, d$score, d$tumor_purity)
})
p_partial <- (sum(perm >= par_bs) + 1) / (length(perm) + 1) # one-sided (Mid higher)

message("\n(D) Rank correlations:")
message(sprintf("  band-score  %+.3f | band-purity %+.3f | score-purity %+.3f",
  raw_bs, raw_bp, raw_sp))
message(sprintf("  partial band-score after purity-rank adjustment = %+.3f ; descriptive permutation p(one-sided) = %.4f",
  par_bs, p_partial))
message("  (diagnostic only; not conditional inference at fixed purity)")

# --- (A) Descriptive Mid > Low comparisons within median-purity strata ------
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
    "  %-9s Low median %.1f (n=%d) | Mid median %.1f (n=%d) | Mid>Low p=%.4f, Pr(Low<Mid)=%.3f",
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
out_dir <- file.path(paths$root, "diagnostics", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_rds <- file.path(out_dir, "reo_confound.rds")
saveRDS(out, out_rds)
message("\nSaved: ", out_rds)
