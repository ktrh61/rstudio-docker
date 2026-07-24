# 120_evaluate_reo_panel.R
# Apply the finalized REO panel (110) to the intermediate-exposure RET tumours
# it was never trained on -- R_Low (assigned share 0-33.3%) and R_Mid (33.3-
# 66.6%) -- as a graded, out-of-sample check. If the panel captures a radiation-
# attributable signal, reversal scores should grade with assigned share:
# Sporadic (train R0) < R_Low < R_Mid < High (train R1). This is exploratory and
# descriptive; it does not alter the panel or its boundary.
# Input : processed/thyr_reo_panel.rds             (from 110; panel + boundary)
#         processed/thyr_clinical.rds              (from 001; driver)
#         processed/thyr_case_assigned_share.rds   (from 041; dose, AS)
#         processed/thyr_se_raw.rds                (from 021; stranded_second)
#         processed/gene_lengths.rds               (from 005)
#         utils/reo_tpm.R
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

source(file.path(paths$root, "utils", "reo_tpm.R"))

# AS bands (percent). High (>= 66.6) is the training arm; these are below it.
AS_LOW_MAX <- 33.3
AS_MID_MAX <- 66.6

# --- Load inputs -----------------------------------------------------------
panel_path <- file.path(paths$processed, "thyr_reo_panel.rds")
clin_path <- file.path(paths$processed, "thyr_clinical.rds")
as_path <- file.path(paths$processed, "thyr_case_assigned_share.rds")
se_path <- file.path(paths$processed, "thyr_se_raw.rds")
len_path <- file.path(paths$processed, "gene_lengths.rds")
for (p in c(panel_path, clin_path, as_path, se_path, len_path)) {
  if (!file.exists(p)) stop("missing input: ", p)
}
reo <- readRDS(panel_path)
panel <- reo$panel
boundary <- reo$boundary
dead_zone <- reo$config$dead_zone
clinical <- readRDS(clin_path)
assigned_share <- as.data.frame(readRDS(as_path))
se <- readRDS(se_path)
gene_lengths <- readRDS(len_path)

# --- Identify R_Low / R_Mid RET tumours ------------------------------------
ret_values <- c("CCDC6-RET", "NCOA4-RET", "RET-OTHER")
ret_cases <- as.character(clinical$REBC_ID)[clinical$Designated_Driver %in% ret_values]

as_tbl <- data.frame(
  case_submitter_id = as.character(assigned_share$REBC_ID),
  dose_mgy = as.numeric(assigned_share$dose_mgy),
  assigned_share = as.numeric(assigned_share$assigned_share_approx),
  stringsAsFactors = FALSE
)
as_tbl <- as_tbl[as_tbl$case_submitter_id %in% ret_cases &
  is.finite(as_tbl$dose_mgy) & as_tbl$dose_mgy > 0 &
  is.finite(as_tbl$assigned_share), , drop = FALSE]
band <- ifelse(as_tbl$assigned_share <= AS_LOW_MAX, "R_Low",
  ifelse(as_tbl$assigned_share <= AS_MID_MAX, "R_Mid", NA_character_))
as_tbl$band <- band
as_tbl <- as_tbl[!is.na(as_tbl$band), , drop = FALSE]

# Paired _merged Primary Tumor sample per case.
cd <- as.data.frame(colData(se))
m <- cd[grepl("_merged", cd$sample_submitter_id), , drop = FALSE]
tumor_of <- setNames(
  m$sample_submitter_id[m$sample_type == "Primary Tumor"],
  m$case_submitter_id[m$sample_type == "Primary Tumor"]
)
as_tbl$tumor_id <- unname(tumor_of[as_tbl$case_submitter_id])
as_tbl <- as_tbl[!is.na(as_tbl$tumor_id), , drop = FALSE]
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

classify <- function(score, b) {
  ifelse(score <= b$negative_max & (b$has_gap | score < b$positive_min), "negative",
    ifelse(score >= b$positive_min & (b$has_gap | score > b$negative_max), "positive",
      "undetermined"))
}
as_tbl$class <- classify(as_tbl$score, boundary)

# --- Report: graded validation ---------------------------------------------
n_pairs <- nrow(panel)
message(sprintf("\nPanel size %d ; boundary negative<=%d positive>=%d",
  n_pairs, boundary$negative_max, boundary$positive_min))
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
  config = list(as_low_max = AS_LOW_MAX, as_mid_max = AS_MID_MAX,
    panel_size = n_pairs, note = "R_Low/R_Mid unfiltered for purity/outliers"),
  samples = as_tbl[, c("case_submitter_id", "band", "assigned_share",
    "dose_mgy", "tumor_id", "score", "class")],
  summary = summary_tbl
)
out_rds <- file.path(paths$processed, "thyr_reo_evaluation.rds")
saveRDS(thyr_reo_evaluation, out_rds)
message("Saved: ", out_rds)

if (dir.exists(paths$output)) {
  out_csv <- file.path(paths$output, "reo_evaluation_samples.csv")
  utils::write.csv(thyr_reo_evaluation$samples, out_csv, row.names = FALSE)
  message("Saved: ", out_csv)
}
