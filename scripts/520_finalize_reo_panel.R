# 520_finalize_reo_panel.R
# Build the final REO panel from the candidate pairs (510) by greedy, non-
# redundant selection, then set a dose-zero-anchored classification boundary.
# Input : processed/thyr_reo_candidate_pairs.rds  (from 510)
#         processed/thyr_normalized_counts.rds     (from 310; R_Tumor sample set)
#         processed/thyr_se_raw.rds                (from 120; single count assay)
#         processed/gene_lengths.rds               (from 020)
#         lib/reo.R
# Output: processed/thyr_reo_panel.rds, output/reo_panel.csv
#
# Greedy selection walks candidate pairs in descending median_diff and keeps a
# pair only if it reuses no gene already in the panel and its per-sample r
# vector correlates (Spearman) below CORRELATION_THRESHOLD with every kept pair,
# stopping at TARGET_PANEL_SIZE. A sample's panel score is the number of pairs
# that reverse (|r| >= dead zone and opposite sign to the R0 reference). The
# classification cutoff is the maximum R0 score, and scores above that cutoff
# are positive. R1 scores describe construction fit but do not determine the
# cutoff. Script 530 applies the resulting rule to R_Low/R_Mid unchanged.

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "reo.R"))

# --- Configuration ---------------------------------------------------------
TARGET_PANEL_SIZE <- 10L
CORRELATION_THRESHOLD <- 0.75 # Spearman; drop a pair correlated above this
# DEAD_ZONE comes from config.R via setup.R (shared with 510).

# --- Load inputs -----------------------------------------------------------
cand_path <- file.path(paths$processed, "thyr_reo_candidate_pairs.rds")
norm_path <- file.path(paths$processed, "thyr_normalized_counts.rds")
se_path <- file.path(paths$processed, "thyr_se_raw.rds")
len_path <- file.path(paths$processed, "gene_lengths.rds")
for (p in c(cand_path, norm_path, se_path, len_path)) {
  if (!file.exists(p)) stop("missing input: ", p)
}
candidates <- readRDS(cand_path)$pairs
dge <- readRDS(norm_path)$units$R_Tumor$dgelist
se <- readRDS(se_path)
gene_lengths <- readRDS(len_path)

arms <- unit_arms(dge$samples$group, "R_Tumor")
r0_samples <- colnames(dge)[arms$sporadic]
r1_samples <- colnames(dge)[arms$high]
all_samples <- c(r0_samples, r1_samples)
log2_tpm <- reo_log2_tpm(se, gene_lengths, all_samples)
message("Candidates: ", nrow(candidates), " ; R0 ", length(r0_samples), " R1 ", length(r1_samples))

# --- Greedy non-redundant selection ----------------------------------------
candidates <- candidates[order(candidates$median_diff, decreasing = TRUE), , drop = FALSE]
used_genes <- character(0)
selected_idx <- integer(0)
selected_r <- NULL # samples x selected-pairs, for correlation checks

for (k in seq_len(nrow(candidates))) {
  up_g <- candidates$up[k]
  down_g <- candidates$down[k]
  if (up_g %in% used_genes || down_g %in% used_genes) next

  new_r <- log2_tpm[up_g, all_samples] - log2_tpm[down_g, all_samples]
  if (!is.null(selected_r)) {
    cors <- apply(selected_r, 2, function(x) {
      abs(stats::cor(x, new_r, method = "spearman"))
    })
    cors[is.na(cors)] <- 0
    if (any(cors >= CORRELATION_THRESHOLD)) next
  }

  selected_idx <- c(selected_idx, k)
  used_genes <- c(used_genes, up_g, down_g)
  selected_r <- cbind(selected_r, new_r)
  if (length(selected_idx) >= TARGET_PANEL_SIZE) break
}

panel <- candidates[selected_idx, , drop = FALSE]
panel$pair_id <- paste0("P", seq_len(nrow(panel)))
rownames(panel) <- NULL
message("Selected panel pairs: ", nrow(panel),
  if (nrow(panel) < TARGET_PANEL_SIZE) paste0(" (below target ", TARGET_PANEL_SIZE, ")") else "")

# --- Per-sample reversal score on the training arms (lib/reo.R) ------------
r0_score <- reversal_score(log2_tpm, r0_samples, panel, DEAD_ZONE)
r1_score <- reversal_score(log2_tpm, r1_samples, panel, DEAD_ZONE)
message(sprintf("R0 score range [%d, %d] ; R1 score range [%d, %d]",
  min(r0_score), max(r0_score), min(r1_score), max(r1_score)))

# --- Boundary: threshold from the Sporadic (R0) arm only -------------------
# The training High arm (R1) is a mixture and contains Sporadic-like cases with
# zero reversals, so a min(R1)-based positive threshold collapses. Instead a
# sample scores positive when it shows MORE pair reversals than any Sporadic
# case: threshold = max(R0 score). This is mixture-robust and interpretable, and
# is the classification-based read (A) of the panel; the graded permutation read
# (B) is done out-of-sample in 530.
max_r0 <- max(r0_score)
boundary <- list(negative_max = max_r0, positive_min = max_r0 + 1L)
message(sprintf("Boundary (R0-based): positive if score > %d", max_r0))
message("Training classification (R0 Sporadic / R1 High):")
print(table(arm = c(rep("R0", length(r0_score)), rep("R1", length(r1_score))),
  class = c(classify_reversal(r0_score, boundary), classify_reversal(r1_score, boundary))))

# --- Assemble and save -----------------------------------------------------
thyr_reo_panel <- list(
  date = Sys.Date(),
  config = list(
    target_panel_size = TARGET_PANEL_SIZE,
    correlation_threshold = CORRELATION_THRESHOLD,
    dead_zone = DEAD_ZONE
  ),
  panel = panel,
  boundary = boundary,
  training = list(
    r0_score = r0_score, r1_score = r1_score,
    r0_samples = r0_samples, r1_samples = r1_samples
  )
)
out_rds <- file.path(paths$processed, "thyr_reo_panel.rds")
saveRDS(thyr_reo_panel, out_rds)
message("Saved: ", out_rds, " (", nrow(panel), " pairs)")

if (dir.exists(paths$output)) {
  out_csv <- file.path(paths$output, "reo_panel.csv")
  utils::write.csv(
    panel[, c("pair_id", "up", "up_name", "down", "down_name",
      "median_diff", "reversal_rate", "r0_q10")],
    out_csv, row.names = FALSE
  )
  message("Saved: ", out_csv)
}
