# gsea_spikein_control.R
# Spike-in positive control for the spec-B gene-set inference (reorg plan v2,
# phase 4 step 15). The historical control (a 15% shift over one Hallmark set
# planted in a unit with no gene-level signal, caught at Westfall-Young fwer
# 0.003 under the old spec) was an ad-hoc measurement; this script makes it a
# recorded diagnostic and re-runs it under the adopted inference (tie-averaged
# normal scores, tie-block ES, per-set permutation p + BH within collection;
# reorg plan v2 appendix B).
#
# The point of the control is that an empty result on real data is a property
# of the data, not of a dead pipeline: a planted signal of modest, biologically
# plausible size must be recovered. Run after the D6 null calibration; reads
# only 310's normalized counts plus a synthetic perturbation.
# Input : processed/thyr_normalized_counts.rds  (from 310)
#         lib/stat_brunnermunzel.R, lib/gsea_permutation.R,
#         lib/gsea_collections.R
# Output: diagnostics/output/gsea_spikein_control.rds

source("setup.R")

suppressPackageStartupMessages({
  library(edgeR)
  library(msigdbr)
  library(parallel)
})

source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))
source(file.path(paths$root, "lib", "gsea_permutation.R"))
source(file.path(paths$root, "lib", "gsea_collections.R"))
source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "annotation.R"))

# --- Configuration ---------------------------------------------------------
SPECIES <- "Homo sapiens"
SPIKE_UNIT <- "B_Tumor" # no gene-level signal against its own null
SPIKE_SET <- "HALLMARK_ADIPOGENESIS" # one mid-size Hallmark set
SPIKE_FACTOR <- 1.15 # 15% multiplicative shift on the High arm
SPIKE_SEED <- 20260807L # shuffle seed; independent of the inference seed.
# Coincides with the calibration's CALIB_SEED (both are the D6 decision
# date), so for a given unit the two diagnostics share a shuffle stream.
# Harmless: they are separate measurements, never combined.

pin_blas_threads()

# --- Load inputs -----------------------------------------------------------
norm_path <- file.path(paths$processed, "thyr_normalized_counts.rds")
if (!file.exists(norm_path)) stop("thyr_normalized_counts.rds not found (310)")
normalized <- readRDS(norm_path)

collections <- load_gene_set_collections(SPECIES)
gene_sets <- collections$gene_sets

out_dir <- file.path(paths$root, "diagnostics", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Build the spiked matrix ------------------------------------------------
dgelist <- normalized$units[[SPIKE_UNIT]]$dgelist
arms <- unit_arms(dgelist$samples$group, SPIKE_UNIT)
cpm_matrix <- edgeR::cpm(
  dgelist,
  normalized.lib.sizes = TRUE, prior.count = 0, log = FALSE
)[, c(arms$sporadic, arms$high), drop = FALSE]
storage.mode(cpm_matrix) <- "double"
ensembl <- strip_ensembl_version(rownames(cpm_matrix))
cpm_matrix <- cpm_matrix[!duplicated(ensembl), , drop = FALSE]
rownames(cpm_matrix) <- ensembl[!duplicated(ensembl)]
nx <- length(arms$sporadic)
n <- ncol(cpm_matrix)

spike_members <- intersect(gene_sets[["H"]][[SPIKE_SET]], rownames(cpm_matrix))
if (length(spike_members) == 0L) stop("Spike set has no genes in the data.")
high_columns <- (nx + 1L):n
cpm_matrix[spike_members, high_columns] <-
  cpm_matrix[spike_members, high_columns] * SPIKE_FACTOR
message(sprintf(
  "Spiked %s (%d genes present) by %.2fx in %s High arm (%d samples)",
  SPIKE_SET, length(spike_members), SPIKE_FACTOR, SPIKE_UNIT,
  length(high_columns)
))

# --- Run the spec-B inference on the spiked unit ----------------------------
index <- lapply(gene_sets, gsea_pathway_index,
  gene_ids = rownames(cpm_matrix),
  min_size = GSEA_MIN_SET_SIZE, max_size = GSEA_MAX_SET_SIZE
)
collection_of <- rep(names(index), lengths(index))
flat <- unlist(unname(index), recursive = FALSE)

scores <- gsea_normal_scores(brunnermunzel_statistics(cpm_matrix, nx))
observed <- gsea_block_scores(scores, flat)

set.seed(SPIKE_SEED)
perm_index <- vapply(seq_len(N_PERM), function(i) sample(n), integer(n))
null_columns <- mclapply(seq_len(N_PERM), function(i) {
  gsea_block_scores(
    gsea_normal_scores(brunnermunzel_statistics(
      cpm_matrix[, perm_index[, i], drop = FALSE], nx
    )),
    flat
  )
}, mc.cores = WORKERS)
null <- gsea_bind_null_columns(null_columns, length(flat))
rownames(null) <- names(flat)

sets <- do.call(rbind, lapply(names(index), function(collection) {
  keep <- collection_of == collection
  standardized <- gsea_nes(observed[keep], null[keep, , drop = FALSE])
  pval <- gsea_pathway_pvalues(standardized$nes, standardized$nes_null)
  data.frame(
    collection = collection,
    pathway = names(observed)[keep],
    size = lengths(flat[keep]),
    ES = unname(observed[keep]),
    NES = unname(standardized$nes),
    pval = pval,
    q_bh = p.adjust(pval, method = "BH"),
    stringsAsFactors = FALSE
  )
}))
rownames(sets) <- NULL

spiked_row <- sets[sets$pathway == SPIKE_SET, , drop = FALSE]
# A control that never entered the tested universe must not read as a pass.
if (nrow(spiked_row) != 1L) {
  stop("Spiked set ", SPIKE_SET, " is not in the tested universe.")
}
hallmark <- sets[sets$collection == "H", , drop = FALSE]
hallmark_rank <- match(SPIKE_SET, hallmark$pathway[order(hallmark$q_bh)])
message(sprintf(
  "Spiked set: NES %.2f ; pval %.4f ; q_bh %.4f (rank %d/%d in H) ; recovered at q<%.2f: %s",
  spiked_row$NES, spiked_row$pval, spiked_row$q_bh, hallmark_rank,
  nrow(hallmark), FDR_CUT, spiked_row$q_bh < FDR_CUT
))
message(sprintf(
  "Hallmark sets besides the spike at q<%.2f: %d",
  FDR_CUT, sum(hallmark$q_bh < FDR_CUT, na.rm = TRUE) -
    as.integer(spiked_row$q_bh < FDR_CUT)
))

# --- Save ------------------------------------------------------------------
gsea_spikein_control <- list(
  date = Sys.Date(),
  config = list(
    unit = SPIKE_UNIT, set = SPIKE_SET, factor = SPIKE_FACTOR,
    n_members_present = length(spike_members),
    seed = SPIKE_SEED, n_perm = N_PERM, q_threshold = FDR_CUT,
    msigdbr_version = as.character(packageVersion("msigdbr"))
  ),
  spiked_row = spiked_row,
  hallmark_rank = hallmark_rank,
  sets = sets
)
out <- file.path(out_dir, "gsea_spikein_control.rds")
saveRDS(gsea_spikein_control, out)
message("Saved: ", out)
