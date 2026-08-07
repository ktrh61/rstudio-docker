# extract_baseline_20260807.R
# Phase-2 baseline extraction (dev QA; see reorg plan v2 §3.0-5).
# Reads the pre-migration processed/ artifacts and writes small machine-readable
# reference CSVs to tests/reference/baseline_20260807/.
#   - Case-ID files are the INVARIANTS the phase-3 refactor must preserve.
#   - Statistical files are the AUDIT RECORD for the intentional phase-4 changes
#     (Storey migration); they are NOT regression targets.
# Known mixed state, recorded as-is: thyr_normalized_counts.rds was regenerated
# without Cook's removal (2026-07-24), while thyr_expression_test.rds and
# thyr_enrichment_test.rds still rest on the earlier with-Cook's normalization.
# Run from the repository root: Rscript tests/reference/extract_baseline_20260807.R

source("setup.R")

out_dir <- file.path(paths$root, "tests", "reference", "baseline_20260807")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

outliers <- readRDS(file.path(paths$processed, "thyr_case_outliers.rds"))
purity <- readRDS(file.path(paths$processed, "thyr_case_purity.rds"))
normalized <- readRDS(file.path(paths$processed, "thyr_normalized_counts.rds"))
expression_test <- readRDS(file.path(paths$processed, "thyr_expression_test.rds"))
enrichment_test <- readRDS(file.path(paths$processed, "thyr_enrichment_test.rds"))
reo_panel <- readRDS(file.path(paths$processed, "thyr_reo_panel.rds"))
reo_eval <- readRDS(file.path(paths$processed, "thyr_reo_evaluation.rds"))

PURITY_THRESHOLD <- 0.6 # convention fixed in the current pipeline (310)

# --- Case-ID invariants ----------------------------------------------------
cases <- merge(
  outliers,
  purity[, c("case_submitter_id", "tumor_purity")],
  by = "case_submitter_id", all.x = TRUE, sort = TRUE
)
cases$purity_pass_0.6 <- as.integer(
  !is.na(cases$tumor_purity) & cases$tumor_purity >= PURITY_THRESHOLD
)
write.csv(cases, file.path(out_dir, "cases_main_bm.csv"), row.names = FALSE)

unit_samples <- do.call(rbind, lapply(names(normalized$units), function(unit) {
  s <- normalized$units[[unit]]$dgelist$samples
  data.frame(
    unit = unit,
    sample_id = rownames(s),
    group = as.character(s$group),
    lib.size = s$lib.size,
    norm.factors = s$norm.factors,
    scaling_coeff = s$scaling_coeff,
    stringsAsFactors = FALSE
  )
}))
write.csv(unit_samples, file.path(out_dir, "unit_samples.csv"), row.names = FALSE)

write.csv(reo_eval$samples, file.path(out_dir, "reo_eval_cases.csv"),
  row.names = FALSE
)
write.csv(reo_panel$panel, file.path(out_dir, "reo_panel.csv"), row.names = FALSE)

# --- Statistical audit baselines -------------------------------------------
gene_summary <- do.call(rbind, lapply(names(expression_test$units), function(unit) {
  g <- expression_test$units[[unit]]$genes
  data.frame(
    unit = unit,
    n_genes = nrow(g),
    min_p_exact = min(g$p_exact),
    n_p_lt_1e2 = sum(g$p_exact < 1e-2),
    n_p_lt_1e3 = sum(g$p_exact < 1e-3),
    n_p_lt_1e4 = sum(g$p_exact < 1e-4),
    n_fdr_perm_lt_0.10 = sum(g$fdr_perm < 0.10),
    stringsAsFactors = FALSE
  )
}))
write.csv(gene_summary, file.path(out_dir, "gene_level_summary.csv"),
  row.names = FALSE
)

omnibus <- do.call(rbind, lapply(names(expression_test$units), function(unit) {
  cbind(unit = unit, expression_test$units[[unit]]$omnibus)
}))
write.csv(omnibus, file.path(out_dir, "gene_level_omnibus.csv"), row.names = FALSE)

top50 <- do.call(rbind, lapply(names(expression_test$units), function(unit) {
  g <- expression_test$units[[unit]]$genes
  g <- g[order(-abs(g$statistic)), ][seq_len(50), ]
  cbind(unit = unit, g)
}))
write.csv(top50, file.path(out_dir, "gene_level_top50.csv"), row.names = FALSE)

gsea_cells <- do.call(rbind, lapply(names(enrichment_test$units), function(unit) {
  sets <- enrichment_test$units[[unit]]
  do.call(rbind, lapply(split(sets, sets$collection), function(s) {
    data.frame(
      unit = unit,
      collection = s$collection[1],
      n_sets = nrow(s),
      min_fwer = min(s$fwer, na.rm = TRUE),
      n_fwer_lt_0.05 = sum(s$fwer < 0.05, na.rm = TRUE),
      min_fdr_subramanian = min(s$fdr_subramanian, na.rm = TRUE),
      n_fdr_lt_0.05 = sum(s$fdr_subramanian < 0.05, na.rm = TRUE),
      n_fdr_lt_0.25 = sum(s$fdr_subramanian < 0.25, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(gsea_cells) <- NULL
write.csv(gsea_cells, file.path(out_dir, "gsea_12cell.csv"), row.names = FALSE)

reo_flat <- unlist(list(
  panel_config = reo_panel$config,
  boundary = reo_panel$boundary,
  eval_config = reo_eval$config,
  read_A_threshold = reo_eval$read_A_threshold,
  read_B = reo_eval$read_B
))
write.csv(
  data.frame(key = names(reo_flat), value = as.character(reo_flat)),
  file.path(out_dir, "reo_summary.csv"),
  row.names = FALSE
)

# --- Checksums -------------------------------------------------------------
snapshot_files <- sort(list.files(
  file.path(paths$processed, "_baseline_20260807"),
  full.names = TRUE
))
rds_md5 <- data.frame(
  file = basename(snapshot_files),
  md5 = unname(tools::md5sum(snapshot_files)),
  stringsAsFactors = FALSE
)
stat_md5 <- do.call(rbind, lapply(names(expression_test$units), function(unit) {
  g <- expression_test$units[[unit]]$genes
  tmp <- tempfile()
  writeLines(sprintf("%s %.15g", g$gene_id, g$statistic), tmp)
  data.frame(
    file = paste0("statistic_vector:", unit),
    md5 = unname(tools::md5sum(tmp)),
    stringsAsFactors = FALSE
  )
}))
write.csv(rbind(rds_md5, stat_md5), file.path(out_dir, "checksums.csv"),
  row.names = FALSE
)

message("Baseline written: ", out_dir)
