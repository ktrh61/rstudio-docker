# normal_signature_agreement.R
# Pre-specified descriptive comparison (reorg plan v2 s0.5, 6th confirmation,
# ratified 2026-08-07): are the exposure traces in RET-case normals and
# BRAF-case normals the same trace? Measured as the Spearman agreement of the
# two units' signed Brunner-Munzel statistic vectors over shared genes --
# threshold-free, no dose assumption, same form as the purity verification.
#
# Readings were fixed before any run (s0.5): high agreement -> consistent
# with a shared glandular trace; no agreement -> identity is UNDECIDABLE
# ("different biology" and "one side has no signal" cannot be distinguished).
# Descriptive only; feeds the s0.6 pattern rules, changes nothing else.
#
# The point estimate travels with a permutation-calibrated reference (same
# convention as the pi0 reporting, v2 D1): labels are shuffled independently
# within each unit (the saved perm_index of each, by reference) and the
# shuffle-b R vector is correlated with the shuffle-b B vector, giving the
# null spread of rho when NEITHER cohort carries label-aligned structure.
# Inter-gene correlation makes this spread far wider than 1/sqrt(n_genes),
# which is why the bare coefficient is uninterpretable without it. An
# observed rho outside this spread says the two contrasts share
# label-aligned structure; it does NOT by itself say the structure is the
# exposure trace rather than a shared covariate (age structure, purity) --
# that separation belongs to the covariate diagnostics (s0.5 2nd).
# Input : processed/thyr_expression_test.rds    (from 410; effect, statistic,
#         perm_index)
#         processed/thyr_normalized_counts.rds  (from 310; per-unit DGEList)
# Output: diagnostics/output/normal_signature_agreement.rds

source("setup.R")

suppressPackageStartupMessages({
  library(edgeR)
  library(parallel)
})

source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))
source(file.path(paths$root, "lib", "units.R"))

pin_blas_threads()

expression_test <- readRDS(
  file.path(paths$processed, "thyr_expression_test.rds")
)
normalized <- readRDS(
  file.path(paths$processed, "thyr_normalized_counts.rds")
)

signed_vector <- function(unit) {
  g <- expression_test$units[[unit]]$genes
  setNames(g$statistic * sign(g$effect - 0.5), g$gene_id)
}

unit_cpm <- function(unit) {
  dgelist <- normalized$units[[unit]]$dgelist
  arms <- unit_arms(dgelist$samples$group, unit)
  m <- edgeR::cpm(
    dgelist,
    normalized.lib.sizes = TRUE, prior.count = 0, log = FALSE
  )[, c(arms$sporadic, arms$high), drop = FALSE]
  storage.mode(m) <- "double"
  list(cpm = m, nx = length(arms$sporadic))
}

r_normal <- signed_vector("R_Normal")
b_normal <- signed_vector("B_Normal")
shared <- intersect(names(r_normal), names(b_normal))
rho <- stats::cor(
  r_normal[shared], b_normal[shared],
  method = "spearman"
)

# Null spread: shuffle-b vs shuffle-b signed statistics on the shared genes.
r_data <- unit_cpm("R_Normal")
b_data <- unit_cpm("B_Normal")
r_perm <- expression_test$units$R_Normal$perm_index
b_perm <- expression_test$units$B_Normal$perm_index
n_perm <- min(ncol(r_perm), ncol(b_perm))
rho_null <- unlist(mclapply(seq_len(n_perm), function(i) {
  rv <- brunnermunzel_statistics(
    r_data$cpm[shared, r_perm[, i], drop = FALSE], r_data$nx
  )
  bv <- brunnermunzel_statistics(
    b_data$cpm[shared, b_perm[, i], drop = FALSE], b_data$nx
  )
  stats::cor(rv, bv, method = "spearman")
}, mc.cores = WORKERS))

p_two <- (sum(abs(rho_null) >= abs(rho)) + 1) / (n_perm + 1)
message(sprintf(
  "Normal-side signature agreement (signed BM, Spearman): rho = %+.4f over %d shared genes",
  rho, length(shared)
))
message(sprintf(
  "  shuffle null: median %+.4f, 95%% band [%+.4f, %+.4f], max |rho| %.4f ; two-sided p %.4f",
  stats::median(rho_null),
  stats::quantile(rho_null, 0.025), stats::quantile(rho_null, 0.975),
  max(abs(rho_null)), p_two
))
message(sprintf(
  "Gene-level context: R_Normal %d DEG (hc p %.3f), B_Normal %d DEG (hc p %.3f)",
  sum(expression_test$units$R_Normal$genes$q_storey < FDR_CUT),
  expression_test$units$R_Normal$omnibus$p[
    expression_test$units$R_Normal$omnibus$test == "hc"
  ],
  sum(expression_test$units$B_Normal$genes$q_storey < FDR_CUT),
  expression_test$units$B_Normal$omnibus$p[
    expression_test$units$B_Normal$omnibus$test == "hc"
  ]
))

out_dir <- file.path(paths$root, "diagnostics", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out <- file.path(out_dir, "normal_signature_agreement.rds")
saveRDS(
  list(
    date = Sys.Date(),
    config = list(
      units = c("R_Normal", "B_Normal"),
      metric = "signed BM statistic (statistic * sign(effect - 0.5))",
      method = "spearman",
      null = "paired within-unit label shuffles (each unit's perm_index)",
      n_perm = n_perm,
      reading = "fixed in reorg plan v2 s0.5 6th before running"
    ),
    n_shared = length(shared),
    rho = rho,
    rho_null = rho_null,
    p_two_sided = p_two
  ),
  out
)
message("Saved: ", out)
