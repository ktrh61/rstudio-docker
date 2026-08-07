# 420_test_gene_sets.R
# Gene-set enrichment for each analysis unit, on the ranking produced by 410.
# Genes are ranked by the signed Brunner-Munzel statistic and the null comes
# from the same label shuffles 410 used, so the gene-level and set-level
# results rest on one permutation null.
# Input : processed/thyr_normalized_counts.rds  (from 310; per-unit DGEList)
#         processed/thyr_expression_test.rds    (from 410; n_perm and seed)
#         lib/stat_brunnermunzel.R, lib/gsea_permutation.R
# Output: processed/thyr_enrichment_test.rds
#           list(date, config, collections, units)
#           units = { R_Tumor, R_Normal, B_Tumor, B_Normal }; each a data frame
#             sets : collection, pathway, size, ES, NES, pval, fwer,
#                    fdr_subramanian, redundant_with, leading_edge
#
# No differentially expressed gene list is used or produced: GSEA consumes the
# whole ranking, which is what makes it the right test when the signal is
# diffuse and no single gene survives correction.
#
# Only process-definition collections are used. Experiment-derived collections
# (C2:CGP, C6, C7) are excluded because their sets are built from perturbation
# experiments -- acute high-dose irradiation, oncogene activation -- whose
# confounding this design exists to avoid. NES and fwer are computed within each
# collection, so adding an exploratory collection cannot dilute the primary one.
#
# Inference is the Westfall-Young family-wise error rate (fwer), and nothing
# else: its null is the largest |NES| anywhere in the collection per shuffle, so
# inter-pathway correlation is carried exactly. The choice matches the level's
# sparsity -- a collection carries signal in one or two sets, so the maximum is
# the right statistic, just as the diffuse gene level in 410 called for a count.
# A spike-in positive control (a 15% shift over one Hallmark set, planted in a
# unit with no signal) is caught at fwer 0.003, so an empty result is a property
# of the data. The claim is only that no set survives fwer correction; no test
# of a diffuse pathway excess is attempted, and none is claimed. pval and
# fdr_subramanian are descriptive: pval is the uncorrected per-set value used
# for ordering, and fdr_subramanian (Subramanian et al. 2005) is carried to
# document why standard GSEA fails here -- it divides by the mean of the pooled
# null, which under inter-pathway correlation follows whichever way the
# collection drifts in one realization, and on this data calls 74% of Hallmark
# in a unit with no gene-level signal.
#
# The size window drops both noise-prone tiny sets and the vague giant ones;
# redundant_with flags a set whose leading edge is largely contained in a
# better-ranked set, which is how an ontology's parent and child terms both
# surface on a single signal.

source("setup.R")

suppressPackageStartupMessages({
  library(edgeR)
  library(msigdbr)
})

source(file.path(paths$root, "lib", "stat_brunnermunzel.R"))
source(file.path(paths$root, "lib", "gsea_permutation.R"))
source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "annotation.R"))

# --- Configuration ---------------------------------------------------------
SPECIES <- "Homo sapiens"
COLLECTIONS <- list(
  "H" = list(collection = "H", subcollection = NA_character_),
  "C2:CP" = list(
    collection = "C2",
    subcollection = c(
      "CP:REACTOME", "CP:WIKIPATHWAYS", "CP:KEGG_MEDICUS",
      "CP:BIOCARTA", "CP:PID"
    )
  ),
  "C5:GO:BP" = list(collection = "C5", subcollection = "GO:BP")
)
PRIMARY_COLLECTION <- "H" # the pre-specified inferential collection
# Unweighted Kolmogorov-Smirnov enrichment score. The studentized Brunner-
# Munzel statistic can reach +/-Inf (its denominator vanishes for some
# allocations), so a magnitude weighting would let a single such gene dominate
# a pathway's score and, in permutations, inflate whichever pathway's max-|NES|
# null happens to contain it -- corrupting the Westfall-Young null. Rank-only
# scoring is bounded and avoids this. The choice follows from the metric's tail,
# known a priori, not from the results: at gseaParam = 1 the primary conclusion
# is the same (R_Tumor Hallmark min fwer 0.64 vs 0.83 here, both null).
GSEA_PARAM <- 0
MIN_SET_SIZE <- 15L
MAX_SET_SIZE <- 500L
REDUNDANT_JACCARD <- 0.5
REDUNDANT_MAX_PVAL <- 0.05 # sets below this get a redundancy annotation
# EXACT_THREADS / BM_EXACT_MAX come from config.R via setup.R. 420 computes BM
# statistics only (no p-value path), so exact.max.allocations is inert here;
# it is set for uniformity with 310/410.

options(
  brunnermunzel.exact.max.allocations = BM_EXACT_MAX,
  brunnermunzel.exact.threads = EXACT_THREADS
)

pin_blas_threads()

# --- Load inputs -----------------------------------------------------------
norm_path <- file.path(paths$processed, "thyr_normalized_counts.rds")
test_path <- file.path(paths$processed, "thyr_expression_test.rds")
if (!file.exists(norm_path)) stop("thyr_normalized_counts.rds not found (310)")
if (!file.exists(test_path)) stop("thyr_expression_test.rds not found (410)")
normalized <- readRDS(norm_path)
expression_test <- readRDS(test_path)

N_PERM <- expression_test$config$n_perm
PERM_SEED <- expression_test$config$perm_seed
message("Reusing 410's permutation null: ", N_PERM, " shuffles, seed ", PERM_SEED)

# --- Gene sets -------------------------------------------------------------
collect_sets <- function(spec) {
  frames <- lapply(spec$subcollection, function(sub) {
    if (is.na(sub)) {
      msigdbr(species = SPECIES, collection = spec$collection)
    } else {
      msigdbr(
        species = SPECIES, collection = spec$collection, subcollection = sub
      )
    }
  })
  frame <- do.call(rbind, frames)
  split(frame$ensembl_gene, frame$gs_name)
}

gene_sets <- lapply(COLLECTIONS, collect_sets)
message("Gene sets before filtering: ", paste(
  names(gene_sets), lengths(gene_sets), sep = "=", collapse = " "
))

# --- Per-unit enrichment ---------------------------------------------------
test_unit <- function(dgelist, unit) {
  arms <- unit_arms(dgelist$samples$group, unit)
  sporadic <- arms$sporadic
  high <- arms$high
  cpm_matrix <- edgeR::cpm(
    dgelist,
    normalized.lib.sizes = TRUE, prior.count = 0, log = FALSE
  )[, c(sporadic, high), drop = FALSE]
  storage.mode(cpm_matrix) <- "double"
  # MSigDB carries unversioned Ensembl identifiers.
  ensembl <- strip_ensembl_version(rownames(cpm_matrix))
  cpm_matrix <- cpm_matrix[!duplicated(ensembl), , drop = FALSE]
  rownames(cpm_matrix) <- ensembl[!duplicated(ensembl)]
  nx <- length(sporadic)
  n <- ncol(cpm_matrix)

  index <- lapply(gene_sets, gsea_pathway_index,
    gene_ids = rownames(cpm_matrix),
    min_size = MIN_SET_SIZE, max_size = MAX_SET_SIZE
  )
  collection_of <- rep(names(index), lengths(index))
  flat <- unlist(unname(index), recursive = FALSE)

  metric <- brunnermunzel_statistics(cpm_matrix, nx)
  observed <- gsea_scores(metric, flat, GSEA_PARAM)

  set.seed(PERM_SEED)
  null <- vapply(
    seq_len(N_PERM),
    function(i) {
      gsea_scores(
        brunnermunzel_statistics(
          cpm_matrix[, sample(n), drop = FALSE], nx
        ),
        flat, GSEA_PARAM
      )
    },
    numeric(length(flat))
  )

  leading_edge <- gsea_leading_edge(metric, flat, rownames(cpm_matrix), GSEA_PARAM)
  per_collection <- lapply(names(index), function(collection) {
    keep <- collection_of == collection
    normalized <- gsea_nes(observed[keep], null[keep, , drop = FALSE])
    result <- data.frame(
      collection = collection,
      pathway = names(observed)[keep],
      size = lengths(flat[keep]),
      ES = unname(observed[keep]),
      NES = unname(normalized$nes),
      pval = gsea_pathway_pvalues(normalized$nes, normalized$nes_null),
      fwer = gsea_westfall_young(normalized$nes, normalized$nes_null),
      fdr_subramanian = gsea_subramanian_fdr(
        normalized$nes, normalized$nes_null
      ),
      stringsAsFactors = FALSE
    )
    result$redundant_with <- gsea_redundancy(
      result, leading_edge, REDUNDANT_JACCARD,
      candidate = !is.na(result$pval) & result$pval < REDUNDANT_MAX_PVAL
    )
    result$leading_edge <- I(unname(leading_edge[result$pathway]))
    result[order(result$fwer, result$pval, -abs(result$NES)), , drop = FALSE]
  })

  sets <- do.call(rbind, per_collection)
  rownames(sets) <- NULL

  primary <- sets[sets$collection == PRIMARY_COLLECTION, , drop = FALSE]
  message(sprintf(
    "  %-9s %d sets | %s min fwer %.3f (fwer<0.05 %d) | top: %s | [ref] fdr_subramanian<0.25 %d",
    unit, nrow(sets), PRIMARY_COLLECTION,
    min(primary$fwer, na.rm = TRUE), sum(primary$fwer < 0.05, na.rm = TRUE),
    sub("^HALLMARK_", "", primary$pathway[1]),
    sum(sets$fdr_subramanian < 0.25, na.rm = TRUE)
  ))
  sets
}

message("Testing units:")
units <- lapply(names(normalized$units), function(unit) {
  test_unit(normalized$units[[unit]]$dgelist, unit)
})
names(units) <- names(normalized$units)

# --- Assemble and save -----------------------------------------------------
thyr_enrichment_test <- list(
  date = Sys.Date(),
  config = list(
    ranking = "signed Brunner-Munzel statistic",
    permutation = "sample labels, shared with 410",
    n_perm = N_PERM,
    perm_seed = PERM_SEED,
    primary_collection = PRIMARY_COLLECTION,
    excluded_collections = c("C2:CGP", "C6", "C7"),
    min_set_size = MIN_SET_SIZE,
    max_set_size = MAX_SET_SIZE,
    redundant_jaccard = REDUNDANT_JACCARD,
    redundant_max_pval = REDUNDANT_MAX_PVAL,
    gsea_param = GSEA_PARAM,
    inference = paste(
      "Westfall-Young fwer, within collection (primary and only);",
      "pval and fdr_subramanian are descriptive"
    ),
    msigdbr_version = as.character(packageVersion("msigdbr")),
    fgsea_version = as.character(packageVersion("fgsea"))
  ),
  collections = names(COLLECTIONS),
  units = units
)

out <- file.path(paths$processed, "thyr_enrichment_test.rds")
saveRDS(thyr_enrichment_test, out)
message("Saved: ", out, " (", length(units), " units)")
