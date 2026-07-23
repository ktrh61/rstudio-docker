# 050_estimate_tumor_purity.R
# Estimate relative tumor purity per case with ContamDE (MUREN normalization),
# for the four analysis groups R_Sporadic / R_High / B_Sporadic / B_High.
# Input : processed/thyr_se_raw.rds              (from 021; single count assay)
#         processed/thyr_clinical.rds            (from 001; driver columns)
#         processed/thyr_case_assigned_share.rds (from 041; dose_mgy, AS)
#         utils/utils_improved.R, utils/norm_improved.R,
#         utils/contamde_purity_functions.R
# Output: processed/thyr_case_purity.rds
#           columns: case_submitter_id, group, tumor_purity
#
# ContamDE outputs a relative purity score normalized within the set of pairs
# passed together (the highest-purity pair is scaled to ~1), so estimation is
# run per group and the score is comparable only within a group.
#
# Group definition (this dataset):
#   driver  : BRAF  = Designated_Driver "BRAF.MutV600E" with no co-mutation
#                     (WGS/RNA CandidateDriverMutation is BRAF-only or empty/NA)
#             RET   = Designated_Driver in {CCDC6-RET, NCOA4-RET, RET-OTHER}
#   exposure: from 041 dose_mgy. dose_mgy == 0 -> Sporadic ;
#             dose_mgy > 0 with assigned_share_approx >= 66.6 -> High.
#             Exposed cases with AS < 66.6 (Low/Mid) are out of scope here.
#   group   = {R,B}_{Sporadic,High}
#
# Pairs are the _merged Primary Tumor / Solid Tissue Normal sample of each case
# (one of each; verified unique for this SE).

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
})

source(file.path(paths$root, "utils", "utils_improved.R"))
source(file.path(paths$root, "utils", "norm_improved.R"))
source(file.path(paths$root, "utils", "contamde_purity_functions.R"))

# ContamDE / MUREN group-internal parallel workers.
workers <- 3L

# AS threshold (percent) separating High from Low/Mid among exposed cases.
as_high_threshold <- 66.6

# BLAS/OMP single-threaded (MUREN parallelizes the outer loop).
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}

# --- Load inputs -----------------------------------------------------------
se_path <- file.path(paths$processed, "thyr_se_raw.rds")
if (!file.exists(se_path)) stop("thyr_se_raw.rds not found (run 021 first)")
se <- readRDS(se_path)
message(
  "SE: ", nrow(se), " genes x ", ncol(se), " samples ; assay: ",
  assayNames(se)
)

clinical_path <- file.path(paths$processed, "thyr_clinical.rds")
if (!file.exists(clinical_path)) stop("thyr_clinical.rds not found (run 001 first)")
clinical <- readRDS(clinical_path)
message("Clinical: ", nrow(clinical), " cases")

as_path <- file.path(paths$processed, "thyr_case_assigned_share.rds")
if (!file.exists(as_path)) stop("thyr_case_assigned_share.rds not found (run 041 first)")
assigned_share <- as.data.frame(readRDS(as_path))
message("Assigned Share: ", nrow(assigned_share), " cases")

# --- Driver classification (001) -------------------------------------------
# RET: exact match on the three RET fusion values of Designated_Driver.
# BRAF: BRAF.MutV600E with no co-mutation (CandidateDriverMutation BRAF-only or
# empty/NA in both WGS and RNA). The two masks are exclusive for this dataset.
ret_values <- c("CCDC6-RET", "NCOA4-RET", "RET-OTHER")
ret_mask <- clinical$Designated_Driver %in% ret_values

mut_ok <- function(x) is.na(x) | x == "" | x == "BRAF"
braf_mask <- !is.na(clinical$Designated_Driver) &
  clinical$Designated_Driver == "BRAF.MutV600E" &
  mut_ok(clinical$WGS_CandidateDriverMutation) &
  mut_ok(clinical$RNA_CandidateDriverMutation)

if (any(braf_mask & ret_mask)) {
  stop("A case matched both BRAF and RET driver masks (unexpected).")
}

driver <- rep(NA_character_, nrow(clinical))
driver[braf_mask] <- "BRAF"
driver[ret_mask] <- "RET"
message("Driver: BRAF ", sum(braf_mask), " ; RET ", sum(ret_mask))

driver_tbl <- data.frame(
  case_submitter_id = as.character(clinical$REBC_ID),
  driver = driver,
  stringsAsFactors = FALSE
)
driver_tbl <- driver_tbl[!is.na(driver_tbl$driver), , drop = FALSE]

# --- Exposure classification (041) -----------------------------------------
# dose_mgy == 0 -> Sporadic ; dose_mgy > 0 with AS >= threshold -> High.
# Other cases (dose_mgy NA, or exposed with AS < threshold) are not assigned an
# exposure band and drop out of scope.
as_tbl <- data.frame(
  case_submitter_id = as.character(assigned_share$REBC_ID),
  dose_mgy = as.numeric(assigned_share$dose_mgy),
  assigned_share = as.numeric(assigned_share$assigned_share_approx),
  stringsAsFactors = FALSE
)

exposure <- rep(NA_character_, nrow(as_tbl))
is_sporadic <- is.finite(as_tbl$dose_mgy) & as_tbl$dose_mgy == 0
is_high <- is.finite(as_tbl$dose_mgy) & as_tbl$dose_mgy > 0 &
  is.finite(as_tbl$assigned_share) & as_tbl$assigned_share >= as_high_threshold
exposure[is_sporadic] <- "Sporadic"
exposure[is_high] <- "High"
as_tbl$exposure <- exposure

# --- Pair resolution (021 colData) -----------------------------------------
# One _merged Primary Tumor and one _merged Solid Tissue Normal per case.
cd <- as.data.frame(colData(se))
is_merged <- grepl("_merged", cd$sample_submitter_id)
m <- cd[is_merged, , drop = FALSE]

t_rows <- m[m$sample_type == "Primary Tumor", , drop = FALSE]
n_rows <- m[m$sample_type == "Solid Tissue Normal", , drop = FALSE]

dup_t <- names(which(table(t_rows$case_submitter_id) > 1))
dup_n <- names(which(table(n_rows$case_submitter_id) > 1))
if (length(dup_t) || length(dup_n)) {
  stop(
    "A case has >1 _merged sample of the same type (pairing not unique): ",
    paste(c(dup_t, dup_n), collapse = ", ")
  )
}

tumor_of <- setNames(t_rows$sample_submitter_id, t_rows$case_submitter_id)
normal_of <- setNames(n_rows$sample_submitter_id, n_rows$case_submitter_id)
paired_cases <- intersect(names(tumor_of), names(normal_of))
message("Paired cases (_merged tumor & normal): ", length(paired_cases))

pair_tbl <- data.frame(
  case_submitter_id = paired_cases,
  tumor_id = unname(tumor_of[paired_cases]),
  normal_id = unname(normal_of[paired_cases]),
  stringsAsFactors = FALSE
)

# --- Assemble target cases -------------------------------------------------
# Driver-classified AND exposure-banded AND paired.
targets <- merge(driver_tbl, as_tbl[, c("case_submitter_id", "exposure")],
  by = "case_submitter_id"
)
targets <- targets[!is.na(targets$exposure), , drop = FALSE]
targets <- merge(targets, pair_tbl, by = "case_submitter_id")

prefix <- ifelse(targets$driver == "RET", "R", "B")
targets$group <- paste0(prefix, "_", targets$exposure)

message("Target cases by group:")
print(table(targets$group))

# --- Gene filter (protein_coding + filterByExpr) ---------------------------
gene_info <- as.data.frame(rowData(se))
is_protein_coding <- gene_info$gene_type == "protein_coding"
counts_all <- assay(se)

filter_genes <- function(normal_ids, tumor_ids) {
  cg <- counts_all[, c(normal_ids, tumor_ids), drop = FALSE]
  grp <- factor(c(
    rep("Normal", length(normal_ids)),
    rep("Tumor", length(tumor_ids))
  ))
  keep_expr <- edgeR::filterByExpr(cg, group = grp)
  is_protein_coding & keep_expr
}

# --- Per-group purity estimation -------------------------------------------
groups <- c("R_Sporadic", "R_High", "B_Sporadic", "B_High")
purity_by_case <- list()

for (g in groups) {
  gi <- targets[targets$group == g, , drop = FALSE]
  if (nrow(gi) == 0) {
    message("Group ", g, ": no target pairs; skipped")
    next
  }
  message("Group ", g, ": ", nrow(gi), " pairs")

  keep <- filter_genes(gi$normal_id, gi$tumor_id)
  counts <- cbind(
    counts_all[keep, gi$normal_id, drop = FALSE],
    counts_all[keep, gi$tumor_id, drop = FALSE]
  )
  np <- nrow(gi)
  colnames(counts) <- c(
    paste0("Normal_", seq_len(np)),
    paste0("Tumor_", seq_len(np))
  )
  message(
    "  genes after filter: ", sum(keep), " ; matrix ",
    nrow(counts), " x ", ncol(counts)
  )

  res <- tryCatch(
    contamde_purity(
      counts = counts,
      subtype = NULL,
      covariate = NULL,
      contaminated = TRUE,
      pairwise_method = "lts",
      workers = workers,
      verbose = FALSE
    ),
    error = function(e) {
      warning("Group ", g, " purity estimation failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(res)) {
    purity_by_case[[g]] <- data.frame(
      case_submitter_id = gi$case_submitter_id,
      group = g,
      tumor_purity = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }

  purity_by_case[[g]] <- data.frame(
    case_submitter_id = gi$case_submitter_id,
    group = g,
    tumor_purity = as.numeric(res$proportion),
    stringsAsFactors = FALSE
  )
  message(sprintf(
    "  purity: mean=%.3f sd=%.3f range=[%.3f, %.3f]",
    mean(res$proportion), sd(res$proportion),
    min(res$proportion), max(res$proportion)
  ))
}

# --- Assemble and save -----------------------------------------------------
thyr_case_purity <- do.call(rbind, purity_by_case)
rownames(thyr_case_purity) <- NULL
thyr_case_purity <- thyr_case_purity[
  order(thyr_case_purity$group, thyr_case_purity$case_submitter_id), ,
  drop = FALSE
]

out <- file.path(paths$processed, "thyr_case_purity.rds")
saveRDS(thyr_case_purity, out)
message(
  "Saved: ", out, " (", nrow(thyr_case_purity), " cases x ",
  ncol(thyr_case_purity), " columns)"
)
