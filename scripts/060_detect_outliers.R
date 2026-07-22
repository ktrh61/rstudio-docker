# 060_detect_outliers.R
# Detect sample outliers with PC-OD (principal-component outlier detection),
# separately for Tumor and Normal within each of the four analysis groups, over
# the high-purity cases retained by 050.
# Input : processed/thyr_case_purity.rds   (from 050; case, group, tumor_purity)
#         processed/thyr_se_raw.rds         (from 021; single count assay)
#         utils/PC-OD_improved.R            (PC_OD)
# Output: processed/thyr_case_outliers.rds
#           columns: case_submitter_id, group, tumor_purity,
#                    has_outlier_tumor, has_outlier_normal   (0/1)
#
# Scope: cases with tumor_purity >= 0.6. Each such case contributes its _merged
# Primary Tumor and Solid Tissue Normal sample. PC-OD is run on eight groups:
# {R_Sporadic, R_High, B_Sporadic, B_High} x {Tumor, Normal}. Within each group
# the count sub-matrix is filtered by edgeR::filterByExpr and turned into log-CPM
# without normalization (normalized.lib.sizes = FALSE, prior.count = 2), then
# PC_OD returns the column indices flagged as outliers.

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
})

source(file.path(paths$root, "utils", "PC-OD_improved.R"))

purity_threshold <- 0.6

# --- Load inputs -----------------------------------------------------------
purity_path <- file.path(paths$processed, "thyr_case_purity.rds")
if (!file.exists(purity_path)) stop("thyr_case_purity.rds not found (run 050 first)")
purity <- readRDS(purity_path)
message("Purity table: ", nrow(purity), " cases")

se_path <- file.path(paths$processed, "thyr_se_raw.rds")
if (!file.exists(se_path)) stop("thyr_se_raw.rds not found (run 021 first)")
se <- readRDS(se_path)
message(
  "SE: ", nrow(se), " genes x ", ncol(se), " samples ; assay: ",
  assayNames(se)
)

# --- Select high-purity cases ----------------------------------------------
target <- purity[!is.na(purity$tumor_purity) &
  purity$tumor_purity >= purity_threshold, , drop = FALSE]
message("High-purity cases (>= ", purity_threshold, "): ", nrow(target))

# --- Resolve _merged Tumor / Normal sample per case ------------------------
cd <- as.data.frame(colData(se))
is_merged <- grepl("_merged", cd$sample_submitter_id)
m <- cd[is_merged, , drop = FALSE]
t_rows <- m[m$sample_type == "Primary Tumor", , drop = FALSE]
n_rows <- m[m$sample_type == "Solid Tissue Normal", , drop = FALSE]
tumor_of <- setNames(t_rows$sample_submitter_id, t_rows$case_submitter_id)
normal_of <- setNames(n_rows$sample_submitter_id, n_rows$case_submitter_id)

target$tumor_id <- unname(tumor_of[target$case_submitter_id])
target$normal_id <- unname(normal_of[target$case_submitter_id])
if (any(is.na(target$tumor_id)) || any(is.na(target$normal_id))) {
  stop("A high-purity case is missing its _merged Tumor or Normal sample.")
}

# --- PC-OD per (group x tissue) --------------------------------------------
counts_all <- assay(se)

# log-CPM without normalization for a set of sample columns, then PC-OD.
# Returns the outlier column positions (indices into sample_ids).
run_pcod <- function(sample_ids) {
  cts <- counts_all[, sample_ids, drop = FALSE]
  y <- edgeR::DGEList(counts = cts)
  keep <- edgeR::filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = TRUE]
  logCPM <- edgeR::cpm(
    y,
    log = TRUE,
    normalized.lib.sizes = FALSE,
    prior.count = 2
  )
  PC_OD(logCPM)
}

target$has_outlier_tumor <- 0L
target$has_outlier_normal <- 0L

groups <- c("R_Sporadic", "R_High", "B_Sporadic", "B_High")
for (g in groups) {
  gi_idx <- which(target$group == g)
  if (length(gi_idx) == 0) {
    message("Group ", g, ": no cases; skipped")
    next
  }
  gi <- target[gi_idx, , drop = FALSE]
  message("Group ", g, ": ", nrow(gi), " cases")

  # Tumor
  out_t <- run_pcod(gi$tumor_id)
  target$has_outlier_tumor[gi_idx[out_t]] <- 1L
  message(
    "  Tumor : ", length(out_t), " outlier(s)",
    if (length(out_t)) {
      paste0(" (", paste(gi$case_submitter_id[out_t],
        collapse = ", "
      ), ")")
    } else {
      ""
    }
  )

  # Normal
  out_n <- run_pcod(gi$normal_id)
  target$has_outlier_normal[gi_idx[out_n]] <- 1L
  message(
    "  Normal: ", length(out_n), " outlier(s)",
    if (length(out_n)) {
      paste0(" (", paste(gi$case_submitter_id[out_n],
        collapse = ", "
      ), ")")
    } else {
      ""
    }
  )
}

# --- Assemble and save -----------------------------------------------------
out_tbl <- target[, c(
  "case_submitter_id", "group", "tumor_purity",
  "has_outlier_tumor", "has_outlier_normal"
), drop = FALSE]
out_tbl <- out_tbl[order(out_tbl$group, out_tbl$case_submitter_id), , drop = FALSE]
rownames(out_tbl) <- NULL

message("Cases with a tumor outlier : ", sum(out_tbl$has_outlier_tumor))
message("Cases with a normal outlier: ", sum(out_tbl$has_outlier_normal))

out <- file.path(paths$processed, "thyr_case_outliers.rds")
saveRDS(out_tbl, out)
message(
  "Saved: ", out, " (", nrow(out_tbl), " cases x ",
  ncol(out_tbl), " columns)"
)
