# 040_build_clinical_master.R
# Build the clinical master (analysis normal form) from Science Data S1 and
# summarize driver information across all cases. This is a summarization layer:
# it records what driver information each case carries and does NOT decide which
# cases enter the analysis. Analysis-target selection (e.g. BRAF single vs
# co-mutation, RET fusion pooling) is deferred to downstream (041/042).
#
# Input : raw/clinical/abg2538-data-s1.txt   (Science Data S1, 440 cases; manual)
#         processed/thyr_se_raw.rds          (from 021; for case-presence pairing)
# Output: processed/thyr_clinical_master.rds (clinical master, analysis normal form)
#         output/thyr_clinical_master.csv    (human-readable full master)
#         output/driver_summary.csv          (case counts by driver classification)
#
# Case ID convention: S1 (REBC-<id>) is the full 440-case set and is the master
# key. SE case IDs from GDC carry a TSS-augmented barcode (REBC-YQ-<id>) for 12
# cases; these denote the same individuals. Pairing strips `REBC-YQ-` to `REBC-`
# on the SE side only, in memory. Neither S1 nor the SE file is overwritten; the
# SE key on disk keeps its YQ form. The master key is the S1 form (no YQ).
#
# Driver aggregation: RNA and WGS candidate columns each carry independent
# information (co-mutation partners, method-specific detections), so the four
# raw candidate columns are kept verbatim. An aggregated set column is derived
# per fusion/mutation family by removing parenthetical aliases (e.g. RET-PTC1),
# splitting comma-separated events, normalizing two known truncated gene symbols
# (PIK3C->PIK3CA, PIK3R->PIK3R1; researcher decision), and taking the RNA/WGS
# union. A mismatch flag records where the RNA and WGS element sets disagree.
# Raw columns are preserved alongside the aggregated columns.

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Column selection
# ---------------------------------------------------------------------------
# Driver columns retained from S1. Designated_DriverType is kept as the raw
# column; a Mut/Fusion split is also derived below (redundant with Pathway, but
# cheap and losslessly recoverable). WGS-derived bulk columns (SSV/SV/SCNA/
# signature/germline, several hundred columns) describe WGS samples that may not
# match the RNA-seq material and are not used; they are dropped here.
driver_cols <- c(
  "Designated_Driver",
  "Designated_DriverGroup",
  "Designated_DriverPathway",
  "Designated_DriverType",
  "DriverGroup_CoMutFig",
  "TwoHitTSCorRNAExpression_Driver",
  "RNA_CandidateDriverFusion",
  "WGS_CandidateDriverFusion",
  "RNA_CandidateDriverMutation",
  "WGS_CandidateDriverMutation"
)

# Clinical variables retained for downstream (POC calculation in 041 needs age
# and dose; demographics kept for description). This set is provisional and
# should be revisited once 041's POC input spec is fixed.
clinical_cols <- c(
  "SEX", "AGE_SURGERY", "AGE_EXPOSURE", "DOSE"
)

id_col <- "REBC_ID"

# ---------------------------------------------------------------------------
# Load Science Data S1
# ---------------------------------------------------------------------------
s1_path <- file.path(paths$raw, "clinical", "abg2538-data-s1.txt")
if (!file.exists(s1_path)) {
  stop("Clinical S1 not found: ", s1_path)
}
s1 <- fread(s1_path)
message("S1 loaded: ", nrow(s1), " cases x ", ncol(s1), " variables")

missing_cols <- setdiff(c(id_col, driver_cols, clinical_cols), names(s1))
if (length(missing_cols) > 0) {
  stop("S1 missing expected columns: ", paste(missing_cols, collapse = ", "))
}

# Master key column (S1 form, no YQ).
master <- s1[, c(id_col, clinical_cols, driver_cols), with = FALSE]
setnames(master, id_col, "case_id")

# ---------------------------------------------------------------------------
# Driver aggregation helpers
# ---------------------------------------------------------------------------
norm_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

# Known truncated gene symbols in S1 candidate columns (researcher decision:
# these are data-entry truncations, normalized for summarization; raw columns
# above keep the original strings).
gene_fixes <- c("PIK3C" = "PIK3CA", "PIK3R" = "PIK3R1")

apply_gene_fixes <- function(elems) {
  hit <- elems %in% names(gene_fixes)
  elems[hit] <- gene_fixes[elems[hit]]
  elems
}

# Turn a candidate string into a sorted, de-duplicated element set: strip
# parenthetical aliases, split on commas, trim, drop empties, normalize symbols.
to_set <- function(s) {
  s <- gsub("\\([^)]*\\)", "", s) # remove "(RET-PTC1)" style aliases
  el <- unlist(strsplit(s, ","))
  el <- trimws(el)
  el <- el[el != ""]
  el <- apply_gene_fixes(el)
  sort(unique(el))
}

# Aggregate one RNA/WGS pair of columns into a union string and a mismatch flag.
# Mismatch flag: "" when the two sets are identical (including both empty),
# "RNA_only" / "WGS_only" when one side is empty, "diff" when both non-empty and
# unequal.
aggregate_pair <- function(rna_vec, wgs_vec) {
  n <- length(rna_vec)
  union_str <- character(n)
  flag <- character(n)
  for (i in seq_len(n)) {
    sr <- to_set(rna_vec[i])
    sw <- to_set(wgs_vec[i])
    u <- sort(unique(c(sr, sw)))
    union_str[i] <- paste(u, collapse = ",")
    if (identical(sr, sw)) {
      flag[i] <- ""
    } else if (length(sr) == 0) {
      flag[i] <- "WGS_only"
    } else if (length(sw) == 0) {
      flag[i] <- "RNA_only"
    } else {
      flag[i] <- "diff"
    }
  }
  list(union = union_str, flag = flag)
}

# ---------------------------------------------------------------------------
# Build aggregated driver columns (raw columns preserved above)
# ---------------------------------------------------------------------------
fus <- aggregate_pair(
  norm_chr(master$RNA_CandidateDriverFusion),
  norm_chr(master$WGS_CandidateDriverFusion)
)
mut <- aggregate_pair(
  norm_chr(master$RNA_CandidateDriverMutation),
  norm_chr(master$WGS_CandidateDriverMutation)
)

master[, fusion_union := fus$union]
master[, fusion_mismatch := fus$flag]
master[, mutation_union := mut$union]
master[, mutation_mismatch := mut$flag]

# Mut/Fusion split derived from Designated_DriverType (raw column kept). Empty
# type (na cases) yields "".
master[, driver_class := norm_chr(Designated_DriverType)]

message(
  "Fusion mismatches (diff/one-sided): ",
  sum(fus$flag != ""), " / ", nrow(master)
)
message(
  "Mutation mismatches (diff/one-sided): ",
  sum(mut$flag != ""), " / ", nrow(master)
)

# ---------------------------------------------------------------------------
# Pair with SE for case presence (YQ normalization on the SE side, in memory)
# ---------------------------------------------------------------------------
se_path <- file.path(paths$processed, "thyr_se_raw.rds")
if (!file.exists(se_path)) {
  stop("SE not found: ", se_path, " (run 021 first)")
}
se <- readRDS(se_path)
se_cases_raw <- unique(as.character(colData(se)$case_submitter_id))

# Strip the TSS-augmented prefix so SE case IDs match the S1 master key.
normalize_case_id <- function(x) sub("^REBC-YQ-", "REBC-", x)
se_cases_norm <- normalize_case_id(se_cases_raw)

# Map normalized -> original SE form, to record which cases carried the YQ form.
se_map <- data.table(
  case_id     = se_cases_norm,
  se_case_id  = se_cases_raw
)
setkey(se_map, case_id)

master[, in_se := case_id %in% se_cases_norm]
master[, se_case_id := NA_character_]
master[se_map, se_case_id := i.se_case_id, on = "case_id"]

message(
  "Cases present in SE (after YQ normalization): ",
  sum(master$in_se), " / ", nrow(master)
)
n_yq <- sum(!is.na(master$se_case_id) &
  master$se_case_id != master$case_id)
message("SE cases carrying the YQ form: ", n_yq)

# ---------------------------------------------------------------------------
# Coarse summary: case counts by driver classification
# ---------------------------------------------------------------------------
# Three nested granularities, all covering the 440-case total with na shown as
# its own category: Group (16, primary), Driver (32, finest), Pathway (9,
# auxiliary). Counts are over all cases, before any analysis-target selection.
summarize_col <- function(col, level_name) {
  v <- norm_chr(master[[col]])
  v[v == ""] <- "na"
  tb <- as.data.table(table(v))
  setnames(tb, c("category", "n"))
  tb[, level := level_name]
  setcolorder(tb, c("level", "category", "n"))
  tb[order(-n)]
}

driver_summary <- rbindlist(list(
  summarize_col("Designated_DriverGroup", "Group"),
  summarize_col("Designated_Driver", "Driver"),
  summarize_col("Designated_DriverPathway", "Pathway")
))

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
saveRDS(master, file.path(paths$processed, "thyr_clinical_master.rds"))
message(
  "Saved master: processed/thyr_clinical_master.rds (",
  nrow(master), " cases x ", ncol(master), " columns)"
)

fwrite(master, file.path(paths$output, "thyr_clinical_master.csv"), na = "NA")
message("Saved readable master: output/thyr_clinical_master.csv")

fwrite(driver_summary, file.path(paths$output, "driver_summary.csv"), na = "NA")
message("Saved driver summary: output/driver_summary.csv")

message(
  "040 complete. Driver classification is recorded, not selected; ",
  "analysis-target judgment is deferred to downstream."
)
