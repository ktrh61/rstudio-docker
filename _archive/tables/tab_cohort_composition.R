# tab_cohort_composition.R
# Summarize case counts by driver classification crossed with Assigned Share
# (AS) band and sample-pair status. This is a summary only: it does not select
# analysis targets or define analysis groups; driver classification is taken
# verbatim from the Designated_* columns. Not part of the analysis pipeline; no
# inference is performed. The observed composition is written as a
# publication-formatting CSV as well as printed.
#
# Output: output/tables/supp_tab_cohort_composition.csv (Table S1: rows = analysis
#         strata and the cases outside them, columns = dose-zero / AS bands /
#         exposed outside strata / total, cell = all (paired); AS shown only for
#         the strata) and *_long.csv (designated-driver x band provenance grid).
# Input : processed/thyr_clinical.rds               (from 030; key REBC_ID)
#         processed/thyr_case_assigned_share.rds     (from 130; AS per case)
#         processed/thyr_se_raw.rds                  (from 120; for pair status)
#
# AS band (5 categories, covering all cases):
#   non_exposed  dose == 0 (status "not_required_sporadic"); AS undefined.
#   no_reference exposed but outside the IREP reference (status
#                "not_in_reference" -- driver-unclassified cases; v2 B.11).
#   (0,33.3)     IREP AS (percent).
#   [33.3,66.6)  IREP AS.
#   [66.6,100]   IREP AS.
#
# Pair status: derived on the fly from SE colData. A case is "paired" when it
# carries both a Primary Tumor and a Solid Tissue Normal sample (sample_type
# exact match). SE case ids are already YQ-normalized to the REBC form (done in
# 110), so they match the REBC_ID key directly. Cases absent from the SE are
# unpaired.
#
# Classification axes: two levels (na kept as its own category).
#   Group  = Designated_DriverGroup
#   Driver = Designated_Driver
# Driver resolves specific fusions that Group collapses; Group is retained
# alongside so low-frequency Driver categories remain interpretable in aggregate.
# (Plain per-category counts, as an earlier Group/Driver summary produced, are the
# band/pair marginals of this table.)

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
})

# --- Load inputs -------------------------------------------------------------
clinical_path <- file.path(paths$processed, "thyr_clinical.rds")
as_path <- file.path(paths$processed, "thyr_case_assigned_share.rds")
se_path <- file.path(paths$processed, "thyr_se_raw.rds")

for (p in c(clinical_path, as_path, se_path)) {
  if (!file.exists(p)) stop("Required input not found: ", p)
}

clinical <- readRDS(clinical_path)
as_dt <- readRDS(as_path)
se <- readRDS(se_path)
setDT(clinical)
setDT(as_dt)

message("Clinical: ", nrow(clinical), " cases ; AS table: ", nrow(as_dt), " cases")

# --- Join AS onto the clinical table and assign the AS band ------------------
dt <- merge(
  clinical[, .(REBC_ID, Designated_DriverGroup, Designated_Driver)],
  as_dt[, .(REBC_ID, assigned_share, assigned_share_status)],
  by = "REBC_ID", all.x = TRUE
)
if (nrow(dt) != nrow(clinical)) {
  stop("Row count changed after AS join: ", nrow(dt), " vs ", nrow(clinical))
}

assign_band <- function(v, st) {
  b <- rep(NA_character_, length(v))
  b[st == "not_required_sporadic"] <- "non_exposed"
  b[st == "not_in_reference"] <- "no_reference"
  ok <- st == "irep"
  b[ok & v > 0 & v < 33.3] <- "(0,33.3)"
  b[ok & v >= 33.3 & v < 66.6] <- "[33.3,66.6)"
  b[ok & v >= 66.6] <- "[66.6,100]"
  b
}
dt[, band := assign_band(assigned_share, assigned_share_status)]

band_levels <- c(
  "non_exposed", "no_reference", "(0,33.3)", "[33.3,66.6)", "[66.6,100]"
)
if (any(is.na(dt$band))) {
  bad <- dt[is.na(band), unique(assigned_share_status)]
  stop("Unassigned AS band for status: ", paste(bad, collapse = ", "))
}
dt[, band := factor(band, levels = band_levels)]

message(
  "AS band counts: ",
  paste(sprintf("%s=%d", band_levels, as.integer(table(dt$band))),
    collapse = " ; "
  )
)

# --- Derive pair status from SE (on the fly; ids already YQ-normalized) ------
cd <- as.data.frame(colData(se))

tumor_type <- "Primary Tumor"
normal_type <- "Solid Tissue Normal"
present <- unique(cd$sample_type)
if (!(tumor_type %in% present)) stop("SE lacks sample_type: ", tumor_type)
if (!(normal_type %in% present)) stop("SE lacks sample_type: ", normal_type)

has_tumor <- unique(cd$case_submitter_id[cd$sample_type == tumor_type])
has_normal <- unique(cd$case_submitter_id[cd$sample_type == normal_type])
paired_cases <- intersect(has_tumor, has_normal)

dt[, pair_class := ifelse(REBC_ID %in% paired_cases, "paired", "unpaired")]
message(
  "Paired cases: ", sum(dt$pair_class == "paired"),
  " / unpaired: ", sum(dt$pair_class == "unpaired")
)

# --- Build the long-format summary -------------------------------------------
count_level <- function(data, col, level_name) {
  v <- as.character(data[[col]])
  v[is.na(v) | v == ""] <- "na"
  tb <- data.table(category = v, band = data$band)
  out <- tb[, .(n = .N), by = .(category, band)]
  out[, level := level_name]
  out[]
}

summarize_subset <- function(data, pair_label) {
  parts <- rbindlist(list(
    count_level(data, "Designated_DriverGroup", "Group"),
    count_level(data, "Designated_Driver", "Driver")
  ))
  parts[, pair_class := pair_label]
  parts[]
}

summary_long <- rbindlist(list(
  summarize_subset(dt, "all"),
  summarize_subset(dt[pair_class == "paired"], "paired"),
  summarize_subset(dt[pair_class == "unpaired"], "unpaired")
))

cat_tot <- summary_long[, .(cat_n = sum(n)), by = .(level, category)]
summary_long <- merge(summary_long, cat_tot, by = c("level", "category"))
setcolorder(summary_long, c("level", "category", "pair_class", "band", "n"))
setorder(summary_long, level, -cat_n, category, pair_class, band)
summary_long[, cat_n := NULL]

# --- Integrity checks ---------------------------------------------------------
wide <- dcast(summary_long, level + category + band ~ pair_class,
  value.var = "n", fill = 0
)
if (!all(c("all", "paired", "unpaired") %in% names(wide))) {
  stop("Missing a pair_class in the summary")
}
mism <- wide[all != paired + unpaired]
if (nrow(mism) > 0) {
  print(mism)
  stop("all != paired + unpaired in ", nrow(mism), " cells")
}

lvl_tot <- summary_long[pair_class == "all", .(n = sum(n)), by = level]
if (any(lvl_tot$n != nrow(clinical))) {
  print(lvl_tot)
  stop("A level's all-total does not equal ", nrow(clinical))
}
message(
  "Integrity checks passed (all = paired + unpaired ; level totals = ",
  nrow(clinical), ")"
)

# --- Print the long summary (provenance log) ---------------------------------
message("Driver x AS band x pair summary (", nrow(summary_long), " rows):")
print(summary_long, nrow = nrow(summary_long))

# Complete the category-by-band grid so that absent combinations are counted
# as zero rather than omitted (long form, kept for provenance).
categories <- unique(summary_long[, .(level, category)])
grid <- categories[, .(band = factor(band_levels, levels = band_levels)),
  by = .(level, category)
]
summary_wide <- dcast(summary_long, level + category + band ~ pair_class,
  value.var = "n", fill = 0
)
summary_wide <- merge(grid, summary_wide,
  by = c("level", "category", "band"), all.x = TRUE
)
for (nm in c("all", "paired", "unpaired")) {
  set(summary_wide, which(is.na(summary_wide[[nm]])), nm, 0L)
}
setcolorder(summary_wide,
  c("level", "category", "band", "all", "paired", "unpaired")
)
setorder(summary_wide, level, category, band)

out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
long_path <- file.path(out_dir, "supp_tab_cohort_composition_long.csv")
utils::write.csv(summary_wide, long_path, row.names = FALSE)
message("Saved (long, provenance): ", long_path)

# --- Publication table (Table S1; analysis-frame rows) -----------------------
# Rows follow the analysis design rather than the source cohort's driver
# catalogue: the two analysis strata (RET fusion-positive by partner; BRAF
# V600E without a co-occurring candidate driver mutation) carry AS-band cells,
# and every case outside the strata is counted in one exposed column without
# AS. The versioned IREP file also holds values for exposed cases outside the
# strata; they enter no analysis and are not displayed (researcher decision
# 2026-08-29: the table shows exactly the AS the analysis used). An em dash
# marks a cell that does not apply; 0 is a counted band with no case.
# Cell = all cases (paired cases).
source(file.path(paths$root, "lib", "cohort_design.R"))
design <- classify_driver(clinical)
dt <- merge(dt,
  data.table(REBC_ID = design$case_submitter_id, driver = design$driver),
  by = "REBC_ID", all.x = TRUE
)
if (nrow(dt) != nrow(clinical)) stop("Row count changed after driver join.")
drv <- as.character(dt$Designated_Driver)
drv[is.na(drv) | drv == ""] <- "na"
is_ret <- !is.na(dt$driver) & dt$driver == "RET"
is_braf <- !is.na(dt$driver) & dt$driver == "BRAF"
dt[, in_stratum := is_ret | is_braf]
dt[, exposed := assigned_share_status != "not_required_sporadic"]
# Guard (mirrors 140): every exposed stratum case carries an IREP value.
if (any(dt[in_stratum & exposed, assigned_share_status] != "irep")) {
  stop("An exposed stratum case lacks an IREP value.")
}
dt[, row := fcase(
  is_ret & drv == "CCDC6-RET", "CCDC6-RET",
  is_ret & drv == "NCOA4-RET", "NCOA4-RET",
  is_ret, "Other RET fusion partner",
  is_braf, "BRAF V600E stratum (no co-occurring candidate driver mutation)",
  drv == "BRAF.MutV600E", "BRAF V600E with a co-occurring candidate driver mutation",
  drv == "na", "No designated driver",
  default = "Other designated driver"
)]
band_pub <- c(`(0,33.3)` = "Low-AS", `[33.3,66.6)` = "Mid-AS", `[66.6,100]` = "High-AS")
cell <- function(d) sprintf("%d (%d)", nrow(d), sum(d$pair_class == "paired"))
make_row <- function(d, label) {
  any_str <- any(d$in_stratum)
  any_out <- any(!d$in_stratum)
  out <- list(category = label, `dose-zero` = cell(d[exposed == FALSE]))
  for (b in names(band_pub)) {
    out[[band_pub[[b]]]] <- if (any_str) {
      cell(d[in_stratum & exposed & as.character(band) == b])
    } else "—"
  }
  out[["exposed, outside strata"]] <- if (any_out) cell(d[!in_stratum & exposed]) else "—"
  out[["total"]] <- cell(d)
  as.data.table(out)
}
pub <- rbindlist(list(
  make_row(dt[row == "CCDC6-RET"], "CCDC6-RET"),
  make_row(dt[row == "NCOA4-RET"], "NCOA4-RET"),
  make_row(dt[row == "Other RET fusion partner"], "Other RET fusion partner"),
  make_row(dt[is_ret], "RET fusion-positive stratum, subtotal"),
  make_row(dt[is_braf], "BRAF V600E stratum (no co-occurring candidate driver mutation)"),
  make_row(dt[row == "BRAF V600E with a co-occurring candidate driver mutation"],
           "BRAF V600E with a co-occurring candidate driver mutation"),
  make_row(dt[row == "Other designated driver"], "Other designated driver"),
  make_row(dt[row == "No designated driver"], "No designated driver"),
  make_row(dt, "All cases")
))
message(sprintf(
  "Strata: RET %d | BRAF %d (of %d designated V600E) | exposed outside strata %d | unused IREP values outside strata %d",
  sum(is_ret), sum(is_braf), sum(drv == "BRAF.MutV600E"),
  sum(!dt$in_stratum & dt$exposed),
  sum(!dt$in_stratum & dt$assigned_share_status == "irep")))

out_path <- file.path(out_dir, "supp_tab_cohort_composition.csv")
utils::write.csv(pub, out_path, row.names = FALSE)
message("Saved (Table S1): ", out_path)
