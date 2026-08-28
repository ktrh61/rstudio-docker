# tab_cohort_composition.R
# Summarize case counts by driver classification crossed with Assigned Share
# (AS) band and sample-pair status. This is a summary only: it does not select
# analysis targets or define analysis groups; driver classification is taken
# verbatim from the Designated_* columns. Not part of the analysis pipeline; no
# inference is performed. The observed composition is written as a
# publication-formatting CSV as well as printed.
#
# Output: output/tables/supp_tab_cohort_composition.csv (Table S1, wide:
#         one row per category, one column per AS band, cell = all (paired);
#         unevaluated band cells outside the IREP reference groups are shown as
#         an em dash, never 0) and *_long.csv (provenance grid, 225 rows).
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

# --- Publication table (wide; Table S1) --------------------------------------
# One row per category, one column per AS band, cell = all cases (paired cases).
# AS was calculated only for the IREP reference set: exposed cases in the
# RET-fusion and BRAF-mutation designated driver groups. For categories outside
# those groups the Low/Mid/High cells were never evaluated and are shown as
# "—", not 0 -- a printed 0 is reserved for an evaluated band with no case
# (researcher decision 2026-08-28: no cell may read as a count that was never
# taken). Column labels follow the manuscript band names.
ref_groups <- c("Fusion.RET", "Mut.BRAF")
na_fill <- function(v) {
  v <- as.character(v)
  v[is.na(v) | v == ""] <- "na"
  v
}
cat_ref <- unique(rbindlist(list(
  dt[, .(level = "Group", category = na_fill(Designated_DriverGroup),
         in_ref = Designated_DriverGroup %in% ref_groups)],
  dt[, .(level = "Driver", category = na_fill(Designated_Driver),
         in_ref = Designated_DriverGroup %in% ref_groups)]
)))
if (anyDuplicated(cat_ref[, .(level, category)]) > 0) {
  stop("A category spans reference and non-reference driver groups.")
}
band_pub <- c(
  non_exposed = "dose-zero", no_reference = "exposed, AS not calculated",
  `(0,33.3)` = "Low-AS", `[33.3,66.6)` = "Mid-AS", `[66.6,100]` = "High-AS"
)
pub <- merge(summary_wide, cat_ref, by = c("level", "category"))
pub[, evaluated := in_ref | band %in% c("non_exposed", "no_reference")]
# Semantics check: an unevaluated cell must carry no case.
if (any(pub[evaluated == FALSE, all] != 0)) {
  print(pub[evaluated == FALSE & all != 0])
  stop("A case carries an AS band outside the IREP reference groups.")
}
pub[, cell := fifelse(evaluated, sprintf("%d (%d)", all, paired), "—")]
pub[, band_lab := factor(band_pub[as.character(band)], levels = unname(band_pub))]
pub_wide <- dcast(pub, level + category ~ band_lab, value.var = "cell")
tot <- summary_wide[, .(cat_n = sum(all),
                        total = sprintf("%d (%d)", sum(all), sum(paired))),
                    by = .(level, category)]
pub_wide <- merge(pub_wide, tot, by = c("level", "category"))
setorder(pub_wide, level, -cat_n, category)
pub_wide[, cat_n := NULL]
setcolorder(pub_wide, c("level", "category", unname(band_pub), "total"))
message("Reference-set cases (AS calculated): ",
        sum(dt$assigned_share_status == "irep"), " = ",
        paste(sprintf("%s %d", ref_groups,
          vapply(ref_groups, function(g)
            sum(dt$assigned_share_status == "irep" & dt$Designated_DriverGroup == g),
            integer(1))), collapse = " + "))

out_path <- file.path(out_dir, "supp_tab_cohort_composition.csv")
utils::write.csv(pub_wide, out_path, row.names = FALSE)
message("Saved (Table S1, wide): ", out_path)
