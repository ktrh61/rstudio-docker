# 042_summarize_driver_as.R
# Summarize case counts by driver classification crossed with Assigned Share
# (AS) band and sample-pair status. This is a summarization layer only: it does
# NOT select analysis targets or define analysis groups (no BRAF V600E
# restriction, no co-mutation exclusion, no RET grepl selection); driver
# classification is taken verbatim from the Designated_* columns, and target
# selection is deferred to downstream, which specifies cases statistically from
# the clinical master, the AS table, and the SE.
#
# Input : processed/thyr_clinical_master.rds       (from 040; master key case_id)
#         processed/thyr_case_assigned_share.rds    (from 041; AS per case)
#         processed/thyr_se_raw.rds                 (from 021; for pair status)
# Output: output/driver_as_summary.csv              (long-format summary)
#
# AS band (4 categories, covering all cases):
#   non_exposed  status == "missing_input"; AS not computed.
#   (0,33.3)     computed AS (percent).
#   [33.3,66.6)  computed AS.
#   [66.6,100]   computed AS.
#
# Pair status: derived on the fly from SE colData, not stored as a flag. A case
# is "paired" when it carries both a Primary Tumor and a Solid Tissue Normal
# sample (sample_type exact match). SE case ids in the TSS-augmented YQ form
# (REBC-YQ-<id>) are normalized to the master key form (REBC-<id>) on the SE
# side only, in memory, matching 040. Cases absent from the SE are unpaired.
#
# Classification axes: two levels (na kept as its own category).
#   Group  = Designated_DriverGroup
#   Driver = Designated_Driver
# Driver resolves specific fusions (e.g. ETV6-NTRK3) that Group collapses;
# Group is retained alongside so low-frequency Driver categories remain
# interpretable in aggregate.
#
# Output columns: level, category, pair_class, band, n
#   level      "Group" | "Driver"
#   pair_class "all" | "paired" | "unpaired"   (all = paired + unpaired)
#   band       one of the four AS bands above
# Only non-empty cells are emitted (following 040's driver_summary style).

source("setup.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
})

# --- Load inputs -------------------------------------------------------------
master_path <- file.path(paths$processed, "thyr_clinical_master.rds")
as_path <- file.path(paths$processed, "thyr_case_assigned_share.rds")
se_path <- file.path(paths$processed, "thyr_se_raw.rds")

for (p in c(master_path, as_path, se_path)) {
  if (!file.exists(p)) stop("Required input not found: ", p)
}

master <- readRDS(master_path)
as_dt <- readRDS(as_path)
se <- readRDS(se_path)
setDT(master)
setDT(as_dt)

message("Master: ", nrow(master), " cases ; AS table: ", nrow(as_dt), " cases")

# --- Join AS onto the master and assign the AS band --------------------------
dt <- merge(
  master[, .(case_id, Designated_DriverGroup, Designated_Driver)],
  as_dt[, .(case_id, assigned_share_approx, assigned_share_approx_status)],
  by = "case_id", all.x = TRUE
)
if (nrow(dt) != nrow(master)) {
  stop("Row count changed after AS join: ", nrow(dt), " vs ", nrow(master))
}

assign_band <- function(v, st) {
  b <- rep(NA_character_, length(v))
  b[st == "missing_input"] <- "non_exposed"
  ok <- st == "computed"
  b[ok & v > 0 & v < 33.3] <- "(0,33.3)"
  b[ok & v >= 33.3 & v < 66.6] <- "[33.3,66.6)"
  b[ok & v >= 66.6] <- "[66.6,100]"
  b
}
dt[, band := assign_band(assigned_share_approx, assigned_share_approx_status)]

band_levels <- c("non_exposed", "(0,33.3)", "[33.3,66.6)", "[66.6,100]")
if (any(is.na(dt$band))) {
  bad <- dt[is.na(band), unique(assigned_share_approx_status)]
  stop("Unassigned AS band for status: ", paste(bad, collapse = ", "))
}
dt[, band := factor(band, levels = band_levels)]

message(
  "AS band counts: ",
  paste(sprintf("%s=%d", band_levels, as.integer(table(dt$band))),
    collapse = " ; "
  )
)

# --- Derive pair status from SE (on the fly; YQ normalized on the SE side) --
cd <- as.data.frame(colData(se))
normalize_case_id <- function(x) sub("^REBC-YQ-", "REBC-", as.character(x))
cd$case_norm <- normalize_case_id(cd$case_submitter_id)

tumor_type <- "Primary Tumor"
normal_type <- "Solid Tissue Normal"
present <- unique(cd$sample_type)
if (!(tumor_type %in% present)) stop("SE lacks sample_type: ", tumor_type)
if (!(normal_type %in% present)) stop("SE lacks sample_type: ", normal_type)

has_tumor <- unique(cd$case_norm[cd$sample_type == tumor_type])
has_normal <- unique(cd$case_norm[cd$sample_type == normal_type])
paired_cases <- intersect(has_tumor, has_normal)

dt[, pair_class := ifelse(case_id %in% paired_cases, "paired", "unpaired")]
message(
  "Paired cases: ", sum(dt$pair_class == "paired"),
  " / unpaired: ", sum(dt$pair_class == "unpaired")
)

# --- Build the long-format summary -------------------------------------------
# For a given classification column and a subset of cases, count (category x
# band) cells.
count_level <- function(data, col, level_name) {
  v <- as.character(data[[col]])
  v[is.na(v) | v == ""] <- "na"
  tb <- data.table(category = v, band = data$band)
  out <- tb[, .(n = .N), by = .(category, band)]
  out[, level := level_name]
  out[]
}

# Build the (level, category, band, n) counts for one pair_class subset.
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

# Column order and sort: level, then within level by descending total n of the
# category, then band order, mirroring 040's descending-count style.
cat_tot <- summary_long[, .(cat_n = sum(n)), by = .(level, category)]
summary_long <- merge(summary_long, cat_tot, by = c("level", "category"))
setcolorder(summary_long, c("level", "category", "pair_class", "band", "n"))
setorder(summary_long, level, -cat_n, category, pair_class, band)
summary_long[, cat_n := NULL]

# --- Integrity checks ---------------------------------------------------------
# all == paired + unpaired for every (level, category, band) cell.
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

# Each level's "all" total must equal the full case count.
lvl_tot <- summary_long[pair_class == "all", .(n = sum(n)), by = level]
if (any(lvl_tot$n != nrow(master))) {
  print(lvl_tot)
  stop("A level's all-total does not equal ", nrow(master))
}
message(
  "Integrity checks passed (all = paired + unpaired ; level totals = ",
  nrow(master), ")"
)

# --- Save ---------------------------------------------------------------------
out_path <- file.path(paths$output, "driver_as_summary.csv")
fwrite(summary_long, out_path, na = "NA")
message("Saved: output/driver_as_summary.csv (", nrow(summary_long), " rows)")

message(
  "042 complete. Driver classification is summarized against AS band and pair ",
  "status; no target selection or grouping is performed here."
)
