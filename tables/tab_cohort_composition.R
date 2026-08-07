# tab_cohort_composition.R
# Summarize case counts by driver classification crossed with Assigned Share
# (AS) band and sample-pair status. This is a summary only: it does not select
# analysis targets or define analysis groups; driver classification is taken
# verbatim from the Designated_* columns. Not part of the analysis pipeline; no
# files are written (results are printed).
#
# Input : processed/thyr_clinical.rds               (from 030; key REBC_ID)
#         processed/thyr_case_assigned_share.rds     (from 130; AS per case)
#         processed/thyr_se_raw.rds                  (from 120; for pair status)
#
# AS band (4 categories, covering all cases):
#   non_exposed  status == "missing_input"; AS not computed.
#   (0,33.3)     computed AS (percent).
#   [33.3,66.6)  computed AS.
#   [66.6,100]   computed AS.
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
  as_dt[, .(REBC_ID, assigned_share_approx, assigned_share_approx_status)],
  by = "REBC_ID", all.x = TRUE
)
if (nrow(dt) != nrow(clinical)) {
  stop("Row count changed after AS join: ", nrow(dt), " vs ", nrow(clinical))
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

# --- Print the summary (no file output) --------------------------------------
message("Driver x AS band x pair summary (", nrow(summary_long), " rows):")
print(summary_long, nrow = nrow(summary_long))
