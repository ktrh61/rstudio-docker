# 001_prepare_clinical.R
# Prepare the clinical table from Science Data S1 (a TSV named .txt) for
# downstream use: read all columns, coerce missing markers to NA, and type each
# column (numeric where the column parses as numeric, character otherwise). No
# columns are dropped and no values are edited; candidate driver columns are
# kept verbatim. Column pruning and any gene-symbol normalization are deferred
# downstream.
#
# Input : raw/clinical/abg2538-data-s1.txt   (Science Data S1, 440 cases; manual)
# Output: processed/thyr_clinical.rds        (prepared clinical table)
#
# Missing values: the source uses "NA", "na", and empty cells as missing; all are
# read as NA. Missing is held as R's NA (detectable by is.na), not as a string
# such as "NA" nor as "".
#
# Typing: all columns are read as character first, then a column is coerced to
# numeric only when every non-missing value parses as a number. Values that fail
# to parse in a coerced column become NA and are reported as a QC signal.

source("setup.R")

# --- Load: all columns as character, source missing markers to NA ------------
s1_path <- file.path(paths$raw, "clinical", "abg2538-data-s1.txt")
if (!file.exists(s1_path)) {
  stop("Clinical S1 not found: ", s1_path)
}
clinical <- data.table::fread(
  s1_path,
  na.strings = c("NA", "na", ""),
  colClasses = "character"
)
message("S1 loaded: ", nrow(clinical), " cases x ", ncol(clinical), " variables")

# --- Type each column: numeric where every non-missing value parses ----------
is_numeric_col <- function(x) {
  nonmiss <- x[!is.na(x)]
  if (length(nonmiss) == 0) {
    return(FALSE)
  } # all-missing: no basis to type as numeric
  parsed <- suppressWarnings(as.numeric(nonmiss))
  !any(is.na(parsed))
}

num_cols <- names(clinical)[vapply(clinical, is_numeric_col, logical(1))]
for (col in num_cols) {
  before_na <- is.na(clinical[[col]])
  coerced <- suppressWarnings(as.numeric(clinical[[col]]))
  new_na <- is.na(coerced) & !before_na
  if (any(new_na)) {
    message(
      "QC: column ", col, " coerced to numeric but ", sum(new_na),
      " non-missing value(s) failed to parse and became NA"
    )
  }
  data.table::set(clinical, j = col, value = coerced)
}
message("Numeric columns: ", length(num_cols), " / ", ncol(clinical))

# --- Save --------------------------------------------------------------------
saveRDS(clinical, file.path(paths$processed, "thyr_clinical.rds"))
message(
  "Saved: processed/thyr_clinical.rds (",
  nrow(clinical), " cases x ", ncol(clinical), " columns)"
)

message("001 complete. Clinical table prepared; all columns kept verbatim.")
