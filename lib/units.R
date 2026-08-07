# Analysis-unit conventions: canonical unit order and the Sporadic/High arm
# split of a unit's group labels. Single definition of the group-label parsing
# that 410/420/510/520 and the figures previously each re-implemented.

UNIT_ORDER <- c("R_Tumor", "R_Normal", "B_Tumor", "B_Normal")

# Split a unit's per-sample group labels (e.g. "R_Sporadic"/"R_High") into the
# two arms, as integer indices. Stops if any sample falls outside the two arms
# (this cannot happen for units built by 310; the check guards silent drift).
unit_arms <- function(group, unit = "") {
  group <- as.character(group)
  sporadic <- which(grepl("Sporadic", group, fixed = TRUE))
  high <- which(grepl("High", group, fixed = TRUE))
  if (length(sporadic) + length(high) != length(group)) {
    stop("Unit ", unit, " has samples outside Sporadic / High.")
  }
  list(sporadic = sporadic, high = high)
}

# The SE's single count assay, verified. 120 stores exactly one assay named by
# the data-driven strand selection ("stranded_second" on this dataset); callers
# must not hardcode that name -- this accessor stops if the invariant breaks.
se_single_assay <- function(se) {
  an <- SummarizedExperiment::assayNames(se)
  if (length(an) != 1L) {
    stop("Expected a single assay in the SE, found: ", paste(an, collapse = ", "))
  }
  SummarizedExperiment::assay(se, 1L)
}
