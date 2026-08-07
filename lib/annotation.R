# Gene-identifier helpers: Ensembl version stripping and gene_id -> gene_name
# mapping from a SummarizedExperiment's rowData.

strip_ensembl_version <- function(x) {
  sub("\\..*$", "", x)
}

# Named vector gene_id -> gene_name. strip_version = TRUE keys the map by
# unversioned Ensembl IDs (for callers that strip versions first, e.g. 510);
# FALSE keys by the versioned rownames (the figures).
gene_name_map <- function(se, strip_version = FALSE) {
  ids <- rownames(se)
  if (strip_version) ids <- strip_ensembl_version(ids)
  setNames(as.data.frame(SummarizedExperiment::rowData(se))$gene_name, ids)
}

# Display label with fallback: the gene name where present, else the
# version-stripped ID. (510 deliberately does NOT use the fallback; it keeps
# NA for unnamed genes.)
gene_label <- function(ids, name_of) {
  nm <- name_of[ids]
  ifelse(is.na(nm) | nm == "", strip_ensembl_version(ids), unname(nm))
}
