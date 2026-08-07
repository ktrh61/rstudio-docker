# MSigDB collections for the gene-set level (420) and its null calibration
# (diagnostics/gsea_null_calibration.R). One definition so the calibration
# exercises exactly the families the analysis uses.
#
# Families: H, C2:CP (canonical pathways), C5:GO:BP, and C2:CGP:radiation --
# the radiation-related subset of C2:CGP, added as an independent exploratory
# family (reorg plan v2 D7, decided 2026-08-07). Its curation rule is fixed
# here, before the new inference ever touches real data: gs_name must match
# RADIATION_INCLUDE and not RADIATION_EXCLUDE. UV-specific sets are excluded
# because the exposure under study is ionizing radiation. The selected set
# names are returned (and saved by callers) as supplement material.
# Construct caveat, recorded ahead of results: these sets describe acute,
# mostly high-dose or in vitro radiation responses, not a decades-later
# exposure memory; a null reads "the trace does not resemble the acute
# response", a positive reads "it does", and both are informative.
#
# C6 and C7 stay excluded on relevance grounds (driver-conditioned design;
# immune signatures outside the hypothesis). NES and the collection-internal
# FDR are computed within each collection, so adding a family cannot dilute
# another family's q-values.

RADIATION_INCLUDE <- "IONIZING|IRRADIAT|RADIATION|(^|_)IR(_|$)"
RADIATION_EXCLUDE <- "(^|_)UV[ABC]?(_|$)|ULTRAVIOLET"

# Size window applied to every family, in 420 and in the null calibration
# alike (one definition so the calibration tests the analyzed universe).
GSEA_MIN_SET_SIZE <- 15L
GSEA_MAX_SET_SIZE <- 500L

GSEA_COLLECTION_SPECS <- list(
  "H" = list(collection = "H", subcollection = NA_character_),
  "C2:CP" = list(
    collection = "C2",
    subcollection = c(
      "CP:REACTOME", "CP:WIKIPATHWAYS", "CP:KEGG_MEDICUS",
      "CP:BIOCARTA", "CP:PID"
    )
  ),
  "C5:GO:BP" = list(collection = "C5", subcollection = "GO:BP"),
  "C2:CGP:radiation" = list(
    collection = "C2", subcollection = "CGP",
    name_include = RADIATION_INCLUDE, name_exclude = RADIATION_EXCLUDE
  )
)

# Returns list(gene_sets, curation): gene_sets is one named list per family
# (gs_name -> Ensembl members), curation records the radiation regexes and
# the set names they selected.
load_gene_set_collections <- function(species = "Homo sapiens") {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package 'msigdbr' is required to load the gene-set collections.")
  }
  collect <- function(spec) {
    frames <- lapply(spec$subcollection, function(sub) {
      if (is.na(sub)) {
        msigdbr::msigdbr(species = species, collection = spec$collection)
      } else {
        msigdbr::msigdbr(
          species = species, collection = spec$collection,
          subcollection = sub
        )
      }
    })
    frame <- do.call(rbind, frames)
    if (!is.null(spec$name_include)) {
      frame <- frame[
        grepl(spec$name_include, frame$gs_name) &
          !grepl(spec$name_exclude, frame$gs_name), ,
        drop = FALSE
      ]
    }
    split(frame$ensembl_gene, frame$gs_name)
  }
  gene_sets <- lapply(GSEA_COLLECTION_SPECS, collect)
  list(
    gene_sets = gene_sets,
    curation = list(
      radiation_include = RADIATION_INCLUDE,
      radiation_exclude = RADIATION_EXCLUDE,
      radiation_sets = names(gene_sets[["C2:CGP:radiation"]])
    )
  )
}
