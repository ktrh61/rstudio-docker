# 010_download_expression.R
# Download REBC-THYR Gene Expression Quantification files (STAR - Counts) from GDC.
# Input : GDC files query (REBC-THYR, RNA-Seq, STAR - Counts, open access)
# Output: raw/expression/<file_id>/<file>.rna_seq.augmented_star_gene_counts.tsv
#         meta/manifest_gene_counts_<timestamp>.tsv       (download manifest)
#         meta/manifest_retry_<timestamp>.tsv             (only if files remain missing)
#
# Each file is verified against the manifest md5sum after download. Partial
# failures are allowed and collected into a retry manifest; re-running the
# script downloads only what is still missing (md5-matching files are skipped).
# The script completes successfully only when every manifest file_id is present
# and md5-verified; otherwise it writes a retry manifest and exits non-zero.

source("setup.R")

library(GenomicDataCommons)

# Per-connection safety valve for individual downloads (not a total-runtime cap).
options(timeout = 180)

expr_dir <- file.path(paths$raw, "expression")
dir.create(expr_dir, recursive = TRUE, showWarnings = FALSE)

# --- Build manifest --------------------------------------------------------
q <- files() |>
  GenomicDataCommons::filter(
    ~ cases.project.project_id == "REBC-THYR" &
      data_category == "Transcriptome Profiling" &
      experimental_strategy == "RNA-Seq" &
      data_type == "Gene Expression Quantification" &
      analysis.workflow_type == "STAR - Counts" &
      access == "open"
  )

manifest <- manifest(q)

if (nrow(manifest) == 0) {
  stop("No files found in manifest")
}
if (!all(c("id", "md5sum", "file_name") %in% names(manifest))) {
  stop("Manifest missing required columns (id, md5sum, file_name)")
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
manifest_file <- file.path(
  paths$meta, paste0("manifest_gene_counts_", timestamp, ".tsv")
)
data.table::fwrite(manifest, manifest_file, sep = "\t", quote = FALSE)

message("Manifest files: ", nrow(manifest))

# --- Helper: is a file_id already present and md5-valid? -------------------
# Expected layout: raw/expression/<file_id>/<file_name>. A file_id counts as
# present only when the expected file exists and its md5 matches the manifest.
is_valid <- function(file_id, file_name, md5_expected) {
  target <- file.path(expr_dir, file_id, file_name)
  if (!file.exists(target)) {
    return(FALSE)
  }
  tools::md5sum(target)[[1]] == md5_expected
}

# --- Download loop ---------------------------------------------------------
file_ids <- manifest$id
n_total <- length(file_ids)

failed <- character()

message("Downloading up to ", n_total, " files to ", expr_dir, " ...")

for (i in seq_len(n_total)) {
  file_id <- file_ids[i]
  file_name <- manifest$file_name[i]
  md5_expected <- manifest$md5sum[i]

  if (is_valid(file_id, file_name, md5_expected)) {
    next
  }

  result <- try(
    gdcdata(file_id, progress = FALSE),
    silent = TRUE
  )

  if (inherits(result, "try-error")) {
    failed <- c(failed, file_id)
  } else {
    cached_path <- as.character(result)[1]
    target_dir <- file.path(expr_dir, file_id)
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    target_file <- file.path(target_dir, basename(cached_path))
    copied <- file.copy(cached_path, target_file, overwrite = TRUE)

    if (!copied ||
      tools::md5sum(target_file)[[1]] != md5_expected) {
      if (file.exists(target_file)) {
        file.remove(target_file)
      }
      failed <- c(failed, file_id)
    }
  }

  if (i %% 100 == 0) {
    message(i, " / ", n_total, " processed")
  }
}

# --- Completeness check ----------------------------------------------------
# Re-verify every manifest file_id against md5, independent of this run's
# success/failure bookkeeping, so a prior partial run is also accounted for.
valid <- vapply(
  seq_len(n_total),
  function(i) {
    is_valid(
      manifest$id[i], manifest$file_name[i], manifest$md5sum[i]
    )
  },
  logical(1)
)
missing_ids <- manifest$id[!valid]

message("Verified present: ", sum(valid), " / ", n_total)

if (length(missing_ids) > 0) {
  retry_manifest <- manifest[manifest$id %in% missing_ids, ]
  retry_file <- file.path(
    paths$meta, paste0("manifest_retry_", timestamp, ".tsv")
  )
  data.table::fwrite(retry_manifest, retry_file, sep = "\t", quote = FALSE)
  message("Missing: ", length(missing_ids), " (retry manifest written)")
  quit(status = 1L)
}

message("All ", n_total, " files present and md5-verified")
