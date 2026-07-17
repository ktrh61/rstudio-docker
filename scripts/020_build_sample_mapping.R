# 020_build_sample_mapping.R
# Build a file-sample mapping for REBC-THYR STAR-Counts files by querying the
# GDC API for case- and sample-level metadata.
# Input : meta/manifest_gene_counts_<timestamp>.tsv   (from 010, latest is used)
# Output: processed/file_sample_mapping.rds           (consumed by 021)
#
# The manifest lists one row per expression file (file_id). For each file this
# script resolves the owning case and biospecimen sample via the GDC API and
# records the identifiers 021 needs to label SummarizedExperiment columns. API
# queries are sequential (rate limits); result parsing is parallelised.

source("setup.R")

library(GenomicDataCommons)
library(parallel)

# --- Load manifest ---------------------------------------------------------
# 010 writes manifests to meta/. Use the most recent one.
manifest_files <- list.files(
  paths$meta,
  pattern = "^manifest_gene_counts_.*\\.tsv$",
  full.names = TRUE
)

if (length(manifest_files) == 0) {
  stop("No manifest found in meta/ (run 010 first)")
}

manifest_file <- tail(sort(manifest_files), 1)
manifest <- data.table::fread(manifest_file)

message("Using manifest: ", basename(manifest_file))
message("Files in manifest: ", nrow(manifest))

# --- Query GDC for case and sample metadata --------------------------------
# Sequential batched queries respect API rate limits. Each batch expands the
# case and sample records attached to a set of file_ids.
file_ids <- manifest$id
n_files <- length(file_ids)
batch_size <- 50L
n_batches <- ceiling(n_files / batch_size)

message("Querying GDC API for ", n_files, " files in ", n_batches, " batches ...")

batch_results <- vector("list", n_batches)

for (b in seq_len(n_batches)) {
  start_idx <- (b - 1L) * batch_size + 1L
  end_idx <- min(b * batch_size, n_files)
  batch_ids <- file_ids[start_idx:end_idx]

  batch_results[[b]] <- tryCatch(
    files() |>
      GenomicDataCommons::filter(~ file_id %in% batch_ids) |>
      GenomicDataCommons::expand(c("cases", "cases.samples")) |>
      results_all(),
    error = function(e) {
      warning("Batch ", b, " failed: ", conditionMessage(e))
      NULL
    }
  )

  if (b %% 5L == 0L) {
    message(b, " / ", n_batches, " batches queried")
  }
  Sys.sleep(0.1)
}

# --- Parse API results in parallel -----------------------------------------
# Each batch result is independent, so parsing distributes across cores. The
# GDC nested structure is flattened into one row per file_id.
parse_batch <- function(batch_result) {
  if (is.null(batch_result)) {
    return(NULL)
  }
  data.table::setDTthreads(1L)

  file_id_to_name <- setNames(batch_result$file_name, batch_result$id)
  rows <- list()

  for (j in seq_along(batch_result$cases)) {
    file_id <- names(batch_result$cases)[j]
    case_data <- batch_result$cases[[j]]
    case_submitter_id <- case_data$submitter_id[1]

    if (is.null(case_data$samples) || length(case_data$samples) == 0) {
      next
    }
    sample_data <- case_data$samples[[1]]

    if (length(sample_data$sample_type) > 0) {
      sample_type_value <- as.character(sample_data$sample_type)[1]
    } else if (length(sample_data$tumor_descriptor) > 0 &&
      length(sample_data$specimen_type) > 0) {
      sample_type_value <- paste(
        sample_data$tumor_descriptor[1],
        sample_data$specimen_type[1],
        "Tumor"
      )
    } else {
      sample_type_value <- "Unknown"
    }

    rows[[file_id]] <- data.frame(
      file_id = file_id,
      file_name = as.character(file_id_to_name[file_id])[1],
      case_submitter_id = as.character(case_submitter_id)[1],
      sample_submitter_id = as.character(sample_data$submitter_id)[1],
      sample_type = sample_type_value,
      sample_id = as.character(sample_data$sample_id)[1],
      stringsAsFactors = FALSE
    )
  }
  rows
}

n_cores <- max(1L, min(detectCores() - 1L, 4L))
message("Parsing results with ", n_cores, " cores ...")

parsed <- mclapply(
  batch_results,
  parse_batch,
  mc.cores = n_cores,
  mc.preschedule = FALSE
)

# Flatten nested per-batch lists into one row set.
mapping_rows <- unlist(parsed, recursive = FALSE)
mapping_rows <- mapping_rows[!vapply(mapping_rows, is.null, logical(1))]
file_sample_mapping <- data.table::rbindlist(mapping_rows, fill = TRUE)

message("Mapped files: ", nrow(file_sample_mapping))

# --- Completeness check ----------------------------------------------------
# Every manifest file_id must resolve to a mapping row.
unmapped <- setdiff(manifest$id, file_sample_mapping$file_id)
if (length(unmapped) > 0) {
  warning("Unmapped file_ids: ", length(unmapped))
}

message(
  "Unique cases: ", length(unique(file_sample_mapping$case_submitter_id)),
  " ; unique samples: ",
  length(unique(file_sample_mapping$sample_submitter_id))
)

# --- Save ------------------------------------------------------------------
mapping_out <- file.path(paths$processed, "file_sample_mapping.rds")
saveRDS(file_sample_mapping, mapping_out)
message("Saved mapping: ", mapping_out)
