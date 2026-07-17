# 021_build_se.R
# Load STAR-Counts TSV files into a SummarizedExperiment with all six assays.
# Input : processed/file_sample_mapping.rds       (from 020)
#         raw/expression/<file_id>/<file>.tsv      (from 010)
# Output: processed/thyr_se_raw.rds                (consumed by 030)
#         meta/loading_metadata_<timestamp>.rds    (run provenance)
#
# Columns are labelled by sample_submitter_id. All six assays (three strand
# counts plus TPM / FPKM / FPKM-UQ) are stored; assay selection is left to 030.
# TSV reading is parallelised; the large per-file result list is released after
# the matrices are populated to cap the memory peak before the SE is built.

source("setup.R")

library(SummarizedExperiment)
library(parallel)

# --- Load mapping ----------------------------------------------------------
mapping_file <- file.path(paths$processed, "file_sample_mapping.rds")
if (!file.exists(mapping_file)) {
  stop("file_sample_mapping.rds not found in processed/ (run 020 first)")
}
file_sample_mapping <- readRDS(mapping_file)

# --- Locate TSV files and attach mapping -----------------------------------
# Layout: raw/expression/<file_id>/<file_name>. file_id is the parent
# directory name.
expr_dir <- file.path(paths$raw, "expression")
tsv_files <- list.files(
  expr_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(tsv_files) == 0) {
  stop("No TSV files found under raw/expression/ (run 010 first)")
}
message("TSV files found: ", length(tsv_files))

file_ids_from_path <- basename(dirname(tsv_files))
tsv_map <- data.frame(
  path = tsv_files,
  file_id = file_ids_from_path,
  stringsAsFactors = FALSE
)

tsv_map <- merge(
  tsv_map,
  file_sample_mapping[, c(
    "file_id", "sample_submitter_id", "case_submitter_id", "sample_type"
  )],
  by = "file_id",
  all.x = TRUE
)

missing_mapping <- is.na(tsv_map$sample_submitter_id)
if (any(missing_mapping)) {
  warning("TSV files without mapping: ", sum(missing_mapping))
  tsv_map <- tsv_map[!missing_mapping, ]
}
message("TSV files with valid mapping: ", nrow(tsv_map))

# Stable column ordering by sample id.
tsv_map <- tsv_map[order(tsv_map$sample_submitter_id), ]

# --- Gene structure from the first file ------------------------------------
# The STAR TSV has a comment line, a header, then four summary rows
# (N_unmapped, N_multimapping, N_noFeature, N_ambiguous) before gene rows.
first_data <- data.table::fread(
  tsv_map$path[1],
  skip = 1L,
  showProgress = FALSE,
  nThread = 1L
)
gene_info <- first_data[5:nrow(first_data), 1:3]
colnames(gene_info) <- c("gene_id", "gene_name", "gene_type")
n_genes <- nrow(gene_info)
n_samples <- nrow(tsv_map)
rm(first_data)

message("Genes: ", n_genes, " ; samples: ", n_samples)

# --- Read all count files in parallel --------------------------------------
# Lightweight vector arguments keep the fork footprint small. Columns are
# extracted by position (four count/expression triplets in columns 4-9).
sample_ids <- tsv_map$sample_submitter_id
paths_vec <- tsv_map$path

read_star_counts <- function(i, paths_vec, ids_vec) {
  data.table::setDTthreads(1L)
  path <- paths_vec[i]

  tryCatch(
    {
      dt <- data.table::fread(
        path,
        skip = 1L,
        showProgress = FALSE,
        nThread = 1L
      )
      dt <- dt[5:nrow(dt), ]
      list(
        idx = i,
        unstranded = as.numeric(dt[[4]]),
        stranded_first = as.numeric(dt[[5]]),
        stranded_second = as.numeric(dt[[6]]),
        tpm_unstranded = as.numeric(dt[[7]]),
        fpkm_unstranded = as.numeric(dt[[8]]),
        fpkm_uq_unstranded = as.numeric(dt[[9]]),
        success = TRUE
      )
    },
    error = function(e) {
      warning("Failed to read ", basename(path), ": ", conditionMessage(e))
      list(idx = i, success = FALSE)
    }
  )
}

n_cores <- max(1L, min(detectCores() - 1L, 8L))
message("Reading count files with ", n_cores, " cores ...")

results <- mclapply(
  seq_len(n_samples),
  read_star_counts,
  paths_vec = paths_vec,
  ids_vec = sample_ids,
  mc.cores = n_cores,
  mc.preschedule = FALSE
)

failed <- !vapply(results, function(x) isTRUE(x$success), logical(1))
if (any(failed)) {
  warning("Files that failed to read: ", sum(failed))
}

# --- Assemble assay matrices -----------------------------------------------
message("Assembling assay matrices ...")

assay_names <- c(
  "unstranded", "stranded_first", "stranded_second",
  "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded"
)

make_matrix <- function() {
  m <- matrix(0, nrow = n_genes, ncol = n_samples)
  rownames(m) <- gene_info$gene_id
  colnames(m) <- sample_ids
  m
}
assay_list <- setNames(
  lapply(assay_names, function(x) make_matrix()),
  assay_names
)

for (r in results) {
  if (!isTRUE(r$success)) {
    next
  }
  col_idx <- r$idx
  for (a in assay_names) {
    assay_list[[a]][, col_idx] <- r[[a]]
  }
}

# Release the large per-file result list before building the SE to cap the
# memory peak (results holds every file's six numeric vectors).
rm(results)
gc()

# --- Build SummarizedExperiment --------------------------------------------
col_data <- data.frame(
  sample_submitter_id = tsv_map$sample_submitter_id,
  case_submitter_id = tsv_map$case_submitter_id,
  sample_type = tsv_map$sample_type,
  file_id = tsv_map$file_id,
  stringsAsFactors = FALSE
)
rownames(col_data) <- col_data$sample_submitter_id

row_data <- data.frame(
  gene_id = gene_info$gene_id,
  gene_name = gene_info$gene_name,
  gene_type = gene_info$gene_type,
  stringsAsFactors = FALSE
)

thyr_se <- SummarizedExperiment(
  assays = assay_list,
  colData = col_data,
  rowData = row_data
)

message(
  "SummarizedExperiment: ", nrow(thyr_se), " genes x ",
  ncol(thyr_se), " samples ; assays: ",
  paste(assayNames(thyr_se), collapse = ", ")
)

# --- Save ------------------------------------------------------------------
se_out <- file.path(paths$processed, "thyr_se_raw.rds")
saveRDS(thyr_se, se_out)
message("Saved SE: ", se_out)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
metadata_out <- file.path(
  paths$meta, paste0("loading_metadata_", timestamp, ".rds")
)
saveRDS(
  list(
    loading_date = Sys.Date(),
    n_files_loaded = n_samples,
    n_genes = n_genes,
    parallel_cores = n_cores
  ),
  metadata_out
)
message("Saved metadata: ", metadata_out)
