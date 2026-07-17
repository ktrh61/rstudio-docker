# 021_build_se.R
# Load STAR-Counts TSV files into a SummarizedExperiment holding the single
# count assay selected by library strandedness.
# Input : processed/file_sample_mapping.rds       (from 020)
#         raw/expression/<file_id>/<file>.tsv      (from 010)
# Output: processed/thyr_se_raw.rds                (consumed by 030)
#         meta/strand_selection_<timestamp>.tsv    (per-sample strand metrics)
#         meta/loading_metadata_<timestamp>.rds    (run provenance)
#
# Strandedness is decided per sample from the STAR gene-count totals of the two
# stranded columns (Dobin's general rule, ratio form): if the smaller total is
# at most half the larger (ratio <= 0.5) the library is stranded and the larger
# column is used; otherwise it is unstranded and the unstranded column is used.
# N_noFeature and N_ambiguous are recorded for audit but not used to decide.
# Columns are labelled by sample_submitter_id. Only the selected count assay is
# stored; TPM and FPKM are dropped. TSV reading is parallelised; the large
# per-file result list is released before the SE is built to cap the memory peak.

source("setup.R")

library(SummarizedExperiment)
library(parallel)

# Ratio threshold for the stranded/unstranded decision (smaller/larger of the
# two stranded gene-count totals). At or below this the library is stranded.
strand_ratio_threshold <- 0.5

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
# Each worker returns the three strand count vectors (for later selection), the
# two meta rows used for audit, and the stranded gene-count totals used to
# decide strandedness. Lightweight vector arguments keep the fork footprint
# small. Column positions: 4 unstranded, 5 stranded_first, 6 stranded_second.
sample_ids <- tsv_map$sample_submitter_id
paths_vec <- tsv_map$path

read_star_counts <- function(i, paths_vec) {
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
      key <- dt[[1]]
      nf <- as.numeric(unlist(dt[key == "N_noFeature", 4:6]))
      am <- as.numeric(unlist(dt[key == "N_ambiguous", 4:6]))

      # Gene rows only (drop the four leading meta rows).
      genes <- dt[5:nrow(dt), ]
      unstranded <- as.numeric(genes[[4]])
      stranded_first <- as.numeric(genes[[5]])
      stranded_second <- as.numeric(genes[[6]])

      list(
        idx = i,
        unstranded = unstranded,
        stranded_first = stranded_first,
        stranded_second = stranded_second,
        sum_unstranded = sum(unstranded),
        sum_first = sum(stranded_first),
        sum_second = sum(stranded_second),
        nf_unstranded = nf[1], nf_first = nf[2], nf_second = nf[3],
        am_unstranded = am[1], am_first = am[2], am_second = am[3],
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
  mc.cores = n_cores,
  mc.preschedule = FALSE
)

failed <- !vapply(results, function(x) isTRUE(x$success), logical(1))
if (any(failed)) {
  stop("Files that failed to read: ", sum(failed))
}

# --- Decide strandedness per sample ----------------------------------------
# Ratio = smaller / larger of the two stranded gene-count totals. ratio <=
# threshold => stranded (take the larger stranded column); otherwise unstranded
# (take the unstranded column).
message(
  "Deciding strandedness (ratio threshold ", strand_ratio_threshold,
  ") ..."
)

decide_strand <- function(r) {
  larger <- max(r$sum_first, r$sum_second)
  smaller <- min(r$sum_first, r$sum_second)
  ratio <- smaller / larger
  if (ratio <= strand_ratio_threshold) {
    selected <- if (r$sum_first > r$sum_second) {
      "stranded_first"
    } else {
      "stranded_second"
    }
  } else {
    selected <- "unstranded"
  }
  list(ratio = ratio, selected = selected)
}

decisions <- lapply(results, decide_strand)
selected_per_sample <- vapply(decisions, function(d) d$selected, character(1))
selected_levels <- unique(selected_per_sample)

message(
  "Selected column(s): ",
  paste(
    sprintf(
      "%s=%d", names(table(selected_per_sample)),
      as.integer(table(selected_per_sample))
    ),
    collapse = " ; "
  )
)

# A single count assay requires one consistent column across all samples. If
# samples disagree, stop rather than mix columns of different strandedness.
if (length(selected_levels) > 1) {
  stop(
    "Samples disagree on strand selection: ",
    paste(selected_levels, collapse = ", "),
    " (see strand_selection metrics)"
  )
}
selected_column <- selected_levels

# --- Strand selection record (meta) ----------------------------------------
strand_selection <- data.frame(
  sample_submitter_id = sample_ids,
  selected = selected_per_sample,
  ratio = vapply(decisions, function(d) d$ratio, numeric(1)),
  sum_unstranded = vapply(results, function(r) r$sum_unstranded, numeric(1)),
  sum_first = vapply(results, function(r) r$sum_first, numeric(1)),
  sum_second = vapply(results, function(r) r$sum_second, numeric(1)),
  nf_unstranded = vapply(results, function(r) r$nf_unstranded, numeric(1)),
  nf_first = vapply(results, function(r) r$nf_first, numeric(1)),
  nf_second = vapply(results, function(r) r$nf_second, numeric(1)),
  am_unstranded = vapply(results, function(r) r$am_unstranded, numeric(1)),
  am_first = vapply(results, function(r) r$am_first, numeric(1)),
  am_second = vapply(results, function(r) r$am_second, numeric(1)),
  stringsAsFactors = FALSE
)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
strand_file <- file.path(
  paths$meta, paste0("strand_selection_", timestamp, ".tsv")
)
data.table::fwrite(strand_selection, strand_file, sep = "\t")
message("Saved strand selection: ", strand_file)

# --- Assemble the selected count matrix ------------------------------------
message("Assembling count matrix for '", selected_column, "' ...")

count_matrix <- matrix(0, nrow = n_genes, ncol = n_samples)
rownames(count_matrix) <- gene_info$gene_id
colnames(count_matrix) <- sample_ids

for (r in results) {
  count_matrix[, r$idx] <- r[[selected_column]]
}

# Release the large per-file result list before building the SE to cap the
# memory peak (results holds every file's three strand count vectors).
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
  assays = setNames(list(count_matrix), selected_column),
  colData = col_data,
  rowData = row_data
)

message(
  "SummarizedExperiment: ", nrow(thyr_se), " genes x ",
  ncol(thyr_se), " samples ; assay: ", assayNames(thyr_se)
)

# --- Save ------------------------------------------------------------------
se_out <- file.path(paths$processed, "thyr_se_raw.rds")
saveRDS(thyr_se, se_out)
message("Saved SE: ", se_out)

metadata_out <- file.path(
  paths$meta, paste0("loading_metadata_", timestamp, ".rds")
)
saveRDS(
  list(
    loading_date = Sys.Date(),
    n_files_loaded = n_samples,
    n_genes = n_genes,
    selected_column = selected_column,
    strand_ratio_threshold = strand_ratio_threshold,
    parallel_cores = n_cores
  ),
  metadata_out
)
message("Saved metadata: ", metadata_out)
