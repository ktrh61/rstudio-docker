# 020_prepare_gene_lengths.R
# Prepare gene lengths (exon union length per gene) from the GDC GENCODE v36 GTF.
# Input : GDC reference GTF (downloaded; md5-verified)
# Output: processed/gene_lengths.rds  (named numeric vector: ENSG id -> length in bp)

source("setup.R")

library(rtracklayer)
library(GenomicRanges)

# --- Reference source (GDC reference files) --------------------------------
# GENCODE v36 GTF as distributed by the GDC. Version rationale is in README.
gtf_url <- "https://api.gdc.cancer.gov/data/be002a2c-3b27-43f3-9e0f-fd47db92a6b5"
gtf_md5 <- "c03931958d4572148650d62eb6dec41a"
gtf_gz <- file.path(paths$raw, "reference", "gencode.v36.annotation.gtf.gz")
gtf_file <- file.path(paths$raw, "reference", "gencode.v36.annotation.gtf")

out_rds <- file.path(paths$processed, "gene_lengths.rds")

dir.create(dirname(gtf_gz), recursive = TRUE, showWarnings = FALSE)

# --- Download --------------------------------------------------------------
# Skip download only if a file with the expected md5 is already present.
have_valid <- file.exists(gtf_gz) && tools::md5sum(gtf_gz)[[1]] == gtf_md5

if (!have_valid) {
  options(timeout = 600)
  download.file(gtf_url, gtf_gz, method = "libcurl", quiet = FALSE)

  md5_observed <- tools::md5sum(gtf_gz)[[1]]
  if (md5_observed != gtf_md5) {
    file.remove(gtf_gz)
    stop(
      "GTF md5 mismatch: expected ", gtf_md5, ", got ", md5_observed,
      " (removed incomplete file)"
    )
  }
}

# --- Decompress ------------------------------------------------------------
if (!file.exists(gtf_file)) {
  R.utils::gunzip(gtf_gz, gtf_file, remove = FALSE, overwrite = TRUE)
}

# --- Compute exon union length per gene ------------------------------------
message("Importing GTF ...")
gtf <- import(gtf_file)
exons <- gtf[gtf$type == "exon"]
exons$gene_id_clean <- sub("\\..*", "", exons$gene_id)

message("Computing gene lengths ...")
exons_by_gene <- split(exons, exons$gene_id_clean)
gene_lengths <- as.numeric(sum(width(reduce(exons_by_gene))))
names(gene_lengths) <- names(exons_by_gene)

# --- Drop zero-length genes ------------------------------------------------
n_zero <- sum(gene_lengths == 0)
if (n_zero > 0) {
  gene_lengths <- gene_lengths[gene_lengths > 0]
}

# --- Save ------------------------------------------------------------------
saveRDS(gene_lengths, out_rds)

message("Genes with length: ", length(gene_lengths))
message("Zero-length genes dropped: ", n_zero)
message("Saved: ", out_rds)
