# sex_chromosome_annotation.R
# Sex-chromosome membership of the discovered gene lists (descriptive
# annotation; researcher Go 2026-08-16).
#
# Purpose: Methods (Candidate confounders) discloses sex as a candidate
# confounder and states that discovered genes are annotated for
# sex-chromosome membership. This diagnostic computes that annotation:
# for each analysis contrast, the number of chrX / chrY genes among the
# tested genes and among the q < FDR_CUT discoveries, for the discoveries
# also split by direction (effect > 0.5 = higher in exposed). Counts
# only -- no enrichment statistic, no test, no p-value; the reading was
# ruled not to need pre-fixing (descriptive diagnostic).
#
# Chromosome source: gene rows of the GDC GENCODE v36 GTF -- the
# reference already used for the gene lengths (020) and the same source
# that verified BHLHB9/chrX by hand (numbers ledger N-80). That lookup
# is re-derived here and asserted as a cross-check.
#
# Input : raw/reference/gencode.v36.annotation.gtf
#         processed/thyr_expression_test.rds  (410: gene_id, effect, q_storey)
#         processed/thyr_se_raw.rds           (gene_id -> gene_name)
# Output: diagnostics/output/sex_chromosome_annotation.rds

source("setup.R")
suppressPackageStartupMessages({
  library(rtracklayer)
})

source(file.path(paths$root, "lib", "units.R"))
source(file.path(paths$root, "lib", "annotation.R"))

gtf_file <- file.path(paths$raw, "reference", "gencode.v36.annotation.gtf")
test <- readRDS(file.path(paths$processed, "thyr_expression_test.rds"))
se <- readRDS(file.path(paths$processed, "thyr_se_raw.rds"))
name_of <- gene_name_map(se)

message("importing GTF gene rows (chromosome per gene_id) ...")
genes_gr <- rtracklayer::import(gtf_file, feature.type = "gene")
chrom_of <- setNames(as.character(GenomicRanges::seqnames(genes_gr)),
                     genes_gr$gene_id)

# --- Per-contrast counts ----------------------------------------------------
summarize_contrast <- function(u) {
  g <- test$units[[u]]$genes
  chrom <- chrom_of[g$gene_id]
  stopifnot(!anyNA(chrom)) # every tested gene must map to the reference
  disc <- g$q_storey < FDR_CUT
  up <- disc & g$effect > 0.5
  down <- disc & g$effect < 0.5
  stopifnot(sum(up) + sum(down) == sum(disc)) # no discovery sits at effect 0.5
  data.frame(
    contrast = u,
    tested = nrow(g),
    tested_chrX = sum(chrom == "chrX"), tested_chrY = sum(chrom == "chrY"),
    discovered = sum(disc),
    disc_chrX = sum(disc & chrom == "chrX"),
    disc_chrY = sum(disc & chrom == "chrY"),
    up_chrX = sum(up & chrom == "chrX"), up_chrY = sum(up & chrom == "chrY"),
    down_chrX = sum(down & chrom == "chrX"),
    down_chrY = sum(down & chrom == "chrY"),
    stringsAsFactors = FALSE
  )
}
tab <- do.call(rbind, lapply(UNIT_ORDER, summarize_contrast))
print(tab, row.names = FALSE)

# --- Cross-check against N-80: the single B_Normal discovery is chrX --------
bn <- test$units[["B_Normal"]]$genes
bn_disc <- bn[bn$q_storey < FDR_CUT, ]
stopifnot(nrow(bn_disc) == 1L,
          chrom_of[bn_disc$gene_id] == "chrX",
          gene_label(bn_disc$gene_id, name_of) == "BHLHB9")
message("N-80 cross-check passed (B_Normal single discovery = BHLHB9, chrX)")

# --- The R_Tumor discovered chrX/chrY genes, listed --------------------------
rt <- test$units[["R_Tumor"]]$genes
rt$chrom <- chrom_of[rt$gene_id]
xy <- rt[rt$q_storey < FDR_CUT & rt$chrom %in% c("chrX", "chrY"), ]
xy$gene_name <- gene_label(xy$gene_id, name_of)
xy$direction <- ifelse(xy$effect > 0.5, "higher in exposed", "lower in exposed")
xy <- xy[order(xy$chrom, xy$p_exact),
         c("gene_id", "gene_name", "chrom", "direction",
           "effect", "p_exact", "q_storey")]
print(xy, row.names = FALSE)

out <- file.path(paths$root, "diagnostics", "output",
                 "sex_chromosome_annotation.rds")
saveRDS(list(summary = tab, r_tumor_xy = xy), out)
message("Saved: ", out)
