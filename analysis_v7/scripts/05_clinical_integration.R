# 05_clinical_integration.R - Integrate clinical data with RNA-seq
# Purpose: Extract BRAF/RET samples and link with clinical data

source("analysis_v7/setup.R")

cat("\n=== Clinical Data Integration (v7) ===\n")
cat("Date:", as.character(Sys.Date()), "\n")

# Load data
se <- readRDS(paste0(paths$processed, "thyr_se_annotated.rds"))
clinical <- data.table::fread(paste0(paths$raw, "abg2538-data-s1.txt"))

cat("\nInput data:\n")
cat("  Total RNA-seq samples:", ncol(se), "\n")
cat("  Total clinical cases:", nrow(clinical), "\n")

# Identify BRAF/RET cases
braf_cases <- clinical[
  (clinical$WGS_CandidateDriverMutation == "" | clinical$WGS_CandidateDriverMutation == "BRAF") &
  (clinical$RNA_CandidateDriverMutation == "" | clinical$RNA_CandidateDriverMutation == "BRAF") &
  clinical$Designated_Driver == "BRAF.MutV600E", 
]

ret_cases <- clinical[clinical$DriverGroup_CoMutFig == "RET-Fusion", ]

braf_ids <- braf_cases$REBC_ID
ret_ids <- ret_cases$REBC_ID
target_ids <- unique(c(braf_ids, ret_ids))

cat("\nTarget cases:\n")
cat("  BRAF V600E:", length(braf_ids), "\n")
cat("  RET fusion:", length(ret_ids), "\n")
cat("  Total unique:", length(target_ids), "\n")

# RET fusion subtypes
cat("\nRET fusion subtypes:\n")
print(table(ret_cases$Designated_Driver))

# Filter SE to target samples
coldata_df <- as.data.frame(colData(se))
target_mask <- coldata_df$case_submitter_id %in% target_ids
se_filtered <- se[, target_mask]

# Create metadata with clear case/sample distinction
metadata <- data.frame(
  sample_id = colData(se_filtered)$sample_submitter_id,
  case_id = colData(se_filtered)$case_submitter_id,
  sample_type = colData(se_filtered)$sample_type,
  is_tumor = colData(se_filtered)$is_tumor,
  is_normal = colData(se_filtered)$is_normal,
  file_id = colData(se_filtered)$file_id,
  stringsAsFactors = FALSE
)

# Add case-level information
clinical_info <- clinical[clinical$REBC_ID %in% target_ids, .(
  REBC_ID, 
  DOSE, 
  AGE_SURGERY, 
  AGE_EXPOSURE,
  SEX,
  Designated_Driver
)]

# Calculate Birth Year and Diagnosis Year
# Chernobyl accident: 1986
clinical_info$birth_year <- 1986 - clinical_info$AGE_EXPOSURE
clinical_info$diagnosis_year <- clinical_info$birth_year + clinical_info$AGE_SURGERY

metadata <- merge(metadata, clinical_info, 
                  by.x = "case_id", by.y = "REBC_ID", all.x = TRUE)

# Add driver status (case-level and sample-level)
# Case-level: the driver mutation of the case
metadata$case_driver <- NA
metadata$case_driver[metadata$case_id %in% braf_ids] <- "BRAF"
metadata$case_driver[metadata$case_id %in% ret_ids] <- "RET"

# 列の整理
metadata <- metadata[, c(
  "sample_id", "case_id", "sample_type", "is_tumor", "is_normal",
  "case_driver", "Designated_Driver",  # case_driverで全て管理
  "SEX", "DOSE", "AGE_SURGERY", "AGE_EXPOSURE", 
  "birth_year", "diagnosis_year", "file_id"
)]

# Summary
cat("\nSamples by case driver and type:\n")
print(table(CaseDriver = metadata$case_driver, SampleType = metadata$sample_type))

cat("\nDesignated drivers (tumors only):\n")
tumor_drivers <- metadata[metadata$is_tumor, c("case_driver", "Designated_Driver")]
print(table(tumor_drivers$Designated_Driver))

cat("\nDose availability by driver:\n")
unique_cases <- metadata[!duplicated(metadata$case_id), ]
print(table(Driver = unique_cases$case_driver, HasDose = !is.na(unique_cases$DOSE)))

# Update SE
colData(se_filtered) <- DataFrame(metadata)
rownames(colData(se_filtered)) <- metadata$sample_id

# Save
saveRDS(se_filtered, paste0(paths$processed, "thyr_se_braf_ret.rds"))
save(se_filtered, file = paste0(paths$processed, "thyr_se_braf_ret.RData"))

data.table::fwrite(metadata, 
                   paste0(paths$processed, "braf_ret_metadata.tsv"), 
                   sep = "\t")

cat("\n=== Clinical Integration Completed ===\n")
cat("Retained ", ncol(se_filtered), " samples from ", length(target_ids), " BRAF/RET cases\n")
cat("Files created:\n")
cat("  - thyr_se_braf_ret.rds\n")
cat("  - braf_ret_metadata.tsv\n")
cat("\nNext: Calculate ERR for exposure classification\n")
