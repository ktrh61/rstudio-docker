# ==============================================================================
# REBC-THYR Radiation-Induced Classification: Complete Pipeline
# final_classification_pipeline.R
# ==============================================================================

# Required libraries
library(SummarizedExperiment)
library(readxl)
library(dplyr)

cat("=== REBC-THYR Radiation-Induced vs Sporadic Classification ===\n")
cat("Target: Identify radiation-induced cases in RET-fusion thyroid cancer\n\n")

# ==============================================================================
# 1. Load Data
# ==============================================================================

cat("Step 1: Loading data...\n")

# Load TPM data
load("./data/processed/thyr_tpm.rda")

# Load high-purity sample lists
load("./data/processed/final_high_purity_sample_lists_v6.rda")

# Load DEG results
load("./data/processed/final_deg_v6_results.rda")

# Get R0 and R1 samples
r0_samples <- final_high_purity_sample_lists_v6$R0$tumor
r1_samples <- final_high_purity_sample_lists_v6$R1$tumor

cat(sprintf("  R0 samples (ERR=0): %d\n", length(r0_samples)))
cat(sprintf("  R1 samples (ERR≥66.6%%): %d\n", length(r1_samples)))

# ==============================================================================
# 2. Define Gene Pairs
# ==============================================================================

cat("\nStep 2: Defining classification model...\n")

# Gene pair definitions based on comprehensive analysis
pair1 <- list(
  up_gene = "ENSG00000198028.4",    # ZNF560
  down_gene = "ENSG00000145423.5",  # SFRP2
  up_symbol = "ZNF560",
  down_symbol = "SFRP2"
)

pair2 <- list(
  up_gene = "ENSG00000167434.10",   # CA4
  down_gene = "ENSG00000143546.10", # S100A8
  up_symbol = "CA4",
  down_symbol = "S100A8"
)

cat("  Pair 1: ZNF560/SFRP2 (Transcription/WNT signaling)\n")
cat("  Pair 2: CA4/S100A8 (Metabolism/Inflammation)\n")

# ==============================================================================
# 3. Calculate Ratios and Thresholds
# ==============================================================================

cat("\nStep 3: Calculating ratios and thresholds...\n")

pc <- 0.1  # pseudocount

# Pair 1 ratios
ratio1_r0 <- log2((tpm[pair1$up_gene, r0_samples] + pc) / 
                    (tpm[pair1$down_gene, r0_samples] + pc))
ratio1_r1 <- log2((tpm[pair1$up_gene, r1_samples] + pc) / 
                    (tpm[pair1$down_gene, r1_samples] + pc))
threshold1 <- (median(ratio1_r0) + median(ratio1_r1)) / 2

# Pair 2 ratios
ratio2_r0 <- log2((tpm[pair2$up_gene, r0_samples] + pc) / 
                    (tpm[pair2$down_gene, r0_samples] + pc))
ratio2_r1 <- log2((tpm[pair2$up_gene, r1_samples] + pc) / 
                    (tpm[pair2$down_gene, r1_samples] + pc))
threshold2 <- (median(ratio2_r0) + median(ratio2_r1)) / 2

cat(sprintf("  Pair 1 threshold: %.2f\n", threshold1))
cat(sprintf("  Pair 2 threshold: %.2f\n", threshold2))

# ==============================================================================
# 4. Classification Function
# ==============================================================================

classify_sample <- function(sample_id, tpm_data, pair1, pair2, 
                            threshold1, threshold2, pc = 0.1) {
  # Calculate ratios
  ratio1 <- log2((tpm_data[pair1$up_gene, sample_id] + pc) / 
                   (tpm_data[pair1$down_gene, sample_id] + pc))
  ratio2 <- log2((tpm_data[pair2$up_gene, sample_id] + pc) / 
                   (tpm_data[pair2$down_gene, sample_id] + pc))
  
  # Apply conservative rule: both must exceed threshold
  if (ratio1 > threshold1 && ratio2 > threshold2) {
    return("radiation-induced")
  } else {
    return("sporadic")
  }
}

# ==============================================================================
# 5. Classify R1 Samples
# ==============================================================================

cat("\nStep 4: Classifying R1 samples...\n")

r1_classifications <- sapply(r1_samples, function(s) {
  classify_sample(s, tpm, pair1, pair2, threshold1, threshold2, pc)
})

n_rad <- sum(r1_classifications == "radiation-induced")
cat(sprintf("  R1 results: %d/%d (%.1f%%) radiation-induced\n", 
            n_rad, length(r1_samples), n_rad/length(r1_samples)*100))

# ==============================================================================
# 6. Test on Intermediate ERR (Optional)
# ==============================================================================

cat("\nStep 5: Testing on intermediate ERR samples...\n")

# Load ERR data
ret_err <- read_xlsx("./data/raw/thyr_poc_ret_fusion.xlsx")
colnames(ret_err)[1] <- "REBC_ID"

# Get intermediate ERR cases
intermediate_cases <- ret_err[!is.na(ret_err$ERR) & 
                                ret_err$ERR > 0 & 
                                ret_err$ERR < 66.6, ]

if (nrow(intermediate_cases) > 0) {
  # Load sample metadata
  sample_data <- read.delim("./data/raw/sample.tsv", stringsAsFactors = FALSE)
  
  # Find matching samples
  intermediate_samples <- c()
  intermediate_err_values <- c()
  
  for (i in 1:nrow(intermediate_cases)) {
    case_id <- intermediate_cases$REBC_ID[i]
    err_val <- intermediate_cases$ERR[i]
    
    # Find tumor samples for this case
    case_samples <- sample_data[sample_data$case_submitter_id == case_id & 
                                  sample_data$sample_type == "Primary Tumor", ]
    
    if (nrow(case_samples) > 0) {
      # Check with and without _merged suffix
      for (sample_id in case_samples$sample_submitter_id) {
        if (sample_id %in% colnames(tpm)) {
          intermediate_samples <- c(intermediate_samples, sample_id)
          intermediate_err_values <- c(intermediate_err_values, err_val)
          break
        } else if (paste0(sample_id, "_merged") %in% colnames(tpm)) {
          intermediate_samples <- c(intermediate_samples, paste0(sample_id, "_merged"))
          intermediate_err_values <- c(intermediate_err_values, err_val)
          break
        }
      }
    }
  }
  
  if (length(intermediate_samples) > 0) {
    # Classify intermediate samples
    intermediate_classifications <- sapply(intermediate_samples, function(s) {
      classify_sample(s, tpm, pair1, pair2, threshold1, threshold2, pc)
    })
    
    # Summarize by ERR group
    err_groups <- cut(intermediate_err_values,
                      breaks = c(0, 33.3, 66.6),
                      labels = c("0-33.3%", "33.3-66.6%"),
                      include.lowest = TRUE)
    
    cat(sprintf("  Found %d intermediate ERR samples\n", length(intermediate_samples)))
    
    for (group in levels(err_groups)) {
      group_samples <- err_groups == group
      if (sum(group_samples) > 0) {
        n_rad_group <- sum(intermediate_classifications[group_samples] == "radiation-induced")
        n_total_group <- sum(group_samples)
        cat(sprintf("  ERR %s: %d/%d (%.1f%%) radiation-induced\n",
                    group, n_rad_group, n_total_group, 
                    n_rad_group/n_total_group*100))
      }
    }
  }
}

# ==============================================================================
# 7. Save Model
# ==============================================================================

cat("\nStep 6: Saving classification model...\n")

classification_model <- list(
  pair1 = c(pair1, threshold = threshold1),
  pair2 = c(pair2, threshold = threshold2),
  classification_rule = "Both ratios must exceed thresholds",
  r0_samples = r0_samples,
  r1_samples = r1_samples,
  r1_classifications = r1_classifications,
  model_date = Sys.time()
)

save(classification_model, file = "./data/processed/final_radiation_classification_model.rda")

cat("\n=== Classification Model Complete ===\n")
cat("Model saved to: ./data/processed/final_radiation_classification_model.rda\n")
cat("\nUsage: For new sample, calculate both ratios and apply conservative rule\n")
cat("Expected performance: ~67% sensitivity in high ERR group\n")