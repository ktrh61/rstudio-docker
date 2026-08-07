# Resolve RET-cohort cases by radiation exposure band, with paired _merged
# tumour / normal sample IDs. Shared by the provisional REO QC scripts
# (diagnostics/reo_lowmid_outliers.R, reo_lowmid_purity.R) so their R_Low / R_Mid definitions match the committed pipeline:
# Sporadic dose 0; among exposed (dose > 0) High assigned share >= 66.6, Mid
# (33.3, 66.6), Low (0, 33.3]. Only RET (R cohort) is returned.
resolve_ret_exposure_cases <- function(clinical, assigned_share, se) {
  ret_values <- c("CCDC6-RET", "NCOA4-RET", "RET-OTHER")
  ret_cases <- as.character(clinical$REBC_ID)[
    clinical$Designated_Driver %in% ret_values
  ]

  a <- data.frame(
    case_submitter_id = as.character(assigned_share$REBC_ID),
    dose_mgy = as.numeric(assigned_share$dose_mgy),
    assigned_share = as.numeric(assigned_share$assigned_share_approx),
    stringsAsFactors = FALSE
  )
  a <- a[a$case_submitter_id %in% ret_cases & is.finite(a$dose_mgy), , drop = FALSE]

  band <- rep(NA_character_, nrow(a))
  band[a$dose_mgy == 0] <- "R_Sporadic"
  exposed <- a$dose_mgy > 0 & is.finite(a$assigned_share)
  band[exposed & a$assigned_share >= 66.6] <- "R_High"
  band[exposed & a$assigned_share > 33.3 & a$assigned_share < 66.6] <- "R_Mid"
  band[exposed & a$assigned_share > 0 & a$assigned_share <= 33.3] <- "R_Low"
  a$band <- band
  a <- a[!is.na(a$band), , drop = FALSE]

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  m <- cd[grepl("_merged", cd$sample_submitter_id), , drop = FALSE]
  tumor_of <- setNames(
    m$sample_submitter_id[m$sample_type == "Primary Tumor"],
    m$case_submitter_id[m$sample_type == "Primary Tumor"]
  )
  normal_of <- setNames(
    m$sample_submitter_id[m$sample_type == "Solid Tissue Normal"],
    m$case_submitter_id[m$sample_type == "Solid Tissue Normal"]
  )
  a$tumor_id <- unname(tumor_of[a$case_submitter_id])
  a$normal_id <- unname(normal_of[a$case_submitter_id])
  a
}
