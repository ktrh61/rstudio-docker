# analysis_v7/setup.R
# Session setup for v7 analysis

# Set working directory
setwd("/home/rstudio/project")
cat("Working directory:", getwd(), "\n")

# Load libraries (pre-installed)
suppressPackageStartupMessages({
  library(GenomicDataCommons)
  library(SummarizedExperiment)
  library(data.table)
  library(dplyr)
})

# Define paths
paths <- list(
  root = "analysis_v7/",
  scripts = "analysis_v7/scripts/",
  raw = "analysis_v7/data/raw/",
  gdc = "analysis_v7/data/raw/gdc/",
  processed = "analysis_v7/data/processed/",
  output = "analysis_v7/output/"
)

# Display environment
cat("\n=== Environment ===\n")
cat("R version:", R.version$version.string, "\n")
cat("Bioconductor:", as.character(BiocManager::version()), "\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("\nReady for v7 analysis!\n")
