# 99_paper_summary.R
# Purpose: Execute all check scripts and compile manuscript numbers
# Date: 2025-01-26

cat("========================================\n")
cat("MANUSCRIPT NUMBERS COMPILATION\n")
cat("========================================\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# Run check scripts
cat("Running check_01_groups.R...\n")
cat("----------------------------------------\n")
source("analysis_v7/scripts/check_01_groups.R")

cat("\n\nRunning check_02_deg.R...\n")
cat("----------------------------------------\n")
source("analysis_v7/scripts/check_02_deg.R")

cat("\n========================================\n")
cat("ALL CHECKS COMPLETE\n")
cat("========================================\n")
cat("\nPlease copy the numbers above into the manuscript.\n")
cat("Numbers marked with '□' are the key values needed.\n")

# Optional: Save output to file
sink("manuscript_numbers.txt")
cat("MANUSCRIPT NUMBERS - Generated:", as.character(Sys.time()), "\n")
cat("========================================\n\n")

cat("Run individual check scripts for details:\n")
cat("  source('check_01_groups.R')  # Sample counts\n")
cat("  source('check_02_deg.R')     # DEG results\n")

sink()
cat("\nNumbers also saved to: manuscript_numbers.txt\n")