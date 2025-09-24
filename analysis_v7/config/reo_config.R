# config/reo_config.R - REO Panel Construction Configuration
# Purpose: Centralized configuration for REO analysis parameters
# Note: This file is read ONLY by Step 1. Steps 2+ use snapshot from RDS.
# Version: v1.0
# Date: 2025-01-26

# ============================================================================
# REO-specific Configuration (from REO_Panel_Strategy.md)
# ============================================================================

CONFIG <- list(
  # Basic parameters
  tau_strength = log2(1.5),        # Strength criterion
  dead_zone = log2(1.2),           # Dead zone threshold
  trim_low = 0.07,                 # Lower trim rate (ADJUSTED from 0.10)
  trim_high = 0.10,                # Upper trim rate
  M = 12,                          # Maximum pairs per gene (ADJUSTED from 10)
  N = 10,                          # Panel size (default)
  T_ratio = 0.4,                   # Decision threshold ratio
  CI_method = "Wilson",            # CI calculation method (no continuity correction)
  
  # Stopping condition parameters
  min_candidates_abs = 100,        # Minimum candidate count (absolute)
  min_candidates_ratio = 0.02,     # Minimum candidate ratio (relative to initial)
  
  # Processing parameters
  VERBOSE = TRUE                   # Verbose output
)

# Log configuration at load time (for transparency)
if (CONFIG$VERBOSE) {
  cat("\n=== REO Configuration Loaded ===\n")
  cat("  tau_strength:", CONFIG$tau_strength, "\n")
  cat("  dead_zone:", CONFIG$dead_zone, "\n")
  cat("  trim_low:", CONFIG$trim_low, "\n")
  cat("  trim_high:", CONFIG$trim_high, "\n")
  cat("  M:", CONFIG$M, "\n")
  cat("  NOTE: trim_low and M have been adjusted from defaults\n")
  cat("================================\n\n")
}