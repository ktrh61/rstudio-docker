# reo_config.R - REO Panel Construction Configuration
# Purpose: Centralized configuration for REO analysis parameters
# Note: This file is read ONLY by Step 1. Steps 2+ use snapshot from RDS.
# Version: v7.0
# Date: 2025-01-27

# ============================================================================
# Configuration Version Information
# ============================================================================

CONFIG <- list(
  # Version tracking
  config_version = "v7.0",
  config_date = "2025-01-27",
  
  # ============================================================================
  # Basic Parameters (Referenced across all steps)
  # ============================================================================
  
  # Core thresholds
  tau_strength = log2(1.5),        # Strength criterion (Step 3a)
  dead_zone = log2(1.2),           # Dead zone threshold (Steps 2/3/5/6) - FIXED, DO NOT ADJUST
  n_eff_threshold = 0.7,           # Minimum effective rate for sample judgment (Steps 5/6)
  
  # Expression trimming
  trim_low = 0.07,                 # Lower trim rate (ADJUSTED from 0.10 due to candidate shortage)
  trim_high = 0.10,                # Upper trim rate
  
  # Pairing constraints
  M = 12,                          # Maximum pairs per gene (ADJUSTED from 10)
  
  # CI calculation
  ci_method = "wilson",            # CI calculation method (no continuity correction)
  ci_level = 0.95,                 # Confidence level for intervals
  
  # Panel size (display default; Step5 overrides via PANEL_CONFIG)
  N = 10,                          # (display default; Step5 overrides via PANEL_CONFIG)
  T_ratio = 0.4,                   # Decision threshold ratio
  
  # Stopping conditions
  min_candidates_abs = 100,        # Minimum candidate count (absolute)
  min_candidates_ratio = 0.02,     # Minimum candidate ratio (relative to initial)
  
  # POC handling - CRITICAL NOTE
  # POC is used ONLY for stratification in Step 6 (R_test)
  # POC is NOT used in REO judgment logic (Steps 2-5)
  poc_cutpoints = c(0, 1/3, 2/3),  # For display and stratification (Step 6)
  poc_note = "POC is for stratification ONLY, not for REO judgment",
  
  # Processing parameters
  VERBOSE = TRUE,                  # Verbose output
  
  # ============================================================================
  # Step 3 Specific Parameters
  # ============================================================================
  
  # Order statistics rule
  a2_order_rule = "nR0_minus_1",   # Use (n_R0-1)th order statistic
  
  # Jonckheere-Terpstra test
  emit_jt_diagnostic = TRUE,       # Record JT Z-values for POC tertiles
  jt_alternative = "increasing",   # Alternative hypothesis for JT test
  
  # ============================================================================
  # Step 4 Specific Parameters (Deduplication)
  # ============================================================================
  
  # Correlation clustering
  cor_threshold = 0.9,             # Spearman correlation threshold for clustering
  cor_method = "spearman",         # Correlation method
  
  # Matching priority (lexicographic ordering)
  matching_priority = c("a2", "r1_reversal", "missing", "alpha"),
  
  # ============================================================================
  # Step 5 Specific Parameters (Panel)
  # ============================================================================
  
  # Panel configuration
  panel_defaults = list(
    size_target = 10,             # Target panel size
    size_min = 8,                 # Minimum acceptable size
    size_max = 12,                # Maximum panel size
    t_ratio = 0.4,                # Threshold ratio for T calculation
    max_missing_rate = 0.20       # Maximum missing rate for panel pairs
  ),
  
  # Classification boundaries
  boundary_band = c(0.3, 0.5),    # REO-score range for "Undetermined"
  
  # QC anchor configuration
  qc_anchor = list(
    expr_pct_low = 0.2,           # Lower expression percentile (20th)
    expr_pct_high = 0.8,          # Upper expression percentile (80th)
    delta_max = log2(1.2),        # Maximum group difference (same as dead_zone)
    n = 300,                      # Number of anchor genes
    method = "MAD"                # Selection method (by MAD in log2TPM scale)
  ),
  
  # ============================================================================
  # Step 6 Specific Parameters (R_test)
  # ============================================================================
  
  # R_test configuration
  rtest = list(
    # Main analysis
    purity_threshold = 0.6,       # Tumor purity threshold (contamDE), not POC
    main_layers = 2,              # Primary analysis: 2 layers (≤33.3%, >33.3%)
    
    # Statistical thresholds
    min_n_stratum = 5,            # Minimum samples per stratum
    ci_method = "wilson",         # CI method (redundant but explicit)
    ci_level = 0.95,              # Confidence level
    
    # Diagnostic outputs
    emit_symmetry = TRUE,         # Output symmetry analysis
    emit_neff = TRUE,             # Output effective sample counts
    emit_jackknife = FALSE        # Jackknife analysis (set TRUE for robustness check)
  ),
  
  # Non-RET exploration
  exploratory_nonRET = TRUE,      # Analyze non-RET samples (interpretation limited)
  nonRET_note = "Non-RET analysis is exploratory only",
  
  # ============================================================================
  # Legacy Parameters (unused in v7 pipeline)
  # ============================================================================
  
  # These are kept for backward compatibility but not used
  # legacy_param1 = NULL,         # legacy (unused in v7 pipeline)
  
  # ============================================================================
  # Output Control
  # ============================================================================
  
  # File naming
  output_prefix = "reo_v7",       # Prefix for output files
  
  # CSV export settings
  export_csv = TRUE,              # Export results as CSV
  csv_digits = 4,                 # Decimal places in CSV output
  
  # Plotting
  make_plots = TRUE,              # Generate diagnostic plots
  plot_format = "pdf",            # Plot file format
  plot_width = 10,                # Plot width in inches
  plot_height = 8                 # Plot height in inches
)

# ============================================================================
# Validation and Logging
# ============================================================================

# Validate critical parameters
stopifnot(
  CONFIG$dead_zone == log2(1.2),  # Dead zone must be fixed
  CONFIG$n_eff_threshold >= 0.5 && CONFIG$n_eff_threshold <= 1.0,
  CONFIG$ci_level > 0 && CONFIG$ci_level < 1,
  CONFIG$cor_threshold >= 0 && CONFIG$cor_threshold <= 1,
  CONFIG$panel_defaults$size_min <= CONFIG$panel_defaults$size_target,
  CONFIG$panel_defaults$size_target <= CONFIG$panel_defaults$size_max
)

# Log configuration at load time (for transparency)
if (CONFIG$VERBOSE) {
  cat("\n=== REO Configuration Loaded (v7.0) ===\n")
  cat(sprintf("  Date: %s\n", CONFIG$config_date))
  cat("\n--- Core Parameters ---\n")
  cat(sprintf("  tau_strength: %.4f\n", CONFIG$tau_strength))
  cat(sprintf("  dead_zone: %.4f (FIXED)\n", CONFIG$dead_zone))
  cat(sprintf("  n_eff_threshold: %.2f\n", CONFIG$n_eff_threshold))
  cat(sprintf("  trim_low: %.2f (adjusted from 0.10)\n", CONFIG$trim_low))
  cat(sprintf("  trim_high: %.2f\n", CONFIG$trim_high))
  cat(sprintf("  M: %d (adjusted from 10)\n", CONFIG$M))
  
  cat("\n--- CI Settings ---\n")
  cat(sprintf("  ci_method: %s\n", CONFIG$ci_method))
  cat(sprintf("  ci_level: %.2f\n", CONFIG$ci_level))
  
  cat("\n--- POC Handling ---\n")
  cat(sprintf("  %s\n", CONFIG$poc_note))
  cat(sprintf("  Cutpoints: %s\n", paste(CONFIG$poc_cutpoints, collapse=", ")))
  
  cat("\n--- Step-Specific Settings ---\n")
  cat(sprintf("  Step 3: JT diagnostics = %s\n", CONFIG$emit_jt_diagnostic))
  cat(sprintf("  Step 4: Correlation threshold = %.2f\n", CONFIG$cor_threshold))
  cat(sprintf("  Step 5: Panel size = %d-%d (target %d)\n", 
              CONFIG$panel_defaults$size_min, 
              CONFIG$panel_defaults$size_max,
              CONFIG$panel_defaults$size_target))
  cat(sprintf("  Step 6: Main analysis = %d layers\n", CONFIG$rtest$main_layers))
  
  cat("========================================\n\n")
}

# ============================================================================
# Helper Functions
# ============================================================================

#' Get parameter with fallback
#' @param param_name Parameter name
#' @param default Default value if not found
get_config_param <- function(param_name, default = NULL) {
  if (param_name %in% names(CONFIG)) {
    return(CONFIG[[param_name]])
  } else {
    if (!is.null(default)) {
      warning(sprintf("Parameter '%s' not found in CONFIG, using default: %s", 
                      param_name, as.character(default)))
      return(default)
    } else {
      stop(sprintf("Required parameter '%s' not found in CONFIG", param_name))
    }
  }
}

#' Export CONFIG to character vector for logging
export_config_log <- function() {
  lines <- c(
    "=== CONFIG SNAPSHOT ===",
    sprintf("Version: %s", CONFIG$config_version),
    sprintf("Date: %s", CONFIG$config_date),
    "",
    "Core Parameters:",
    sprintf("  tau_strength: %.4f", CONFIG$tau_strength),
    sprintf("  dead_zone: %.4f", CONFIG$dead_zone),
    sprintf("  n_eff_threshold: %.2f", CONFIG$n_eff_threshold),
    sprintf("  ci_method: %s", CONFIG$ci_method),
    sprintf("  ci_level: %.2f", CONFIG$ci_level),
    "",
    "Processing:",
    sprintf("  trim_low: %.2f", CONFIG$trim_low),
    sprintf("  trim_high: %.2f", CONFIG$trim_high),
    sprintf("  M: %d", CONFIG$M),
    sprintf("  cor_threshold: %.2f", CONFIG$cor_threshold),
    "",
    "Panel:",
    sprintf("  size: %d-%d (target %d)", 
            CONFIG$panel_defaults$size_min,
            CONFIG$panel_defaults$size_max,
            CONFIG$panel_defaults$size_target),
    sprintf("  QC anchors: %d", CONFIG$qc_anchor$n),
    "",
    "R_test:",
    sprintf("  layers: %d", CONFIG$rtest$main_layers),
    sprintf("  min_n: %d", CONFIG$rtest$min_n_stratum),
    "======================="
  )
  return(lines)
}

# ============================================================================
# Notes for Implementation
# ============================================================================

# 1. Step 1 reads this file and saves CONFIG to RDS
# 2. Steps 2-6 use CONFIG from RDS (not this file)
# 3. Each step's RDS includes CONFIG snapshot for reproducibility
# 4. Modifications to CONFIG should bump config_version
# 5. Dead zone is FIXED at log2(1.2) - do not adjust
# 6. POC is for stratification only, not for REO logic