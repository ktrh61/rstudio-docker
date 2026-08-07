# =============================================================================
# dose_thyroid_as.R
#
# Standalone Monte Carlo estimator of the thyroid Assigned Share (AS) associated
# with the expected value of the excess relative risk (ERR), approximating the
# NCI/NIH Interactive RadioEpidemiological Program (IREP) for a single case.
#
# Scope (fixed):
#   - cancer site : thyroid
#   - exposure    : acute
#   - radiation   : electrons > 15 keV (treated as reference radiation, REF = 1)
#   - dose input  : constant (point estimate), no dose-distribution sampling
#
# The estimator reproduces IREP's "Assigned Share associated with the expected
# value of ERR" (AS = mean(ERR) / (mean(ERR) + 1)). It is an approximation of
# the same algorithm, not a bit-for-bit reproduction of IREP output.
#
# Sources:
#   Kocher DC, et al. (2008) Health Phys 95(1):119-147.
#   NCI-CDC Working Group (2003) NIH Publication 03-5387 (Land et al.).
#
# Dependencies: base R only (stats::approxfun, stats::splinefun). No packages,
# no project files, no global state. Reads nothing from disk. Pure function of
# its arguments; the global RNG state is saved and restored on exit so that
# calling the function has no side effect on the caller's random stream.
# =============================================================================


# --- Model constants ---------------------------------------------------------

# Thyroid ERR/Sv statistical uncertainty, as a lognormal distribution
# (geometric mean GM, geometric standard deviation GSD) indexed by age at
# exposure. Values are Table IV.D.8 of NIH (2003), verbatim.
.THYROID_ERRSV <- data.frame(
  age = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
  gm = c(
    9.463, 6.262, 4.136, 2.732, 1.804, 1.192,
    0.788, 0.521, 0.345, 0.228, 0.151
  ),
  gsd = c(
    2.183, 1.924, 1.976, 2.160, 2.301, 2.367,
    2.365, 2.379, 2.732, 3.140, 3.611
  )
)

# Discrete DDREF distribution for breast and thyroid under chronic low-LET
# exposure (Kocher 2008, Fig. 5b; NIH 2003, Fig. IV.F.2). Mean = 1.6.
.DDREF_VALUES <- c(0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0)
.DDREF_PROBS <- c(0.01, 0.04, 0.35, 0.23, 0.23, 0.10, 0.04)

# Reference dose D_L (Gy) below which an acute-exposure DDREF is phased in,
# log-uniform on this interval (Kocher 2008, Eq. 25).
.DL_RANGE_GY <- c(0.03, 0.20)

# Latency: sigmoid midpoint tau ~ Triangular(a, mode, b), thyroid NIOSH-IREP.
.LATENCY_TAU <- c(a = 4, mode = 4.5, b = 5.5)

# Age-at-exposure interpolators for the ERR/Sv lognormal parameters.
# log(GM) is linear in age; log(GSD) is curvilinear, so a monotone spline is
# used between tabulated points. Ages above 50 are clipped to 50.
.gm_logfun <- stats::approxfun(.THYROID_ERRSV$age, log(.THYROID_ERRSV$gm),
  rule = 2
)
.gsd_logfun <- stats::splinefun(.THYROID_ERRSV$age, log(.THYROID_ERRSV$gsd),
  method = "monoH.FC"
)


# --- Internal helpers --------------------------------------------------------

# Lognormal ERR/Sv parameters (GM, GSD) at a given age at exposure.
.thyroid_errsv_params <- function(age_exposure) {
  e <- min(age_exposure, 50)
  list(gm = exp(.gm_logfun(e)), gsd = exp(.gsd_logfun(e)))
}

# Sample a lognormal parameterised by geometric mean and geometric SD.
.rlnorm_gmgsd <- function(n, gm, gsd) {
  stats::rlnorm(n, meanlog = log(gm), sdlog = log(gsd))
}

# Sample the effective acute-exposure DDREF for a thyroid case at a given dose.
# A DDREF is applied only below an uncertain reference dose D_L; at or above
# D_L it is 1. Below D_L the chronic DDREF is phased in via a logistic function
# of dose (Kocher 2008, Eq. 25).
.sample_ddref_acute <- function(n, dose_Gy) {
  ddref_chronic <- sample(.DDREF_VALUES,
    size = n, replace = TRUE,
    prob = .DDREF_PROBS
  )
  dl <- exp(stats::runif(n, log(.DL_RANGE_GY[1]), log(.DL_RANGE_GY[2])))
  di <- 0.5 * dl
  s <- di / log(500)
  phased <- 1 + (1 - 1 / ddref_chronic) / (1 + exp((dose_Gy - di) / s))
  ddref <- 1 / phased
  ddref[dose_Gy >= dl] <- 1
  ddref
}

# Sample from Triangular(a, mode, b).
.rtriangle <- function(n, a, mode, b) {
  u <- stats::runif(n)
  fc <- (mode - a) / (b - a)
  ifelse(u < fc,
    a + sqrt(u * (b - a) * (mode - a)),
    b - sqrt((1 - u) * (b - a) * (b - mode))
  )
}

# Latency adjustment: sigmoid in time since exposure, uncertain midpoint tau.
# Approaches 1 for times well beyond the minimum latency period.
.latency_adjust <- function(n, years_since_exposure) {
  tau <- .rtriangle(n, .LATENCY_TAU["a"], .LATENCY_TAU["mode"], .LATENCY_TAU["b"])
  shape <- (7 - 2) / (2 * log(0.99 / 0.01))
  1 / (1 + exp((tau - years_since_exposure) / shape))
}


# --- Public function ---------------------------------------------------------

#' Thyroid Assigned Share at the expected value of ERR (single case)
#'
#' Estimates the thyroid Assigned Share associated with the expected value of
#' the excess relative risk for one exposed individual, under acute exposure to
#' electrons > 15 keV with a constant (point-estimate) dose. Approximates the
#' corresponding IREP output.
#'
#' @param dose_mGy       Thyroid absorbed dose, in mGy (point estimate).
#' @param age_exposure   Age at exposure, in years (clipped at 50 internally).
#' @param age_surgery    Age at diagnosis/surgery, in years. Used with
#'                       age_exposure to obtain time since exposure for the
#'                       latency adjustment.
#' @param n_iter         Number of Monte Carlo iterations (default 10000).
#' @param seed           RNG seed for reproducibility (default 99, IREP's
#'                       default). The global RNG state is restored on exit.
#'
#' @return A single numeric value: the Assigned Share associated with the
#'   expected value of ERR, expressed as a percentage (0-100).
compute_thyroid_as <- function(dose_mGy, age_exposure, age_surgery,
                               n_iter = 10000, seed = 99) {
  # Save and restore the caller's RNG state so this function has no side effect.
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)

  dose_Gy <- dose_mGy / 1000 # electrons: wR = 1, so Sv = Gy
  years_since_exposure <- age_surgery - age_exposure

  params <- .thyroid_errsv_params(age_exposure)
  errsv <- .rlnorm_gmgsd(n_iter, params$gm, params$gsd)
  ddref <- .sample_ddref_acute(n_iter, dose_Gy)
  lat <- .latency_adjust(n_iter, years_since_exposure)

  # REF_L = 1 (electrons as reference radiation); no A-bomb dosimetry-bias
  # factor (that correction is specific to DS86 dose reconstruction and does
  # not apply to the Ron et al. 1995 pooled thyroid model).
  err <- (errsv / ddref) * dose_Gy * lat
  err[err < 0] <- 0

  err_mean <- mean(err)
  100 * err_mean / (err_mean + 1)
}
