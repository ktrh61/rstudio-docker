# Monte Carlo studentized permutation Brunner-Munzel test.
#
# This file is intentionally repository-local rather than an R package. Source
# it from any working directory; the companion C++ file is resolved relative to
# this file's own location.

.bm_find_source_file <- function() {
  frames <- sys.frames()
  candidates <- unlist(
    lapply(frames, function(frame) {
      ofile <- frame$ofile
      if (is.null(ofile) || length(ofile) != 1L || !nzchar(ofile)) {
        return(character())
      }
      ofile
    }),
    use.names = FALSE
  )
  if (length(candidates) == 0L) {
    stop(
      "Cannot determine the location of utils/brunnermunzel_mc.R; ",
      "load it with source().",
      call. = FALSE
    )
  }
  normalizePath(tail(candidates, 1L), mustWork = TRUE)
}

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop(
    "Package 'Rcpp' is required to source utils/brunnermunzel_mc.R.",
    call. = FALSE
  )
}

.bm_source_file <- .bm_find_source_file()
.bm_cpp_file <- file.path(dirname(.bm_source_file), "brunnermunzel_mc.cpp")
if (!file.exists(.bm_cpp_file)) {
  stop("Companion C++ file not found: ", .bm_cpp_file, call. = FALSE)
}
.bm_cpp_file <- normalizePath(.bm_cpp_file, mustWork = TRUE)
.bm_cpp_signature <- unname(tools::md5sum(.bm_cpp_file))

if (
  !exists(".brunnermunzel_mc_state", inherits = FALSE) ||
    !is.environment(.brunnermunzel_mc_state)
) {
  .brunnermunzel_mc_state <- new.env(parent = baseenv())
}

.bm_required_cpp <- c(
  "bm_observed_cpp",
  "bm_mc_null_cpp",
  "bm_mc_extreme_cpp"
)
.bm_needs_compile <-
  !all(vapply(
    .bm_required_cpp,
    exists,
    logical(1),
    envir = .brunnermunzel_mc_state,
    inherits = FALSE
  )) ||
  !identical(.brunnermunzel_mc_state$cpp_file, .bm_cpp_file) ||
  !identical(.brunnermunzel_mc_state$cpp_signature, .bm_cpp_signature)

if (.bm_needs_compile) {
  Rcpp::sourceCpp(
    file = .bm_cpp_file,
    env = .brunnermunzel_mc_state,
    rebuild = FALSE,
    showOutput = FALSE,
    verbose = FALSE
  )
  .brunnermunzel_mc_state$cpp_file <- .bm_cpp_file
  .brunnermunzel_mc_state$cpp_signature <- .bm_cpp_signature
  .brunnermunzel_mc_state$null_cache <- new.env(parent = emptyenv())
  .brunnermunzel_mc_state$cache_bytes <- 0
  .brunnermunzel_mc_state$cache_clock <- 0
  .brunnermunzel_mc_state$null_generation_count <- 0
}

if (
  !exists("null_cache", envir = .brunnermunzel_mc_state, inherits = FALSE) ||
    !is.environment(.brunnermunzel_mc_state$null_cache)
) {
  .brunnermunzel_mc_state$null_cache <- new.env(parent = emptyenv())
  .brunnermunzel_mc_state$cache_bytes <- 0
  .brunnermunzel_mc_state$cache_clock <- 0
  .brunnermunzel_mc_state$null_generation_count <- 0
}

.bm_validate_samples <- function(x, y) {
  is_plain_numeric_vector <- function(value) {
    is.numeric(value) &&
      !is.object(value) &&
      typeof(value) %in% c("integer", "double") &&
      is.null(dim(value))
  }
  if (!is_plain_numeric_vector(x) || !is_plain_numeric_vector(y)) {
    stop("x and y must be numeric vectors.", call. = FALSE)
  }
  if (length(x) == 0L || length(y) == 0L) {
    stop("x and y must not be empty.", call. = FALSE)
  }
  if (length(x) < 2L || length(y) < 2L) {
    stop("not enough observations", call. = FALSE)
  }
  if (anyNA(x) || anyNA(y) || !all(is.finite(x)) || !all(is.finite(y))) {
    stop("x and y must contain only finite values (no NA, NaN, or Inf).",
         call. = FALSE)
  }
  invisible(NULL)
}

.bm_validate_integer <- function(value, name, positive = FALSE) {
  if (
    length(value) != 1L ||
      is.logical(value) ||
      is.object(value) ||
      !typeof(value) %in% c("integer", "double") ||
      !is.finite(value) ||
      value != trunc(value)
  ) {
    stop(name, " must be a finite integer of length one.", call. = FALSE)
  }
  if (positive && value <= 0) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  if (value > .Machine$integer.max || value < -.Machine$integer.max) {
    stop(name, " is outside the range supported by R integers.", call. = FALSE)
  }
  as.integer(value)
}

.bm_effect <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  pooled_rank <- rank(c(x, y), ties.method = "average")
  (mean(pooled_rank[nx + seq_len(ny)]) - (ny + 1) / 2) / nx
}

.bm_statistic <- function(x, y) {
  .bm_validate_samples(x, y)
  nx <- length(x)
  ny <- length(y)
  n <- nx + ny
  pooled_rank <- rank(c(x, y), ties.method = "average")
  x_index <- seq_len(nx)
  y_index <- nx + seq_len(ny)
  mean_rank_x <- mean(pooled_rank[x_index])
  mean_rank_y <- mean(pooled_rank[y_index])
  d_x <- pooled_rank[x_index] - rank(x, ties.method = "average")
  d_y <- pooled_rank[y_index] - rank(y, ties.method = "average")
  s_x <- sum((d_x - mean(d_x))^2) / (nx - 1)
  s_y <- sum((d_y - mean(d_y))^2) / (ny - 1)
  numerator <- nx * ny * (mean_rank_y - mean_rank_x)
  denominator <- n * sqrt(nx * s_x + ny * s_y)

  if (denominator > 0) {
    return(numerator / denominator)
  }
  if (numerator > 0) {
    return(Inf)
  }
  if (numerator < 0) {
    return(-Inf)
  }
  0
}

.bm_with_local_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(code())
  }

  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  code()
}

.bm_cache_max_bytes <- function() {
  value <- getOption(
    "brunnermunzel.mc.cache.max.bytes",
    256 * 1024^2
  )
  if (
    length(value) != 1L ||
      is.logical(value) ||
      is.object(value) ||
      !typeof(value) %in% c("integer", "double") ||
      !is.finite(value) ||
      value < 0
  ) {
    stop(
      "Option 'brunnermunzel.mc.cache.max.bytes' must be a ",
      "non-negative finite number.",
      call. = FALSE
    )
  }
  as.double(value)
}

.bm_cache_key <- function(block_sizes, nx, ny, B, seed) {
  paste(
    .brunnermunzel_mc_state$cpp_signature,
    nx,
    ny,
    B,
    seed,
    paste(RNGkind(), collapse = ","),
    paste(block_sizes, collapse = ","),
    sep = "|"
  )
}

.bm_cache_get <- function(key) {
  cache <- .brunnermunzel_mc_state$null_cache
  if (!exists(key, envir = cache, inherits = FALSE)) {
    return(NULL)
  }
  entry <- get(key, envir = cache, inherits = FALSE)
  .brunnermunzel_mc_state$cache_clock <-
    .brunnermunzel_mc_state$cache_clock + 1
  entry$last_used <- .brunnermunzel_mc_state$cache_clock
  assign(key, entry, envir = cache)
  entry$null
}

.bm_cache_put <- function(key, null_distribution) {
  max_bytes <- .bm_cache_max_bytes()
  entry <- list(
    null = null_distribution,
    bytes = 0,
    last_used = .brunnermunzel_mc_state$cache_clock + 1
  )
  bytes <- as.double(object.size(entry)) + as.double(object.size(key))
  entry$bytes <- bytes
  if (max_bytes == 0 || bytes > max_bytes) {
    return(invisible(FALSE))
  }

  cache <- .brunnermunzel_mc_state$null_cache
  if (exists(key, envir = cache, inherits = FALSE)) {
    old_entry <- get(key, envir = cache, inherits = FALSE)
    rm(list = key, envir = cache)
    .brunnermunzel_mc_state$cache_bytes <-
      .brunnermunzel_mc_state$cache_bytes - old_entry$bytes
  }

  while (
    .brunnermunzel_mc_state$cache_bytes + bytes > max_bytes &&
      length(cache_keys <- ls(cache, all.names = TRUE)) > 0L
  ) {
    last_used <- vapply(cache_keys, function(cache_key) {
      get(cache_key, envir = cache, inherits = FALSE)$last_used
    }, numeric(1))
    victim <- cache_keys[[which.min(last_used)]]
    victim_entry <- get(victim, envir = cache, inherits = FALSE)
    rm(list = victim, envir = cache)
    .brunnermunzel_mc_state$cache_bytes <-
      .brunnermunzel_mc_state$cache_bytes - victim_entry$bytes
  }

  .brunnermunzel_mc_state$cache_clock <-
    .brunnermunzel_mc_state$cache_clock + 1
  entry$last_used <- .brunnermunzel_mc_state$cache_clock
  assign(key, entry, envir = cache)
  .brunnermunzel_mc_state$cache_bytes <-
    .brunnermunzel_mc_state$cache_bytes + bytes
  invisible(TRUE)
}

.bm_generate_null <- function(block_sizes, nx, B) {
  null_distribution <- .brunnermunzel_mc_state$bm_mc_null_cpp(
    as.integer(block_sizes),
    as.integer(nx),
    B
  )
  .brunnermunzel_mc_state$null_generation_count <-
    .brunnermunzel_mc_state$null_generation_count + 1
  null_distribution
}

brunnermunzel_mc_test <- function(
    x,
    y,
    alternative = c("two.sided", "greater", "less"),
    force = FALSE,
    est = c("original", "difference"),
    B = 999999L,
    seed = NULL,
    ...) {
  alternative <- match.arg(alternative)
  est <- match.arg(est)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  .bm_validate_samples(x, y)
  B <- .bm_validate_integer(B, "B", positive = TRUE)
  if (!is.null(seed)) {
    seed <- .bm_validate_integer(seed, "seed")
  }

  # `force` is deliberately not evaluated. It is retained solely for call
  # compatibility with brunnermunzel 2.0; Monte Carlo execution never falls
  # back to the asymptotic test.

  alternative_code <- switch(
    alternative,
    two.sided = 1L,
    greater = 2L,
    less = 3L
  )
  nx <- length(x)
  ny <- length(y)
  observed <- .brunnermunzel_mc_state$bm_observed_cpp(
    as.numeric(c(x, y)),
    as.integer(nx)
  )

  if (is.null(seed)) {
    null_distribution <- .bm_generate_null(observed$block_sizes, nx, B)
  } else {
    cache_key <- .bm_cache_key(observed$block_sizes, nx, ny, B, seed)
    null_distribution <- .bm_cache_get(cache_key)
    if (is.null(null_distribution)) {
      null_distribution <- .bm_with_local_seed(seed, function() {
        .bm_generate_null(observed$block_sizes, nx, B)
      })
      .bm_cache_put(cache_key, null_distribution)
    }
  }
  n_extreme <- .brunnermunzel_mc_state$bm_mc_extreme_cpp(
    null_distribution,
    observed$statistic,
    as.integer(nx),
    as.integer(ny),
    alternative_code
  )

  effect <- .bm_effect(x, y)
  if (est == "original") {
    ESTIMATE <- effect
    names(ESTIMATE) <- "P(X<Y)+.5*P(X=Y)"
  } else {
    ESTIMATE <- 2 * effect - 1
    names(ESTIMATE) <- "P(X<Y)-P(X>Y)"
  }
  PVALUE <- (n_extreme + 1) / (as.double(B) + 1)

  result <- structure(
    list(
      method = "permuted Brunner-Munzel Test",
      data.name = DNAME,
      p.value = PVALUE,
      estimate = ESTIMATE
    ),
    class = "htest"
  )
  attr(result, "mc") <- list(
    B = B,
    seed = seed,
    n.extreme = n_extreme,
    statistic = observed$statistic
  )
  result
}

rm(
  .bm_source_file,
  .bm_cpp_file,
  .bm_cpp_signature,
  .bm_required_cpp,
  .bm_needs_compile,
  .bm_find_source_file
)
