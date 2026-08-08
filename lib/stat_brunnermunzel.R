# Studentized permutation Brunner-Munzel test, by exact enumeration where that
# is affordable and by Monte Carlo otherwise.
#
# The permutation null depends on the data only through the tie pattern, so
# whenever C(n, nx) is small enough to enumerate, every allocation can be
# visited once and the p-value is exact -- no sampling error and no 1/(B+1)
# floor. `method = "auto"` (the default) takes that route when
# C(n, nx) <= getOption("brunnermunzel.exact.max.allocations", 1e8) and falls
# back to Monte Carlo above it.
#
# `brunnermunzel_pvalues()` applies the same machinery to a whole matrix of
# features at once, which is the efficient way to use it: features sharing a
# tie pattern share one enumeration or one Monte Carlo null.
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
      "Cannot determine the location of lib/stat_brunnermunzel.R; ",
      "load it with source().",
      call. = FALSE
    )
  }
  normalizePath(tail(candidates, 1L), mustWork = TRUE)
}

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop(
    "Package 'Rcpp' is required to source lib/stat_brunnermunzel.R.",
    call. = FALSE
  )
}

.bm_source_file <- .bm_find_source_file()
.bm_cpp_file <- file.path(dirname(.bm_source_file), "stat_brunnermunzel.cpp")
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
  "bm_mc_extreme_cpp",
  "bm_n_allocations_cpp",
  "bm_observed_matrix_cpp",
  "bm_exact_extreme_cpp"
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

.bm_exact_max_allocations <- function() {
  value <- getOption("brunnermunzel.exact.max.allocations", 1e8)
  if (
    length(value) != 1L ||
      is.logical(value) ||
      is.object(value) ||
      !typeof(value) %in% c("integer", "double") ||
      is.na(value) ||
      value < 0
  ) {
    stop(
      "Option 'brunnermunzel.exact.max.allocations' must be a ",
      "non-negative number.",
      call. = FALSE
    )
  }
  as.double(value)
}

.bm_exact_threads <- function() {
  value <- getOption("brunnermunzel.exact.threads", NULL)
  if (is.null(value)) {
    cores <- tryCatch(
      parallel::detectCores(logical = FALSE),
      error = function(condition) NA_integer_
    )
    if (length(cores) != 1L || is.na(cores) || cores < 1L) {
      cores <- 1L
    }
    # Capped so that sourcing this file never oversubscribes a shared machine
    # by default; raise the option deliberately if more threads are wanted.
    return(max(1L, min(16L, as.integer(cores))))
  }
  .bm_validate_integer(
    value,
    "Option 'brunnermunzel.exact.threads'",
    positive = TRUE
  )
}

.bm_alternative_code <- function(alternative) {
  switch(alternative, two.sided = 1L, greater = 2L, less = 3L)
}

# Decides between enumeration and sampling, and reports the number of
# allocations so callers can record how the p-value was obtained.
.bm_resolve_method <- function(method, n, nx) {
  total <- .brunnermunzel_mc_state$bm_n_allocations_cpp(
    as.integer(n),
    as.integer(nx)
  )
  if (method == "auto") {
    method <- if (is.finite(total) && total <= .bm_exact_max_allocations()) {
      "exact"
    } else {
      "mc"
    }
  } else if (method == "exact" && !is.finite(total)) {
    stop(
      "Exact enumeration requires C(n, nx) to be exactly representable, ",
      "which it is not for n = ", n, " and nx = ", nx, ".",
      call. = FALSE
    )
  }
  list(method = method, total = total)
}

.bm_exact_counts <- function(block_sizes, nx, observed, alternative_code) {
  .brunnermunzel_mc_state$bm_exact_extreme_cpp(
    as.integer(block_sizes),
    as.integer(nx),
    as.numeric(observed),
    alternative_code,
    .bm_exact_threads()
  )
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

# A null for one tie pattern, from the cache when an explicit seed makes it
# reproducible. Without a seed every call draws afresh and bypasses the cache.
.bm_null_for <- function(block_sizes, nx, ny, B, seed) {
  if (is.null(seed)) {
    return(.bm_generate_null(block_sizes, nx, B))
  }
  cache_key <- .bm_cache_key(block_sizes, nx, ny, B, seed)
  null_distribution <- .bm_cache_get(cache_key)
  if (is.null(null_distribution)) {
    null_distribution <- .bm_with_local_seed(seed, function() {
      .bm_generate_null(block_sizes, nx, B)
    })
    .bm_cache_put(cache_key, null_distribution)
  }
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
    method = c("auto", "exact", "mc"),
    ...) {
  alternative <- match.arg(alternative)
  est <- match.arg(est)
  method <- match.arg(method)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  .bm_validate_samples(x, y)
  B <- .bm_validate_integer(B, "B", positive = TRUE)
  if (!is.null(seed)) {
    seed <- .bm_validate_integer(seed, "seed")
  }

  # `force` is deliberately not evaluated. It is retained solely for call
  # compatibility with brunnermunzel 2.0; neither execution path ever falls
  # back to the asymptotic test.

  alternative_code <- .bm_alternative_code(alternative)
  nx <- length(x)
  ny <- length(y)
  observed <- .brunnermunzel_mc_state$bm_observed_cpp(
    as.numeric(c(x, y)),
    as.integer(nx)
  )
  resolved <- .bm_resolve_method(method, nx + ny, nx)

  if (resolved$method == "exact") {
    counted <- .bm_exact_counts(
      observed$block_sizes,
      nx,
      observed$statistic,
      alternative_code
    )
    n_extreme <- counted$count
    n_permutations <- counted$total
    # The observed allocation is itself enumerated and is always extreme
    # against itself, so no plus-one correction is needed or wanted here.
    PVALUE <- n_extreme / n_permutations
  } else {
    null_distribution <- .bm_null_for(observed$block_sizes, nx, ny, B, seed)
    n_extreme <- .brunnermunzel_mc_state$bm_mc_extreme_cpp(
      null_distribution,
      observed$statistic,
      as.integer(nx),
      as.integer(ny),
      alternative_code
    )
    n_permutations <- as.double(B)
    PVALUE <- (n_extreme + 1) / (n_permutations + 1)
  }

  effect <- .bm_effect(x, y)
  if (est == "original") {
    ESTIMATE <- effect
    names(ESTIMATE) <- "P(X<Y)+.5*P(X=Y)"
  } else {
    ESTIMATE <- 2 * effect - 1
    names(ESTIMATE) <- "P(X<Y)-P(X>Y)"
  }

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
    method = resolved$method,
    n.permutations = n_permutations,
    B = if (resolved$method == "mc") B else NA_integer_,
    seed = if (resolved$method == "mc") seed else NULL,
    n.extreme = n_extreme,
    statistic = observed$statistic
  )
  result
}

.bm_validate_matrix <- function(data, nx) {
  if (
    !is.matrix(data) ||
      !is.numeric(data) ||
      is.object(data) ||
      !typeof(data) %in% c("integer", "double")
  ) {
    stop("data must be a plain numeric matrix.", call. = FALSE)
  }
  if (anyNA(data) || !all(is.finite(data))) {
    stop(
      "data must contain only finite values (no NA, NaN, or Inf).",
      call. = FALSE
    )
  }
  nx <- .bm_validate_integer(nx, "nx", positive = TRUE)
  if (nx < 2L || ncol(data) - nx < 2L) {
    stop("each group must contain at least two observations", call. = FALSE)
  }
  nx
}

# Standard Brunner-Munzel statistic for every row of `data`. Columns 1..nx are
# group x. Cheap relative to brunnermunzel_pvalues() because no permutation
# null is built, which is what makes label-shuffling loops affordable.
brunnermunzel_statistics <- function(data, nx) {
  nx <- .bm_validate_matrix(data, nx)
  if (nrow(data) == 0L) {
    return(stats::setNames(numeric(0), rownames(data)))
  }
  statistics <- .brunnermunzel_mc_state$bm_observed_matrix_cpp(
    matrix(as.double(data), nrow = nrow(data)),
    as.integer(nx)
  )$statistic
  names(statistics) <- rownames(data)
  statistics
}

# The Brunner-Munzel effect P(X<Y) + 0.5 P(X=Y) for every row of `data`, on the
# same column convention as brunnermunzel_statistics().
brunnermunzel_effects <- function(data, nx) {
  nx <- .bm_validate_matrix(data, nx)
  if (nrow(data) == 0L) {
    return(stats::setNames(numeric(0), rownames(data)))
  }
  ny <- ncol(data) - nx
  pooled <- apply(data, 1L, rank)
  effects <- (colMeans(pooled[nx + seq_len(ny), , drop = FALSE]) -
                (ny + 1) / 2) / nx
  names(effects) <- rownames(data)
  effects
}

# Two-group p-values for every row of `data` at once. Columns 1..nx are group
# x and the remainder are group y. Rows sharing a tie pattern share a single
# enumeration (or a single Monte Carlo null), which is what makes this far
# cheaper than looping brunnermunzel_mc_test() over the rows.
brunnermunzel_pvalues <- function(
    data,
    nx,
    alternative = c("two.sided", "greater", "less"),
    B = 999999L,
    seed = NULL,
    method = c("auto", "exact", "mc")) {
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  nx <- .bm_validate_matrix(data, nx)
  n <- ncol(data)
  ny <- n - nx
  if (!is.null(seed)) {
    seed <- .bm_validate_integer(seed, "seed")
  }
  if (nrow(data) == 0L) {
    return(stats::setNames(numeric(0), rownames(data)))
  }

  alternative_code <- .bm_alternative_code(alternative)
  resolved <- .bm_resolve_method(method, n, nx)
  # B parameterizes only the Monte Carlo null; on the exact path it is unused
  # and deliberately not validated, so callers on that path need no dummy.
  if (resolved$method == "mc") {
    B <- .bm_validate_integer(B, "B", positive = TRUE)
  }
  observed <- .brunnermunzel_mc_state$bm_observed_matrix_cpp(
    matrix(as.double(data), nrow = nrow(data)),
    as.integer(nx)
  )

  pvalues <- numeric(nrow(data))
  for (rows in split(seq_len(nrow(data)), observed$pattern)) {
    block_sizes <- as.integer(
      strsplit(observed$pattern[[rows[1L]]], ",", fixed = TRUE)[[1L]]
    )
    statistics <- observed$statistic[rows]
    if (resolved$method == "exact") {
      counted <- .bm_exact_counts(
        block_sizes, nx, statistics, alternative_code
      )
      pvalues[rows] <- counted$count / counted$total
    } else {
      null_distribution <- .bm_null_for(block_sizes, nx, ny, B, seed)
      n_extreme <- vapply(
        statistics,
        function(statistic) {
          .brunnermunzel_mc_state$bm_mc_extreme_cpp(
            null_distribution,
            statistic,
            as.integer(nx),
            as.integer(ny),
            alternative_code
          )
        },
        numeric(1)
      )
      pvalues[rows] <- (n_extreme + 1) / (as.double(B) + 1)
    }
  }

  attr(pvalues, "mc") <- list(
    method = resolved$method,
    n.permutations = if (resolved$method == "exact") resolved$total else as.double(B),
    statistic = observed$statistic
  )
  names(pvalues) <- rownames(data)
  pvalues
}

rm(
  .bm_source_file,
  .bm_cpp_file,
  .bm_cpp_signature,
  .bm_required_cpp,
  .bm_needs_compile,
  .bm_find_source_file
)
