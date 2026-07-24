.repo_root <- normalizePath(
  file.path(testthat::test_path(), "..", ".."),
  mustWork = TRUE
)
.implementation_file <- file.path(.repo_root, "utils", "brunnermunzel_mc.R")
source(.implementation_file)

.reset_null_cache <- function() {
  cache <- .brunnermunzel_mc_state$null_cache
  cache_keys <- ls(cache, all.names = TRUE)
  if (length(cache_keys) > 0L) {
    rm(list = cache_keys, envir = cache)
  }
  .brunnermunzel_mc_state$cache_bytes <- 0
  .brunnermunzel_mc_state$cache_clock <- 0
  .brunnermunzel_mc_state$null_generation_count <- 0
  invisible(NULL)
}

.reference_bm <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  n <- nx + ny
  pooled <- rank(c(x, y), ties.method = "average")
  within_x <- rank(x, ties.method = "average")
  within_y <- rank(y, ties.method = "average")
  x_index <- seq_len(nx)
  y_index <- nx + seq_len(ny)
  mean_x <- mean(pooled[x_index])
  mean_y <- mean(pooled[y_index])
  d_x <- pooled[x_index] - within_x
  d_y <- pooled[y_index] - within_y
  s_x <- sum((d_x - mean(d_x))^2) / (nx - 1)
  s_y <- sum((d_y - mean(d_y))^2) / (ny - 1)
  numerator <- nx * ny * (mean_y - mean_x)
  denominator <- n * sqrt(nx * s_x + ny * s_y)
  statistic <- if (denominator > 0) {
    numerator / denominator
  } else if (numerator > 0) {
    Inf
  } else if (numerator < 0) {
    -Inf
  } else {
    0
  }

  list(
    pooled.rank = pooled,
    within.rank.x = within_x,
    within.rank.y = within_y,
    variance.x = s_x,
    variance.y = s_y,
    statistic = statistic,
    effect = (mean_y - (ny + 1) / 2) / nx
  )
}

.enumerate_bm <- function(x, y, alternative) {
  nx <- length(x)
  combined <- c(x, y)
  allocations <- combn(seq_along(combined), nx, simplify = FALSE)
  statistics <- vapply(allocations, function(index) {
    .reference_bm(combined[index], combined[-index])$statistic
  }, numeric(1))
  observed <- .reference_bm(x, y)$statistic
  extreme <- switch(
    alternative,
    two.sided = abs(statistics) >= abs(observed),
    greater = statistics <= observed,
    less = statistics >= observed
  )
  list(
    statistic = observed,
    permutation.statistics = statistics,
    n.extreme = sum(extreme),
    p.value = mean(extreme)
  )
}

testthat::test_that("source loads only the implementation and preserves cwd", {
  old_wd <- getwd()
  root_files <- list.files(.repo_root, all.files = TRUE)
  env <- new.env(parent = globalenv())

  source(.implementation_file, local = env)

  testthat::expect_identical(getwd(), old_wd)
  testthat::expect_identical(list.files(.repo_root, all.files = TRUE), root_files)
  testthat::expect_true(is.function(env$brunnermunzel_mc_test))
  testthat::expect_true(is.environment(env$.brunnermunzel_mc_state))
  testthat::expect_true(is.function(
    env$.brunnermunzel_mc_state$bm_observed_cpp
  ))
  testthat::expect_true(is.function(
    env$.brunnermunzel_mc_state$bm_mc_null_cpp
  ))
  testthat::expect_true(is.function(
    env$.brunnermunzel_mc_state$bm_mc_extreme_cpp
  ))
  testthat::expect_true(is.function(
    env$.brunnermunzel_mc_state$bm_n_allocations_cpp
  ))
  testthat::expect_true(is.function(
    env$.brunnermunzel_mc_state$bm_observed_matrix_cpp
  ))
  testthat::expect_true(is.function(
    env$.brunnermunzel_mc_state$bm_exact_extreme_cpp
  ))
  testthat::expect_true(is.function(env$brunnermunzel_pvalues))
})

testthat::test_that("the C++ path is resolved from the R file, not cwd", {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tempdir())
  temporary_wd <- getwd()
  env <- new.env(parent = globalenv())

  source(.implementation_file, local = env)

  testthat::expect_identical(getwd(), temporary_wd)
  testthat::expect_true(is.function(env$brunnermunzel_mc_test))
  testthat::expect_s3_class(
    env$brunnermunzel_mc_test(1:3, 4:6, B = 9L, seed = 1L),
    "htest"
  )
})

testthat::test_that("return value matches the principal CRAN structure", {
  x <- c(1, 2, 2)
  y <- c(2, 3, 4)
  result <- brunnermunzel_mc_test(x, y, B = 199L, seed = 7L)

  testthat::expect_s3_class(result, "htest")
  testthat::expect_identical(
    names(result),
    c("method", "data.name", "p.value", "estimate")
  )
  testthat::expect_identical(result$method, "permuted Brunner-Munzel Test")
  testthat::expect_identical(result$data.name, "x and y")
  testthat::expect_identical(
    names(result$estimate),
    "P(X<Y)+.5*P(X=Y)"
  )
  testthat::expect_type(result$p.value, "double")
  testthat::expect_silent(capture.output(print(result)))
})

testthat::test_that("effect estimates and names match brunnermunzel 2.0", {
  testthat::skip_if_not_installed("brunnermunzel")
  data_sets <- list(
    list(x = c(1, 2, 2), y = c(2, 3, 4)),
    list(x = c(1, 1, 3, 4), y = c(1, 2, 2, 5)),
    list(x = c(0, 0, 1), y = c(0, 1, 1))
  )

  for (data in data_sets) {
    for (estimate_type in c("original", "difference")) {
      expected <- brunnermunzel::brunnermunzel.permutation.test(
        data$x,
        data$y,
        force = TRUE,
        est = estimate_type
      )$estimate
      actual <- brunnermunzel_mc_test(
        data$x,
        data$y,
        est = estimate_type,
        B = 19L,
        seed = 1L
      )$estimate
      testthat::expect_equal(actual, expected, tolerance = 1e-14)
      testthat::expect_identical(names(actual), names(expected))
    }
  }
})

testthat::test_that("observed studentized statistics match readable R reference", {
  data_sets <- list(
    list(x = c(1, 2, 3), y = c(2, 4, 6)),
    list(x = c(1, 1, 3, 4), y = c(1, 2, 2, 5)),
    list(x = c(0, 0, 0), y = c(1, 1, 1))
  )

  for (data in data_sets) {
    reference <- .reference_bm(data$x, data$y)
    result <- brunnermunzel_mc_test(
      data$x,
      data$y,
      B = 1L,
      seed = 9L
    )
    testthat::expect_equal(
      attr(result, "mc")$statistic,
      reference$statistic,
      tolerance = 1e-14
    )
    testthat::expect_equal(
      .bm_statistic(data$x, data$y),
      reference$statistic,
      tolerance = 1e-14
    )
    testthat::expect_equal(
      unname(result$estimate),
      reference$effect,
      tolerance = 1e-14
    )
  }
})

testthat::test_that("standard statistics independently match CRAN BM statistics", {
  testthat::skip_if_not_installed("brunnermunzel")
  data_sets <- list(
    list(
      x = c(1, 2, 3),
      y = c(2, 4, 6)
    ),
    list(
      x = c(1, 1, 3, 4),
      y = c(1, 2, 2, 5)
    ),
    list(
      x = c(1, 2, 2, 3, 5),
      y = c(2, 4, 4)
    )
  )

  # attr(result, "mc")$statistic is the conventional BM statistic used by
  # brunnermunzel.test(). The CRAN permutation implementation internally
  # multiplies it by the positive constant (nx + ny) / (nx * ny), which does
  # not alter the ordering of statistics for fixed group sizes.
  for (data in data_sets) {
    expected <- unname(
      brunnermunzel::brunnermunzel.test(data$x, data$y)$statistic
    )
    actual <- attr(
      brunnermunzel_mc_test(
        data$x,
        data$y,
        B = 1L,
        seed = 1L
      ),
      "mc"
    )$statistic
    testthat::expect_equal(actual, expected, tolerance = 1e-12)
  }
})

testthat::test_that("all small allocations match the studentized reference", {
  combined <- c(1, 1, 2, 2, 3, 4)
  allocations <- combn(seq_along(combined), 3L, simplify = FALSE)

  for (index in allocations) {
    x <- combined[index]
    y <- combined[-index]
    reference <- .reference_bm(x, y)
    result <- brunnermunzel_mc_test(x, y, B = 1L, seed = 2L)

    testthat::expect_equal(
      attr(result, "mc")$statistic,
      reference$statistic,
      tolerance = 1e-14
    )
    testthat::expect_equal(
      reference$pooled.rank,
      rank(c(x, y), ties.method = "average")
    )
    testthat::expect_equal(
      reference$within.rank.x,
      rank(x, ties.method = "average")
    )
    testthat::expect_equal(
      reference$within.rank.y,
      rank(y, ties.method = "average")
    )
  }
})

testthat::test_that("all alternatives follow exact-enumeration directions", {
  x <- c(1, 2, 2)
  y <- c(2, 3, 4)
  exact <- lapply(
    c("two.sided", "greater", "less"),
    function(alternative) .enumerate_bm(x, y, alternative)
  )
  names(exact) <- c("two.sided", "greater", "less")

  testthat::expect_equal(exact$two.sided$p.value, 0.30)
  testthat::expect_equal(exact$greater$p.value, 1.00)
  testthat::expect_equal(exact$less$p.value, 0.15)

  for (alternative in names(exact)) {
    result <- brunnermunzel_mc_test(
      x,
      y,
      alternative = alternative,
      B = 19999L,
      seed = 104L,
      method = "mc"
    )
    standard_error <- sqrt(
      exact[[alternative]]$p.value *
        (1 - exact[[alternative]]$p.value) /
        20000
    )
    testthat::expect_equal(
      result$p.value,
      exact[[alternative]]$p.value,
      tolerance = max(0.005, 5 * standard_error)
    )
  }
})

testthat::test_that("Monte Carlo values approach CRAN exact values stably", {
  testthat::skip_if_not_installed("brunnermunzel")
  x <- c(1, 1, 3, 4)
  y <- c(1, 2, 2, 5)

  for (alternative in c("two.sided", "greater", "less")) {
    exact <- brunnermunzel::brunnermunzel.permutation.test(
      x,
      y,
      alternative = alternative,
      force = TRUE
    )
    monte_carlo <- brunnermunzel_mc_test(
      x,
      y,
      alternative = alternative,
      B = 29999L,
      seed = 812L,
      method = "mc"
    )
    standard_error <- sqrt(
      exact$p.value * (1 - exact$p.value) / 30000
    )

    testthat::expect_equal(monte_carlo$estimate, exact$estimate)
    testthat::expect_equal(
      monte_carlo$p.value,
      exact$p.value,
      tolerance = max(0.005, 5 * standard_error)
    )
  }
})

testthat::test_that("plus-one correction and inclusive boundaries are used", {
  result <- brunnermunzel_mc_test(
    c(1, 2, 2),
    c(2, 3, 4),
    B = 1000L,
    seed = 31L,
    method = "mc"
  )
  details <- attr(result, "mc")
  testthat::expect_equal(
    result$p.value,
    (details$n.extreme + 1) / (details$B + 1)
  )
  testthat::expect_gt(result$p.value, 0)

  tied <- brunnermunzel_mc_test(
    rep(0, 3),
    rep(0, 3),
    B = 100L,
    seed = 31L,
    method = "mc"
  )
  testthat::expect_equal(attr(tied, "mc")$n.extreme, 100)
  testthat::expect_equal(tied$p.value, 1)
})

testthat::test_that("CRAN-scale comparisons include near-equal boundaries", {
  count_extreme <- function(null, observed, alternative) {
    .brunnermunzel_mc_state$bm_mc_extreme_cpp(
      sort(null),
      observed,
      2L,
      2L,
      alternative
    )
  }

  testthat::expect_equal(
    count_extreme(c(-1 + 5e-15, 0, 1 - 2e-14, 1 - 5e-15), 1, 1L),
    2
  )
  testthat::expect_equal(
    count_extreme(c(-1.1, -1, -1 + 5e-15, -1 + 2e-14), -1, 2L),
    3
  )
  testthat::expect_equal(
    count_extreme(c(1 - 2e-14, 1 - 5e-15, 1, 1.1), 1, 3L),
    3
  )

  testthat::expect_equal(count_extreme(c(-Inf, 0, Inf), Inf, 1L), 2)
  testthat::expect_equal(count_extreme(c(-Inf, 0, Inf), -Inf, 2L), 1)
  testthat::expect_equal(count_extreme(c(-Inf, 0, Inf), Inf, 3L), 1)
  testthat::expect_equal(count_extreme(c(-1, 0, 1), 0, 1L), 3)
})

testthat::test_that("explicit seeds share nulls for identical tie patterns", {
  .reset_null_cache()
  old_options <- options(
    brunnermunzel.mc.cache.max.bytes = 256 * 1024^2
  )
  on.exit({
    options(old_options)
    .reset_null_cache()
  }, add = TRUE)

  first <- brunnermunzel_mc_test(
    c(1, 2, 3),
    c(4, 5, 6),
    B = 499L,
    seed = 101L,
    method = "mc"
  )
  transformed <- brunnermunzel_mc_test(
    c(10, 20, 30),
    c(40, 50, 60),
    B = 499L,
    seed = 101L,
    method = "mc"
  )
  testthat::expect_equal(.brunnermunzel_mc_state$null_generation_count, 1)
  testthat::expect_equal(first$p.value, transformed$p.value)

  tied <- brunnermunzel_mc_test(
    c(1, 1, 3),
    c(2, 4, 4),
    B = 499L,
    seed = 101L,
    method = "mc"
  )
  tied_transformed <- brunnermunzel_mc_test(
    c(10, 10, 30),
    c(20, 40, 40),
    B = 499L,
    seed = 101L,
    method = "mc"
  )
  testthat::expect_equal(.brunnermunzel_mc_state$null_generation_count, 2)
  testthat::expect_equal(tied$p.value, tied_transformed$p.value)
})

testthat::test_that("different cache key components do not share nulls", {
  .reset_null_cache()
  on.exit(.reset_null_cache(), add = TRUE)

  pattern_a_x <- c(1, 1, 3)
  pattern_a_y <- c(2, 4, 4)
  pattern_b_x <- c(1, 2, 2)
  pattern_b_y <- c(3, 3, 4)

  invisible(brunnermunzel_mc_test(
    pattern_a_x, pattern_a_y, B = 99L, seed = 1L, method = "mc"
  ))
  invisible(brunnermunzel_mc_test(
    pattern_b_x, pattern_b_y, B = 99L, seed = 1L, method = "mc"
  ))
  invisible(brunnermunzel_mc_test(
    pattern_a_x, pattern_a_y, B = 100L, seed = 1L, method = "mc"
  ))
  invisible(brunnermunzel_mc_test(
    pattern_a_x, pattern_a_y, B = 99L, seed = 2L, method = "mc"
  ))
  invisible(brunnermunzel_mc_test(
    c(1, 1, 2, 3), c(4, 4), B = 99L, seed = 1L, method = "mc"
  ))

  # pattern A has block sizes 2,1,1,2 and pattern B has 1,2,2,1:
  # the ordered sequence, rather than its sorted multiset, is part of the key.
  testthat::expect_equal(.brunnermunzel_mc_state$null_generation_count, 5)
  testthat::expect_length(
    ls(.brunnermunzel_mc_state$null_cache, all.names = TRUE),
    5
  )
})

testthat::test_that("explicit-seed p-values do not depend on gene call order", {
  .reset_null_cache()
  on.exit(.reset_null_cache(), add = TRUE)
  arguments_a <- list(
    x = c(1, 2, 3),
    y = c(4, 5, 6),
    B = 999L,
    seed = 321L,
    method = "mc"
  )
  arguments_b <- list(
    x = c(1, 1, 3),
    y = c(2, 4, 4),
    B = 999L,
    seed = 321L,
    method = "mc"
  )

  a_then_b <- list(
    a = do.call(brunnermunzel_mc_test, arguments_a),
    b = do.call(brunnermunzel_mc_test, arguments_b)
  )
  .reset_null_cache()
  b_then_a <- list(
    b = do.call(brunnermunzel_mc_test, arguments_b),
    a = do.call(brunnermunzel_mc_test, arguments_a)
  )

  testthat::expect_equal(a_then_b$a$p.value, b_then_a$a$p.value)
  testthat::expect_equal(a_then_b$b$p.value, b_then_a$b$p.value)
  testthat::expect_equal(.brunnermunzel_mc_state$null_generation_count, 2)
})

testthat::test_that("seed NULL consumes RNG and bypasses the null cache", {
  .reset_null_cache()
  on.exit(.reset_null_cache(), add = TRUE)
  x <- c(1, 1, 2, 4)
  y <- c(1, 3, 3, 5)

  set.seed(123L)
  first <- brunnermunzel_mc_test(x, y, B = 499L, method = "mc")
  state_after_first <- .Random.seed
  invisible(brunnermunzel_mc_test(x, y, B = 499L, method = "mc"))
  state_after_second <- .Random.seed

  testthat::expect_equal(.brunnermunzel_mc_state$null_generation_count, 2)
  testthat::expect_length(
    ls(.brunnermunzel_mc_state$null_cache, all.names = TRUE),
    0
  )
  testthat::expect_false(identical(state_after_first, state_after_second))

  set.seed(123L)
  repeated <- brunnermunzel_mc_test(x, y, B = 499L, method = "mc")
  testthat::expect_equal(first$p.value, repeated$p.value)
})

testthat::test_that("the LRU cache respects its configured byte limit", {
  .reset_null_cache()
  B <- 199L
  old_options <- options(
    brunnermunzel.mc.cache.max.bytes = 1024^2
  )
  on.exit({
    options(old_options)
    .reset_null_cache()
  }, add = TRUE)

  invisible(brunnermunzel_mc_test(
    c(1, 2, 3), c(4, 5, 6), B = B, seed = 1L, method = "mc"
  ))
  first_key <- ls(
    .brunnermunzel_mc_state$null_cache,
    all.names = TRUE
  )
  first_entry <- get(
    first_key,
    envir = .brunnermunzel_mc_state$null_cache,
    inherits = FALSE
  )
  max_bytes <- first_entry$bytes * 1.5
  .reset_null_cache()
  options(brunnermunzel.mc.cache.max.bytes = max_bytes)

  invisible(brunnermunzel_mc_test(
    c(1, 2, 3), c(4, 5, 6), B = B, seed = 1L, method = "mc"
  ))
  invisible(brunnermunzel_mc_test(
    c(1, 1, 3), c(2, 4, 4), B = B, seed = 1L, method = "mc"
  ))
  testthat::expect_lte(.brunnermunzel_mc_state$cache_bytes, max_bytes)
  testthat::expect_length(
    ls(.brunnermunzel_mc_state$null_cache, all.names = TRUE),
    1
  )

  invisible(brunnermunzel_mc_test(
    c(1, 2, 3), c(4, 5, 6), B = B, seed = 1L, method = "mc"
  ))
  testthat::expect_equal(.brunnermunzel_mc_state$null_generation_count, 3)
})

testthat::test_that("explicit and ambient seeds are reproducible", {
  .reset_null_cache()
  x <- c(1, 1, 2, 4)
  y <- c(1, 3, 3, 5)

  first <- brunnermunzel_mc_test(x, y, B = 999L, seed = 123L, method = "mc")
  second <- brunnermunzel_mc_test(x, y, B = 999L, seed = 123L, method = "mc")
  testthat::expect_equal(first$p.value, second$p.value)
  testthat::expect_equal(first$estimate, second$estimate)

  set.seed(456L)
  ambient_first <- brunnermunzel_mc_test(x, y, B = 999L, method = "mc")
  set.seed(456L)
  ambient_second <- brunnermunzel_mc_test(x, y, B = 999L, method = "mc")
  testthat::expect_equal(ambient_first$p.value, ambient_second$p.value)

  set.seed(789L)
  seed_before <- .Random.seed
  invisible(brunnermunzel_mc_test(x, y, B = 99L, seed = 123L, method = "mc"))
  testthat::expect_identical(.Random.seed, seed_before)
})

testthat::test_that("degenerate all-identical data are neutral and finite", {
  original <- brunnermunzel_mc_test(
    rep(0, 3),
    rep(0, 3),
    est = "original",
    B = 99L,
    seed = 1L
  )
  difference <- brunnermunzel_mc_test(
    rep(0, 3),
    rep(0, 3),
    est = "difference",
    B = 99L,
    seed = 1L
  )

  testthat::expect_equal(unname(original$estimate), 0.5)
  testthat::expect_equal(unname(difference$estimate), 0)
  testthat::expect_equal(original$p.value, 1)
  testthat::expect_false(is.nan(original$p.value))
})

testthat::test_that("invalid inputs fail explicitly", {
  valid <- c(1, 2, 3)
  testthat::expect_error(brunnermunzel_mc_test(numeric(), valid))
  testthat::expect_error(brunnermunzel_mc_test(valid, numeric()))
  testthat::expect_error(brunnermunzel_mc_test(1, valid))
  testthat::expect_error(brunnermunzel_mc_test(valid, 1))
  testthat::expect_error(brunnermunzel_mc_test(c("1", "2"), valid))
  testthat::expect_error(brunnermunzel_mc_test(factor(c("1", "2")), valid))
  testthat::expect_error(brunnermunzel_mc_test(as.Date("2020-01-01") + 0:1, valid))
  testthat::expect_error(brunnermunzel_mc_test(c(1 + 1i, 2 + 0i), valid))
  testthat::expect_error(brunnermunzel_mc_test(matrix(valid), valid))
  testthat::expect_error(brunnermunzel_mc_test(c(1, NA, 2), valid))
  testthat::expect_error(brunnermunzel_mc_test(c(1, NaN, 2), valid))
  testthat::expect_error(brunnermunzel_mc_test(c(1, Inf, 2), valid))
  testthat::expect_error(brunnermunzel_mc_test(valid, valid, B = 0))
  testthat::expect_error(brunnermunzel_mc_test(valid, valid, B = -1))
  testthat::expect_error(brunnermunzel_mc_test(valid, valid, B = 1.5))
  testthat::expect_error(brunnermunzel_mc_test(
    valid, valid, alternative = "invalid"
  ))
  testthat::expect_error(brunnermunzel_mc_test(
    valid, valid, est = "invalid"
  ))
  testthat::expect_error(brunnermunzel_mc_test(
    valid, valid, seed = NA_integer_
  ))
  testthat::expect_error(brunnermunzel_mc_test(
    valid, valid, seed = 1.5
  ))
})

testthat::test_that("partial matching, force, and downstream extraction work", {
  x <- c(1, 2, 2)
  y <- c(2, 3, 4)
  result <- brunnermunzel_mc_test(
    x,
    y,
    alternative = "t",
    force = NA,
    est = "o",
    B = 999L,
    seed = 123L
  )
  out <- c(
    pvalue = result$p.value,
    effect = unname(result$estimate)
  )

  testthat::expect_identical(names(out), c("pvalue", "effect"))
  testthat::expect_type(out, "double")
  testthat::expect_true(all(is.finite(out)))
})

# --- Exact enumeration -----------------------------------------------------

.exact_cases <- list(
  list(x = c(1, 2, 2), y = c(2, 3, 4)),
  list(x = c(1, 1, 3, 4), y = c(1, 2, 2, 5)),
  list(x = c(1, 2, 3), y = c(4, 5, 6)),
  list(x = rep(0, 3), y = rep(0, 3)),
  list(x = c(1, 1, 1, 2), y = c(1, 1, 2, 2)),
  list(x = c(5, 3, 9, 1, 7), y = c(2, 8, 4, 6)),
  list(x = c(-1, -1, 0, 0, 1, 1), y = c(0, 0, 1, 1, 2, 2))
)

testthat::test_that("the exact path reproduces full enumeration", {
  for (case in .exact_cases) {
    for (alternative in c("two.sided", "greater", "less")) {
      reference <- .enumerate_bm(case$x, case$y, alternative)
      result <- brunnermunzel_mc_test(
        case$x, case$y, alternative = alternative, method = "exact"
      )
      details <- attr(result, "mc")

      testthat::expect_identical(details$method, "exact")
      testthat::expect_equal(details$n.extreme, reference$n.extreme)
      testthat::expect_equal(
        details$n.permutations,
        choose(length(case$x) + length(case$y), length(case$x))
      )
      testthat::expect_equal(result$p.value, reference$p.value)
    }
  }
})

testthat::test_that("threaded enumeration visits every allocation once", {
  # C(24, 12) is large enough that the worker split is actually taken; the
  # small cases above would collapse to a single worker and never exercise it.
  nx <- 12L
  block_sizes <- rep(1L, 24L)
  total <- choose(24, 12)
  # A two-sided key of zero is met by every allocation and a `less` key of
  # -Inf likewise, so either count landing anywhere but `total` means the
  # worker ranges overlapped or left a gap.
  probes <- c(0, -Inf, Inf, 1, 2, 5)

  for (alternative in c(1L, 3L)) {
    counted <- lapply(c(1L, 3L, 16L), function(threads) {
      old_options <- options(brunnermunzel.exact.threads = threads)
      on.exit(options(old_options), add = TRUE)
      .bm_exact_counts(block_sizes, nx, probes, alternative)
    })

    for (result in counted) {
      testthat::expect_equal(result$total, total)
      testthat::expect_equal(
        result$count[[if (alternative == 1L) 1L else 2L]],
        total
      )
      # Counts must fall as the key gets more extreme.
      finite_probe_counts <- result$count[4:6]
      testthat::expect_true(
        all(diff(finite_probe_counts) <= 0)
      )
      testthat::expect_true(all(result$count >= 0 & result$count <= total))
    }
    testthat::expect_identical(counted[[1]]$count, counted[[2]]$count)
    testthat::expect_identical(counted[[1]]$count, counted[[3]]$count)
  }
})

testthat::test_that("exact counts are unaffected by the thread count", {
  x <- c(4, 8, 15, 16, 23, 42, 7, 9)
  y <- c(3, 11, 5, 20, 1, 33, 6, 14, 2)
  reference <- .enumerate_bm(x, y, "two.sided")$p.value
  for (threads in c(1L, 2L, 5L, 16L)) {
    old_options <- options(brunnermunzel.exact.threads = threads)
    single <- brunnermunzel_mc_test(x, y, method = "exact")
    options(old_options)
    testthat::expect_equal(single$p.value, reference)
  }
})

testthat::test_that("exact p-values fall below the Monte Carlo floor", {
  # Perfect separation at n = 24: the exact two-sided p is far under the
  # 1/(B + 1) floor that a Monte Carlo run of this size could ever report.
  x <- seq_len(12)
  y <- 100 + seq_len(12)
  exact <- brunnermunzel_mc_test(x, y, method = "exact")
  monte_carlo <- brunnermunzel_mc_test(x, y, B = 9999L, seed = 5L, method = "mc")

  testthat::expect_lt(exact$p.value, 1 / 10000)
  testthat::expect_equal(monte_carlo$p.value, 1 / 10000)
  testthat::expect_gt(exact$p.value, 0)
  testthat::expect_equal(
    attr(exact, "mc")$n.permutations,
    choose(24, 12)
  )
})

testthat::test_that("auto switches on the allocation-count budget", {
  x <- c(1, 2, 3, 4)
  y <- c(5, 6, 7, 8)
  n_allocations <- choose(8, 4)

  old_options <- options(
    brunnermunzel.exact.max.allocations = n_allocations
  )
  testthat::expect_identical(
    attr(brunnermunzel_mc_test(x, y, B = 99L, seed = 1L), "mc")$method,
    "exact"
  )
  options(brunnermunzel.exact.max.allocations = n_allocations - 1)
  testthat::expect_identical(
    attr(brunnermunzel_mc_test(x, y, B = 99L, seed = 1L), "mc")$method,
    "mc"
  )
  options(old_options)

  testthat::expect_error(
    brunnermunzel_mc_test(x, y, method = "exact", B = 99L),
    NA
  )
  # The budget is validated where it is read, not where it is set.
  options(brunnermunzel.exact.max.allocations = -1)
  testthat::expect_error(brunnermunzel_mc_test(x, y))
  options(brunnermunzel.exact.max.allocations = "many")
  testthat::expect_error(brunnermunzel_mc_test(x, y))
  options(old_options)
})

testthat::test_that("the exact path leaves the RNG and the null cache alone", {
  .reset_null_cache()
  on.exit(.reset_null_cache(), add = TRUE)
  x <- c(1, 1, 2, 4)
  y <- c(1, 3, 3, 5)

  set.seed(404L)
  seed_before <- .Random.seed
  first <- brunnermunzel_mc_test(x, y, method = "exact")
  testthat::expect_identical(.Random.seed, seed_before)
  testthat::expect_equal(.brunnermunzel_mc_state$null_generation_count, 0)
  testthat::expect_length(
    ls(.brunnermunzel_mc_state$null_cache, all.names = TRUE),
    0
  )

  # A seed is accepted for call compatibility but cannot change the answer.
  second <- brunnermunzel_mc_test(x, y, seed = 7L, method = "exact")
  testthat::expect_equal(first$p.value, second$p.value)
  testthat::expect_true(is.na(attr(first, "mc")$B))
})

testthat::test_that("allocation counts are exact and saturate to infinity", {
  count <- .brunnermunzel_mc_state$bm_n_allocations_cpp
  testthat::expect_equal(count(6L, 3L), 20)
  testthat::expect_equal(count(28L, 12L), choose(28, 12))
  testthat::expect_equal(count(36L, 9L), choose(36, 9))
  testthat::expect_equal(count(10L, 0L), 1)
  testthat::expect_equal(count(10L, 10L), 1)
  testthat::expect_true(is.infinite(count(200L, 100L)))
  testthat::expect_error(count(5L, 6L))
})

# --- Matrix interface ------------------------------------------------------

testthat::test_that("batched p-values match the one-at-a-time calls", {
  set.seed(2024L)
  nx <- 4L
  ny <- 5L
  rows <- list(
    c(c(1, 2, 3, 4), c(5, 6, 7, 8, 9)),
    c(c(9, 8, 7, 6), c(5, 4, 3, 2, 1)),
    c(c(1, 1, 2, 2), c(1, 2, 2, 3, 3)),
    c(rep(4, 4), rep(4, 5)),
    round(rnorm(nx + ny), 3),
    round(rnorm(nx + ny), 3)
  )
  data <- do.call(rbind, rows)
  rownames(data) <- paste0("gene", seq_len(nrow(data)))

  for (alternative in c("two.sided", "greater", "less")) {
    for (method in c("exact", "mc")) {
      batched <- brunnermunzel_pvalues(
        data, nx,
        alternative = alternative, B = 999L, seed = 11L, method = method
      )
      one_at_a_time <- vapply(
        seq_len(nrow(data)),
        function(i) {
          brunnermunzel_mc_test(
            data[i, seq_len(nx)],
            data[i, nx + seq_len(ny)],
            alternative = alternative,
            B = 999L, seed = 11L, method = method
          )$p.value
        },
        numeric(1)
      )
      testthat::expect_equal(as.numeric(batched), one_at_a_time)
      testthat::expect_identical(names(batched), rownames(data))
      testthat::expect_identical(attr(batched, "mc")$method, method)
    }
  }
})

testthat::test_that("batched exact p-values match full enumeration", {
  nx <- 4L
  data <- rbind(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    c(1, 1, 2, 2, 1, 2, 2, 3, 3),
    c(9, 9, 9, 9, 1, 1, 1, 1, 1)
  )
  for (alternative in c("two.sided", "greater", "less")) {
    batched <- brunnermunzel_pvalues(
      data, nx, alternative = alternative, method = "exact"
    )
    reference <- vapply(
      seq_len(nrow(data)),
      function(i) {
        .enumerate_bm(
          data[i, seq_len(nx)], data[i, -seq_len(nx)], alternative
        )$p.value
      },
      numeric(1)
    )
    testthat::expect_equal(as.numeric(batched), reference)
  }
})

testthat::test_that("the matrix interface validates its inputs", {
  data <- rbind(c(1, 2, 3, 4, 5, 6), c(2, 3, 4, 5, 6, 7))
  testthat::expect_error(brunnermunzel_pvalues(as.vector(data), 3L))
  testthat::expect_error(brunnermunzel_pvalues(data, 1L))
  testthat::expect_error(brunnermunzel_pvalues(data, 5L))
  testthat::expect_error(brunnermunzel_pvalues(data, 3L, B = 0))
  testthat::expect_error(brunnermunzel_pvalues(data, 3L, seed = 1.5))
  testthat::expect_error(
    brunnermunzel_pvalues(rbind(c(1, NA, 3, 4, 5, 6)), 3L)
  )
  testthat::expect_error(
    brunnermunzel_pvalues(rbind(c(1, Inf, 3, 4, 5, 6)), 3L)
  )
  testthat::expect_error(
    brunnermunzel_pvalues(matrix(letters[1:12], nrow = 2), 3L)
  )
  testthat::expect_length(
    brunnermunzel_pvalues(data[0, , drop = FALSE], 3L),
    0
  )
})

testthat::test_that("tie patterns group genes without changing answers", {
  nx <- 3L
  # Rows 1 and 2 share a tie pattern; row 3 has its own.
  data <- rbind(
    c(1, 2, 3, 4, 5, 6),
    c(10, 20, 30, 40, 50, 60),
    c(1, 1, 3, 4, 4, 6)
  )
  observed <- .brunnermunzel_mc_state$bm_observed_matrix_cpp(data, nx)
  testthat::expect_identical(observed$pattern[1], observed$pattern[2])
  testthat::expect_false(identical(observed$pattern[1], observed$pattern[3]))
  testthat::expect_equal(observed$statistic[1], observed$statistic[2])

  .reset_null_cache()
  on.exit(.reset_null_cache(), add = TRUE)
  invisible(brunnermunzel_pvalues(data, nx, B = 99L, seed = 3L, method = "mc"))
  testthat::expect_equal(.brunnermunzel_mc_state$null_generation_count, 2)
})
