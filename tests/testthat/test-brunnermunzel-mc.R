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
      seed = 104L
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
      seed = 812L
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
    seed = 31L
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
    seed = 31L
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
    seed = 101L
  )
  transformed <- brunnermunzel_mc_test(
    c(10, 20, 30),
    c(40, 50, 60),
    B = 499L,
    seed = 101L
  )
  testthat::expect_equal(.brunnermunzel_mc_state$null_generation_count, 1)
  testthat::expect_equal(first$p.value, transformed$p.value)

  tied <- brunnermunzel_mc_test(
    c(1, 1, 3),
    c(2, 4, 4),
    B = 499L,
    seed = 101L
  )
  tied_transformed <- brunnermunzel_mc_test(
    c(10, 10, 30),
    c(20, 40, 40),
    B = 499L,
    seed = 101L
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
    pattern_a_x, pattern_a_y, B = 99L, seed = 1L
  ))
  invisible(brunnermunzel_mc_test(
    pattern_b_x, pattern_b_y, B = 99L, seed = 1L
  ))
  invisible(brunnermunzel_mc_test(
    pattern_a_x, pattern_a_y, B = 100L, seed = 1L
  ))
  invisible(brunnermunzel_mc_test(
    pattern_a_x, pattern_a_y, B = 99L, seed = 2L
  ))
  invisible(brunnermunzel_mc_test(
    c(1, 1, 2, 3), c(4, 4), B = 99L, seed = 1L
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
    seed = 321L
  )
  arguments_b <- list(
    x = c(1, 1, 3),
    y = c(2, 4, 4),
    B = 999L,
    seed = 321L
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
  first <- brunnermunzel_mc_test(x, y, B = 499L)
  state_after_first <- .Random.seed
  invisible(brunnermunzel_mc_test(x, y, B = 499L))
  state_after_second <- .Random.seed

  testthat::expect_equal(.brunnermunzel_mc_state$null_generation_count, 2)
  testthat::expect_length(
    ls(.brunnermunzel_mc_state$null_cache, all.names = TRUE),
    0
  )
  testthat::expect_false(identical(state_after_first, state_after_second))

  set.seed(123L)
  repeated <- brunnermunzel_mc_test(x, y, B = 499L)
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
    c(1, 2, 3), c(4, 5, 6), B = B, seed = 1L
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
    c(1, 2, 3), c(4, 5, 6), B = B, seed = 1L
  ))
  invisible(brunnermunzel_mc_test(
    c(1, 1, 3), c(2, 4, 4), B = B, seed = 1L
  ))
  testthat::expect_lte(.brunnermunzel_mc_state$cache_bytes, max_bytes)
  testthat::expect_length(
    ls(.brunnermunzel_mc_state$null_cache, all.names = TRUE),
    1
  )

  invisible(brunnermunzel_mc_test(
    c(1, 2, 3), c(4, 5, 6), B = B, seed = 1L
  ))
  testthat::expect_equal(.brunnermunzel_mc_state$null_generation_count, 3)
})

testthat::test_that("explicit and ambient seeds are reproducible", {
  .reset_null_cache()
  x <- c(1, 1, 2, 4)
  y <- c(1, 3, 3, 5)

  first <- brunnermunzel_mc_test(x, y, B = 999L, seed = 123L)
  second <- brunnermunzel_mc_test(x, y, B = 999L, seed = 123L)
  testthat::expect_equal(first$p.value, second$p.value)
  testthat::expect_equal(first$estimate, second$estimate)

  set.seed(456L)
  ambient_first <- brunnermunzel_mc_test(x, y, B = 999L)
  set.seed(456L)
  ambient_second <- brunnermunzel_mc_test(x, y, B = 999L)
  testthat::expect_equal(ambient_first$p.value, ambient_second$p.value)

  set.seed(789L)
  seed_before <- .Random.seed
  invisible(brunnermunzel_mc_test(x, y, B = 99L, seed = 123L))
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
