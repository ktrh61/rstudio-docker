.repo_root <- normalizePath(
  file.path(testthat::test_path(), "..", ".."),
  mustWork = TRUE
)
source(file.path(.repo_root, "lib", "gsea_permutation.R"))

# --- Tie-free agreement with standard GSEA (fgsea kernel, gseaParam = 1) ----
# On tie-free input every block is a singleton, so the block ES must equal the
# standard weighted GSEA statistic exactly (reorg plan v2 D3).

test_that("block ES equals fgsea calcGseaStat on tie-free input", {
  skip_if_not_installed("fgsea")
  set.seed(42)
  n <- 2000L
  for (rep in 1:5) {
    scores <- sample(rnorm(n)) # continuous: tie-free almost surely
    stopifnot(!anyDuplicated(scores))
    index <- lapply(c(5L, 20L, 100L, 400L), function(k) sample.int(n, k))
    names(index) <- paste0("set", seq_along(index))
    block_es <- gsea_block_scores(scores, index)

    ordering <- order(scores, decreasing = TRUE)
    sorted <- scores[ordering]
    position <- integer(n)
    position[ordering] <- seq_len(n)
    fgsea_es <- vapply(index, function(k) {
      fgsea::calcGseaStat(sorted, sort.int(position[k]), gseaParam = 1)
    }, numeric(1))
    expect_equal(block_es, fgsea_es, tolerance = 1e-12)
  }
})

test_that("block leading edge equals fgsea leading edge on tie-free input", {
  skip_if_not_installed("fgsea")
  set.seed(7)
  n <- 500L
  scores <- sample(rnorm(n))
  gene_ids <- sprintf("g%04d", seq_len(n))
  index <- lapply(c(10L, 50L, 200L), function(k) sample.int(n, k))
  names(index) <- paste0("set", seq_along(index))
  block_edge <- gsea_block_leading_edge(scores, index, gene_ids)

  ordering <- order(scores, decreasing = TRUE)
  sorted <- scores[ordering]
  position <- integer(n)
  position[ordering] <- seq_len(n)
  for (name in names(index)) {
    ref <- fgsea::calcGseaStat(
      sorted, sort.int(position[index[[name]]]),
      gseaParam = 1, returnLeadingEdge = TRUE
    )$leadingEdge
    expect_setequal(block_edge[[name]], gene_ids[ordering[ref]])
  }
})

# --- Small reference vectors (hand-computed) --------------------------------

test_that("block ES matches the hand-computed tied reference", {
  # Sorted scores 3, 2, 2, 1, -1, -2; blocks [1], [2,3], [4], [5], [6].
  # Hits: the score-3 gene and one of the tied score-2 genes.
  # N_R = 5, miss step 1/4.
  #   end of block 1: 3/5 - 0 = 0.6
  #   end of block 2: 5/5 - 1/4 = 0.75  <- ES (tie-break-free extreme)
  # (A within-block evaluation would give 1.0 or 0.75 depending on the
  # arbitrary tie-break; the block end is the invariant value.)
  scores <- c(3, 2, 2, 1, -1, -2)
  expect_equal(
    unname(gsea_block_scores(scores, list(s = c(1L, 2L)))), 0.75
  )
  # Same set, but the other tied gene: identical by construction.
  expect_equal(
    unname(gsea_block_scores(scores, list(s = c(1L, 3L)))), 0.75
  )
  # Same data with the genes stored in a different order: sets are gene
  # indices, not sorted positions, so the ES must not change.
  shuffled <- c(-1, 2, 3, 2, 1, -2) # gene 3 has score 3, genes 2 and 4 tie
  expect_equal(
    unname(gsea_block_scores(shuffled, list(s = c(3L, 2L)))), 0.75
  )
})

test_that("block ES matches the hand-computed negative reference", {
  # Tie-free scores 5..1, hits at the bottom two. Running difference reaches
  # -1 just before the first hit, which is the max deviation.
  scores <- c(5, 4, 3, 2, 1)
  expect_equal(
    unname(gsea_block_scores(scores, list(s = c(4L, 5L)))), -1
  )
  edge <- gsea_block_leading_edge(
    scores, list(s = c(4L, 5L)), sprintf("g%d", 1:5)
  )
  expect_setequal(edge$s, c("g4", "g5"))
})

test_that("a set with zero total hit weight scores 0 with an empty edge", {
  scores <- c(1, 0, 0, -1)
  expect_equal(unname(gsea_block_scores(scores, list(s = c(2L, 3L)))), 0)
  edge <- gsea_block_leading_edge(
    scores, list(s = c(2L, 3L)), sprintf("g%d", 1:4)
  )
  expect_length(edge$s, 0L)
})

# --- Normal scores ----------------------------------------------------------

test_that("normal scores are tie-shared, monotone, and bound infinities", {
  metric <- c(-Inf, -1, 0, 1, Inf, Inf)
  scores <- gsea_normal_scores(metric)
  expect_true(all(is.finite(scores)))
  expect_equal(scores[5], scores[6])
  expect_true(all(diff(scores[1:5]) > 0))
  expect_equal(scores, qnorm((c(1, 2, 3, 4, 5.5, 5.5) - 0.5) / 6))
})

# --- Tail-ratio q (hand-computed) -------------------------------------------

test_that("tail-ratio q matches the hand-computed reference", {
  nes <- c(2, 1, -1)
  nes_null <- rbind(c(1.5, 0.5), c(0.5, -0.5), c(-1.5, 0.5))
  # Pooled positive null {1.5, 0.5, 0.5, 0.5}, negative {-1.5, -0.5};
  # observed positive {2, 1}, negative {-1}.
  #   q(2)  = [(0+1)/(4+1)] / (1/2) = 0.4
  #   q(1)  = [(1+1)/(4+1)] / (2/2) = 0.4
  #   q(-1) = [(1+1)/(2+1)] / (1/1) = 2/3
  expect_equal(
    gsea_tail_ratio_q(nes, nes_null), c(0.4, 0.4, 2 / 3),
    tolerance = 1e-12
  )
})

test_that("tail-ratio q is plus-one smoothed and capped at 1", {
  nes <- c(3, 0.1)
  nes_null <- rbind(c(0.5, 0.4), c(0.3, 0.2))
  q <- gsea_tail_ratio_q(nes, nes_null)
  expect_true(all(q > 0)) # never exactly zero
  expect_true(all(q <= 1))
  # NES = 3 above every null value: q = (0+1)/(4+1) / (1/2) = 0.4
  expect_equal(q[1], 0.4, tolerance = 1e-12)
})
