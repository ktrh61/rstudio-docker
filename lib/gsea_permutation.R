# Phenotype-permutation GSEA, spec B (reorg plan v2 D2/D3).
#
# The permutation is over sample labels, never over genes: gene permutation
# treats genes as exchangeable and is anti-conservative whenever they are
# correlated, which across a transcriptome they are.
#
# Ranking metric (gsea_normal_scores): tie-averaged normal scores of the
# signed Brunner-Munzel statistic, qnorm((rank - 0.5) / G) with average ranks
# on ties. The studentized statistic can reach +/-Inf (its denominator
# vanishes for some allocations); the normal score bounds every gene's weight
# while preserving the ordering, which is the only information the protocol's
# rank-based commitments treat as real.
#
# Enrichment score (gsea_block_scores): weighted running sum (Subramanian
# et al. 2005, gseaParam = 1) evaluated only at tie-block boundaries. Inside a
# block of tied scores no ordering exists, so any within-block running-sum
# value would depend on an arbitrary tie-break; block ends are invariant. On
# tie-free input this is exactly the standard GSEA statistic (verified against
# fgsea in tests/testthat/test-gsea-block-es.R). Boundary conventions, fixed
# here: the positive extreme is scanned at the end of every block containing a
# hit, the negative extreme just before every such block; when the two
# extremes tie in magnitude the positive one is returned; a set whose hit
# weights sum to zero scores 0. Leading edges are reported at block
# granularity (whole tie-blocks, never a fraction of one).
#
# Inference: per-set permutation p (gsea_pathway_pvalues -- each set against
# its own null row, sign-conditional) with BH within the collection, applied
# by the caller. BH is the Storey q-value machinery with pi0 held at its
# conservative bound 1; the protocol-wide q framework is unchanged, only the
# pi0 plug-in is not estimated at this level (24-3,900 dependent p-values
# moving on one collective-drift mode leave the estimator variance-dominated,
# and its errors point anti-conservative exactly in the drift realizations).
# This choice is the one the held-out null calibration passed: the pooled
# tail-ratio FDR (Subramanian 2005) and a per-realization restandardized
# variant were both measured miscalibrated under this dependence structure
# (global-null P(>=1 discovery) 0.14 and 0.22 against nominal 0.10; per-set
# p + BH measured 0.045). Measurements and the mechanism are recorded in
# diagnostics/output/ and reorg plan v2 appendix B.
# Callers evaluate one collection at a time so that each collection's NES
# normalization is built from that collection alone.

# --- Ranking metric ---------------------------------------------------------
# Tie-averaged normal scores: monotone in the metric, tied metrics share one
# score, and +/-Inf map to finite extremes.
gsea_normal_scores <- function(metric) {
  r <- rank(metric, ties.method = "average")
  qnorm((r - 0.5) / length(metric))
}

# Pathways as positions into `gene_ids`, dropping members the data does not
# carry and sets left outside the size window.
gsea_pathway_index <- function(pathways, gene_ids, min_size, max_size) {
  index <- lapply(pathways, function(members) {
    positions <- match(unique(members), gene_ids)
    positions[!is.na(positions)]
  })
  index[lengths(index) >= min_size & lengths(index) <= max_size]
}

# --- Tie-block invariant enrichment score -----------------------------------
# Per-metric preprocessing shared by every set: descending order, each gene's
# position in that order, and the tie-block structure of the sorted scores.
gsea_block_ranking <- function(scores) {
  ordering <- order(scores, decreasing = TRUE)
  sorted <- scores[ordering]
  n <- length(sorted)
  new_block <- c(TRUE, sorted[-1L] != sorted[-n])
  block_id <- cumsum(new_block)
  block_start <- which(new_block)
  block_end <- c(block_start[-1L] - 1L, n)
  position <- integer(n)
  position[ordering] <- seq_len(n)
  list(
    ordering = ordering, sorted = sorted, n = n,
    position = position, block_id = block_id,
    block_start = block_start, block_end = block_end
  )
}

# Core running-sum evaluation for one set. P_hit weights each hit by |score|
# (gseaParam = 1) and P_miss steps by 1/(n - k); both are cumulated to the
# ends of the tie-blocks that contain hits (positive candidates) and to the
# positions just before those blocks (negative candidates), which are the only
# order-invariant extremes of the running difference.
.gsea_block_es <- function(ranking, set_genes, leading_edge = FALSE) {
  h <- sort.int(ranking$position[set_genes])
  k <- length(h)
  w <- abs(ranking$sorted[h])
  total_w <- sum(w)
  miss_scale <- 1 / (ranking$n - k)
  if (total_w == 0) {
    if (leading_edge) {
      return(list(es = 0, edge = integer(0)))
    }
    return(0)
  }
  bid <- ranking$block_id[h]
  group_last <- which(c(bid[-1L] != bid[-k], TRUE))
  cum_w <- cumsum(w)[group_last]
  cum_h <- group_last
  bids <- bid[group_last]
  n_groups <- length(group_last)
  v_pos <- cum_w / total_w -
    (ranking$block_end[bids] - cum_h) * miss_scale
  prev_w <- c(0, cum_w[-n_groups])
  prev_h <- c(0L, cum_h[-n_groups])
  v_neg <- prev_w / total_w -
    (ranking$block_start[bids] - 1L - prev_h) * miss_scale
  v_max <- max(v_pos)
  v_min <- min(v_neg)
  es <- if (v_max >= -v_min) v_max else v_min
  if (!leading_edge) {
    return(es)
  }
  edge <- if (es >= 0) {
    h[seq_len(group_last[which.max(v_pos)])]
  } else {
    first <- prev_h[which.min(v_neg)] + 1L
    if (first > k) integer(0) else h[first:k]
  }
  list(es = es, edge = edge)
}

gsea_block_scores <- function(scores, index) {
  ranking <- gsea_block_ranking(scores)
  vapply(index, function(k) .gsea_block_es(ranking, k), numeric(1))
}

gsea_block_leading_edge <- function(scores, index, gene_ids) {
  ranking <- gsea_block_ranking(scores)
  lapply(index, function(k) {
    gene_ids[ranking$ordering[.gsea_block_es(ranking, k, TRUE)$edge]]
  })
}

# --- NES --------------------------------------------------------------------
# `null` is pathways x permutations, from the same label shuffles that
# produced `observed`. Scores are normalized within pathway and sign, so a
# pathway's observed NES and its own null NES stay on one scale.
gsea_nes <- function(observed, null) {
  positive_mean <- apply(null, 1L, function(z) {
    z <- z[is.finite(z) & z > 0]
    if (length(z) == 0L) NA_real_ else mean(z)
  })
  negative_mean <- apply(null, 1L, function(z) {
    z <- z[is.finite(z) & z < 0]
    if (length(z) == 0L) NA_real_ else mean(abs(z))
  })
  list(
    nes = observed / ifelse(observed >= 0, positive_mean, negative_mean),
    nes_null = null / ifelse(null >= 0, positive_mean, negative_mean)
  )
}

# --- Per-set permutation p (the inference carrier; BH is applied on top) ----
# Each pathway against its own null row, conditional on the observed sign.
# The plus-one counting keeps p > 0 at permutation resolution.
gsea_pathway_pvalues <- function(nes, nes_null) {
  vapply(seq_along(nes), function(i) {
    if (!is.finite(nes[i])) {
      return(NA_real_)
    }
    z <- nes_null[i, ]
    z <- z[is.finite(z)]
    if (nes[i] >= 0) {
      (sum(z >= nes[i]) + 1) / (sum(z >= 0) + 1)
    } else {
      (sum(z <= nes[i]) + 1) / (sum(z <= 0) + 1)
    }
  }, numeric(1))
}

# --- Redundancy annotation --------------------------------------------------
# Marks a pathway redundant when its leading edge is largely contained in that
# of a better-ranked pathway, which is what makes an ontology's parent and
# child terms both surface on one signal. Only `candidate` rows are compared:
# the comparison is quadratic, and a set nobody would report needs no parent.
gsea_redundancy <- function(result, leading_edge, threshold, candidate) {
  parent <- rep(NA_character_, nrow(result))
  ordering <- which(candidate)
  if (length(ordering) == 0L) {
    return(parent)
  }
  ordering <- ordering[
    order(result$pval[ordering], -abs(result$NES[ordering]))
  ]
  kept <- integer(0)
  for (i in ordering) {
    edge <- leading_edge[[result$pathway[i]]]
    if (length(edge) == 0L) {
      next
    }
    for (j in kept) {
      other <- leading_edge[[result$pathway[j]]]
      shared <- length(intersect(edge, other))
      if (shared / length(union(edge, other)) >= threshold) {
        parent[i] <- result$pathway[j]
        break
      }
    }
    if (is.na(parent[i])) {
      kept <- c(kept, i)
    }
  }
  parent
}
