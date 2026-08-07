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
# Inference (gsea_tail_ratio_q): the Subramanian tail-ratio FDR computed
# within one collection -- the ratio of the pooled-null tail fraction to the
# observed tail fraction at each set's NES, by sign. The pooled null
# (collection sets x shuffles) is what gives the estimate its m x B
# resolution; NES standardization makes the sets exchangeable under the null,
# which matches a hypothesis that carries no set-level prediction. Fixed
# conventions: signs are handled separately against sign-restricted
# references; the null tail fraction is plus-one smoothed like a permutation
# p-value (so q > 0 always); the observed tail fraction is exact (never zero
# at a set's own NES); an empty null side yields the smoothed ceiling rather
# than NA; q is capped at 1; no monotonization is applied -- the raw ratio is
# reported.
#
# The Westfall-Young FWER (gsea_westfall_young) is retained as a sensitivity
# column only: its null is the largest |NES| anywhere in the collection per
# shuffle, so inter-pathway correlation is carried exactly, but it answers a
# sparse-signal question the protocol no longer poses at this level.
# Callers evaluate one collection at a time so that each collection's pooled
# null and NES normalization are built from that collection alone.

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

# --- Per-set permutation p (descriptive; ordering only) ---------------------
# Each pathway against its own null. The collection-wide drift that inter-
# pathway correlation produces sits in both terms, so it cancels here.
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

# --- Collection-internal tail-ratio FDR (inference) -------------------------
# Subramanian et al. (2005) within one collection; conventions fixed in the
# file header. Sorted sign-restricted references turn each tail count into a
# binary search.
gsea_tail_ratio_q <- function(nes, nes_null) {
  pooled <- nes_null[is.finite(nes_null)]
  finite <- nes[is.finite(nes)]
  pooled_positive <- sort(pooled[pooled >= 0])
  pooled_negative <- sort(pooled[pooled < 0])
  observed_positive <- sort(finite[finite >= 0])
  observed_negative <- sort(finite[finite < 0])
  count_at_least <- function(x, reference) {
    length(reference) - findInterval(x, reference, left.open = TRUE)
  }
  count_at_most <- function(x, reference) {
    findInterval(x, reference)
  }

  q <- rep(NA_real_, length(nes))
  positive <- which(is.finite(nes) & nes >= 0)
  negative <- which(is.finite(nes) & nes < 0)
  if (length(positive)) {
    null_tail <- (count_at_least(nes[positive], pooled_positive) + 1) /
      (length(pooled_positive) + 1)
    observed_tail <- count_at_least(nes[positive], observed_positive) /
      length(observed_positive)
    q[positive] <- pmin(1, null_tail / observed_tail)
  }
  if (length(negative)) {
    null_tail <- (count_at_most(nes[negative], pooled_negative) + 1) /
      (length(pooled_negative) + 1)
    observed_tail <- count_at_most(nes[negative], observed_negative) /
      length(observed_negative)
    q[negative] <- pmin(1, null_tail / observed_tail)
  }
  q
}

# --- Westfall-Young FWER (sensitivity column) -------------------------------
# Single-step: the null is the largest |NES| anywhere in the collection per
# shuffle, so the correlation among pathways is carried exactly rather than
# assumed away.
gsea_westfall_young <- function(nes, nes_null) {
  null_max <- apply(abs(nes_null), 2L, max, na.rm = TRUE)
  n_perm <- length(null_max)
  vapply(seq_along(nes), function(i) {
    if (!is.finite(nes[i])) {
      return(NA_real_)
    }
    (sum(null_max >= abs(nes[i])) + 1) / (n_perm + 1)
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
