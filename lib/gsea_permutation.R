# Phenotype-permutation GSEA.
#
# The permutation is over sample labels, never over genes: gene permutation
# treats genes as exchangeable and is anti-conservative whenever they are
# correlated, which across a transcriptome they are. Enrichment scores come
# from fgsea's kernel and are normalized to NES per pathway and sign
# (Subramanian et al. 2005).
#
# Inference here is the Westfall-Young family-wise error rate
# (gsea_westfall_young), not an FDR: its null is the largest |NES| anywhere in
# the collection per shuffle, so inter-pathway correlation is carried exactly.
# The Subramanian FDR (gsea_subramanian_fdr) is provided only as a descriptive
# contrast -- it divides by the mean of the pooled null, which under
# inter-pathway correlation tracks whichever way the whole collection drifts in
# one realization, and was measured on the study data to call 32/50 Hallmark
# sets in a unit with no gene-level signal. Callers evaluate one collection at
# a time so that each collection's max-|NES| null (and its NES normalization)
# is built from that collection alone.

if (!requireNamespace("fgsea", quietly = TRUE)) {
  stop(
    "Package 'fgsea' is required to source lib/gsea_permutation.R.",
    call. = FALSE
  )
}

# fgsea's batch kernel scores every pathway in a single call, is far faster
# than looping the exported one, and agrees with it exactly. Fall back to the
# exported function if a future fgsea stops shipping it.
.gsea_score_batch <- local({
  namespace <- asNamespace("fgsea")
  if (exists("calcGseaStatBatchCpp", envir = namespace, inherits = FALSE)) {
    batch <- get("calcGseaStatBatchCpp", envir = namespace)
    function(sorted, index) {
      as.numeric(batch(
        stats = sorted, selectedGenes = index, geneRanks = seq_along(sorted)
      ))
    }
  } else {
    function(sorted, index) {
      vapply(index, function(k) fgsea::calcGseaStat(sorted, k), numeric(1))
    }
  }
})

.gsea_ranking <- function(metric) {
  ordering <- order(metric, decreasing = TRUE)
  ranks <- integer(length(metric))
  ranks[ordering] <- seq_along(metric)
  list(ordering = ordering, sorted = metric[ordering], ranks = ranks)
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

# fgsea applies the exponent by transforming the statistics, which is also the
# only way to reach the batch kernel. gsea_param = 0 is the unweighted
# Kolmogorov-Smirnov form: every hit counts once and the magnitudes drop out.
.gsea_weight <- function(sorted, gsea_param) {
  if (gsea_param == 1) {
    return(sorted)
  }
  sign(sorted) * abs(sorted)^gsea_param
}

gsea_scores <- function(metric, index, gsea_param) {
  ranking <- .gsea_ranking(metric)
  scores <- .gsea_score_batch(
    .gsea_weight(ranking$sorted, gsea_param),
    lapply(index, function(k) sort(ranking$ranks[k]))
  )
  names(scores) <- names(index)
  scores
}

gsea_leading_edge <- function(metric, index, gene_ids, gsea_param) {
  ranking <- .gsea_ranking(metric)
  weighted <- .gsea_weight(ranking$sorted, gsea_param)
  lapply(index, function(k) {
    edge <- fgsea::calcGseaStat(
      weighted, sort(ranking$ranks[k]), returnLeadingEdge = TRUE
    )$leadingEdge
    gene_ids[ranking$ordering[edge]]
  })
}

# Enforce the step-up shape the Broad procedure gives its q-values: a more
# extreme NES may not carry a larger FDR than a less extreme one.
.gsea_monotone_fdr <- function(nes, fdr) {
  out <- fdr
  positive <- which(is.finite(nes) & nes >= 0)
  negative <- which(is.finite(nes) & nes < 0)
  if (length(positive)) {
    o <- positive[order(nes[positive], decreasing = TRUE)]
    out[o] <- cummin(fdr[o])
  }
  if (length(negative)) {
    o <- negative[order(nes[negative])]
    out[o] <- cummin(fdr[o])
  }
  out
}

# `null` is pathways x permutations, from the same label shuffles that produced
# `observed`. Scores are normalized within pathway and sign, so a pathway's
# observed NES and its own null NES stay on one scale.
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

# Single-step Westfall-Young: the null is the largest |NES| anywhere in the
# collection per shuffle, so the correlation among pathways is carried exactly
# rather than assumed away. Controls the family-wise error rate.
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

# Subramanian et al. (2005). Kept for comparability with standard GSEA output
# and deliberately NOT used for inference here: it compares the observed NES
# distribution with the mean of the pooled null, which under inter-pathway
# correlation is driven by the direction the whole collection happens to drift
# in one realization. Verified on this data to call 74% of Hallmark in a unit
# with no signal.
gsea_subramanian_fdr <- function(nes, nes_null) {
  pooled <- nes_null[is.finite(nes_null)]
  finite <- nes[is.finite(nes)]
  # Sorted references so each tail fraction is a binary search rather than a
  # scan of the whole pooled null, which has n_pathways x n_permutations entries.
  pooled_positive <- sort(pooled[pooled >= 0])
  pooled_negative <- sort(pooled[pooled < 0])
  observed_positive <- sort(finite[finite >= 0])
  observed_negative <- sort(finite[finite < 0])
  at_least <- function(x, reference) {
    if (length(reference) == 0L) {
      return(rep(0, length(x)))
    }
    (length(reference) -
      findInterval(x, reference, left.open = TRUE)) / length(reference)
  }
  at_most <- function(x, reference) {
    if (length(reference) == 0L) {
      return(rep(0, length(x)))
    }
    findInterval(x, reference) / length(reference)
  }

  fdr <- rep(NA_real_, length(nes))
  positive <- which(is.finite(nes) & nes >= 0)
  negative <- which(is.finite(nes) & nes < 0)
  if (length(positive)) {
    denominator <- at_least(nes[positive], observed_positive)
    fdr[positive] <- ifelse(
      denominator == 0, 1,
      pmin(1, at_least(nes[positive], pooled_positive) / denominator)
    )
  }
  if (length(negative)) {
    denominator <- at_most(nes[negative], observed_negative)
    fdr[negative] <- ifelse(
      denominator == 0, 1,
      pmin(1, at_most(nes[negative], pooled_negative) / denominator)
    )
  }

  .gsea_monotone_fdr(nes, fdr)
}

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
