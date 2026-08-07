# ==============================================================================
# MUREN core
# Statistical model follows the original MUREN implementation. Pairwise work is
# evaluated directly in deterministic chunks to reduce parallel overhead.
# Requires lib/norm_muren_helpers.R to be sourced first.
# ==============================================================================

if (!requireNamespace("foreach", quietly = TRUE)) {
  stop("Package 'foreach' is required.")
}
if (!requireNamespace("doSNOW", quietly = TRUE)) {
  stop("Package 'doSNOW' is required.")
}
if (!requireNamespace("iterators", quietly = TRUE)) {
  stop("Package 'iterators' is required.")
}
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required.")
}
if (!requireNamespace("parallel", quietly = TRUE)) {
  stop("Package 'parallel' is required.")
}

`%dopar%` <- foreach::`%dopar%`

muren_norm <- function(reads,
                       refs = "saturated",
                       pairwise_method = "lts",
                       single_param = TRUE,
                       res_return = "counts",
                       filter_gene = TRUE,
                       trim = 10,
                       maxiter = 70,
                       workers,
                       ...) {
  if (missing(workers)) {
    stop("'workers' must be explicitly specified.")
  }

  is_flag <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }

  if (!is_flag(single_param)) {
    stop("'single_param' must be TRUE or FALSE.")
  }
  if (!is_flag(filter_gene)) {
    stop("'filter_gene' must be TRUE or FALSE.")
  }
  if (!is.character(pairwise_method) || length(pairwise_method) != 1L ||
      is.na(pairwise_method) || !pairwise_method %in% c("lts", "mode")) {
    stop("'pairwise_method' must be either 'lts' or 'mode'.")
  }
  if (!single_param && pairwise_method != "lts") {
    stop("The double-parameter form is defined only for pairwise_method='lts'.")
  }
  if (!is.character(res_return) || length(res_return) != 1L ||
      is.na(res_return) ||
      !res_return %in% c("counts", "log_counts", "scaling_coeff")) {
    stop("'res_return' must be 'counts', 'log_counts', or 'scaling_coeff'.")
  }
  if (!single_param && res_return == "scaling_coeff") {
    stop("res_return='scaling_coeff' is valid only when single_param=TRUE.")
  }
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) || trim < 0) {
    stop("'trim' must be one finite non-negative number.")
  }
  if (!is.numeric(maxiter) || length(maxiter) != 1L || !is.finite(maxiter) ||
      maxiter <= 0 || !is.wholenumber(maxiter)) {
    stop("'maxiter' must be a positive integer.")
  }
  if (!is.numeric(workers) || length(workers) != 1L || !is.finite(workers) ||
      workers <= 0 || !is.wholenumber(workers)) {
    stop("'workers' must be a fixed positive integer.")
  }
  if (!is.data.frame(reads) && !is.matrix(reads)) {
    stop("'reads' must be a data.frame or matrix.")
  }
  if (nrow(reads) <= 1L) {
    stop("At least two genes are required.")
  }

  workers <- as.integer(workers)
  maxiter <- as.integer(maxiter)
  pair_args <- list(...)
  raw_reads <- reads

  i_sample <- rep(TRUE, ncol(raw_reads))
  if (is.data.frame(raw_reads)) {
    i_sample <- vapply(raw_reads, is.numeric, logical(1))
  }
  if (!any(i_sample)) {
    stop("No numeric sample columns were found in 'reads'.")
  }

  reads <- as.matrix(raw_reads[, i_sample, drop = FALSE])
  if (any(!is.finite(reads)) || any(reads < 0)) {
    stop("Numeric sample columns must contain finite non-negative counts.")
  }
  if (ncol(reads) <= 1L) {
    stop("At least two samples are required.")
  }

  i_gene <- filter_gene_l(reads, trim)
  kept_gene <- if (filter_gene) i_gene else rep(TRUE, nrow(reads))
  reads <- reads[kept_gene, , drop = FALSE]
  if (nrow(reads) <= 1L) {
    stop("At least two genes must remain after filtering.")
  }

  log_raw_reads_mx <- lg(reads)
  n_exp <- ncol(reads)
  n_gene <- nrow(reads)

  # Resolve the reference set once. The target sample is never removed from a
  # reference set that contains it, matching the original MUREN protocol.
  if (is.character(refs)) {
    if (length(refs) == 1L && identical(refs, "saturated")) {
      ref_idx <- seq_len(n_exp)
    } else {
      if (length(refs) == 0L || anyNA(refs) || anyDuplicated(refs) ||
          is.null(colnames(reads)) || !all(refs %in% colnames(reads))) {
        stop("Bad specification of references !")
      }
      ref_idx <- which(colnames(reads) %in% refs)
    }
  } else if (is.numeric(refs)) {
    if (length(refs) == 0L || any(!is.finite(refs)) ||
        !all(is.wholenumber(refs)) || anyDuplicated(refs) ||
        any(refs < 1L) || any(refs > n_exp)) {
      stop("Bad specification of references !")
    }
    ref_idx <- as.integer(refs)
  } else {
    stop("Bad specification of references !")
  }

  used_refs <- sort(unique(ref_idx))
  unused_refs <- setdiff(seq_len(n_exp), used_refs)

  pairs <- cbind(
    target = rep(seq_len(n_exp), each = length(ref_idx)),
    reference = rep(ref_idx, times = n_exp)
  )
  locations <- pairs[, "reference"] +
    (pairs[, "target"] - 1L) * n_exp

  # For single-parameter LTS, d(i,j) = -d(j,i) and d(i,i) = 0.
  # Compute each unordered non-diagonal pair once, then reconstruct the full
  # directed matrix required by median polish.
  optimized_sp_lts <- single_param && pairwise_method == "lts"
  if (optimized_sp_lts) {
    diagonal_pair <- pairs[, "target"] == pairs[, "reference"]
    lower_sample <- pmin(pairs[, "target"], pairs[, "reference"])
    upper_sample <- pmax(pairs[, "target"], pairs[, "reference"])
    pair_key <- paste(lower_sample, upper_sample, sep = ":")
    calculation_keys <- unique(pair_key[!diagonal_pair])
    calculation_rows <- match(calculation_keys, pair_key)
    work_pairs <- cbind(
      target = lower_sample[calculation_rows],
      reference = upper_sample[calculation_rows]
    )
    result_index <- match(pair_key, calculation_keys)
    result_sign <- ifelse(pairs[, "target"] == lower_sample, 1, -1)
  } else {
    work_pairs <- pairs
  }

  required_helpers <- c(
    "lg", "ep", "TOL",
    if (single_param) "polish_coeff" else "polish_one_gene",
    if (!single_param) "reg_dp" else if (pairwise_method == "lts") "reg_sp" else "mode_sp"
  )
  missing_helpers <- required_helpers[
    !vapply(required_helpers, exists, logical(1), envir = .GlobalEnv, inherits = TRUE)
  ]
  if (length(missing_helpers) > 0L) {
    stop(
      "Required MUREN helpers not found: ",
      paste(missing_helpers, collapse = ", ")
    )
  }

  cl <- parallel::makeCluster(workers, type = "PSOCK")
  doSNOW::registerDoSNOW(cl)
  parallel::clusterSetRNGStream(cl, 12345L)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  parallel::clusterExport(cl, varlist = required_helpers, envir = .GlobalEnv)

  n_work <- nrow(work_pairs)
  n_parts <- min(workers, n_work)
  split_idx <- split(
    seq_len(n_work),
    rep(seq_len(n_parts), each = ceiling(n_work / n_parts), length.out = n_work)
  )

  if (single_param) {
    pair_results <- foreach::foreach(
      idx = iterators::iter(split_idx),
      .packages = if (pairwise_method == "lts") "MASS" else character(0),
      .export = c(
        if (pairwise_method == "lts") "reg_sp" else "mode_sp",
        "log_raw_reads_mx", "work_pairs", "pair_args", "pairwise_method"
      ),
      .combine = "c",
      .inorder = TRUE
    ) %dopar% {
      out <- numeric(length(idx))
      for (ii in seq_along(idx)) {
        p <- idx[ii]
        i <- work_pairs[p, "target"]
        j <- work_pairs[p, "reference"]
        call_args <- c(
          list(s_k = log_raw_reads_mx[, i], s_r = log_raw_reads_mx[, j]),
          pair_args
        )
        if (pairwise_method == "lts") {
          out[ii] <- do.call(reg_sp, call_args)
        } else {
          out[ii] <- do.call(mode_sp, call_args)
        }
      }
      out
    }
    pair_results <- as.numeric(pair_results)

    if (optimized_sp_lts) {
      res_pairwise <- numeric(nrow(pairs))
      res_pairwise[!diagonal_pair] <-
        pair_results[result_index[!diagonal_pair]] *
        result_sign[!diagonal_pair]
    } else {
      res_pairwise <- pair_results
    }
  } else {
    res_chunks <- foreach::foreach(
      idx = iterators::iter(split_idx),
      .packages = "MASS",
      .export = c(
        "reg_dp", "log_raw_reads_mx", "work_pairs", "pair_args", "n_gene"
      ),
      .combine = "cbind",
      .inorder = TRUE
    ) %dopar% {
      out <- matrix(NA_real_, nrow = n_gene, ncol = length(idx))
      for (ii in seq_along(idx)) {
        p <- idx[ii]
        i <- work_pairs[p, "target"]
        j <- work_pairs[p, "reference"]
        out[, ii] <- do.call(
          reg_dp,
          c(
            list(s_k = log_raw_reads_mx[, i], s_r = log_raw_reads_mx[, j]),
            pair_args
          )
        )
      }
      out
    }
    res_pairwise <- as.matrix(res_chunks)
  }

  if (single_param) {
    if (length(res_pairwise) != length(locations)) {
      stop("Pairwise result length mismatch (single_param).")
    }

    coef_sp <- polish_coeff(
      fitted_n = res_pairwise,
      n_exp = n_exp,
      locations = locations,
      unused_refs = unused_refs,
      maxiter = maxiter
    )
    coef_sp <- 2^as.vector(coef_sp)
    names(coef_sp) <- colnames(reads)

    if (res_return == "scaling_coeff") {
      return(1 / coef_sp)
    }

    polished_mx <- sweep(
      as.matrix(raw_reads[, i_sample, drop = FALSE]),
      2,
      coef_sp,
      `*`
    )
    if (res_return == "log_counts") {
      polished_mx <- lg(polished_mx)
    }
  } else {
    if (!is.matrix(res_pairwise) || nrow(res_pairwise) != n_gene) {
      stop("Pairwise result shape mismatch (double_param).")
    }

    polished_mx <- foreach::foreach(
      n = iterators::iter(seq_len(n_gene)),
      .combine = rbind,
      .export = c("polish_one_gene"),
      .inorder = TRUE
    ) %dopar% {
      polish_one_gene(
        fitted_n = res_pairwise[n, ],
        n_exp = n_exp,
        locations = locations,
        unused_refs = unused_refs,
        maxiter = maxiter
      )
    }

    polished_mx[log_raw_reads_mx < TOL | polished_mx < 0] <- 0
    if (res_return == "counts") {
      polished_mx <- ep(polished_mx)
    }
  }

  if (!is.null(rownames(raw_reads))) {
    rownames(polished_mx) <- if (single_param) {
      rownames(raw_reads)
    } else {
      rownames(raw_reads)[kept_gene]
    }
  }
  colnames(polished_mx) <- colnames(reads)

  if (is.data.frame(raw_reads)) {
    output_rows <- if (single_param) rep(TRUE, nrow(raw_reads)) else kept_gene
    result <- raw_reads[output_rows, , drop = FALSE]
    result[, i_sample] <- polished_mx
    return(result)
  }

  as.matrix(polished_mx)
}
