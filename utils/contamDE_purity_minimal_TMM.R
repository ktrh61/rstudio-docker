# ==============================================================================
# contamDE_purity_minimal_TMM.R
# Minimal contamDE purity with TMM normalization (original contamDE approach)
# Dependencies: edgeR, limma, qvalue
# ==============================================================================

# --- voom with TMM normalization ---
limma_voom_tmm <- function(counts) {
  counts <- as.matrix(counts)
  N <- ncol(counts) / 2
  stopifnot(N == floor(ncol(counts) / 2))
  
  # DGEList with TMM normalization
  d <- edgeR::DGEList(counts = counts)
  d <- edgeR::calcNormFactors(d, method = "TMM")
  
  # Design matrix (Normal/Tumor)
  condition <- factor(c(rep("Normal", N), rep("Tumor", N)))
  design <- stats::model.matrix(~ 0 + condition)
  colnames(design) <- levels(condition)
  
  # voom transformation
  v <- limma::voom(d, design, normalize.method = "none")
  
  # Linear model and contrasts
  fit <- limma::lmFit(v, design)
  ctr <- limma::makeContrasts(contrasts = "Tumor - Normal", levels = design)
  fit2 <- limma::contrasts.fit(fit, ctr)
  
  # eBayes with original contamDE settings
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  
  # Effective library sizes
  size <- d$samples$lib.size * d$samples$norm.factors
  size <- size / mean(size)
  
  list(
    p.limma = fit2$p.value[, 1],
    log2FC.limma = fit2$coefficients[, 1],
    size = size
  )
}

# --- Main contamDE purity estimation function ---
contamDE_purity_minimal_TMM <- function(counts,
                                        subtype = NULL,
                                        covariate = NULL,
                                        is.contaminated = TRUE,
                                        prior.count = 1,      # Original contamDE default
                                        verbose = TRUE) {
  
  if (verbose) cat("Starting tumor purity estimation (TMM version)...\n")
  
  counts <- as.matrix(counts)
  N <- ncol(counts) / 2
  G <- nrow(counts)
  stopifnot(N == floor(ncol(counts) / 2))
  
  # Initial ranking with limma+voom (TMM)
  d <- limma_voom_tmm(counts)
  p.limma <- d$p.limma
  log2FC.limma <- d$log2FC.limma
  size <- d$size
  
  # Design matrix (original contamDE structure)
  if (is.null(subtype) || length(unique(subtype)) == 1) {
    subtype <- rep(1, N)
    if (is.null(covariate)) {
      design <- matrix(subtype, ncol = 1)
    } else {
      covariate <- matrix(covariate, ncol = N)
      design <- stats::model.matrix(~ covariate)
    }
  } else {
    subtype <- factor(subtype)
    if (is.null(covariate)) {
      design <- stats::model.matrix(~ 0 + subtype)
    } else {
      covariate0 <- matrix(covariate, nrow = N)
      design <- stats::model.matrix(~ 0 + subtype + covariate0)
    }
  }
  
  # Normalize counts by TMM size factors
  count.norm <- t(t(counts) / size)
  count.normal <- count.norm[, (1:N)] + prior.count
  count.tumor <- count.norm[, -(1:N)] + prior.count
  
  # Calculate log2 ratios
  y <- log2(count.tumor) - log2(count.normal)
  
  # Initialize purity
  w_hat <- rep(1, N)
  up <- down <- integer(0)
  
  if (is.contaminated) {
    # Multiple testing correction with qvalue
    p_adj <- tryCatch(
      qvalue::qvalue(p.limma, pi0.method = "bootstrap")$qvalues,
      error = function(e) {
        if (verbose) cat("  qvalue failed, using BH adjustment\n")
        stats::p.adjust(p.limma, method = "BH")
      }
    )
    
    # Limit to top 1000 genes (original contamDE approach)
    n_sig <- sum(p_adj < 0.1, na.rm = TRUE)
    if (n_sig > 1000L) {
      idx_sig <- which(is.finite(p_adj) & p_adj < 0.1)
      idx_sorted <- idx_sig[order(p_adj[idx_sig])]
      keep <- idx_sorted[seq_len(1000L)]
      p_adj[-keep] <- 1
      if (verbose) {
        cat(sprintf("  Limited to top 1000 genes (from %d with q < 0.1)\n", n_sig))
      }
    }
    
    # Define informative genes
    up <- which(p_adj < 0.1 & log2FC.limma > log2(1.5))
    down <- which(p_adj < 0.1 & log2FC.limma < -log2(1.5))
    
    if (verbose) {
      cat(sprintf("  Informative genes - Up: %d, Down: %d\n", 
                  length(up), length(down)))
    }
    
    # Calculate purity estimates
    if (length(up) > 0) {
      y_up <- y[up, , drop = FALSE]
    } else {
      y_up <- matrix(0, 0, N)
    }
    
    if (length(down) > 0) {
      y_down <- y[down, , drop = FALSE]
    } else {
      y_down <- matrix(0, 0, N)
    }
    
    sumup <- colSums(y_up)
    sumdown <- colSums(y_down)
    sum.max <- max(sumup - sumdown)
    
    if (!is.finite(sum.max) || sum.max <= 0) {
      if (verbose) {
        cat("  Warning: No informative genes found. Setting purity to 1.0.\n")
      }
      w_hat <- rep(1.0, N)
    } else {
      w_hat <- (sumup - sumdown) / sum.max
      w_hat <- pmax(0, pmin(1, w_hat))
    }
    
    if (verbose) {
      cat(sprintf("  Purity summary: mean=%.3f, sd=%.3f\n", 
                  mean(w_hat), sd(w_hat)))
    }
  } else {
    if (verbose) {
      cat("  Contamination correction disabled. Setting purity to 1.0.\n")
    }
  }
  
  # Return results
  list(
    proportion = w_hat,
    size = size,
    design = design,
    y = y,
    n_pairs = N,
    n_genes = G,
    informative_genes = list(up = up, down = down),
    p_init = p.limma,
    log2FC_init = log2FC.limma,
    normalization_method = "TMM",
    estimation_date = Sys.time()
  )
}

cat("contamDE TMM version loaded successfully!\n")