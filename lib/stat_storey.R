# stat_storey.R
# The protocol-wide multiple-testing correction (reorg plan v2 D1): Storey's
# q-value with the plug-in pi0 estimate at a fixed lambda,
#   pi0_hat = min(1, mean(p > lambda) / (1 - lambda)),   lambda = 0.5
#   q       = min(1, pi0_hat * BH(p))
# Monotone because BH is monotone. One definition, so the DEGES screen
# (310 via lib/norm_deges.R) and the gene-level inference (410) cannot drift
# apart: both headers claim they apply the same correction, and this file is
# what makes that claim structural rather than a hand-synced comment.
# At the set level the protocol holds pi0 at its conservative bound 1, where
# the q-value reduces to plain BH (reorg plan v2 appendix B.1) -- callers use
# p.adjust directly there.

storey_pi0 <- function(p, lambda = 0.5) {
  min(1, mean(p > lambda) / (1 - lambda))
}

storey_q <- function(p, lambda = 0.5) {
  pi0_hat <- storey_pi0(p, lambda)
  list(pi0 = pi0_hat, q = pmin(1, pi0_hat * p.adjust(p, method = "BH")))
}
