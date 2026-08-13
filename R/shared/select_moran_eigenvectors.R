# =============================================================================
# Shared — Moran's-I eigenvector selection (genetics' PVR::PVR(method="moran")
# analogue)
#
# Hand-rolled equivalent of PVR::PVR(method = "moran") for a plain eigenvector
# matrix: PVR() requires a `phylo` object, so it can't run on genetics'
# population-structure PCA eigenvectors. Adds eigenvectors one at a time and
# stops once residual Moran's I (tested against Euclidean distance in
# eigenvector space) is no longer significant. Used by both
# [4]_GENETIC_PVR.R and [4.1]_GENETIC_PVR_sPCA.R so they select the same subset.
# =============================================================================

library(spdep)

#' Select the smallest leading eigenvector subset that leaves no significant
#' residual spatial (here: population-structure) autocorrelation.
#'
#' @param y      Response vector.
#' @param x      Environmental covariate vector (e.g. geodist_H1_span, km),
#'                same length as y; always included alongside the candidate
#'                eigenvectors, mirroring PVR(..., envVar = x, ...).
#' @param E      Matrix of candidate eigenvectors, columns ordered by
#'                importance (e.g. PC1, PC2, ...); nrow(E) == length(y).
#' @param sig.t  Significance threshold. Default 0.05, matching PVR's default.
#' @return E[, 1:k] for the smallest k (possibly 0) whose residual Moran's I
#'          is non-significant; the full E with a warning if none qualify.
select_moran_eigenvectors <- function(y, x, E, sig.t = 0.05) {
  stopifnot(
    "nrow(E) must match length(y)" = nrow(E) == length(y),
    "length(x) must match length(y)" = length(x) == length(y)
  )

  # Fixed connectivity from the full eigenvector set, tested against every candidate k.
  D <- as.matrix(dist(E))
  W <- ifelse(D == 0, 0, 1 / D^2)
  diag(W) <- 0
  listw <- mat2listw(W, style = "W", zero.policy = TRUE)

  moran_p <- function(k) {
    m <- if (k == 0) lm(y ~ x) else lm(y ~ E[, seq_len(k), drop = FALSE] + x)
    lm.morantest(m, listw, zero.policy = TRUE)$p.value
  }

  for (k in 0:ncol(E)) {
    p <- moran_p(k)
    if (p > sig.t) {
      message(sprintf(
        "select_moran_eigenvectors: stopped at k = %d eigenvector%s (Moran's I p = %.4f > %.3f).",
        k, if (k == 1) "" else "s", p, sig.t))
      return(if (k == 0) E[, 0, drop = FALSE] else E[, seq_len(k), drop = FALSE])
    }
  }

  warning("select_moran_eigenvectors: residual autocorrelation stayed significant through all ",
          ncol(E), " eigenvectors; returning the full set.")
  E
}
