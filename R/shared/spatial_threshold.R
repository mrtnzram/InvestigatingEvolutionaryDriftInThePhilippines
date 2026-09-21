# =============================================================================
# Shared — geography-only distance threshold for spatial weights
#
# The [4.1] / [4.2] / [4.3] sPCA/MEM scripts connect two languages in the
# spatial weights matrix only when their network distance is <= tau. tau is
# fixed from geography alone: the longest edge of the minimum spanning tree
# over the distance matrix, i.e. the smallest cutoff at which the neighbour
# graph is connected and every language has at least one neighbour.
#
# This replaces the earlier data-driven search (20 quantile candidates, picking
# the one that minimized |Moran's I| or maximized the sPCA eigenvalue). That
# search chose tau from the response being tested, which the permutation tests
# did not account for, and at small N it jumped between distant candidates
# (Mindanao: 369 -> 899 km after a minor change to the response).
# =============================================================================

#' Longest minimum-spanning-tree edge of a distance matrix (Prim's algorithm).
#'
#' Zero off-diagonal distances (languages sharing coordinates) are valid MST
#' edges, so duplicate-coordinate groups are handled without special-casing.
#'
#' @param D Symmetric, finite distance matrix (km).
#' @return The MST's longest edge length (km), used as the threshold tau.
mst_threshold <- function(D) {
  D <- unname(as.matrix(D))
  stopifnot(
    "Distance matrix must be square." = nrow(D) == ncol(D),
    "Distance matrix must be finite (no NA/Inf)." = all(is.finite(D)),
    "Distance matrix must be symmetric." = isSymmetric(D)
  )
  n <- nrow(D)
  if (n < 2) return(NA_real_)

  in_tree <- c(TRUE, rep(FALSE, n - 1))
  reach   <- D[1, ]   # cheapest edge from the tree to each node
  longest <- 0
  for (step in seq_len(n - 1)) {
    reach[in_tree] <- Inf
    j <- which.min(reach)
    longest    <- max(longest, reach[j])
    in_tree[j] <- TRUE
    reach      <- pmin(reach, D[j, ])
  }
  longest
}

#' Inverse-squared-distance weights, truncated at tau (not row-standardized;
#' mat2listw(style = "W") or spca(matWeight = ) standardize downstream).
#'
#' @param D   Distance matrix (km).
#' @param tau Cutoff (km); pairs farther apart get weight 0.
#' @return Weight matrix with zero diagonal. Zero-distance pairs get 0, not Inf.
threshold_weights <- function(D, tau) {
  W <- 1 / D^2
  W[!is.finite(W)] <- 0
  W[D > tau] <- 0
  diag(W) <- 0
  W
}
