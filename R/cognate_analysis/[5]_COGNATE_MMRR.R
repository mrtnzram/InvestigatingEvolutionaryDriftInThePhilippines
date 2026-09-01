# =============================================================================
# [5] Cognate Analysis — Multiple Matrix Regression with Randomization (MMRR)
#
# Isolation-by-distance test for cognate-class structure. Regresses pairwise
# cognate-coding similarity (cosine, over the NEXUS matrix's binary characters)
# jointly on terrain-penalized migration distance and patristic phylogenetic
# distance, with permutation p-values, so geography and shared ancestry —
# collinear here — are assessed together rather than by a single-predictor
# Mantel test. Mirrors [5]_PHONEME_MMRR.R / [5]_GRAMMAR_MMRR.R's joint
# two-predictor design (cognate has a real phylogeny, unlike
# [5]_GENETIC_MMRR.R's single-predictor case).
#
# Cognate has no [1]_*_cosine_similarity.R step (unlike phoneme/grammar), so
# the similarity matrix is built fresh here, from the NEXUS cognate-coding
# matrix [0]_COGNATE_LOANS.R already extracted, rather than read back from an
# earlier script.
#
# Inputs:  data/cognate/initial_datasets/COGNATE_NEXUS_matrix.csv,
#          data/cognate/network_distance/COGNATE_dist_matrix.csv       (written by [3])
#          data/cognate/initial_datasets/COGNATE_phylo_dist_matrix.csv (written by [0]_Phylogenetic_Tree.R)
# Outputs: data/cognate/MMRR/COGNATE_sim_matrix.csv,
#          data/cognate/MMRR/COGNATE_mmrr_results.csv        (joint two-predictor model)
#          data/cognate/MMRR/COGNATE_mmrr_single_results.csv (single-predictor models)
#          figures/cognate/mmrr/cognate_mmrr_single_*.png  (3 single-predictor)
#          figures/cognate/mmrr/cognate_mmrr_pairplot.png,
#          figures/cognate/mmrr/cognate_mmrr_partial_regression.png
# =============================================================================

library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)
library(conflicted)

# Project-wide conflict preferences (kept identical across the pipeline).
conflicts_prefer(purrr::map)
conflicts_prefer(stats::sd)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(tidyr::extract)

# ---- 0. Load data -----------------------------------------------------------
nexus_mat <- read.csv(here("data", "cognate", "initial_datasets", "COGNATE_NEXUS_matrix.csv"),
                      row.names = "glottocode", check.names = FALSE) |>
  as.matrix()

# ---- 1. Similarity matrix (cosine), built fresh from the NEXUS matrix -------
# No IDF weighting: unlike phoneme/grammar's per-column IDF, there's no natural
# frequency table for cognate-coding characters, matching the "native scaling,
# no manual weighting" choice already made for this same matrix in
# cognate_feems.ipynb. Plain cosine similarity, vectorized rather than
# [1]_PHONEME_cosine_similarity.R's pairwise loop (that loop exists there to
# apply IDF weights per pair; with no weighting, matrix algebra is equivalent
# and avoids an O(n^2) loop over ~7800 columns).
calculate_cosine_similarity <- function(X, epsilon = 1e-9) {
  norms <- sqrt(rowSums(X^2))
  (X %*% t(X)) / (outer(norms, norms) + epsilon)
}
COGNATE_sim_matrix <- calculate_cosine_similarity(nexus_mat)

# ---- 2. Pairwise terrain-penalized distance matrix (X_geo) ------------------
# Built once by [3]_COGNATE_network_distance.R and read back here instead of
# rebuilt, matching phoneme/grammar's convention.
COGNATE_dist_matrix <- read.csv(here("data", "cognate", "network_distance", "COGNATE_dist_matrix.csv"),
                                row.names = 1, check.names = FALSE) |>
  as.matrix()

# ---- 7. Phylogenetic distance matrix (X_phylo) + alignment ------------------
# Written by [0]_Phylogenetic_Tree.R "for parity with the phoneme/grammar
# pipelines" but never consumed until now.
COGNATE_phylo_dist_matrix <- read.csv(
  here("data", "cognate", "initial_datasets", "COGNATE_phylo_dist_matrix.csv"),
  row.names = 1, check.names = FALSE
) |>
  as.matrix()

# Align all three matrices to a common glottocode set and ordering before
# vectorizing; the intersection is defensive against upstream tree-set changes.
common <- Reduce(intersect, list(
  rownames(COGNATE_sim_matrix),
  rownames(COGNATE_dist_matrix),
  rownames(COGNATE_phylo_dist_matrix)
))

dropped <- setdiff(
  union(rownames(COGNATE_sim_matrix), rownames(COGNATE_phylo_dist_matrix)),
  common
)
if (length(dropped) > 0) {
  message("MMRR alignment dropped ", length(dropped),
          " language(s) not shared across all matrices: ",
          paste(dropped, collapse = ", "))
}

Y_ling  <- COGNATE_sim_matrix[common, common]
X_geo   <- COGNATE_dist_matrix[common, common]
X_phylo <- COGNATE_phylo_dist_matrix[common, common]

# X diagonals are zeroed for the guards below; Y_ling's diagonal is left at its
# natural self-similarity of ~1 (unfold() below only reads the strict lower
# triangle, so the diagonal value never enters the fit either way).
diag(X_geo) <- 0
diag(X_phylo) <- 0

# Guards: MMRR assumes symmetric matrices sharing one label order; X's are also
# zero-diagonal, Y_ling is bounded cosine similarity.
stopifnot(
  "Y/X_geo labels differ"    = identical(dimnames(Y_ling), dimnames(X_geo)),
  "Y/X_phylo labels differ"  = identical(dimnames(Y_ling), dimnames(X_phylo)),
  "Y not square"             = nrow(Y_ling) == ncol(Y_ling),
  "Y_ling not symmetric"     = isSymmetric(unname(Y_ling)),
  "X_geo not symmetric"      = isSymmetric(unname(X_geo)),
  "X_phylo not symmetric"    = isSymmetric(unname(X_phylo)),
  "Y_ling out of [-1, 1]"    = all(Y_ling >= -1 & Y_ling <= 1),
  "X_geo diagonal nonzero"   = all(diag(X_geo)    == 0),
  "X_phylo diagonal nonzero" = all(diag(X_phylo)  == 0)
)
message("MMRR input: ", nrow(Y_ling), " languages, ",
        choose(nrow(Y_ling), 2), " pairwise comparisons.")

# ---- 8. Matrix heatmaps (quick visual QC) -----------------------------------
melt_matrix <- function(m) {
  as_tibble(m, rownames = "Var1") |>
    pivot_longer(-Var1, names_to = "Var2", values_to = "value")
}

heatmap_p <- function(m, title) {
  ggplot(melt_matrix(m), aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "yellow", high = "red") +
    labs(title = title, x = "", y = "") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_fixed()
}

print(
  heatmap_p(X_geo,   "Migration distance") +
  heatmap_p(Y_ling,  "Cognate-coding similarity") +
  heatmap_p(X_phylo, "Phylogenetic distance")
)

# ---- 9. MMRR (Wang 2013) ----------------------------------------------------
# Wang Lab's canonical MMRR (Dryad doi:10.5061/dryad.kt71r/1, packaged as
# algatr::mmrr), reproduced here so the script has no external dependency —
# unfold()/MMRR() copied verbatim from [5]_PHONEME_MMRR.R, same convention
# [5]_GENETIC_MMRR.R already follows. Reference: Wang, I.J. (2013) "Examining
# the full effects of landscape heterogeneity on spatial genetic variation."
# Evolution 67:3403.
#
# Pairwise cells are not independent, so significance comes from relabeling Y's
# rows and columns with one permutation vector — preserving symmetry and the zero
# diagonal — and refitting with the predictors held fixed.

# unfold(): lower-triangle entries of a matrix as a vector; scale = TRUE
# standardizes them (Wang's default), so coefficients are standardized betas.
unfold <- function(X, scale = TRUE) {
  x <- vector()
  for (i in 2:nrow(X)) x <- c(x, X[i, 1:(i - 1)])
  if (scale) x <- scale(x, center = TRUE, scale = TRUE)
  x
}

MMRR <- function(Y, X, nperm = 9999, scale = TRUE) {
  nrowsY <- nrow(Y)
  y <- unfold(Y, scale = scale)
  if (is.null(names(X))) names(X) <- paste0("X", seq_along(X))
  Xmats <- sapply(X, unfold, scale = scale)
  colnames(Xmats) <- names(X)

  fit  <- lm(y ~ Xmats)
  summ <- summary(fit)
  coeffs    <- fit$coefficients
  r.squared <- summ$r.squared
  tstat     <- summ$coefficients[, "t value"]
  Fstat     <- summ$fstatistic[1]

  # Permutation null: relabel languages (rows + cols of Y jointly), refit.
  tprob <- rep(1, length(tstat))
  Fprob <- 1
  for (i in seq_len(nperm)) {
    rand  <- sample(seq_len(nrowsY))
    Yperm <- Y[rand, rand]
    yperm <- unfold(Yperm, scale = scale)
    summp <- summary(lm(yperm ~ Xmats))
    Fprob <- Fprob + as.numeric(summp$fstatistic[1] >= Fstat)
    tprob <- tprob + as.numeric(abs(summp$coefficients[, "t value"]) >= abs(tstat))
  }

  tp <- tprob / (nperm + 1)
  Fp <- Fprob / (nperm + 1)
  names(coeffs) <- names(tstat) <- names(tp) <- c("Intercept", names(X))

  list(coefficients = coeffs, tstatistic = tstat, tpvalue = tp,
       Fstatistic = Fstat, Fpvalue = Fp, r.squared = r.squared)
}

dir.create(here("figures", "cognate", "mmrr"), recursive = TRUE, showWarnings = FALSE)

# ---- 9a. Single-predictor MMRRs (run BEFORE the joint model) ----------------
# Same MMRR() above, just handed a one-element predictor list — the permutation
# scheme is identical, so with one predictor this reduces to a permutation-tested
# simple regression on the vectorized lower triangles, i.e. the Mantel test this
# design replaces. Because unfold() standardizes, the single-predictor beta IS
# the Pearson r and R^2 = beta^2.
#
# These are reported alongside the joint model on purpose. Geography and
# phylogeny are collinear, so each can be strongly significant alone while
# NEITHER clears the threshold once the other is in the model — the joint fit
# tests only each predictor's UNIQUE contribution. Reporting the joint model by
# itself would understate the result; reporting only these would overstate how
# separable the two effects are. The third fit (geo ~ phylo) quantifies exactly
# that overlap.
single_mmrr <- function(Ymat, Xmat, y_lab, x_lab, title, file, nperm = 9999) {
  set.seed(1)
  fit <- MMRR(Ymat, list(x = Xmat), nperm = nperm)
  beta <- unname(fit$coefficients["x"])
  pval <- unname(fit$tpvalue["x"])
  r2   <- fit$r.squared

  # Permutation p has a floor of 1/(nperm+1); show it as "<" rather than an
  # exact value it cannot resolve.
  p_txt <- if (pval <= 1 / (nperm + 1)) {
    sprintf("p < %.4f", 1 / (nperm + 1))
  } else {
    sprintf("p = %.4f", pval)
  }

  # Axes are plotted in RAW units (km, patristic, cosine) so they read cleanly;
  # standardizing only rescales the axes, leaving the fit identical.
  df <- tibble(x = unfold(Xmat, scale = FALSE), y = unfold(Ymat, scale = FALSE))

  p <- ggplot(df, aes(x, y)) +
    geom_point(alpha = 0.25, size = 0.7, colour = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, colour = "firebrick", linewidth = 0.9) +
    labs(
      title    = title,
      subtitle = sprintf("beta = %+.3f (standardized)   %s   R² = %.4f",
                         beta, p_txt, r2),
      x = x_lab, y = y_lab
    ) +
    theme_bw() +
    theme(plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 10, colour = "grey25"))

  print(p)
  ggsave(here("figures", "cognate", "mmrr", file), p,
         width = 6.5, height = 4.5, units = "in", dpi = 300)

  tibble(model = title, beta = beta, p_value = pval, r_squared = r2)
}

COGNATE_mmrr_single <- bind_rows(
  single_mmrr(Y_ling, X_geo,
              y_lab = "Cognate-coding similarity (cosine)",
              x_lab = "Relative migration distance (km)",
              title = "Cognate similarity vs Geographic distance",
              file  = "cognate_mmrr_single_similarity_vs_geography.png"),
  single_mmrr(Y_ling, X_phylo,
              y_lab = "Cognate-coding similarity (cosine)",
              x_lab = "Patristic phylogenetic distance",
              title = "Cognate similarity vs Phylogenetic distance",
              file  = "cognate_mmrr_single_similarity_vs_phylogeny.png"),
  # Collinearity check: the shared structure that flattens both joint betas.
  single_mmrr(X_geo, X_phylo,
              y_lab = "Relative migration distance (km)",
              x_lab = "Patristic phylogenetic distance",
              title = "Geographic distance vs Phylogenetic distance",
              file  = "cognate_mmrr_single_geography_vs_phylogeny.png")
)
print(COGNATE_mmrr_single)
write.csv(COGNATE_mmrr_single,
          file = here("data", "cognate", "MMRR", "COGNATE_mmrr_single_results.csv"), row.names = FALSE)


# ---- 9b. Joint two-predictor MMRR -------------------------------------------
set.seed(1)  # reproducible permutation p-values
mmrr_fit <- MMRR(Y_ling, list(geo = X_geo, phylo = X_phylo), nperm = 9999)
print(mmrr_fit)

# ---- 10. Tidy results tibble ------------------------------------------------
COGNATE_mmrr_results <- tibble(
  dataset    = "COGNATE",
  beta_geo   = unname(mmrr_fit$coefficients["geo"]),
  p_geo      = unname(mmrr_fit$tpvalue["geo"]),
  beta_phylo = unname(mmrr_fit$coefficients["phylo"]),
  p_phylo    = unname(mmrr_fit$tpvalue["phylo"]),
  r_squared  = mmrr_fit$r.squared,
  p_model    = unname(mmrr_fit$Fpvalue)
)
print(COGNATE_mmrr_results)
write.csv(COGNATE_mmrr_results,
          file = here("data", "cognate", "MMRR", "COGNATE_mmrr_results.csv"), row.names = FALSE)

# ---- 11. Visualization ------------------------------------------------------
# Both figures are built on the standardized unfolded lower triangles, so they
# read on the same scale as the fitted (standardized) MMRR coefficients.
# (figures/cognate/mmrr already created in §9a, ahead of the first ggsave())

mmrr_df <- tibble(
  ling  = as.numeric(unfold(Y_ling)),
  geo   = as.numeric(unfold(X_geo)),
  phylo = as.numeric(unfold(X_phylo))
)

# (a) Pairplot as a single 3x3 facet_grid so every cell shares its column (x) and
# row (y) scale. Lower triangle: scatter + linear fit; upper: Pearson r;
# diagonal: marginal density. Built directly to avoid a GGally dependency.
pair_labs <- c(ling = "Similarity", geo = "Geo distance", phylo = "Phylo distance")
pair_vars <- names(pair_labs)
pair_idx  <- setNames(seq_along(pair_vars), pair_vars)
as_pair_factor <- function(v) factor(pair_labs[v], levels = pair_labs)

# Cell text is centred on each variable's midpoint, not 0: the standardized
# inputs are mildly right-skewed, so the panel middle is not exactly 0.
pair_mid <- colMeans(sapply(mmrr_df[pair_vars], range))

# Lower triangle (row index > col index): scatter of col-var (x) vs row-var (y).
scatter_df <- expand_grid(row = pair_vars, col = pair_vars) |>
  filter(pair_idx[row] > pair_idx[col]) |>
  mutate(dat = map2(row, col, ~ tibble(x = mmrr_df[[.y]], y = mmrr_df[[.x]]))) |>
  unnest(dat) |>
  mutate(row = as_pair_factor(row), col = as_pair_factor(col))

# Upper triangle (row index < col index): Pearson r, centred in the panel.
cor_df <- expand_grid(row = pair_vars, col = pair_vars) |>
  filter(pair_idx[row] < pair_idx[col]) |>
  mutate(
    lab = sprintf("r = %.2f", map2_dbl(row, col, ~ cor(mmrr_df[[.x]], mmrr_df[[.y]]))),
    x   = pair_mid[col], y = pair_mid[row],
    row = as_pair_factor(row), col = as_pair_factor(col)
  )

# Anchor every panel's free scale to the col-var (x) and row-var (y) ranges;
# without it the text/density-only cells collapse to a degenerate axis.
anchor_df <- expand_grid(row = pair_vars, col = pair_vars) |>
  mutate(dat = map2(row, col, ~ tibble(x = range(mmrr_df[[.y]]), y = range(mmrr_df[[.x]])))) |>
  unnest(dat) |>
  mutate(row = as_pair_factor(row), col = as_pair_factor(col))

# Diagonal densities are rescaled into the lower 85% of the panel because the
# row's y-scale is the variable's value range, not a 0-1 density axis; from/to
# clip the curve so its tails stay inside the anchored panel.
density_df <- map_dfr(pair_vars, function(v) {
  rng <- range(mmrr_df[[v]])
  d   <- density(mmrr_df[[v]], from = rng[1], to = rng[2])
  tibble(v = v, x = d$x,
         ymin = rng[1],
         ymax = rng[1] + (d$y / max(d$y)) * 0.85 * diff(rng))
}) |>
  mutate(row = as_pair_factor(v), col = as_pair_factor(v))

pairplot <- ggplot() +
  geom_blank(data = anchor_df, aes(x, y)) +
  geom_ribbon(data = density_df, aes(x = x, ymin = ymin, ymax = ymax),
              fill = "grey80", color = "grey40", linewidth = 0.3) +
  geom_point(data = scatter_df, aes(x, y),
             alpha = 0.4, size = 0.6, color = "steelblue") +
  geom_smooth(data = scatter_df, aes(x, y), method = "lm", se = FALSE,
              color = "firebrick", linewidth = 0.7) +
  geom_text(data = cor_df, aes(x, y, label = lab), size = 4.5) +
  facet_grid(row ~ col, scales = "free") +
  theme_bw() +
  theme(
    axis.title       = element_blank(),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 9),
    panel.spacing    = unit(4, "pt")
  ) +
  labs(title = "MMRR inputs")
print(pairplot)
ggsave(here("figures", "cognate", "mmrr", "cognate_mmrr_pairplot.png"),
       pairplot, width = 6.5, height = 6, units = "in", dpi = 300)

# (b) Added-variable (partial regression) plots: each panel shows the
# isolation-by-distance relationship holding the other predictor constant, and
# its fitted slope equals the corresponding MMRR beta.
#   geo | phylo : resid(Y ~ X_phylo) vs resid(X_geo ~ X_phylo), slope = beta_geo
#   phylo | geo : resid(Y ~ X_geo)   vs resid(X_phylo ~ X_geo), slope = beta_phylo
av_plot <- function(y, x, other, xlab, ylab, beta, pval) {
  ry <- resid(lm(y ~ other))
  rx <- resid(lm(x ~ other))
  ggplot(tibble(rx = rx, ry = ry), aes(rx, ry)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "firebrick", linewidth = 1) +
    theme_bw() +
    labs(
      x = xlab, y = ylab,
      subtitle = sprintf("beta = %.3f,  p = %.4f", beta, pval)
    )
}

p_geo_av <- av_plot(
  mmrr_df$ling, mmrr_df$geo, mmrr_df$phylo,
  xlab = "Geo distance | phylo (residuals)",
  ylab = "Similarity | phylo (residuals)",
  beta = COGNATE_mmrr_results$beta_geo, pval = COGNATE_mmrr_results$p_geo
)

p_phylo_av <- av_plot(
  mmrr_df$ling, mmrr_df$phylo, mmrr_df$geo,
  xlab = "Phylo distance | geo (residuals)",
  ylab = "Similarity | geo (residuals)",
  beta = COGNATE_mmrr_results$beta_phylo, pval = COGNATE_mmrr_results$p_phylo
)

partial_plot <- (p_geo_av + p_phylo_av) +
  plot_annotation(
    title = "MMRR partial regression: cognate-coding similarity vs. geography + phylogeny",
    subtitle = sprintf("Model R^2 = %.3f,  permutation p = %.4f  (9,999 permutations)",
                       COGNATE_mmrr_results$r_squared, COGNATE_mmrr_results$p_model)
  )
print(partial_plot)
ggsave(here("figures", "cognate", "mmrr", "cognate_mmrr_partial_regression.png"),
       partial_plot, width = 10, height = 4.5, units = "in", dpi = 300)

# ---- 12. Persist matrices ---------------------------------------------------
write.csv(COGNATE_sim_matrix, file = here("data", "cognate", "MMRR", "COGNATE_sim_matrix.csv"), row.names = TRUE)
