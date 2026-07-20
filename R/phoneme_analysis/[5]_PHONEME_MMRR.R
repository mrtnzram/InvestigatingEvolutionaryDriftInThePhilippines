# =============================================================================
# [5] Phoneme Analysis — Multiple Matrix Regression with Randomization (MMRR)
#
# Control analysis for internal linguistic diffusion. Regresses pairwise
# GENERAL phonemic dissimilarity between the Philippine languages jointly on
# two pairwise predictors:
#
#     Y_ling = b0 + b1 * X_geo + b2 * X_phylo + e
#
#   Y_ling  : 1 - cosine similarity (the same general dissimilarity the old
#             Mantel test used — NOT Spanish-specific, NOT a delta measure)
#   X_geo   : terrain-penalized migration distance through the waypoint network
#   X_phylo : patristic (cophenetic) phylogenetic distance from the ABVD tree
#
# This isolates internal diffusion — geographic and/or genealogical — from the
# external (colonial-contact) signal that is this study's primary focus.
# Replacing the single-predictor Mantel test with MMRR lets geography and shared
# ancestry be assessed jointly, since geographically close languages are often
# also close relatives (see the geo/phylo collinearity flagged in [4] §9).
#
# This file keeps the original Dijkstra routing (shortest_path_trace) to build
# the full pairwise migration-distance matrix — it is intentionally NOT the
# simplified "distance to nearest node" used in [3]_PHONEME_network_distance.R,
# because MMRR needs genuine language-to-language distances.
#
# Inputs:  data/PHONEME_cosine_matrix.csv, data/RUHLENdf_PH.csv,
#          data/nodes.csv, data/edges.csv,
#          data/PHONEME_phylo_dist_matrix.csv  (written by [0]_Phylogenetic_Tree.R)
# Outputs: data/PHONEME_dist_matrix.csv, data/PHONEME_diss_matrix.csv,
#          data/PHONEME_mmrr_results.csv,
#          figures/phoneme/mmrr/phoneme_mmrr_pairplot.png,
#          figures/phoneme/mmrr/phoneme_mmrr_partial_regression.png
# =============================================================================

library(readr)
library(tidyverse)
library(dplyr)
library(geosphere)
library(sf)
library(maps)
library(sfheaders)
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
RUHLENdf <- read_csv(here("data", "RUHLENdf_PH.csv"))

ph_lang <- RUHLENdf |>
  filter(Language_type == "Philippine Language") |>
  pull(language)

# Cosine-similarity matrix written by [1]_PHONEME_cosine_similarity.R
cosine_matrix <- read.csv(here("data", "PHONEME_cosine_matrix.csv"),
                          row.names = 1, check.names = FALSE) |>
  as.matrix()

nodes <- read.csv(here("data", "nodes.csv"))
edges <- read.csv(here("data", "edges.csv"))

nodes$id   <- as.character(nodes$id)
edges$from <- as.character(edges$from)
edges$to   <- as.character(edges$to)

# Add reverse edges for bidirectional routing
reverse_edges <- edges |> rename(from = to, to = from)
edges <- bind_rows(edges, reverse_edges) |> distinct()

# ---- 1. Dissimilarity matrix (1 - cosine), Philippine languages -------------
cosine_matrix_phil <- cosine_matrix[ph_lang, ph_lang]
PHONEME_diss_matrix <- 1 - cosine_matrix_phil

# ---- 2. Land mask (Philippines + Malaysia) ----------------------------------
world_map <- map_data("world") |> filter(region %in% c("Philippines", "Malaysia"))

land_sf <- sf_polygon(obj = world_map, polygon_id = "group", x = "long", y = "lat") |>
  st_union() |>
  st_sf(geometry = _) |>
  st_set_crs(4326)

# ---- 3. Terrain-penalized network edges + weighted graph --------------------
# NOTE: edges are penalized at 44.18 (matching the original [1] route weighting),
# while the language connectors below use 4.44. This penalty mismatch is carried
# over verbatim from the monolith to preserve results — flag for review.
edge_land_penalty <- 44.18

edges <- edges |>
  rowwise() |>
  mutate(weight = {
    from_coords <- nodes |> filter(id == from)
    to_coords   <- nodes |> filter(id == to)
    if (nrow(from_coords) == 0 || nrow(to_coords) == 0) NA_real_
    else distHaversine(c(from_coords$longitude, from_coords$latitude),
                       c(to_coords$longitude, to_coords$latitude))
  }) |>
  ungroup()

edge_lines <- edges |>
  rowwise() |>
  mutate(
    geometry = list(st_linestring(matrix(c(
      nodes$longitude[nodes$id == from], nodes$latitude[nodes$id == from],
      nodes$longitude[nodes$id == to],   nodes$latitude[nodes$id == to]
    ), ncol = 2, byrow = TRUE)))
  ) |>
  ungroup() |>
  st_as_sf(crs = 4326)

edge_lines <- edge_lines |>
  rowwise() |>
  mutate(
    land_part = list(st_intersection(geometry, land_sf)),
    sea_part  = list(st_difference(geometry, land_sf)),

    land_len = as.numeric(if (!is.null(land_part) && length(land_part) > 0) st_length(land_part) else 0),
    sea_len  = as.numeric(if (!is.null(sea_part)  && length(sea_part)  > 0) st_length(sea_part)  else 0),

    weighted_cost = land_len * edge_land_penalty + sea_len,
    crosses_land  = land_len > 0
  ) |>
  ungroup()

# Adjacency list keyed by node id, weighted by terrain-penalized cost
all_ids <- unique(c(edge_lines$from, edge_lines$to))
graph <- lapply(all_ids, function(id) {
  neighbors <- edge_lines |> filter(from == id) |> select(to, weighted_cost)
  if (nrow(neighbors) == 0) tibble(to = character(), weight = numeric())
  else rename(neighbors, weight = weighted_cost)
})
names(graph) <- all_ids

# ---- 4. Routing helpers -----------------------------------------------------
# Every node in nodes.csv participates in the edge graph (no orphans), so a
# language's geographically nearest node is always a routable entry point onto
# the network. shortest_path_trace() is a hand-rolled Dijkstra kept verbatim
# from the original Mantel script so the X_geo matrix is reproducible
# bit-for-bit. For a faster equivalent, igraph::distances() on the same weighted
# graph returns identical shortest-path costs (see [3] for the igraph pattern);
# left as-is to avoid any numeric drift in the published results.
find_nearest_node <- function(coords) {
  distances <- distHaversine(matrix(c(nodes$longitude, nodes$latitude), ncol = 2),
                             coords)
  nodes$id[which.min(distances)]
}

shortest_path_trace <- function(start_id, end_id, graph) {
  visited <- setNames(rep(FALSE, length(graph)), names(graph))
  dist <- setNames(rep(Inf, length(graph)), names(graph))
  prev <- setNames(rep(NA_character_, length(graph)), names(graph))
  dist[start_id] <- 0
  queue <- data.frame(id = start_id, dist = 0)

  while (nrow(queue) > 0) {
    queue <- queue[order(queue$dist), ]
    current <- queue$id[1]
    current_dist <- queue$dist[1]
    queue <- queue[-1, ]

    if (visited[current]) next
    visited[current] <- TRUE

    neighbors <- graph[[current]]
    if (is.null(neighbors)) next

    for (i in seq_len(nrow(neighbors))) {
      neighbor <- neighbors$to[i]
      weight <- neighbors$weight[i]
      if (is.na(weight)) next
      if (dist[neighbor] > current_dist + weight) {
        dist[neighbor] <- current_dist + weight
        prev[neighbor] <- current
        queue <- rbind(queue, data.frame(id = neighbor, dist = dist[neighbor]))
      }
    }
  }

  if (!is.finite(dist[end_id])) return(list(distance = NA_real_, path = NULL))

  # Reconstruct path
  path <- end_id
  while (!is.na(prev[path[1]])) {
    path <- c(prev[path[1]], path)
  }
  return(list(distance = dist[end_id], path = path))
}

# ---- 5. Per-language connector penalties (coord -> nearest node) -------------
land_penalty <- 4.44

phil_df <- RUHLENdf |>
  filter(Language_type == "Philippine Language") |>
  mutate(
    start_coords = map2(longitude, latitude, ~ c(.x, .y)),
    nearest_node = map_chr(start_coords, find_nearest_node)
  )

connector_df <- phil_df |>
  mutate(
    connector_geom = map2(start_coords, nearest_node, ~ st_linestring(rbind(
      .x,
      c(nodes$longitude[nodes$id == .y], nodes$latitude[nodes$id == .y])
    )))
  ) |>
  mutate(connector_geom_sfc = st_sfc(connector_geom, crs = 4326)) |>
  rowwise() |>
  mutate(
    land_part = list(st_intersection(connector_geom_sfc, land_sf)),
    sea_part  = list(st_difference(connector_geom_sfc, land_sf)),

    land_len = as.numeric(if (!is.null(land_part) && length(land_part) > 0) st_length(land_part) else 0),
    sea_len  = as.numeric(if (!is.null(sea_part)  && length(sea_part)  > 0) st_length(sea_part)  else 0),

    connector_penalty = land_len * land_penalty + sea_len
  ) |>
  ungroup()

# ---- 6. Pairwise terrain-penalized distance matrix (X_geo) ------------------
phil_pairs <- expand_grid(lang1 = phil_df$language, lang2 = phil_df$language) |>
  filter(lang1 != lang2) |>
  left_join(phil_df |> select(language, node1 = nearest_node), by = c("lang1" = "language")) |>
  left_join(phil_df |> select(language, node2 = nearest_node), by = c("lang2" = "language")) |>
  left_join(connector_df |> select(language, penalty1 = connector_penalty), by = c("lang1" = "language")) |>
  left_join(connector_df |> select(language, penalty2 = connector_penalty), by = c("lang2" = "language"))

phil_pairs <- phil_pairs |>
  rowwise() |>
  mutate(
    trace = list(shortest_path_trace(node1, node2, graph)),
    tree_dist = trace$distance,
    geodist_H1_span = if (is.na(tree_dist)) NA_real_ else
      (penalty1 + tree_dist + penalty2) / 1000
  ) |>
  ungroup()

dist_matrix <- phil_pairs |>
  select(lang1, lang2, geodist_H1_span) |>
  pivot_wider(names_from = lang2, values_from = geodist_H1_span) |>
  column_to_rownames("lang1") |>
  as.matrix()

PHONEME_dist_matrix <- dist_matrix[ph_lang, ph_lang]
# phil_pairs excludes lang1 == lang2, so pivot_wider() leaves the diagonal NA;
# set those self-distances to 0. With a fully connected node graph every
# off-diagonal pair is reachable, so no other NAs arise here.
PHONEME_dist_matrix[is.na(PHONEME_dist_matrix)] <- 0

# Dijkstra routes i->j and j->i sum identical edge weights but in a different
# order, so the raw matrix can carry ~1e-12 floating-point asymmetry. Symmetrize
# so the MMRR symmetry guard below is exact (unfold() reads only the lower
# triangle regardless, so this does not change the fit).
PHONEME_dist_matrix <- (PHONEME_dist_matrix + t(PHONEME_dist_matrix)) / 2

# ---- 7. Phylogenetic distance matrix (X_phylo) + alignment ------------------
# Patristic (cophenetic) distances written by [0]_Phylogenetic_Tree.R, keyed by
# language name. check.names = FALSE preserves names-with-spaces so they match
# the cosine/dist matrices exactly.
PHONEME_phylo_dist_matrix <- read.csv(
  here("data", "PHONEME_phylo_dist_matrix.csv"),
  row.names = 1, check.names = FALSE
) |>
  as.matrix()

# Align all three matrices to a common language set + identical ordering before
# vectorizing. The intersection is defensive: with the current data all 58
# Philippine languages are shared, but this keeps the script correct if the tree
# set ever changes upstream.
common <- Reduce(intersect, list(
  rownames(PHONEME_diss_matrix),
  rownames(PHONEME_dist_matrix),
  rownames(PHONEME_phylo_dist_matrix)
))

dropped <- setdiff(
  union(rownames(PHONEME_diss_matrix), rownames(PHONEME_phylo_dist_matrix)),
  common
)
if (length(dropped) > 0) {
  message("MMRR alignment dropped ", length(dropped),
          " language(s) not shared across all matrices: ",
          paste(dropped, collapse = ", "))
}

Y_ling  <- PHONEME_diss_matrix[common, common]
X_geo   <- PHONEME_dist_matrix[common, common]
X_phylo <- PHONEME_phylo_dist_matrix[common, common]

# Self-distances are exactly zero by definition. 1 - cosine leaves ~1e-10 of
# floating-point noise on Y_ling's diagonal (cosine self-similarity is 1 up to
# rounding); zero all three diagonals so the guards below are exact. unfold()
# only reads the strict lower triangle, so this never affects the fit.
diag(Y_ling) <- 0
diag(X_geo) <- 0
diag(X_phylo) <- 0

# Guards: MMRR assumes symmetric, zero-diagonal matrices sharing one label order.
stopifnot(
  "Y/X_geo labels differ"    = identical(dimnames(Y_ling), dimnames(X_geo)),
  "Y/X_phylo labels differ"  = identical(dimnames(Y_ling), dimnames(X_phylo)),
  "Y not square"             = nrow(Y_ling) == ncol(Y_ling),
  "Y_ling not symmetric"     = isSymmetric(unname(Y_ling)),
  "X_geo not symmetric"      = isSymmetric(unname(X_geo)),
  "X_phylo not symmetric"    = isSymmetric(unname(X_phylo)),
  "Y_ling diagonal nonzero"  = all(diag(Y_ling)   == 0),
  "X_geo diagonal nonzero"   = all(diag(X_geo)    == 0),
  "X_phylo diagonal nonzero" = all(diag(X_phylo)  == 0)
)
message("MMRR input: ", nrow(Y_ling), " languages, ",
        choose(nrow(Y_ling), 2), " pairwise comparisons.")

# ---- 8. Matrix heatmaps (quick visual QC) -----------------------------------
# Modernized off reshape2::melt -> as_tibble(rownames) + pivot_longer.
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
  heatmap_p(Y_ling,  "Phonemic dissimilarity") +
  heatmap_p(X_phylo, "Phylogenetic distance")
)

# ---- 9. MMRR (Wang 2013) ----------------------------------------------------
# Canonical Multiple Matrix Regression with Randomization from the Wang Lab
# (landscapegenetics.org; Dryad doi:10.5061/dryad.kt71r/1; packaged as
# algatr::mmrr). Reproduced here, lightly adapted, so the analysis has no
# external dependency. Reference: Wang, I.J. (2013) "Examining the full effects
# of landscape heterogeneity on spatial genetic variation." Evolution 67:3403.
#
# Ordinary regression p-values are invalid because pairwise matrix cells are not
# independent (each language participates in many pairs). Significance is instead
# assessed by relabeling the rows AND columns of Y with a single permutation
# vector (Y[rand, rand]) — this preserves symmetry (Y_ij = Y_ji) and the zero
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

set.seed(1)  # reproducible permutation p-values
mmrr_fit <- MMRR(Y_ling, list(geo = X_geo, phylo = X_phylo), nperm = 9999)
print(mmrr_fit)

# ---- 10. Tidy results tibble ------------------------------------------------
PHONEME_mmrr_results <- tibble(
  dataset    = "PHONEME",
  beta_geo   = unname(mmrr_fit$coefficients["geo"]),
  p_geo      = unname(mmrr_fit$tpvalue["geo"]),
  beta_phylo = unname(mmrr_fit$coefficients["phylo"]),
  p_phylo    = unname(mmrr_fit$tpvalue["phylo"]),
  r_squared  = mmrr_fit$r.squared,
  p_model    = unname(mmrr_fit$Fpvalue)
)
print(PHONEME_mmrr_results)
write.csv(PHONEME_mmrr_results,
          file = here("data", "PHONEME_mmrr_results.csv"), row.names = FALSE)

# ---- 11. Visualization ------------------------------------------------------
# Both figures are built on the standardized unfolded lower triangles, so they
# read on the same scale as the fitted (standardized) MMRR coefficients.
dir.create(here("figures", "phoneme", "mmrr"), recursive = TRUE, showWarnings = FALSE)

mmrr_df <- tibble(
  ling  = as.numeric(unfold(Y_ling)),
  geo   = as.numeric(unfold(X_geo)),
  phylo = as.numeric(unfold(X_phylo))
)

# (a) Pairplot / scatterplot matrix, laid out as a single 3x3 facet_grid so every
# cell shares its column (x) and row (y) scale and the panel edges line up exactly
# — unlike a patchwork of nine independent plots, whose differing axes never
# aligned. Variable names are the facet strips (top = column var, right = row
# var). Lower triangle: scatter + linear fit. Upper triangle: Pearson r (the
# geo<->phylo collinearity that motivates the joint model, mirroring the check in
# [4] §9). Diagonal: the variable's marginal density. Built directly rather than
# with GGally::ggpairs to keep the script dependency-free.
pair_labs <- c(ling = "Dissimilarity", geo = "Geo distance", phylo = "Phylo distance")
pair_vars <- names(pair_labs)
pair_idx  <- setNames(seq_along(pair_vars), pair_vars)
as_pair_factor <- function(v) factor(pair_labs[v], levels = pair_labs)

# Cell text is centred on each variable's midpoint so it sits inside every
# free-scaled panel (the standardized inputs centre near 0 but are mildly
# right-skewed, so the panel middle is not exactly 0).
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

# Anchor every panel's free scale to the col-var (x) and row-var (y) ranges. The
# text/density-only cells — the first row and last column carry no scatter —
# would otherwise collapse to a degenerate axis; this pins all panels to the true
# variable ranges, as a proper scatterplot matrix does.
anchor_df <- expand_grid(row = pair_vars, col = pair_vars) |>
  mutate(dat = map2(row, col, ~ tibble(x = range(mmrr_df[[.y]]), y = range(mmrr_df[[.x]])))) |>
  unnest(dat) |>
  mutate(row = as_pair_factor(row), col = as_pair_factor(col))

# Diagonal: marginal density of each variable. The row's y-scale is the
# variable's value range (not a 0-1 density axis), so the density height is
# rescaled into the lower ~85% of the panel (baseline at the variable's min) —
# the shape shows without breaking the shared scales. from/to clip the curve to
# the data range so its tails stay inside the anchored panel.
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
ggsave(here("figures", "phoneme", "mmrr", "phoneme_mmrr_pairplot.png"),
       pairplot, width = 6.5, height = 6, units = "in", dpi = 300)

# (b) Added-variable (partial regression) plots. Chosen over a coefficient plot:
# with two predictors and a single dataset these show the ACTUAL isolation-by-
# distance relationship holding the other predictor constant (the scientific
# claim), and each panel's fitted slope equals the corresponding MMRR beta. A
# 2-point coefficient plot conveys less, and permutation yields a null
# distribution rather than a natural CI around the estimate.
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
  ylab = "Dissimilarity | phylo (residuals)",
  beta = PHONEME_mmrr_results$beta_geo, pval = PHONEME_mmrr_results$p_geo
)

p_phylo_av <- av_plot(
  mmrr_df$ling, mmrr_df$phylo, mmrr_df$geo,
  xlab = "Phylo distance | geo (residuals)",
  ylab = "Dissimilarity | geo (residuals)",
  beta = PHONEME_mmrr_results$beta_phylo, pval = PHONEME_mmrr_results$p_phylo
)

partial_plot <- (p_geo_av + p_phylo_av) +
  plot_annotation(
    title = "MMRR partial regression: phonemic dissimilarity vs. geography + phylogeny",
    subtitle = sprintf("Model R^2 = %.3f,  permutation p = %.4f  (9,999 permutations)",
                       PHONEME_mmrr_results$r_squared, PHONEME_mmrr_results$p_model)
  )
print(partial_plot)
ggsave(here("figures", "phoneme", "mmrr", "phoneme_mmrr_partial_regression.png"),
       partial_plot, width = 10, height = 4.5, units = "in", dpi = 300)

# ---- 12. Persist matrices ---------------------------------------------------
write.csv(PHONEME_dist_matrix, file = here("data", "PHONEME_dist_matrix.csv"), row.names = TRUE)
write.csv(PHONEME_diss_matrix, file = here("data", "PHONEME_diss_matrix.csv"), row.names = TRUE)
