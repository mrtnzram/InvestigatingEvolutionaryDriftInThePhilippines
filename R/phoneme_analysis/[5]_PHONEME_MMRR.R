# =============================================================================
# [5] Phoneme Analysis — Multiple Matrix Regression with Randomization (MMRR)
#
# Control analysis for internal linguistic diffusion. Regresses pairwise general
# phonemic similarity (cosine, not Spanish-specific) jointly on terrain-penalized
# migration distance and patristic phylogenetic distance, with permutation
# p-values, so geography and shared ancestry — collinear here — are assessed
# together rather than by a single-predictor Mantel test.
#
# Note: the pairwise distance matrix is built by Dijkstra routing over the full
# language-to-language network, not the nearest-node distance used in [3]. Each
# pair takes the cheaper of the via-network route and a direct terrain-penalized
# line (§6), under a single land penalty shared by every leg (§3).
#
# Inputs:  data/PHONEME_cosine_matrix.csv, data/RUHLENdf_PH.csv,
#          data/nodes.csv, data/edges.csv,
#          data/PHONEME_phylo_dist_matrix.csv  (written by [0]_Phylogenetic_Tree.R)
# Outputs: data/PHONEME_dist_matrix.csv, data/PHONEME_sim_matrix.csv,
#          data/PHONEME_mmrr_results.csv        (joint two-predictor model)
#          data/PHONEME_mmrr_single_results.csv (single-predictor models)
#          figures/phoneme/mmrr/phoneme_mmrr_single_*.png  (3 single-predictor)
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

# ---- 1. Similarity matrix (cosine), Philippine languages --------------------
cosine_matrix_phil <- cosine_matrix[ph_lang, ph_lang]
PHONEME_sim_matrix <- cosine_matrix_phil

# ---- 2. Land mask (Philippines + Malaysia) ----------------------------------
world_map <- map_data("world") |> filter(region %in% c("Philippines", "Malaysia"))

land_sf <- sf_polygon(obj = world_map, polygon_id = "group", x = "long", y = "lat") |>
  st_union() |>
  st_sf(geometry = _) |>
  st_set_crs(4326)

# ---- 3. Terrain-penalized network edges + weighted graph --------------------
# ONE terrain penalty for every leg of every route — network edges, language
# connectors, and the direct language-to-language line in §6. Land travel costs
# LAND_PENALTY x sea travel per metre, which is a property of the terrain, not
# of which layer a segment belongs to. This resolves the mismatch inherited from
# the original monolith (edges at 44.18 vs connectors at 4.44), which made a
# metre of land 10x dearer on the network than on a connector and left the
# routing internally inconsistent. 4.44 is the value [3] already uses.
LAND_PENALTY <- 4.44
edge_land_penalty <- LAND_PENALTY

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
# shortest_path_trace() is a hand-rolled Dijkstra kept verbatim from the original
# Mantel script so X_geo stays bit-for-bit reproducible; igraph::distances() on
# the same weighted graph (the [3] pattern) is the faster equivalent.
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
land_penalty <- LAND_PENALTY   # same terrain penalty as the network edges (§3)

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
# Each pair takes the cheaper of two candidate routes (a true least-cost choice,
# mirroring the pmin() [3] uses for its H2 measure):
#
#   via-network : connector1 + shortest path(node1 -> node2) + connector2
#   direct      : the terrain-penalized straight line between the two languages
#
# Without the direct option every pair was forced out to the network and back,
# which is indefensible for nearby languages: the 58 languages share only 17
# nearest nodes, so 184 pairs (11%) have node1 == node2 and the network leg
# contributes 0 — leaving a pure out-and-back detour. Ivatanen and Itbayaten sit
# at identical coordinates yet scored 257 km; South Kalinga and Central Kalinga
# are 10 km apart and scored 725 km. Both legs are penalized with the same
# LAND_PENALTY (§3), so the comparison is like-for-like.
phil_pairs <- expand_grid(lang1 = phil_df$language, lang2 = phil_df$language) |>
  filter(lang1 != lang2) |>
  left_join(phil_df |> select(language, node1 = nearest_node), by = c("lang1" = "language")) |>
  left_join(phil_df |> select(language, node2 = nearest_node), by = c("lang2" = "language")) |>
  left_join(connector_df |> select(language, penalty1 = connector_penalty), by = c("lang1" = "language")) |>
  left_join(connector_df |> select(language, penalty2 = connector_penalty), by = c("lang2" = "language")) |>
  left_join(phil_df |> select(language, lon1 = longitude, lat1 = latitude), by = c("lang1" = "language")) |>
  left_join(phil_df |> select(language, lon2 = longitude, lat2 = latitude), by = c("lang2" = "language"))

# Direct-line cost, split into land/sea exactly like the connectors in §5.
# Degenerate zero-length lines (languages recorded at identical coordinates)
# would make st_linestring/st_intersection misbehave, so they short-circuit to 0.
direct_penalized_cost <- function(lon1, lat1, lon2, lat2) {
  if (isTRUE(all.equal(c(lon1, lat1), c(lon2, lat2)))) return(0)
  geom <- st_sfc(st_linestring(rbind(c(lon1, lat1), c(lon2, lat2))), crs = 4326)
  land_part <- st_intersection(geom, land_sf)
  sea_part  <- st_difference(geom, land_sf)
  land_part <- land_part[!st_is_empty(land_part)]
  sea_part  <- sea_part[!st_is_empty(sea_part)]
  land_len <- if (length(land_part) > 0) as.numeric(sum(st_length(land_part))) else 0
  sea_len  <- if (length(sea_part)  > 0) as.numeric(sum(st_length(sea_part)))  else 0
  land_len * LAND_PENALTY + sea_len
}

phil_pairs <- phil_pairs |>
  rowwise() |>
  mutate(
    trace = list(shortest_path_trace(node1, node2, graph)),
    tree_dist = trace$distance,
    via_network_cost = if (is.na(tree_dist)) NA_real_ else penalty1 + tree_dist + penalty2,
    direct_cost      = direct_penalized_cost(lon1, lat1, lon2, lat2),
    # na.rm: an unreachable network route must not veto the direct one.
    geodist_H1_span  = pmin(via_network_cost, direct_cost, na.rm = TRUE) / 1000,
    used_direct      = !is.na(direct_cost) &
      (is.na(via_network_cost) | direct_cost < via_network_cost)
  ) |>
  ungroup()

message(sprintf(
  "X_geo routing: %d of %d ordered pairs (%.0f%%) took the direct line rather than the network.",
  sum(phil_pairs$used_direct), nrow(phil_pairs), 100 * mean(phil_pairs$used_direct)
))

dist_matrix <- phil_pairs |>
  select(lang1, lang2, geodist_H1_span) |>
  pivot_wider(names_from = lang2, values_from = geodist_H1_span) |>
  column_to_rownames("lang1") |>
  as.matrix()

PHONEME_dist_matrix <- dist_matrix[ph_lang, ph_lang]
# phil_pairs excludes lang1 == lang2, so the diagonal comes back NA from
# pivot_wider(); the graph is connected, so no other NAs arise.
PHONEME_dist_matrix[is.na(PHONEME_dist_matrix)] <- 0

# i->j and j->i sum the same weights in a different order, leaving ~1e-12 of
# asymmetry; symmetrize so the MMRR guard below is exact.
PHONEME_dist_matrix <- (PHONEME_dist_matrix + t(PHONEME_dist_matrix)) / 2

# ---- 7. Phylogenetic distance matrix (X_phylo) + alignment ------------------
# check.names = FALSE preserves names-with-spaces so the labels still match the
# cosine/dist matrices exactly.
PHONEME_phylo_dist_matrix <- read.csv(
  here("data", "PHONEME_phylo_dist_matrix.csv"),
  row.names = 1, check.names = FALSE
) |>
  as.matrix()

# Align all three matrices to a common language set and ordering before
# vectorizing; the intersection is defensive against upstream tree-set changes.
common <- Reduce(intersect, list(
  rownames(PHONEME_sim_matrix),
  rownames(PHONEME_dist_matrix),
  rownames(PHONEME_phylo_dist_matrix)
))

dropped <- setdiff(
  union(rownames(PHONEME_sim_matrix), rownames(PHONEME_phylo_dist_matrix)),
  common
)
if (length(dropped) > 0) {
  message("MMRR alignment dropped ", length(dropped),
          " language(s) not shared across all matrices: ",
          paste(dropped, collapse = ", "))
}

Y_ling  <- PHONEME_sim_matrix[common, common]
X_geo   <- PHONEME_dist_matrix[common, common]
X_phylo <- PHONEME_phylo_dist_matrix[common, common]

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
  heatmap_p(Y_ling,  "Phonemic similarity") +
  heatmap_p(X_phylo, "Phylogenetic distance")
)

# ---- 9. MMRR (Wang 2013) ----------------------------------------------------
# Wang Lab's canonical MMRR (Dryad doi:10.5061/dryad.kt71r/1, packaged as
# algatr::mmrr), reproduced here so the script has no external dependency.
# Reference: Wang, I.J. (2013) "Examining the full effects of landscape
# heterogeneity on spatial genetic variation." Evolution 67:3403.
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

# dir.create() hoisted here, before the first ggsave() in §9a — it used to sit
# at §11 and left a fresh checkout failing on the very first single-predictor plot.
dir.create(here("figures", "phoneme", "mmrr"), recursive = TRUE, showWarnings = FALSE)

# ---- 9a. Single-predictor MMRRs (run BEFORE the joint model) ----------------
# Same MMRR() above, just handed a one-element predictor list — the permutation
# scheme is identical, so with one predictor this reduces to a permutation-tested
# simple regression on the vectorized lower triangles, i.e. exactly the Mantel
# test this script replaced. Because unfold() standardizes, the single-predictor
# beta IS the Pearson r and R^2 = beta^2.
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

  # Axes are plotted in RAW units (km, patristic, 1 - cosine) so they read
  # cleanly; standardizing only rescales the axes, leaving the fit identical.
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
  ggsave(here("figures", "phoneme", "mmrr", file), p,
         width = 6.5, height = 4.5, units = "in", dpi = 300)

  tibble(model = title, beta = beta, p_value = pval, r_squared = r2)
}

PHONEME_mmrr_single <- bind_rows(
  single_mmrr(Y_ling, X_geo,
              y_lab = "Phonemic similarity (cosine)",
              x_lab = "Relative migration distance (km)",
              title = "Linguistic similarity vs Geographic distance",
              file  = "phoneme_mmrr_single_similarity_vs_geography.png"),
  single_mmrr(Y_ling, X_phylo,
              y_lab = "Phonemic similarity (cosine)",
              x_lab = "Patristic phylogenetic distance",
              title = "Linguistic similarity vs Phylogenetic distance",
              file  = "phoneme_mmrr_single_similarity_vs_phylogeny.png"),
  # Collinearity check: the shared structure that flattens both joint betas.
  single_mmrr(X_geo, X_phylo,
              y_lab = "Relative migration distance (km)",
              x_lab = "Patristic phylogenetic distance",
              title = "Geographic distance vs Phylogenetic distance",
              file  = "phoneme_mmrr_single_geography_vs_phylogeny.png")
)
print(PHONEME_mmrr_single)
write.csv(PHONEME_mmrr_single,
          file = here("data", "PHONEME_mmrr_single_results.csv"), row.names = FALSE)


# ---- 9b. Joint two-predictor MMRR -------------------------------------------
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
# (figures/phoneme/mmrr already created in §9a, ahead of the first ggsave())

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
ggsave(here("figures", "phoneme", "mmrr", "phoneme_mmrr_pairplot.png"),
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
  beta = PHONEME_mmrr_results$beta_geo, pval = PHONEME_mmrr_results$p_geo
)

p_phylo_av <- av_plot(
  mmrr_df$ling, mmrr_df$phylo, mmrr_df$geo,
  xlab = "Phylo distance | geo (residuals)",
  ylab = "Similarity | geo (residuals)",
  beta = PHONEME_mmrr_results$beta_phylo, pval = PHONEME_mmrr_results$p_phylo
)

partial_plot <- (p_geo_av + p_phylo_av) +
  plot_annotation(
    title = "MMRR partial regression: phonemic similarity vs. geography + phylogeny",
    subtitle = sprintf("Model R^2 = %.3f,  permutation p = %.4f  (9,999 permutations)",
                       PHONEME_mmrr_results$r_squared, PHONEME_mmrr_results$p_model)
  )
print(partial_plot)
ggsave(here("figures", "phoneme", "mmrr", "phoneme_mmrr_partial_regression.png"),
       partial_plot, width = 10, height = 4.5, units = "in", dpi = 300)

# ---- 12. Persist matrices ---------------------------------------------------
write.csv(PHONEME_dist_matrix, file = here("data", "PHONEME_dist_matrix.csv"), row.names = TRUE)
write.csv(PHONEME_sim_matrix, file = here("data", "PHONEME_sim_matrix.csv"), row.names = TRUE)
