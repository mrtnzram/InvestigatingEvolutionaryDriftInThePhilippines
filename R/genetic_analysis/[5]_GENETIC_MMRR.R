# ── [5]_GENETIC_MMRR.R ───────────────────────────────────────────────────────
# Multiple Matrix Regression with Randomization (MMRR), genetic analogue of
# [5]_PHONEME_MMRR.R / [5]_GRAMMAR_MMRR.R. Regresses pairwise population IBS
# similarity (from PLINK --distance 1-ibs) on terrain-penalized migration
# distance, with permutation p-values. Single predictor only — there is no
# genetic phylogeny to pair with geography here, so the joint two-predictor
# model of the language scripts collapses to one fit, and the added-variable
# (partial regression) figure is skipped: with one predictor it would just
# reproduce the plain scatter.
#
# Input:   data/genetic_dist.mdist, data/genetic_dist.mdist.id (PLINK, individual-
#          level), data/GENETIC_final.csv (from [3], population connector costs),
#          data/nodes.csv, data/edges.csv
# Outputs: data/GENETIC_dist_matrix.csv, data/GENETIC_sim_matrix.csv,
#          data/GENETIC_mmrr_results.csv,
#          figures/genetic/mmrr/genetic_mmrr_single_similarity_vs_geography.png,
#          figures/genetic/mmrr/genetic_mmrr_pairplot.png
# ──────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(geosphere)
library(sf)
library(maps)
library(sfheaders)
library(patchwork)
library(here)
library(conflicted)

conflicts_prefer(purrr::map)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

# ── 1. Y: individual PLINK matrix -> population IBS similarity ──────────────
# .mdist is a square, headerless, tab-delimited 1-IBS distance matrix over
# individuals; .id gives each row/col's (FID = population, IID = sample).
mdist_ids <- read_tsv(here("data", "genetic_dist.mdist.id"),
                      col_names = c("FID", "IID"), col_types = "cc")
M <- as.matrix(read_tsv(here("data", "genetic_dist.mdist"),
                        col_names = FALSE, col_types = cols(.default = "d")))
dimnames(M) <- NULL
stopifnot("mdist is not square" = nrow(M) == ncol(M),
          "mdist/.id row count differs" = nrow(M) == nrow(mdist_ids))

# Same reconciliation [0]_GENETIC_ADMX_MALDER.R applies to the qpAdm table, so
# labels here join GENETIC_final.csv's `population` column.
fid <- ifelse(mdist_ids$FID == "ManoboRajahKabunsuwan", "ManoboRK", mdist_ids$FID)

# Block-average individuals -> populations. Two rowsum() passes give the block
# SUMS (Z'MZ for the group-indicator matrix Z, without ever materializing Z);
# dividing by outer(n, n) turns each block sum into its mean. Off-diagonal
# cells are unbiased pairwise means; the diagonal is not (it includes each
# individual's own zero self-distance, understating within-population
# distance by a factor of 1 - 1/n_k) but that is moot below, since the
# diagonal is overwritten and unfold() only reads the strict lower triangle.
n_pop   <- table(fid)
pop_sum <- rowsum(t(rowsum(M, fid)), fid)
pop_dist <- pop_sum / outer(as.numeric(n_pop), as.numeric(n_pop))
dimnames(pop_dist) <- list(rownames(pop_sum), rownames(pop_sum))

stopifnot(
  "population count != 115"   = nrow(pop_dist) == 115,
  "pop_dist not symmetric"    = isSymmetric(unname(pop_dist)),
  "pop_dist out of [0, 1]"    = all(pop_dist >= 0 & pop_dist <= 1)
)
message(sprintf(
  "Y_gen: %d populations, %d-%d individuals each (Arta/Batak thinnest at 3).",
  nrow(pop_dist), min(n_pop), max(n_pop)
))

diag(pop_dist) <- 0
GENETIC_sim_matrix <- 1 - pop_dist   # IBS proportion; diagonal = 1 like the cosine matrices

# ── 2. X_geo: terrain-penalized pairwise migration distance ─────────────────
GENETIC_final <- read.csv(here("data", "GENETIC_final.csv"))
stopifnot("GENETIC_final.csv population set != mdist" =
            setequal(GENETIC_final$population, rownames(GENETIC_sim_matrix)))

nodes <- read.csv(here("data", "nodes.csv")) |> mutate(id = as.character(id))
edges <- read.csv(here("data", "edges.csv")) |>
  mutate(from = as.character(from), to = as.character(to))

world_map <- map_data("world") |> filter(region %in% c("Philippines", "Malaysia"))
land_sf <- sf_polygon(obj = world_map, polygon_id = "group", x = "long", y = "lat") |>
  st_union() |>
  st_sf(geometry = _) |>
  st_set_crs(4326)

LAND_PENALTY <- 4.44   # matches [3] and the corrected phoneme/grammar [5] scripts

# split_and_penalize() / build_network_edges(): copied from
# [3]_GENETIC_network_distance.R so the node-to-node network geometry here is
# built exactly the way [3] built the connector costs this script reuses below.
split_and_penalize <- function(geom, land_sf, land_penalty) {
  land_part <- st_intersection(geom, st_geometry(land_sf))
  sea_part  <- st_difference(geom, st_geometry(land_sf))
  land_part <- land_part[!st_is_empty(land_part)]
  sea_part  <- sea_part[!st_is_empty(sea_part)]
  land_len <- if (length(land_part) > 0) as.numeric(sum(st_length(land_part))) else 0
  sea_len  <- if (length(sea_part)  > 0) as.numeric(sum(st_length(sea_part)))  else 0
  list(land_len = land_len, sea_len = sea_len,
       weighted_cost = land_len * land_penalty + sea_len)
}

build_network_edges <- function(nodes, edges, land_sf, land_penalty) {
  edges_bi <- bind_rows(edges, edges |> rename(from = to, to = from)) |> distinct()
  edges_sf <- edges_bi |>
    rowwise() |>
    mutate(geometry = list(st_linestring(matrix(c(
      nodes$longitude[nodes$id == from], nodes$latitude[nodes$id == from],
      nodes$longitude[nodes$id == to],   nodes$latitude[nodes$id == to]
    ), ncol = 2, byrow = TRUE)))) |>
    ungroup() |>
    st_as_sf(crs = 4326)
  edges_sf |>
    rowwise() |>
    mutate(.split = list(split_and_penalize(geometry, land_sf, land_penalty)),
           weighted_cost = .split$weighted_cost) |>
    ungroup() |>
    select(-.split) |>
    st_as_sf()
}

network_edges <- build_network_edges(nodes, edges, land_sf, LAND_PENALTY)

# All-pairs node-to-node cost on the 39-node waypoint network, one igraph call
# (the [3] pattern) rather than a per-pair Dijkstra loop.
g <- igraph::graph_from_data_frame(
  d = network_edges |> st_drop_geometry() |> select(from, to, weighted_cost),
  directed = FALSE, vertices = nodes
)
node_dist_m <- igraph::distances(g, weights = igraph::E(g)$weighted_cost)

# Per-population connector cost (point -> nearest network node) and node id:
# already computed by [3] as geodist_H1_span (km) / nearest_node — reused
# rather than rebuilt, so this script doesn't re-derive [3]'s connector geometry.
pop_geo <- GENETIC_final |>
  select(population, nearest_node, connector_km = geodist_H1_span, longitude, latitude) |>
  mutate(nearest_node = as.character(nearest_node))

# Direct-line cost, split into land/sea exactly like [3]'s connectors and the
# language [5] scripts' direct option. Degenerate zero-length lines (the 4
# coordinate sites shared by >1 population) short-circuit to 0.
direct_penalized_cost <- function(lon1, lat1, lon2, lat2) {
  if (isTRUE(all.equal(c(lon1, lat1), c(lon2, lat2)))) return(0)
  geom <- st_sfc(st_linestring(rbind(c(lon1, lat1), c(lon2, lat2))), crs = 4326)
  split_and_penalize(geom, land_sf, LAND_PENALTY)$weighted_cost
}

# Unordered pairs only (Y and X are both symmetric) — 6,555 pairs vs. 13,110
# for the full ordered grid the language scripts use at n ~ 58.
pop_pairs <- as.data.frame(t(combn(pop_geo$population, 2)), stringsAsFactors = FALSE) |>
  set_names(c("pop1", "pop2")) |>
  as_tibble() |>
  left_join(pop_geo |> select(population, node1 = nearest_node, conn1 = connector_km,
                              lon1 = longitude, lat1 = latitude),
            by = c("pop1" = "population")) |>
  left_join(pop_geo |> select(population, node2 = nearest_node, conn2 = connector_km,
                              lon2 = longitude, lat2 = latitude),
            by = c("pop2" = "population"))

pop_pairs <- pop_pairs |>
  rowwise() |>
  mutate(
    node_to_node_km  = node_dist_m[node1, node2] / 1000,
    via_network_cost = conn1 + node_to_node_km + conn2,
    direct_cost      = direct_penalized_cost(lon1, lat1, lon2, lat2) / 1000,
    geodist_H1_span  = pmin(via_network_cost, direct_cost, na.rm = TRUE),
    used_direct      = !is.na(direct_cost) &
      (is.na(via_network_cost) | direct_cost < via_network_cost)
  ) |>
  ungroup()

message(sprintf(
  "X_geo routing: %d of %d population pairs (%.0f%%) took the direct line rather than the network.",
  sum(pop_pairs$used_direct), nrow(pop_pairs), 100 * mean(pop_pairs$used_direct)
))

pop_order <- rownames(GENETIC_sim_matrix)
GENETIC_dist_matrix <- matrix(0, nrow = length(pop_order), ncol = length(pop_order),
                              dimnames = list(pop_order, pop_order))
GENETIC_dist_matrix[cbind(pop_pairs$pop1, pop_pairs$pop2)] <- pop_pairs$geodist_H1_span
GENETIC_dist_matrix[cbind(pop_pairs$pop2, pop_pairs$pop1)] <- pop_pairs$geodist_H1_span

# ── 3. Alignment + guards ────────────────────────────────────────────────────
common <- intersect(rownames(GENETIC_sim_matrix), rownames(GENETIC_dist_matrix))
Y_gen <- GENETIC_sim_matrix[common, common]
X_geo <- GENETIC_dist_matrix[common, common]
diag(X_geo) <- 0

stopifnot(
  "Y/X_geo labels differ" = identical(dimnames(Y_gen), dimnames(X_geo)),
  "Y not square"          = nrow(Y_gen) == ncol(Y_gen),
  "Y_gen not symmetric"   = isSymmetric(unname(Y_gen)),
  "X_geo not symmetric"   = isSymmetric(unname(X_geo)),
  "Y_gen out of [0, 1]"   = all(Y_gen >= 0 & Y_gen <= 1),
  "X_geo diagonal nonzero" = all(diag(X_geo) == 0)
)
message("MMRR input: ", nrow(Y_gen), " populations, ",
        choose(nrow(Y_gen), 2), " pairwise comparisons.")

# ── 4. Matrix heatmaps (quick visual QC) ─────────────────────────────────────
melt_matrix <- function(m) {
  as_tibble(m, rownames = "Var1") |>
    pivot_longer(-Var1, names_to = "Var2", values_to = "value")
}
heatmap_p <- function(m, title) {
  ggplot(melt_matrix(m), aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "yellow", high = "red") +
    labs(title = title, x = "", y = "") +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
    coord_fixed()
}
print(heatmap_p(X_geo, "Migration distance") + heatmap_p(Y_gen, "IBS similarity"))

# ── 5. MMRR (Wang 2013), single predictor ────────────────────────────────────
# unfold()/MMRR() copied verbatim from [5]_PHONEME_MMRR.R (Wang Lab's canonical
# MMRR, reproduced here so the script has no external dependency).
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
mmrr_fit <- MMRR(Y_gen, list(geo = X_geo), nperm = 9999)
print(mmrr_fit)

GENETIC_mmrr_results <- tibble(
  dataset   = "GENETIC",
  beta_geo  = unname(mmrr_fit$coefficients["geo"]),
  p_geo     = unname(mmrr_fit$tpvalue["geo"]),
  r_squared = mmrr_fit$r.squared,
  p_model   = unname(mmrr_fit$Fpvalue)
)
print(GENETIC_mmrr_results)
write.csv(GENETIC_mmrr_results,
          file = here("data", "GENETIC_mmrr_results.csv"), row.names = FALSE)

# ── 6. Visualization ─────────────────────────────────────────────────────────
dir.create(here("figures", "genetic", "mmrr"), recursive = TRUE, showWarnings = FALSE)

# (a) Single-predictor scatter, raw axes (km, IBS proportion) so it reads
# cleanly; beta/R^2 come from the standardized fit above, not refit here.
p_txt <- if (GENETIC_mmrr_results$p_geo <= 1e-4) "p < 0.0001" else
  sprintf("p = %.4f", GENETIC_mmrr_results$p_geo)

mmrr_df <- tibble(
  gen = as.numeric(unfold(Y_gen)),
  geo = as.numeric(unfold(X_geo))
)
raw_df <- tibble(x = unfold(X_geo, scale = FALSE), y = unfold(Y_gen, scale = FALSE))

single_plot <- ggplot(raw_df, aes(x, y)) +
  geom_point(alpha = 0.15, size = 0.5, colour = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, colour = "firebrick", linewidth = 0.9) +
  labs(
    title    = "Genetic similarity vs Geographic distance",
    subtitle = sprintf("beta = %+.3f (standardized)   %s   R² = %.4f",
                       GENETIC_mmrr_results$beta_geo, p_txt, GENETIC_mmrr_results$r_squared),
    x = "Relative migration distance (km)", y = "IBS similarity"
  ) +
  theme_bw() +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10, colour = "grey25"))
print(single_plot)
ggsave(here("figures", "genetic", "mmrr", "genetic_mmrr_single_similarity_vs_geography.png"),
       single_plot, width = 6.5, height = 4.5, units = "in", dpi = 300)

# (b) 2x2 pairplot — same facet_grid machinery as the 3x3 language pairplots,
# generic over pair_vars, so a 2-element pair_labs collapses it to one
# scatter cell, one Pearson-r cell, and two marginal densities.
pair_labs <- c(gen = "IBS similarity", geo = "Geo distance")
pair_vars <- names(pair_labs)
pair_idx  <- setNames(seq_along(pair_vars), pair_vars)
as_pair_factor <- function(v) factor(pair_labs[v], levels = pair_labs)

pair_mid <- colMeans(sapply(mmrr_df[pair_vars], range))

scatter_df <- expand_grid(row = pair_vars, col = pair_vars) |>
  filter(pair_idx[row] > pair_idx[col]) |>
  mutate(dat = map2(row, col, ~ tibble(x = mmrr_df[[.y]], y = mmrr_df[[.x]]))) |>
  unnest(dat) |>
  mutate(row = as_pair_factor(row), col = as_pair_factor(col))

cor_df <- expand_grid(row = pair_vars, col = pair_vars) |>
  filter(pair_idx[row] < pair_idx[col]) |>
  mutate(
    lab = sprintf("r = %.2f", map2_dbl(row, col, ~ cor(mmrr_df[[.x]], mmrr_df[[.y]]))),
    x   = pair_mid[col], y = pair_mid[row],
    row = as_pair_factor(row), col = as_pair_factor(col)
  )

anchor_df <- expand_grid(row = pair_vars, col = pair_vars) |>
  mutate(dat = map2(row, col, ~ tibble(x = range(mmrr_df[[.y]]), y = range(mmrr_df[[.x]])))) |>
  unnest(dat) |>
  mutate(row = as_pair_factor(row), col = as_pair_factor(col))

density_df <- map_dfr(pair_vars, function(v) {
  rng <- range(mmrr_df[[v]])
  d   <- density(mmrr_df[[v]], from = rng[1], to = rng[2])
  tibble(v = v, x = d$x, ymin = rng[1],
         ymax = rng[1] + (d$y / max(d$y)) * 0.85 * diff(rng))
}) |>
  mutate(row = as_pair_factor(v), col = as_pair_factor(v))

pairplot <- ggplot() +
  geom_blank(data = anchor_df, aes(x, y)) +
  geom_ribbon(data = density_df, aes(x = x, ymin = ymin, ymax = ymax),
              fill = "grey80", color = "grey40", linewidth = 0.3) +
  geom_point(data = scatter_df, aes(x, y),
             alpha = 0.3, size = 0.5, color = "steelblue") +
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
ggsave(here("figures", "genetic", "mmrr", "genetic_mmrr_pairplot.png"),
       pairplot, width = 6.5, height = 6, units = "in", dpi = 300)

# ── 7. Persist matrices ──────────────────────────────────────────────────────
write.csv(GENETIC_dist_matrix, file = here("data", "GENETIC_dist_matrix.csv"), row.names = TRUE)
write.csv(GENETIC_sim_matrix, file = here("data", "GENETIC_sim_matrix.csv"), row.names = TRUE)
