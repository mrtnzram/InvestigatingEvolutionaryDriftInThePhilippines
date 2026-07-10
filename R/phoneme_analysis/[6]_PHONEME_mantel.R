# =============================================================================
# [5] Phoneme Analysis — Mantel test (isolation by distance)
# Tests whether pairwise phonemic dissimilarity between Philippine languages
# tracks pairwise terrain-penalized travel distance through the waypoint network.
#
# This file keeps the original Dijkstra routing (shortest_path_trace) to build
# the full pairwise distance matrix — it is intentionally NOT the simplified
# "distance to nearest node" used in [3]_PHONEME_network_distance.R, because a
# Mantel test needs genuine language-to-language distances.
#
# Inputs:  data/PHONEME_cosine_matrix.csv, data/RUHLENdf_PH.csv,
#          data/nodes.csv, data/edges.csv
# Outputs: data/PHONEME_dist_matrix.csv, data/PHONEME_diss_matrix.csv
# =============================================================================

library(readr)
library(tidyverse)
library(dplyr)
library(geosphere)
library(sf)
library(maps)
library(sfheaders)
library(reshape2)
library(vegan)
library(ggplot2)
library(here)

# ---- 0. Load data -----------------------------------------------------------
RUHLENdf <- read_csv(here("data", "RUHLENdf_PH.csv"))

ph_lang <- RUHLENdf %>%
  filter(Language_type == "Philippine Language") %>%
  pull(language)

# Cosine-similarity matrix written by [1]_PHONEME_cosine_similarity.R
cosine_matrix <- read.csv(here("data", "PHONEME_cosine_matrix.csv"),
                          row.names = 1, check.names = FALSE) %>%
  as.matrix()

nodes <- read.csv(here("data", "nodes.csv"))
edges <- read.csv(here("data", "edges.csv"))

nodes$id   <- as.character(nodes$id)
edges$from <- as.character(edges$from)
edges$to   <- as.character(edges$to)

# Add reverse edges for bidirectional routing
reverse_edges <- edges %>% rename(from = to, to = from)
edges <- bind_rows(edges, reverse_edges) %>% distinct()

# ---- 1. Dissimilarity matrix (1 - cosine), Philippine languages -------------
cosine_matrix_phil <- cosine_matrix[ph_lang, ph_lang]
PHONEME_diss_matrix <- 1 - cosine_matrix_phil

# ---- 2. Land mask (Philippines + Malaysia) ----------------------------------
world_map <- map_data("world") %>% filter(region %in% c("Philippines", "Malaysia"))

land_sf <- sf_polygon(obj = world_map, polygon_id = "group", x = "long", y = "lat") %>%
  st_union() %>%
  st_sf(geometry = .) %>%
  st_set_crs(4326)

# ---- 3. Terrain-penalized network edges + weighted graph --------------------
# NOTE: edges are penalized at 44.18 (matching the original [1] route weighting),
# while the language connectors below use 4.44. This penalty mismatch is carried
# over verbatim from the monolith to preserve results — flag for review.
edge_land_penalty <- 44.18

edges <- edges %>%
  rowwise() %>%
  mutate(weight = {
    from_coords <- nodes %>% filter(id == from)
    to_coords   <- nodes %>% filter(id == to)
    if (nrow(from_coords) == 0 || nrow(to_coords) == 0) NA_real_
    else distHaversine(c(from_coords$longitude, from_coords$latitude),
                       c(to_coords$longitude, to_coords$latitude))
  }) %>%
  ungroup()

edge_lines <- edges %>%
  rowwise() %>%
  mutate(
    geometry = list(st_linestring(matrix(c(
      nodes$longitude[nodes$id == from], nodes$latitude[nodes$id == from],
      nodes$longitude[nodes$id == to],   nodes$latitude[nodes$id == to]
    ), ncol = 2, byrow = TRUE)))
  ) %>%
  ungroup() %>%
  st_as_sf(crs = 4326)

edge_lines <- edge_lines %>%
  rowwise() %>%
  mutate(
    land_part = list(st_intersection(geometry, land_sf)),
    sea_part  = list(st_difference(geometry, land_sf)),

    land_len = as.numeric(if (!is.null(land_part) && length(land_part) > 0) st_length(land_part) else 0),
    sea_len  = as.numeric(if (!is.null(sea_part)  && length(sea_part)  > 0) st_length(sea_part)  else 0),

    weighted_cost = land_len * edge_land_penalty + sea_len,
    crosses_land  = land_len > 0
  ) %>%
  ungroup()

# Adjacency list keyed by node id, weighted by terrain-penalized cost
all_ids <- unique(c(edge_lines$from, edge_lines$to))
graph <- lapply(all_ids, function(id) {
  neighbors <- edge_lines %>% filter(from == id) %>% select(to, weighted_cost)
  if (nrow(neighbors) == 0) tibble(to = character(), weight = numeric())
  else rename(neighbors, weight = weighted_cost)
})
names(graph) <- all_ids

# ---- 4. Routing helpers -----------------------------------------------------
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

phil_df <- RUHLENdf %>%
  filter(Language_type == "Philippine Language") %>%
  mutate(
    start_coords = map2(longitude, latitude, ~ c(.x, .y)),
    nearest_node = map_chr(start_coords, find_nearest_node)
  )

connector_df <- phil_df %>%
  mutate(
    connector_geom = map2(start_coords, nearest_node, ~ st_linestring(rbind(
      .x,
      c(nodes$longitude[nodes$id == .y], nodes$latitude[nodes$id == .y])
    )))
  ) %>%
  mutate(connector_geom_sfc = st_sfc(connector_geom, crs = 4326)) %>%
  rowwise() %>%
  mutate(
    land_part = list(st_intersection(connector_geom_sfc, land_sf)),
    sea_part  = list(st_difference(connector_geom_sfc, land_sf)),

    land_len = as.numeric(if (!is.null(land_part) && length(land_part) > 0) st_length(land_part) else 0),
    sea_len  = as.numeric(if (!is.null(sea_part)  && length(sea_part)  > 0) st_length(sea_part)  else 0),

    connector_penalty = land_len * land_penalty + sea_len
  ) %>%
  ungroup()

# ---- 6. Pairwise terrain-penalized distance matrix --------------------------
phil_pairs <- expand_grid(lang1 = phil_df$language, lang2 = phil_df$language) %>%
  filter(lang1 != lang2) %>%
  left_join(phil_df %>% select(language, node1 = nearest_node), by = c("lang1" = "language")) %>%
  left_join(phil_df %>% select(language, node2 = nearest_node), by = c("lang2" = "language")) %>%
  left_join(connector_df %>% select(language, penalty1 = connector_penalty), by = c("lang1" = "language")) %>%
  left_join(connector_df %>% select(language, penalty2 = connector_penalty), by = c("lang2" = "language"))

phil_pairs <- phil_pairs %>%
  rowwise() %>%
  mutate(
    trace = list(shortest_path_trace(node1, node2, graph)),
    tree_dist = trace$distance,
    geodist_H1_span = if (is.na(tree_dist)) NA_real_ else
      (penalty1 + tree_dist + penalty2) / 1000
  ) %>%
  ungroup()

dist_matrix <- phil_pairs %>%
  select(lang1, lang2, geodist_H1_span) %>%
  pivot_wider(names_from = lang2, values_from = geodist_H1_span) %>%
  column_to_rownames("lang1") %>%
  as.matrix()

PHONEME_dist_matrix <- dist_matrix[ph_lang, ph_lang]
PHONEME_dist_matrix[is.na(PHONEME_dist_matrix)] <- 0

# ---- 7. Plot the two matrices -----------------------------------------------
melt_phoneme_dist_matrix <- melt(PHONEME_dist_matrix)
melt_phoneme_diss_matrix <- melt(PHONEME_diss_matrix)

dist_matrix_p <- ggplot(melt_phoneme_dist_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_fixed()

diss_matrix_p <- ggplot(melt_phoneme_diss_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_fixed()

dist_matrix_p + diss_matrix_p

# ---- 8. Mantel test ---------------------------------------------------------
# Convert to distance objects
x_dist <- as.dist(PHONEME_dist_matrix)
y_dist <- as.dist(PHONEME_diss_matrix)

# Run Mantel test
mantel_result <- mantel(x_dist, y_dist, method = "spearman", permutations = 999)
print(mantel_result)

# Convert matrices to vectors
x_vec <- as.vector(as.dist(PHONEME_dist_matrix))
y_vec <- as.vector(as.dist(PHONEME_diss_matrix))

ggplot(data.frame(x = x_vec, y = y_vec), aes(x = x, y = y)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  theme_bw() +
  labs(
    title = "Mantel Test: Phoneme Distance vs. Dissimilarity",
    x = "Relative Migration Pairwise Distance (km)",
    y = "Phonemic Dissimilarity"
  )

# ---- 9. Persist matrices ----------------------------------------------------
write.csv(PHONEME_dist_matrix, file = here("data", "PHONEME_dist_matrix.csv"), row.names = TRUE)
write.csv(PHONEME_diss_matrix, file = here("data", "PHONEME_diss_matrix.csv"), row.names = TRUE)
