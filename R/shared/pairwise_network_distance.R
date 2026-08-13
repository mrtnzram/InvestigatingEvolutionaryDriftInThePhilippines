# =============================================================================
# Shared — Pairwise terrain-penalized migration-cost matrix
#
# Builds the language/population pairwise cost matrix, called once per domain
# from [3]_*_network_distance.R. Direct line wins only when
# direct_cost <= direct_ratio * via_network_cost (default 0.5), not plain pmin().
#
# Input:   points_df (id_col, longitude, latitude, nearest_node,
#          geodist_H1_span), nodes (data/nodes.csv), network_edges (sf, from
#          [3]'s build_network_edges()), land_sf (land mask).
# Output:  symmetric, zero-diagonal, unit-named matrix (km).
# =============================================================================

library(dplyr)
library(sf)

# igraph:: kept namespaced rather than library()'d — attaching it adds a third
# union()/setdiff() candidate that breaks [5]'s conflicts_prefer() setup when
# [3] and [5] run in the same session.

# Own copy of [3]'s land/sea split so this file sources standalone.
.pairwise_split_and_penalize <- function(geom, land_sf, land_penalty) {
  land_part <- st_intersection(geom, st_geometry(land_sf))
  sea_part  <- st_difference(geom, st_geometry(land_sf))

  land_part <- land_part[!st_is_empty(land_part)]
  sea_part  <- sea_part[!st_is_empty(sea_part)]

  land_len <- if (length(land_part) > 0) as.numeric(sum(st_length(land_part))) else 0
  sea_len  <- if (length(sea_part)  > 0) as.numeric(sum(st_length(sea_part)))  else 0

  list(weighted_cost = land_len * land_penalty + sea_len)
}

# Direct penalized cost (metres); zero-length pairs short-circuit to 0.
.pairwise_direct_cost <- function(lon1, lat1, lon2, lat2, land_sf, land_penalty) {
  if (isTRUE(all.equal(c(lon1, lat1), c(lon2, lat2)))) return(0)
  geom <- st_sfc(st_linestring(rbind(c(lon1, lat1), c(lon2, lat2))), crs = 4326)
  .pairwise_split_and_penalize(geom, land_sf, land_penalty)$weighted_cost
}

#' Build a symmetric pairwise terrain-penalized migration-cost matrix.
#'
#' @param points_df    One row per unit; must contain `id_col`, `longitude`,
#'                      `latitude`, `nearest_node`, `geodist_H1_span` (the
#'                      unit -> nearest network node connector cost, km).
#' @param id_col        Name of the unique identifier column in `points_df`
#'                      (e.g. "language" or "population").
#' @param nodes          Waypoint network nodes (data/nodes.csv), `id` as character.
#' @param network_edges  sf object of network edges with a `weighted_cost` column
#'                      (from [3]'s build_network_edges()).
#' @param land_sf        Land mask (Philippines + Malaysia).
#' @param land_penalty   Terrain penalty applied to any land-crossing segment.
#' @param direct_ratio   Direct line wins only when direct_cost <= direct_ratio *
#'                      via_network_cost. Default 0.5.
#' @return A symmetric, zero-diagonal matrix (km), dimnamed by `points_df[[id_col]]`.
build_pairwise_distance_matrix <- function(points_df, id_col, nodes, network_edges, land_sf,
                                            land_penalty = 4.44, direct_ratio = 0.5) {

  stopifnot(
    "direct_ratio must be in (0, 1]" = direct_ratio > 0 && direct_ratio <= 1,
    "id_col must be unique in points_df" = !anyDuplicated(points_df[[id_col]])
  )

  ids <- points_df[[id_col]]

  # All-pairs node-to-node cost, one igraph call over the shared waypoint mesh.
  g <- igraph::graph_from_data_frame(
    d = network_edges |> st_drop_geometry() |> select(from, to, weighted_cost),
    directed = FALSE, vertices = nodes
  )
  node_dist_m <- igraph::distances(g, weights = igraph::E(g)$weighted_cost)

  pairs <- as.data.frame(t(combn(ids, 2)), stringsAsFactors = FALSE) |>
    setNames(c("id1", "id2")) |>
    as_tibble() |>
    left_join(points_df |>
                select(id1 = all_of(id_col), node1 = nearest_node,
                       conn1 = geodist_H1_span, lon1 = longitude, lat1 = latitude) |>
                mutate(node1 = as.character(node1)),
              by = "id1") |>
    left_join(points_df |>
                select(id2 = all_of(id_col), node2 = nearest_node,
                       conn2 = geodist_H1_span, lon2 = longitude, lat2 = latitude) |>
                mutate(node2 = as.character(node2)),
              by = "id2")

  pairs <- pairs |>
    rowwise() |>
    mutate(
      node_to_node_km  = node_dist_m[node1, node2] / 1000,
      via_network_cost = conn1 + node_to_node_km + conn2,
      direct_cost      = .pairwise_direct_cost(lon1, lat1, lon2, lat2, land_sf, land_penalty) / 1000,
      used_direct      = direct_cost <= direct_ratio * via_network_cost,
      pair_cost        = if_else(used_direct, direct_cost, via_network_cost)
    ) |>
    ungroup()

  message(sprintf(
    "pairwise distance: %d of %d pairs (%.0f%%) used the direct line (<= %.0f%% of the network route).",
    sum(pairs$used_direct), nrow(pairs), 100 * mean(pairs$used_direct), 100 * direct_ratio
  ))

  mat <- matrix(0, nrow = length(ids), ncol = length(ids), dimnames = list(ids, ids))
  mat[cbind(pairs$id1, pairs$id2)] <- pairs$pair_cost
  mat[cbind(pairs$id2, pairs$id1)] <- pairs$pair_cost
  mat
}
