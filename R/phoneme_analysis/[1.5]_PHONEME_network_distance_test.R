# ── network_distance_pipeline.R ──────────────────────────────────────────────
# Replaces compute_shortest_path_df(). Instead of routing each language to a
# single reference point (ref_coords1) via Dijkstra, this computes each
# language's terrain-penalized cost to ENTER the navigable network — i.e.
# the connector from the language's coordinate to its nearest graph node.
#
# Dropped vs. the old version:
#   - ref_coords1 / find-nearest-node-to-Manila
#   - shortest_path_trace() (Dijkstra) — no longer needed; "distance to the
#     network" terminates the moment a language reaches any node
#   - per-row plot_path() — replaced by one overview arrow plot
#
# Kept:
#   - land/sea terrain penalty (land_penalty), now factored into one helper
#     used by both the network edges AND the language connectors, removing
#     the 3x duplicated land/sea-split logic from the original script
#   - the final arrow plot (full_tree_sf, main vs. connector, directional
#     arrowheads), saved as an .rds exactly as before
# ──────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(geosphere)
library(sf)
library(maps)
library(sfheaders)
library(here)

# ── 0. Load data ─────────────────────────────────────────────────────────────
PHONEME_cossim <- read.csv(here("data", "PHONEME_cossim.csv"))

nodes <- read.csv(here("data", "nodes.csv")) |> mutate(id = as.character(id))
edges <- read.csv(here("data", "edges.csv")) |>
  mutate(from = as.character(from), to = as.character(to))


# ── 1. Land mask (Philippines + Malaysia) ────────────────────────────────────
world_map <- map_data("world") |> filter(region %in% c("Philippines", "Malaysia"))

land_sf <- sf_polygon(obj = world_map, polygon_id = "group", x = "long", y = "lat") |>
  st_union() |>
  st_sf(geometry = _) |>
  st_set_crs(4326)


# ── 2. Shared helper: split a linestring into land/sea + terrain-penalize ────
# Used for BOTH network edges and language→node connectors, so the land/sea
# logic only exists once instead of being copy-pasted per use case.
split_and_penalize <- function(geom, land_sf, land_penalty) {
  
  land_part <- st_intersection(geom, st_geometry(land_sf))
  sea_part  <- st_difference(geom, st_geometry(land_sf))
  
  # st_intersection/st_difference at the bare sfc level RETAIN empty results
  # (unlike the sf data.frame method, which drops them) — filter explicitly
  # rather than relying on length() to tell empty from non-empty.
  land_part <- land_part[!st_is_empty(land_part)]
  sea_part  <- sea_part[!st_is_empty(sea_part)]
  
  land_len <- if (length(land_part) > 0) as.numeric(sum(st_length(land_part))) else 0
  sea_len  <- if (length(sea_part)  > 0) as.numeric(sum(st_length(sea_part)))  else 0
  
  list(
    land_len      = land_len,
    sea_len       = sea_len,
    crosses_land  = land_len > 0,
    weighted_cost = land_len * land_penalty + sea_len,
    land_geom     = if (length(land_part) > 0) land_part else NULL,  # sfc (may hold >1 piece)
    sea_geom      = if (length(sea_part)  > 0) sea_part  else NULL   # sfc (may hold >1 piece)
  )
}


# ── 3. Build the network edges (the "main" path) ─────────────────────────────
# NOTE: geometry is built and converted to a real sfc column (via st_as_sf)
# in its own stage BEFORE the land/sea split. Under rowwise(), a column
# built inline as list(st_linestring(...)) is already unwrapped to the bare
# sfg when referenced later in the SAME mutate() call — so `geometry[[1]]`
# would index into the geometry's underlying matrix instead of extracting
# it, handing split_and_penalize() a raw numeric. Once geometry is a proper
# sfc (post st_as_sf), a second rowwise() stage extracts each row's geometry
# correctly.
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
    mutate(
      .split        = list(split_and_penalize(geometry, land_sf, land_penalty)),
      land_len      = .split$land_len,
      sea_len       = .split$sea_len,
      crosses_land  = .split$crosses_land,
      weighted_cost = .split$weighted_cost,
      land_geom     = list(.split$land_geom),
      sea_geom      = list(.split$sea_geom)
    ) |>
    ungroup() |>
    select(-.split) |>
    st_as_sf()
}


# ── 4. Nearest-node finder ────────────────────────────────────────────────────
find_nearest_node <- function(coords, nodes) {
  d <- distHaversine(matrix(c(nodes$longitude, nodes$latitude), ncol = 2), coords)
  nodes$id[which.min(d)]
}


# ── 5. Per-language: terrain-penalized distance to the network ───────────────
compute_network_distance_df <- function(df, nodes, edges, land_sf, land_penalty = 4.44) {
  
  stopifnot(land_penalty >= 1)
  
  nodes <- nodes |> mutate(id = as.character(id))
  edges <- edges |> mutate(from = as.character(from), to = as.character(to))
  
  df_geom <- df |>
    rowwise() |>
    mutate(
      start_coords  = list(c(longitude, latitude)),
      nearest_node  = find_nearest_node(start_coords, nodes),
      connector_geom = list(st_linestring(rbind(
        start_coords,
        c(nodes$longitude[nodes$id == nearest_node],
          nodes$latitude[nodes$id == nearest_node])
      )))
    ) |>
    ungroup() |>
    mutate(connector_geom = st_sfc(connector_geom, crs = 4326))   # list-col -> real sfc
  
  df_geom |>
    rowwise() |>
    mutate(
      .split          = list(split_and_penalize(connector_geom, land_sf, land_penalty)),
      land_len        = .split$land_len,
      sea_len         = .split$sea_len,
      crosses_land    = .split$crosses_land,
      geodist_H1_span = .split$weighted_cost / 1000,     # metres -> km
      land_geom       = list(.split$land_geom),
      sea_geom        = list(.split$sea_geom)
    ) |>
    ungroup() |>
    select(-.split)
}


# ── 6. Run ─────────────────────────────────────────────────────────────────
edge_lines <- build_network_edges(nodes, edges, land_sf, land_penalty = 44.18)

PHONEME_cossim <- compute_network_distance_df(
  PHONEME_cossim, nodes, edges, land_sf, land_penalty = 44.18
)

PHONEME_cossim |>
  select(any_of(c("language", "nearest_node", "geodist_H1_span",
                  "land_len", "sea_len", "crosses_land"))) |>
  arrange(geodist_H1_span) |>
  print(n = 20)


# ── 7. Assemble full_tree_sf (main edges + per-language connectors) ─────────
# Mirrors the original land_segments_sf / sea_segments_sf / connector_sf
# assembly, but built once via map()/list_rbind() instead of manual loops,
# and with a single connector per language (no "end" connector, since there
# is no longer a ref_coords1 target).
geom_rows <- function(geom_list, crosses_land_value, source_value) {
  # NOTE: purrr::map is explicitly namespaced — library(maps) is loaded after
  # library(tidyverse), so the unqualified map() resolves to maps::map()
  # (cartographic plotting), not purrr::map() (list mapping), and silently
  # fails with a cryptic "database type not supported" error.
  #
  # to_sfg_list() normalizes whatever split_and_penalize() stored for a row
  # — a bare sfg, a (possibly multi-piece) sfc, or NULL — into a flat list of
  # bare sfg objects, which is the only thing st_sfc() accepts. as.list() on
  # an sfc is the one operation sf guarantees always drops correctly to sfg;
  # relying on `[[1]]` further upstream wasn't safe across this rowwise chain.
  to_sfg_list <- function(x) {
    if (is.null(x))         return(list())
    if (inherits(x, "sfc")) return(as.list(x))
    if (inherits(x, "sfg")) return(list(x))
    list()
  }
  
  unwrapped <- purrr::map(geom_list, 1)                     # undo rowwise's list(x) wrapper
  sfg_list  <- do.call(c, purrr::map(unwrapped, to_sfg_list))  # flatten to one list of sfg
  
  if (length(sfg_list) == 0) return(NULL)
  st_sf(
    geometry     = st_sfc(sfg_list, crs = 4326),
    crosses_land = crosses_land_value,
    source       = source_value
  )
}

main_sf <- bind_rows(
  geom_rows(edge_lines$land_geom, TRUE,  "main"),
  geom_rows(edge_lines$sea_geom,  FALSE, "main")
)

connector_sf <- bind_rows(
  geom_rows(PHONEME_cossim$land_geom, TRUE,  "connector"),
  geom_rows(PHONEME_cossim$sea_geom,  FALSE, "connector")
)

full_tree_sf <- bind_rows(main_sf, connector_sf)


# ── 8. Simplify each main edge to a straight start→end segment, for arrows ──
# (Land/sea-split edges can be multi-vertex; arrows look cleaner on the
# overall edge direction rather than every sub-segment.)
simplify_to_straight <- function(sf_lines) {
  sf_lines |>
    mutate(.row = row_number()) |>
    group_by(.row) |>
    group_modify(~ {
      coords <- st_coordinates(.x)
      if (nrow(coords) < 2) return(.x)
      new_geom <- st_sfc(
        st_linestring(rbind(coords[1, c("X", "Y")], coords[nrow(coords), c("X", "Y")])),
        crs = st_crs(.x)
      )
      st_set_geometry(.x, new_geom)
    }) |>
    ungroup() |>
    select(-.row) |>
    st_as_sf()
}

arrow_main <- full_tree_sf |>
  filter(source == "main") |>
  simplify_to_straight()

arrow_connectors <- full_tree_sf |>
  filter(source == "connector")


# ── 9. Overview plot: network + language points + directional arrows ────────
plot_network <- function(full_tree_sf, arrow_main, arrow_connectors,
                         points_df, lon_range = c(116, 127), lat_range = c(4, 21)) {
  ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
                 fill = "gray95", color = "gray70") +
    geom_sf(data = full_tree_sf, linewidth = 1, color = "black") +
    geom_sf(data = arrow_connectors,
            arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
            color = "black") +
    geom_sf(data = arrow_main,
            arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
            color = "black") +
    geom_point(data = points_df, aes(x = longitude, y = latitude),
               size = 3, shape = 21) +
    coord_sf(xlim = lon_range, ylim = lat_range) +
    theme_minimal() +
    labs(title = "Phonemic Historical Routes Waypoint System",
         x = "Longitude", y = "Latitude")
}

print(plot_network(full_tree_sf, arrow_main, arrow_connectors, PHONEME_cossim))


# ── 10. Arrow-only plot (connectors get arrows, main path stays plain) ──────
arrow_plot <- ggplot() +
  geom_sf(data = full_tree_sf, color = "black", linewidth = 1) +
  geom_sf(data = arrow_connectors,
          arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
          color = "black", linewidth = 1) +
  coord_sf(xlim = c(116, 127), ylim = c(4, 21)) +
  theme_minimal()

print(arrow_plot)
saveRDS(arrow_plot, file = here("data", "phoneme_waypoint_plot.rds"))